#!/usr/bin/env python3
"""
TackleBox: FlyForge v1.2.0

Comprehensive bait/probe designer for hybridization capture enrichment
with integrated support for in-house RNA probe synthesis via oligo pools.

FlyForge designs custom bait panels for single-species targets or large
metagenomic projects with thousands of reference genomes. It combines
flexible tiling and multi-stage quality filtering with oligo pool design
for cost-effective in-house probe synthesis.

Features:
    - Flexible tiling with iterative density adjustment to hit target bait count
    - Preprocessing (U->T, N-run replacement, self-repeat softmasking)
    - Complementarity removal in input sequences (CARPDM-style)
    - Multi-stage filtering (ambiguous bases, repeat-masked, Tm, LguI sites)
    - BLAST specificity filtering against user-supplied databases
    - Self-BLAST redundancy and complementarity collapsing
    - cd-hit-est clustering
    - O-pool design: T7 promoter + amplification primer + BspQI/LguI cut site
      (following CARPDM / Hackenberger et al. 2025)
    - BLAST validation against input targets with coverage analysis
    - Comprehensive plots and summary statistics
    - Outputs both bare probes and full oligo pool sequences for ordering

Oligo pool synthesis approach adapted from:
    Hackenberger et al. 2025. CARPDM: cost-effective antibiotic resistome
    profiling of metagenomic samples using targeted enrichment.
    Applied and Environmental Microbiology 91(3):e01876-24.

Author:  Tyler J. Murchie
Created with assistance from Claude (Anthropic) and ChatGPT 5.4 (OpenAI).
License: MIT
"""

import argparse
import os
import re
import sys
import tempfile
import subprocess
import time
import shlex
import shutil
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Set, Optional
from hashlib import sha256
from io import StringIO
from itertools import combinations_with_replacement, permutations

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO, SeqUtils
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Blast import NCBIXML
import primer3 as primer3_mod
from tqdm import tqdm

__version__ = "1.2.0"

AMBIGUOUS_BASES = set("NRYWSMKHBVDnrywsmkhbvd")
AMBIGUOUS_RE = re.compile(r'[^ATGC]')

T7_PROMOTER = "GCTAATACGACTCACTATAGGG"
T7_PROMOTER_LEN = len(T7_PROMOTER)
BSPQI_SITE_FWD = "GAAGAGC"
BSPQI_SITE_REV = "GCTCTTC"
DEFAULT_OPRIMER_SUFFIX = "GCTCTTCG"

# ============================================================================
# Data structures
# ============================================================================

@dataclass
class Bait:
    bait_id: str
    seq: str
    ref_id: str
    ref_start: int   # 1-based inclusive
    ref_end: int      # 1-based inclusive
    gc_frac: float
    masked_frac: float
    ambiguous_count: int
    tm: float = 0.0


@dataclass
class RefStats:
    ref_id: str
    length_original: int
    length_preprocessed: int
    modified_bases: int = 0
    frac_modified: float = 0.0
    dropped_for_length: bool = False
    n_baits_tiled: int = 0
    n_baits_after_amb: int = 0
    n_baits_after_mask: int = 0
    n_baits_after_tm: int = 0
    n_baits_after_blast: int = 0
    n_baits_after_cluster: int = 0
    n_baits_after_redundancy: int = 0
    n_baits_final: int = 0
    coverage_mean: float = 0.0
    coverage_min: int = 0
    coverage_max: int = 0
    coverage_fraction_covered: float = 0.0
    proportion_of_final_baits: float = 0.0


# ============================================================================
# Logging & Progress Display
# ============================================================================

def format_duration(seconds: float) -> str:
    """Format seconds into human-readable duration string."""
    if seconds < 60:
        return f"{seconds:.1f}s"
    elif seconds < 3600:
        m, s = divmod(int(seconds), 60)
        return f"{m}m {s}s"
    else:
        h, remainder = divmod(int(seconds), 3600)
        m, s = divmod(remainder, 60)
        return f"{h}h {m}m {s}s"


def make_progress_bar(fraction: float, width: int = 40) -> str:
    """Create a text-based progress bar."""
    filled = int(width * fraction)
    bar = "█" * filled + "░" * (width - filled)
    pct = fraction * 100
    return f"[{bar}] {pct:5.1f}%"


def make_logger(progress_path):
    """Create a logger that writes to stderr and optionally to a file."""
    log_handle = open(progress_path, "w") if progress_path else None

    def log(msg: str):
        timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
        line = f"[{timestamp}] {msg}"
        print(line, file=sys.stderr)
        if log_handle is not None:
            print(line, file=log_handle)
            log_handle.flush()

    return log, log_handle


class ProgressTracker:
    """Visual progress tracker for the FlyForge pipeline."""

    # Display names for each step
    STEP_NAMES = {
        "preprocessing":         "Preprocessing sequences",
        "complement_removal":    "Removing complementary regions",
        "self_masking":          "Self-repeat softmasking",
        "tiling":                "Tiling baits",
        "ambiguous_filter":      "Filtering ambiguous bases",
        "masked_filter":         "Filtering repeat-masked baits",
        "tm_filter":             "Filtering by melting temperature",
        "lgui_filter":           "Filtering LguI/BspQI sites",
        "complement_bait_filter":"Removing complement bait pairs",
        "blast_filter":          "BLAST percent-identity filter",
        "specificity_filter":    "BLAST specificity filter (nt)",
        "self_blast_redundancy": "Self-BLAST redundancy filter",
        "clustering":            "cd-hit-est clustering",
        "opool_design":          "Designing oligo pool (T7 + primer)",
        "validation":            "BLAST validation & analysis",
        "write_output":          "Writing output files",
    }

    def __init__(self, steps: List[str], log_fn, t_start: float):
        self.steps = steps
        self.total = len(steps)
        self.completed = 0
        self.log_fn = log_fn
        self.t_start = t_start
        self.step_times: List[float] = []
        self.last_time = t_start
        self.status: Dict[str, str] = {}  # step -> "pending"|"running"|"done"
        for s in steps:
            self.status[s] = "pending"

    def _print_dashboard(self, current_step: str = None,
                         details: str = "", step_dt: float = 0):
        """Print the full progress dashboard to stderr."""
        now = time.time()
        elapsed = now - self.t_start
        fraction = self.completed / self.total if self.total > 0 else 0

        # ETA calculation
        if self.step_times:
            avg_step = sum(self.step_times) / len(self.step_times)
            remaining_steps = self.total - self.completed
            eta = avg_step * remaining_steps
        else:
            eta = 0

        lines = []
        lines.append("")
        lines.append("=" * 72)
        lines.append(f"  TackleBox: FlyForge — Pipeline Progress")
        lines.append("=" * 72)
        lines.append(f"  {make_progress_bar(fraction, 50)}")
        lines.append(f"  Elapsed: {format_duration(elapsed)}  |  "
                     f"ETA: {format_duration(eta) if self.completed < self.total else 'done'}")
        lines.append("-" * 72)

        # Step checklist
        for i, step in enumerate(self.steps, 1):
            display = self.STEP_NAMES.get(step, step)
            st = self.status[step]
            if st == "done":
                icon = "  [x]"
            elif st == "running":
                icon = "  [>]"
            else:
                icon = "  [ ]"
            line = f"{icon} {i:2d}. {display}"
            if st == "done" and step == current_step:
                line += f"  ({format_duration(step_dt)}"
                if details:
                    line += f" — {details}"
                line += ")"
            lines.append(line)

        lines.append("-" * 72)
        lines.append(f"  Steps completed: {self.completed}/{self.total}")
        lines.append("=" * 72)
        lines.append("")

        print("\n".join(lines), file=sys.stderr)

    def start_step(self, step: str):
        """Mark a step as running and log it."""
        self.status[step] = "running"
        display = self.STEP_NAMES.get(step, step)
        self.log_fn(f"[Step {self.completed + 1}/{self.total}] "
                    f"Starting: {display}...")

    def finish_step(self, step: str, details: str = ""):
        """Mark a step as done and print the dashboard."""
        now = time.time()
        dt = now - self.last_time
        self.last_time = now
        self.step_times.append(dt)
        self.completed += 1
        self.status[step] = "done"

        display = self.STEP_NAMES.get(step, step)
        elapsed = now - self.t_start
        msg = (f"[Step {self.completed}/{self.total}] "
               f"Finished: {display} ({format_duration(dt)})")
        if details:
            msg += f" — {details}"
        msg += f"  [Total elapsed: {format_duration(elapsed)}]"
        self.log_fn(msg)

        self._print_dashboard(current_step=step, details=details,
                              step_dt=dt)


# ============================================================================
# FASTA I/O helpers
# ============================================================================

def read_fasta(path: str) -> Dict[str, str]:
    """Read a FASTA file and return dict of {id: sequence}."""
    seqs: Dict[str, str] = {}
    current_id = None
    chunks: List[str] = []
    with (sys.stdin if path == "-" else open(path, "r")) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    seqs[current_id] = "".join(chunks)
                current_id = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if current_id is not None:
            seqs[current_id] = "".join(chunks)
    return seqs


def write_fasta(path: str, records: List[Tuple[str, str]], wrap: int = 80) -> None:
    """Write a FASTA file from list of (id, sequence) tuples."""
    with open(path, "w") as out:
        for rid, seq in records:
            out.write(f">{rid}\n")
            for i in range(0, len(seq), wrap):
                out.write(seq[i:i+wrap] + "\n")


def write_seqrecords_fasta(records: List[SeqRecord], path: str) -> None:
    """Write SeqRecord list to FASTA."""
    SeqIO.write(records, path, 'fasta')


# ============================================================================
# Preprocessing
# ============================================================================

def preprocess_sequence(seq: str) -> Tuple[str, int]:
    """
    Preprocess sequence:
    - uppercase
    - U->T
    - N-runs (1-10) replaced with T-runs
    Returns (processed_seq, n_modified_bases)
    """
    seq_proc = seq.upper().replace("U", "T")

    def repl(m):
        return "T" * len(m.group(0))

    seq_proc = re.sub(r"N{1,10}", repl, seq_proc)
    n_modified = sum(
        1 for a, b in zip(seq.upper(), seq_proc) if a != b
    )
    return seq_proc, n_modified


def self_repeat_softmask(seqs: Dict[str, str], k: int = 15, threshold: int = 3) -> Dict[str, str]:
    """
    Simple self-repeat masker:
    - Count k-mers across all sequences.
    - Any k-mer occurring >= threshold times is considered repetitive.
    - All positions covered by such k-mers are converted to lowercase.
    """
    kmer_counts = Counter()
    for sid, seq in seqs.items():
        s = seq.upper()
        for i in range(0, len(s) - k + 1):
            kmer = s[i:i+k]
            if "N" in kmer:
                continue
            kmer_counts[kmer] += 1

    repeat_kmers = {km for km, c in kmer_counts.items() if c >= threshold}

    masked: Dict[str, str] = {}
    for sid, seq in seqs.items():
        s = list(seq.upper())
        mask = [False] * len(s)
        for i in range(0, len(s) - k + 1):
            kmer = "".join(s[i:i+k])
            if kmer in repeat_kmers:
                for j in range(i, i + k):
                    mask[j] = True
        for i, m in enumerate(mask):
            if m:
                s[i] = s[i].lower()
        masked[sid] = "".join(s)
    return masked


# ============================================================================
# Complementary region removal (from CARPDM)
# ============================================================================

def collapse_slices(slice_lst: List[Tuple[int, int]],
                    result_lst: List = None) -> List[Tuple[int, int]]:
    """Recursively collapse overlapping slices."""
    if result_lst is None:
        result_lst = []
    if len(slice_lst) == 1:
        result_lst.append(slice_lst[0])
        return result_lst
    elif slice_lst[0][-1] >= slice_lst[1][0] and len(slice_lst) == 2:
        result_lst.append((slice_lst[0][0], max(slice_lst[0][1], slice_lst[1][1])))
        return result_lst
    elif slice_lst[0][-1] >= slice_lst[1][0]:
        result_lst.append((slice_lst[0][0], max(slice_lst[0][1], slice_lst[1][1])))
        return collapse_slices(slice_lst[2:], result_lst)
    else:
        result_lst.append(slice_lst[0])
        return collapse_slices(slice_lst[1:], result_lst)


def chop_record(input_record: SeqRecord, slices_to_rc: List[Tuple[int, int]],
                probe_length: int) -> List[SeqRecord]:
    """Chop a record to remove complementary regions."""
    all_slices = [position for s in slices_to_rc for position in s]
    all_slices.insert(0, 0)
    all_slices.append(len(input_record.seq) + 1)
    slices_to_keep = [(max(0, all_slices[i]), all_slices[i + 1])
                      for i in range(0, len(all_slices), 2)]
    records = []
    for s in slices_to_keep:
        working_record = input_record[s[0]:s[1]]
        if len(working_record.seq) >= probe_length:
            working_record.id = working_record.id + f'_[{s[0]}:{s[1]}]'
            working_record.description = ''
            records.append(working_record)
    for s in slices_to_rc:
        working_record = input_record[s[0]:s[1]]
        if len(working_record.seq) >= probe_length:
            working_record.id = working_record.id + f'_[{s[0]}:{s[1]}]'
            working_record.description = ''
            working_record.seq = working_record.seq.reverse_complement()
            records.append(working_record)
    return records


def cut_seqs(input_fasta: List[SeqRecord],
             slice_dict: Dict[str, Set[Tuple[int, int]]],
             probe_length: int) -> List[SeqRecord]:
    """Return list of SeqRecords with complementary regions removed."""
    slice_seqs = slice_dict.keys()
    new_fasta = []
    for record in input_fasta:
        if record.id in slice_seqs:
            slices = slice_dict[record.id]
            old_slice_lst = sorted(list(slices), key=lambda x: (x[0], x[1]))
            if len(old_slice_lst) > 1:
                while True:
                    new_slice_lst = collapse_slices(old_slice_lst)
                    if new_slice_lst == old_slice_lst or len(new_slice_lst) == 1:
                        final_slice_lst = new_slice_lst
                        break
                    old_slice_lst = new_slice_lst
            else:
                final_slice_lst = old_slice_lst
            record_lst = chop_record(record, final_slice_lst, probe_length)
            new_fasta.extend(record_lst)
        else:
            new_fasta.append(record)
    return new_fasta


def remove_complementary_targets(input_fasta: List[SeqRecord],
                                 output_dir: str, prefix: str,
                                 num_threads: int, probe_length: int,
                                 log_fn) -> str:
    """BLAST input against itself and remove complementary regions.
    Returns path to output FASTA."""
    header = ('qseqid sseqid length pident nident staxid sstrand '
              'qstart qend sstart send qlen slen').split()
    output_path = f'{prefix}_input_no_comp.fna'

    while True:
        str_fasta = ''.join([record.format('fasta') for record in input_fasta])
        database = f'{output_dir}/self_db_temp'
        cmd = shlex.split(
            f'makeblastdb -out {database} -dbtype nucl -title self_db_temp')
        subprocess.run(cmd, input=str_fasta, text=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        outfmt = f'6 {" ".join(header)}'
        cmd = shlex.split(
            f'blastn -outfmt "{outfmt}" -num_threads {num_threads} '
            f'-db {database} -strand minus')
        blast_output = subprocess.run(cmd, text=True, input=str_fasta,
                                      capture_output=True).stdout
        if not blast_output.strip():
            SeqIO.write(input_fasta, output_path, 'fasta')
            break

        self_blast_df = pd.read_csv(StringIO(blast_output), sep='\t', names=header)
        working_df = self_blast_df.loc[self_blast_df['nident'] >= 30]
        if working_df.empty:
            SeqIO.write(input_fasta, output_path, 'fasta')
            break

        top_hit = working_df['qseqid'].mode().squeeze()
        if isinstance(top_hit, pd.Series):
            top_hit = top_hit.sort_values(ascending=False).iloc[0]
        top_hit_df = working_df.loc[working_df['qseqid'] == top_hit]

        slice_dict = defaultdict(set)
        for t in top_hit_df[['sseqid', 'sstart', 'send']].itertuples():
            slice_dict[t[1]].add((min(t[2], t[3]), max(t[2], t[3])))
        input_fasta = cut_seqs(input_fasta, slice_dict, probe_length)

    # Clean up temp db files
    for f in os.listdir(output_dir):
        if f.startswith('self_db_temp'):
            os.remove(os.path.join(output_dir, f))

    log_fn(f"Complementary region removal complete. "
           f"{len(input_fasta)} segments remain.")
    return output_path


# ============================================================================
# Bait metrics computation
# ============================================================================

def compute_bait_metrics(seq: str) -> Tuple[float, float, int]:
    """Compute GC fraction, masked fraction, and ambiguous base count."""
    length = len(seq)
    if length == 0:
        return 0.0, 0.0, 0
    upper = seq.upper()
    gc = upper.count("G") + upper.count("C")
    gc_frac = gc / length
    masked = sum(1 for b in seq if b.islower())
    masked_frac = masked / length
    ambiguous = sum(1 for b in seq if b in AMBIGUOUS_BASES)
    return gc_frac, masked_frac, ambiguous


def compute_tm(seq: str) -> float:
    """Compute melting temperature using nearest-neighbor for RNA/DNA."""
    try:
        return mt.Tm_NN(Seq(seq.upper()), nn_table=mt.R_DNA_NN1)
    except Exception:
        return 0.0


# ============================================================================
# Tiling
# ============================================================================

def tile_sequence(ref_id: str, seq: str, bait_len: int = 80,
                  tiling_density: float = 3.0, omit_short_leq: int = 70,
                  pad_min: int = 71, circular: bool = False) -> List[Bait]:
    """Generate tiled baits from a (softmasked) reference sequence.

    When ``circular`` is True, windows are allowed to wrap across the
    end/start boundary so circular genomes (for example mitochondria) do not
    lose enrichment at the linearized termini.
    """
    baits: List[Bait] = []
    L = len(seq)
    if L <= omit_short_leq:
        return baits

    if L < bait_len and L >= pad_min:
        padded = seq + "T" * (bait_len - L)
        gc_frac, masked_frac, ambiguous = compute_bait_metrics(padded)
        bait_id = f"{ref_id}|b1|pos1-{bait_len}"
        baits.append(Bait(
            bait_id=bait_id, seq=padded, ref_id=ref_id,
            ref_start=1, ref_end=L, gc_frac=gc_frac,
            masked_frac=masked_frac, ambiguous_count=ambiguous,
            tm=compute_tm(padded)))
        return baits

    if L < bait_len:
        return baits

    step = max(1, int(round(bait_len / tiling_density)))
    if circular:
        starts = list(range(0, L, step))
    else:
        starts = list(range(0, L - bait_len + 1, step))
        last_start = L - bait_len
        if starts[-1] != last_start:
            starts.append(last_start)
        starts = sorted(set(starts))

    for idx, start in enumerate(starts, start=1):
        if circular and start + bait_len > L:
            wrap = (start + bait_len) - L
            frag = seq[start:] + seq[:wrap]
        else:
            frag = seq[start:start + bait_len]
        gc_frac, masked_frac, ambiguous = compute_bait_metrics(frag)
        bait_id = f"{ref_id}|b{idx}|pos{start+1}-{start + bait_len}"
        baits.append(Bait(
            bait_id=bait_id, seq=frag, ref_id=ref_id,
            ref_start=start + 1, ref_end=start + bait_len,
            gc_frac=gc_frac, masked_frac=masked_frac,
            ambiguous_count=ambiguous, tm=compute_tm(frag)))
    return baits


def estimate_baits_for_density(lengths: Dict[str, int], bait_len: int,
                               omit_short_leq: int, pad_min: int,
                               tiling_density: float) -> int:
    """Estimate total bait count for a given tiling density (pre-filter)."""
    total = 0
    for L in lengths.values():
        if L <= omit_short_leq:
            continue
        if pad_min <= L < bait_len:
            total += 1
        elif L >= bait_len:
            step = max(1, int(round(bait_len / tiling_density)))
            n_starts = ((L - bait_len) // step) + 1
            last_start = (n_starts - 1) * step
            if last_start != L - bait_len:
                n_starts += 1
            total += n_starts
    return total


def choose_tiling_density(lengths: Dict[str, int], bait_len: int,
                          omit_short_leq: int, pad_min: int,
                          max_baits: int, min_density: float,
                          requested_density: float) -> Tuple[float, int]:
    """Find tiling density giving <= max_baits (estimated), or use min_density."""
    if max_baits <= 0:
        return requested_density, estimate_baits_for_density(
            lengths, bait_len, omit_short_leq, pad_min, requested_density)

    est_req = estimate_baits_for_density(
        lengths, bait_len, omit_short_leq, pad_min, requested_density)
    if est_req <= max_baits:
        return requested_density, est_req

    est_min = estimate_baits_for_density(
        lengths, bait_len, omit_short_leq, pad_min, min_density)
    if est_min > max_baits:
        return min_density, est_min

    low = min_density
    high = requested_density
    best_d = min_density
    best_est = est_min

    for _ in range(25):
        mid = 0.5 * (low + high)
        est_mid = estimate_baits_for_density(
            lengths, bait_len, omit_short_leq, pad_min, mid)
        if est_mid <= max_baits:
            best_d = mid
            best_est = est_mid
            low = mid
        else:
            high = mid
        if high - low < 1e-3:
            break

    return best_d, best_est


# ============================================================================
# Filtering
# ============================================================================

def filter_ambiguous(baits: List[Bait], ambiguous_cutoff: int = 10) -> Tuple[List[Bait], int]:
    """Remove baits with >= ambiguous_cutoff ambiguous bases."""
    kept = [b for b in baits if b.ambiguous_count < ambiguous_cutoff]
    return kept, len(baits) - len(kept)


def filter_masked_fraction(baits: List[Bait], max_masked_frac: float = 0.25) -> Tuple[List[Bait], int]:
    """Remove baits with masked fraction > max_masked_frac."""
    kept = [b for b in baits if b.masked_frac <= max_masked_frac]
    return kept, len(baits) - len(kept)


def filter_melting_temp(baits: List[Bait], min_tm: float = 50.0) -> Tuple[List[Bait], int]:
    """Remove baits with Tm < min_tm."""
    kept = [b for b in baits if b.tm >= min_tm]
    return kept, len(baits) - len(kept)


def filter_lgui_sites(baits: List[Bait]) -> Tuple[List[Bait], int]:
    """Remove baits containing LguI/BspQI cut sites (GAAGAGC or GCTCTTC)
    which would interfere with oligo pool synthesis."""
    kept = []
    for b in baits:
        s = b.seq.upper()
        if 'GAAGAGC' not in s and 'GCTCTTC' not in s:
            kept.append(b)
    return kept, len(baits) - len(kept)


def filter_complementary_baits(baits: List[Bait]) -> Tuple[List[Bait], int]:
    """Remove baits whose reverse complement is already in the set."""
    seq_set = {b.seq.upper() for b in baits}
    kept = []
    removed_ids = set()
    for b in baits:
        rc = str(Seq(b.seq.upper()).reverse_complement())
        if rc in seq_set and b.bait_id not in removed_ids:
            # Keep this one, but mark its complement for removal if found later
            # Actually: remove one of each pair. Use sorted order.
            pass
        kept.append(b)

    # Better approach: iterate sorted, remove if RC exists and hasn't been kept
    kept = []
    seen_seqs = set()
    removed = 0
    for b in sorted(baits, key=lambda x: x.bait_id):
        s = b.seq.upper()
        rc = str(Seq(s).reverse_complement())
        if rc in seen_seqs:
            removed += 1
        else:
            kept.append(b)
            seen_seqs.add(s)
    return kept, removed


# ============================================================================
# BLAST-based filtering
# ============================================================================

def blast_filter_baits(baits: List[Bait], blast_db: str,
                       evalue: float = 1e-5, min_pident: float = 90.0,
                       max_hits_per_query: int = 1, threads: int = 1,
                       ) -> Tuple[List[Bait], int]:
    """Filter baits by BLASTn hits against a user-supplied DB.
    Baits with any hit meeting thresholds are removed."""
    if not baits:
        return baits, 0

    with tempfile.TemporaryDirectory() as tmpdir:
        query_fa = os.path.join(tmpdir, "baits.fasta")
        write_fasta(query_fa, [(b.bait_id, b.seq) for b in baits])

        out_tsv = os.path.join(tmpdir, "blast.tsv")
        cmd = [
            "blastn", "-task", "blastn", "-query", query_fa,
            "-db", blast_db,
            "-outfmt", "6 qseqid pident length evalue bitscore",
            "-evalue", str(evalue),
            "-max_target_seqs", str(max_hits_per_query),
            "-num_threads", str(threads),
            "-dust", "no", "-soft_masking", "false",
        ]

        # Set BLASTDB environment for nt database
        env = os.environ.copy()
        db_dir = os.path.dirname(blast_db)
        if db_dir:
            env['BLASTDB'] = db_dir

        with open(out_tsv, "w") as out:
            subprocess.run(cmd, check=True, stdout=out, stderr=subprocess.DEVNULL,
                           env=env)

        bad_ids = set()
        with open(out_tsv) as fh:
            for line in fh:
                parts = line.strip().split()
                if len(parts) >= 5:
                    qid, pident_val = parts[0], float(parts[1])
                    if pident_val >= min_pident:
                        bad_ids.add(qid)

    kept = [b for b in baits if b.bait_id not in bad_ids]
    return kept, len(baits) - len(kept)


def blast_specificity_filter(baits: List[Bait], blast_db: str,
                             probe_length: int, threads: int = 1,
                             is_nt: bool = False) -> Tuple[List[Bait], int]:
    """CARPDM-style identity-based specificity filter.
    Uses nident (number of identical residues) rather than percent identity.
    For nt database: only removes non-bacterial hits within identity range.
    For other DBs: removes any hit with > probe_length * 0.625 identities."""
    if not baits:
        return baits, 0

    with tempfile.TemporaryDirectory() as tmpdir:
        query_fa = os.path.join(tmpdir, "baits.fasta")
        write_fasta(query_fa, [(b.bait_id, b.seq.upper()) for b in baits])

        out_tsv = os.path.join(tmpdir, "blast.tsv")
        env = os.environ.copy()
        db_dir = os.path.dirname(blast_db)
        if db_dir:
            env['BLASTDB'] = db_dir
        db_name = os.path.basename(blast_db)

        if is_nt:
            header = 'qseqid sseqid length nident staxid sstrand sskingdom'
        else:
            header = 'qseqid sseqid length nident'

        cmd = shlex.split(
            f'blastn -query {query_fa} -out {out_tsv} '
            f'-outfmt "6 {header}" -num_threads {threads} '
            f'-db {db_name}')
        subprocess.run(cmd, env=env, stderr=subprocess.DEVNULL)

        header_list = header.split()
        try:
            blast_df = pd.read_csv(out_tsv, sep='\t', names=header_list)
        except Exception:
            return baits, 0

        if blast_df.empty:
            return baits, 0

        cutoff = probe_length * 0.625
        if is_nt:
            non_bac_high = probe_length * 0.975
            filter_set = {'Archaea', 'Eukaryota', 'Viruses'}
            probes_to_drop = set(blast_df.loc[
                (blast_df['sskingdom'].isin(filter_set))
                & (blast_df['nident'] < non_bac_high)
                & (blast_df['nident'] > cutoff),
                'qseqid'])
        else:
            probes_to_drop = set(blast_df.loc[
                blast_df['nident'] > cutoff, 'qseqid'])

    kept = [b for b in baits if b.bait_id not in probes_to_drop]
    return kept, len(baits) - len(kept)


# ============================================================================
# Self-BLAST redundancy and complementarity filtering (from CARPDM)
# ============================================================================

def self_blast_filter(baits: List[Bait], output_dir: str, prefix: str,
                      threads: int, probe_num_cutoff: int,
                      log_fn) -> Tuple[List[Bait], int, int]:
    """Perform self-BLAST to remove complementary and redundant probes.
    Returns (kept_baits, n_complementary_removed, n_redundancy_removed)."""
    if not baits:
        return baits, 0, 0

    # Write current baits to temp file
    temp_fa = os.path.join(output_dir, f'{prefix}_pre_self_blast.fna')
    write_fasta(temp_fa, [(b.bait_id, b.seq.upper()) for b in baits])

    # Build temp BLAST DB and run self-BLAST
    header = 'qseqid sseqid length nident sstrand'.split()
    database = os.path.join(output_dir, 'self_db_temp')
    cmd = shlex.split(
        f'makeblastdb -in {temp_fa} -out {database} -dbtype nucl')
    subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    blast_out = os.path.join(output_dir, f'{prefix}_self_blast.txt')
    outfmt = f'6 {" ".join(header)}'
    cmd = shlex.split(
        f'blastn -query {temp_fa} -out {blast_out} -outfmt "{outfmt}" '
        f'-num_threads {threads} -db {database}')
    subprocess.run(cmd, stderr=subprocess.DEVNULL)

    # Clean up temp DB
    for f in os.listdir(output_dir):
        if f.startswith('self_db_temp'):
            os.remove(os.path.join(output_dir, f))

    # Parse results
    try:
        self_blast_df = pd.read_csv(blast_out, sep='\t', names=header)
    except Exception:
        return baits, 0, 0

    if self_blast_df.empty:
        return baits, 0, 0

    # Step 1: Remove complementary probes
    rev_comp_df = self_blast_df.loc[
        (self_blast_df['qseqid'] != self_blast_df['sseqid'])
        & (self_blast_df['nident'] >= 30)
        & (self_blast_df['sstrand'] == 'minus')]

    comp_probes_to_drop = set()
    if not rev_comp_df.empty:
        query_set = set(rev_comp_df['qseqid'].to_list())
        hit_set = set(rev_comp_df['sseqid'].to_list())
        comp_probes_to_drop.update(hit_set ^ query_set)
        rev_comp_df2 = rev_comp_df.loc[
            (~rev_comp_df['qseqid'].isin(comp_probes_to_drop))
            & (~rev_comp_df['sseqid'].isin(comp_probes_to_drop))]
        hit_counts = rev_comp_df2['qseqid'].value_counts()
        queries = [k for k, v in sorted(hit_counts.items(),
                                        key=lambda item: (item[1], item[0]))]
        for query in queries:
            if query not in comp_probes_to_drop:
                hits = list(rev_comp_df2.loc[
                    rev_comp_df2['qseqid'] == query, 'sseqid'].unique())
                comp_probes_to_drop.update(hits)

    n_comp_removed = len(comp_probes_to_drop)
    log_fn(f"Self-BLAST complementarity filter: {n_comp_removed} probes removed")

    # Step 2: Redundancy filter
    probes_to_drop = set(comp_probes_to_drop)
    n_total_probes = len(baits)
    blast_df = self_blast_df.loc[
        (self_blast_df['qseqid'] != self_blast_df['sseqid'])
        & (self_blast_df['sstrand'] == 'plus')]

    for id_cutoff in reversed(range(int(baits[0].seq.__len__() if baits else 80))):
        remaining = n_total_probes - len(probes_to_drop)
        if remaining < probe_num_cutoff:
            break
        working_df = blast_df.loc[
            (blast_df['nident'] == id_cutoff)
            & (~blast_df['qseqid'].isin(probes_to_drop))
            & (~blast_df['sseqid'].isin(probes_to_drop))]
        if working_df.empty:
            continue
        hit_set = set(working_df['sseqid'].to_list())
        query_set = set(working_df['qseqid'].to_list())
        probes_to_drop.update(hit_set ^ query_set)
        working_df = working_df.loc[
            (~working_df['qseqid'].isin(probes_to_drop))
            & (~working_df['sseqid'].isin(probes_to_drop))]
        hit_counts = working_df['qseqid'].value_counts()
        queries = [k for k, v in sorted(hit_counts.items(), reverse=True,
                                        key=lambda item: (item[1], item[0]))]
        for query in queries:
            if query not in probes_to_drop:
                hits = list(working_df.loc[
                    working_df['qseqid'] == query, 'sseqid'].unique())
                probes_to_drop.update(hits)

    n_red_removed = len(probes_to_drop) - n_comp_removed
    log_fn(f"Self-BLAST redundancy filter: {n_red_removed} additional probes removed")

    # Build kept list preserving order
    kept = [b for b in baits if b.bait_id not in probes_to_drop]

    # Clean up temp file
    if os.path.exists(temp_fa):
        os.remove(temp_fa)

    return kept, n_comp_removed, n_red_removed


# ============================================================================
# cd-hit-est clustering
# ============================================================================

def cluster_baits_cd_hit(baits: List[Bait], identity: float = 0.95,
                         overlap: float = 0.83,
                         threads: int = 1) -> Tuple[List[Bait], int]:
    """Cluster baits using cd-hit-est, keeping one representative per cluster."""
    if not baits:
        return baits, 0

    with tempfile.TemporaryDirectory() as tmpdir:
        in_fa = os.path.join(tmpdir, "baits_in.fasta")
        out_fa = os.path.join(tmpdir, "baits_cdhit.fasta")
        write_fasta(in_fa, [(b.bait_id, b.seq) for b in baits])

        cmd = [
            "cd-hit-est", "-i", in_fa, "-o", out_fa,
            "-c", str(identity), "-aS", str(overlap),
            "-M", "0", "-T", str(threads), "-d", "0",
        ]

        try:
            subprocess.run(cmd, check=True,
                           stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
        except FileNotFoundError:
            sys.stderr.write("WARNING: cd-hit-est not found; skipping clustering.\n")
            return baits, 0

        kept_ids = set()
        with open(out_fa) as fh:
            for line in fh:
                if line.startswith(">"):
                    kept_ids.add(line[1:].strip().split()[0])

    id_to_bait = {b.bait_id: b for b in baits}
    kept = [id_to_bait[i] for i in kept_ids if i in id_to_bait]
    return kept, len(baits) - len(kept)


# ============================================================================
# O-Pool design: T7 promoter + primer design (from CARPDM)
# ============================================================================

def design_opool(baits: List[Bait], output_dir: str, prefix: str,
                 log_fn) -> Tuple[str, str, str, str]:
    """Design an order-ready oligo pool for in-house RNA bait synthesis.

    Returns (oligo_pool_path, bare_probes_path, primer_fasta_path, primer_seq).

    Final oligo structure:
    5'-GCTAATACGACTCACTATAGGG [probe] [reverse_complement(amplification_primer)]-3'

    The selected amplification primer always ends in ``GCTCTTCG`` so that its
    reverse complement appends ``CGAAGAGC`` at the 3' end of the oligo, yielding
    the CARPDM-style BspQI/LguI-compatible construct.
    """
    if not baits:
        raise RuntimeError("No baits available for oligo-pool design.")

    if shutil.which('blastn') is None:
        raise RuntimeError(
            'blastn was not found in PATH. BLAST+ is required for oligo-pool '
            'primer selection. Install BLAST+ or rerun with --no-opool.'
        )

    bare_probes_path = os.path.join(output_dir, f'{prefix}_probes.fna')
    sorted_baits = sorted(baits, key=lambda x: x.bait_id)
    write_fasta(bare_probes_path,
                [(b.bait_id, b.seq.upper()) for b in sorted_baits])

    log_fn('Generating candidate amplification primers...')
    primer2 = T7_PROMOTER
    primer2_tm = primer3_mod.calc_tm(primer2)
    primer_length = 18
    first_nts = ''
    last_nts = DEFAULT_OPRIMER_SUFFIX
    len_remaining = primer_length - len(first_nts) - len(last_nts)
    nts = 'ATGC'

    order_set = set()
    seqs = combinations_with_replacement(nts, len_remaining)
    for seq in seqs:
        orders = set([''.join(order) for order in permutations(seq)])
        order_set.update(orders)

    candidate_records = []
    for seq in sorted(order_set):
        primer_seq = first_nts + seq + last_nts
        if any([
            'AAA' in primer_seq,
            'TTT' in primer_seq,
            'GGG' in primer_seq,
            'CCC' in primer_seq,
            BSPQI_SITE_FWD in primer_seq,
            BSPQI_SITE_REV in primer_seq[:-2]
        ]):
            continue

        primer_tm = primer3_mod.calc_tm(primer_seq)
        primer_hairpin = primer3_mod.calc_hairpin(primer_seq)
        primer_homodimer = primer3_mod.calc_homodimer(primer_seq)
        primer_heterodimer = primer3_mod.calc_heterodimer(primer_seq, primer2)

        if all([
            (primer2_tm - 0.1) < primer_tm < (primer2_tm + 0.1),
            primer_hairpin.tm < 20,
            primer_homodimer.tm < 20,
            primer_heterodimer.tm < 20
        ]):
            pid = f'primer_{sha256(primer_seq.encode("utf-8")).hexdigest()}'
            candidate_records.append(
                SeqRecord(Seq(primer_seq), id=pid, description=''))

    log_fn(f'Found {len(candidate_records)} candidate primers')
    if not candidate_records:
        raise RuntimeError(
            'No acceptable amplification primers were found for oligo-pool '
            'construction. This run should not be used for synthesis ordering.'
        )

    cand_primer_fa = os.path.join(output_dir, f'{prefix}_candidate_primers_temp.fasta')
    SeqIO.write(candidate_records, cand_primer_fa, 'fasta')

    t7_probes_fa = os.path.join(output_dir, f'{prefix}_T7_probes_temp.fasta')
    t7_records = [
        SeqRecord(Seq(T7_PROMOTER + b.seq.upper()), id=b.bait_id, description='')
        for b in sorted_baits
    ]
    SeqIO.write(t7_records, t7_probes_fa, 'fasta')

    header = 'qseqid sseqid length nident sstrand'.split()
    primer_blast_out = os.path.join(output_dir, f'{prefix}_primer_blast_temp.txt')
    outfmt = f'6 {" ".join(header)}'
    cmd = [
        'blastn', '-task', 'blastn-short',
        '-query', cand_primer_fa,
        '-out', primer_blast_out,
        '-outfmt', outfmt,
        '-subject', t7_probes_fa,
        '-word_size', '4',
        '-evalue', '50',
        '-dust', 'no',
        '-soft_masking', 'false',
    ]
    subprocess.run(cmd, check=True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)

    try:
        primer_blast_df = pd.read_csv(primer_blast_out, sep='	', names=header)
        if primer_blast_df.empty:
            raise ValueError('No primer BLAST hits were produced.')
        id_dict = {}
        for p in primer_blast_df['qseqid'].unique():
            id_dict[p] = int(primer_blast_df.loc[
                primer_blast_df['qseqid'] == p, 'nident'].sum())
        final_primer_id = min(id_dict, key=id_dict.get)
        log_fn(f'Selected primer {final_primer_id} with '
               f'{id_dict[final_primer_id]} total identities to probe set')
    except Exception as exc:
        raise RuntimeError(
            'Failed to parse primer-selection BLAST results. This run should '
            'not be used for oligo synthesis ordering.'
        ) from exc

    final_primer_seq = None
    for rec in candidate_records:
        if rec.id == final_primer_id:
            final_primer_seq = str(rec.seq)
            break
    if final_primer_seq is None:
        raise RuntimeError('Selected amplification primer could not be recovered.')

    primer_rc = str(Seq(final_primer_seq).reverse_complement())
    if not primer_rc.startswith('CGAAGAGC'):
        raise RuntimeError(
            "Selected primer does not yield the expected 3' BspQI/LguI motif after reverse-complement appending."
        )

    oligo_path = os.path.join(output_dir, f'{prefix}_oligo_pool.fna')
    oligo_records = []
    for b in sorted_baits:
        full_seq = T7_PROMOTER + b.seq.upper() + primer_rc
        rec = SeqRecord(Seq(full_seq), id=b.bait_id, description='')
        oligo_records.append(rec)
    SeqIO.write(oligo_records, oligo_path, 'fasta')

    primer_fasta_path = os.path.join(output_dir, f'{prefix}_amplification_primers.fna')
    with open(primer_fasta_path, 'w') as f:
        f.write(f'>T7_primer\n{T7_PROMOTER}\n')
        f.write(f'>{final_primer_id}\n{final_primer_seq}\n')


    for tmp in [cand_primer_fa, t7_probes_fa, primer_blast_out]:
        if os.path.exists(tmp):
            os.remove(tmp)

    log_fn(f'O-pool design complete. Oligo pool written to {oligo_path}')
    log_fn(f'T7 promoter: {T7_PROMOTER} ({T7_PROMOTER_LEN} nt)')
    log_fn(f'Amplification primer: {final_primer_seq}')
    log_fn(f"Primer RC (appended to 3'): {primer_rc}")
    log_fn(f"Full oligo structure: 5'-{T7_PROMOTER}[probe]{primer_rc}-3'")

    return oligo_path, bare_probes_path, primer_fasta_path, final_primer_seq


# ============================================================================
# Validation metadata helpers
# ============================================================================

def parse_probe_header_metadata(probe_id: str, seqlen: int) -> Dict[str, object]:
    """Extract positional metadata from probe identifiers when present."""
    m = re.match(r"^(?P<ref>.+)\|[^|]+\|pos(?P<start>\d+)-(?P<end>\d+)$", probe_id)
    if m:
        return {
            "probe_id": probe_id,
            "ref_hint": m.group("ref"),
            "start": int(m.group("start")),
            "end": int(m.group("end")),
            "scheme": "explicit_pos",
            "trusted": True,
        }

    m = re.match(r"^(?P<ref>.+)_pos(?P<start>\d+)-(?P<end>\d+)$", probe_id)
    if m:
        return {
            "probe_id": probe_id,
            "ref_hint": m.group("ref"),
            "start": int(m.group("start")),
            "end": int(m.group("end")),
            "scheme": "suffix_pos",
            "trusted": True,
        }

    m = re.match(r"^(?P<ref>.+)_(?P<start0>\d+)$", probe_id)
    if m:
        start0 = int(m.group("start0"))
        return {
            "probe_id": probe_id,
            "ref_hint": m.group("ref"),
            "start": start0 + 1,
            "end": start0 + seqlen,
            "scheme": "suffix_start0",
            "trusted": False,
            "suffix_start0": start0,
        }

    return {
        "probe_id": probe_id,
        "ref_hint": None,
        "start": None,
        "end": None,
        "scheme": "none",
        "trusted": False,
    }


def _normalize_ref_token(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "", str(name)).lower()


def resolve_ref_hint_to_target(ref_hint: Optional[str], target_ids: List[str]) -> Optional[str]:
    if not ref_hint:
        return target_ids[0] if len(target_ids) == 1 else None
    if ref_hint in target_ids:
        return ref_hint

    norm_hint = _normalize_ref_token(ref_hint)
    norm_targets = {tid: _normalize_ref_token(tid) for tid in target_ids}
    exact = [tid for tid, norm in norm_targets.items() if norm == norm_hint]
    if len(exact) == 1:
        return exact[0]

    contains = [
        tid for tid, norm in norm_targets.items()
        if norm_hint and (norm_hint in norm or norm in norm_hint)
    ]
    if len(contains) == 1:
        return contains[0]

    acc_tokens = re.findall(r"[A-Z]{1,4}[_-]?\d+(?:\.\d+)?", str(ref_hint))
    for token in acc_tokens:
        token_norm = _normalize_ref_token(token)
        matched = [tid for tid, norm in norm_targets.items() if token_norm and token_norm in norm]
        if len(matched) == 1:
            return matched[0]

    return target_ids[0] if len(target_ids) == 1 else None


def parse_circular_id_set(circular_ids_arg: Optional[str], available_ids: List[str]) -> Set[str]:
    """Resolve a comma-delimited set of reference IDs to treat as circular."""
    if not circular_ids_arg:
        return set()
    requested = [tok.strip() for tok in str(circular_ids_arg).split(',') if tok.strip()]
    resolved: Set[str] = set()
    for token in requested:
        match = resolve_ref_hint_to_target(token, available_ids)
        if match is None:
            raise RuntimeError(f"Could not resolve circular reference identifier: {token}")
        resolved.add(match)
    return resolved


def infer_trusted_probe_metadata(
    probe_records: List[SeqRecord],
    target_lengths: Dict[str, int],
) -> Dict[str, Dict[str, object]]:
    trusted: Dict[str, Dict[str, object]] = {}
    suffix_candidates: Dict[str, Dict[str, object]] = {}

    for rec in probe_records:
        meta = parse_probe_header_metadata(rec.id, len(rec.seq))
        if meta["scheme"] in {"explicit_pos", "suffix_pos"}:
            trusted[rec.id] = meta
        elif meta["scheme"] == "suffix_start0":
            suffix_candidates[rec.id] = meta

    if suffix_candidates and len(target_lengths) == 1:
        target_len = next(iter(target_lengths.values()))
        vals = sorted(int(meta["start"]) for meta in suffix_candidates.values() if meta["start"] is not None)
        if vals:
            span = vals[-1] - vals[0]
            if span >= max(50, int(0.5 * target_len)) and span <= int(2.0 * target_len):
                for probe_id, meta in suffix_candidates.items():
                    promoted = dict(meta)
                    promoted["trusted"] = True
                    trusted[probe_id] = promoted

    return trusted


def apply_interval_coverage(cov_array: np.ndarray, start_1based: int, end_1based: int) -> None:
    if cov_array.size == 0 or start_1based is None or end_1based is None:
        return
    span = max(0, int(end_1based) - int(start_1based) + 1)
    if span <= 0:
        return
    L = cov_array.size
    start0 = (int(start_1based) - 1) % L
    if start0 + span <= L:
        cov_array[start0:start0 + span] += 1
    else:
        first = L - start0
        cov_array[start0:] += 1
        cov_array[:span - first] += 1


def hsp_to_match_array(hsp) -> np.ndarray:
    matches = ' '.join(hsp.match.replace('|', '1').replace(' ', '0'))
    match_array = np.fromstring(matches, dtype=int, sep=' ')
    if '-' in hsp.sbjct:
        sbjct_gaps = ' '.join(re.sub(r'[ATGC]', '0', hsp.sbjct.replace('-', '1')))
        gap_arr = np.fromstring(sbjct_gaps, dtype=int, sep=' ')
        gap_pos = np.where(gap_arr == 1)
        match_array = np.delete(match_array, gap_pos)
    return match_array


def hsp_subject_start(hsp) -> int:
    return hsp.sbjct_start if hsp.strand[0] == hsp.strand[1] else hsp.sbjct_end


def pick_primary_hit(valid_hits: List[Dict[str, object]], probe_meta: Optional[Dict[str, object]] = None) -> Optional[Dict[str, object]]:
    if not valid_hits:
        return None

    probe_meta = probe_meta or {}
    expected_target = probe_meta.get('target_id')
    expected_start = probe_meta.get('start')

    def key(hit: Dict[str, object]):
        target_penalty = 0
        if expected_target is not None:
            target_penalty = 0 if hit['target'] == expected_target else 1
        start_penalty = 0
        if expected_start is not None:
            start_penalty = abs(int(hit['subject_start']) - int(expected_start))
        return (
            target_penalty,
            start_penalty,
            -int(hit['identities']),
            -int(hit['align_length']),
            -float(hit['bits']),
            int(hit['subject_start']),
        )

    return min(valid_hits, key=key)


# ============================================================================
# BLAST validation and analysis (adapted from CARPDM)
# ============================================================================

def blast_validation(probe_fasta: str, design_fasta: str, output_dir: str,
                     prefix: str, probe_length: int,
                     log_fn, circular_refs: Optional[Set[str]] = None) -> Dict[str, np.ndarray]:
    """BLAST probes against design targets and write analysis tables/plots.

    Coverage assignment is conservative: each probe contributes at most one
    primary on-target placement. When probe identifiers encode positional
    provenance (the standard FlyForge ``|posstart-end`` format), that metadata
    is trusted for coverage placement so repetitive regions do not inflate
    depth by counting every valid BLAST HSP.
    """
    if shutil.which('blastn') is None:
        raise RuntimeError(
            'blastn was not found in PATH. BLAST+ is required for validation and '
            'coverage analysis.'
        )

    circular_refs = set(circular_refs or set())
    design_records = list(SeqIO.parse(design_fasta, 'fasta'))
    original_lengths = {record.id: len(record.seq) for record in design_records}

    blast_subject = design_fasta
    if circular_refs:
        extension = max(0, probe_length - 1)
        extended_records = []
        for record in design_records:
            seq = str(record.seq)
            if record.id in circular_refs and extension > 0 and len(seq) > 0:
                seq = seq + seq[:extension]
            extended_records.append((record.id, seq))
        blast_subject = os.path.join(output_dir, f'{prefix}_validation_targets_for_blast.fna')
        write_fasta(blast_subject, extended_records)

    blast_xml = os.path.join(output_dir, f'{prefix}_final_blast.xml')
    cmd = [
        'blastn',
        '-query', probe_fasta,
        '-subject', blast_subject,
        '-out', blast_xml,
        '-outfmt', '5',
        '-dust', 'no',
        '-soft_masking', 'false',
    ]
    subprocess.run(cmd, check=True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)

    coverage_dict: Dict[str, np.ndarray] = {}
    gc_dict: Dict[str, float] = {}
    target_seq_dict: Dict[str, str] = {}
    probe_target_sets: Dict[str, Set[str]] = defaultdict(set)
    target_count_dict: Dict[str, int] = {}

    probe_records = list(SeqIO.parse(probe_fasta, 'fasta'))
    target_ids = [record.id for record in design_records]
    target_lengths = dict(original_lengths)
    trusted_meta = infer_trusted_probe_metadata(probe_records, target_lengths)

    for record in design_records:
        name = record.id
        coverage_dict[name] = np.zeros(original_lengths[name], dtype=float)
        gc_dict[name] = SeqUtils.gc_fraction(record.seq)
        target_seq_dict[name] = str(record.seq)
        target_count_dict[name] = 0

    identity_cutoff = int(probe_length * 0.625)
    pair_list = []
    seen_pairs: Set[Tuple[str, str]] = set()

    with open(blast_xml) as infile:
        for record in NCBIXML.parse(infile):
            probe = record.query.split()[0]
            meta = trusted_meta.get(probe)
            if meta is not None:
                meta = dict(meta)
                meta['target_id'] = resolve_ref_hint_to_target(meta.get('ref_hint'), target_ids)

            valid_hits: List[Dict[str, object]] = []
            for alignment in record.alignments:
                target = getattr(alignment, 'hit_id', None)
                if not target:
                    hit_def = getattr(alignment, 'hit_def', alignment.title)
                    target = hit_def.split()[0]
                if target not in coverage_dict:
                    continue
                for hsp in alignment.hsps:
                    if hsp.identities <= identity_cutoff:
                        continue
                    if re.search(r'--+', hsp.sbjct) is not None:
                        continue
                    if re.search(r'--+', hsp.query) is not None:
                        continue
                    valid_hits.append({
                        'target': target,
                        'hsp': hsp,
                        'identities': int(hsp.identities),
                        'align_length': int(hsp.align_length),
                        'bits': float(hsp.bits),
                        'subject_start': int(hsp_subject_start(hsp)),
                    })

            assigned_target = None
            if meta is not None and meta.get('target_id') in coverage_dict and meta.get('start') is not None and meta.get('end') is not None:
                assigned_target = str(meta['target_id'])
                apply_interval_coverage(coverage_dict[assigned_target], int(meta['start']), int(meta['end']))
            else:
                primary = pick_primary_hit(valid_hits, meta)
                if primary is not None:
                    assigned_target = str(primary['target'])
                    match_array = hsp_to_match_array(primary['hsp'])
                    target_len = coverage_dict[assigned_target].size
                    if assigned_target in circular_refs:
                        start1 = int(primary['subject_start'])
                        end1 = start1 + int(match_array.size) - 1
                        apply_interval_coverage(coverage_dict[assigned_target], start1, end1)
                    else:
                        pad_left = int(primary['subject_start']) - 1
                        pad_right = target_len - (pad_left + match_array.size)
                        if pad_right >= 0:
                            padded = np.pad(
                                match_array,
                                (pad_left, pad_right),
                                mode='constant',
                                constant_values=(0, 0),
                            )
                            coverage_dict[assigned_target] += padded

            if assigned_target is not None:
                pair = (assigned_target, probe)
                if pair not in seen_pairs:
                    seen_pairs.add(pair)
                    pair_list.append(pair)
                    target_count_dict[assigned_target] += 1
                probe_target_sets[probe].add(assigned_target)

    pair_df = pd.DataFrame(pair_list, columns=['target', 'probe'])
    pair_prefix = os.path.join(output_dir, prefix)
    pair_df.to_csv(f'{pair_prefix}_target_probe_pairs.csv', index=False)

    target_rows = []
    for target in coverage_dict:
        array = coverage_dict[target]
        target_rows.append({
            'target_name': target,
            'length': array.size,
            'gc': gc_dict.get(target, 0),
            'sequence': target_seq_dict.get(target, ''),
            'probe_count': target_count_dict.get(target, 0),
            'proportion_covered': np.count_nonzero(array) / array.size if array.size > 0 else 0,
            'mean_depth': float(np.mean(array)) if array.size > 0 else 0.0,
            'std_dev': float(np.std(array)) if array.size > 0 else 0.0,
            'min_depth': int(np.min(array)) if array.size > 0 else 0,
            'max_depth': int(np.max(array)) if array.size > 0 else 0,
        })
    target_info = pd.DataFrame(target_rows)
    target_info.to_csv(f'{pair_prefix}_target_info.csv', index=False)

    probe_rows = []
    for record in probe_records:
        seq = str(record.seq)
        probe_rows.append({
            'probe_id': record.id,
            'gc': SeqUtils.gc_fraction(seq),
            'tm': mt.Tm_NN(seq, nn_table=mt.R_DNA_NN1),
            'num_targets': len(probe_target_sets.get(record.id, set())),
            'sequence': seq,
            'contains_lgui_site': int(BSPQI_SITE_FWD in seq.upper() or BSPQI_SITE_REV in seq.upper()),
        })
    probe_info = pd.DataFrame(probe_rows)
    probe_info.to_csv(f'{pair_prefix}_probe_info.csv', index=False)

    plot_dir = os.path.join(output_dir, f'{prefix}_plots')
    os.makedirs(plot_dir, exist_ok=True)
    indiv_dir = os.path.join(plot_dir, 'individual_target_coverages')
    os.makedirs(indiv_dir, exist_ok=True)

    for target in coverage_dict:
        cov_array = coverage_dict[target]
        fig, ax = plt.subplots(figsize=(12, 4))
        sns.lineplot(x=np.arange(len(cov_array)), y=cov_array, ax=ax)
        ax.set_title(f'{target} Coverage')
        ax.set_xlabel('Position (bp)')
        ax.set_ylabel('Coverage Depth')
        safe_name = target.replace('/', '|').replace(' ', '_')
        plt.tight_layout()
        plt.savefig(os.path.join(indiv_dir, f'{safe_name}_coverage.png'), dpi=150)
        plt.savefig(os.path.join(indiv_dir, f'{safe_name}_coverage.svg'))
        plt.close(fig)

    if not probe_info.empty and 'gc' in probe_info.columns:
        fig, ax = plt.subplots()
        sns.violinplot(y=probe_info['gc'].astype(float), ax=ax)
        ax.set_title('Probe GC Content Distribution')
        plt.tight_layout()
        plt.savefig(os.path.join(plot_dir, 'probe_gc.png'), dpi=300)
        plt.savefig(os.path.join(plot_dir, 'probe_gc.svg'))
        plt.close(fig)

    if not probe_info.empty and 'tm' in probe_info.columns:
        fig, ax = plt.subplots()
        sns.violinplot(y=probe_info['tm'].astype(float), ax=ax)
        ax.set_title('Probe Melting Temperature Distribution')
        plt.tight_layout()
        plt.savefig(os.path.join(plot_dir, 'probe_tm.png'), dpi=300)
        plt.savefig(os.path.join(plot_dir, 'probe_tm.svg'))
        plt.close(fig)

    if not probe_info.empty and 'num_targets' in probe_info.columns:
        fig, ax = plt.subplots()
        sns.violinplot(y=probe_info['num_targets'].astype(float), ax=ax)
        ax.set_title('Number of Targets per Probe')
        plt.tight_layout()
        plt.savefig(os.path.join(plot_dir, 'probe_num_targets.png'), dpi=300)
        plt.savefig(os.path.join(plot_dir, 'probe_num_targets.svg'))
        plt.close(fig)

    if not target_info.empty:
        fig, ax = plt.subplots()
        sns.violinplot(y=target_info['length'].astype(float), ax=ax)
        ax.set_title('Target Length Distribution')
        plt.tight_layout()
        plt.savefig(os.path.join(plot_dir, 'target_len.png'), dpi=300)
        plt.savefig(os.path.join(plot_dir, 'target_len.svg'))
        plt.close(fig)

        fig, ax = plt.subplots()
        sns.violinplot(y=target_info['gc'].astype(float), ax=ax)
        ax.set_title('Target GC Content Distribution')
        plt.tight_layout()
        plt.savefig(os.path.join(plot_dir, 'target_gc.png'), dpi=300)
        plt.savefig(os.path.join(plot_dir, 'target_gc.svg'))
        plt.close(fig)

        fig, ax = plt.subplots()
        sns.violinplot(y=target_info['probe_count'].astype(float), ax=ax)
        ax.set_title('Target Probe Count Distribution')
        plt.tight_layout()
        plt.savefig(os.path.join(plot_dir, 'target_probe_count.png'), dpi=300)
        plt.savefig(os.path.join(plot_dir, 'target_probe_count.svg'))
        plt.close(fig)

        fig, ax = plt.subplots()
        sns.violinplot(y=target_info['proportion_covered'].astype(float), ax=ax)
        ax.set_title('Target Coverage Proportion Distribution')
        plt.tight_layout()
        plt.savefig(os.path.join(plot_dir, 'target_coverage_prop.png'), dpi=300)
        plt.savefig(os.path.join(plot_dir, 'target_coverage_prop.svg'))
        plt.close(fig)

        fig, ax = plt.subplots()
        sns.violinplot(y=target_info['mean_depth'].astype(float), ax=ax)
        ax.set_title('Target Coverage Depth Distribution')
        plt.tight_layout()
        plt.savefig(os.path.join(plot_dir, 'target_coverage_depth.png'), dpi=300)
        plt.savefig(os.path.join(plot_dir, 'target_coverage_depth.svg'))
        plt.close(fig)

        fig, ax = plt.subplots()
        sns.violinplot(y=target_info['std_dev'].astype(float), ax=ax)
        ax.set_title('Target Coverage Std Deviation Distribution')
        plt.tight_layout()
        plt.savefig(os.path.join(plot_dir, 'target_coverage_stdev.png'), dpi=300)
        plt.savefig(os.path.join(plot_dir, 'target_coverage_stdev.svg'))
        plt.close(fig)

    log_fn(f'BLAST validation and analysis complete. Plots in {plot_dir}')
    return coverage_dict


# ============================================================================
# Summary statistics
# ============================================================================

def summarize_baits(step_name: str, baits: List[Bait]) -> Dict[str, float]:
    if not baits:
        return {
            "step": step_name, "n_baits": 0,
            "gc_mean": 0.0, "gc_min": 0.0, "gc_max": 0.0,
            "tm_mean": 0.0, "tm_min": 0.0, "tm_max": 0.0,
            "masked_mean": 0.0, "ambiguous_mean": 0.0,
        }
    gcs = [b.gc_frac for b in baits]
    tms = [b.tm for b in baits]
    masks = [b.masked_frac for b in baits]
    ambs = [b.ambiguous_count for b in baits]
    return {
        "step": step_name,
        "n_baits": len(baits),
        "gc_mean": sum(gcs) / len(baits),
        "gc_min": min(gcs),
        "gc_max": max(gcs),
        "tm_mean": sum(tms) / len(baits),
        "tm_min": min(tms),
        "tm_max": max(tms),
        "masked_mean": sum(masks) / len(baits),
        "ambiguous_mean": sum(ambs) / len(baits),
    }


# ============================================================================
# Main pipeline
# ============================================================================

def run_pipeline(args):
    # Setup output directory
    os.makedirs(args.output_dir, exist_ok=True)

    log, log_handle = make_logger(args.progress_log)
    t_start = time.time()

    # Track steps for progress
    steps_planned = ["preprocessing"]
    if args.remove_complements:
        steps_planned.append("complement_removal")
    if not args.skip_self_mask:
        steps_planned.append("self_masking")
    steps_planned.append("tiling")
    steps_planned.append("ambiguous_filter")
    steps_planned.append("masked_filter")
    if args.min_tm > 0:
        steps_planned.append("tm_filter")
    if not args.no_opool:
        steps_planned.append("lgui_filter")
    steps_planned.append("complement_bait_filter")
    if args.blast_db:
        steps_planned.append("blast_filter")
    if args.specificity_db:
        steps_planned.append("specificity_filter")
    if not args.no_redundancy:
        steps_planned.append("self_blast_redundancy")
    if not args.no_cluster:
        steps_planned.append("clustering")
    if not args.no_opool:
        steps_planned.append("opool_design")
    steps_planned.append("validation")
    steps_planned.append("write_output")

    progress = ProgressTracker(steps_planned, log, t_start)

    log(f"TackleBox: FlyForge v{__version__}")
    log(f"Starting pipeline with {len(args.input)} input FASTA file(s).")

    # ===== 1. Read and preprocess sequences =====
    progress.start_step("preprocessing")
    all_seqs: Dict[str, str] = {}
    ref_stats: Dict[str, RefStats] = {}

    for path in args.input:
        seqs = read_fasta(path)
        for sid, seq in seqs.items():
            processed, n_mod = preprocess_sequence(seq)
            all_seqs[sid] = processed
            ref_stats[sid] = RefStats(
                ref_id=sid,
                length_original=len(seq),
                length_preprocessed=len(processed),
                modified_bases=n_mod,
                frac_modified=(n_mod / len(seq)) if len(seq) > 0 else 0.0,
            )

    progress.finish_step("preprocessing",
                         details=f"{len(all_seqs)} reference sequences")

    lengths_dict = {sid: len(seq) for sid, seq in all_seqs.items()}
    total_bp = sum(lengths_dict.values())
    log(f"Total reference length: {total_bp:,} bp")

    if args.max_total_bp and total_bp > args.max_total_bp:
        log(f"ERROR: total reference length {total_bp} bp exceeds "
            f"max_total_bp={args.max_total_bp}. Aborting.")
        if log_handle is not None:
            log_handle.close()
        sys.exit(1)

    # ===== 2. Optional complementary region removal =====
    if args.remove_complements:
        progress.start_step("complement_removal")
        input_records = []
        for sid, seq in all_seqs.items():
            input_records.append(SeqRecord(Seq(seq), id=sid, description=''))
        comp_out = remove_complementary_targets(
            input_records, args.output_dir,
            os.path.join(args.output_dir, args.prefix),
            args.threads, args.bait_length, log)
        # Re-read processed sequences
        all_seqs = read_fasta(comp_out)
        lengths_dict = {sid: len(seq) for sid, seq in all_seqs.items()}
        progress.finish_step("complement_removal",
                             details=f"{len(all_seqs)} segments after complement removal")

    # ===== 3. Determine tiling density =====
    used_tiling_density = args.tiling_density
    estimated_baits_for_used = None
    if args.max_baits is not None:
        used_tiling_density, est_baits = choose_tiling_density(
            lengths_dict, args.bait_length, args.omit_short_leq,
            args.pad_min, args.max_baits, args.min_tiling_density,
            args.tiling_density)
        estimated_baits_for_used = est_baits
        if est_baits > args.max_baits:
            log(f"WARNING: even at min tiling density "
                f"{used_tiling_density:.3f}, estimated bait count "
                f"{est_baits} exceeds max_baits={args.max_baits}. "
                "Using min density.")
        else:
            log(f"Auto-adjusted tiling density from "
                f"{args.tiling_density} to {used_tiling_density:.3f} "
                f"to target max_baits={args.max_baits} "
                f"(estimated ~{est_baits}).")
    else:
        log(f"Using tiling density {used_tiling_density}.")

    circular_ids = set(lengths_dict.keys()) if args.circular else parse_circular_id_set(args.circular_ids, list(lengths_dict.keys()))
    if circular_ids:
        log("Circular tiling enabled for: " + ", ".join(sorted(circular_ids)))

    # ===== 4. Self-repeat softmasking =====
    if args.skip_self_mask:
        masked_seqs = all_seqs
    else:
        progress.start_step("self_masking")
        masked_seqs = self_repeat_softmask(
            all_seqs, k=args.repeat_k, threshold=args.repeat_threshold)
        progress.finish_step("self_masking")

    # ===== 5. Tiling =====
    progress.start_step("tiling")
    all_baits: List[Bait] = []
    ref_baits_after_tiling: Dict[str, List[Bait]] = {}

    for ref_id, seq in masked_seqs.items():
        baits = tile_sequence(
            ref_id, seq, bait_len=args.bait_length,
            tiling_density=used_tiling_density,
            omit_short_leq=args.omit_short_leq, pad_min=args.pad_min,
            circular=(ref_id in circular_ids))
        all_baits.extend(baits)
        ref_baits_after_tiling[ref_id] = baits
        rs = ref_stats.get(ref_id)
        if rs:
            rs.n_baits_tiled = len(baits)
            if len(seq) <= args.omit_short_leq:
                rs.dropped_for_length = True

    # Per-reference coverage from tiled baits
    if not args.no_coverage:
        for ref_id, seq in masked_seqs.items():
            rs = ref_stats.get(ref_id)
            if rs is None:
                continue
            b_list = ref_baits_after_tiling.get(ref_id, [])
            L = len(seq)
            if L == 0 or not b_list:
                continue
            cov = [0] * L
            for b in b_list:
                start = max(0, b.ref_start - 1)
                end = min(b.ref_end, L)
                for i in range(start, end):
                    cov[i] += 1
            rs.coverage_mean = sum(cov) / L
            rs.coverage_min = min(cov)
            rs.coverage_max = max(cov)
            rs.coverage_fraction_covered = sum(1 for c in cov if c > 0) / L

    stats: List[Dict] = []
    stats.append(summarize_baits("after_tiling", all_baits))
    progress.finish_step("tiling", details=f"{len(all_baits)} baits after tiling")

    # ===== 6. Ambiguous-base filter =====
    progress.start_step("ambiguous_filter")
    all_baits, removed_amb = filter_ambiguous(all_baits, args.ambiguous_cutoff)
    s = summarize_baits("after_ambiguous_filter", all_baits)
    s["removed"] = removed_amb
    stats.append(s)
    _update_ref_counts(ref_stats, all_baits, 'n_baits_after_amb')
    progress.finish_step("ambiguous_filter",
                         details=f"{removed_amb} removed; {len(all_baits)} remain")

    # ===== 7. Masked-fraction filter =====
    progress.start_step("masked_filter")
    all_baits, removed_mask = filter_masked_fraction(all_baits, args.max_masked_frac)
    s = summarize_baits("after_masked_filter", all_baits)
    s["removed"] = removed_mask
    stats.append(s)
    _update_ref_counts(ref_stats, all_baits, 'n_baits_after_mask')
    progress.finish_step("masked_filter",
                         details=f"{removed_mask} removed; {len(all_baits)} remain")

    # ===== 8. Melting temperature filter =====
    if args.min_tm > 0:
        progress.start_step("tm_filter")
        all_baits, removed_tm = filter_melting_temp(all_baits, args.min_tm)
        s = summarize_baits("after_tm_filter", all_baits)
        s["removed"] = removed_tm
        stats.append(s)
        _update_ref_counts(ref_stats, all_baits, 'n_baits_after_tm')
        progress.finish_step("tm_filter",
                             details=f"{removed_tm} removed (Tm < {args.min_tm}); "
                                     f"{len(all_baits)} remain")

    # ===== 9. LguI site filter (if o-pool mode) =====
    if not args.no_opool:
        progress.start_step("lgui_filter")
        all_baits, removed_lgui = filter_lgui_sites(all_baits)
        s = summarize_baits("after_lgui_filter", all_baits)
        s["removed"] = removed_lgui
        stats.append(s)
        progress.finish_step("lgui_filter",
                             details=f"{removed_lgui} removed (LguI/BspQI sites); "
                                     f"{len(all_baits)} remain")

    # ===== 10. Complementary bait filter =====
    progress.start_step("complement_bait_filter")
    all_baits, removed_comp = filter_complementary_baits(all_baits)
    s = summarize_baits("after_complement_filter", all_baits)
    s["removed"] = removed_comp
    stats.append(s)
    progress.finish_step("complement_bait_filter",
                         details=f"{removed_comp} removed (perfect complements); "
                                 f"{len(all_baits)} remain")

    # ===== 11. Optional BLAST filter (pident-based) =====
    if args.blast_db:
        progress.start_step("blast_filter")
        all_baits, removed_blast = blast_filter_baits(
            all_baits, blast_db=args.blast_db, evalue=args.blast_evalue,
            min_pident=args.blast_min_pident,
            max_hits_per_query=args.blast_max_hits, threads=args.threads)
        s = summarize_baits("after_blast_filter", all_baits)
        s["removed"] = removed_blast
        stats.append(s)
        _update_ref_counts(ref_stats, all_baits, 'n_baits_after_blast')
        progress.finish_step("blast_filter",
                             details=f"{removed_blast} removed; {len(all_baits)} remain")

    # ===== 12. Optional specificity filter (nident-based, from CARPDM) =====
    if args.specificity_db:
        progress.start_step("specificity_filter")
        is_nt = (os.path.basename(args.specificity_db) == 'nt')
        all_baits, removed_spec = blast_specificity_filter(
            all_baits, blast_db=args.specificity_db,
            probe_length=args.bait_length, threads=args.threads,
            is_nt=is_nt)
        s = summarize_baits("after_specificity_filter", all_baits)
        s["removed"] = removed_spec
        stats.append(s)
        progress.finish_step("specificity_filter",
                             details=f"{removed_spec} removed; {len(all_baits)} remain")

    # ===== 13. Self-BLAST redundancy filter =====
    if not args.no_redundancy:
        progress.start_step("self_blast_redundancy")
        all_baits, n_comp_rm, n_red_rm = self_blast_filter(
            all_baits, args.output_dir, args.prefix, args.threads,
            args.probe_num_cutoff, log)
        s = summarize_baits("after_self_blast_filter", all_baits)
        s["removed_complementary"] = n_comp_rm
        s["removed_redundant"] = n_red_rm
        stats.append(s)
        _update_ref_counts(ref_stats, all_baits, 'n_baits_after_redundancy')
        progress.finish_step("self_blast_redundancy",
                             details=f"{n_comp_rm + n_red_rm} removed; "
                                     f"{len(all_baits)} remain")

    # ===== 14. cd-hit-est clustering =====
    if not args.no_cluster:
        progress.start_step("clustering")
        all_baits, removed_cluster = cluster_baits_cd_hit(
            all_baits, identity=args.cluster_identity,
            overlap=args.cluster_overlap, threads=args.threads)
        s = summarize_baits("after_clustering", all_baits)
        s["removed"] = removed_cluster
        stats.append(s)
        _update_ref_counts(ref_stats, all_baits, 'n_baits_after_cluster')
        progress.finish_step("clustering",
                             details=f"{removed_cluster} removed; {len(all_baits)} remain")

    # Update final counts
    _update_ref_counts(ref_stats, all_baits, 'n_baits_final')
    total_final = len(all_baits)
    for rs in ref_stats.values():
        rs.proportion_of_final_baits = (
            rs.n_baits_final / total_final if total_final > 0 else 0)

    # ===== 15. O-pool design =====
    oligo_path = None
    bare_probes_path = None
    primer_fasta_path = None
    primer_seq = ''
    if not args.no_opool:
        progress.start_step("opool_design")
        oligo_path, bare_probes_path, primer_fasta_path, primer_seq = design_opool(
            all_baits, args.output_dir, args.prefix, log)
        progress.finish_step("opool_design",
                             details=f"Oligo pool: {oligo_path}")
    else:
        # Write bare probes even without o-pool
        bare_probes_path = os.path.join(args.output_dir,
                                        f'{args.prefix}_probes.fna')
        write_fasta(bare_probes_path,
                    [(b.bait_id, b.seq.upper())
                     for b in sorted(all_baits, key=lambda x: x.bait_id)])

    # ===== 16. BLAST validation =====
    progress.start_step("validation")
    input_fasta_path = args.input[0]  # Use first input for validation
    blast_validation(bare_probes_path, input_fasta_path,
                     args.output_dir, args.prefix,
                     args.bait_length, log,
                     circular_refs=circular_ids)

    target_csv = os.path.join(args.output_dir, f'{args.prefix}_target_info.csv')
    try:
        target_df = pd.read_csv(target_csv)
        target_df = target_df.set_index('target_name')
        for ref_id, rs in ref_stats.items():
            if ref_id in target_df.index:
                row = target_df.loc[ref_id]
                rs.coverage_mean = float(row.get('mean_depth', 0.0))
                rs.coverage_min = int(row.get('min_depth', 0))
                rs.coverage_max = int(row.get('max_depth', 0))
                rs.coverage_fraction_covered = float(row.get('proportion_covered', 0.0))
    except Exception as exc:
        log(f'WARNING: unable to refresh per-reference final coverage from validation output: {exc}')

    progress.finish_step("validation", details="BLAST validation and plotting complete")

    # ===== 17. Write final outputs =====
    progress.start_step("write_output")
    # Write main bait FASTA (bare probes)
    final_baits_path = os.path.join(args.output_dir,
                                     f'{args.prefix}_final_baits.fa')
    write_fasta(final_baits_path,
                [(b.bait_id, b.seq.upper())
                 for b in sorted(all_baits, key=lambda x: x.bait_id)])

    # Write per-reference summary
    per_ref_path = os.path.join(args.output_dir,
                                f'{args.prefix}_per_ref_stats.tsv')
    with open(per_ref_path, "w") as out:
        header = [
            "ref_id", "length_original", "length_preprocessed",
            "modified_bases", "frac_modified", "dropped_for_length",
            "n_baits_tiled", "n_baits_after_amb", "n_baits_after_mask",
            "n_baits_after_tm", "n_baits_after_blast",
            "n_baits_after_cluster", "n_baits_after_redundancy",
            "n_baits_final",
            "coverage_mean", "coverage_min", "coverage_max",
            "coverage_fraction_covered", "proportion_of_final_baits",
        ]
        out.write("\t".join(header) + "\n")
        for rs in ref_stats.values():
            row = [
                rs.ref_id, str(rs.length_original), str(rs.length_preprocessed),
                str(rs.modified_bases), f"{rs.frac_modified:.6f}",
                str(int(rs.dropped_for_length)),
                str(rs.n_baits_tiled), str(rs.n_baits_after_amb),
                str(rs.n_baits_after_mask), str(rs.n_baits_after_tm),
                str(rs.n_baits_after_blast), str(rs.n_baits_after_cluster),
                str(rs.n_baits_after_redundancy), str(rs.n_baits_final),
                f"{rs.coverage_mean:.6f}", str(rs.coverage_min),
                str(rs.coverage_max), f"{rs.coverage_fraction_covered:.6f}",
                f"{rs.proportion_of_final_baits:.6f}",
            ]
            out.write("\t".join(row) + "\n")

    # Write comprehensive summary
    summary_path = os.path.join(args.output_dir,
                                f'{args.prefix}_summary.tsv')
    with open(summary_path, "w") as out:
        out.write("# TackleBox: FlyForge v{}\n".format(__version__))
        out.write(f"# Date: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        out.write("# ===== Parameters =====\n")
        param_info = {
            "input_files": ", ".join(args.input),
            "output_dir": args.output_dir,
            "prefix": args.prefix,
            "bait_length": args.bait_length,
            "tiling_density_requested": args.tiling_density,
            "tiling_density_used": f"{used_tiling_density:.4f}",
            "estimated_pre_filter_baits": estimated_baits_for_used or "N/A",
            "max_baits": args.max_baits or "N/A",
            "min_tiling_density": args.min_tiling_density,
            "omit_short_leq": args.omit_short_leq,
            "pad_min": args.pad_min,
            "ambiguous_cutoff": args.ambiguous_cutoff,
            "max_masked_frac": args.max_masked_frac,
            "min_tm": args.min_tm,
            "circular": args.circular,
            "circular_ids": ",".join(sorted(circular_ids)) if circular_ids else "N/A",
            "repeat_k": args.repeat_k,
            "repeat_threshold": args.repeat_threshold,
            "skip_self_mask": args.skip_self_mask,
            "remove_complements": args.remove_complements,
            "blast_db": args.blast_db or "N/A",
            "blast_evalue": args.blast_evalue,
            "blast_min_pident": args.blast_min_pident,
            "specificity_db": args.specificity_db or "N/A",
            "cluster_identity": args.cluster_identity,
            "cluster_overlap": args.cluster_overlap,
            "no_cluster": args.no_cluster,
            "no_redundancy": args.no_redundancy,
            "probe_num_cutoff": args.probe_num_cutoff,
            "no_opool": args.no_opool,
            "threads": args.threads,
            "n_input_sequences": len(ref_stats),
            "total_input_bp": total_bp,
            "final_bait_count": len(all_baits),
        }
        if primer_seq:
            param_info["T7_primer"] = T7_PROMOTER
            param_info["amplification_primer"] = primer_seq
            param_info["oligo_length"] = (
                len(T7_PROMOTER) + args.bait_length +
                len(str(Seq(primer_seq).reverse_complement())))

        for k, v in param_info.items():
            out.write(f"# {k}\t{v}\n")

        out.write("# ===== Step-wise Statistics =====\n")
        header_written = False
        for s in stats:
            if not header_written:
                out.write("\t".join(s.keys()) + "\n")
                header_written = True
            out.write("\t".join(str(v) for v in s.values()) + "\n")

    progress.finish_step("write_output", details=f"Final baits: {len(all_baits)}")

    total_time = time.time() - t_start

    # ===== Comprehensive end-of-run summary =====
    # Compute final probe-level stats
    final_gcs = [b.gc_frac for b in all_baits] if all_baits else [0]
    final_tms = [b.tm for b in all_baits] if all_baits else [0]

    # Read target info for coverage summary
    target_csv = os.path.join(args.output_dir, f'{args.prefix}_target_info.csv')
    try:
        target_df = pd.read_csv(target_csv)
        mean_cov_prop = target_df['proportion_covered'].mean()
        mean_cov_depth = target_df['mean_depth'].mean()
        n_targets = len(target_df)
    except Exception:
        mean_cov_prop = 0
        mean_cov_depth = 0
        n_targets = len(ref_stats)

    # Gather step-wise bait counts for the filter summary
    filter_summary_lines = []
    for s in stats:
        name = s.get("step", "")
        n = s.get("n_baits", 0)
        removed = s.get("removed", s.get("removed_complementary", 0) + s.get("removed_redundant", 0))
        if removed:
            filter_summary_lines.append(f"    {name:<35s} {n:>8,d} baits  ({removed:,d} removed)")
        else:
            filter_summary_lines.append(f"    {name:<35s} {n:>8,d} baits")

    # Build summary text
    sep = "=" * 72
    lines = [
        "",
        sep,
        "  TackleBox: FlyForge v{} — Run Summary".format(__version__),
        sep,
        "",
        "  INPUT",
        f"    Input file(s):           {', '.join(os.path.basename(f) for f in args.input)}",
        f"    Reference sequences:     {n_targets:,d}",
        f"    Total reference bp:      {total_bp:,d}",
        "",
        "  DESIGN PARAMETERS",
        f"    Bait length:             {args.bait_length} bp",
        f"    Tiling density:          {used_tiling_density:.2f}x  (step size: {max(1, int(round(args.bait_length / used_tiling_density)))} bp)",
        f"    Min melting temp:        {args.min_tm} C" if args.min_tm > 0 else f"    Min melting temp:        off",
        f"    Circular tiling:         {'yes' if circular_ids else 'no'}" + (f" ({', '.join(sorted(circular_ids))})" if circular_ids else ""),
        f"    Complement removal:      {'yes' if args.remove_complements else 'no'}",
        f"    Self-repeat masking:     {'yes (k={}, threshold={})'.format(args.repeat_k, args.repeat_threshold) if not args.skip_self_mask else 'no (skipped)'}",
        "",
        "  FILTERING PIPELINE",
    ] + filter_summary_lines + [
        "",
        "  FINAL BAIT SET",
        f"    Total probes:            {len(all_baits):,d}",
        f"    GC content:              mean {sum(final_gcs)/len(final_gcs):.1%}, range {min(final_gcs):.1%}–{max(final_gcs):.1%}",
        f"    Melting temp:            mean {sum(final_tms)/len(final_tms):.1f} C, range {min(final_tms):.1f}–{max(final_tms):.1f} C",
        f"    Target coverage:         {mean_cov_prop:.1%} of reference bases covered",
        f"    Mean coverage depth:     {mean_cov_depth:.2f}x",
    ]

    if not args.no_opool and primer_seq:
        oligo_len = len(T7_PROMOTER) + args.bait_length + len(str(Seq(primer_seq).reverse_complement()))
        lines += [
            "",
            "  OLIGO POOL SYNTHESIS",
            f"    T7 promoter:             {T7_PROMOTER} ({T7_PROMOTER_LEN} nt)",
            f"    Amplification primer:    {primer_seq} ({len(primer_seq)} nt)",
            f"    Primer RC (3' append):   {str(Seq(primer_seq).reverse_complement())}",
            f"    Total oligo length:      {oligo_len} nt",
            f"    Structure:               5'-T7-[{args.bait_length}bp probe]-RC_primer-3'",
        ]

    if args.blast_db:
        lines.append(f"    BLAST filter DB:         {args.blast_db}")
    if args.specificity_db:
        lines.append(f"    Specificity filter DB:   {args.specificity_db}")

    lines += [
        "",
        "  OUTPUT FILES",
        f"    {args.prefix}_final_baits.fa           Bare probe sequences (for QC / mapping)",
    ]
    if oligo_path:
        lines.append(f"    {os.path.basename(oligo_path):<40s} Full oligo sequences (upload to Twist)")
    if primer_fasta_path:
        lines.append(f"    {os.path.basename(primer_fasta_path):<40s} Amplification primers")
    lines += [
        f"    {args.prefix}_summary.tsv              Parameters and step-wise statistics",
        f"    {args.prefix}_per_ref_stats.tsv         Per-reference coverage and bait counts",
        f"    {args.prefix}_target_info.csv           Per-target analysis",
        f"    {args.prefix}_probe_info.csv            Per-probe analysis",
        f"    {args.prefix}_target_probe_pairs.csv    Target:probe pairings",
        f"    {args.prefix}_plots/                    Violin plots and coverage plots",
        "",
        f"  Total runtime: {format_duration(total_time)}",
        sep,
        "",
    ]

    summary_text = "\n".join(lines)
    print(summary_text, file=sys.stderr)
    if log_handle is not None:
        print(summary_text, file=log_handle)
        log_handle.flush()

    log(f"Pipeline complete. All outputs in: {args.output_dir}/")

    if log_handle is not None:
        log_handle.close()


def _update_ref_counts(ref_stats: Dict[str, RefStats], baits: List[Bait],
                       attr: str) -> None:
    """Update per-reference bait counts for a given attribute."""
    counts = defaultdict(int)
    for b in baits:
        counts[b.ref_id] += 1
    for ref_id, rs in ref_stats.items():
        setattr(rs, attr, counts.get(ref_id, 0))


# ============================================================================
# Argument parsing
# ============================================================================

def main():
    ap = argparse.ArgumentParser(
        description=(
            "TackleBox: FlyForge - Comprehensive bait/probe designer for "
            "hybridization capture enrichment with in-house synthesis support."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic bait design with 4x tiling
  FlyForge -i target.fasta --prefix my_baits --tiling-density 4

  # With BLAST specificity filter and o-pool synthesis
  FlyForge -i target.fasta --prefix my_baits \\
    --specificity-db /path/to/nt --threads 8

  # Match MyBaits Expert panel (80nt, ~4x tiling)
  FlyForge -i NC_003427.fasta --prefix bear_mito \\
    --bait-length 80 --tiling-density 4 --threads 4
        """
    )

    # Core I/O
    ap.add_argument("-i", "--input", nargs="+", required=True,
                    help="Input FASTA file(s).")
    ap.add_argument("--prefix", required=True,
                    help="Prefix for all output files.")
    ap.add_argument("--output-dir", default="flyforge_output",
                    help="Output directory (default: flyforge_output).")
    ap.add_argument("--progress-log", default=None,
                    help="Progress log file path.")

    # Design parameters
    ap.add_argument("--bait-length", type=int, default=80,
                    help="Bait/probe length (default: 80).")
    ap.add_argument("--tiling-density", type=float, default=3.0,
                    help="Tiling density; step = bait_length/density (default: 3.0).")
    ap.add_argument("--omit-short-leq", type=int, default=70,
                    help="Drop refs with length <= this (default: 70).")
    ap.add_argument("--pad-min", type=int, default=71,
                    help="Pad refs >= this and < bait_length to full bait (default: 71).")
    ap.add_argument("--min-tm", type=float, default=0.0,
                    help="Minimum melting temperature filter (default: 0 = off). "
                         "Recommended: 50 for stringent filtering.")
    ap.add_argument("--circular", action="store_true",
                    help="Treat all input references as circular when tiling and validation.")
    ap.add_argument("--circular-ids", default=None,
                    help="Comma-delimited subset of reference IDs to treat as circular.")

    # Filtering
    ap.add_argument("--ambiguous-cutoff", type=int, default=10,
                    help="Remove baits with >= this many ambiguous bases (default: 10).")
    ap.add_argument("--max-masked-frac", type=float, default=0.25,
                    help="Max lowercase (repeat-masked) fraction (default: 0.25).")

    # Self-repeat masking
    ap.add_argument("--repeat-k", type=int, default=15,
                    help="k-mer size for self-repeat masking (default: 15).")
    ap.add_argument("--repeat-threshold", type=int, default=3,
                    help="Min k-mer count to flag as repetitive (default: 3).")
    ap.add_argument("--skip-self-mask", action="store_true",
                    help="Skip internal self-repeat masking.")

    # Complementary region removal
    ap.add_argument("--remove-complements", action="store_true",
                    help="Remove complementary regions from input sequences "
                         "(CARPDM-style, recommended for multi-gene inputs).")

    # BLAST filter (pident-based)
    ap.add_argument("--blast-db",
                    help="BLASTn DB for percent-identity filtering.")
    ap.add_argument("--blast-evalue", type=float, default=1e-5,
                    help="BLAST e-value cutoff (default: 1e-5).")
    ap.add_argument("--blast-min-pident", type=float, default=90.0,
                    help="BLAST min percent identity to remove bait (default: 90).")
    ap.add_argument("--blast-max-hits", type=int, default=1,
                    help="BLAST max_target_seqs per bait (default: 1).")

    # Specificity filter (nident-based, from CARPDM)
    ap.add_argument("--specificity-db",
                    help="BLAST DB for nident-based specificity filtering. "
                         "If basename is 'nt', uses taxonomy-aware filtering.")

    # Clustering
    ap.add_argument("--no-cluster", action="store_true",
                    help="Skip cd-hit-est clustering.")
    ap.add_argument("--cluster-identity", type=float, default=0.95,
                    help="cd-hit-est identity threshold (default: 0.95).")
    ap.add_argument("--cluster-overlap", type=float, default=0.83,
                    help="cd-hit-est overlap threshold (default: 0.83).")

    # Self-BLAST redundancy
    ap.add_argument("--no-redundancy", action="store_true",
                    help="Skip self-BLAST redundancy/complementarity filter.")
    ap.add_argument("--probe-num-cutoff", type=int, default=100000,
                    help="Max probe count target for redundancy filter (default: 100000).")

    # O-pool synthesis
    ap.add_argument("--no-opool", action="store_true",
                    help="Skip oligo pool design (T7 + primer). "
                         "Only use if ordering RNA probes directly.")

    # Bait count control
    ap.add_argument("--max-baits", type=int, default=None,
                    help="Target max bait count; auto-adjusts tiling density.")
    ap.add_argument("--min-tiling-density", type=float, default=1.0,
                    help="Min tiling density when auto-adjusting (default: 1.0).")
    ap.add_argument("--max-total-bp", type=int, default=0,
                    help="Max total reference bp; aborts if exceeded (default: 0 = off).")
    ap.add_argument("--no-coverage", action="store_true",
                    help="Skip per-base coverage calculation.")

    # Performance
    ap.add_argument("--threads", type=int, default=1,
                    help="Threads for BLAST/cd-hit-est (default: 1).")

    # Version
    ap.add_argument("-v", "--version", action="version",
                    version=f"TackleBox: FlyForge v{__version__}")

    args = ap.parse_args()

    # Set up progress log default
    if args.progress_log is None:
        args.progress_log = os.path.join(args.output_dir,
                                         f'{args.prefix}_progress.log')

    run_pipeline(args)


if __name__ == "__main__":
    main()
