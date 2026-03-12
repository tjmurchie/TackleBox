#!/usr/bin/env python3
"""
TackleBox: FlyForgeAudit v1.0.0

Companion post-design auditing and panel-augmentation module for FlyForge.

Modes
-----
audit
    Evaluate an existing bait set against one or more reference FASTA files,
    reproduce FlyForge-style coverage statistics/plots, and optionally screen
    the bait panel against an avoidance BLAST database.

augment
    Start from an existing bait set plus new target references, identify regions
    of the new targets that are not adequately covered by the current panel,
    design the minimal additional bait set needed to meet the requested coverage
    threshold, screen those new baits against the existing panel, and emit both
    the extra bare-bait FASTA and an order-ready oligo-pool FASTA.

Author: Tyler J. Murchie
Created with assistance from Claude (Anthropic) and ChatGPT 5.4 (OpenAI).
License: MIT
"""

import argparse
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from Bio import SeqIO, SeqUtils
from Bio.Blast import NCBIXML
from Bio.Seq import Seq

import FlyForge as ff

__version__ = "1.1.0"


# ============================================================================
# Small utilities
# ============================================================================


def log_factory(log_path: str):
    handle = open(log_path, "w")

    def log(msg: str):
        ts = time.strftime("%Y-%m-%d %H:%M:%S")
        line = f"[{ts}] {msg}"
        print(line, file=sys.stderr)
        print(line, file=handle)
        handle.flush()

    return log, handle


@dataclass
class AnalysisResult:
    coverage_dict: Dict[str, np.ndarray]
    target_info: pd.DataFrame
    probe_info: pd.DataFrame
    pair_df: pd.DataFrame
    blast_xml: str


@dataclass
class PreparedTargets:
    raw_seqs: Dict[str, str]
    processed_seqs: Dict[str, str]
    masked_seqs: Dict[str, str]
    fasta_path: str
    blast_fasta_path: str
    circular_ids: set
    modified_bases: Dict[str, int]


def summarize_recommendations(recommendations: pd.DataFrame) -> Tuple[int, int]:
    if recommendations is None or recommendations.empty:
        return 0, 0
    severities = recommendations.get("severity", pd.Series(dtype=str)).astype(str).str.lower()
    informational = int((severities == "info").sum())
    actionable = int(len(recommendations) - informational)
    return actionable, informational


def format_recommendation_lines(recommendations: pd.DataFrame) -> List[str]:
    if recommendations is None or recommendations.empty:
        return ["    - No recommendations generated."]
    lines = []
    for row in recommendations.itertuples(index=False):
        sev = str(getattr(row, "severity", "info")).upper()
        category = getattr(row, "category", "note")
        target = getattr(row, "target", "panel")
        rec = getattr(row, "recommendation", "")
        evidence = getattr(row, "evidence", "")
        lines.append(f"    - [{sev}] {category} | {target}: {rec}")
        if evidence:
            lines.append(f"      Evidence: {evidence}")
    return lines


def build_extended_target_fasta(raw_seqs: Dict[str, str], bait_length: int, circular_ids: set, output_path: str) -> str:
    records = []
    extension = max(0, bait_length - 1)
    for ref_id, seq in raw_seqs.items():
        if ref_id in circular_ids and extension > 0 and len(seq) > 0:
            seq = seq + seq[:extension]
        records.append((ref_id, seq))
    ff.write_fasta(output_path, records)
    return output_path


def ensure_in_path(programs: Iterable[str]) -> None:
    missing = [p for p in programs if shutil.which(p) is None]
    if missing:
        raise RuntimeError(
            "Required external program(s) not found in PATH: " + ", ".join(missing)
        )


def read_fasta_paths(paths: List[str]) -> Dict[str, str]:
    seqs: Dict[str, str] = {}
    for path in paths:
        for rec in SeqIO.parse(path, "fasta"):
            if rec.id in seqs:
                raise RuntimeError(
                    f"Duplicate FASTA identifier detected across inputs: {rec.id}"
                )
            seqs[rec.id] = str(rec.seq)
    return seqs


def parse_bait_id_metadata(bait_id: str, seqlen: int) -> Tuple[str, int, int]:
    meta = ff.parse_probe_header_metadata(bait_id, seqlen)
    ref_hint = str(meta.get("ref_hint") or bait_id)
    start = int(meta.get("start") or 1)
    end = int(meta.get("end") or seqlen)
    return ref_hint, start, end


def read_baits_as_objects(fasta_path: str) -> List[ff.Bait]:
    baits: List[ff.Bait] = []
    for rec in SeqIO.parse(fasta_path, "fasta"):
        seq = str(rec.seq).upper()
        ref_id, start, end = parse_bait_id_metadata(rec.id, len(seq))
        gc_frac, masked_frac, ambiguous = ff.compute_bait_metrics(seq)
        baits.append(
            ff.Bait(
                bait_id=rec.id,
                seq=seq,
                ref_id=ref_id,
                ref_start=start,
                ref_end=end,
                gc_frac=gc_frac,
                masked_frac=masked_frac,
                ambiguous_count=ambiguous,
                tm=ff.compute_tm(seq),
            )
        )
    return baits


def write_baits_fasta(path: str, baits: List[ff.Bait]) -> None:
    ff.write_fasta(path, [(b.bait_id, b.seq.upper()) for b in baits])


def infer_bait_length_from_panel(baits: List[ff.Bait]) -> int:
    lengths = sorted({len(b.seq) for b in baits})
    if not lengths:
        raise RuntimeError("No bait sequences were found.")
    if len(lengths) != 1:
        raise RuntimeError(
            "FlyForgeAudit requires a uniform bait length. Observed lengths: "
            + ", ".join(str(x) for x in lengths)
        )
    return lengths[0]


def prepare_targets(
    input_paths: List[str],
    output_dir: str,
    prefix: str,
    bait_length: int,
    remove_complements: bool,
    skip_self_mask: bool,
    repeat_k: int,
    repeat_threshold: int,
    circular_all: bool,
    circular_ids_arg: Optional[str],
    log_fn,
) -> PreparedTargets:
    raw_seqs = read_fasta_paths(input_paths)
    modified_bases: Dict[str, int] = {}
    processed = {}
    for ref_id, seq in raw_seqs.items():
        proc, n_mod = ff.preprocess_sequence(seq)
        processed[ref_id] = proc
        modified_bases[ref_id] = n_mod

    work_fasta = os.path.join(output_dir, f"{prefix}_processed_targets.fna")
    ff.write_fasta(work_fasta, list(processed.items()))

    if remove_complements:
        records = list(SeqIO.parse(work_fasta, "fasta"))
        comp_prefix = os.path.join(output_dir, prefix)
        work_fasta = ff.remove_complementary_targets(
            records,
            output_dir=output_dir,
            prefix=comp_prefix,
            num_threads=1,
            probe_length=bait_length,
            log_fn=log_fn,
        )
        processed = {rec.id: str(rec.seq) for rec in SeqIO.parse(work_fasta, "fasta")}

    masked = processed if skip_self_mask else ff.self_repeat_softmask(
        processed, k=repeat_k, threshold=repeat_threshold
    )
    circular_ids = set(masked.keys()) if circular_all else ff.parse_circular_id_set(circular_ids_arg, list(masked.keys()))
    masked_fasta = os.path.join(output_dir, f"{prefix}_analysis_targets.fna")
    ff.write_fasta(masked_fasta, list(masked.items()))
    blast_fasta = os.path.join(output_dir, f"{prefix}_analysis_targets_for_blast.fna")
    build_extended_target_fasta(masked, bait_length, circular_ids, blast_fasta)

    return PreparedTargets(
        raw_seqs=raw_seqs,
        processed_seqs=processed,
        masked_seqs=masked,
        fasta_path=masked_fasta,
        blast_fasta_path=blast_fasta,
        circular_ids=circular_ids,
        modified_bases=modified_bases,
    )


# ============================================================================
# BLAST-backed panel analysis
# ============================================================================


def analyze_probe_panel(
    probe_fasta: str,
    target_fasta: str,
    output_dir: str,
    prefix: str,
    probe_length: int,
    log_fn,
    circular_ids: Optional[set] = None,
    target_lengths_override: Optional[Dict[str, int]] = None,
    write_outputs: bool = True,
    write_plots: bool = True,
) -> AnalysisResult:
    ensure_in_path(["blastn"])

    blast_xml = os.path.join(output_dir, f"{prefix}_final_blast.xml")
    cmd = [
        "blastn",
        "-query",
        probe_fasta,
        "-subject",
        target_fasta,
        "-out",
        blast_xml,
        "-outfmt",
        "5",
        "-dust",
        "no",
        "-soft_masking",
        "false",
    ]
    subprocess.run(cmd, check=True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)

    coverage_dict: Dict[str, np.ndarray] = {}
    gc_dict: Dict[str, float] = {}
    target_seq_dict: Dict[str, str] = {}
    target_count_dict: Dict[str, int] = {}
    probe_target_sets: Dict[str, set] = defaultdict(set)

    target_records = list(SeqIO.parse(target_fasta, "fasta"))
    probe_records = list(SeqIO.parse(probe_fasta, "fasta"))
    target_ids = [record.id for record in target_records]
    target_lengths = dict(target_lengths_override or {record.id: len(record.seq) for record in target_records})
    circular_ids = set(circular_ids or set())
    trusted_meta = ff.infer_trusted_probe_metadata(probe_records, target_lengths)

    for record in target_records:
        orig_len = int(target_lengths.get(record.id, len(record.seq)))
        coverage_dict[record.id] = np.zeros(orig_len, dtype=float)
        gc_dict[record.id] = SeqUtils.gc_fraction(str(record.seq)[:orig_len])
        target_seq_dict[record.id] = str(record.seq)[:orig_len]
        target_count_dict[record.id] = 0

    identity_cutoff = int(probe_length * 0.625)
    pair_list: List[Tuple[str, str]] = []
    seen_pairs = set()

    with open(blast_xml) as infile:
        for record in NCBIXML.parse(infile):
            probe = record.query.split()[0]
            meta = trusted_meta.get(probe)
            if meta is not None:
                meta = dict(meta)
                meta["target_id"] = ff.resolve_ref_hint_to_target(meta.get("ref_hint"), target_ids)

            valid_hits: List[Dict[str, object]] = []
            for alignment in record.alignments:
                target = getattr(alignment, "hit_id", None)
                if not target:
                    target = getattr(alignment, "hit_def", alignment.title).split()[0]
                if target not in coverage_dict:
                    continue
                for hsp in alignment.hsps:
                    if hsp.identities <= identity_cutoff:
                        continue
                    if re.search(r"--+", hsp.sbjct) is not None:
                        continue
                    if re.search(r"--+", hsp.query) is not None:
                        continue
                    valid_hits.append(
                        {
                            "target": target,
                            "hsp": hsp,
                            "identities": int(hsp.identities),
                            "align_length": int(hsp.align_length),
                            "bits": float(hsp.bits),
                            "subject_start": int(ff.hsp_subject_start(hsp)),
                        }
                    )

            assigned_target = None
            if meta is not None and meta.get("target_id") in coverage_dict and meta.get("start") is not None and meta.get("end") is not None:
                assigned_target = str(meta["target_id"])
                ff.apply_interval_coverage(coverage_dict[assigned_target], int(meta["start"]), int(meta["end"]))
            else:
                primary = ff.pick_primary_hit(valid_hits, meta)
                if primary is not None:
                    assigned_target = str(primary["target"])
                    match_array = ff.hsp_to_match_array(primary["hsp"])
                    target_len = coverage_dict[assigned_target].size
                    if assigned_target in circular_ids:
                        start1 = int(primary["subject_start"])
                        end1 = start1 + int(match_array.size) - 1
                        ff.apply_interval_coverage(coverage_dict[assigned_target], start1, end1)
                    else:
                        pad_left = int(primary["subject_start"]) - 1
                        pad_right = target_len - (pad_left + match_array.size)
                        if pad_right >= 0:
                            padded = np.pad(
                                match_array,
                                (pad_left, pad_right),
                                mode="constant",
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

    pair_df = pd.DataFrame(pair_list, columns=["target", "probe"])

    target_rows = []
    for target_name, cov_array in coverage_dict.items():
        target_rows.append(
            {
                "target_name": target_name,
                "length": cov_array.size,
                "gc": gc_dict[target_name],
                "sequence": target_seq_dict[target_name],
                "probe_count": target_count_dict[target_name],
                "proportion_covered": np.count_nonzero(cov_array) / cov_array.size if cov_array.size else 0,
                "mean_depth": float(np.mean(cov_array)) if cov_array.size else 0.0,
                "std_dev": float(np.std(cov_array)) if cov_array.size else 0.0,
                "min_depth": int(np.min(cov_array)) if cov_array.size else 0,
                "max_depth": int(np.max(cov_array)) if cov_array.size else 0,
            }
        )
    target_info = pd.DataFrame(target_rows)

    probe_rows = []
    for record in probe_records:
        seq = str(record.seq).upper()
        probe_rows.append(
            {
                "probe_id": record.id,
                "gc": SeqUtils.gc_fraction(seq),
                "tm": ff.compute_tm(seq),
                "num_targets": len(probe_target_sets.get(record.id, set())),
                "sequence": seq,
                "contains_lgui_site": int(ff.BSPQI_SITE_FWD in seq or ff.BSPQI_SITE_REV in seq),
            }
        )
    probe_info = pd.DataFrame(probe_rows)

    if write_outputs:
        pair_prefix = os.path.join(output_dir, prefix)
        pair_df.to_csv(f"{pair_prefix}_target_probe_pairs.csv", index=False)
        target_info.to_csv(f"{pair_prefix}_target_info.csv", index=False)
        probe_info.to_csv(f"{pair_prefix}_probe_info.csv", index=False)

    if write_plots:
        plot_dir = os.path.join(output_dir, f"{prefix}_plots")
        os.makedirs(plot_dir, exist_ok=True)
        indiv_dir = os.path.join(plot_dir, "individual_target_coverages")
        os.makedirs(indiv_dir, exist_ok=True)

        for target in coverage_dict:
            cov_array = coverage_dict[target]
            fig, ax = plt.subplots(figsize=(12, 4))
            sns.lineplot(x=np.arange(len(cov_array)), y=cov_array, ax=ax)
            ax.set_title(f"{target} Coverage")
            ax.set_xlabel("Position (bp)")
            ax.set_ylabel("Coverage Depth")
            safe_name = target.replace("/", "|").replace(" ", "_")
            plt.tight_layout()
            plt.savefig(os.path.join(indiv_dir, f"{safe_name}_coverage.png"), dpi=150)
            plt.savefig(os.path.join(indiv_dir, f"{safe_name}_coverage.svg"))
            plt.close(fig)

        def violin(df: pd.DataFrame, col: str, title: str, stem: str):
            if df.empty or col not in df.columns:
                return
            fig, ax = plt.subplots()
            sns.violinplot(y=df[col].astype(float), ax=ax)
            ax.set_title(title)
            plt.tight_layout()
            plt.savefig(os.path.join(plot_dir, f"{stem}.png"), dpi=300)
            plt.savefig(os.path.join(plot_dir, f"{stem}.svg"))
            plt.close(fig)

        violin(probe_info, "gc", "Probe GC Content Distribution", "probe_gc")
        violin(probe_info, "tm", "Probe Melting Temperature Distribution", "probe_tm")
        violin(probe_info, "num_targets", "Number of Targets per Probe", "probe_num_targets")
        violin(target_info, "length", "Target Length Distribution", "target_len")
        violin(target_info, "gc", "Target GC Content Distribution", "target_gc")
        violin(target_info, "probe_count", "Target Probe Count Distribution", "target_probe_count")
        violin(target_info, "proportion_covered", "Target Coverage Proportion Distribution", "target_coverage_prop")
        violin(target_info, "mean_depth", "Target Coverage Depth Distribution", "target_coverage_depth")
        violin(target_info, "std_dev", "Target Coverage Std Deviation Distribution", "target_coverage_stdev")

    return AnalysisResult(
        coverage_dict=coverage_dict,
        target_info=target_info,
        probe_info=probe_info,
        pair_df=pair_df,
        blast_xml=blast_xml,
    )


# ============================================================================
# Avoid-database screening and recommendations
# ============================================================================


def blast_against_avoid_db(
    probe_fasta: str,
    blast_db: str,
    output_tsv: str,
    min_pident: float,
    max_hits: int,
    threads: int,
) -> pd.DataFrame:
    ensure_in_path(["blastn"])
    env = os.environ.copy()
    db_dir = os.path.dirname(blast_db)
    if db_dir:
        env["BLASTDB"] = db_dir
    db_name = os.path.basename(blast_db)
    outfmt = "6 qseqid sseqid pident length nident qlen bitscore evalue sstrand"
    cmd = [
        "blastn",
        "-query",
        probe_fasta,
        "-db",
        db_name,
        "-out",
        output_tsv,
        "-outfmt",
        outfmt,
        "-num_threads",
        str(threads),
        "-max_target_seqs",
        str(max_hits),
        "-dust",
        "no",
        "-soft_masking",
        "false",
    ]
    subprocess.run(cmd, check=True, env=env, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)

    cols = outfmt.split()[1:]
    try:
        df = pd.read_csv(output_tsv, sep="\t", names=cols)
    except Exception:
        return pd.DataFrame(columns=cols)
    if df.empty:
        return df
    df = df.loc[df["pident"] >= min_pident].copy()
    df["identity_fraction"] = df["nident"] / df["qlen"]
    df.to_csv(output_tsv, sep="\t", index=False)
    return df


def build_recommendations(
    target_info: pd.DataFrame,
    probe_info: pd.DataFrame,
    bait_count: int,
    output_tsv: str,
    output_txt: str,
    desired_coverage: float,
    coverage_fraction_goal: float,
    avoid_hits: Optional[pd.DataFrame] = None,
) -> pd.DataFrame:
    recs: List[dict] = []

    for row in target_info.itertuples(index=False):
        if row.proportion_covered < coverage_fraction_goal:
            recs.append(
                {
                    "severity": "high",
                    "category": "coverage_gap",
                    "target": row.target_name,
                    "recommendation": "Add spike-in baits for uncovered or weakly covered regions.",
                    "evidence": f"Coverage fraction {row.proportion_covered:.3f} < goal {coverage_fraction_goal:.3f}",
                }
            )
        if row.mean_depth < desired_coverage:
            recs.append(
                {
                    "severity": "medium",
                    "category": "low_depth",
                    "target": row.target_name,
                    "recommendation": "Increase tiling density or add supplemental baits to raise depth.",
                    "evidence": f"Mean depth {row.mean_depth:.3f} < goal {desired_coverage:.3f}",
                }
            )
        if row.std_dev > max(row.mean_depth, 1.0):
            recs.append(
                {
                    "severity": "medium",
                    "category": "uneven_coverage",
                    "target": row.target_name,
                    "recommendation": "Inspect individual coverage plots and smooth sparse intervals with additional baits.",
                    "evidence": f"Coverage SD {row.std_dev:.3f} exceeds mean depth {row.mean_depth:.3f}",
                }
            )

    if not probe_info.empty:
        low_tm = probe_info.loc[probe_info["tm"] < 50]
        high_gc = probe_info.loc[(probe_info["gc"] > 0.8) | (probe_info["gc"] < 0.2)]
        lgui = probe_info.loc[probe_info.get("contains_lgui_site", 0) == 1]
        if len(low_tm) > 0:
            recs.append(
                {
                    "severity": "medium",
                    "category": "thermodynamics",
                    "target": "panel",
                    "recommendation": "Review low-Tm baits; these may hybridize weakly under stringent capture conditions.",
                    "evidence": f"{len(low_tm)} / {bait_count} probes have Tm < 50 C",
                }
            )
        if len(high_gc) > 0:
            recs.append(
                {
                    "severity": "low",
                    "category": "gc_outliers",
                    "target": "panel",
                    "recommendation": "Inspect extreme GC-content probes for capture bias or synthesis difficulty.",
                    "evidence": f"{len(high_gc)} / {bait_count} probes fall outside 20%-80% GC",
                }
            )
        if len(lgui) > 0:
            recs.append(
                {
                    "severity": "high",
                    "category": "opool_compatibility",
                    "target": "panel",
                    "recommendation": "Remove or redesign probes containing internal BspQI/LguI motifs before oligo-pool ordering.",
                    "evidence": f"{len(lgui)} probes contain {ff.BSPQI_SITE_FWD} or {ff.BSPQI_SITE_REV}",
                }
            )

    if avoid_hits is not None and not avoid_hits.empty:
        n_probes = avoid_hits["qseqid"].nunique()
        max_ident = avoid_hits["pident"].max()
        recs.append(
            {
                "severity": "high",
                "category": "off_target_hits",
                "target": "panel",
                "recommendation": "Review or replace probes with strong hits to the avoidance database.",
                "evidence": f"{n_probes} probes hit the avoid DB at >= configured % identity (max {max_ident:.2f}%)",
            }
        )

    if not recs:
        recs.append(
            {
                "severity": "info",
                "category": "no_action",
                "target": "panel",
                "recommendation": "No major audit flags were detected under the requested thresholds.",
                "evidence": "Coverage and probe QC metrics passed the configured review rules.",
            }
        )

    df = pd.DataFrame(recs)
    df.to_csv(output_tsv, sep="\t", index=False)
    with open(output_txt, "w") as out:
        out.write("# FlyForgeAudit recommendations\n")
        for row in df.itertuples(index=False):
            out.write(
                f"- [{row.severity}] {row.category} | {row.target} | {row.recommendation} | {row.evidence}\n"
            )
    return df


def write_per_ref_stats(
    output_path: str,
    target_info: pd.DataFrame,
    pair_df: pd.DataFrame,
    total_panel_baits: int,
    modified_bases: Optional[Dict[str, int]] = None,
) -> None:
    modified_bases = modified_bases or {}
    unique_probe_counts = (
        pair_df.groupby("target")["probe"].nunique().to_dict() if not pair_df.empty else {}
    )
    rows = []
    for row in target_info.itertuples(index=False):
        matched = int(unique_probe_counts.get(row.target_name, 0))
        rows.append(
            {
                "ref_id": row.target_name,
                "length_original": int(row.length),
                "length_preprocessed": int(row.length),
                "modified_bases": int(modified_bases.get(row.target_name, 0)),
                "frac_modified": float(modified_bases.get(row.target_name, 0) / row.length) if row.length else 0.0,
                "dropped_for_length": 0,
                "n_baits_tiled": matched,
                "n_baits_after_amb": matched,
                "n_baits_after_mask": matched,
                "n_baits_after_tm": matched,
                "n_baits_after_blast": matched,
                "n_baits_after_cluster": matched,
                "n_baits_after_redundancy": matched,
                "n_baits_final": matched,
                "coverage_mean": float(row.mean_depth),
                "coverage_min": int(row.min_depth),
                "coverage_max": int(row.max_depth),
                "coverage_fraction_covered": float(row.proportion_covered),
                "proportion_of_final_baits": matched / total_panel_baits if total_panel_baits else 0.0,
            }
        )
    pd.DataFrame(rows).to_csv(output_path, sep="\t", index=False)


def write_summary(
    output_path: str,
    mode: str,
    params: Dict[str, object],
    target_info: pd.DataFrame,
    probe_info: pd.DataFrame,
    recommendations: pd.DataFrame,
    extra_lines: Optional[List[str]] = None,
) -> None:
    actionable_flags, informational_notes = summarize_recommendations(recommendations)
    with open(output_path, "w") as out:
        out.write(f"# TackleBox: FlyForgeAudit v{__version__}\n")
        out.write(f"# Date: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        out.write(f"# Mode: {mode}\n")
        out.write("# ===== Parameters =====\n")
        for key, value in params.items():
            out.write(f"# {key}\t{value}\n")
        out.write("# ===== Headline metrics =====\n")
        out.write(
            f"targets\t{len(target_info)}\n"
            f"panel_baits\t{len(probe_info)}\n"
            f"mean_prop_covered\t{target_info['proportion_covered'].mean() if not target_info.empty else 0.0:.6f}\n"
            f"mean_depth\t{target_info['mean_depth'].mean() if not target_info.empty else 0.0:.6f}\n"
            f"n_recommendations\t{len(recommendations)}\n"
            f"actionable_flags\t{actionable_flags}\n"
            f"informational_notes\t{informational_notes}\n"
        )
        if extra_lines:
            out.write("# ===== Extra =====\n")
            for line in extra_lines:
                out.write(line + "\n")

def greedy_cover_candidates(
    ref_id: str,
    seq: str,
    current_cov: np.ndarray,
    target_depth: int,
    bait_len: int,
    omit_short_leq: int,
    pad_min: int,
    counter: Dict[str, int],
) -> List[ff.Bait]:
    seq = seq
    L = len(seq)
    candidates: List[ff.Bait] = []

    if L <= omit_short_leq:
        return candidates

    if L < bait_len:
        deficit_total = int(np.maximum(0, target_depth - current_cov).sum())
        if deficit_total > 0 and L >= pad_min:
            counter[ref_id] += 1
            padded = seq + ("T" * (bait_len - L))
            gc_frac, masked_frac, ambiguous = ff.compute_bait_metrics(padded)
            candidates.append(
                ff.Bait(
                    bait_id=f"{ref_id}|aug{counter[ref_id]}|pos1-{bait_len}",
                    seq=padded,
                    ref_id=ref_id,
                    ref_start=1,
                    ref_end=L,
                    gc_frac=gc_frac,
                    masked_frac=masked_frac,
                    ambiguous_count=ambiguous,
                    tm=ff.compute_tm(padded),
                )
            )
        return candidates

    deficits = np.maximum(0, target_depth - current_cov.astype(int)).astype(int)
    if deficits.sum() == 0:
        return candidates

    while deficits.sum() > 0:
        prefix = np.concatenate([[0], np.cumsum(deficits)])
        best_start = None
        best_score = 0
        for start in range(0, L - bait_len + 1):
            score = int(prefix[start + bait_len] - prefix[start])
            if score > best_score:
                best_start = start
                best_score = score
        if best_start is None or best_score <= 0:
            break
        frag = seq[best_start:best_start + bait_len]
        counter[ref_id] += 1
        gc_frac, masked_frac, ambiguous = ff.compute_bait_metrics(frag)
        candidates.append(
            ff.Bait(
                bait_id=f"{ref_id}|aug{counter[ref_id]}|pos{best_start + 1}-{best_start + bait_len}",
                seq=frag,
                ref_id=ref_id,
                ref_start=best_start + 1,
                ref_end=best_start + bait_len,
                gc_frac=gc_frac,
                masked_frac=masked_frac,
                ambiguous_count=ambiguous,
                tm=ff.compute_tm(frag),
            )
        )
        deficits[best_start:best_start + bait_len] = np.maximum(
            0, deficits[best_start:best_start + bait_len] - 1
        )

    return candidates


def deduplicate_baits_by_sequence(baits: List[ff.Bait]) -> List[ff.Bait]:
    seen = set()
    kept: List[ff.Bait] = []
    for bait in baits:
        seq = bait.seq.upper()
        if seq not in seen:
            seen.add(seq)
            kept.append(bait)
    return kept


def screen_against_existing_panel(
    candidates: List[ff.Bait],
    existing_baits: List[ff.Bait],
    output_dir: str,
    prefix: str,
    cluster_identity: float,
    cluster_overlap: float,
    threads: int,
) -> Tuple[List[ff.Bait], int]:
    if not candidates or not existing_baits:
        return candidates, 0

    existing_seqs = {b.seq.upper() for b in existing_baits}
    existing_rcs = {str(Seq(b.seq.upper()).reverse_complement()) for b in existing_baits}

    quick_kept = []
    dropped_ids = set()
    for bait in candidates:
        seq = bait.seq.upper()
        if seq in existing_seqs or seq in existing_rcs:
            dropped_ids.add(bait.bait_id)
        else:
            quick_kept.append(bait)

    if not quick_kept:
        return quick_kept, len(candidates)

    ensure_in_path(["blastn"])
    with tempfile.TemporaryDirectory() as tmpdir:
        cand_fa = os.path.join(tmpdir, "candidates.fna")
        existing_fa = os.path.join(tmpdir, "existing.fna")
        out_tsv = os.path.join(tmpdir, "candidate_vs_existing.tsv")
        write_baits_fasta(cand_fa, quick_kept)
        write_baits_fasta(existing_fa, existing_baits)
        outfmt = "6 qseqid length nident qlen sstrand"
        cmd = [
            "blastn",
            "-query",
            cand_fa,
            "-subject",
            existing_fa,
            "-out",
            out_tsv,
            "-outfmt",
            outfmt,
            "-num_threads",
            str(threads),
            "-dust",
            "no",
            "-soft_masking",
            "false",
        ]
        subprocess.run(cmd, check=True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
        try:
            df = pd.read_csv(out_tsv, sep="\t", names=outfmt.split()[1:])
        except Exception:
            df = pd.DataFrame(columns=outfmt.split()[1:])

    if not df.empty:
        plus_cut = df.loc[
            (df["sstrand"] == "plus")
            & (df["nident"] >= (df["qlen"] * cluster_identity))
            & (df["length"] >= (df["qlen"] * cluster_overlap)),
            "qseqid",
        ]
        minus_cut = df.loc[(df["sstrand"] == "minus") & (df["nident"] >= 30), "qseqid"]
        dropped_ids.update(set(plus_cut).union(set(minus_cut)))

    kept = [bait for bait in quick_kept if bait.bait_id not in dropped_ids]
    return kept, len(candidates) - len(kept)


def filter_new_candidates(
    candidates: List[ff.Bait],
    existing_baits: List[ff.Bait],
    args,
    output_dir: str,
    prefix: str,
    log_fn,
) -> Tuple[List[ff.Bait], List[str]]:
    messages: List[str] = []
    working = deduplicate_baits_by_sequence(candidates)
    messages.append(f"Initial proposed candidates: {len(working)}")

    working, removed = ff.filter_ambiguous(working, args.ambiguous_cutoff)
    messages.append(f"After ambiguous filter: {len(working)} ({removed} removed)")

    working, removed = ff.filter_masked_fraction(working, args.max_masked_frac)
    messages.append(f"After masked-fraction filter: {len(working)} ({removed} removed)")

    if args.min_tm > 0:
        working, removed = ff.filter_melting_temp(working, args.min_tm)
        messages.append(f"After Tm filter: {len(working)} ({removed} removed)")

    if not args.no_opool:
        working, removed = ff.filter_lgui_sites(working)
        messages.append(f"After LguI/BspQI filter: {len(working)} ({removed} removed)")

    working, removed = ff.filter_complementary_baits(working)
    messages.append(f"After perfect-complement filter: {len(working)} ({removed} removed)")

    working, removed = screen_against_existing_panel(
        working,
        existing_baits,
        output_dir=output_dir,
        prefix=prefix,
        cluster_identity=args.cluster_identity,
        cluster_overlap=args.cluster_overlap,
        threads=args.threads,
    )
    messages.append(f"After existing-panel cross-screen: {len(working)} ({removed} removed)")

    if args.blast_db and working:
        working, removed = ff.blast_filter_baits(
            working,
            blast_db=args.blast_db,
            evalue=args.blast_evalue,
            min_pident=args.blast_min_pident,
            max_hits_per_query=args.blast_max_hits,
            threads=args.threads,
        )
        messages.append(f"After BLAST exclusion filter: {len(working)} ({removed} removed)")

    if args.specificity_db and working:
        is_nt = os.path.basename(args.specificity_db) == "nt"
        working, removed = ff.blast_specificity_filter(
            working,
            blast_db=args.specificity_db,
            probe_length=args.bait_length,
            threads=args.threads,
            is_nt=is_nt,
        )
        messages.append(f"After specificity filter: {len(working)} ({removed} removed)")

    if not args.no_redundancy and working:
        working, comp_rm, red_rm = ff.self_blast_filter(
            working,
            output_dir,
            f"{prefix}_candidate_screen",
            args.threads,
            args.probe_num_cutoff,
            log_fn,
        )
        messages.append(
            f"After self-BLAST filter: {len(working)} ({comp_rm + red_rm} removed)"
        )

    if not args.no_cluster and working:
        working, removed = ff.cluster_baits_cd_hit(
            working,
            identity=args.cluster_identity,
            overlap=args.cluster_overlap,
            threads=args.threads,
        )
        messages.append(f"After clustering: {len(working)} ({removed} removed)")

    return working, messages


# ============================================================================
# Mode implementations
# ============================================================================


def run_audit(args) -> None:
    os.makedirs(args.output_dir, exist_ok=True)
    log_path = os.path.join(args.output_dir, f"{args.prefix}_progress.log")
    log, log_handle = log_factory(log_path)
    t_start = time.time()

    steps = ["preprocessing"]
    if args.remove_complements:
        steps.append("complement_removal")
    if not args.skip_self_mask:
        steps.append("self_masking")
    steps.append("validation")
    if args.avoid_db:
        steps.append("specificity_filter")
    steps.append("write_output")
    progress = ff.ProgressTracker(steps, log, t_start)

    log(f"TackleBox: FlyForgeAudit v{__version__} — audit mode")
    bait_objs = read_baits_as_objects(args.baits)
    bait_length = infer_bait_length_from_panel(bait_objs)
    panel_fasta = os.path.join(args.output_dir, f"{args.prefix}_final_baits.fa")
    write_baits_fasta(panel_fasta, bait_objs)

    progress.start_step("preprocessing")
    targets = prepare_targets(
        args.reference,
        output_dir=args.output_dir,
        prefix=args.prefix,
        bait_length=bait_length,
        remove_complements=args.remove_complements,
        skip_self_mask=args.skip_self_mask,
        repeat_k=args.repeat_k,
        repeat_threshold=args.repeat_threshold,
        circular_all=args.circular,
        circular_ids_arg=args.circular_ids,
        log_fn=log,
    )
    if args.remove_complements:
        progress.status["complement_removal"] = "done"
    if not args.skip_self_mask:
        progress.status["self_masking"] = "done"
    circular_detail = f"; circular: {', '.join(sorted(targets.circular_ids))}" if targets.circular_ids else ""
    progress.finish_step("preprocessing", details=f"{len(targets.masked_seqs)} targets prepared{circular_detail}")

    progress.start_step("validation")
    result = analyze_probe_panel(
        probe_fasta=panel_fasta,
        target_fasta=targets.blast_fasta_path,
        output_dir=args.output_dir,
        prefix=args.prefix,
        probe_length=bait_length,
        log_fn=log,
        circular_ids=targets.circular_ids,
        target_lengths_override={k: len(v) for k, v in targets.masked_seqs.items()},
        write_outputs=True,
        write_plots=True,
    )
    progress.finish_step("validation", details="Panel analysis and plotting complete")

    per_ref_path = os.path.join(args.output_dir, f"{args.prefix}_per_ref_stats.tsv")
    write_per_ref_stats(
        per_ref_path,
        result.target_info,
        result.pair_df,
        total_panel_baits=len(bait_objs),
        modified_bases=targets.modified_bases,
    )

    avoid_hits = None
    if args.avoid_db:
        progress.start_step("specificity_filter")
        avoid_path = os.path.join(args.output_dir, f"{args.prefix}_avoid_hits.tsv")
        avoid_hits = blast_against_avoid_db(
            panel_fasta,
            args.avoid_db,
            avoid_path,
            args.avoid_min_pident,
            args.avoid_max_hits,
            args.threads,
        )
        flagged = avoid_hits["qseqid"].nunique() if not avoid_hits.empty else 0
        log(f"Avoid-database screen complete: {flagged} probes flagged")
        progress.finish_step("specificity_filter", details=f"{flagged} probes flagged")

    rec_tsv = os.path.join(args.output_dir, f"{args.prefix}_recommendations.tsv")
    rec_txt = os.path.join(args.output_dir, f"{args.prefix}_recommendations.txt")
    recommendations = build_recommendations(
        result.target_info,
        result.probe_info,
        bait_count=len(bait_objs),
        output_tsv=rec_tsv,
        output_txt=rec_txt,
        desired_coverage=args.desired_coverage_depth,
        coverage_fraction_goal=args.coverage_fraction_goal,
        avoid_hits=avoid_hits,
    )

    progress.start_step("write_output")
    summary_path = os.path.join(args.output_dir, f"{args.prefix}_summary.tsv")
    params = {
        "mode": "audit",
        "baits": args.baits,
        "reference": ", ".join(args.reference),
        "bait_length": bait_length,
        "desired_coverage_depth": args.desired_coverage_depth,
        "coverage_fraction_goal": args.coverage_fraction_goal,
        "avoid_db": args.avoid_db or "N/A",
        "avoid_min_pident": args.avoid_min_pident,
        "avoid_max_hits": args.avoid_max_hits,
        "remove_complements": args.remove_complements,
        "skip_self_mask": args.skip_self_mask,
        "circular": args.circular,
        "circular_ids": ",".join(sorted(targets.circular_ids)) if targets.circular_ids else "N/A",
    }
    write_summary(summary_path, "audit", params, result.target_info, result.probe_info, recommendations)

    actionable_flags, informational_notes = summarize_recommendations(recommendations)
    rec_lines = format_recommendation_lines(recommendations)
    progress.finish_step("write_output", details="Reports written")

    elapsed = ff.format_duration(time.time() - t_start)
    mean_prop = result.target_info["proportion_covered"].mean() if not result.target_info.empty else 0.0
    mean_depth = result.target_info["mean_depth"].mean() if not result.target_info.empty else 0.0
    print(
        "\n" + "=" * 72 + "\n"
        f"  TackleBox: FlyForgeAudit v{__version__} — Audit Summary\n"
        + "=" * 72 + "\n"
        f"  Panel baits:             {len(bait_objs):,d}\n"
        f"  Targets audited:         {len(result.target_info):,d}\n"
        f"  Mean coverage fraction:  {mean_prop:.1%}\n"
        f"  Mean coverage depth:     {mean_depth:.2f}x\n"
        f"  Actionable flags:        {actionable_flags:,d}\n"
        f"  Informational notes:     {informational_notes:,d}\n"
        + "-" * 72 + "\n"
        + "  Recommendations\n"
        + "\n".join(rec_lines) + "\n"
        + "-" * 72 + "\n"
        + f"  Output directory:        {args.output_dir}\n"
        + f"  Runtime:                 {elapsed}\n"
        + "=" * 72,
        file=sys.stderr,
    )
    log_handle.close()

def run_augment(args) -> None:
    os.makedirs(args.output_dir, exist_ok=True)
    log_path = os.path.join(args.output_dir, f"{args.prefix}_progress.log")
    log, log_handle = log_factory(log_path)
    t_start = time.time()

    steps = ["preprocessing"]
    if args.remove_complements:
        steps.append("complement_removal")
    if not args.skip_self_mask:
        steps.append("self_masking")
    steps.append("tiling")
    if not args.no_opool:
        steps.append("opool_design")
    steps.append("validation")
    if args.avoid_db:
        steps.append("specificity_filter")
    steps.append("write_output")
    progress = ff.ProgressTracker(steps, log, t_start)

    log(f"TackleBox: FlyForgeAudit v{__version__} — augment mode")
    existing_baits = read_baits_as_objects(args.existing_baits)
    bait_length = infer_bait_length_from_panel(existing_baits)
    if args.bait_length is not None and args.bait_length != bait_length:
        raise RuntimeError(
            f"Requested --bait-length {args.bait_length} does not match existing panel length {bait_length}."
        )
    args.bait_length = bait_length

    existing_panel_fa = os.path.join(args.output_dir, f"{args.prefix}_existing_final_baits.fa")
    write_baits_fasta(existing_panel_fa, existing_baits)

    progress.start_step("preprocessing")
    targets = prepare_targets(
        args.new_targets,
        output_dir=args.output_dir,
        prefix=args.prefix,
        bait_length=bait_length,
        remove_complements=args.remove_complements,
        skip_self_mask=args.skip_self_mask,
        repeat_k=args.repeat_k,
        repeat_threshold=args.repeat_threshold,
        circular_all=args.circular,
        circular_ids_arg=args.circular_ids,
        log_fn=log,
    )
    if args.remove_complements:
        progress.status["complement_removal"] = "done"
    if not args.skip_self_mask:
        progress.status["self_masking"] = "done"
    circular_detail = f"; circular: {', '.join(sorted(targets.circular_ids))}" if targets.circular_ids else ""
    progress.finish_step("preprocessing", details=f"{len(targets.masked_seqs)} targets prepared{circular_detail}")

    counter = defaultdict(int)
    new_baits: List[ff.Bait] = []
    deficit_snapshots: List[str] = []

    progress.start_step("tiling")
    for iteration in range(1, args.max_augment_iterations + 1):
        combined = existing_baits + new_baits
        combined_tmp = os.path.join(args.output_dir, f"{args.prefix}_iter{iteration}_combined_tmp.fna")
        write_baits_fasta(combined_tmp, combined)

        iter_result = analyze_probe_panel(
            probe_fasta=combined_tmp,
            target_fasta=targets.blast_fasta_path,
            output_dir=args.output_dir,
            prefix=f"{args.prefix}_iter{iteration}",
            probe_length=bait_length,
            log_fn=log,
            circular_ids=targets.circular_ids,
            target_lengths_override={k: len(v) for k, v in targets.masked_seqs.items()},
            write_outputs=False,
            write_plots=False,
        )

        deficits_total = int(
            sum(np.maximum(0, args.min_existing_coverage - arr.astype(int)).sum() for arr in iter_result.coverage_dict.values())
        )
        mean_prop = iter_result.target_info["proportion_covered"].mean() if not iter_result.target_info.empty else 0.0
        deficit_snapshots.append(
            f"iteration_{iteration}\tpanel_baits={len(combined)}\tdeficit_bases={deficits_total}\tmean_prop_covered={mean_prop:.6f}"
        )
        log(f"Iteration {iteration}: deficit bases remaining = {deficits_total}")
        if deficits_total == 0:
            break

        proposals: List[ff.Bait] = []
        for ref_id, seq in targets.masked_seqs.items():
            cov = iter_result.coverage_dict.get(ref_id)
            if cov is None:
                cov = np.zeros(len(seq), dtype=float)
            proposals.extend(
                greedy_cover_candidates(
                    ref_id=ref_id,
                    seq=seq,
                    current_cov=cov,
                    target_depth=args.min_existing_coverage,
                    bait_len=bait_length,
                    omit_short_leq=args.omit_short_leq,
                    pad_min=args.pad_min,
                    counter=counter,
                )
            )
        if not proposals:
            log("No additional candidate windows could be proposed; stopping augmentation.")
            break

        filtered, messages = filter_new_candidates(
            proposals,
            existing_baits + new_baits,
            args,
            args.output_dir,
            f"{args.prefix}_iter{iteration}",
            log,
        )
        for msg in messages:
            log(f"Iteration {iteration}: {msg}")
        if not filtered:
            log("No candidate baits survived filtering; stopping augmentation.")
            break

        prior_n = len(new_baits)
        new_baits.extend(filtered)
        new_baits = deduplicate_baits_by_sequence(new_baits)
        if len(new_baits) == prior_n:
            log("No net new baits were added after deduplication; stopping augmentation.")
            break
    progress.finish_step("tiling", details=f"{len(new_baits)} extra baits proposed")

    extra_fasta = os.path.join(args.output_dir, f"{args.prefix}_extra_final_baits.fa")
    write_baits_fasta(extra_fasta, new_baits)

    extra_oligo_path = None
    extra_primer_path = None
    if not args.no_opool and new_baits:
        progress.start_step("opool_design")
        extra_oligo_path, _, extra_primer_path, _ = ff.design_opool(
            new_baits,
            args.output_dir,
            f"{args.prefix}_extra",
            log,
        )
        progress.finish_step("opool_design", details=f"Extra oligo pool: {os.path.basename(extra_oligo_path)}")

    merged_baits = existing_baits + new_baits
    merged_fasta = os.path.join(args.output_dir, f"{args.prefix}_final_baits.fa")
    write_baits_fasta(merged_fasta, merged_baits)

    progress.start_step("validation")
    final_result = analyze_probe_panel(
        probe_fasta=merged_fasta,
        target_fasta=targets.blast_fasta_path,
        output_dir=args.output_dir,
        prefix=args.prefix,
        probe_length=bait_length,
        log_fn=log,
        circular_ids=targets.circular_ids,
        target_lengths_override={k: len(v) for k, v in targets.masked_seqs.items()},
        write_outputs=True,
        write_plots=True,
    )
    progress.finish_step("validation", details="Merged panel analysis complete")

    per_ref_path = os.path.join(args.output_dir, f"{args.prefix}_per_ref_stats.tsv")
    write_per_ref_stats(
        per_ref_path,
        final_result.target_info,
        final_result.pair_df,
        total_panel_baits=len(merged_baits),
        modified_bases=targets.modified_bases,
    )

    avoid_hits = None
    if args.avoid_db:
        progress.start_step("specificity_filter")
        avoid_path = os.path.join(args.output_dir, f"{args.prefix}_avoid_hits.tsv")
        avoid_hits = blast_against_avoid_db(
            merged_fasta,
            args.avoid_db,
            avoid_path,
            args.avoid_min_pident,
            args.avoid_max_hits,
            args.threads,
        )
        if new_baits:
            blast_against_avoid_db(
                extra_fasta,
                args.avoid_db,
                os.path.join(args.output_dir, f"{args.prefix}_extra_avoid_hits.tsv"),
                args.avoid_min_pident,
                args.avoid_max_hits,
                args.threads,
            )
        flagged = avoid_hits["qseqid"].nunique() if not avoid_hits.empty else 0
        progress.finish_step("specificity_filter", details=f"{flagged} merged-panel probes flagged")

    rec_tsv = os.path.join(args.output_dir, f"{args.prefix}_recommendations.tsv")
    rec_txt = os.path.join(args.output_dir, f"{args.prefix}_recommendations.txt")
    recommendations = build_recommendations(
        final_result.target_info,
        final_result.probe_info,
        bait_count=len(merged_baits),
        output_tsv=rec_tsv,
        output_txt=rec_txt,
        desired_coverage=args.min_existing_coverage,
        coverage_fraction_goal=args.coverage_fraction_goal,
        avoid_hits=avoid_hits,
    )

    progress.start_step("write_output")
    summary_path = os.path.join(args.output_dir, f"{args.prefix}_summary.tsv")
    params = {
        "mode": "augment",
        "existing_baits": args.existing_baits,
        "new_targets": ", ".join(args.new_targets),
        "bait_length": bait_length,
        "min_existing_coverage": args.min_existing_coverage,
        "coverage_fraction_goal": args.coverage_fraction_goal,
        "max_augment_iterations": args.max_augment_iterations,
        "avoid_db": args.avoid_db or "N/A",
        "avoid_min_pident": args.avoid_min_pident,
        "avoid_max_hits": args.avoid_max_hits,
        "no_opool": args.no_opool,
        "remove_complements": args.remove_complements,
        "skip_self_mask": args.skip_self_mask,
        "circular": args.circular,
        "circular_ids": ",".join(sorted(targets.circular_ids)) if targets.circular_ids else "N/A",
        "extra_baits": len(new_baits),
        "merged_panel_baits": len(merged_baits),
    }
    extra_lines = deficit_snapshots + [
        f"extra_oligo_pool\t{extra_oligo_path or 'N/A'}",
        f"extra_amplification_primers\t{extra_primer_path or 'N/A'}",
    ]
    write_summary(summary_path, "augment", params, final_result.target_info, final_result.probe_info, recommendations, extra_lines=extra_lines)

    actionable_flags, informational_notes = summarize_recommendations(recommendations)
    rec_lines = format_recommendation_lines(recommendations)
    progress.finish_step("write_output", details="Augment reports written")

    elapsed = ff.format_duration(time.time() - t_start)
    mean_prop = final_result.target_info["proportion_covered"].mean() if not final_result.target_info.empty else 0.0
    mean_depth = final_result.target_info["mean_depth"].mean() if not final_result.target_info.empty else 0.0
    print(
        "\n" + "=" * 72 + "\n"
        f"  TackleBox: FlyForgeAudit v{__version__} — Augment Summary\n"
        + "=" * 72 + "\n"
        f"  Existing baits:          {len(existing_baits):,d}\n"
        f"  Extra baits designed:    {len(new_baits):,d}\n"
        f"  Merged panel size:       {len(merged_baits):,d}\n"
        f"  Targets audited:         {len(final_result.target_info):,d}\n"
        f"  Mean coverage fraction:  {mean_prop:.1%}\n"
        f"  Mean coverage depth:     {mean_depth:.2f}x\n"
        f"  Actionable flags:        {actionable_flags:,d}\n"
        f"  Informational notes:     {informational_notes:,d}\n"
        + "-" * 72 + "\n"
        + "  Recommendations\n"
        + "\n".join(rec_lines) + "\n"
        + "-" * 72 + "\n"
        + f"  Output directory:        {args.output_dir}\n"
        + f"  Runtime:                 {elapsed}\n"
        + "=" * 72,
        file=sys.stderr,
    )
    log_handle.close()

def add_shared_filtering_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--bait-length", type=int, default=None,
                        help="Expected bait length. In audit mode this is inferred from the panel if omitted.")
    parser.add_argument("--omit-short-leq", type=int, default=70,
                        help="Drop references with length <= this when proposing spike-in baits (default: 70).")
    parser.add_argument("--pad-min", type=int, default=71,
                        help="Pad references >= this and < bait_length to full bait length (default: 71).")
    parser.add_argument("--min-tm", type=float, default=50.0,
                        help="Minimum Tm filter for newly designed augment baits (default: 50).")
    parser.add_argument("--ambiguous-cutoff", type=int, default=10,
                        help="Remove baits with >= this many ambiguous bases (default: 10).")
    parser.add_argument("--max-masked-frac", type=float, default=0.25,
                        help="Maximum lowercase repeat-masked fraction (default: 0.25).")
    parser.add_argument("--repeat-k", type=int, default=15,
                        help="k-mer size for self-repeat masking (default: 15).")
    parser.add_argument("--repeat-threshold", type=int, default=3,
                        help="Minimum repeat-kmer count to trigger masking (default: 3).")
    parser.add_argument("--skip-self-mask", action="store_true",
                        help="Skip internal repeat masking of the references under review.")
    parser.add_argument("--remove-complements", action="store_true",
                        help="Run CARPDM/FlyForge-style complementary target cleanup before analysis.")
    parser.add_argument("--blast-db",
                        help="BLASTn DB for percent-identity off-target filtering when proposing augment baits.")
    parser.add_argument("--blast-evalue", type=float, default=1e-5,
                        help="BLAST e-value cutoff for augment bait filtering (default: 1e-5).")
    parser.add_argument("--blast-min-pident", type=float, default=90.0,
                        help="BLAST min percent identity to remove a proposed augment bait (default: 90).")
    parser.add_argument("--blast-max-hits", type=int, default=5,
                        help="BLAST max_target_seqs for proposed augment baits (default: 5).")
    parser.add_argument("--specificity-db",
                        help="CARPDM-style nident-based specificity filter DB for proposed augment baits.")
    parser.add_argument("--no-cluster", action="store_true",
                        help="Skip cd-hit-est clustering of proposed augment baits.")
    parser.add_argument("--cluster-identity", type=float, default=0.95,
                        help="cd-hit-est identity threshold for augment baits (default: 0.95).")
    parser.add_argument("--cluster-overlap", type=float, default=0.83,
                        help="cd-hit-est overlap threshold for augment baits (default: 0.83).")
    parser.add_argument("--no-redundancy", action="store_true",
                        help="Skip self-BLAST redundancy/complementarity filtering for augment baits.")
    parser.add_argument("--probe-num-cutoff", type=int, default=100000,
                        help="Max probe-count target for self-BLAST redundancy filter (default: 100000).")
    parser.add_argument("--avoid-db",
                        help="Curated avoidance BLAST database to screen the audited or augmented panel against.")
    parser.add_argument("--avoid-min-pident", type=float, default=80.0,
                        help="Minimum percent identity to flag an avoid-database hit (default: 80).")
    parser.add_argument("--avoid-max-hits", type=int, default=10,
                        help="Maximum BLAST hits per bait reported for avoid-database screening (default: 10).")
    parser.add_argument("--coverage-fraction-goal", type=float, default=0.95,
                        help="Coverage-fraction goal used in recommendations (default: 0.95).")
    parser.add_argument("--circular", action="store_true",
                        help="Treat all reference targets as circular for analysis/augmentation.")
    parser.add_argument("--circular-ids", default=None,
                        help="Comma-delimited subset of reference IDs to treat as circular.")
    parser.add_argument("--threads", type=int, default=1,
                        help="Threads for BLAST/cd-hit-est (default: 1).")


def main() -> None:
    ap = argparse.ArgumentParser(
        description="TackleBox: FlyForgeAudit - audit and augment existing bait panels.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Audit an existing bait set against one or more references
  FlyForgeAudit audit --baits existing_panel.fa --reference targets/*.fasta \
    --prefix panel_audit --output-dir audit_out --avoid-db /path/to/exclude_db

  # Add the minimal extra bait set needed to cover new organisms
  FlyForgeAudit augment --existing-baits existing_panel.fa --new-targets new_refs/*.fasta \
    --prefix panel_spikein --output-dir spikein_out --min-existing-coverage 1 --threads 8
        """,
    )
    sub = ap.add_subparsers(dest="mode", required=True)

    audit = sub.add_parser("audit", help="Audit an existing bait panel against reference FASTA files.")
    audit.add_argument("--baits", required=True, help="Existing bare-bait FASTA (no oligo-pool primers).")
    audit.add_argument("--reference", nargs="+", required=True, help="Reference FASTA file(s) to audit against.")
    audit.add_argument("--prefix", required=True, help="Prefix for output files.")
    audit.add_argument("--output-dir", default="flyforge_audit_output", help="Output directory.")
    audit.add_argument("--desired-coverage-depth", type=float, default=1.0,
                       help="Coverage depth goal used in recommendations (default: 1.0).")
    add_shared_filtering_args(audit)

    augment = sub.add_parser("augment", help="Design the minimal spike-in bait set for new targets.")
    augment.add_argument("--existing-baits", required=True, help="Existing bare-bait FASTA (no oligo-pool primers).")
    augment.add_argument("--new-targets", nargs="+", required=True, help="New target FASTA file(s) to add to the panel.")
    augment.add_argument("--prefix", required=True, help="Prefix for output files.")
    augment.add_argument("--output-dir", default="flyforge_augment_output", help="Output directory.")
    augment.add_argument("--min-existing-coverage", type=int, default=1,
                         help="Minimum desired coverage depth on the new targets after augmentation (default: 1).")
    augment.add_argument("--max-augment-iterations", type=int, default=5,
                         help="Maximum iterative augment rounds (default: 5).")
    augment.add_argument("--no-opool", action="store_true",
                         help="Skip order-ready extra oligo-pool generation.")
    add_shared_filtering_args(augment)

    ap.add_argument("-v", "--version", action="version",
                    version=f"TackleBox: FlyForgeAudit v{__version__}")
    args = ap.parse_args()

    if args.mode == "audit":
        run_audit(args)
    elif args.mode == "augment":
        run_augment(args)
    else:
        raise RuntimeError(f"Unsupported mode: {args.mode}")


if __name__ == "__main__":
    main()
