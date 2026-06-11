"""External tool invocation helpers: BLAST, MMSeqs2, generic hooks."""
from __future__ import annotations

import csv
import os
import shutil
import subprocess
import tempfile
from typing import Dict, Optional

from .utils import info, warn


def run(cmd: list, log: Optional[str] = None, verbose: bool = True) -> str:
    """Run *cmd* as a subprocess, stream stdout/stderr to *log* if given.

    Raises ``RuntimeError`` on non-zero exit.
    """
    if verbose:
        info("Running: " + " ".join(str(x) for x in cmd))
    p = subprocess.run(cmd, text=True, capture_output=True)
    if log:
        with open(log, "a", encoding="utf-8") as fh:
            fh.write("$ " + " ".join(str(x) for x in cmd) + "\n")
            fh.write(p.stdout)
            fh.write(p.stderr)
            fh.write("\n")
    if p.returncode:
        raise RuntimeError(
            f"Command failed (exit {p.returncode}): {' '.join(str(x) for x in cmd)}\n{p.stderr}"
        )
    return p.stdout


def run_blast(query: str, db: str, out: str, cfgsec: dict) -> None:
    """Run blastn with parameters drawn from *cfgsec* (a config sub-dict)."""
    if not shutil.which("blastn"):
        raise RuntimeError(
            "blastn not found in PATH. Install NCBI BLAST+ and ensure blastn is executable."
        )
    cmd = [
        "blastn",
        "-query", query,
        "-db", db,
        "-task", cfgsec.get("blast_task", "megablast"),
        "-outfmt", cfgsec.get("outfmt", "6 qseqid sseqid pident length qstart qend evalue bitscore"),
        "-out", out,
    ]
    if "max_target_seqs" in cfgsec:
        cmd += ["-max_target_seqs", str(cfgsec["max_target_seqs"])]
    if "max_hsps" in cfgsec:
        cmd += ["-max_hsps", str(cfgsec["max_hsps"])]
    if "evalue" in cfgsec:
        cmd += ["-evalue", str(cfgsec["evalue"])]
    if cfgsec.get("num_threads", 1) > 1:
        cmd += ["-num_threads", str(int(cfgsec["num_threads"]))]
    run(cmd, out + ".log")


def _split_fasta_batches(path: str, batch_size: int):
    """Yield lists of (header, seq) pairs from *path* in chunks of *batch_size*."""
    batch: list = []
    with open(path, encoding="utf-8") as fh:
        header = None
        seq_parts: list = []
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    batch.append((header, "".join(seq_parts)))
                    if len(batch) >= batch_size:
                        yield batch
                        batch = []
                header = line
                seq_parts = []
            else:
                seq_parts.append(line)
        if header is not None:
            batch.append((header, "".join(seq_parts)))
    if batch:
        yield batch


def run_mmseqs(query: str, db: str, out: str, cfgsec: dict,
               batch_info: Optional[list] = None) -> None:
    """Run ``mmseqs easy-search`` with output columns identical to parse_tax_blast().

    The format string produces tab-separated columns in the same order as the
    BLAST outfmt used by run_blast():
      query  target  pident  alnlen  qlen  qstart  qend  evalue  bits  taxid  taxname

    Taxonomy columns (taxid, taxname) require the MMSeqs2 database to have been
    indexed with taxonomy (``mmseqs createtaxdb``).  Without taxonomy indexing
    these columns will be empty, and parse_tax_blast() will fall back to
    string-matching mode.

    Input is processed in batches of *batch_size* sequences (default 2000) to
    guard against prefilter segfaults seen in MMSeqs2 v16 when translating large
    nucleotide inputs against a protein database (search-type 2).  MMSeqs2 v18+
    handles arbitrarily large inputs without issue; raise batch_size or set it
    to 0 (unlimited) if you are on a fixed version that is not affected.

    Each completed batch is checkpointed to *out*.mmseqs_batches/ so that a
    failed run can resume from the last successful batch rather than restarting
    from scratch.  Pass a mutable ``[completed, total]`` list as *batch_info*
    to share batch progress with a BlastTicker running in another thread.

    Set ``mmseqs_binary`` in *cfgsec* to the full path of the mmseqs executable
    if it is not on PATH (e.g. when using a separate conda environment).
    """
    mmseqs_bin = cfgsec.get("mmseqs_binary") or shutil.which("mmseqs")
    if not mmseqs_bin:
        raise RuntimeError(
            "mmseqs not found in PATH. Install MMSeqs2 and ensure it is executable,"
            " or set mmseqs_binary in the taxonomy_blast config section."
        )
    # Use a local temp directory (defaults to /tmp) rather than a path next to
    # the output file, which may be on a network filesystem.  NFS + MMSeqs2's
    # multithreaded prefilter I/O causes consistent segfaults.
    local_tmp_root = cfgsec.get("tmp_dir", "") or tempfile.gettempdir()
    batch_size = int(cfgsec.get("batch_size", 2000))
    min_seq_id = float(cfgsec.get("min_pident", 70.0)) / 100.0
    pid = os.getpid()

    def _cmd(batch_query: str, batch_out: str, batch_tmp: str) -> list:
        c = [
            mmseqs_bin, "easy-search",
            batch_query, db, batch_out, batch_tmp,
            "--format-output",
            "query,target,pident,alnlen,qlen,qstart,qend,evalue,bits,taxid,taxname",
            "--min-seq-id", f"{min_seq_id:.4f}",
            "-e", str(cfgsec.get("evalue", "1e-10")),
            "--threads", str(int(cfgsec.get("num_threads", 1))),
        ]
        if cfgsec.get("max_target_seqs"):
            c += ["--max-seqs", str(cfgsec["max_target_seqs"])]
        if cfgsec.get("search_type") is not None:
            c += ["--search-type", str(int(cfgsec["search_type"]))]
        return c

    # Load all batches upfront to know total count for checkpointing and progress.
    all_batches = list(_split_fasta_batches(query, batch_size))
    n_batches = len(all_batches)
    total_seqs = sum(len(b) for b in all_batches)

    # Per-batch checkpoint dir — survives failures so runs can resume mid-step.
    ckpt_dir = out + ".mmseqs_batches"
    os.makedirs(ckpt_dir, exist_ok=True)

    def _ckpt(idx: int) -> str:
        return os.path.join(ckpt_dir, f"batch_{idx:04d}.tsv")

    completed = {i for i in range(n_batches) if os.path.exists(_ckpt(i))}
    n_done = len(completed)

    if n_done > 0:
        info(f"  Resuming mmseqs2: {n_done}/{n_batches} batches already complete,"
             f" continuing from batch {n_done + 1}")
    else:
        info(f"  Running {n_batches} mmseqs2 batches × {batch_size} seqs ({total_seqs:,} total)")

    if batch_info is not None:
        batch_info[0] = n_done
        batch_info[1] = n_batches

    batch_dirs: list = []
    success = False
    try:
        with open(out, "w", encoding="utf-8") as out_fh:
            # Pre-fill output with results from already-completed batches.
            for i in sorted(completed):
                with open(_ckpt(i), encoding="utf-8") as fh:
                    shutil.copyfileobj(fh, out_fh)
            out_fh.flush()

            for batch_idx, batch in enumerate(all_batches):
                if batch_idx in completed:
                    continue

                batch_tmp = os.path.join(local_tmp_root, f"mmseqs_spinner_{pid}_b{batch_idx}")
                os.makedirs(batch_tmp, exist_ok=True)
                batch_dirs.append(batch_tmp)
                batch_fasta = os.path.join(batch_tmp, "query.fasta")
                batch_tsv = os.path.join(batch_tmp, "results.tsv")

                with open(batch_fasta, "w", encoding="utf-8") as fh:
                    for header, seq in batch:
                        fh.write(header + "\n")
                        fh.write(seq + "\n")

                run(_cmd(batch_fasta, batch_tsv, batch_tmp),
                    out + f".batch{batch_idx}.log", verbose=False)

                # Atomic checkpoint: write to .tmp then rename so partial writes
                # are never mistaken for a completed batch on resume.
                cp = _ckpt(batch_idx)
                cp_tmp = cp + ".tmp"
                with open(cp_tmp, "w", encoding="utf-8") as tmp_fh:
                    if os.path.exists(batch_tsv):
                        with open(batch_tsv, encoding="utf-8") as in_fh:
                            shutil.copyfileobj(in_fh, tmp_fh)
                os.replace(cp_tmp, cp)

                with open(cp, encoding="utf-8") as fh:
                    shutil.copyfileobj(fh, out_fh)
                out_fh.flush()

                shutil.rmtree(batch_tmp, ignore_errors=True)
                batch_dirs = [d for d in batch_dirs if d != batch_tmp]

                if batch_info is not None:
                    batch_info[0] = batch_idx + 1

        success = True
    finally:
        for d in batch_dirs:
            shutil.rmtree(d, ignore_errors=True)
        if success:
            shutil.rmtree(ckpt_dir, ignore_errors=True)


def run_template_command(template: str, input_fasta: str, outprefix: str, label: str) -> Optional[str]:
    """Run an external hook command defined as a format string.

    Placeholders: ``{input_fasta}``, ``{outprefix}``, ``{label}``.
    Returns the path to the log file, or ``None`` if *template* is empty.
    """
    if not template:
        return None
    cmd = template.format(input_fasta=input_fasta, outprefix=outprefix, label=label)
    info("Running external hook: " + cmd)
    # shell=True is required here because the template may include pipes/redirects.
    proc = subprocess.run(cmd, shell=True, text=True, capture_output=True)  # noqa: S602
    log_path = f"{outprefix}.{label}.log"
    with open(log_path, "w", encoding="utf-8") as log:
        log.write("$ " + cmd + "\n")
        log.write(proc.stdout)
        log.write(proc.stderr)
    if proc.returncode:
        raise RuntimeError(
            f"{label} hook failed with exit code {proc.returncode}; see {log_path}"
        )
    return log_path


def parse_generic_hook_table(path: str, ann: Dict, hit_reason: str, review_reason: str) -> None:
    """Parse a simple TSV/CSV produced by an FCS or VecScreen postprocessor.

    The table must contain one of the columns ``accession``, ``qseqid``,
    ``query``, or ``record_key`` to identify the record.  Any cell containing
    the words reject / contaminant / adapter / vector / foreign triggers
    *hit_reason*; otherwise *review_reason* is added.

    This generic parser lets Spinner interoperate with changing external-tool
    output formats by adding a minimal postprocessing step that produces this
    normalised TSV.
    """
    if not path:
        return
    try:
        delimiter = "\t" if path.endswith(".tsv") else ","
        with open(path, encoding="utf-8", errors="replace") as f:
            reader = csv.DictReader(f, delimiter=delimiter)
            for row in reader:
                q = (
                    row.get("accession")
                    or row.get("qseqid")
                    or row.get("query")
                    or row.get("record_key")
                    or ""
                ).strip()
                if q not in ann:
                    continue
                text = " ".join(str(v).lower() for v in row.values())
                if any(x in text for x in ("reject", "contaminant", "adapter", "vector", "foreign")):
                    ann[q].add_reason(hit_reason)
                else:
                    ann[q].add_reason(review_reason)
    except FileNotFoundError:
        warn(f"Hook table not found: {path}")
