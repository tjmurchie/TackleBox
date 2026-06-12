"""Fillet batch pipeline — run the full metagenomics pipeline on a sequencing batch.

Reads a metadata Excel/TSV/CSV file describing all libraries, then runs each
through the complete pipeline (preprocess → screen → build DB → align → classify)
with numbered output directories, resumable progress tracking, and a live-updating
stats spreadsheet.

Two-phase shared-DB design (default, use_shared_db=True)
---------------------------------------------------------
1. Preprocess + Kraken2-screen ALL samples first.
2. Union the taxid hits across ALL samples, expand via taxonomy, run blastdbcmd ONCE.
3. Align + classify each sample against the single shared LAST database.

This means blastdbcmd (the main bottleneck for large taxid lists) runs exactly once
per batch rather than once per library.  Pass --no-shared-db to fall back to the
old per-library DB behaviour.

Usage::

    fillet batch \\
      --metadata samples.xlsx \\
      --outdir /results/batch1 \\
      --kraken-db /data/krakenuniq_db \\
      --blast-db-dir /data/NCBI_NT/2026-03-31 \\
      --nodes /data/taxonomy/nodes.dmp \\
      --names /data/taxonomy/names.dmp \\
      --threads 32 --resume -y

Directory layout::

    outdir/
      batch_progress.json       — per-sample status (pending/screened/complete/failed)
      batch_stats.tsv           — live-updated stats for all samples
      batch_stats.xlsx          — same, Excel format (if openpyxl installed)
      batch_metadata_copy.tsv   — copy of input metadata for audit trail
      shared_last_db/           — single shared LAST DB (two-phase mode)
        batch_shared.combined.fa
        batch_shared.suf  ...
        batch_union_taxids.txt
      01_NHB-7-Lib1/            — outputs for library 1
        preprocessing/
        screening/
        NHB-7-Lib1.last.maf
        NHB-7-Lib1.viewer.sqlite
      02_Blank-01/
      ...
      batch.merged.viewer.sqlite   — all samples merged (created at end)
      batch.merged_samples.tsv     — sample manifest for merged DB
"""
from __future__ import annotations

import json
import shutil
import sys
import time
import traceback
from contextlib import nullcontext
from pathlib import Path
from typing import Any, Dict, List, Optional, Set

from .batchmeta import load_batch_metadata, validate_batch_metadata, is_two_color_platform
from .batchstats import BatchStats
from .io import read_fossil_table_structured, filter_fossil_taxa_for_sample
from .metagenome import run_metagenome_pipeline, run_screening_phase
from .progress import PipelineDisplay, hms


# ---------------------------------------------------------------------------
# Progress state machine (JSON-backed)
# ---------------------------------------------------------------------------

class BatchProgress:
    """Simple JSON-backed state machine for per-library pipeline status.

    Valid statuses: pending → screened → complete | failed
    """

    def __init__(self, path: Path):
        self.path = path
        self._state: Dict[str, Any] = {}
        if path.exists():
            try:
                self._state = json.loads(path.read_text(encoding="utf-8"))
            except Exception:
                self._state = {}

    def status(self, lib_id: str) -> str:
        return self._state.get(lib_id, {}).get("status", "pending")

    def set_status(self, lib_id: str, status: str, **kwargs: Any) -> None:
        if lib_id not in self._state:
            self._state[lib_id] = {}
        self._state[lib_id]["status"] = status
        self._state[lib_id].update(kwargs)
        self._save()

    def sample_outdir(self, lib_id: str) -> Optional[str]:
        return self._state.get(lib_id, {}).get("outdir")

    def get_field(self, lib_id: str, field: str) -> Any:
        return self._state.get(lib_id, {}).get(field)

    def _save(self) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.path.write_text(
            json.dumps(self._state, indent=2, default=str), encoding="utf-8"
        )

    def summary(self) -> Dict[str, int]:
        counts: Dict[str, int] = {"pending": 0, "screened": 0, "running": 0, "complete": 0, "failed": 0}
        for s in self._state.values():
            key = s.get("status", "pending")
            counts[key] = counts.get(key, 0) + 1
        return counts


# ---------------------------------------------------------------------------
# Numbered directory helpers
# ---------------------------------------------------------------------------

def _numbered_dir(n: int, total: int, lib_id: str) -> str:
    """Return a zero-padded directory name like '007_NHB-7-Lib1'."""
    width = max(2, len(str(total)))
    safe = lib_id.replace("/", "_").replace("\\", "_").replace(" ", "_")
    return f"{n:0{width}d}_{safe}"


# ---------------------------------------------------------------------------
# Shared DB build helper
# ---------------------------------------------------------------------------

def _build_shared_db(
    *,
    all_raw_taxids: Set[str],
    nodes_path: str | Path,
    blast_db_dir: Optional[str | Path],
    custom_ref_fastas: Optional[List[str | Path]],
    shared_db_prefix: Path,
    threads: int,
    verbose: bool,
    max_seed_rank: str = "family",
    max_seqs_per_taxon: Optional[int] = None,
) -> Dict[str, Any]:
    """Expand taxids from all samples, run blastdbcmd once, build shared lastdb."""
    from .builddb import expand_taxids, filter_seeds_by_rank, extract_sequences_blastdbcmd, build_last_db
    import io as _io, logging

    shared_db_prefix.parent.mkdir(parents=True, exist_ok=True)
    fasta_path = shared_db_prefix.parent / (shared_db_prefix.name + ".combined.fa")
    union_taxid_file = shared_db_prefix.parent / "batch_union_taxids.txt"

    if verbose:
        print()
        print("  ── Building shared LAST database ──────────────────────────────")
        print(f"     Seed taxids (union): {len(all_raw_taxids):,}")
        print(f"     Output prefix:       {shared_db_prefix}", flush=True)

    log_buf = _io.StringIO()
    handler = logging.StreamHandler(log_buf)
    logger = logging.getLogger("fillet.builddb")
    logger.addHandler(handler)

    n_expanded = n_nt_seqs = n_custom_seqs = 0
    try:
        filtered_seeds, n_rank_dropped = filter_seeds_by_rank(
            all_raw_taxids, str(nodes_path), max_seed_rank
        )
        if verbose and n_rank_dropped:
            print(f"     Rank filter ({max_seed_rank}): dropped {n_rank_dropped:,} high-rank seeds, "
                  f"{len(filtered_seeds):,} remaining", flush=True)
        expanded = expand_taxids(filtered_seeds, str(nodes_path))
        n_expanded = len(expanded)
        if verbose:
            print(f"     After expansion:     {n_expanded:,}", flush=True)

        if blast_db_dir and expanded:
            if verbose:
                print("     Running blastdbcmd (this is the slow step) ...", flush=True)
            n_nt_seqs = extract_sequences_blastdbcmd(
                taxids=expanded,
                blast_db_dir=str(blast_db_dir),
                output_fasta=fasta_path,
                max_seqs_per_taxon=max_seqs_per_taxon,
                verbose=False,
            )
        else:
            fasta_path.write_text("", encoding="utf-8")

        if custom_ref_fastas:
            from .metagenome import _append_fasta
            for ref in custom_ref_fastas:
                n_custom_seqs += _append_fasta(ref, fasta_path)

        if verbose:
            print(f"     NT sequences:        {n_nt_seqs:,}", flush=True)
            print("     Building lastdb index ...", flush=True)

        build_last_db(
            input_fasta=fasta_path,
            output_prefix=shared_db_prefix,
            threads=threads,
            verbose=False,
        )
    finally:
        logger.removeHandler(handler)

    # Write the union taxid file so later runs can inspect it
    union_taxid_file.write_text("\n".join(sorted(all_raw_taxids)) + "\n", encoding="utf-8")

    if verbose:
        print(f"     ✓  Shared LAST DB ready: {shared_db_prefix}", flush=True)

    return {
        "taxids_expanded": n_expanded,
        "n_sequences": n_nt_seqs + n_custom_seqs,
        "shared_db_prefix": str(shared_db_prefix),
    }


# ---------------------------------------------------------------------------
# Main batch runner
# ---------------------------------------------------------------------------

def run_batch_pipeline(
    *,
    metadata_path: str | Path,
    outdir: str | Path,
    # Screening
    kraken_db: Optional[str | Path] = None,
    screener: str = "auto",
    min_kraken_reads: int = 5,
    min_kraken_unique_kmers: int = 50,
    max_kraken_dup_rate: Optional[float] = None,
    # DB build
    nodes_path: str | Path,
    blast_db_dir: Optional[str | Path] = None,
    last_db_prefix_base: Optional[str | Path] = None,
    skip_build_if_exists: bool = True,
    custom_ref_fastas: Optional[List[str | Path]] = None,
    use_shared_db: bool = True,
    max_seed_rank: str = "family",
    max_seqs_per_taxon: Optional[int] = None,
    chunked_db: bool = True,
    chunked_db_rank: Optional[str] = None,
    parallel_chunks: int = 1,
    chunk_timeout_sec: int = 1800,
    keep_chunk_mafs: bool = True,
    # Alignment
    last_train_file: Optional[str | Path] = None,
    evalue: str = "1e-5",
    max_target_seqs: int = 2000,
    chunk_size: int = 50000,
    threads: int = 16,
    # Classify
    taxonomy_tsv: Optional[str] = None,
    names_path: Optional[str] = None,
    root_taxid: str = "1",
    sample_sheet: Optional[str] = None,
    regional_taxa: Optional[str] = None,
    config: Optional[str] = None,
    em_iterations: int = 1,
    max_viewer_reads: int = 5000,
    sqlite_viewer: bool = True,
    # Preprocessing defaults (overridden per-sample by platform)
    min_length: int = 24,
    dedup: bool = True,
    adapter_check: bool = True,
    preprocess_tool: str = "auto",
    # Batch control
    resume: bool = True,
    merge_at_end: bool = True,
    verbose: bool = True,
    stop_on_error: bool = False,
    # Age/site-resolved fossil support
    fossil_table: Optional[str] = None,
) -> Dict[str, Any]:
    """Run the full metagenomics pipeline on every library in the metadata file.

    Parameters
    ----------
    use_shared_db : build one shared LAST database from the union of all Kraken2
        hits across the whole batch, then align every library against it.
        This runs blastdbcmd exactly once per batch instead of once per library.
        Default: True.  Pass False to build a separate DB per library (old behaviour).
    last_db_prefix_base : override the base directory for the shared (or per-sample)
        LAST DB.  Default: OUTDIR/shared_last_db/
    stop_on_error : if True, halt the batch on first failure. Default: continue.
    """
    metadata_path = Path(metadata_path)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # ── Load and validate metadata ───────────────────────────────────────────
    if verbose:
        print(f"\n  Loading metadata from {metadata_path} ...")
    samples = load_batch_metadata(metadata_path)
    if not samples:
        raise ValueError(f"No samples found in {metadata_path}")

    errors = validate_batch_metadata(samples)
    if errors:
        print("\n  Metadata validation errors:", file=sys.stderr)
        for e in errors:
            print(f"    ✗  {e}", file=sys.stderr)
        raise ValueError(f"{len(errors)} metadata error(s) — fix before running batch")

    # ── Load fossil records (age/site-resolved) ──────────────────────────────
    _fossil_records = read_fossil_table_structured(fossil_table) if fossil_table else []
    if _fossil_records and verbose:
        print(f"  Loaded {len(_fossil_records)} fossil occurrence record(s) from {fossil_table}")

    shutil.copy2(metadata_path, outdir / f"batch_metadata_copy{metadata_path.suffix}")

    n_samples = len(samples)
    progress = BatchProgress(outdir / "batch_progress.json")
    stats = BatchStats(samples)

    stats_json = outdir / "batch_stats.json"
    if resume and stats_json.exists():
        try:
            saved = BatchStats.load_from_json(stats_json)
            for lib_id in [s["library_id"] for s in samples]:
                if progress.status(lib_id) in ("complete", "screened"):
                    row = saved.get(lib_id)
                    if row:
                        stats.update(lib_id, row)
        except Exception:
            pass

    # ── Banner ────────────────────────────────────────────────────────────────
    if verbose:
        n_complete = sum(1 for s in samples if progress.status(s["library_id"]) == "complete")
        n_screened = sum(1 for s in samples if progress.status(s["library_id"]) == "screened")
        n_pending = n_samples - n_complete - n_screened
        db_mode = "shared (one DB per batch)" if use_shared_db else "per-library"
        print()
        print("╔" + "═" * 70 + "╗")
        print(f"║  Fillet Batch Metagenomics Pipeline {'':30}  ║")
        print(f"║  Metadata:  {str(metadata_path)[:55]:<55}  ║")
        print(f"║  Libraries: {n_samples:<4}  |  Pending: {n_pending:<4}  |  Complete: {n_complete:<4}  |  Threads: {threads:<4}  ║")
        print(f"║  DB mode:   {db_mode:<60}  ║")
        print("╚" + "═" * 70 + "╝")
        print(flush=True)

    batch_start = time.time()
    results: List[Dict[str, Any]] = []
    n_complete = 0
    n_failed = 0

    # ═══════════════════════════════════════════════════════════════════════════
    # TWO-PHASE SHARED-DB MODE
    # ═══════════════════════════════════════════════════════════════════════════
    if use_shared_db:
        db_base = Path(last_db_prefix_base) if last_db_prefix_base else outdir / "shared_last_db"
        shared_db_prefix = db_base / "batch_shared"

        # ── Phase 1: Preprocess + Screen all samples ───────────────────────────
        if verbose:
            print("  Phase 1: Preprocessing and screening all samples")
            print("  " + "─" * 60, flush=True)

        all_raw_taxids: Set[str] = set()
        screening_cache: Dict[str, Dict[str, Any]] = {}  # lib_id → {query, raw_taxids, ...}

        for idx, sample in enumerate(samples, start=1):
            lib_id = sample["library_id"]
            dir_name = _numbered_dir(idx, n_samples, lib_id)
            sample_outdir = outdir / dir_name
            cur_status = progress.status(lib_id)

            if resume and cur_status in ("screened", "complete"):
                # Already screened — reload cached taxids and query path
                cached_query = progress.get_field(lib_id, "preprocessed_query")
                cached_taxids_file = progress.get_field(lib_id, "taxids_file")
                if cached_query and cached_taxids_file and Path(cached_taxids_file).exists():
                    cached_taxids = set(Path(cached_taxids_file).read_text().splitlines())
                    cached_taxids.discard("")
                    all_raw_taxids.update(cached_taxids)
                    screening_cache[lib_id] = {
                        "query": cached_query,
                        "raw_taxids": cached_taxids,
                        "preprocess_stats": {},
                    }
                    if verbose:
                        mark = "✓" if cur_status == "complete" else "~"
                        print(f"  {mark}  [{idx}/{n_samples}]  {lib_id}  — screening already done "
                              f"({len(cached_taxids):,} taxids)")
                    continue
                # Fall through to re-run screening if cache is stale

            if cur_status == "complete":
                n_complete += 1
                if verbose:
                    print(f"  ✓  [{idx}/{n_samples}]  {lib_id}  — already complete")
                continue

            if verbose:
                print(f"\n  [{idx}/{n_samples}]  {lib_id}  — screening ...", flush=True)

            try:
                sample_poly_g = sample.get("two_color_chemistry", False)
                screen_result = run_screening_phase(
                    fastq_r1=sample.get("r1_path") or None,
                    fastq_r2=sample.get("r2_path") or None,
                    bam=sample.get("bam_path") or None,
                    preprocess_tool=preprocess_tool,
                    min_length=min_length,
                    dedup=dedup,
                    adapter_check=adapter_check,
                    poly_g_trim=sample_poly_g,
                    sample_id=lib_id,
                    outdir=sample_outdir,
                    kraken_db=kraken_db,
                    screener=screener,
                    min_kraken_reads=min_kraken_reads,
                    min_kraken_unique_kmers=min_kraken_unique_kmers,
                    max_kraken_dup_rate=max_kraken_dup_rate,
                    threads=threads,
                    verbose=verbose,
                )

                raw_taxids = screen_result["raw_taxids"]
                all_raw_taxids.update(raw_taxids)

                # Cache the taxid list for this sample so we can reload on resume
                sample_taxids_file = sample_outdir / "screening" / f"{lib_id}.screen_taxids.txt"
                sample_taxids_file.parent.mkdir(parents=True, exist_ok=True)
                sample_taxids_file.write_text("\n".join(sorted(raw_taxids)) + "\n", encoding="utf-8")

                screening_cache[lib_id] = screen_result
                progress.set_status(
                    lib_id, "screened",
                    outdir=str(sample_outdir),
                    preprocessed_query=screen_result["query"],
                    taxids_file=str(sample_taxids_file),
                    screened_at=time.strftime("%Y-%m-%dT%H:%M:%S"),
                )

                prep = screen_result.get("preprocess_stats", {})
                stats.update(lib_id, {
                    "raw_reads": prep.get("raw_reads"),
                    "after_polyg_trim": prep.get("after_polyg_trim"),
                    "after_adapter_trim": prep.get("after_trim"),
                    "after_length_filter": prep.get("after_length_filter"),
                    "after_adapter_check": prep.get("after_adapter_check"),
                    "after_dedup": prep.get("after_dedup"),
                    "pct_dedup": prep.get("pct_dedup"),
                    "mean_read_length": prep.get("mean_read_length"),
                    "reads_into_screening": prep.get("after_dedup"),
                    "taxids_raw": screen_result.get("taxids_from_screener"),
                })
                _write_batch_stats(stats, outdir, stats_json)

            except Exception as exc:
                n_failed += 1
                tb = traceback.format_exc()
                progress.set_status(lib_id, "failed",
                                    error=str(exc),
                                    failed_at=time.strftime("%Y-%m-%dT%H:%M:%S"))
                err_log = sample_outdir / "error.log"
                sample_outdir.mkdir(parents=True, exist_ok=True)
                err_log.write_text(tb, encoding="utf-8")
                if verbose:
                    print(f"  ✗  {lib_id}  screening FAILED: {exc}", file=sys.stderr)
                if stop_on_error:
                    raise

        if not all_raw_taxids:
            raise RuntimeError("No taxids found across any sample — check Kraken2 results.")

        # ── Build shared LAST DB (once for the whole batch) ────────────────────
        from .metagenome import _db_exists

        if skip_build_if_exists and _db_exists(shared_db_prefix):
            if verbose:
                print(f"\n  Shared LAST DB already exists — skipping build: {shared_db_prefix}")
            db_build_result: Dict[str, Any] = {"shared_db_prefix": str(shared_db_prefix)}
        else:
            if not blast_db_dir and not custom_ref_fastas:
                raise ValueError("--blast-db-dir or --custom-ref-fasta is required to build the LAST database")
            db_build_result = _build_shared_db(
                all_raw_taxids=all_raw_taxids,
                nodes_path=nodes_path,
                blast_db_dir=blast_db_dir,
                custom_ref_fastas=custom_ref_fastas,
                shared_db_prefix=shared_db_prefix,
                threads=threads,
                verbose=verbose,
                max_seed_rank=max_seed_rank,
                max_seqs_per_taxon=max_seqs_per_taxon,
            )

        # Write union taxid list for reference
        union_taxid_file = db_base / "batch_union_taxids.txt"
        if not union_taxid_file.exists():
            union_taxid_file.write_text("\n".join(sorted(all_raw_taxids)) + "\n", encoding="utf-8")

        # ── Phase 2: Align + Classify each sample ──────────────────────────────
        if verbose:
            print()
            print("  Phase 2: Aligning and classifying all samples")
            print("  " + "─" * 60, flush=True)

        for idx, sample in enumerate(samples, start=1):
            lib_id = sample["library_id"]
            dir_name = _numbered_dir(idx, n_samples, lib_id)
            sample_outdir = outdir / dir_name
            cur_status = progress.status(lib_id)

            if resume and cur_status == "complete":
                n_complete += 1
                if verbose:
                    print(f"  ✓  [{idx}/{n_samples}]  {lib_id}  — already complete")
                continue

            if cur_status == "failed" and lib_id not in screening_cache:
                if verbose:
                    print(f"  ✗  [{idx}/{n_samples}]  {lib_id}  — failed in Phase 1, skipping")
                continue

            screen_data = screening_cache.get(lib_id)
            if screen_data is None:
                if verbose:
                    print(f"  !  [{idx}/{n_samples}]  {lib_id}  — no screening data, skipping")
                continue

            preprocessed_query = screen_data["query"]

            if verbose:
                print(f"\n  ┌─ [{idx}/{n_samples}]  {lib_id}")
                print(f"  │   Output: {sample_outdir}")
                print("  │", flush=True)

            progress.set_status(lib_id, "running", started_at=time.strftime("%Y-%m-%dT%H:%M:%S"))
            stats.start(lib_id)

            try:
                sample_poly_g = sample.get("two_color_chemistry", False)

                _fos_taxa, _fos_ev = filter_fossil_taxa_for_sample(
                    _fossil_records, sample.get("site_name"), sample.get("age_BP")
                ) if _fossil_records else (None, [])

                sample_result = run_metagenome_pipeline(
                    # Supply preprocessed FASTA directly (skip Stage 0)
                    query=preprocessed_query,
                    sample_id=lib_id,
                    outdir=sample_outdir,
                    # Supply union taxid file (skip Stage 1 Kraken2)
                    taxid_file=str(union_taxid_file),
                    kraken_db=None,
                    # Point at pre-built shared DB (skip Stage 2)
                    nodes_path=nodes_path,
                    blast_db_dir=blast_db_dir,
                    last_db_prefix=shared_db_prefix,
                    skip_build_if_exists=True,
                    custom_ref_fastas=custom_ref_fastas or [],
                    max_seed_rank=max_seed_rank,
                    max_seqs_per_taxon=max_seqs_per_taxon,
                    chunked_db=chunked_db,
                    chunked_db_rank=chunked_db_rank,
                    parallel_chunks=parallel_chunks,
                    chunk_timeout_sec=chunk_timeout_sec,
                    keep_chunk_mafs=keep_chunk_mafs,
                    # Alignment + classification
                    last_train_file=last_train_file,
                    evalue=evalue,
                    max_target_seqs=max_target_seqs,
                    chunk_size=chunk_size,
                    threads=threads,
                    resume=resume,
                    taxonomy_tsv=taxonomy_tsv,
                    names_path=names_path,
                    root_taxid=root_taxid,
                    sample_sheet=sample_sheet,
                    regional_taxa=regional_taxa,
                    config=config,
                    em_iterations=em_iterations,
                    max_viewer_reads=max_viewer_reads,
                    no_viewer=False,
                    sqlite_viewer=sqlite_viewer,
                    fossil_taxa_override=_fos_taxa if _fos_taxa else None,
                    fos_evidence_texts_override=_fos_ev if _fos_ev else None,
                    verbose=verbose,
                )

                prep = screen_data.get("preprocess_stats", {})
                stats.update(lib_id, {
                    "raw_reads": prep.get("raw_reads"),
                    "after_polyg_trim": prep.get("after_polyg_trim"),
                    "after_adapter_trim": prep.get("after_trim"),
                    "after_length_filter": prep.get("after_length_filter"),
                    "after_adapter_check": prep.get("after_adapter_check"),
                    "after_dedup": prep.get("after_dedup"),
                    "pct_dedup": prep.get("pct_dedup"),
                    "mean_read_length": prep.get("mean_read_length"),
                    "reads_into_screening": prep.get("after_dedup"),
                    "taxids_raw": screen_data.get("taxids_from_screener"),
                    "taxids_expanded": db_build_result.get("taxids_expanded"),
                    "nt_sequences": db_build_result.get("n_sequences"),
                    "align_chunks": sample_result.get("align_chunks"),
                    "align_resumed": sample_result.get("align_resumed"),
                    "reads_assigned": sample_result.get("reads_assigned"),
                    "reads_unassigned": sample_result.get("reads_unassigned"),
                    "pct_assigned": sample_result.get("pct_assigned"),
                    "taxa_with_reads": sample_result.get("taxa_with_reads"),
                    "reads_with_damage": sample_result.get("reads_with_damage"),
                    "mean_damage_score": sample_result.get("mean_damage_score"),
                })
                stats.complete(lib_id)
                progress.set_status(lib_id, "complete",
                                    completed_at=time.strftime("%Y-%m-%dT%H:%M:%S"),
                                    sqlite=sample_result.get("sqlite_path"))
                n_complete += 1
                if verbose:
                    print(f"  └─ ✓  {lib_id}  complete")
                results.append({"library_id": lib_id, "status": "complete", **sample_result})

            except Exception as exc:
                n_failed += 1
                tb = traceback.format_exc()
                stats.complete(lib_id, failed=True)
                progress.set_status(lib_id, "failed",
                                    error=str(exc),
                                    failed_at=time.strftime("%Y-%m-%dT%H:%M:%S"))
                err_log = sample_outdir / "error.log"
                sample_outdir.mkdir(parents=True, exist_ok=True)
                err_log.write_text(tb, encoding="utf-8")
                if verbose:
                    print(f"  └─ ✗  {lib_id}  FAILED: {exc}", file=sys.stderr)
                if stop_on_error:
                    raise

            _write_batch_stats(stats, outdir, stats_json)

    # ═══════════════════════════════════════════════════════════════════════════
    # PER-LIBRARY DB MODE (legacy / --no-shared-db)
    # ═══════════════════════════════════════════════════════════════════════════
    else:
        for idx, sample in enumerate(samples, start=1):
            lib_id = sample["library_id"]
            dir_name = _numbered_dir(idx, n_samples, lib_id)
            sample_outdir = outdir / dir_name

            if resume and progress.status(lib_id) == "complete":
                n_complete += 1
                if verbose:
                    print(f"  ✓  [{idx:>{len(str(n_samples))}}/{n_samples}]  {lib_id}  — already complete")
                continue

            if last_db_prefix_base:
                db_prefix = Path(last_db_prefix_base) / lib_id / lib_id
            else:
                db_prefix = outdir / "shared_last_db" / lib_id / lib_id

            if verbose:
                print()
                print(f"  ┌─ [{idx}/{n_samples}]  {lib_id}")
                print(f"  │   Output: {sample_outdir}")
                print("  │", flush=True)

            progress.set_status(lib_id, "running", outdir=str(sample_outdir),
                                started_at=time.strftime("%Y-%m-%dT%H:%M:%S"))
            stats.start(lib_id)

            try:
                sample_poly_g = sample.get("two_color_chemistry", False)
                _fos_taxa, _fos_ev = filter_fossil_taxa_for_sample(
                    _fossil_records, sample.get("site_name"), sample.get("age_BP")
                ) if _fossil_records else (None, [])
                sample_result = run_metagenome_pipeline(
                    fastq_r1=sample.get("r1_path") or None,
                    fastq_r2=sample.get("r2_path") or None,
                    bam=sample.get("bam_path") or None,
                    preprocess_tool=preprocess_tool,
                    min_length=min_length,
                    dedup=dedup,
                    adapter_check=adapter_check,
                    poly_g_trim=sample_poly_g,
                    sample_id=lib_id,
                    outdir=sample_outdir,
                    kraken_db=kraken_db,
                    screener=screener,
                    min_kraken_reads=min_kraken_reads,
                    min_kraken_unique_kmers=min_kraken_unique_kmers,
                    max_kraken_dup_rate=max_kraken_dup_rate,
                    nodes_path=nodes_path,
                    blast_db_dir=blast_db_dir,
                    last_db_prefix=db_prefix,
                    skip_build_if_exists=skip_build_if_exists,
                    custom_ref_fastas=custom_ref_fastas or [],
                    max_seed_rank=max_seed_rank,
                    max_seqs_per_taxon=max_seqs_per_taxon,
                    chunked_db=chunked_db,
                    chunked_db_rank=chunked_db_rank,
                    parallel_chunks=parallel_chunks,
                    chunk_timeout_sec=chunk_timeout_sec,
                    keep_chunk_mafs=keep_chunk_mafs,
                    last_train_file=last_train_file,
                    evalue=evalue,
                    max_target_seqs=max_target_seqs,
                    chunk_size=chunk_size,
                    threads=threads,
                    resume=resume,
                    taxonomy_tsv=taxonomy_tsv,
                    names_path=names_path,
                    root_taxid=root_taxid,
                    sample_sheet=sample_sheet,
                    regional_taxa=regional_taxa,
                    config=config,
                    em_iterations=em_iterations,
                    max_viewer_reads=max_viewer_reads,
                    no_viewer=False,
                    sqlite_viewer=sqlite_viewer,
                    fossil_taxa_override=_fos_taxa if _fos_taxa else None,
                    fos_evidence_texts_override=_fos_ev if _fos_ev else None,
                    verbose=verbose,
                )

                prep = sample_result.get("preprocess_stats", {})
                stats.update(lib_id, {
                    "raw_reads": prep.get("raw_reads"),
                    "after_polyg_trim": prep.get("after_polyg_trim"),
                    "after_adapter_trim": prep.get("after_trim"),
                    "after_length_filter": prep.get("after_length_filter"),
                    "after_adapter_check": prep.get("after_adapter_check"),
                    "after_dedup": prep.get("after_dedup"),
                    "pct_dedup": prep.get("pct_dedup"),
                    "mean_read_length": prep.get("mean_read_length"),
                    "reads_into_screening": prep.get("after_dedup"),
                    "taxids_raw": sample_result.get("taxids_from_screener"),
                    "taxids_expanded": sample_result.get("taxids_expanded")
                        if isinstance(sample_result.get("taxids_expanded"), int) else None,
                    "nt_sequences": sample_result.get("n_sequences")
                        if isinstance(sample_result.get("n_sequences"), int) else None,
                    "align_chunks": sample_result.get("align_chunks"),
                    "align_resumed": sample_result.get("align_resumed"),
                    "reads_assigned": sample_result.get("reads_assigned"),
                    "reads_unassigned": sample_result.get("reads_unassigned"),
                    "pct_assigned": sample_result.get("pct_assigned"),
                    "taxa_with_reads": sample_result.get("taxa_with_reads"),
                    "reads_with_damage": sample_result.get("reads_with_damage"),
                    "mean_damage_score": sample_result.get("mean_damage_score"),
                })
                stats.complete(lib_id)
                progress.set_status(lib_id, "complete",
                                    completed_at=time.strftime("%Y-%m-%dT%H:%M:%S"),
                                    sqlite=sample_result.get("sqlite_path"))
                n_complete += 1
                if verbose:
                    print(f"  └─ ✓  {lib_id}  complete")
                results.append({"library_id": lib_id, "status": "complete", **sample_result})

            except Exception as exc:
                n_failed += 1
                tb = traceback.format_exc()
                stats.complete(lib_id, failed=True)
                progress.set_status(lib_id, "failed",
                                    error=str(exc),
                                    failed_at=time.strftime("%Y-%m-%dT%H:%M:%S"))
                err_log = sample_outdir / "error.log"
                sample_outdir.mkdir(parents=True, exist_ok=True)
                err_log.write_text(tb, encoding="utf-8")
                if verbose:
                    print(f"  └─ ✗  {lib_id}  FAILED: {exc}", file=sys.stderr)
                if stop_on_error:
                    raise

            _write_batch_stats(stats, outdir, stats_json)

    # ── Batch-level merge ────────────────────────────────────────────────────
    merged_sqlite = None
    if merge_at_end and n_complete > 0:
        sqlite_paths = []
        for s in samples:
            sdir_name = _numbered_dir(
                next(i for i, x in enumerate(samples, 1) if x["library_id"] == s["library_id"]),
                n_samples, s["library_id"]
            )
            p = outdir / sdir_name / f"{s['library_id']}.viewer.sqlite"
            if p.exists():
                sqlite_paths.append(str(p))

        if len(sqlite_paths) >= 2:
            from .viewer_sqlite import merge_viewer_dbs
            merged_sqlite = outdir / "batch.merged.viewer.sqlite"
            if verbose:
                print()
                print(f"  Merging {len(sqlite_paths)} SQLite databases → {merged_sqlite}")
            merge_result = merge_viewer_dbs(sqlite_paths, merged_sqlite)
            if verbose:
                samples_str = ", ".join(merge_result["samples"][:5])
                if len(merge_result["samples"]) > 5:
                    samples_str += f" + {len(merge_result['samples'])-5} more"
                print(f"  ✓  Merged: {samples_str}")

    # ── Final summary ─────────────────────────────────────────────────────────
    elapsed = time.time() - batch_start
    if verbose:
        print()
        print("  " + "━" * 68)
        print(f"  Batch complete  [{hms(elapsed)}]")
        print(f"  Libraries:  {n_samples}  |  Complete: {n_complete}  |  Failed: {n_failed}")
        print(f"  Stats:      {outdir / 'batch_stats.tsv'}")
        if merged_sqlite:
            print(f"  Merged DB:  {merged_sqlite}")
            print(f"  Open with:  fillet gui {merged_sqlite}")
        print("  " + "━" * 68)
        print()
        stats.print_summary()

    return {
        "n_samples": n_samples,
        "n_complete": n_complete,
        "n_failed": n_failed,
        "merged_sqlite": str(merged_sqlite) if merged_sqlite else None,
        "stats_tsv": str(outdir / "batch_stats.tsv"),
        "outdir": str(outdir),
    }


def _write_batch_stats(stats: BatchStats, outdir: Path, json_path: Path) -> None:
    stats.write_tsv(outdir / "batch_stats.tsv")
    stats.write_json(json_path)
    try:
        stats.write_excel(outdir / "batch_stats.xlsx")
    except Exception:
        pass
