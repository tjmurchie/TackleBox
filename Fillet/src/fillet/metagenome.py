"""Fillet metagenomics pipeline: screen → build subset DB → align (LAST) → classify.

Typical usage::

    from fillet.metagenome import run_metagenome_pipeline

    run_metagenome_pipeline(
        query="reads.fasta",
        sample_id="NHB-7",
        outdir="/results/NHB-7",
        kraken_db="/data/krakenuniq_db",
        nodes_path="/data/taxonomy/nodes.dmp",
        names_path="/data/taxonomy/names.dmp",
        blast_db_dir="/data/NCBI_NT/2026-03-31",
        last_train_file="/data/LAST/nt_adna.train",
        custom_ref_fasta="/data/flyguide/local_taxa.fasta",  # optional
        threads=32,
    )

Why single-volume LAST is fast
-------------------------------
The full NT LAST database spans 106 volumes (≈110 GB each, ≈11.6 TB total).
lastal searches each volume in sequence even with -P32, so the NFS I/O of loading
11 TB of suffix arrays dominates runtime — not the alignment itself.

A subset DB extracted for KrakenUniq-identified taxa is a single volume that fits in
RAM, so lastal -P{threads} genuinely parallelises across all query chunks within
that one volume. In practice this gives ~30× speedup over the full NT DB on the same
hardware.
"""
from __future__ import annotations

import argparse
import shutil
import time
from contextlib import nullcontext
from pathlib import Path
from typing import Any, Dict, List, Optional, Set

from .align import count_fasta_records, run_chunked_alignment
from .builddb import expand_taxids, filter_seeds_by_rank, extract_sequences_blastdbcmd, build_last_db
from .progress import PipelineDisplay, hms


def _load_taxids_from_file(path: str | Path) -> Set[str]:
    taxids: Set[str] = set()
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            t = line.strip()
            if t and not t.startswith("#"):
                taxids.add(t)
    return taxids


def _db_exists(prefix: Path) -> bool:
    return (
        (prefix.with_suffix(".suf")).exists()
        or bool(list(prefix.parent.glob(prefix.name + "*.suf")))
        or (prefix.parent / "chunk_plan.json").exists()
    )


def _open_fasta(path: str | Path):
    p = str(path)
    if p.endswith(".gz"):
        import gzip
        return gzip.open(p, "rt", encoding="utf-8", errors="replace")
    return open(p, "r", encoding="utf-8", errors="replace")


def _append_fasta(src: str | Path, dst: str | Path) -> int:
    """Append sequences from src FASTA to dst (already open for appending). Returns record count."""
    n = 0
    with open(dst, "a", encoding="utf-8") as out:
        with _open_fasta(src) as fh:
            for line in fh:
                out.write(line)
                if line.startswith(">"):
                    n += 1
    return n


def run_screening_phase(
    *,
    fastq_r1: Optional[str | Path] = None,
    fastq_r2: Optional[str | Path] = None,
    bam: Optional[str | Path] = None,
    preprocess_tool: str = "auto",
    min_length: int = 24,
    dedup: bool = True,
    adapter_check: bool = True,
    poly_g_trim: bool = False,
    query: Optional[str | Path] = None,
    sample_id: str,
    outdir: str | Path,
    kraken_db: Optional[str | Path] = None,
    screener: str = "auto",
    min_kraken_reads: int = 5,
    min_kraken_unique_kmers: int = 50,
    max_kraken_dup_rate: Optional[float] = None,
    taxid_file: Optional[str | Path] = None,
    nodes_path: Optional[str | Path] = None,
    threads: int = 16,
    verbose: bool = True,
) -> Dict[str, Any]:
    """Run Stage 0 (preprocessing) and Stage 1 (Kraken2/KrakenUniq screening) only.

    Returns a dict with:
      ``query``           — path to the preprocessed FASTA (input to Stage 3)
      ``raw_taxids``      — set of taxid strings from the screener (or taxid_file)
      ``preprocess_stats``— preprocessing statistics dict (empty if query was given)
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    result: Dict[str, Any] = {"preprocess_stats": {}}

    do_preprocess = fastq_r1 is not None or bam is not None

    # ── Stage 0: Preprocessing ────────────────────────────────────────────────
    if do_preprocess:
        prep_outdir = outdir / "preprocessing"
        prep_log = outdir / "preprocessing.log"
        if verbose:
            print(f"  [{sample_id}] Stage 0: preprocessing ...", flush=True)
        from .preprocessing import preprocess_reads, find_preprocessor
        if not bam:
            find_preprocessor(preprocess_tool)
        from .preprocessing import preprocess_reads
        preprocessed_fasta, prep_stats = preprocess_reads(
            outdir=prep_outdir,
            prefix=sample_id,
            r1=fastq_r1,
            r2=fastq_r2,
            bam=bam,
            min_length=min_length,
            threads=threads,
            tool=preprocess_tool,
            dedup=dedup,
            adapter_check=adapter_check,
            poly_g_trim=poly_g_trim,
            log_path=prep_log,
        )
        query = preprocessed_fasta
        result["preprocess_stats"] = prep_stats
        if verbose:
            after = prep_stats.get("after_dedup", "?")
            print(f"  [{sample_id}]   After dedup: {after:,}" if isinstance(after, int)
                  else f"  [{sample_id}]   After dedup: {after}", flush=True)
    else:
        if query is None:
            raise ValueError(f"[{sample_id}] Provide --query, --fastq-r1, or --bam as input")
        query = Path(query)

    result["query"] = str(query)

    # ── Stage 1: Screening ────────────────────────────────────────────────────
    raw_taxids: Set[str] = set()
    if taxid_file is not None:
        raw_taxids = _load_taxids_from_file(taxid_file)
        result["taxids_from_screener"] = len(raw_taxids)
    else:
        if not kraken_db:
            raise ValueError(f"[{sample_id}] --kraken-db is required unless --taxid-file is given")
        screen_outdir = outdir / "screening"
        screen_log = outdir / "screening.log"
        if verbose:
            print(f"  [{sample_id}] Stage 1: Kraken2 screening ...", flush=True)
        from .krakenutils import screen_reads
        raw_taxids, report_path, tool_used = screen_reads(
            query=query,
            db=kraken_db,
            outdir=screen_outdir,
            sample_id=sample_id,
            threads=threads,
            screener=screener,
            min_reads=min_kraken_reads,
            min_unique_kmers=min_kraken_unique_kmers,
            max_dup_rate=max_kraken_dup_rate,
            log_path=screen_log,
        )
        if not raw_taxids:
            raise RuntimeError(
                f"[{sample_id}] No taxa found with ≥{min_kraken_reads} reads in {tool_used} report.\n"
                f"Check {report_path} and consider lowering --min-kraken-reads."
            )
        result["taxids_from_screener"] = len(raw_taxids)
        if verbose:
            print(f"  [{sample_id}]   Taxa found: {len(raw_taxids):,}", flush=True)

    result["raw_taxids"] = raw_taxids
    return result


def run_metagenome_pipeline(
    *,
    # ── FASTQ / BAM input (Stage 0) ──────────────────────────────────────────
    fastq_r1: Optional[str | Path] = None,
    fastq_r2: Optional[str | Path] = None,
    bam: Optional[str | Path] = None,
    preprocess_tool: str = "auto",
    min_length: int = 24,
    dedup: bool = True,
    adapter_check: bool = True,
    # ── Step-into controls ───────────────────────────────────────────────────
    from_stage: Optional[str] = None,   # "screen"|"build-db"|"align"|"classify"
    maf_path: Optional[str | Path] = None,  # explicit MAF when from_stage="classify"
    poly_g_trim: bool = False,
    # ── FASTA input (or output of Stage 0) ───────────────────────────────────
    query: Optional[str | Path] = None,
    sample_id: str,
    outdir: str | Path,
    # ── Stage 1: Taxonomic screening ─────────────────────────────────────────
    kraken_db: Optional[str | Path] = None,
    screener: str = "auto",
    min_kraken_reads: int = 5,
    min_kraken_unique_kmers: int = 50,
    max_kraken_dup_rate: Optional[float] = None,
    taxid_file: Optional[str | Path] = None,
    # ── Stage 2: Subset LAST DB build ────────────────────────────────────────
    nodes_path: str | Path,
    blast_db_dir: Optional[str | Path] = None,
    last_db_prefix: Optional[str | Path] = None,
    skip_build_if_exists: bool = True,
    custom_ref_fastas: Optional[List[str | Path]] = None,
    max_seed_rank: str = "family",
    max_seqs_per_taxon: Optional[int] = None,
    # Chunked DB options
    chunked_db: bool = True,
    chunked_db_rank: Optional[str] = None,  # None = auto-select
    parallel_chunks: int = 1,
    chunk_timeout_sec: int = 1800,
    keep_chunk_mafs: bool = True,
    # ── Stage 3: LAST alignment ───────────────────────────────────────────────
    last_train_file: Optional[str | Path] = None,
    last_min_score: Optional[int] = None,   # -e: min gapped alignment score; None = LAST auto
    last_m: Optional[int] = None,           # -m: max initial matches per query position; None = LAST default (10)
    evalue: str = "1e-5",
    max_target_seqs: int = 2000,
    max_hsps: int = 1,
    chunk_size: int = 50000,
    resume: bool = True,
    threads: int = 16,
    last_threads: Optional[int] = None,
    # ── Stage 4: RELIC-LCA classification ────────────────────────────────────
    taxonomy_tsv: Optional[str] = None,
    names_path: Optional[str] = None,
    root_taxid: str = "1",
    sample_sheet: Optional[str] = None,
    regional_taxa: Optional[str] = None,
    config: Optional[str] = None,
    em_iterations: int = 1,
    dirichlet_alpha: float = 0.0,
    damage_weight: float = 0.0,
    no_deconvolve: bool = False,
    use_mini_db: bool = True,
    max_clade_seqs_per_taxon: int = 5000,
    max_clade_db_mb: float = 50.0,
    mini_db_timeout_min: float = 30.0,
    max_parallel_mini_dbs: int = 1,
    curate_refs: bool = True,
    min_ref_length: int = 300,
    max_ref_n_fraction: float = 0.5,
    ref_scope: str = "all",
    ref_mode: str = "full",
    nt_index_path: Optional[str | Path] = None,
    custom_ref_taxid_map: Optional[str | Path] = None,
    dedup_similarity: float = 0.90,
    max_viewer_reads: int = 5000,
    no_viewer: bool = False,
    sqlite_viewer: bool = True,
    yes: bool = False,
    # ── RELIC-LCA CLI overrides (QUICK_LCA_FLAGS) ────────────────────────────
    # Mapping of config-key → value; None values are ignored.  Allows callers
    # such as cmd_metagenome to pass through --top-bitscore-fraction etc.
    lca_overrides: Optional[Dict[str, Any]] = None,
    # ── Classification mode profile ───────────────────────────────────────────
    mode: str = "auto",
    # ── UniWeight uniqueness index ────────────────────────────────────────────
    uniqueness_index: Optional[str | Path] = None,
    # ── Per-sample fossil support (site/age-resolved, injected by batch) ──────
    fossil_taxa_override: Optional[Set[str]] = None,
    fos_evidence_texts_override: Optional[List[str]] = None,
    # ── Display ───────────────────────────────────────────────────────────────
    verbose: bool = True,
) -> Dict[str, Any]:
    """Full metagenomics pipeline.

    Stages
    ------
    0. (Optional) FASTQ/BAM preprocessing — trim, merge PE reads, adapter check, dedup
    1. Kraken[Uniq] screening — fast first-pass over all reads
    2. Subset LAST DB — NT extraction + optional custom reference merge + lastdb index
    3. LAST alignment — multi-threaded MAF output with aligned qseq/sseq for damage
    4. RELIC-LCA classification — EM-refined LCA, damage detection, SQLite viewer

    Input
    -----
    Provide one of:
    - ``query`` — pre-processed FASTA (skip Stage 0)
    - ``fastq_r1`` / ``fastq_r2`` — raw paired-end FASTQ (fastp/AdapterRemoval)
    - ``fastq_r1`` alone — single-end or pre-merged FASTQ
    - ``bam`` — BAM from eager/BWA; reads extracted by samtools

    The subset DB is single-volume so lastal -P{threads} genuinely parallelises
    all query chunks (unlike the 106-volume full NT DB which forces sequential I/O).

    last_train_file
    ---------------
    If None, the bundled aDNA training file (nt_adna.train) is used automatically.
    It encodes substitution rates and gap penalties appropriate for aDNA
    (including C→T deamination at read ends), generated by last-train on aDNA data.

    Parameters
    ----------
    custom_ref_fastas : list of FASTA paths (from FlyGuide/Spinner or other curated
        sources) that are appended to the NT extraction before running lastdb.
        Useful for local/regional taxa not well-represented in NT.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Auto-select bundled training file if not provided
    if last_train_file is None:
        from .builddb import find_default_train_file
        last_train_file = find_default_train_file()

    do_preprocess = fastq_r1 is not None or bam is not None
    # Pure custom-ref mode: user supplies curated FASTA(s), no KrakenUniq, no NT extraction.
    _pure_custom_ref = (
        bool(custom_ref_fastas)
        and not blast_db_dir
        and not kraken_db
        and not taxid_file
    )
    skip_screening = (taxid_file is not None and kraken_db is None) or _pure_custom_ref
    do_screening = not skip_screening

    # ── from_stage: resolve skip level and validate prereqs ──────────────────
    _FROM_MAP = {"screen": 1, "build-db": 2, "align": 3, "classify": 4}
    _from = _FROM_MAP.get(from_stage, 0) if from_stage else 0

    if _from >= 1:
        # Resuming after preprocessing — query must be given or discoverable
        do_preprocess = False
        if query is None:
            _auto_q = Path(outdir) / "preprocessing" / f"{sample_id}.dedup.fasta"
            if _auto_q.exists():
                query = _auto_q
            else:
                raise ValueError(
                    f"--from-stage {from_stage!r} requires --query (preprocessed FASTA). "
                    f"Auto-discovery failed at {_auto_q}."
                )

    if _from >= 2:
        # Skip screening entirely; raw_taxids loaded from taxid_file if given
        skip_screening = True
        do_screening = False

    if _from >= 3:
        # Assume existing LAST DB; skip the build stage
        skip_build_if_exists = True

    # Pre-resolve last_db_prefix and MAF path (needed even for skipped stages)
    if last_db_prefix is None:
        last_db_prefix = Path(outdir) / "subset_last_db" / sample_id
    last_db_prefix = Path(last_db_prefix)

    _override_maf: Optional[Path] = None
    if _from == 4:
        _override_maf = Path(maf_path) if maf_path else Path(outdir) / f"{sample_id}.last.maf"
        if not _override_maf.exists():
            raise FileNotFoundError(
                f"--from-stage classify requires an existing MAF file.\n"
                f"Expected: {_override_maf}\n"
                f"Use --maf to specify an alternate path."
            )

    if _from == 3 and not _db_exists(last_db_prefix):
        # --from-stage align needs the DB; --from-stage classify only needs the MAF (checked above)
        raise FileNotFoundError(
            f"--from-stage {from_stage!r} requires an existing LAST DB at {last_db_prefix}.\n"
            "Run the full pipeline first or specify --last-db-prefix."
        )

    # Stage count: only count stages that will display
    if _from > 0:
        n_stages = sum([
            1 if _from <= 1 else 0,  # screening (or taxid load)
            1 if _from <= 2 else 0,  # DB build
            1 if _from <= 3 else 0,  # alignment
            1,                        # classification always
        ])
    else:
        n_stages = (1 if do_preprocess else 0) + (4 if do_screening else 3)

    # For banner: count reads if FASTA given, else show "?" until after preprocessing
    n_reads: Optional[int] = None
    if query is not None and not do_preprocess:
        try:
            n_reads = count_fasta_records(query)
        except Exception:
            pass

    reads_str = f"{n_reads:,}" if n_reads else "?"
    subtitle = f"Sample: {sample_id}  |  Reads: {reads_str}  |  Threads: {threads}"
    runtime_log = outdir / "pipeline_runtime.json"

    # Write config snapshot before PipelineDisplay loads the file so the
    # merge wizard can discover database paths and input files later.
    _pipeline_config = {
        "sample_id": sample_id,
        "outdir": str(outdir),
        "blast_db_dir": str(blast_db_dir) if blast_db_dir else None,
        "last_db_prefix": str(last_db_prefix) if last_db_prefix else None,
        "fastq_r1": str(fastq_r1) if fastq_r1 else None,
        "fastq_r2": str(fastq_r2) if fastq_r2 else None,
        "bam": str(bam) if bam else None,
        "query": str(query) if query else None,
        "kraken_db": str(kraken_db) if kraken_db else None,
        "maf_path": str(outdir / f"{sample_id}.last.maf"),
    }
    try:
        _prior = {}
        if runtime_log.exists():
            _prior = json.loads(runtime_log.read_text(encoding="utf-8"))
        _prior["config"] = _pipeline_config
        runtime_log.write_text(json.dumps(_prior, indent=2, default=str), encoding="utf-8")
    except Exception:
        pass

    disp = PipelineDisplay(n_stages=n_stages, title="Fillet Metagenomics Pipeline",
                           subtitle=subtitle, runtime_path=runtime_log)

    if verbose:
        disp.banner()
        if do_preprocess:
            if bam:
                disp.info(f"Input BAM: {bam}")
            else:
                disp.info(f"Input R1:  {fastq_r1}")
                if fastq_r2:
                    disp.info(f"Input R2:  {fastq_r2}")
        else:
            disp.info(f"Query:   {query}")
        disp.info(f"Output:  {outdir}")
        if last_train_file:
            disp.info(f"LAST train: {last_train_file}")
        if custom_ref_fastas:
            for ref in custom_ref_fastas:
                disp.info(f"Custom ref: {ref}")
        print(flush=True)

    result: Dict[str, Any] = {}
    stage = 0

    # ── STAGE 0: FASTQ / BAM preprocessing ───────────────────────────────────
    if do_preprocess:
        stage += 1
        prep_outdir = outdir / "preprocessing"
        prep_log = outdir / "preprocessing.log"
        prep_label = "BAM → FASTA" if bam else ("Paired-end FASTQ preprocessing" if fastq_r2 else "Single-end FASTQ preprocessing")

        if verbose:
            disp.stage_start(stage, prep_label)
            if bam:
                disp.info(f"Input:     {bam}")
            else:
                disp.info(f"R1:        {fastq_r1}")
                if fastq_r2:
                    disp.info(f"R2:        {fastq_r2}")
            disp.info(f"Min length:    {min_length} bp  |  Dedup: {dedup}  |  Adapter check: {adapter_check}")

        from .preprocessing import preprocess_reads, find_preprocessor

        if not bam:
            binary, tool_name = find_preprocessor(preprocess_tool)
            if verbose:
                disp.info(f"Tool:      {tool_name}")

        ctx = disp.spinner("Trimming, merging and deduplicating reads ...") if verbose else nullcontext()
        with ctx:
            from .preprocessing import preprocess_reads
            preprocessed_fasta, prep_stats = preprocess_reads(
                outdir=prep_outdir,
                prefix=sample_id,
                r1=fastq_r1,
                r2=fastq_r2,
                bam=bam,
                min_length=min_length,
                threads=threads,
                tool=preprocess_tool,
                dedup=dedup,
                adapter_check=adapter_check,
                poly_g_trim=poly_g_trim,
                log_path=prep_log,
            )

        query = preprocessed_fasta
        result["preprocess_stats"] = prep_stats
        if verbose:
            after_trim = prep_stats.get("after_trim", "?")
            after_dedup = prep_stats.get("after_dedup", "?")
            disp.stage_done(stage, prep_label, [
                f"After trim/merge:    {after_trim:,}" if isinstance(after_trim, int) else f"After trim:  {after_trim}",
                f"After adapter check: {prep_stats.get('after_adapter_check', '?')}",
                f"After dedup:         {after_dedup:,}" if isinstance(after_dedup, int) else f"After dedup: {after_dedup}",
                f"Output FASTA:        {query}",
                f"Log:                 {prep_log}",
            ])
    else:
        if query is None:
            raise ValueError("Provide --query (FASTA), --fastq-r1, or --bam as input")
        query = Path(query)
        # When using a pre-processed FASTA, apply the same adapter remnant check
        # that runs inside preprocess_reads for FASTQ inputs.
        if adapter_check and _from <= 0:
            from .preprocessing import filter_adapters_fasta
            prep_outdir = outdir / "preprocessing"
            prep_outdir.mkdir(parents=True, exist_ok=True)
            checked_fasta = prep_outdir / f"{sample_id}.adapter-checked.fasta"
            done_flag = prep_outdir / f"{sample_id}.adapter-checked.done"
            _n_removed = 0
            if not done_flag.exists():
                if verbose:
                    disp.info(f"Adapter remnant check on input FASTA ({query.name}) ...")
                _n_in, _n_removed = filter_adapters_fasta(query, checked_fasta)
                done_flag.write_text(
                    f"n_in={_n_in} n_removed={_n_removed}\n", encoding="utf-8"
                )
                if verbose:
                    if _n_removed > 0:
                        disp.info(
                            f"Adapter check: removed {_n_removed:,} reads with remnants "
                            f"({_n_in:,} → {_n_in - _n_removed:,})"
                        )
                    else:
                        disp.info(f"Adapter check: no remnants in {_n_in:,} reads")
            else:
                try:
                    _vals = dict(kv.split("=") for kv in done_flag.read_text(encoding="utf-8").split())
                    _n_removed = int(_vals.get("n_removed", 0))
                except Exception:
                    _n_removed = 0
                if verbose:
                    disp.info(f"Adapter check (resume): n_removed={_n_removed}")
            if _n_removed > 0 and checked_fasta.exists() and checked_fasta.stat().st_size > 0:
                query = checked_fasta

    # ── STAGE: Taxonomic screening / taxid loading ───────────────────────────
    raw_taxids: Set[str] = set()

    if _from <= 1:
        stage += 1

        if do_screening:
            if not kraken_db:
                raise ValueError("--kraken-db is required unless --taxid-file is given")
            screen_outdir = outdir / "screening"
            screen_log = outdir / "screening.log"

            if verbose:
                disp.stage_start(stage, "Taxonomic screening")
                tool_label = screener if screener != "auto" else "auto (krakenuniq preferred)"
                disp.info(f"Screener:  {tool_label}")
                disp.info(f"Database:  {kraken_db}")
                disp.info(f"Min reads: {min_kraken_reads}"
                          + (f"  |  Min unique k-mers: {min_kraken_unique_kmers}" if min_kraken_unique_kmers else "")
                          + (f"  |  Max dup rate: {max_kraken_dup_rate}" if max_kraken_dup_rate is not None else ""))

            from .krakenutils import screen_reads, find_screener
            _, tool_name = find_screener(screener)

            ctx = disp.spinner(f"Running {tool_name} ...") if verbose else nullcontext()
            with ctx:
                raw_taxids, report_path, tool_used = screen_reads(
                    query=query,
                    db=kraken_db,
                    outdir=screen_outdir,
                    sample_id=sample_id,
                    threads=threads,
                    screener=screener,
                    min_reads=min_kraken_reads,
                    min_unique_kmers=min_kraken_unique_kmers,
                    max_dup_rate=max_kraken_dup_rate,
                    log_path=screen_log,
                )

            if not raw_taxids:
                raise RuntimeError(
                    f"No taxa found with ≥{min_kraken_reads} reads in {tool_used} report.\n"
                    f"Check {report_path} and consider lowering --min-kraken-reads."
                )

            result["taxids_from_screener"] = len(raw_taxids)
            from .krakenutils import get_unclassified_fraction
            unclass_frac = get_unclassified_fraction(report_path)
            result["screener_unclassified_fraction"] = unclass_frac
            if verbose:
                disp.stage_done(stage, "Taxonomic screening", [
                    f"Tool:        {tool_used}",
                    f"Report:      {report_path}",
                    f"Taxa found:  {len(raw_taxids):,}  (≥{min_kraken_reads} reads)",
                    f"Unclassified: {unclass_frac:.1%}",
                    f"Log:         {screen_log}",
                ])
            if unclass_frac > 0.50:
                import warnings
                warnings.warn(
                    f"[{sample_id}] {unclass_frac:.1%} of reads were UNCLASSIFIED by {tool_used}. "
                    f"Organisms not covered by the screening database will be ABSENT from the "
                    f"LAST reference database and may generate false-positive hits to related taxa. "
                    f"Consider using a more comprehensive database (e.g. KrakenUniq NT) or "
                    f"supplying --taxid-file with known taxa of interest.",
                    stacklevel=2,
                )
        else:
            if verbose:
                disp.stage_start(stage, "Loading taxid list")
            if not taxid_file:
                raise ValueError("Provide --kraken-db for screening or --taxid-file to bypass it")
            raw_taxids = _load_taxids_from_file(taxid_file)
            result["taxids_from_screener"] = len(raw_taxids)
            if verbose:
                disp.stage_done(stage, "Loading taxid list", [
                    f"Source:  {taxid_file}",
                    f"Taxids:  {len(raw_taxids):,}",
                ])

    elif _from == 2:
        # from_stage=build-db: load taxids from file or auto-discover existing report
        if taxid_file:
            raw_taxids = _load_taxids_from_file(taxid_file)
        else:
            _auto_report = outdir / "screening" / f"{sample_id}.krakenuniq.report.txt"
            if not _auto_report.exists():
                _auto_report = outdir / "screening" / f"{sample_id}.kraken2.report.txt"
            if _auto_report.exists():
                from .krakenutils import parse_krakenuniq_report
                raw_taxids = parse_krakenuniq_report(
                    _auto_report, min_reads=min_kraken_reads,
                    min_unique_kmers=min_kraken_unique_kmers,
                    max_dup_rate=max_kraken_dup_rate,
                )
                if verbose:
                    print(f"  Auto-discovered screening report: {_auto_report.name} "
                          f"({len(raw_taxids):,} taxids)", flush=True)
            else:
                raise ValueError(
                    f"--from-stage build-db requires either --taxid-file or an existing "
                    f"screening report in {outdir}/screening/. No report found at {_auto_report}."
                )
        result["taxids_from_screener"] = len(raw_taxids)
    # else: from_stage >= align/classify — raw_taxids not needed (DB already exists)

    # ── STAGE: Build subset LAST database ────────────────────────────────────
    # last_db_prefix already resolved above (before stage guards)
    if _from <= 2:
        stage += 1

    if _from <= 2 and verbose:
        disp.stage_start(stage, "Building subset LAST database")
        disp.info(f"Output prefix: {last_db_prefix}")
        disp.info(f"BLAST DB dir:  {blast_db_dir or 'N/A'}")
        disp.info(f"Custom refs:   {len(custom_ref_fastas or [])}")
        disp.info(f"Nodes:         {nodes_path}")

    taxid_cache = last_db_prefix.parent / (last_db_prefix.name + ".taxids.txt")

    if _from >= 3:
        # DB existence already validated; skip stage display entirely
        result["taxids_expanded"] = "skipped"
        result["n_sequences"] = "skipped"
    elif skip_build_if_exists and _db_exists(last_db_prefix):
        if verbose:
            disp.info("LAST database already exists — skipping build.")
        cached_expanded = 0
        if taxid_cache.exists():
            cached_expanded = sum(1 for ln in taxid_cache.read_text().splitlines() if ln.strip())
        result["taxids_expanded"] = cached_expanded or "cached"
        result["n_sequences"] = "cached"
        if verbose:
            disp.stage_done(stage, "Building subset LAST database", ["Reused existing database."])
    else:
        if not blast_db_dir and not custom_ref_fastas:
            raise ValueError("--blast-db-dir or --custom-ref-fasta is required to build the LAST database")

        build_log = outdir / "builddb.log"
        fasta_path = last_db_prefix.parent / (last_db_prefix.name + ".combined.fa")
        fasta_path.parent.mkdir(parents=True, exist_ok=True)

        n_expanded = 0
        n_nt_sequences = 0
        n_custom_sequences = 0

        ctx = disp.spinner("Expanding taxonomy + extracting sequences + indexing LAST DB ...") if verbose else nullcontext()
        with ctx:
            import logging, io as _io

            log_buf = _io.StringIO()
            builddb_log = logging.getLogger("fillet.builddb")
            handler = logging.StreamHandler(log_buf)
            builddb_log.addHandler(handler)

            try:
                # Step 1: Filter seeds above max_seed_rank, then BFS expand
                filtered_seeds, n_rank_dropped = filter_seeds_by_rank(
                    raw_taxids, str(nodes_path), max_seed_rank
                )
                if n_rank_dropped and verbose:
                    disp.info(
                        f"Rank filter ({max_seed_rank}): dropped {n_rank_dropped:,} "
                        f"high-rank seeds, {len(filtered_seeds):,} remaining"
                    )
                expanded_taxids = expand_taxids(filtered_seeds, str(nodes_path))
                n_expanded = len(expanded_taxids)

                # Exclude taxids already covered by custom refs from NT extraction.
                # This prevents pulling lower-quality NT copies when a curated version exists.
                _excluded_taxids: Set[str] = set()
                if custom_ref_fastas and blast_db_dir:
                    from .builddb import extract_taxids_from_fasta
                    for _crf in custom_ref_fastas:
                        try:
                            _excluded_taxids.update(extract_taxids_from_fasta(Path(_crf)))
                        except Exception:
                            pass
                if custom_ref_taxid_map:
                    _excluded_taxids.update(_load_taxids_from_file(custom_ref_taxid_map))
                if _excluded_taxids:
                    _n_excl = len(expanded_taxids & _excluded_taxids)
                    expanded_taxids -= _excluded_taxids
                    n_expanded = len(expanded_taxids)
                    if verbose and _n_excl:
                        disp.info(
                            f"Custom-ref dedup: {_n_excl:,} taxid(s) excluded from NT extraction "
                            f"(already covered by custom reference FASTAs)"
                        )

                # Auto-suggest a per-taxon cap if none given, based on available RAM.
                _eff_max_seqs = max_seqs_per_taxon
                if _eff_max_seqs is None and n_expanded > 0:
                    from .builddb import available_memory_gb, suggest_max_seqs_per_taxon
                    _avail = available_memory_gb()
                    _eff_max_seqs = suggest_max_seqs_per_taxon(n_expanded, _avail)
                    if _eff_max_seqs is not None and verbose:
                        disp.info(
                            f"Memory auto-cap: {_avail:.0f} GB available → "
                            f"max {_eff_max_seqs:,} seqs/taxon "
                            f"(override with --max-seqs-per-taxon)"
                        )

                # Step 2: Extract matching NT sequences → combined.fa (non-chunked path only).
                # The chunked path (step 4a) extracts per-chunk directly; extracting to
                # combined.fa first and then ignoring it is ~2× wasted blastdbcmd time.
                _will_chunk = chunked_db and blast_db_dir and n_expanded > 0
                if blast_db_dir and expanded_taxids and not _will_chunk:
                    n_nt_sequences = extract_sequences_blastdbcmd(
                        taxids=expanded_taxids,
                        blast_db_dir=str(blast_db_dir),
                        output_fasta=fasta_path,
                        verbose=False,
                        max_seqs_per_taxon=_eff_max_seqs,
                    )
                elif not _will_chunk:
                    # No NT extraction — seed the file for custom-ref-only mode
                    fasta_path.write_text("", encoding="utf-8")

                # Step 3: Append custom reference FASTAs (FlyGuide/Spinner etc.).
                # For the chunked path, custom refs are injected at the chunk level (each
                # chunk already contains its taxon's NT sequences; custom refs go alongside).
                # If no chunked path: append to combined.fa as usual.
                if custom_ref_fastas and not _will_chunk:
                    for ref in custom_ref_fastas:
                        n_custom_sequences += _append_fasta(ref, fasta_path)

                # Apply reference scope filter (before curation: scope first, then quality)
                if ref_scope != "all" and not _will_chunk and fasta_path.exists() and fasta_path.stat().st_size > 0:
                    from .builddb import filter_fasta_by_scope
                    _scope_stats = filter_fasta_by_scope(fasta_path, ref_scope, verbose=False)
                    if verbose:
                        disp.info(
                            f"Ref scope ({ref_scope}): {_scope_stats['n_in']:,} → "
                            f"{_scope_stats['n_kept']:,} sequences "
                            f"({_scope_stats['n_removed']:,} removed)"
                        )

                # Curate combined FASTA (length filter, N-filter, dedup)
                if curate_refs and not _will_chunk and fasta_path.exists() and fasta_path.stat().st_size > 0:
                    from .builddb import curate_extracted_fasta
                    _curate_stats = curate_extracted_fasta(
                        fasta_path=fasta_path,
                        min_length=min_ref_length,
                        max_n_fraction=max_ref_n_fraction,
                        verbose=verbose,
                    )
                    if verbose:
                        disp.info(
                            f"Reference curation: {_curate_stats['n_in']:,} → "
                            f"{_curate_stats['n_out']:,} sequences "
                            f"({_curate_stats['n_removed_dup']:,} dupes, "
                            f"{_curate_stats['n_removed_length']:,} short, "
                            f"{_curate_stats['n_removed_n']:,} N-heavy removed)"
                        )

                # Step 4a (chunked): partition by rank, build N small DBs.
                # Step 4b (single):  legacy monolithic DB with preflight RAM check.
                if chunked_db and blast_db_dir and n_expanded > 0:
                    from .builddb import (
                        available_memory_gb, build_chunk_plan, save_chunk_plan,
                        build_chunked_db as _build_chunked_db, show_chunk_options_menu,
                        build_taxon_kingdom_map_from_nodes,
                    )
                    _avail_gb = available_memory_gb()
                    _chunks_dir = last_db_prefix.parent / "chunks"
                    chunk_plan_path = last_db_prefix.parent / "chunk_plan.json"

                    if chunk_plan_path.exists() and skip_build_if_exists and _from < 2:
                        from .builddb import load_chunk_plan
                        _chunk_plan = load_chunk_plan(chunk_plan_path)
                        if verbose:
                            disp.info(f"Loaded chunk plan: {len(_chunk_plan)} chunks from {chunk_plan_path}")
                    else:
                        _chunk_plan, _options_info = build_chunk_plan(
                            expanded_taxids=expanded_taxids,
                            nodes_path=nodes_path,
                            chunks_dir=_chunks_dir,
                            available_gb=_avail_gb,
                            threads=threads,
                            names_path=names_path,
                            force_rank=chunked_db_rank,
                        )
                        # Interactive menu when rank not forced and we have options
                        if not chunked_db_rank and _options_info.get("options"):
                            _sel_idx = show_chunk_options_menu(
                                _options_info,
                                yes=yes,
                                timeout_sec=chunk_timeout_sec,
                                n_expanded=n_expanded,
                                n_seeds=len(raw_taxids),
                            )
                            _sel_opt = _options_info["options"][_sel_idx]
                            # Re-build plan if user picked a different rank
                            auto_rank = _chunk_plan[0].chunk_id.split("_")[0] if _chunk_plan else ""
                            if _sel_opt["rank"] != auto_rank:
                                _chunk_plan, _ = build_chunk_plan(
                                    expanded_taxids=expanded_taxids,
                                    nodes_path=nodes_path,
                                    chunks_dir=_chunks_dir,
                                    available_gb=_avail_gb,
                                    threads=threads,
                                    names_path=names_path,
                                    force_rank=_sel_opt["rank"],
                                )
                        save_chunk_plan(_chunk_plan, chunk_plan_path)
                        if verbose:
                            disp.info(f"Chunk plan saved: {len(_chunk_plan)} chunk(s) → {chunk_plan_path}")

                    # Pre-compute per-taxon kingdom map so each taxon gets the
                    # appropriate reference type selection (plants→cp/mt only,
                    # animals→mt only, fungi→ITS-primary, etc.).
                    _taxon_kingdoms: Optional[Dict[str, str]] = None
                    if nodes_path:
                        _all_chunk_taxids: Set[str] = set()
                        for _cs in _chunk_plan:
                            _all_chunk_taxids.update(_cs.taxids)
                        _taxon_kingdoms = build_taxon_kingdom_map_from_nodes(
                            _all_chunk_taxids, nodes_path
                        )
                        if verbose:
                            from collections import Counter as _Counter
                            _kc = _Counter(_taxon_kingdoms.values())
                            disp.info(
                                f"[build-db] Kingdom map: "
                                + ", ".join(f"{k}={v:,}" for k, v in sorted(_kc.items()))
                            )

                    _build_chunked_db(
                        _chunk_plan,
                        blast_db_dir=str(blast_db_dir),
                        threads=threads,
                        parallel_chunks=parallel_chunks,
                        verbose=verbose,
                        max_seqs_per_taxon=_eff_max_seqs,
                        force_rebuild=(_from == 2),
                        custom_ref_fastas=custom_ref_fastas or [],
                        min_seq_length=min_ref_length,
                        ref_mode=ref_mode,
                        nt_index_path=nt_index_path,
                        dark_taxids_tsv=outdir / "dark_taxids.tsv",
                        dedup_similarity=dedup_similarity,
                        taxon_kingdoms=_taxon_kingdoms,
                        mask_cp_ir=(_taxon_kingdoms is not None),
                    )
                    _dark_tsv = outdir / "dark_taxids.tsv"
                    if _dark_tsv.exists() and _dark_tsv.stat().st_size > 0:
                        _n_dark = sum(1 for _l in open(_dark_tsv, encoding="utf-8") if _l.strip() and not _l.startswith("taxid"))
                        if verbose:
                            disp.info(
                                f"[build-db] Dark taxids: {_n_dark} taxid(s) had no qualifying sequences "
                                f"in the reference DB.\n"
                                f"           These reads will be assigned to the nearest represented ancestor.\n"
                                f"           List: {_dark_tsv}"
                            )
                    result["chunk_plan"] = str(chunk_plan_path)
                    result["n_chunks"] = len(_chunk_plan)
                else:
                    # Single monolithic DB path
                    _chunk_plan = None
                    if fasta_path.exists() and _eff_max_seqs is None:
                        from .builddb import (
                            available_memory_gb, memory_preflight_prompt,
                        )
                        _fasta_bytes = fasta_path.stat().st_size
                        _avail_gb = available_memory_gb()
                        _choice = memory_preflight_prompt(
                            fasta_size_bytes=_fasta_bytes,
                            available_gb=_avail_gb,
                            n_taxids=n_expanded,
                            threads=threads,
                            yes=yes,
                        )
                        if _choice == 0:
                            raise RuntimeError("Aborted by user at memory pre-flight check.")
                        elif _choice != -1 and _choice is not None:
                            from .builddb import _downsample_fasta_by_taxon
                            _tmp = fasta_path.with_suffix(".precap.tmp")
                            fasta_path.rename(_tmp)
                            n_nt_sequences = _downsample_fasta_by_taxon(
                                _tmp, fasta_path, _choice, verbose=True
                            )
                            if verbose:
                                disp.info(
                                    f"Applied cap: {_choice:,} seqs/taxon → "
                                    f"{fasta_path.stat().st_size / 1024**3:.1f} GB FASTA"
                                )
                    build_last_db(
                        input_fasta=fasta_path,
                        output_prefix=last_db_prefix,
                        threads=threads,
                        verbose=False,
                    )
            finally:
                builddb_log.removeHandler(handler)

        with open(build_log, "a", encoding="utf-8") as f:
            f.write(log_buf.getvalue())

        taxid_cache.parent.mkdir(parents=True, exist_ok=True)
        taxid_cache.write_text("\n".join(sorted(raw_taxids)) + "\n", encoding="utf-8")

        result["taxids_expanded"] = n_expanded
        result["n_sequences"] = n_nt_sequences + n_custom_sequences
        if verbose:
            lines = [
                f"Seed taxids:      {len(raw_taxids):,}",
                f"After expansion:  {n_expanded:,}",
            ]
            if blast_db_dir:
                lines.append(f"NT sequences:     {n_nt_sequences:,}")
            if custom_ref_fastas:
                lines.append(f"Custom seqs:      {n_custom_sequences:,}")
            lines += [
                f"Total sequences:  {n_nt_sequences + n_custom_sequences:,}",
                f"LAST DB prefix:   {last_db_prefix}",
                f"Log:              {build_log}",
            ]
            disp.stage_done(stage, "Building subset LAST database", lines)

    result["last_db_prefix"] = str(last_db_prefix)

    # ── STAGE: LAST alignment ────────────────────────────────────────────────
    # Pre-resolve MAF path (used by both alignment stage and classification)
    aligned = _override_maf or outdir / f"{sample_id}.last.maf"
    align_workdir = outdir / f"{sample_id}.align_work"

    if _from <= 3:
        stage += 1

    # Chunked alignment: use chunk_plan.json if present and chunked_db enabled.
    _chunk_plan_path = last_db_prefix.parent / "chunk_plan.json"
    _use_chunked_align = chunked_db and _chunk_plan_path.exists()

    if _from == 4:
        # Alignment skipped — MAF already validated above
        result["align_resumed"] = "skipped"
        if verbose:
            disp.info(f"Alignment skipped (--from-stage classify): using {aligned}")
    elif _use_chunked_align:
        from .builddb import load_chunk_plan, align_chunked, merge_chunked_mafs
        _loaded_chunk_plan = load_chunk_plan(_chunk_plan_path)

        if verbose:
            disp.stage_start(stage, f"LAST alignment — chunked ({len(_loaded_chunk_plan)} DB(s))")
            disp.info(f"Query:         {query}")
            disp.info(f"Chunk MAFs:    {align_workdir / 'chunk_mafs'}")
            print(flush=True)

        chunk_maf_paths = align_chunked(
            plan=_loaded_chunk_plan,
            query_fasta=Path(query),
            outdir=align_workdir,
            threads=threads,
            parallel_chunks=parallel_chunks,
            last_train_file=str(last_train_file) if last_train_file else None,
            last_min_score=last_min_score,
            last_m=last_m,
            evalue=evalue,
            max_target_seqs=max_target_seqs,
            chunk_size=chunk_size,
            resume=resume,
            keep_chunk_mafs=keep_chunk_mafs,
            verbose=verbose,
        )

        if verbose:
            disp.info(f"Merging {len(chunk_maf_paths)} chunk MAF(s) -> {aligned}")

        n_blocks = merge_chunked_mafs(
            chunk_maf_paths=chunk_maf_paths,
            output_maf=aligned,
            max_target_seqs=max_target_seqs,
            verbose=verbose,
        )

        result["align_chunks"] = len(chunk_maf_paths)
        result["align_resumed"] = 0
        result["align_blocks"] = n_blocks
        if verbose:
            disp.stage_done(stage, "LAST alignment — chunked", [
                f"Chunk MAFs:  {len(chunk_maf_paths):,}",
                f"Blocks kept: {n_blocks:,}",
                f"Output:      {aligned}",
            ])

    else:
        if verbose:
            disp.stage_start(stage, "LAST alignment (MAF, multi-threaded)")
            disp.info(f"LAST DB:       {last_db_prefix}  (single-volume — true -P{last_threads or threads} threading)")
            if last_train_file:
                disp.info(f"Training file: {last_train_file}")
            disp.info(f"Output:        {aligned}")
            disp.info(f"Chunk size:    {chunk_size:,} reads  |  Resume: {resume}")
            print(flush=True)

        align_result = run_chunked_alignment(
            aligner="last",
            query=query,
            output=aligned,
            workdir=align_workdir,
            last_db_prefix=str(last_db_prefix),
            last_train_file=str(last_train_file) if last_train_file else None,
            last_min_score=last_min_score,
            last_m=last_m,
            last_threads=last_threads or threads,
            threads=threads,
            evalue=evalue,
            max_target_seqs=max_target_seqs,
            max_hsps=max_hsps,
            chunk_size=chunk_size,
            resume=resume,
            last_output_format="maf",
            verbose=verbose,
        )

        result["align_chunks"] = align_result.chunks_total
        result["align_resumed"] = align_result.resumed_chunks
        if verbose:
            disp.stage_done(stage, "LAST alignment (MAF, multi-threaded)", [
                f"Chunks:    {align_result.chunks_total:,}  (resumed: {align_result.resumed_chunks:,})",
                f"Output:    {aligned}",
            ])

    # ── STAGE: RELIC-LCA classification ──────────────────────────────────────
    stage += 1
    if verbose:
        disp.stage_start(stage, "RELIC-LCA classification")
        print(flush=True)

    from .cli import cmd_classify, QUICK_LCA_FLAGS

    classify_args = argparse.Namespace(
        input=[str(aligned)],
        sample_id=[sample_id],
        columns=None,
        query_fasta=[str(query)],
        acc2taxid=None,
        nodes=str(nodes_path) if not taxonomy_tsv else None,
        names=names_path,
        taxonomy_tsv=taxonomy_tsv,
        root_taxid=root_taxid,
        sample_sheet=sample_sheet,
        regional_taxa=regional_taxa,
        contaminants_file=None,
        palynology_table=None,
        fossil_table=None,
        depth_file=None,
        reroute_taxids=[],
        config=config,
        outdir=str(outdir),
        prefix=sample_id,
        max_viewer_reads=max_viewer_reads,
        no_viewer=no_viewer,
        sqlite_viewer=sqlite_viewer,
        verbose=verbose,
        em_iterations=em_iterations,
        dirichlet_alpha=dirichlet_alpha,
        damage_weight=damage_weight,
        no_deconvolve=no_deconvolve,
        no_mini_db=not use_mini_db,
        max_clade_seqs_per_taxon=max_clade_seqs_per_taxon,
        max_clade_db_mb=max_clade_db_mb,
        mini_db_timeout_min=mini_db_timeout_min,
        max_parallel_mini_dbs=max_parallel_mini_dbs,
        blast_db_dir=str(blast_db_dir) if blast_db_dir else None,
        nt_index=str(nt_index_path) if nt_index_path else None,
        mode=mode,
        uniqueness_index=str(uniqueness_index) if uniqueness_index else None,
    )
    # Inject per-sample fossil support so cmd_classify can skip file-based loading.
    if fossil_taxa_override is not None:
        classify_args.__fossil_taxa_override__ = fossil_taxa_override
        classify_args.__fos_evidence_texts_override__ = fos_evidence_texts_override or []
    # Apply LCA overrides from the caller (e.g. CLI flags passed through cmd_metagenome).
    # Values not in lca_overrides (or None values) fall back to config/TOML defaults.
    _overrides = lca_overrides or {}
    for _flag, _section, key, _caster in QUICK_LCA_FLAGS:
        setattr(classify_args, key, _overrides.get(key, None))
    # Boolean flags that live outside QUICK_LCA_FLAGS (handled by apply_config_overrides).
    classify_args.suppress_unclassified_nodes = bool(_overrides.get("suppress_unclassified_nodes", False))
    classify_args.cross_clade_graduated_off = bool(_overrides.get("cross_clade_graduated_off", False))

    cmd_classify(classify_args)

    sqlite_path = outdir / f"{sample_id}.viewer.sqlite"
    result["sqlite_path"] = str(sqlite_path) if sqlite_path.exists() else None

    # Parse classify output for BatchStats
    assign_path = outdir / f"{sample_id}.read_assignments.tsv"
    if assign_path.exists():
        import csv as _csv
        reads_assigned = reads_unassigned = reads_with_damage = 0
        damage_scores: List[float] = []
        assigned_taxids: Set[str] = set()
        with assign_path.open("r", encoding="utf-8") as _fh:
            for row in _csv.DictReader(_fh, delimiter="\t"):
                if row.get("assigned_taxid", "0") not in ("0", ""):
                    reads_assigned += 1
                    assigned_taxids.add(row["assigned_taxid"])
                else:
                    reads_unassigned += 1
                ds = float(row.get("damage_score") or 0)
                if ds > 0:
                    reads_with_damage += 1
                    damage_scores.append(ds)
        total_reads = reads_assigned + reads_unassigned
        result["reads_assigned"] = reads_assigned
        result["reads_unassigned"] = reads_unassigned
        result["pct_assigned"] = 100.0 * reads_assigned / max(1, total_reads)
        result["taxa_with_reads"] = len(assigned_taxids)
        result["reads_with_damage"] = reads_with_damage
        result["mean_damage_score"] = sum(damage_scores) / len(damage_scores) if damage_scores else 0.0

    if verbose:
        disp.stage_done(stage, "RELIC-LCA classification", [
            f"SQLite viewer:  {sqlite_path}" if sqlite_path.exists() else "Classification complete.",
            f"Serve with:     fillet serve-viewer --db {sqlite_path}" if sqlite_path.exists() else "",
        ])

        n_screener = result.get("taxids_from_screener")
        n_expanded = result.get("taxids_expanded")
        n_seq = result.get("n_sequences")
        n_assigned = result.get("reads_assigned")
        n_taxa = result.get("taxa_with_reads")
        mean_dmg = result.get("mean_damage_score")
        disp.pipeline_done([
            f"Screener taxa:    {n_screener:,}" if isinstance(n_screener, int) else "",
            f"Expanded taxids:  {n_expanded:,}" if isinstance(n_expanded, int) else "",
            f"Sequences in DB:  {n_seq:,}" if isinstance(n_seq, int) else "",
            f"Reads assigned:   {n_assigned:,}" if isinstance(n_assigned, int) else "",
            f"Taxa with reads:  {n_taxa:,}" if isinstance(n_taxa, int) else "",
            f"Mean damage:      {mean_dmg:.4f}" if isinstance(mean_dmg, float) else "",
            f"SQLite viewer:    {result.get('sqlite_path') or 'N/A'}",
        ])

    return result


def run_batch_metagenome_pipeline(
    samples: List[Dict[str, Any]],
    *,
    shared_db_outdir: str | Path,
    kraken_db: Optional[str | Path] = None,
    screener: str = "auto",
    min_kraken_reads: int = 5,
    min_kraken_unique_kmers: int = 50,
    max_kraken_dup_rate: Optional[float] = None,
    nodes_path: str | Path,
    names_path: Optional[str | Path] = None,
    blast_db_dir: Optional[str | Path] = None,
    custom_ref_fastas: Optional[List[str | Path]] = None,
    max_seed_rank: str = "family",
    max_seqs_per_taxon: Optional[int] = None,
    chunked_db: bool = True,
    chunked_db_rank: Optional[str] = None,
    parallel_chunks: int = 1,
    chunk_timeout_sec: int = 1800,
    keep_chunk_mafs: bool = True,
    last_train_file: Optional[str | Path] = None,
    last_min_score: Optional[int] = None,
    last_m: Optional[int] = None,
    evalue: str = "1e-5",
    max_target_seqs: int = 2000,
    max_hsps: int = 1,
    chunk_size: int = 50000,
    resume: bool = True,
    threads: int = 16,
    last_threads: Optional[int] = None,
    taxonomy_tsv: Optional[str] = None,
    root_taxid: str = "1",
    sample_sheet: Optional[str] = None,
    regional_taxa: Optional[str] = None,
    config: Optional[str] = None,
    em_iterations: int = 1,
    dirichlet_alpha: float = 0.0,
    damage_weight: float = 0.0,
    no_deconvolve: bool = False,
    use_mini_db: bool = True,
    max_clade_seqs_per_taxon: int = 5000,
    max_clade_db_mb: float = 50.0,
    mini_db_timeout_min: float = 30.0,
    max_parallel_mini_dbs: int = 1,
    curate_refs: bool = True,
    min_ref_length: int = 300,
    max_ref_n_fraction: float = 0.5,
    ref_scope: str = "all",
    ref_mode: str = "full",
    nt_index_path: Optional[str | Path] = None,
    custom_ref_taxid_map: Optional[str | Path] = None,
    dedup_similarity: float = 0.90,
    no_viewer: bool = False,
    sqlite_viewer: bool = True,
    max_viewer_reads: int = 5000,
    lca_overrides: Optional[Dict[str, Any]] = None,
    uniqueness_index: Optional[str | Path] = None,
    yes: bool = True,
    verbose: bool = True,
) -> None:
    """Run the metagenomics pipeline for multiple samples sharing one LAST DB per group.

    Each dict in ``samples`` must contain:
      ``sample_id``  — str
      ``outdir``     — str | Path (per-sample output directory)
      ``group_id``   — str (optional, default 'default'); samples with the same
                       group_id share one LAST DB build, paying the cost only once.

    Plus one of:
      ``query``      — path to a preprocessed FASTA (skips preprocessing)
      ``fastq_r1``   — raw R1 FASTQ (triggers preprocessing + screening)
      ``fastq_r2``   — R2 FASTQ for paired-end (optional)

    Pipeline phases
    ---------------
    1. Screen every sample with KrakenUniq/Kraken2 (or auto-discover existing reports).
    2. For each group: union all taxids, build **one** shared subset LAST DB.
       The first sample in the group is also aligned and classified during this phase.
    3. Align and classify every remaining sample against the group's shared DB.

    The shared DBs are written to ``shared_db_outdir/<group_id>/``.
    """
    shared_db_outdir = Path(shared_db_outdir)
    shared_db_outdir.mkdir(parents=True, exist_ok=True)

    def _sample_lca(s: Dict) -> Dict[str, Any]:
        """Merge global lca_overrides with any per-sample _lca_overrides from the sample sheet."""
        merged = dict(lca_overrides or {})
        merged.update(s.get("_lca_overrides") or {})
        return merged

    # ── Phase 1: Screen (or auto-discover) all samples ────────────────────────
    screened: Dict[str, Dict[str, Any]] = {}
    if verbose:
        print(f"\n[batch] Phase 1 — screening {len(samples)} sample(s)", flush=True)

    for s in samples:
        sid = s["sample_id"]
        outdir = Path(s["outdir"])

        # Check for an existing screening report to avoid re-running KrakenUniq
        report = outdir / "screening" / f"{sid}.krakenuniq.report.txt"
        if not report.exists():
            report = outdir / "screening" / f"{sid}.kraken2.report.txt"

        auto_query = s.get("query") or str(outdir / "preprocessing" / f"{sid}.dedup.fasta")
        if report.exists() and Path(auto_query).exists():
            from .krakenutils import parse_krakenuniq_report
            raw_taxids = parse_krakenuniq_report(
                report,
                min_reads=min_kraken_reads,
                min_unique_kmers=min_kraken_unique_kmers,
                max_dup_rate=max_kraken_dup_rate,
            )
            if verbose:
                print(f"  [{sid}] reusing existing screening report: {len(raw_taxids):,} taxids",
                      flush=True)
            screened[sid] = {"query": auto_query, "raw_taxids": raw_taxids}
        else:
            result = run_screening_phase(
                query=s.get("query"),
                fastq_r1=s.get("fastq_r1"),
                fastq_r2=s.get("fastq_r2"),
                sample_id=sid,
                outdir=outdir,
                kraken_db=kraken_db,
                screener=screener,
                min_kraken_reads=min_kraken_reads,
                min_kraken_unique_kmers=min_kraken_unique_kmers,
                max_kraken_dup_rate=max_kraken_dup_rate,
                nodes_path=nodes_path,
                threads=threads,
                verbose=verbose,
            )
            screened[sid] = result

    # ── Phase 2: Build one shared DB per group ────────────────────────────────
    groups: Dict[str, List[Dict]] = {}
    for s in samples:
        gid = s.get("group_id", "default")
        groups.setdefault(gid, []).append(s)

    group_db_prefixes: Dict[str, Path] = {}
    done_samples: Set[str] = set()

    if verbose:
        print(f"\n[batch] Phase 2 — building {len(groups)} shared DB(s)", flush=True)

    for gid, group_samples in groups.items():
        group_db_dir = shared_db_outdir / gid
        group_db_dir.mkdir(parents=True, exist_ok=True)
        group_db_prefix = group_db_dir / gid
        group_db_prefixes[gid] = group_db_prefix

        if _db_exists(group_db_prefix):
            if verbose:
                print(f"  [group:{gid}] shared DB exists — skipping build", flush=True)
            continue

        # Union taxids from all samples in this group
        all_taxids: Set[str] = set()
        for s in group_samples:
            all_taxids |= screened[s["sample_id"]]["raw_taxids"]

        if verbose:
            n_sids = [s["sample_id"] for s in group_samples]
            print(
                f"  [group:{gid}] unioning {len(all_taxids):,} taxids "
                f"from {len(group_samples)} samples: {', '.join(n_sids)}",
                flush=True,
            )

        taxid_file = group_db_dir / f"{gid}.taxids.txt"
        taxid_file.write_text(
            "\n".join(sorted(all_taxids)) + "\n", encoding="utf-8"
        )

        # Run full pipeline for the first sample using the unioned taxid file.
        # from_stage='build-db' + taxid_file → skip screen → build DB + align + classify.
        # The DB is placed at group_db_prefix (shared location).
        first = group_samples[0]
        fsid = first["sample_id"]
        run_metagenome_pipeline(
            query=screened[fsid]["query"],
            sample_id=fsid,
            outdir=first["outdir"],
            from_stage="build-db",
            taxid_file=str(taxid_file),
            last_db_prefix=group_db_prefix,
            nodes_path=nodes_path,
            names_path=names_path,
            blast_db_dir=blast_db_dir,
            custom_ref_fastas=custom_ref_fastas or [],
            max_seed_rank=max_seed_rank,
            max_seqs_per_taxon=max_seqs_per_taxon,
            chunked_db=chunked_db,
            chunked_db_rank=chunked_db_rank,
            parallel_chunks=parallel_chunks,
            chunk_timeout_sec=chunk_timeout_sec,
            keep_chunk_mafs=keep_chunk_mafs,
            last_train_file=last_train_file,
            last_min_score=last_min_score,
            last_m=last_m,
            evalue=evalue,
            max_target_seqs=max_target_seqs,
            max_hsps=max_hsps,
            chunk_size=chunk_size,
            resume=resume,
            threads=threads,
            last_threads=last_threads,
            taxonomy_tsv=taxonomy_tsv,
            root_taxid=root_taxid,
            sample_sheet=sample_sheet,
            regional_taxa=regional_taxa,
            config=config,
            em_iterations=em_iterations,
            dirichlet_alpha=dirichlet_alpha,
            damage_weight=damage_weight,
            no_deconvolve=no_deconvolve,
            use_mini_db=use_mini_db,
            max_clade_seqs_per_taxon=max_clade_seqs_per_taxon,
            max_clade_db_mb=max_clade_db_mb,
            mini_db_timeout_min=mini_db_timeout_min,
            max_parallel_mini_dbs=max_parallel_mini_dbs,
            curate_refs=curate_refs,
            min_ref_length=min_ref_length,
            max_ref_n_fraction=max_ref_n_fraction,
            ref_scope=ref_scope,
            ref_mode=ref_mode,
            nt_index_path=nt_index_path,
            custom_ref_taxid_map=custom_ref_taxid_map,
            dedup_similarity=dedup_similarity,
            no_viewer=no_viewer,
            sqlite_viewer=sqlite_viewer,
            max_viewer_reads=max_viewer_reads,
            lca_overrides=_sample_lca(first),
            uniqueness_index=uniqueness_index,
            yes=yes,
            verbose=verbose,
        )
        done_samples.add(fsid)

    # ── Phase 3: Align + classify remaining samples ────────────────────────────
    remaining = [s for s in samples if s["sample_id"] not in done_samples]
    if verbose and remaining:
        print(
            f"\n[batch] Phase 3 — aligning/classifying {len(remaining)} remaining sample(s)",
            flush=True,
        )

    for s in remaining:
        sid = s["sample_id"]
        gid = s.get("group_id", "default")
        db_prefix = group_db_prefixes[gid]
        run_metagenome_pipeline(
            query=screened[sid]["query"],
            sample_id=sid,
            outdir=s["outdir"],
            from_stage="align",
            last_db_prefix=db_prefix,
            nodes_path=nodes_path,
            names_path=names_path,
            blast_db_dir=None,
            last_train_file=last_train_file,
            last_min_score=last_min_score,
            last_m=last_m,
            evalue=evalue,
            max_target_seqs=max_target_seqs,
            max_hsps=max_hsps,
            chunk_size=chunk_size,
            resume=resume,
            threads=threads,
            last_threads=last_threads,
            taxonomy_tsv=taxonomy_tsv,
            root_taxid=root_taxid,
            sample_sheet=sample_sheet,
            regional_taxa=regional_taxa,
            config=config,
            em_iterations=em_iterations,
            dirichlet_alpha=dirichlet_alpha,
            damage_weight=damage_weight,
            no_deconvolve=no_deconvolve,
            use_mini_db=use_mini_db,
            max_clade_seqs_per_taxon=max_clade_seqs_per_taxon,
            max_clade_db_mb=max_clade_db_mb,
            mini_db_timeout_min=mini_db_timeout_min,
            max_parallel_mini_dbs=max_parallel_mini_dbs,
            curate_refs=curate_refs,
            min_ref_length=min_ref_length,
            max_ref_n_fraction=max_ref_n_fraction,
            ref_scope=ref_scope,
            ref_mode=ref_mode,
            nt_index_path=nt_index_path,
            custom_ref_taxid_map=custom_ref_taxid_map,
            dedup_similarity=dedup_similarity,
            no_viewer=no_viewer,
            sqlite_viewer=sqlite_viewer,
            max_viewer_reads=max_viewer_reads,
            lca_overrides=_sample_lca(s),
            uniqueness_index=uniqueness_index,
            yes=yes,
            verbose=verbose,
        )
