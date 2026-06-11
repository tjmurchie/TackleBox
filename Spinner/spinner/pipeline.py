"""Main pipeline orchestration: run_pipeline().

This is the central function that ties together annotation, external-tool
hooks, scoring, capping, and output writing.  It is called by cli.py for
both ``audit`` and ``filter`` subcommands (and their screen-* variants).
"""
from __future__ import annotations

import datetime
import os
import sys
import time
from collections import Counter
from pathlib import Path
from typing import List, Optional

from . import VERSION
from .annotation import annotate, rescue_cross_species_duplicates, trim_and_rescreen_adapters

# Bundled config files shipped with Spinner — used as fallback defaults
# when the caller does not supply explicit paths.
_CONFIGS_DIR = Path(__file__).resolve().parent.parent / "configs"
_BUNDLED_DEFAULTS = {
    "regions_config": "regions_config_example.tsv",
    "adapters_tsv":   "adapters_default.tsv",
    "bad_keywords_tsv": "bad_keywords_ancient.tsv",
}
_MODE_CONFIGS = {
    "reference_db": "spinner_reference_db.yml",
    "bait_panel":   "spinner_bait_panel.yml",
}
from .capping import cap_refs, rescue_sole_representatives
from .clustering import run_uchime, run_vsearch
from .decisions import score_decide, write_decisions
from .external import parse_generic_hook_table, run_blast, run_mmseqs, run_template_command
from .fasta import parse_fasta, write_keyed_fasta
from .reporting import write_split_fastas, write_summary_html, write_summary_tsv
from .taxonomy_blast import (
    load_taxdb,
    make_windowed_fasta,
    parse_tax_blast,
    parse_tax_blast_escalation,
    parse_windowed_blast,
)
from .utils import BlastTicker, fmt_seconds, info, section, stage, warn
from .vector_screen import parse_vector_blast


def _write_command_txt(path: str, argv: List[str]) -> None:
    with open(path, "w", encoding="utf-8") as f:
        f.write("# Spinner command\n")
        f.write(" ".join(argv) + "\n")


def _write_resolved_config(cfg: dict, path: str) -> None:
    try:
        import yaml
        with open(path, "w", encoding="utf-8") as f:
            yaml.dump(cfg, f, default_flow_style=False, allow_unicode=True)
    except ImportError:
        import json
        with open(path, "w", encoding="utf-8") as f:
            json.dump(cfg, f, indent=2)


def _cleanup_temp(paths: list) -> None:
    for p in paths:
        if p and os.path.exists(p):
            try:
                os.remove(p)
            except OSError:
                pass


def _save_checkpoint(ann: dict, record_count: int, path: str) -> None:
    """Pickle the post-annotation state so a killed run can resume."""
    import pickle
    try:
        with open(path, "wb") as f:
            pickle.dump({"ann": ann, "record_count": record_count}, f, protocol=4)
        info(f"  Checkpoint saved → {os.path.basename(path)}")
    except Exception as e:
        warn(f"  Checkpoint save failed (non-fatal): {e}")


def _load_checkpoint(path: str, record_count: int) -> Optional[dict]:
    """Return pickled ann dict if the checkpoint matches the current input, else None."""
    import pickle
    if not os.path.exists(path):
        return None
    try:
        with open(path, "rb") as f:
            data = pickle.load(f)
        stored = data.get("record_count")
        if stored != record_count:
            warn(f"  [resume] Checkpoint record count ({stored}) != current ({record_count})"
                 " — ignoring stale checkpoint and starting fresh")
            return None
        return data["ann"]
    except Exception as e:
        warn(f"  [resume] Could not load checkpoint ({e}) — starting fresh")


def _print_step_banner(label: str, enabled: bool, reason: str = "") -> None:
    status = "ENABLED" if enabled else "skipped"
    note = f"  ({reason})" if reason else ""
    info(f"  {label:<30} {status}{note}")


def _print_final_summary(ann: dict, outprefix: str, filter_mode: bool,
                          elapsed: float, steps: dict) -> None:
    """Print a detailed end-of-run summary to stderr."""
    counts = Counter(a.decision for a in ann.values())
    n_total = len(ann)
    keep   = counts.get("KEEP", 0)
    review = counts.get("REVIEW", 0)
    reject = counts.get("REJECT", 0)

    stage(f"Run complete  [{fmt_seconds(elapsed)} total]")

    # Decision counts with percentages
    pct = lambda n: f"{100 * n / n_total:.1f}%" if n_total else "0%"
    info(f"  {'Records processed:':<28} {n_total:,}")
    info(f"  {'KEEP:':<28} {keep:>7,}  ({pct(keep)})")
    info(f"  {'REVIEW:':<28} {review:>7,}  ({pct(review)})")
    info(f"  {'REJECT:':<28} {reject:>7,}  ({pct(reject)})")

    # Class × decision breakdown
    by_class: Counter = Counter()
    for a in ann.values():
        by_class[(a.marker_class, a.decision)] += 1
    classes = sorted({k for k, _ in by_class})
    if classes:
        section("Decisions by marker class")
        header = f"  {'Class':<14}" + "".join(f"  {d:<8}" for d in ("KEEP", "REVIEW", "REJECT"))
        info(header)
        for cls in classes:
            row = (f"  {cls:<14}" +
                   "".join(f"  {by_class.get((cls, d), 0):<8}" for d in ("KEEP", "REVIEW", "REJECT")))
            info(row)

    # Kingdom breakdown
    by_kingdom: Counter = Counter()
    for a in ann.values():
        by_kingdom[(a.kingdom, a.decision)] += 1
    kingdoms = sorted({k for k, _ in by_kingdom})
    if len(kingdoms) > 1:
        section("Decisions by kingdom")
        header = f"  {'Kingdom':<14}" + "".join(f"  {d:<8}" for d in ("KEEP", "REVIEW", "REJECT"))
        info(header)
        for kd in kingdoms:
            row = (f"  {kd:<14}" +
                   "".join(f"  {by_kingdom.get((kd, d), 0):<8}" for d in ("KEEP", "REVIEW", "REJECT")))
            info(row)

    # Top reasons
    reasons: Counter = Counter(r for a in ann.values() for r in a.reasons)
    if reasons:
        section("Top rejection / review reasons")
        for rsn, n in reasons.most_common(10):
            info(f"  {n:>7,}  {rsn}")

    # Adapter summary
    adapter_hits = Counter(a.adapter_name for a in ann.values() if a.adapter_hit)
    if adapter_hits:
        section("Adapter contamination detected")
        for name, n in adapter_hits.most_common():
            info(f"  {n:>7,}  {name}")

    # Taxonomy summary (if run)
    tax_counts = Counter(a.taxonomy_status for a in ann.values())
    if any(st != "NOT_CHECKED" for st in tax_counts):
        section("Taxonomy BLAST status")
        for st, n in tax_counts.most_common():
            info(f"  {n:>7,}  {st}")

    # Windowed BLAST summary (if run)
    wind_counts = Counter(a.windowed_status for a in ann.values())
    if any(st != "NOT_CHECKED" for st in wind_counts):
        section("Windowed BLAST status")
        for st, n in wind_counts.most_common():
            info(f"  {n:>7,}  {st}")

    # Chimera screen summary (if run)
    chimera_count = sum(1 for a in ann.values() if "chimera_detected" in a.reasons)
    borderline_count = sum(1 for a in ann.values() if "chimera_borderline" in a.reasons)
    if chimera_count or borderline_count:
        section("Chimera screen results")
        info(f"  {chimera_count:>7,}  chimera_detected")
        info(f"  {borderline_count:>7,}  chimera_borderline")

    # Sole representative rescues (if enabled)
    rescued_count = sum(1 for a in ann.values() if "sole_representative" in a.reasons)
    if rescued_count:
        section("Sole representative rescues")
        info(f"  {rescued_count:>7,}  records promoted to KEEP as sole species representative")

    # Cross-species duplicate rescues
    dup_rescued = sum(1 for a in ann.values() if "rescued_duplicate" in a.reasons)
    if dup_rescued:
        section("Cross-species duplicate rescues")
        info(f"  {dup_rescued:>7,}  records promoted as species representative (shared haplotype)")

    # Adapter trim rescues
    trim_rescued = sum(1 for a in ann.values() if "adapter_trimmed_rescued" in a.reasons)
    if trim_rescued:
        section("Adapter trim rescues")
        info(f"  {trim_rescued:>7,}  records trimmed at adapter position and re-screened clean")

    # Cross-kingdom escalation rescues
    nr_rescued = sum(1 for a in ann.values() if "taxonomy_rescued_nr_protein" in a.reasons)
    nt_rescued = sum(1 for a in ann.values() if "taxonomy_rescued_nt_blast" in a.reasons)
    if nr_rescued or nt_rescued:
        section("Cross-kingdom escalation rescues")
        if nr_rescued:
            info(f"  {nr_rescued:>7,}  rescued by NR protein (Swiss-Prot coverage gap)")
        if nt_rescued:
            info(f"  {nt_rescued:>7,}  rescued by NT nucleotide database")

    # Output files
    section("Output files")
    outputs = [
        (outprefix + ".decisions.tsv",           "Audit table — one row per input record"),
        (outprefix + ".summary.tsv",              "Summary statistics (machine-readable)"),
        (outprefix + ".summary.html",             "Summary report (open in browser)"),
        (outprefix + ".command.txt",              "Exact command for reproducibility"),
        (outprefix + ".run_config.resolved.yml",  "Full resolved config snapshot"),
    ]
    if filter_mode:
        outputs += [
            (outprefix + ".keep.fasta",   f"High-confidence records ({keep:,})"),
            (outprefix + ".review.fasta", f"Records to inspect manually ({review:,})"),
            (outprefix + ".reject.fasta", f"Rejected records — kept for audit ({reject:,})"),
        ]
    for path, desc in outputs:
        exists = "+" if os.path.exists(path) else "!"
        info(f"  [{exists}] {path}")
        info(f"       {desc}")

    # Next steps hints
    section("What to do next")
    if filter_mode:
        info("  1. Open summary.html in a browser to review counts, classes, and top reasons.")
        info("  2. Inspect review.fasta — these records may be salvageable after manual check.")
        info("  3. Investigate any unexpected rejections with the explain subcommand:")
        info(f"       ./Spinner explain --decisions {outprefix}.decisions.tsv --accession ACCESSION")
        info("  4. Regenerate the HTML report after manual TSV edits:")
        info(f"       ./Spinner report --decisions {outprefix}.decisions.tsv --outprefix <prefix>")
        if not steps.get("taxonomy_blast"):
            info("  5. Add taxonomy sanity checking by re-running with spinner_with_nt_blast.yml")
            info("     (requires blastn + NCBI nt database)")
        if not steps.get("cluster"):
            info("  6. Reduce haplotype redundancy by enabling cluster: true in your config")
            info("     (requires vsearch)")
    else:
        info("  1. Review decisions.tsv — no FASTAs were written (audit mode).")
        info("  2. Re-run with the 'filter' subcommand to produce keep/review/reject FASTAs:")
        info(f"       ./Spinner filter [same args] --outprefix <prefix>")


def run_pipeline(args, filter_mode: bool) -> None:
    """Run the full Spinner annotation and decision pipeline.

    Parameters
    ----------
    args:
        Parsed argparse.Namespace with at minimum fasta, outprefix, config.
    filter_mode:
        When True, write keep / review / reject FASTA files.
    """
    from .config import load_config

    # Resolve --mode to bundled config if --config not explicitly provided
    mode = getattr(args, "mode", "")
    if mode and not getattr(args, "config", ""):
        mode_file = _CONFIGS_DIR / _MODE_CONFIGS[mode]
        if mode_file.exists():
            args.config = str(mode_file)
        else:
            warn(f"Bundled config for --mode '{mode}' not found: {mode_file}; using defaults")

    run_start = time.time()
    cfg = load_config(getattr(args, "config", "") or "")

    # --- CLI overrides for dedicated screen-* subcommands ---
    if getattr(args, "taxonomy_blast_db", ""):
        cfg["steps"]["taxonomy_blast"] = True
        cfg["taxonomy_blast"]["blast_db"] = args.taxonomy_blast_db
    if getattr(args, "vector_blast_db", ""):
        cfg["steps"]["vector_screen"] = True
        cfg["vector_screen"]["blast_db"] = args.vector_blast_db
    if getattr(args, "windowed_blast_db", ""):
        cfg["steps"]["windowed_blast"] = True
        cfg["windowed_blast"]["blast_db"] = args.windowed_blast_db
    if getattr(args, "enable_cluster", False):
        cfg["steps"]["cluster"] = True

    # --- map CLI input-file args into config ---
    if getattr(args, "regions_config", ""):
        cfg["inputs"]["regions_config"] = args.regions_config
    if getattr(args, "species_kingdom", ""):
        cfg["inputs"]["species_kingdom"] = args.species_kingdom
    if getattr(args, "adapters", ""):
        cfg["inputs"]["adapters_tsv"] = args.adapters
    if getattr(args, "bad_keywords", ""):
        cfg["inputs"]["bad_keywords_tsv"] = args.bad_keywords

    # --keep-temp override
    if getattr(args, "keep_temp", False):
        cfg["run"]["keep_temp_files"] = True

    # --- fall back to bundled config files when not explicitly provided ---
    for key, filename in _BUNDLED_DEFAULTS.items():
        if not cfg["inputs"].get(key):
            candidate = _CONFIGS_DIR / filename
            if candidate.exists():
                cfg["inputs"][key] = str(candidate)

    steps = cfg.get("steps", {})

    # --- Build ordered step plan for "Step X/N" display ---
    _step_plan = ["Load FASTA", "Annotate and QC screen"]
    if steps.get("vector_screen", False):  _step_plan.append("Vector screen")
    if steps.get("taxonomy_blast", False): _step_plan.append("Taxonomy BLAST")
    if steps.get("windowed_blast", False): _step_plan.append("Windowed BLAST — chimerism screen")
    if steps.get("chimera_screen", False): _step_plan.append("Chimera screen — vsearch uchime")
    if steps.get("fcs_adaptor", False):    _step_plan.append("FCS-adaptor screen")
    if steps.get("fcs_gx", False):         _step_plan.append("FCS-GX contamination screen")
    if steps.get("cluster", False):        _step_plan.append("Clustering")
    _step_plan.append("Score, cap and write outputs")
    _N = len(_step_plan)
    _s = [0]  # mutable step counter

    def _stage(label: str) -> None:
        _s[0] += 1
        stage(f"[{_s[0]}/{_N}]  {label}")

    # --- reproducibility outputs ---
    outprefix = args.outprefix
    Path(outprefix).parent.mkdir(parents=True, exist_ok=True)
    _write_command_txt(outprefix + ".command.txt", sys.argv)
    _write_resolved_config(cfg, outprefix + ".run_config.resolved.yml")

    temp_files: list = []

    # --- Print run header ---
    stage(f"Spinner {VERSION}  |  {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    mode_label = ("filter (keep/review/reject FASTAs will be written)"
                  if filter_mode else "audit (decisions only, no FASTA output)")
    info(f"  Mode:    {mode_label}")
    info(f"  Config:  {getattr(args, 'config', 'defaults') or 'defaults'}")
    info(f"  Output:  {outprefix}.*")

    section("Steps enabled for this run")
    _print_step_banner("basic_qc",           steps.get("basic_qc", True))
    _print_step_banner("classify_regions",   steps.get("classify_regions", True))
    _print_step_banner("adapter_screen",     steps.get("adapter_screen", True))
    _print_step_banner("bad_keyword_screen", steps.get("bad_keyword_screen", True))
    _print_step_banner("vector_screen",      steps.get("vector_screen", False))
    _print_step_banner("taxonomy_blast",     steps.get("taxonomy_blast", False))
    _print_step_banner("windowed_blast",     steps.get("windowed_blast", False))
    _print_step_banner("chimera_screen",     steps.get("chimera_screen", False))
    _print_step_banner("cluster",            steps.get("cluster", False))
    _print_step_banner("cap_references",     steps.get("cap_references", True))
    rescue = cfg.get("capping", {}).get("rescue_sole_representatives", False)
    _print_step_banner("rescue_sole_reps",   rescue)

    # Checkpoint file — saved after annotation so a killed run can resume.
    checkpoint_file = outprefix + ".checkpoint.pkl"
    temp_files.append(checkpoint_file)

    # -----------------------------------------------------------------------
    _stage("Load FASTA")
    t0 = time.time()
    records = parse_fasta(args.fasta)
    info(f"Loaded {len(records):,} records from {len(args.fasta)} file(s)"
         f"  [{fmt_seconds(time.time() - t0)}]")

    # -----------------------------------------------------------------------
    _stage("Annotate and QC screen")
    _ann_from_ckpt = _load_checkpoint(checkpoint_file, len(records))
    if _ann_from_ckpt is not None:
        ann = _ann_from_ckpt
        info(f"  [resume] Loaded annotation checkpoint — skipping QC scan ({len(ann):,} records)")
    else:
        section("Per-record QC, classification, adapter and keyword scanning")
        t0 = time.time()
        ann = annotate(
            records,
            cfg,
            cfg["inputs"].get("species_kingdom"),
            cfg["inputs"].get("regions_config"),
            cfg["inputs"].get("adapters_tsv"),
            cfg["inputs"].get("bad_keywords_tsv"),
        )
        info(f"Annotation complete  [{fmt_seconds(time.time() - t0)}]")

        # Inline annotation stats (pre-rescue / pre-trim counts)
        dup_seq  = sum(1 for a in ann.values() if a.duplicate_sequence)
        dup_acc  = sum(1 for a in ann.values() if a.duplicate_accession)
        adapt    = sum(1 for a in ann.values() if a.adapter_hit)
        kw_hit   = sum(1 for a in ann.values() if a.bad_keyword_hit)
        hi_n     = sum(1 for a in ann.values() if "n_fraction_high" in a.reasons)
        lo_cplx  = sum(1 for a in ann.values() if "low_complexity" in a.reasons)
        info(f"  Duplicate sequences:   {dup_seq:,}")
        info(f"  Duplicate accessions:  {dup_acc:,}")
        info(f"  Adapter hits:          {adapt:,}")
        info(f"  Bad keyword hits:      {kw_hit:,}")
        info(f"  High N-fraction:       {hi_n:,}")
        info(f"  Low complexity:        {lo_cplx:,}")

        _save_checkpoint(ann, len(records), checkpoint_file)

    # Cross-species duplicate rescue: species whose sequences are all exact copies
    # of another species' sequences get one representative un-flagged so that the
    # species is not left with zero KEEP candidates.  Runs on both fresh and resumed
    # runs (idempotent — already-rescued records are skipped).
    rescue_n = rescue_cross_species_duplicates(ann, cfg)
    if rescue_n:
        info(f"  Cross-species duplicate rescue: {rescue_n:,} records promoted")

    # Adapter trim-and-rescreen: internal adapter hits are trimmed at the detected
    # position and the clean flank is re-checked.  Rescued records have their
    # sequences updated in-place so external tools receive the trimmed version.
    trim_n = trim_and_rescreen_adapters(
        records, ann, cfg, cfg["inputs"].get("adapters_tsv", "")
    )
    if trim_n:
        info(f"  Adapter trim-and-rescreen: {trim_n:,} records trimmed and rescued")

    # Exclude duplicate records from external screening — they produce identical
    # BLAST/chimera results as their primary copy and skipping them saves
    # significant runtime, especially for windowed BLAST.
    dup_keys = {k for k, a in ann.items() if a.duplicate_accession or a.duplicate_sequence}
    screened_records = [r for r in records if r.id not in dup_keys]
    if dup_keys:
        info(f"  {len(dup_keys):,} duplicate records excluded from external screens")

    # Write keyed FASTA for external tool input (duplicates excluded).
    # Always regenerated — fast, and needed by downstream external steps.
    tmp_keyed = outprefix + ".tmp.keyed.fasta"
    write_keyed_fasta(screened_records, tmp_keyed)
    temp_files.append(tmp_keyed)

    # -----------------------------------------------------------------------
    if steps.get("vector_screen", False):
        _stage("Vector screen via BLAST / UniVec")
        db = cfg["vector_screen"].get("blast_db", "")
        if db:
            info(f"  Database: {db}")
            out_v = outprefix + ".vector_blast.tsv"
            t0 = time.time()
            try:
                with BlastTicker("Vector BLAST",
                                 output_file=out_v,
                                 total_queries=len(screened_records),
                                 avg_hsp=3.0):
                    run_blast(tmp_keyed, db, out_v, cfg["vector_screen"])
                parse_vector_blast(out_v, ann, cfg)
                vec_hits = sum(1 for a in ann.values() if a.vector_hit)
                info(f"  Vector hits found: {vec_hits:,}  [{fmt_seconds(time.time() - t0)}]")
            except Exception as e:
                if cfg["run"].get("fail_on_missing_external_tool", False):
                    raise
                warn(f"Vector screen skipped/failed: {e}")
        else:
            warn("vector_screen enabled but vector_screen.blast_db not set — skipping")

    # -----------------------------------------------------------------------
    taxdb = None
    if steps.get("taxonomy_blast", False):
        tax_method = cfg["taxonomy_blast"].get("method", "blastn").lower()
        _stage(f"Taxonomy BLAST  [{tax_method}]")
        db = cfg["taxonomy_blast"].get("blast_db", "")
        taxdump_dir = cfg["taxonomy_blast"].get("taxdump_dir", "")
        if db:
            info(f"  Database: {db}")
            if taxdump_dir:
                info(f"  Taxdump dir:    {taxdump_dir}")
                info("  Loading NCBI taxdump (~30 s, ~2 GB RAM) ...")
            else:
                info("  Taxdump not configured — using string-matching mode only")
            t0 = time.time()
            taxdb = load_taxdb(cfg)
            if taxdump_dir and taxdb:
                info(f"  Taxdump loaded  [{fmt_seconds(time.time() - t0)}]")

            n_threads = cfg["taxonomy_blast"].get("num_threads", 1)
            max_qlen = int(cfg["taxonomy_blast"].get("max_query_length", 0))
            if max_qlen > 0:
                tax_records = [r for r in screened_records if len(r.seq_upper) <= max_qlen]
                n_skipped = len(screened_records) - len(tax_records)
                if n_skipped:
                    info(f"  Skipping {n_skipped:,} sequences >{max_qlen:,} bp"
                         f" (max_query_length) — will be marked taxonomy_not_checked")
                tmp_keyed_tax = outprefix + ".tmp.keyed.tax.fasta"
                write_keyed_fasta(tax_records, tmp_keyed_tax)
                temp_files.append(tmp_keyed_tax)
                tax_query_fasta = tmp_keyed_tax
            else:
                tax_records = screened_records
                tax_query_fasta = tmp_keyed

            out_t = outprefix + ".taxonomy_blast.tsv"
            ckpt_dir = out_t + ".mmseqs_batches"
            already_done = (
                os.path.exists(out_t)
                and os.path.getsize(out_t) > 0
                and not os.path.exists(ckpt_dir)
            )
            if already_done:
                info(f"  [resume] {os.path.basename(out_t)} already exists — reusing cached results")
            else:
                info(f"  Running {tax_method} against {len(tax_records):,} records"
                     f"  ({len(dup_keys):,} duplicates pre-excluded)  threads={n_threads}")
                t0 = time.time()
                try:
                    ticker_label = (f"Taxonomy {tax_method}  ({len(screened_records):,} queries,"
                                    f" {n_threads} threads)")
                    batch_info: list = [0, 0]
                    with BlastTicker(ticker_label,
                                     output_file=out_t,
                                     total_queries=len(tax_records),
                                     avg_hsp=5.0) as ticker:
                        if tax_method == "mmseqs2":
                            ticker.batch_info = batch_info
                            run_mmseqs(tax_query_fasta, db, out_t, cfg["taxonomy_blast"],
                                       batch_info=batch_info)
                        else:
                            run_blast(tax_query_fasta, db, out_t, cfg["taxonomy_blast"])
                    info(f"  Search complete  [{fmt_seconds(time.time() - t0)}]")
                except Exception as e:
                    if cfg["run"].get("fail_on_missing_external_tool", False):
                        raise
                    warn(f"Taxonomy search skipped/failed: {e}")
            if os.path.exists(out_t):
                info("  Parsing results and updating annotations ...")
                parse_tax_blast(out_t, ann, cfg, taxdb)
                tax_counts = Counter(a.taxonomy_status for a in ann.values())
                for status, n in tax_counts.most_common():
                    info(f"    {n:>7,}  {status}")

                # --- Multi-database cross-kingdom escalation ---
                # Sequences rejected as cross-kingdom by Swiss-Prot are re-checked
                # against NR protein (level 1) and NT nucleotide (level 2).  Only
                # the small subset that failed the fast database reaches the slow one.
                if cfg["taxonomy_blast"].get("escalate_cross_kingdom", False):
                    cross_k_keys = {k for k, a in ann.items()
                                    if "taxonomy_cross_kingdom" in a.reasons}
                    if cross_k_keys:
                        info(f"  Cross-kingdom escalation: {len(cross_k_keys):,} records to re-verify")
                        esc_records = [r for r in screened_records if r.id in cross_k_keys]

                        # Level 1 — NR protein (MMSeqs2)
                        nr_db = cfg["taxonomy_blast"].get("nr_protein_db", "")
                        if nr_db and esc_records:
                            out_esc_nr = outprefix + ".escalation_nr.tsv"
                            esc_nr_fasta = outprefix + ".tmp.esc_nr.fasta"
                            temp_files.append(esc_nr_fasta)
                            ck_esc_nr = out_esc_nr + ".mmseqs_batches"
                            nr_already_done = (
                                os.path.exists(out_esc_nr)
                                and os.path.getsize(out_esc_nr) > 0
                                and not os.path.exists(ck_esc_nr)
                            )
                            if nr_already_done:
                                info("  [resume] escalation_nr.tsv already exists — reusing")
                            else:
                                write_keyed_fasta(esc_records, esc_nr_fasta)
                                nr_cfg = dict(cfg["taxonomy_blast"])
                                nr_cfg["blast_db"] = nr_db
                                if cfg["taxonomy_blast"].get("nr_mmseqs_binary"):
                                    nr_cfg["mmseqs_binary"] = cfg["taxonomy_blast"]["nr_mmseqs_binary"]
                                try:
                                    info(f"  Level 1 — NR protein ({len(esc_records):,} seqs): {nr_db}")
                                    esc_batch_info: list = [0, 0]
                                    with BlastTicker(
                                        f"Escalation NR  ({len(esc_records):,} seqs)",
                                        output_file=out_esc_nr,
                                        total_queries=len(esc_records),
                                        avg_hsp=5.0,
                                    ) as esc_ticker:
                                        esc_ticker.batch_info = esc_batch_info
                                        run_mmseqs(esc_nr_fasta, nr_db, out_esc_nr, nr_cfg,
                                                   batch_info=esc_batch_info)
                                except Exception as e:
                                    warn(f"Escalation NR search failed: {e}")
                            if os.path.exists(out_esc_nr):
                                nr_rescued = parse_tax_blast_escalation(
                                    out_esc_nr, ann, cfg, taxdb, "nr_protein"
                                )
                                info(f"    NR protein rescued: {nr_rescued:,} records")

                        # Level 2 — NT nucleotide BLAST (only remaining cross-kingdom)
                        nt_db = cfg["taxonomy_blast"].get("nt_blast_db", "")
                        still_cross_keys = {k for k, a in ann.items()
                                            if "taxonomy_cross_kingdom" in a.reasons}
                        nt_records = [r for r in screened_records
                                      if r.id in still_cross_keys]
                        if nt_db and nt_records:
                            out_esc_nt = outprefix + ".escalation_nt.tsv"
                            esc_nt_fasta = outprefix + ".tmp.esc_nt.fasta"
                            temp_files.append(esc_nt_fasta)
                            nt_already_done = (
                                os.path.exists(out_esc_nt)
                                and os.path.getsize(out_esc_nt) > 0
                            )
                            if nt_already_done:
                                info("  [resume] escalation_nt.tsv already exists — reusing")
                            else:
                                write_keyed_fasta(nt_records, esc_nt_fasta)
                                nt_cfg = {
                                    "blast_task": "megablast",
                                    "max_target_seqs": cfg["taxonomy_blast"].get("max_target_seqs", 10),
                                    "evalue": cfg["taxonomy_blast"].get("evalue", "1e-10"),
                                    "num_threads": cfg["taxonomy_blast"].get("num_threads", 1),
                                    "outfmt": "6 qseqid saccver pident length qlen qstart qend evalue bitscore staxids sscinames",
                                }
                                try:
                                    info(f"  Level 2 — NT megablast ({len(nt_records):,} seqs): {nt_db}")
                                    with BlastTicker(
                                        f"Escalation NT  ({len(nt_records):,} seqs)",
                                        output_file=out_esc_nt,
                                        total_queries=len(nt_records),
                                        avg_hsp=5.0,
                                    ):
                                        run_blast(esc_nt_fasta, nt_db, out_esc_nt, nt_cfg)
                                except Exception as e:
                                    warn(f"Escalation NT search failed: {e}")
                            if os.path.exists(out_esc_nt):
                                nt_rescued = parse_tax_blast_escalation(
                                    out_esc_nt, ann, cfg, taxdb, "nt_blast", nt_mode=True
                                )
                                info(f"    NT megablast rescued: {nt_rescued:,} records")

                        final_cross = sum(1 for a in ann.values()
                                          if "taxonomy_cross_kingdom" in a.reasons)
                        info(f"  Confirmed cross-kingdom after escalation: {final_cross:,}")
                    else:
                        info("  Cross-kingdom escalation: no rejections to re-verify")
        else:
            warn("taxonomy_blast enabled but taxonomy_blast.blast_db not set — skipping")

    # -----------------------------------------------------------------------
    if steps.get("windowed_blast", False):
        _stage("Windowed BLAST — chimerism screen")
        db = cfg["windowed_blast"].get("blast_db") or cfg["taxonomy_blast"].get("blast_db", "")
        if db:
            info(f"  BLAST database: {db}")
            wb_cfg = cfg.get("windowed_blast", {})
            info(f"  Window size: {wb_cfg.get('window_size', 500)} bp"
                 f"  step: {wb_cfg.get('window_step', 250)} bp"
                 f"  min sequence length: {wb_cfg.get('enabled_for_min_length', 1000)} bp")
            win_fasta = outprefix + ".windows.fasta"
            temp_files.append(win_fasta)
            nwin = make_windowed_fasta(screened_records, win_fasta, cfg)
            if nwin:
                n_threads = cfg["windowed_blast"].get("num_threads", 1)
                info(f"  {nwin:,} windows from {sum(1 for r in screened_records if len(r.seq_upper) >= wb_cfg.get('enabled_for_min_length', 1000)):,}"
                     f" long records  (duplicates pre-excluded)  threads={n_threads}")
                out_w = outprefix + ".windowed_blast.tsv"
                if os.path.exists(out_w) and os.path.getsize(out_w) > 0:
                    info(f"  [resume] {os.path.basename(out_w)} already exists — reusing cached results")
                else:
                    t0 = time.time()
                    try:
                        with BlastTicker(f"Windowed BLAST  ({nwin:,} windows, {n_threads} threads)",
                                         output_file=out_w,
                                         total_queries=nwin,
                                         avg_hsp=2.0):
                            run_blast(win_fasta, db, out_w, cfg["windowed_blast"])
                        info(f"  BLAST complete  [{fmt_seconds(time.time() - t0)}]")
                    except Exception as e:
                        if cfg["run"].get("fail_on_missing_external_tool", False):
                            raise
                        warn(f"Windowed BLAST skipped/failed: {e}")
                if os.path.exists(out_w):
                    parse_windowed_blast(out_w, ann, cfg, taxdb)
                    wind_counts = Counter(a.windowed_status for a in ann.values()
                                         if a.windowed_status != "NOT_CHECKED")
                    for status, n in wind_counts.most_common():
                        info(f"    {n:>7,}  {status}")
            else:
                info("  No records meet minimum length for windowed BLAST — skipping")
        else:
            warn("windowed_blast enabled but no blast_db configured — skipping")

    # -----------------------------------------------------------------------
    if steps.get("chimera_screen", False):
        _stage("Chimera screen — vsearch uchime")
        t0 = time.time()
        n_chim, n_border = run_uchime(tmp_keyed, outprefix, cfg, ann,
                                       total_queries=len(screened_records))
        info(f"  Chimeras detected: {n_chim:,}  |  Borderline: {n_border:,}"
             f"  [{fmt_seconds(time.time() - t0)}]")

    # -----------------------------------------------------------------------
    if steps.get("fcs_adaptor", False):
        _stage("FCS-adaptor screen")
        cmd = cfg.get("fcs_adaptor", {}).get("command", "")
        try:
            run_template_command(cmd, tmp_keyed, outprefix, "fcs_adaptor")
            table = cfg.get("fcs_adaptor", {}).get(
                "results_tsv", outprefix + ".fcs_adaptor.tsv"
            )
            parse_generic_hook_table(table, ann, "fcs_adaptor_hit", "fcs_adaptor_review")
        except Exception as e:
            if cfg["run"].get("fail_on_missing_external_tool", False):
                raise
            warn(f"FCS-adaptor hook skipped/failed: {e}")

    # -----------------------------------------------------------------------
    if steps.get("fcs_gx", False):
        _stage("FCS-GX contamination screen")
        cmd = cfg.get("fcs_gx", {}).get("command", "")
        try:
            run_template_command(cmd, tmp_keyed, outprefix, "fcs_gx")
            table = cfg.get("fcs_gx", {}).get(
                "results_tsv", outprefix + ".fcs_gx.tsv"
            )
            parse_generic_hook_table(table, ann, "fcs_gx_contaminant", "fcs_gx_review")
        except Exception as e:
            if cfg["run"].get("fail_on_missing_external_tool", False):
                raise
            warn(f"FCS-GX hook skipped/failed: {e}")

    # -----------------------------------------------------------------------
    if steps.get("cluster", False):
        _stage("Clustering")
        info(f"  Identity threshold: {cfg.get('cluster', {}).get('identity', 0.99)}")
        run_vsearch(tmp_keyed, outprefix, cfg, ann,
                    total_queries=len(screened_records))
        centroids = sum(1 for a in ann.values() if a.cluster_role == "centroid")
        members   = sum(1 for a in ann.values() if a.cluster_role == "member")
        info(f"  Cluster centroids: {centroids:,}  |  non-centroid members: {members:,}")

    # -----------------------------------------------------------------------
    _stage("Score, cap and write outputs")
    t0 = time.time()

    section("First-pass scoring")
    score_decide(ann, cfg)
    pre = Counter(a.decision for a in ann.values())
    info(f"  Before capping — KEEP={pre.get('KEEP', 0):,}  "
         f"REVIEW={pre.get('REVIEW', 0):,}  REJECT={pre.get('REJECT', 0):,}")

    section("Applying species x marker capping")
    cap_refs(ann, cfg)
    cap_hit = sum(1 for a in ann.values() if "cap_exceeded" in a.reasons)
    info(f"  Records exceeding cap: {cap_hit:,}")

    section("Final scoring")
    score_decide(ann, cfg)

    section("Sole representative rescue")
    rescued = rescue_sole_representatives(ann, cfg)
    if rescued:
        info(f"  Rescued {rescued:,} record(s) with reason 'sole_representative'"
             " (only KEEP for their species × class)")
    else:
        info("  Not enabled or no species required rescue"
             if not cfg.get("capping", {}).get("rescue_sole_representatives", False)
             else "  All identified species already have ≥1 KEEP record")

    section("Writing output files")
    write_decisions(ann, outprefix + ".decisions.tsv")
    info(f"  decisions.tsv  ({len(ann):,} rows)")
    write_summary_tsv(ann, outprefix)
    info("  summary.tsv")
    if steps.get("report", True):
        write_summary_html(ann, outprefix)
        info("  summary.html")
    if filter_mode:
        write_split_fastas(records, ann, outprefix, cfg)
        final = Counter(a.decision for a in ann.values())
        info(f"  keep.fasta   ({final.get('KEEP', 0):,} records)")
        info(f"  review.fasta ({final.get('REVIEW', 0):,} records)")
        info(f"  reject.fasta ({final.get('REJECT', 0):,} records)")
    info(f"  Scoring + writing complete  [{fmt_seconds(time.time() - t0)}]")

    # --- temp file cleanup ---
    if not cfg["run"].get("keep_temp_files", False):
        _cleanup_temp(temp_files)

    # --- detailed final summary ---
    _print_final_summary(ann, outprefix, filter_mode, time.time() - run_start, steps)
