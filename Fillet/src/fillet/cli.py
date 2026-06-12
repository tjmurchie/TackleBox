from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any, Dict, Generator, Iterable, List, Optional, Sequence, Tuple
import shutil
import subprocess
import sys
import time

from . import __version__
from .align import run_chunked_alignment, fillet_blast_fields
from .config import load_config
from .io import (
    group_hits_by_read,
    iter_blast_hits,
    iter_last_maf_hits,
    read_acc2taxid,
    read_contaminants,
    read_fasta_lengths,
    read_fasta_records,
    read_regional_taxa,
    read_sample_sheet,
    read_support_table,
    write_tsv,
)
from .relic_lca import RelicLCAClassifier, ReadAssignment, TaxonSummary, assignments_to_rows, summaries_to_rows, summarize_assignments, compute_read_damage
from .taxonomy import Taxonomy
from .viewer import build_viewer_payload, write_viewer_html, write_viewer_json
from .viewer_sqlite import write_viewer_sqlite, serve_viewer
from .report import write_metamerge_tables, write_megan_wide, write_holi_compat, write_read2tax, write_krona_like


QUICK_LCA_FLAGS = [
    ("--min-bitscore", "hit_filters", "min_bitscore", float),
    ("--max-evalue", "hit_filters", "max_evalue", str),
    ("--min-aln-len", "hit_filters", "min_aln_len", int),
    ("--min-pident", "hit_filters", "min_pident", float),
    ("--min-qcov", "hit_filters", "min_qcov", float),
    ("--min-read-len", "hit_filters", "min_read_len", int),
    ("--max-read-len", "hit_filters", "max_read_len", int),
    ("--top-bitscore-fraction", "hit_filters", "top_bitscore_fraction", float),
    ("--max-hits-per-read", "hit_filters", "max_hits_per_read", int),
    ("--max-hits-per-taxon", "hit_filters", "max_hits_per_taxon", int),
    ("--max-hits-per-genus", "hit_filters", "max_hits_per_genus", int),
    ("--min-node-posterior", "assignment", "min_node_posterior", float),
    ("--min-posterior-margin", "assignment", "min_posterior_margin", float),
    ("--min-posterior-family", "assignment", "min_posterior_family", float),
    ("--min-posterior-genus", "assignment", "min_posterior_genus", float),
    ("--min-posterior-species", "assignment", "min_posterior_species", float),
    ("--min-em-coherence-fraction", "assignment", "min_em_coherence_fraction", float),
    ("--em-coherence-max-lca-level", "assignment", "em_coherence_max_lca_level", int),
    ("--max-assignment-rank-level", "assignment", "max_assignment_rank_level", int),
    ("--min-cumulative-reads", "summary", "min_cumulative_reads_for_confident_call", int),
    ("--min-best-reference-breadth", "summary", "min_best_reference_breadth_for_confident_call", float),
    ("--max-top-locus-fraction", "summary", "max_top_locus_fraction", float),
    ("--max-stack-concentration-filter", "summary", "max_stack_concentration_filter", float),
    ("--min-mean-windows-per-ref", "summary", "min_mean_windows_per_ref", float),
    ("--max-per-ref-ses-filter", "summary", "max_per_ref_ses_filter", float),
    ("--ses-filter", "summary", "max_stack_excess_score_filter", float),
    ("--ses-combined-filter", "summary", "ses_combined_filter", float),
    ("--ses-combined-max-windows", "summary", "ses_combined_max_windows_per_ref", float),
    ("--min-reads-fungi", "summary", "min_reads_fungi", int),
    ("--min-composite-authenticity", "summary", "min_composite_authenticity", float),
    ("--composite-damage-weight-fungi", "summary", "composite_damage_weight_fungi", float),
    ("--incongruent-dominant-fraction", "summary", "incongruent_dominant_fraction", float),
    ("--damage-mode", "damage", "damage_mode", str),
]


def add_taxonomy_args(p: argparse.ArgumentParser, required: bool = True) -> None:
    taxgrp = p.add_mutually_exclusive_group(required=required)
    taxgrp.add_argument("--taxonomy-tsv", help="Small/curated taxonomy TSV with columns taxid,parent_taxid,rank,name.")
    taxgrp.add_argument("--nodes", help="NCBI nodes.dmp; requires --names.")
    p.add_argument("--names", help="NCBI names.dmp when using --nodes.")
    p.add_argument("--root-taxid", default="1")


def add_lca_override_args(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("quick RELIC-LCA overrides", "Override common config values from the command line. These are base/floor settings; RELIC-LCA still scores evidence continuously and can stop at a stable ancestor rather than forcing a brittle species call.")
    for flag, _section, key, caster in QUICK_LCA_FLAGS:
        arg_name = flag
        kwargs = {"dest": key, "default": None}
        if caster is int:
            kwargs["type"] = int
        elif caster is float:
            kwargs["type"] = float
        else:
            kwargs["type"] = str
        g.add_argument(arg_name, **kwargs)
    # Boolean flag: graduated cross-clade penalty (on by default)
    g.add_argument(
        "--no-cross-clade-graduated",
        dest="cross_clade_graduated_off",
        action="store_true",
        default=False,
        help=(
            "Disable graduated cross-clade penalty and revert to the original "
            "binary threshold + flat factor (cross_clade_hit_penalty)."
        ),
    )
    g.add_argument(
        "--suppress-unclassified-nodes",
        dest="suppress_unclassified_nodes",
        action="store_true",
        default=False,
        help=(
            "Remove NCBI taxonomy bucket nodes from output: 'X (in: Y)', "
            "'uncultured X', 'environmental samples', 'unclassified sequences'. "
            "These are never informative target taxa."
        ),
    )


def apply_config_overrides(cfg: Dict, args: argparse.Namespace) -> Dict:
    for _flag, section, key, _caster in QUICK_LCA_FLAGS:
        val = getattr(args, key, None)
        if val is not None:
            cfg.setdefault(section, {})[key] = val
    if getattr(args, "cross_clade_graduated_off", False):
        cfg.setdefault("weights", {})["cross_clade_graduated"] = False
    if getattr(args, "suppress_unclassified_nodes", False):
        cfg.setdefault("summary", {})["suppress_unclassified_nodes"] = True
    return cfg


def apply_mode_profile(cfg: Dict, mode: Optional[str]) -> Dict:
    """Merge per-mode default parameters into *cfg* (in-place).

    Profile values only fill keys that are not already set by the user's TOML
    config, so explicit TOML settings always take priority over mode defaults.
    CLI flags applied after this call (via apply_config_overrides) will override
    both profile and TOML values.
    """
    if not mode or mode == "auto":
        return cfg
    from .data.mode_profiles import MODE_PROFILES
    profile = MODE_PROFILES.get(mode, {})
    for section, values in profile.items():
        cfg.setdefault(section, {})
        for k, v in values.items():
            cfg[section].setdefault(k, v)
    return cfg


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="fillet",
        description="Fillet: independent RELIC-LCA classifier + MEGAN-like browser for TackleBox.",
    )
    p.add_argument("--version", action="version", version=f"Fillet {__version__}")
    sub = p.add_subparsers(dest="cmd")

    c = sub.add_parser("classify", help="Classify BLAST/LAST-like tabular alignments with RELIC-LCA.")
    c.add_argument("--input", nargs="+", required=True, help="One or more BLAST outfmt6/LAST blasttab files.")
    c.add_argument("--sample-id", nargs="+", help="Sample IDs matching --input. If omitted, inferred from filenames.")
    c.add_argument("--columns", default=None, help="Column names for tabular input. Default expects Fillet BLAST fields; standard 12-col BLAST also works.")
    c.add_argument("--query-fasta", nargs="+", help="Optional FASTA files matching --input; enables read length, query coverage, and embedded sequence inspection.")
    c.add_argument("--acc2taxid", help="Optional accession-to-taxid map for inputs lacking staxids, e.g. LAST blasttab.")
    add_taxonomy_args(c, required=True)
    c.add_argument("--sample-sheet", help="Optional CSV/TSV/XLSX with sample_id and role: sample, negative, positive, environmental.")
    c.add_argument("--regional-taxa", help="Optional CSV/TSV/XLSX regional taxa/prior table keyed by taxid or name.")
    c.add_argument("--contaminants-file", dest="contaminants_file", default=None,
                   help="Optional list of contaminant taxids/names (one per line or TSV with taxid column). "
                        "Matching taxa get a 'user_contaminant' quality flag in all outputs.")
    c.add_argument("--palynology-table", dest="palynology_table", default=None,
                   help="Optional palynology/pollen taxid table (TSV/CSV/XLSX or plain taxid list). "
                        "Matching taxa get 'pal_support' flag and contribute to authenticity tier.")
    c.add_argument("--fossil-table", dest="fossil_table", default=None,
                   help="Optional fossil occurrence table (TSV/CSV/XLSX or plain taxid list). "
                        "Matching taxa get 'fos_support' flag and contribute to authenticity tier.")
    c.add_argument("--depth-file", dest="depth_file", default=None,
                   help="TSV with columns sample_id and total_reads for RPM normalization. "
                        "If omitted, total reads are inferred from the alignment input.")
    c.add_argument("--reroute-taxid", dest="reroute_taxids", action="append", default=[], metavar="SRC[:DST]",
                   help="After classification, reassign all reads at SRC taxid to DST taxid (or SRC's parent if DST "
                        "omitted). Can be specified multiple times.")
    c.add_argument("--config", help="Optional TOML config overriding defaults.")
    c.add_argument("--outdir", required=True, help="Output directory.")
    c.add_argument("--prefix", default="fillet", help="Output filename prefix.")
    c.add_argument("--max-viewer-reads", type=int, default=5000, help="Max read examples embedded in standalone HTML viewer.")
    c.add_argument("--no-viewer", action="store_true", help="Skip HTML viewer generation.")
    c.add_argument("--sqlite-viewer", action="store_true", default=True, help="Write a SQLite-backed viewer database for large projects (default: on).")
    c.add_argument("--no-sqlite-viewer", dest="sqlite_viewer", action="store_false", help="Do not write the SQLite viewer database.")
    c.add_argument("--em-iterations", type=int, default=1, dest="em_iterations",
                   help="Number of sample-level EM refinement passes (default: 1). "
                        "Set to 0 to disable. EM redistributes ambiguous reads toward "
                        "taxa already supported by other reads in the same sample.")
    c.add_argument("--dirichlet-alpha", type=float, default=0.0, dest="dirichlet_alpha",
                   metavar="ALPHA",
                   help="[Experimental] Dirichlet pseudocount for EM shrinkage. "
                        "Taxa with fewer weighted reads than ALPHA receive no EM boost, "
                        "suppressing false-positive feedback from spurious assignments. "
                        "0.0 = disabled (default). Recommended: 0.5.")
    c.add_argument("--damage-weight", type=float, default=0.0, dest="damage_weight",
                   metavar="W",
                   help="[Experimental] Weight for damage-integrated likelihood scoring. "
                        "Re-scores ambiguous reads using a binomial LR: "
                        "P(damage|ancient_taxon) / P(damage|modern). "
                        "Requires qseq/sseq alignment columns (LAST output). "
                        "0.0 = disabled (default). Recommended: 0.3.")
    c.add_argument("--no-deconvolve", action="store_true", dest="no_deconvolve",
                   help="Skip within-clade Bayesian deconvolution step.")
    c.add_argument("--no-mini-db", action="store_true", dest="no_mini_db",
                   help="Skip per-clade mini-DB enrichment for deconvolution (faster, less accurate).")
    c.add_argument("--max-clade-seqs", type=int, default=5000, dest="max_clade_seqs_per_taxon",
                   metavar="N",
                   help="Max sequences/taxon for deconvolution mini-DBs (default: 5000).")
    c.add_argument("--max-clade-db-mb", type=float, default=50.0, dest="max_clade_db_mb",
                   help="Max FASTA size (MB) per mini-DB before diversity subsampling [default: 50]")
    c.add_argument("--mini-db-timeout-min", type=float, default=30.0, dest="mini_db_timeout_min",
                   help="Max minutes for mini-DB lastdb build; 0=no limit [default: 30]")
    c.add_argument("--max-parallel-mini-dbs", type=int, default=1, dest="max_parallel_mini_dbs",
                   help="Mini-DBs to build in parallel [default: 1]")
    c.add_argument("--no-curate-refs", action="store_true", dest="no_curate_refs",
                   help="Disable reference curation (length, N-content, deduplication) [default: enabled]")
    c.add_argument("--min-ref-length", type=int, default=300, dest="min_ref_length",
                   help="Minimum reference sequence length in bp; sequences shorter than this are excluded during extraction and curation (default: 300; prevents 107 bp predicted ncRNA from consuming alignment DB slots ahead of mitochondrial genomes and mRNA)")
    c.add_argument("--max-ref-n-fraction", type=float, default=0.5, dest="max_ref_n_fraction",
                   help="Maximum fraction of N/ambiguous bases in a reference sequence [default: 0.5]")
    c.add_argument("--verbose", action="store_true")
    c.add_argument(
        "--mode", default="auto", dest="mode",
        choices=["auto", "bacteria-genome", "animal-mito", "plant-chloroplast", "fungi-its", "mixed"],
        help=(
            "Classification mode. Sets evidence-quality filter defaults appropriate for the "
            "reference DB type. 'bacteria-genome': whole-genome bacterial DB — enables stk "
            "filter (0.10), min-windows (2), SES+windows (40K/15), and noise-node suppression. "
            "'animal-mito'/'plant-chloroplast': organelle DB — stk filter (0.25), no SES. "
            "'fungi-its': barcode DB — relaxed stk filter (0.70). "
            "'mixed': multi-kingdom DB — noise-node suppression only. "
            "'auto' (default): no profile, use TOML config and explicit CLI flags only."
        ),
    )
    c.add_argument(
        "--uniqueness-index", dest="uniqueness_index", default=None, metavar="PATH",
        help="SQLite uniqueness index from `fillet build-uniqueness-index`. "
             "When provided, per-alignment uniqueness weights are applied before LCA scoring: "
             "alignments to conserved loci shared by many taxa receive near-zero weight, "
             "reducing false positives from distributed conserved-locus hits (e.g. plant cp genes).",
    )
    add_lca_override_args(c)

    v = sub.add_parser("view", help="Build/rebuild the standalone HTML viewer from Fillet output TSVs.")
    v.add_argument("--assignments", required=True, help="Read assignment TSV from fillet classify.")
    v.add_argument("--summary", required=True, help="Taxon summary TSV from fillet classify.")
    add_taxonomy_args(v, required=True)
    v.add_argument("--out-html", required=True)
    v.add_argument("--max-viewer-reads", type=int, default=5000)
    v.add_argument("--out-sqlite", help="Optional path to also write a SQLite-backed viewer database.")

    i = sub.add_parser("init", help="Create a runnable skeleton directory with config and sample sheets.")
    i.add_argument("--outdir", required=True)

    a = sub.add_parser("run-align", help="Restartable BLAST or LAST alignment wrapper with chunks, ETA, resume, and Fillet-friendly output.")
    a.add_argument("--aligner", choices=["blast", "last"], default="blast")
    a.add_argument("--query", required=True, help="Input FASTA.")
    a.add_argument("--out", required=True, help="Combined output BLAST-like table.")
    a.add_argument("--workdir", help="Working directory for chunks/logs. Default: OUT.parent/fillet_align_work.")
    a.add_argument("--threads", type=int, default=16)
    a.add_argument("--chunk-size", type=int, default=50000)
    a.add_argument("--resume", action="store_true", default=True, help="Resume from completed alignment chunks (default: on).")
    a.add_argument("--no-resume", dest="resume", action="store_false", help="Rerun all chunks from scratch.")
    a.add_argument("--db", help="BLAST DB basename; required for --aligner blast.")
    a.add_argument("--task", default="blastn", choices=["blastn", "megablast"])
    a.add_argument("--evalue", default="1e-5")
    a.add_argument("--max-target-seqs", type=int, default=2000)
    a.add_argument("--max-hsps", type=int, default=1)
    a.add_argument("--with-inspection-fields", type=int, default=1, help="1 includes qlen/qcovhsp/slen/qseq/sseq/stitle/btop/sstrand for the viewer.")
    a.add_argument("--last-db-prefix", help="LAST DB prefix; required for --aligner last.")
    a.add_argument("--last-train-file", help="Optional last-train model file.")
    a.add_argument("--last-threads", type=int, help="LAST -P threads; default uses --threads.")
    a.add_argument("--last-output-format", choices=["tab", "maf"], default="tab",
                   help="'tab' (blasttab via maf-convert, default) or 'maf' (raw MAF — preserves qseq/sseq for damage detection).")

    # Back-compatible alias. It now uses the chunked/resumable BLAST engine.
    b = sub.add_parser("run-blast", help="Restartable convenience wrapper around BLASTn/megablast that writes Fillet-friendly outfmt6.")
    b.add_argument("--query", required=True)
    b.add_argument("--db", required=True)
    b.add_argument("--out", required=True)
    b.add_argument("--workdir", help="Working directory for chunks/logs. Default: OUT.parent/fillet_blast_work.")
    b.add_argument("--task", default="blastn", choices=["blastn", "megablast"])
    b.add_argument("--threads", type=int, default=16)
    b.add_argument("--evalue", default="1e-5")
    b.add_argument("--max-target-seqs", type=int, default=2000)
    b.add_argument("--max-hsps", type=int, default=1)
    b.add_argument("--chunk-size", type=int, default=50000)
    b.add_argument("--resume", action="store_true", default=True)
    b.add_argument("--no-resume", dest="resume", action="store_false")
    b.add_argument("--with-inspection-fields", dest="with_inspection_fields", action="store_true", default=True, help="Include qlen/qcovhsp/slen/qseq/sseq/stitle/btop/sstrand so the GUI can show MEGAN-like pairwise read/reference alignments. Default: on.")
    b.add_argument("--no-inspection-fields", dest="with_inspection_fields", action="store_false", help="Disable qseq/sseq/stitle/BTOP output for smaller BLAST tables.")

    w = sub.add_parser("workflow", help="Run restartable BLAST/LAST on FASTA and immediately classify with RELIC-LCA/viewer outputs.")
    w.add_argument("--aligner", choices=["blast", "last"], default="blast")
    w.add_argument("--query", required=True, help="Input FASTA for this sample/library.")
    w.add_argument("--sample-id", required=True)
    w.add_argument("--outdir", required=True)
    w.add_argument("--prefix", default=None, help="Output prefix. Default: sample-id.")
    w.add_argument("--threads", type=int, default=16)
    w.add_argument("--chunk-size", type=int, default=50000)
    w.add_argument("--resume", action="store_true", default=True)
    w.add_argument("--no-resume", dest="resume", action="store_false")
    w.add_argument("--db", help="BLAST DB basename; required for --aligner blast.")
    w.add_argument("--task", default="blastn", choices=["blastn", "megablast"])
    w.add_argument("--evalue", default="1e-5")
    w.add_argument("--max-target-seqs", type=int, default=2000)
    w.add_argument("--max-hsps", type=int, default=1)
    w.add_argument("--last-db-prefix", help="LAST DB prefix; required for --aligner last.")
    w.add_argument("--last-train-file", help="Optional last-train model file.")
    w.add_argument("--last-threads", type=int)
    w.add_argument("--last-output-format", choices=["tab", "maf"], default="tab",
                   help="'tab' (blasttab via maf-convert, default) or 'maf' (raw MAF — preserves qseq/sseq for damage detection).")
    w.add_argument("--acc2taxid", help="Optional accession-to-taxid map, required/recommended for LAST.")
    add_taxonomy_args(w, required=True)
    w.add_argument("--sample-sheet", help="Optional sample/control sheet.")
    w.add_argument("--regional-taxa", help="Optional regional/ecological prior table.")
    w.add_argument("--config", help="Optional TOML config overriding defaults.")
    w.add_argument("--max-viewer-reads", type=int, default=5000)
    w.add_argument("--no-viewer", action="store_true")
    w.add_argument("--sqlite-viewer", action="store_true", default=True, help="Write a SQLite-backed viewer database (default: on).")
    w.add_argument("--no-sqlite-viewer", dest="sqlite_viewer", action="store_false")
    w.add_argument("--em-iterations", type=int, default=1, dest="em_iterations",
                   help="Number of EM refinement passes (default: 1). Set to 0 to disable.")
    w.add_argument("--verbose", action="store_true")
    w.add_argument("--yes", "-y", action="store_true", help="Skip the settings confirmation prompt and run immediately.")
    w.add_argument(
        "--mode", default="auto", dest="mode",
        choices=["auto", "bacteria-genome", "animal-mito", "plant-chloroplast", "fungi-its", "mixed"],
        help="Classification mode profile. See 'fillet classify --help' for details.",
    )
    add_lca_override_args(w)


    # ── plot ───────────────────────────────────────────────────────────────────
    pl = sub.add_parser("plot", help="Generate publication-quality plots from taxon_summary.tsv files.")
    pl.add_argument("--input", "-i", nargs="+", required=True, metavar="TSV",
                    help="One or more taxon_summary.tsv files (from fillet classify).")
    pl.add_argument("--type", dest="plot_type", default="bar",
                    choices=["bar", "stacked", "bubble", "heatmap", "strat", "strat-climate", "damage"],
                    help="Plot type: bar, stacked (bar), bubble, heatmap, strat (stratigraphic), "
                         "strat-climate (strat + climate proxy), damage. Default: bar.")
    pl.add_argument("--out", "-o", required=True, metavar="FILE",
                    help="Output file path. Extension determines format: .png or .svg.")
    pl.add_argument("--metric", default="cumulative_hard_reads",
                    help="Column to use for plot values. Default: cumulative_hard_reads.")
    pl.add_argument("--rank", default=None, metavar="RANK",
                    help="Filter to a specific rank (e.g. family, genus, species).")
    pl.add_argument("--top-n", type=int, default=20, metavar="N",
                    help="Number of top taxa to show. Default: 20.")
    pl.add_argument("--min-reads", type=int, default=1, metavar="N",
                    help="Exclude taxa with fewer than N total reads across all samples.")
    pl.add_argument("--taxa", nargs="+", metavar="TAXON",
                    help="Explicit list of taxon names to include (overrides --top-n).")
    pl.add_argument("--samples", nargs="+", metavar="SAMPLE",
                    help="Explicit ordered list of sample IDs to include.")
    pl.add_argument("--title", default=None, help="Plot title.")
    pl.add_argument("--normalize", action="store_true", default=False,
                    help="Normalize values (relative abundance or percent composition).")
    pl.add_argument("--log-scale", action="store_true", default=False,
                    help="Use log scale for heatmap colour.")
    pl.add_argument("--dpi", type=int, default=150, help="DPI for PNG output. Default: 150.")
    pl.add_argument("--width", type=float, default=None, help="Figure width in inches.")
    pl.add_argument("--height", type=float, default=None, help="Figure height in inches.")
    pl.add_argument("--climate-file", dest="climate_file", default=None, metavar="CSV",
                    help="User-supplied climate proxy CSV/TSV (columns: depth/age, value).")
    pl.add_argument("--climate-label", dest="climate_label", default="Climate proxy",
                    help="Label for user-supplied climate proxy panel.")
    pl.add_argument("--show-lr04", dest="show_lr04", action="store_true", default=False,
                    help="Overlay LR04 benthic δ¹⁸O stack (Lisiecki & Raymo 2005).")
    pl.add_argument("--show-gisp2-b", dest="show_gisp2_b", action="store_true", default=False,
                    help="Overlay GISP2 B-core annual accumulation (Meese et al. 1997, ~700 cal BP).")
    pl.add_argument("--show-gisp2-d", dest="show_gisp2_d", action="store_true", default=False,
                    help="Overlay GISP2 D-core depth-age model (Meese et al. 1997, 0–163 ka).")
    pl.add_argument("--show-mis", dest="show_mis", action="store_true", default=False,
                    help="Shade MIS stage boundaries on all panels.")
    pl.add_argument("--use-age", dest="use_age", action="store_true", default=False,
                    help="Use age_BP column for y-axis (default: depth_cm).")
    pl.add_argument("--14c-dates", dest="dates_14c", default=None, metavar="FILE",
                    help="TSV/CSV of pre-calibrated ¹⁴C dates to overlay as error bars. "
                         "Columns: name, depth_cm, median_BP, lo68, hi68, lo95, hi95.")

    sv = sub.add_parser("serve-viewer", help="Serve a SQLite-backed Fillet viewer locally in your browser.")
    sv.add_argument("--db", required=True, help="Fillet .viewer.sqlite database from classify/workflow/view.")
    sv.add_argument("--host", default="127.0.0.1")
    sv.add_argument("--port", type=int, default=8765)

    g = sub.add_parser("gui", help="Launch the Fillet interactive GUI viewer (requires PyQt6).")
    g.add_argument("db", nargs="?", help="Optional .viewer.sqlite database to open on startup.")

    # ── build-db ──────────────────────────────────────────────────────────────
    bd = sub.add_parser(
        "build-db",
        help="Build a targeted LAST database for a set of taxids (subset of NCBI NT).",
        formatter_class=__import__("argparse").RawDescriptionHelpFormatter,
        description="""
Build a single-volume LAST database from NCBI NT sequences that belong to the
specified taxa (and optionally all their descendants).

Steps:
  1. Expand seed taxids to all descendants via nodes.dmp (BFS)
  2. Extract matching sequences from NCBI NT (BLAST+ format) via blastdbcmd
  3. Build a LAST database with lastdb -P threads -Q 0

The resulting database can be searched with:
  fillet run-align --aligner last --last-db-prefix OUT_PREFIX ...

Example:
  fillet build-db --taxids 9612 40674 \\
    --nodes /path/to/nodes.dmp \\
    --blast-db-dir /path/to/NCBI_NT/2025-08-05 \\
    --out-prefix /path/to/mydb/mammals_last \\
    --threads 32
""",
    )
    bd.add_argument(
        "--taxids", "-t", nargs="+", metavar="TAXID",
        help="One or more seed NCBI taxids. All descendants are included unless --no-expand.",
    )
    bd.add_argument(
        "--taxid-file", metavar="FILE",
        help="File with one taxid per line (alternative to --taxids).",
    )
    bd.add_argument(
        "--nodes", required=True, metavar="NODES_DMP",
        help="Path to NCBI taxonomy nodes.dmp.",
    )
    bd.add_argument(
        "--blast-db-dir", required=True, metavar="DIR",
        help="Directory containing NCBI NT BLAST+ database volumes (sets BLASTDB env).",
    )
    bd.add_argument(
        "--out-prefix", required=True, metavar="PREFIX",
        help="Output prefix for LAST database files and extracted FASTA.",
    )
    bd.add_argument(
        "--no-expand", action="store_true", default=False,
        help="Do not expand seed taxids to descendants (use seeds only).",
    )
    bd.add_argument(
        "--threads", type=int, default=16, metavar="N",
        help="Threads for lastdb indexing (default: 16).",
    )
    bd.add_argument(
        "--verbose", "-v", action="store_true", default=True,
    )

    # ── metagenome ────────────────────────────────────────────────────────────
    met = sub.add_parser(
        "metagenome",
        help="Full metagenomics pipeline: screen (KrakenUniq/Kraken2) → subset LAST DB → align → classify.",
        formatter_class=__import__("argparse").RawDescriptionHelpFormatter,
        description="""
Full metagenomics pipeline (screen → build DB → align → classify):

  Stage 1  Kraken[Uniq] screening — fast first-pass identification of candidate taxa
  Stage 2  Subset LAST DB build — extract matching NT sequences, build single-volume LAST DB
  Stage 3  LAST alignment — sensitive MAF alignment with damage-preserving output
  Stage 4  RELIC-LCA classification — EM-refined LCA with damage detection

Example:
  fillet metagenome \\
    --query reads.fasta --sample-id NHB-7 --outdir /results/NHB-7 \\
    --kraken-db /data/krakenuniq_db \\
    --blast-db-dir /data/NCBI_NT/2026-03-31 \\
    --nodes /data/taxonomy/nodes.dmp --names /data/taxonomy/names.dmp \\
    --last-train-file /data/LAST/nt_adna.train \\
    --threads 32
""",
    )
    met_input = met.add_argument_group(
        "Input (provide one of: --query, --fastq-r1/r2, or --bam)"
    )
    met_input.add_argument("--query", metavar="FASTA",
                           help="Pre-processed FASTA input. Skip Stage 0 if provided.")
    met_input.add_argument("--fastq-r1", metavar="R1.fastq.gz",
                           help="Raw paired-end R1 FASTQ (gzip OK). Triggers Stage 0 preprocessing.")
    met_input.add_argument("--fastq-r2", metavar="R2.fastq.gz",
                           help="Raw paired-end R2 FASTQ. Omit for single-end / pre-merged.")
    met_input.add_argument("--bam", metavar="FILE.bam",
                           help="BAM file from nf-core/eager or BWA — reads extracted via samtools. "
                                "No re-trimming: BAM reads are assumed to be already processed by eager.")
    met.add_argument("--sample-id", required=True)
    met.add_argument("--outdir", required=True)
    met.add_argument("--threads", type=int, default=16)

    met_resume = met.add_argument_group("Step-into / resume controls")
    met_resume.add_argument(
        "--from-stage",
        choices=["screen", "build-db", "align", "classify"],
        default=None, dest="from_stage", metavar="STAGE",
        help=(
            "Start pipeline at STAGE, skipping earlier stages (prerequisite files must exist). "
            "  classify — re-run RELIC-LCA only; needs --query + existing MAF; "
            "fastest option for iterating on LCA parameters. "
            "  align — re-run LAST alignment + classification; needs --query + existing LAST DB. "
            "  build-db — re-build LAST DB + align + classify; needs --query + --taxid-file. "
            "  screen — skip preprocessing only; needs --query."
        ),
    )
    met_resume.add_argument(
        "--maf", metavar="FILE", default=None,
        help=(
            "Explicit MAF alignment file path for --from-stage classify. "
            "Defaults to {outdir}/{sample-id}.last.maf if omitted."
        ),
    )

    met_prep = met.add_argument_group("Stage 0: Read preprocessing (FASTQ/BAM input only)")
    met_prep.add_argument("--preprocess-tool", choices=["auto", "fastp", "AdapterRemoval"], default="auto",
                          help="Trimming/merging tool (default: auto, prefers fastp).")
    met_prep.add_argument("--platform", metavar="PLATFORM", default=None,
                          help="Sequencing platform (e.g. NextSeq1000, HiSeq2500). "
                               "2-colour Illumina platforms (NextSeq/NovaSeq) automatically enable poly-G trimming.")
    met_prep.add_argument("--poly-g-trim", dest="poly_g_trim", action="store_true", default=False,
                          help="Force poly-G tail trimming (auto-enabled for 2-colour platforms).")
    met_prep.add_argument("--min-length", type=int, default=24, metavar="N",
                          help="Minimum read length after trimming (default: 24 bp).")
    met_prep.add_argument("--no-dedup", dest="dedup", action="store_false", default=True,
                          help="Skip exact-sequence deduplication.")
    met_prep.add_argument("--no-adapter-check", dest="adapter_check", action="store_false", default=True,
                          help="Skip extra Hamming-1 adapter remnant filter.")

    met_screen = met.add_argument_group("Stage 1: Taxonomic screening")
    met_screen.add_argument("--kraken-db", metavar="DIR",
                            help="KrakenUniq or Kraken2 database directory. Required unless --taxid-file is given.")
    met_screen.add_argument("--screener", choices=["auto", "krakenuniq", "kraken2"], default="auto",
                            help="Screening tool (default: auto, prefers krakenuniq).")
    met_screen.add_argument("--min-kraken-reads", type=int, default=5, metavar="N",
                            help="Min reads for a taxon to be included (default: 5).")
    met_screen.add_argument("--min-kraken-unique-kmers", type=int, default=50, metavar="N",
                            help="Min unique k-mers (KrakenUniq only; default: 50). Filters conserved-sequence "
                                 "false positives — taxa found only via a handful of repeated k-mers shared "
                                 "across distant lineages.")
    met_screen.add_argument("--max-kraken-dup-rate", type=float, default=None, metavar="X",
                            help="Max k-mer duplication rate (KrakenUniq only; default: no limit). "
                                 "Rejects taxa whose k-mers are highly repetitive (e.g. >20 means each "
                                 "unique k-mer was seen 20+ times on average).")
    met_screen.add_argument("--taxid-file", metavar="FILE",
                            help="Pre-computed taxid list (one per line); skips Kraken screening.")

    met_db = met.add_argument_group("Stage 2: Subset LAST database")
    met_db.add_argument("--nodes", required=True, metavar="NODES_DMP",
                        help="NCBI taxonomy nodes.dmp — required for taxonomy expansion and classification.")
    met_db.add_argument("--blast-db-dir", metavar="DIR",
                        help="Directory containing NCBI NT BLAST+ database volumes (sets BLASTDB env). "
                             "Required unless --skip-build-db and a pre-built DB exists at --last-db-prefix.")
    met_db.add_argument("--custom-ref-fasta", nargs="+", metavar="FASTA",
                        help="One or more FASTA files to merge into the LAST DB alongside NT sequences. "
                             "Use this for FlyGuide/Spinner curated reference databases, "
                             "local species not in NT, or any supplementary reference sequences.")
    met_db.add_argument("--last-db-prefix", metavar="PREFIX",
                        help="LAST DB prefix for subset database. Default: OUTDIR/subset_last_db/SAMPLE_ID. "
                             "Point multiple samples at the same prefix to reuse a shared DB.")
    met_db.add_argument("--skip-build-db", action="store_true",
                        help="Skip DB build if LAST DB already exists at --last-db-prefix.")
    met_db.add_argument("--max-seed-rank", default="family", metavar="RANK",
                        help="Maximum taxonomic rank allowed as a seed for DB expansion (default: family). "
                             "Seeds above this rank (order, class, phylum, …) are dropped before BFS expansion "
                             "to prevent a single high-rank k-mer hit from pulling in millions of sequences.")

    met_db.add_argument("--no-chunked-db", dest="chunked_db", action="store_false", default=True,
                        help="Disable hierarchical chunked DB; use single monolithic LAST DB (legacy, "
                             "requires --max-seqs-per-taxon to avoid OOM). Default: chunked on.")
    met_db.add_argument("--chunked-db-rank", default=None, dest="chunked_db_rank", metavar="RANK",
                        help="Force partition rank for chunked DB (e.g. family, order, class). "
                             "Default: auto-select based on available RAM.")
    met_db.add_argument("--parallel-chunks", type=int, default=1, dest="parallel_chunks", metavar="N",
                        help="Number of chunk DBs to build / align in parallel (default: 1). "
                             "Each chunk uses threads/N threads. Increase on many-core servers.")
    met_db.add_argument("--chunk-timeout", type=int, default=1800, dest="chunk_timeout_sec", metavar="SEC",
                        help="Seconds to wait for user rank selection before auto-selecting safe option "
                             "(default: 1800 = 30 min).")
    met_db.add_argument("--no-keep-chunk-mafs", dest="keep_chunk_mafs", action="store_false", default=True,
                        help="Delete per-chunk MAF files after merging (saves disk, loses resumability).")

    met_align = met.add_argument_group("Stage 3: LAST alignment")
    met_align.add_argument("--last-train-file", metavar="FILE",
                           help="LAST training file (last-train output) — strongly recommended for aDNA.")
    met_align.add_argument("--last-min-score", type=int, metavar="N", dest="last_min_score",
                           help="Minimum lastal gapped alignment score (-e). Lower values recover "
                                "damaged alignments that fall below LAST's auto threshold. "
                                "Recommended for ancient DNA: 150–180. Default: LAST automatic.")
    met_align.add_argument("--last-m", type=int, metavar="N", dest="last_m",
                           help="Max initial seed matches per query position (-m). LAST default=10. "
                                "Increase to 1000–10000 for small reference DBs (e.g. organelles-only) "
                                "to avoid missing self-hits at conserved positions (12S/16S rRNA, CO1). "
                                "Do NOT increase for large DBs (full NT) — catastrophically slow. "
                                "Default: LAST automatic (10).")
    met_align.add_argument("--last-threads", type=int, metavar="N",
                           help="Threads for lastal (default: --threads).")
    met_align.add_argument("--evalue", default="1e-5")
    met_align.add_argument("--max-target-seqs", type=int, default=2000)
    met_align.add_argument("--chunk-size", type=int, default=50000, metavar="N",
                           help="Reads per alignment chunk (default: 50000).")
    met_align.add_argument("--resume", action="store_true", default=True,
                           help="Resume from completed alignment chunks (default: on).")
    met_align.add_argument("--no-resume", dest="resume", action="store_false",
                           help="Rerun all alignment chunks from scratch.")

    met_clf = met.add_argument_group("Stage 4: RELIC-LCA classification")
    met_clf.add_argument("--names", metavar="NAMES_DMP",
                         help="NCBI names.dmp (used with --nodes for taxonomy). "
                              "Alternative: --taxonomy-tsv.")
    met_clf.add_argument("--taxonomy-tsv", metavar="FILE",
                         help="Curated taxonomy TSV (alternative to --nodes/--names).")
    met_clf.add_argument("--root-taxid", default="1")
    met_clf.add_argument("--sample-sheet", metavar="FILE")
    met_clf.add_argument("--regional-taxa", metavar="FILE")
    met_clf.add_argument("--config", metavar="TOML",
                         help="Optional TOML config overriding RELIC-LCA defaults.")
    met_clf.add_argument("--em-iterations", type=int, default=1)
    met_clf.add_argument("--dirichlet-alpha", type=float, default=0.0, dest="dirichlet_alpha",
                         metavar="ALPHA",
                         help="[Experimental] Dirichlet pseudocount for EM shrinkage (default: 0.0 = off).")
    met_clf.add_argument("--damage-weight", type=float, default=0.0, dest="damage_weight",
                         metavar="W",
                         help="[Experimental] Damage likelihood weight (default: 0.0 = off). Requires LAST output.")
    met_clf.add_argument("--no-deconvolve", action="store_true", dest="no_deconvolve",
                         help="Skip within-clade Bayesian deconvolution step.")
    met_clf.add_argument("--no-mini-db", action="store_true", dest="no_mini_db",
                         help="Skip per-clade mini-DB enrichment for deconvolution (faster, less accurate).")
    met_clf.add_argument("--max-clade-seqs", type=int, default=5000, dest="max_clade_seqs_per_taxon",
                         metavar="N",
                         help="Max sequences/taxon for deconvolution mini-DBs (default: 5000).")
    met_clf.add_argument("--max-clade-db-mb", type=float, default=50.0, dest="max_clade_db_mb",
                         help="Max FASTA size (MB) per mini-DB before diversity subsampling [default: 50]")
    met_clf.add_argument("--mini-db-timeout-min", type=float, default=30.0, dest="mini_db_timeout_min",
                         help="Max minutes for mini-DB lastdb build; 0=no limit [default: 30]")
    met_clf.add_argument("--max-parallel-mini-dbs", type=int, default=1, dest="max_parallel_mini_dbs",
                         help="Mini-DBs to build in parallel [default: 1]")
    met_clf.add_argument("--no-curate-refs", action="store_true", dest="no_curate_refs",
                         help="Disable reference curation (length, N-content, deduplication) [default: enabled]")
    met_clf.add_argument("--min-ref-length", type=int, default=300, dest="min_ref_length",
                         help="Minimum reference sequence length in bp; sequences shorter than this are excluded during extraction and curation (default: 300; prevents 107 bp predicted ncRNA from consuming alignment DB slots ahead of mitochondrial genomes and mRNA)")
    met_clf.add_argument("--max-ref-n-fraction", type=float, default=0.5, dest="max_ref_n_fraction",
                         help="Maximum fraction of N/ambiguous bases in a reference sequence [default: 0.5]")
    met_clf.add_argument(
        "--ref-scope", default="all", dest="ref_scope",
        choices=["all", "organelle", "markers", "organelle+markers"],
        help=(
            "Filter references to only the specified genomic scope after NT extraction "
            "[default: all]. 'organelle': mitochondrial/plastid seqs + NC_ accessions "
            "(bacterial complete genomes). 'markers': COI/16S/18S/ITS/rbcL/matK/cytb/ND2/12S/trnL. "
            "'organelle+markers': union of both. Dramatically shrinks the LAST DB for "
            "eukaryote-biased samples without losing microbial coverage (NC_ are kept). "
            "Superseded by --ref-mode when that option is specified."
        ),
    )
    met_clf.add_argument(
        "--ref-mode", default="full", dest="ref_mode",
        choices=["full", "organelles-markers", "organelles-only", "organelles-nuclear"],
        help=(
            "Two-pass tier-sorted reference selection mode [default: full]. "
            "Sequences are classified into priority tiers before extraction: "
            "(1) organelle genomes (mito/chloroplast), "
            "(2) RefSeq complete chromosomes/bacterial genomes, "
            "(3) RefSeq transcripts (NM_/NR_), "
            "(4) barcoding markers (COI/ITS/16S/rbcL/matK/cytb), "
            "(5) long GenBank (≥1kb), "
            "(6) predicted sequences (XM_/XR_, capped). "
            "organelles-nuclear: like organelles-markers but adds nuclear_chrom_cap=1 "
            "with max 15 Mb — admits complete bacterial/archaeal genomes while blocking "
            "large eukaryote chromosomes. "
            "Scaffold/unplaced contigs are excluded in all modes. "
            "'full': all tiers with per-type sub-caps (best for complex samples). "
            "'organelles-markers': organelle genomes + barcode markers only "
            "(smaller DB, good for aDNA with expected organelle coverage). "
            "'organelles-only': organelle genomes only (minimal DB, fastest build). "
            "When --ref-mode is set, it replaces --ref-scope filtering."
        ),
    )
    met_clf.add_argument(
        "--nt-index", default=None, dest="nt_index", metavar="FILE",
        help=(
            "Path to a pre-built NT accession index (nt_accession_index.tsv.gz) "
            "created by 'fillet build-nt-index'. When supplied, the metadata pass "
            "reads from this file instead of querying blastdbcmd, which is faster "
            "for repeated DB builds against the same NT version. Optional."
        ),
    )
    met_clf.add_argument(
        "--custom-ref-taxid-map", default=None, dest="custom_ref_taxid_map", metavar="FILE",
        help=(
            "File listing taxids present in --custom-ref-fasta (one per line). "
            "Those taxids are excluded from NT extraction to avoid duplicating sequences. "
            "Not needed when custom FASTA headers contain taxid= or |taxid: tags."
        ),
    )
    met_clf.add_argument("--max-seqs-per-taxon", type=int, default=None, dest="max_seqs_per_taxon",
                         metavar="N",
                         help="Cap sequences extracted per taxon when building the subset LAST DB. "
                              "Prevents memory blowout on taxa with thousands of NT entries. "
                              "Default: auto (set based on available RAM). Use 0 to disable the cap.")
    met_clf.add_argument(
        "--dedup-similarity", default=0.90, type=float, dest="dedup_similarity",
        metavar="FLOAT",
        help="Sourmash Jaccard similarity threshold for removing near-duplicate reference "
             "sequences within each taxid (0–1, default 0.90). Only active with --ref-mode. "
             "Set to 1.0 to disable.",
    )
    met_clf.add_argument("--max-viewer-reads", type=int, default=5000)
    met_clf.add_argument("--no-viewer", action="store_true")
    met_clf.add_argument("--no-sqlite-viewer", dest="sqlite_viewer", action="store_false", default=True)
    met_clf.add_argument(
        "--uniqueness-index", dest="uniqueness_index", default=None, metavar="PATH",
        help="SQLite uniqueness index (built with `fillet build-uniqueness-index`). "
             "Applies per-alignment uniqueness weights before LCA scoring — "
             "suppresses conserved-locus FPs (e.g. cp gene reads assigned across plant families). "
             "Pass the path to the .uniqueness.sqlite3 file produced for the Eukaryota chunk.",
    )
    met_clf.add_argument("--verbose", "-v", action="store_true", default=True)
    met_clf.add_argument("--yes", "-y", action="store_true",
                         help="Skip the settings confirmation prompt.")
    met_clf.add_argument(
        "--mode", default="auto", dest="mode",
        choices=["auto", "bacteria-genome", "animal-mito", "plant-chloroplast", "fungi-its", "mixed"],
        help="Classification mode profile. See 'fillet classify --help' for details.",
    )
    add_lca_override_args(met)

    # ── batch-metagenome ──────────────────────────────────────────────────────
    bmet = sub.add_parser(
        "batch-metagenome",
        help="Run the metagenomics pipeline for multiple samples sharing one LAST DB per group.",
        formatter_class=__import__("argparse").RawDescriptionHelpFormatter,
        description="""
Run the full metagenomics pipeline for a batch of samples.

Samples with the same group_id share a single LAST DB build — the screening
results for all group members are unioned before building, so each group only
pays the expensive lastdb cost once.

Sample sheet TSV columns (tab-separated, header required):
  sample_id   — required
  outdir      — required (per-sample output directory)
  group_id    — optional (default: 'default'); samples with the same group_id share a DB
  query       — preprocessed FASTA; if absent, fastq_r1 must be set
  fastq_r1    — R1 FASTQ (triggers preprocessing + screening)
  fastq_r2    — R2 FASTQ (optional, paired-end)

Example:
  fillet batch-metagenome --samples batch.tsv \\
    --shared-db-dir /results/shared_dbs \\
    --kraken-db /data/krakenuniq_nt_db \\
    --blast-db-dir /data/NCBI_NT/2026-03-31 \\
    --nodes /data/taxonomy/nodes.dmp --names /data/taxonomy/names.dmp \\
    --threads 32
""",
    )
    bmet.add_argument("--samples", required=True, metavar="TSV",
                      help="Sample sheet TSV (sample_id, outdir, group_id, query/fastq_r1/fastq_r2).")
    bmet.add_argument("--shared-db-dir", required=True, metavar="DIR", dest="shared_db_dir",
                      help="Directory where shared LAST DBs are written (one sub-dir per group).")
    bmet.add_argument("--threads", type=int, default=16)
    bmet.add_argument("--nodes", required=True, metavar="NODES_DMP")
    bmet.add_argument("--names", metavar="NAMES_DMP")
    bmet.add_argument("--kraken-db", metavar="DIR", dest="kraken_db")
    bmet.add_argument("--screener", choices=["auto", "krakenuniq", "kraken2"], default="auto")
    bmet.add_argument("--min-kraken-reads", type=int, default=5, dest="min_kraken_reads")
    bmet.add_argument("--min-kraken-unique-kmers", type=int, default=50, dest="min_kraken_unique_kmers")
    bmet.add_argument("--max-kraken-dup-rate", type=float, default=None, dest="max_kraken_dup_rate")
    bmet.add_argument("--blast-db-dir", metavar="DIR", dest="blast_db_dir")
    bmet.add_argument("--max-seed-rank", default="family", dest="max_seed_rank")
    bmet.add_argument("--max-seqs-per-taxon", type=int, default=None, dest="max_seqs_per_taxon")
    bmet.add_argument("--em-iterations", type=int, default=1, dest="em_iterations")
    bmet.add_argument("--no-chunked-db", dest="chunked_db", action="store_false", default=True)
    bmet.add_argument("--chunked-db-rank", default=None, dest="chunked_db_rank")
    bmet.add_argument("--parallel-chunks", type=int, default=1, dest="parallel_chunks")
    bmet.add_argument("--last-train-file", metavar="FILE", dest="last_train_file")
    bmet.add_argument("--last-m", type=int, default=None, metavar="N", dest="last_m",
                      help="Max initial matches per query position in LAST (-m). "
                           "Default: LAST default (10). Use 1000 for large DBs, 10000 for small organelle DBs.")
    bmet.add_argument("--last-min-score", type=int, default=None, metavar="N", dest="last_min_score",
                      help="Minimum gapped alignment score for LAST (-e). "
                           "Use 150 with --damage-mode ancient to tighten specificity.")
    bmet.add_argument("--evalue", default="1e-5")
    bmet.add_argument("--max-target-seqs", type=int, default=2000, dest="max_target_seqs")
    bmet.add_argument("--chunk-size", type=int, default=50000, dest="chunk_size")
    bmet.add_argument("--resume", action="store_true", default=True)
    bmet.add_argument("--no-resume", dest="resume", action="store_false")
    bmet.add_argument("--config", metavar="TOML")
    bmet.add_argument("--no-deconvolve", action="store_true", dest="no_deconvolve")
    bmet.add_argument("--no-mini-db", action="store_true", dest="no_mini_db")
    bmet.add_argument("--max-clade-seqs", type=int, default=5000, dest="max_clade_seqs_per_taxon",
                      metavar="N",
                      help="Max sequences/taxon for deconvolution mini-DBs (default: 5000).")
    bmet.add_argument("--max-clade-db-mb", type=float, default=50.0, dest="max_clade_db_mb",
                      help="Max FASTA size (MB) per mini-DB before diversity subsampling [default: 50]")
    bmet.add_argument("--mini-db-timeout-min", type=float, default=30.0, dest="mini_db_timeout_min",
                      help="Max minutes for mini-DB lastdb build; 0=no limit [default: 30]")
    bmet.add_argument("--max-parallel-mini-dbs", type=int, default=1, dest="max_parallel_mini_dbs",
                      help="Mini-DBs to build in parallel [default: 1]")
    bmet.add_argument("--no-curate-refs", action="store_true", dest="no_curate_refs",
                      help="Disable reference curation (length, N-content, deduplication) [default: enabled]")
    bmet.add_argument("--min-ref-length", type=int, default=300, dest="min_ref_length",
                      help="Minimum reference sequence length in bp; sequences shorter than this are excluded during extraction and curation (default: 300; prevents 107 bp predicted ncRNA from consuming alignment DB slots ahead of mitochondrial genomes and mRNA)")
    bmet.add_argument("--max-ref-n-fraction", type=float, default=0.5, dest="max_ref_n_fraction",
                      help="Maximum fraction of N/ambiguous bases in a reference sequence [default: 0.5]")
    bmet.add_argument(
        "--ref-scope", default="all", dest="ref_scope",
        choices=["all", "organelle", "markers", "organelle+markers"],
        help="Filter references to genomic scope after extraction [default: all]. Superseded by --ref-mode.",
    )
    bmet.add_argument(
        "--ref-mode", default="full", dest="ref_mode",
        choices=["full", "organelles-markers", "organelles-only", "organelles-nuclear"],
        help="Two-pass tier-sorted reference selection mode [default: full]. See 'fillet metagenome --help'.",
    )
    bmet.add_argument(
        "--nt-index", default=None, dest="nt_index", metavar="FILE",
        help="Path to pre-built NT accession index (nt_accession_index.tsv.gz). Optional.",
    )
    bmet.add_argument(
        "--custom-ref-taxid-map", default=None, dest="custom_ref_taxid_map", metavar="FILE",
        help="File of taxids in --custom-ref-fasta to exclude from NT extraction.",
    )
    bmet.add_argument(
        "--dedup-similarity", default=0.90, type=float, dest="dedup_similarity",
        metavar="FLOAT",
        help="Sourmash Jaccard similarity threshold for removing near-duplicate reference "
             "sequences within each taxid (0–1, default 0.90). Only active with --ref-mode. "
             "Set to 1.0 to disable.",
    )
    bmet.add_argument("--no-viewer", action="store_true", dest="no_viewer")
    bmet.add_argument("--no-sqlite-viewer", dest="sqlite_viewer", action="store_false", default=True)
    bmet.add_argument(
        "--uniqueness-index", dest="uniqueness_index", default=None, metavar="PATH",
        help="SQLite uniqueness index (from `fillet build-uniqueness-index`). "
             "Passed to all samples in the batch; downweights conserved-locus hits before LCA scoring.",
    )
    bmet.add_argument("--yes", "-y", action="store_true")
    bmet.add_argument("--verbose", "-v", action="store_true", default=True)
    bmet.add_argument(
        "--mode", default="auto", dest="mode",
        choices=["auto", "bacteria-genome", "animal-mito", "plant-chloroplast", "fungi-its", "mixed"],
        help="Classification mode profile. See 'fillet classify --help' for details.",
    )
    add_lca_override_args(bmet)

    mg = sub.add_parser(
        "merge",
        help="Merge two or more Fillet viewer databases (.viewer.sqlite) into one.",
    )
    mg.add_argument(
        "--input", "-i",
        nargs="+",
        required=True,
        metavar="DB",
        help="Two or more .viewer.sqlite files to merge (in order; first DB's taxonomy wins).",
    )
    mg.add_argument(
        "--output", "-o",
        required=True,
        metavar="OUT",
        help="Path for the merged output .viewer.sqlite.",
    )
    mg.add_argument(
        "--select-samples", nargs="+", metavar="SAMPLE",
        help="Include only these sample IDs in the merged output (applied per source).",
    )
    mg.add_argument(
        "--exclude-samples", nargs="+", metavar="SAMPLE",
        help="Exclude these sample IDs from the merged output.",
    )

    # ── batch-status ───────────────────────────────────────────────────────────
    bs = sub.add_parser(
        "batch-status",
        help="Show batch progress summary without re-running.",
    )
    bs.add_argument("--outdir", "-o", required=True, help="Batch output directory (contains batch_progress.json).")

    # ── update-taxonomy ────────────────────────────────────────────────────────
    ut = sub.add_parser(
        "update-taxonomy",
        help="Download a fresh NCBI taxdump and extract nodes.dmp / names.dmp.",
    )
    ut.add_argument(
        "--outdir", "-o", required=True,
        help="Directory to extract into. A dated subdirectory (YYYY-MM-DD) is created automatically.",
    )
    ut.add_argument(
        "--no-date-subdir", action="store_true", default=False,
        help="Write directly into --outdir instead of creating a YYYY-MM-DD subdirectory.",
    )
    ut.add_argument(
        "--keep-tarball", action="store_true", default=False,
        help="Keep the downloaded taxdump.tar.gz after extraction.",
    )

    # ── build-nt-index ────────────────────────────────────────────────────────
    bni = sub.add_parser(
        "build-nt-index",
        help="Build a pre-computed NT accession metadata index for faster tier-sorted DB builds.",
        description=(
            "Runs blastdbcmd once over the entire NT database to extract accession, taxid, "
            "length, and title metadata for every sequence (~2-4 hours for 124M sequences). "
            "The resulting nt_accession_index.tsv.gz can be passed to 'fillet metagenome "
            "--nt-index' to skip the per-build metadata pass and speed up repeated DB builds "
            "against the same NT version."
        ),
    )
    bni.add_argument("--blast-db-dir", required=True, dest="blast_db_dir", metavar="DIR",
                     help="Directory containing the NT BLAST+ database volumes.")
    bni.add_argument("--output", "-o", default=None, dest="output", metavar="FILE",
                     help="Output path for the index (default: {blast_db_dir}/nt_accession_index.tsv.gz).")

    # ── batch ──────────────────────────────────────────────────────────────────
    bt = sub.add_parser(
        "batch",
        help="Run the full metagenomics pipeline on a whole sequencing batch from a metadata sheet.",
        formatter_class=__import__("argparse").RawDescriptionHelpFormatter,
        description="""
Run the full metagenomics pipeline (preprocess → screen → build DB → align → classify)
on every library listed in a metadata Excel/TSV/CSV file.

Each library gets a numbered output directory (01_LibraryID/, 02_LibraryID/, ...).
Progress is checkpointed after every sample so the batch can be resumed if interrupted.
A live stats spreadsheet (batch_stats.tsv / .xlsx) is updated after each sample.
At the end, all SQLite viewer databases are merged into batch.merged.viewer.sqlite.

Metadata file format (one row per library):
  Required columns:  library_id, r1_path, platform
  Key columns:       sample_id, r2_path, batch_id, library_type,
                     sample_age_type, material_type, is_control, control_type,
                     age_BP, site_name, group, include

Generate a template with:
  fillet batch --generate-template samples_template.xlsx

Example run:
  fillet batch \\
    --metadata samples.xlsx \\
    --outdir /results/batch1 \\
    --kraken-db /data/krakenuniq_db \\
    --blast-db-dir /data/NCBI_NT/2026-03-31 \\
    --nodes /data/taxonomy/nodes.dmp \\
    --names /data/taxonomy/names.dmp \\
    --threads 32 --resume -y
""",
    )
    bt.add_argument("--metadata", metavar="FILE",
                    help="Metadata Excel (.xlsx), TSV, or CSV file. One row per library.")
    bt.add_argument("--generate-template", metavar="FILE",
                    help="Write an empty metadata template to FILE and exit.")
    bt.add_argument("--outdir", metavar="DIR",
                    help="Batch output directory (required unless --generate-template).")
    bt.add_argument("--threads", type=int, default=16)
    bt.add_argument("--resume", action="store_true", default=True,
                    help="Skip already-complete samples (default: on).")
    bt.add_argument("--no-resume", dest="resume", action="store_false")
    bt.add_argument("--stop-on-error", action="store_true",
                    help="Halt batch on first failure (default: continue and log error).")
    bt.add_argument("--yes", "-y", action="store_true", help="Skip confirmation prompt.")
    bt.add_argument("--no-shared-db", dest="use_shared_db", action="store_false", default=True,
                    help="Build a separate LAST database per library instead of one shared DB per batch.")

    bt_screen = bt.add_argument_group("Stage 1: Taxonomic screening")
    bt_screen.add_argument("--kraken-db", metavar="DIR")
    bt_screen.add_argument("--screener", choices=["auto", "krakenuniq", "kraken2"], default="auto")
    bt_screen.add_argument("--min-kraken-reads", type=int, default=5, metavar="N")
    bt_screen.add_argument("--min-kraken-unique-kmers", type=int, default=50, metavar="N")
    bt_screen.add_argument("--max-kraken-dup-rate", type=float, default=None, metavar="X")

    bt_db = bt.add_argument_group("Stage 2: Subset LAST database")
    bt_db.add_argument("--nodes", required=False, metavar="NODES_DMP",
                       help="NCBI taxonomy nodes.dmp (required).")
    bt_db.add_argument("--blast-db-dir", metavar="DIR")
    bt_db.add_argument("--custom-ref-fasta", nargs="+", metavar="FASTA",
                       help="Custom reference FASTAs (FlyGuide/Spinner) merged into every sample's LAST DB.")
    bt_db.add_argument("--last-db-prefix-base", metavar="DIR",
                       help="Base directory for per-sample LAST DBs. Default: OUTDIR/shared_last_db/")
    bt_db.add_argument("--skip-build-db", action="store_true")
    bt_db.add_argument("--max-seed-rank", default="family", metavar="RANK",
                       help="Maximum rank allowed as a DB expansion seed (default: family). "
                            "Seeds above this rank are dropped before BFS expansion.")
    bt_db.add_argument("--max-seqs-per-taxon", type=int, default=None, dest="max_seqs_per_taxon",
                       metavar="N",
                       help="Cap sequences extracted per taxon from NT (default: no cap). "
                            "Recommended: 300. Prevents OOM when --no-chunked-db is used.")
    bt_db.add_argument("--no-chunked-db", dest="chunked_db", action="store_false", default=True,
                       help="Disable hierarchical chunked DB; use single monolithic LAST DB.")
    bt_db.add_argument("--chunked-db-rank", default=None, dest="chunked_db_rank", metavar="RANK",
                       help="Force partition rank for chunked DB (e.g. family, order, class).")
    bt_db.add_argument("--parallel-chunks", type=int, default=1, dest="parallel_chunks", metavar="N",
                       help="Number of chunk DBs to build/align in parallel (default: 1).")
    bt_db.add_argument("--chunk-timeout", type=int, default=1800, dest="chunk_timeout_sec", metavar="SEC",
                       help="Seconds to wait for rank selection before auto-selecting (default: 1800).")
    bt_db.add_argument("--no-keep-chunk-mafs", dest="keep_chunk_mafs", action="store_false", default=True,
                       help="Delete per-chunk MAF files after merging.")

    bt_align = bt.add_argument_group("Stage 3: LAST alignment")
    bt_align.add_argument("--last-train-file", metavar="FILE",
                          help="LAST training file. Default: bundled nt_adna.train.")
    bt_align.add_argument("--evalue", default="1e-5")
    bt_align.add_argument("--max-target-seqs", type=int, default=2000)
    bt_align.add_argument("--chunk-size", type=int, default=50000, metavar="N")

    bt_clf = bt.add_argument_group("Stage 4: Classification")
    bt_clf.add_argument("--names", metavar="NAMES_DMP")
    bt_clf.add_argument("--taxonomy-tsv", metavar="FILE")
    bt_clf.add_argument("--root-taxid", default="1")
    bt_clf.add_argument("--sample-sheet", metavar="FILE",
                        help="Fillet sample sheet (roles/groups). Defaults to batch metadata roles.")
    bt_clf.add_argument("--regional-taxa", metavar="FILE")
    bt_clf.add_argument("--config", metavar="TOML")
    bt_clf.add_argument("--em-iterations", type=int, default=1)
    bt_clf.add_argument("--max-viewer-reads", type=int, default=5000)
    bt_clf.add_argument("--no-merge", dest="merge_at_end", action="store_false", default=True,
                        help="Skip merging all SQLite outputs at the end.")

    bt_prep = bt.add_argument_group("Stage 0: Preprocessing")
    bt_prep.add_argument("--min-length", type=int, default=24, metavar="N")
    bt_prep.add_argument("--no-dedup", dest="dedup", action="store_false", default=True)
    bt_prep.add_argument("--no-adapter-check", dest="adapter_check", action="store_false", default=True)
    bt_prep.add_argument("--preprocess-tool", choices=["auto", "fastp", "AdapterRemoval"], default="auto")

    return p


def _progress_groups(
    groups: Iterable[Tuple[str, list]],
    sid: str,
    verbose: bool,
    interval: int = 10_000,
) -> Generator[Tuple[str, list], None, None]:
    """Wrap a read-group iterator and print periodic progress to stderr."""
    if not verbose:
        yield from groups
        return
    t0 = time.monotonic()
    for n, item in enumerate(groups, 1):
        yield item
        if n % interval == 0:
            elapsed = time.monotonic() - t0
            rate = n / elapsed if elapsed > 0 else 0.0
            print(f"[Fillet] {sid}: {n:,} reads @ {rate:.0f} reads/s", file=sys.stderr, flush=True)


def infer_sample_id(path: str) -> str:
    name = Path(path).name
    for suffix in [".outfmt6.gz", ".outfmt6", ".blasttab.gz", ".blasttab", ".tab.gz", ".tsv.gz", ".tab", ".tsv", ".gz"]:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return Path(name).stem


def _read_depth_file(path: Optional[str]) -> Dict[str, int]:
    """Read a TSV with sample_id and total_reads columns; returns {sample_id: int}."""
    if not path:
        return {}
    import csv as _csv
    out: Dict[str, int] = {}
    with open(path, "r", encoding="utf-8", newline="") as fh:
        for row in _csv.DictReader(fh, delimiter="\t"):
            sid = row.get("sample_id") or row.get("sample") or row.get("library_id")
            depth = row.get("total_reads") or row.get("reads") or row.get("depth")
            if sid and depth:
                try:
                    out[str(sid).strip()] = int(float(str(depth).strip()))
                except ValueError:
                    pass
    return out


def _sample_roles(sheet: Dict[str, Dict[str, str]]) -> Dict[str, str]:
    out = {}
    for sid, row in sheet.items():
        role = row.get("role") or row.get("type") or row.get("control_type") or row.get("sample_type") or "sample"
        out[sid] = role
    return out


def _load_tax(args: argparse.Namespace) -> Taxonomy:
    return Taxonomy.from_paths(nodes=args.nodes, names=args.names, taxonomy_tsv=args.taxonomy_tsv, root_taxid=args.root_taxid)


def cmd_classify(args: argparse.Namespace) -> int:
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    cfg = load_config(args.config)
    cfg = apply_mode_profile(cfg, getattr(args, "mode", "auto"))
    cfg = apply_config_overrides(cfg, args)
    tax = _load_tax(args)
    acc2taxid = read_acc2taxid(args.acc2taxid)
    sample_sheet = read_sample_sheet(args.sample_sheet)
    sample_roles = _sample_roles(sample_sheet)
    regional = read_regional_taxa(args.regional_taxa)
    contaminants = read_contaminants(getattr(args, "contaminants_file", None))
    palynology_taxa = read_support_table(getattr(args, "palynology_table", None))
    # Per-sample override (injected by batch for site/age-resolved fossil support)
    _fos_override = getattr(args, "__fossil_taxa_override__", None)
    if _fos_override is not None:
        fossil_taxa = None  # batch-wide set not used; per-sample dict takes over
    else:
        fossil_taxa = read_support_table(getattr(args, "fossil_table", None))
    sample_depths = _read_depth_file(getattr(args, "depth_file", None))
    uq_con = None
    uq_path = getattr(args, "uniqueness_index", None)
    if uq_path:
        from .builddb import load_uniqueness_index
        uq_con = load_uniqueness_index(Path(uq_path))
        if args.verbose:
            print(f"[Fillet] Uniqueness index loaded: {uq_path}", file=sys.stderr)
    dirichlet_alpha = float(getattr(args, "dirichlet_alpha", 0.0))
    damage_weight = float(getattr(args, "damage_weight", 0.0))
    if dirichlet_alpha > 0.0 or damage_weight > 0.0:
        from .relic_lca_bayes import BayesianRelicLCAClassifier
        clf = BayesianRelicLCAClassifier(
            tax, cfg, regional_taxa=regional,
            dirichlet_alpha=dirichlet_alpha,
            damage_likelihood_weight=damage_weight,
        )
        if args.verbose:
            print(
                f"[Fillet] Experimental Bayesian mode: "
                f"dirichlet_alpha={dirichlet_alpha}, damage_weight={damage_weight}",
                file=sys.stderr,
            )
    else:
        clf = RelicLCAClassifier(tax, cfg, regional_taxa=regional, uniqueness_index=uq_con)

    inputs = list(args.input)
    sample_ids = args.sample_id or [infer_sample_id(x) for x in inputs]
    if len(sample_ids) != len(inputs):
        raise SystemExit("--sample-id must have same number of values as --input")
    fastas = args.query_fasta or [None] * len(inputs)
    if len(fastas) != len(inputs):
        raise SystemExit("--query-fasta must have same number of values as --input")

    em_iterations = int(getattr(args, "em_iterations", 1))
    assignments: List[ReadAssignment] = []
    for path, sid, fasta in zip(inputs, sample_ids, fastas):
        if args.verbose:
            print(f"[Fillet] Reading {path} as sample {sid}", file=sys.stderr)
        qlens = read_fasta_lengths(fasta) if fasta else {}
        seqs = read_fasta_records(fasta) if fasta else {}
        is_maf = str(path).endswith('.maf') or str(path).endswith('.maf.gz')
        if is_maf:
            hits = iter_last_maf_hits(path, acc2taxid=acc2taxid, query_lengths=qlens, sample_id=sid)
        else:
            hits = iter_blast_hits(path, columns=args.columns, acc2taxid=acc2taxid, query_lengths=qlens, sample_id=sid)
        groups = _progress_groups(group_hits_by_read(hits), sid, args.verbose)
        sample_assignments = clf.classify_groups_with_em(
            groups, sample_id=sid, sequences=seqs, em_iterations=em_iterations
        )
        assignments.extend(sample_assignments)
        if args.verbose:
            n_assigned = sum(1 for a in sample_assignments if a.assigned_taxid != "0")
            n_em = sum(1 for a in sample_assignments if a.em_refined)
            n_dmg = sum(1 for a in sample_assignments if a.damage_score > 0)
            print(f"[Fillet] {sid}: {n_assigned}/{len(sample_assignments)} reads assigned", file=sys.stderr)
            if em_iterations > 0:
                print(f"[Fillet] {sid}: {n_em} reads refined by EM pass", file=sys.stderr)
            if n_dmg:
                dmg_vals = [a.damage_score for a in sample_assignments if a.damage_score > 0]
                print(f"[Fillet] {sid}: {n_dmg} reads have damage metrics (mean={sum(dmg_vals)/len(dmg_vals):.3f})", file=sys.stderr)

    # Apply manual taxon re-routing (Item 12): reassign reads at SRC taxid to DST (or parent).
    reroute_map: Dict[str, str] = {}
    for spec in getattr(args, "reroute_taxids", []):
        if ":" in spec:
            src, dst = spec.split(":", 1)
        else:
            src = spec
            dst = tax.parent(src) or "1"
        reroute_map[src.strip()] = dst.strip()
    if reroute_map:
        for a in assignments:
            if a.assigned_taxid in reroute_map:
                a.assigned_taxid = reroute_map[a.assigned_taxid]

    # Build per-sample fossil dicts — either from batch override or from flat file
    _fos_override = getattr(args, "__fossil_taxa_override__", None)
    if _fos_override is not None:
        _fos_ev = getattr(args, "__fos_evidence_texts_override__", []) or []
        _fossil_taxa_by_sample = {sid: _fos_override for sid in sample_ids}
        _fos_evidence_by_sample = {sid: _fos_ev for sid in sample_ids}
    else:
        _fossil_taxa_by_sample = None
        _fos_evidence_by_sample = None

    summaries = summarize_assignments(
        tax, assignments,
        sample_roles=sample_roles,
        config=cfg,
        contaminants=contaminants,
        regional_taxa=regional,
        palynology_taxa=palynology_taxa,
        fossil_taxa=fossil_taxa,
        sample_depths=sample_depths,
        fossil_taxa_by_sample=_fossil_taxa_by_sample,
        fos_evidence_by_sample=_fos_evidence_by_sample,
    )
    if dirichlet_alpha > 0.0 or damage_weight > 0.0:
        from .relic_lca_bayes import annotate_dirichlet_ci, estimate_taxon_damage_rates
        dmg_rates = {}
        for sid in {a.sample_id for a in assignments}:
            dmg_rates.update(estimate_taxon_damage_rates(assignments, sid))
        summaries = annotate_dirichlet_ci(summaries, dirichlet_alpha, dmg_rates)

    # Within-clade deconvolution (default on; skip with --no-deconvolve).
    deconv_rows: list = []
    if not getattr(args, "no_deconvolve", False):
        from .deconvolution import run_deconvolution, deconvolution_to_rows
        deconv_results = run_deconvolution(
            tax, summaries, assignments,
            nodes_path=getattr(args, 'nodes', None),
            blast_db_dir=getattr(args, 'blast_db_dir', None),
            query_fasta=Path(args.query_fasta[0]) if getattr(args, 'query_fasta', None) else None,
            deconv_db_dir=None,  # auto-derives from query_fasta.parent
            threads=getattr(args, 'threads', 8),
            max_clade_seqs_per_taxon=getattr(args, 'max_clade_seqs_per_taxon', 5000),
            use_mini_db=not getattr(args, 'no_mini_db', False),
            max_clade_db_mb=getattr(args, 'max_clade_db_mb', 50.0),
            mini_db_timeout_sec=int(getattr(args, 'mini_db_timeout_min', 30.0) * 60) if getattr(args, 'mini_db_timeout_min', 30.0) > 0 else None,
            max_parallel_mini_dbs=getattr(args, 'max_parallel_mini_dbs', 1),
            mini_db_index_path=getattr(args, 'nt_index', None),
        )
        if deconv_results:
            deconv_rows = deconvolution_to_rows(deconv_results)
            deconv_path = outdir / f"{args.prefix}.deconvolution.tsv"
            deconv_fields = list(deconv_rows[0].keys()) if deconv_rows else []
            write_tsv(deconv_path, deconv_rows, deconv_fields)
            print(f"[Fillet] Wrote: {deconv_path} ({len(deconv_results)} clade(s))")

    assignment_fields = list(assignments_to_rows(assignments[:1])[0].keys()) if assignments else list(ReadAssignment.__dataclass_fields__.keys())
    summary_fields = list(summaries_to_rows(summaries[:1])[0].keys()) if summaries else list(TaxonSummary.__dataclass_fields__.keys())
    assign_path = outdir / f"{args.prefix}.read_assignments.tsv"
    summary_path = outdir / f"{args.prefix}.taxon_summary.tsv"
    write_tsv(assign_path, assignments_to_rows(assignments), assignment_fields)
    write_tsv(summary_path, summaries_to_rows(summaries), summary_fields)
    write_metamerge_tables(outdir, summaries)
    write_megan_wide(outdir, summaries)
    write_holi_compat(outdir, summaries, taxonomy=tax)
    write_read2tax(outdir, assignments)
    write_krona_like(outdir, assignments, tax)

    if not args.no_viewer:
        payload = build_viewer_payload(tax, assignments, summaries, max_reads=args.max_viewer_reads)
        write_viewer_json(payload, outdir / f"{args.prefix}.viewer.json")
        write_viewer_html(payload, outdir / f"{args.prefix}.viewer.html")
    if getattr(args, "sqlite_viewer", True):
        sqlite_path = outdir / f"{args.prefix}.viewer.sqlite"
        _dark_tsv = outdir / "dark_taxids.tsv"
        write_viewer_sqlite(sqlite_path, tax, assignments, summaries,
                            nodes_dmp_path=getattr(args, "nodes", None),
                            deconvolution_rows=deconv_rows or None,
                            dark_taxids_tsv=_dark_tsv if _dark_tsv.exists() else None)

    print(f"[Fillet] Wrote: {assign_path}")
    print(f"[Fillet] Wrote: {summary_path}")
    if not args.no_viewer:
        print(f"[Fillet] Wrote: {outdir / (args.prefix + '.viewer.html')}")
    if getattr(args, "sqlite_viewer", True):
        print(f"[Fillet] Wrote: {outdir / (args.prefix + '.viewer.sqlite')}")
        print(f"[Fillet] Serve large viewer with: fillet serve-viewer --db {outdir / (args.prefix + '.viewer.sqlite')} --port 8765")
    return 0


def _read_tsv(path: str) -> List[Dict[str, str]]:
    import csv
    with open(path, "r", encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _coerce_assignment_row(r: Dict[str, str]) -> ReadAssignment:
    fields = ReadAssignment.__dataclass_fields__
    out = {k: r.get(k, "") for k in fields}
    for k in ["posterior", "entropy", "best_bitscore", "best_pident", "second_path_posterior",
              "damage_ct_5p", "damage_ga_3p", "damage_score"]:
        out[k] = float(out.get(k) or 0)
    for k in ["n_hits_raw", "n_hits_used", "best_aln_len", "em_refined"]:
        out[k] = int(float(out.get(k) or 0))
    for k in ["best_qcov"]:
        out[k] = None if out.get(k) in {"", "None", None} else float(out.get(k))
    for k in ["read_len", "best_qstart", "best_qend", "best_sstart", "best_send", "best_slen"]:
        out[k] = None if out.get(k) in {"", "None", None} else int(float(out.get(k)))
    return ReadAssignment(**out)


def _coerce_summary_row(r: Dict[str, str]) -> TaxonSummary:
    fields = TaxonSummary.__dataclass_fields__
    out = {k: r.get(k, "") for k in fields}
    for k in ["direct_hard_reads", "cumulative_hard_reads", "conflicted_reads", "evidence_reads", "unique_references", "max_stack_depth"]:
        out[k] = int(float(out.get(k) or 0))
    for k in ["direct_weighted_reads", "cumulative_weighted_reads", "mean_posterior", "negative_weighted_reads", "blank_fraction", "best_reference_breadth", "mean_reference_breadth", "top_locus_fraction", "mean_damage_score", "max_damage_score"]:
        out[k] = float(out.get(k) or 0)
    return TaxonSummary(**out)


def cmd_view(args: argparse.Namespace) -> int:
    tax = _load_tax(args)
    assignments = [_coerce_assignment_row(r) for r in _read_tsv(args.assignments)]
    summaries = [_coerce_summary_row(r) for r in _read_tsv(args.summary)]
    payload = build_viewer_payload(tax, assignments, summaries, max_reads=args.max_viewer_reads)
    write_viewer_html(payload, args.out_html)
    if getattr(args, "out_sqlite", None):
        write_viewer_sqlite(args.out_sqlite, tax, assignments, summaries)
        print(f"[Fillet] Wrote: {args.out_sqlite}")
    print(f"[Fillet] Wrote: {args.out_html}")
    return 0


def cmd_init(args: argparse.Namespace) -> int:
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    from importlib import resources
    cfg_text = resources.files("fillet.data").joinpath("default_config.toml").read_text(encoding="utf-8")
    (out / "fillet_config.toml").write_text(cfg_text, encoding="utf-8")
    (out / "samplesheet.tsv").write_text(
        "sample_id\trole\tgroup\tnotes\n"
        "Sample_01\tsample\tLayer_A\tregular sediment sample\n"
        "Blank_01\tnegative\tcontrols\textraction or library blank\n"
        "Positive_01\tpositive\tcontrols\toptional positive control\n",
        encoding="utf-8",
    )
    (out / "regional_taxa.tsv").write_text(
        "taxid\tscientific_name\tstatus\thabitat\tweight\tnotes\n"
        "9903\tBos taurus\tintroduced_modern\tterrestrial\t0.80\tExample only\n"
        "9901\tBison bison\tpresent_past_or_regional\tterrestrial\t1.35\tExample only\n",
        encoding="utf-8",
    )
    print(f"[Fillet] Initialized skeleton in {out}")
    return 0


def cmd_run_align(args: argparse.Namespace) -> int:
    out = Path(args.out)
    workdir = Path(args.workdir) if getattr(args, 'workdir', None) else out.parent / 'fillet_align_work'
    run_chunked_alignment(
        aligner=args.aligner,
        query=args.query,
        output=out,
        workdir=workdir,
        blast_db=args.db,
        blast_task=getattr(args, 'task', 'blastn'),
        evalue=getattr(args, 'evalue', '1e-5'),
        max_target_seqs=args.max_target_seqs,
        max_hsps=args.max_hsps,
        threads=args.threads,
        chunk_size=args.chunk_size,
        resume=bool(args.resume),
        with_inspection_fields=bool(args.with_inspection_fields),
        last_db_prefix=args.last_db_prefix,
        last_train_file=args.last_train_file,
        last_threads=args.last_threads,
        last_output_format=getattr(args, 'last_output_format', 'tab'),
    )
    return 0


def cmd_run_blast(args: argparse.Namespace) -> int:
    # Back-compatible wrapper. If --with-inspection-fields was not supplied, v0.3 behavior is preserved;
    # workflow/run-align default to inspection fields because the viewer benefits greatly from qseq/sseq/stitle.
    out = Path(args.out)
    workdir = Path(args.workdir) if args.workdir else out.parent / 'fillet_blast_work'
    run_chunked_alignment(
        aligner='blast',
        query=args.query,
        output=out,
        workdir=workdir,
        blast_db=args.db,
        blast_task=args.task,
        evalue=args.evalue,
        max_target_seqs=args.max_target_seqs,
        max_hsps=args.max_hsps,
        threads=args.threads,
        chunk_size=args.chunk_size,
        resume=bool(args.resume),
        with_inspection_fields=bool(args.with_inspection_fields),
    )
    return 0


def _estimate_runtime(n_reads: int, threads: int, with_inspection: bool, db_path: str | None = None) -> str:
    """Return a human-readable BLAST + classify runtime estimate."""
    import os
    db_gb = 0.0
    if db_path:
        for ext in [".nal", ".nin", ".nsq", ""]:
            p = Path(str(db_path) + ext)
            if p.exists():
                db_gb += p.stat().st_size / 1e9
        if db_gb == 0:
            base = Path(db_path)
            for p in base.parent.glob(base.name + "*"):
                db_gb += p.stat().st_size / 1e9
    rate_per_thread = 30 if with_inspection else 60
    blast_lo = n_reads / max(1, rate_per_thread * 1.5 * threads)
    blast_hi = n_reads / max(1, rate_per_thread * 0.5 * threads)
    classify_est = n_reads / 500

    def _hms(s: float) -> str:
        if s < 120:
            return f"{int(s)}s"
        if s < 7200:
            return f"{int(s // 60)}m"
        h = int(s // 3600)
        m = int((s % 3600) // 60)
        return f"{h}h {m}m"

    db_note = f"  DB: ~{db_gb:.0f} GB" if db_gb > 0 else ""
    total_lo = blast_lo + classify_est
    total_hi = blast_hi + classify_est
    return (
        f"FASTA: {n_reads:,} reads  Threads: {threads}{db_note}\n"
        f"  BLAST: {_hms(blast_lo)} – {_hms(blast_hi)}"
        + (" (inspection mode — includes qseq/sseq)" if with_inspection else "") + "\n"
        f"  Classify: ~{_hms(classify_est)}\n"
        f"  TOTAL: ~{_hms(total_lo)} – {_hms(total_hi)}"
    )


def _display_workflow_settings(args: argparse.Namespace, cfg: dict, n_reads: int | None = None) -> None:
    """Print a MEGAN7-style settings summary to stdout before a workflow run."""
    W = 76  # total width including ║ borders
    lca = cfg.get("hit_filters", {})
    asgn = cfg.get("assignment", {})
    # ║  {content}  ║ — 3 chars on each side → inner = W - 6
    inner = W - 6

    def _line(label: str, val: str) -> str:
        avail = inner - len(label) - 1
        val_trunc = val if len(val) <= avail else val[:avail - 1] + "…"
        gap = inner - len(label) - len(val_trunc)
        return f"║  {label}{' ' * max(1, gap)}{val_trunc}  ║"

    def _note(text: str) -> str:
        text_trunc = text if len(text) <= inner else text[:inner - 1] + "…"
        gap = inner - len(text_trunc)
        return f"║  {text_trunc}{' ' * max(0, gap)}  ║"

    border_top = "╔" + "═" * (W - 2) + "╗"
    border_mid = "╠" + "═" * (W - 2) + "╣"
    border_bot = "╚" + "═" * (W - 2) + "╝"

    def _section(title: str) -> str:
        pad = inner - len(title)
        return f"║  {title}{'─' * max(1, pad)}  ║"

    lines = [
        border_top,
        _line("Fillet workflow", f"v{__import__('fillet').__version__}"),
        border_mid,
        _section("Input"),
        _line("Query FASTA", str(getattr(args, "query", ""))),
        _line("Sample ID", str(getattr(args, "sample_id", ""))),
        _line("Output directory", str(getattr(args, "outdir", ""))),
    ]
    if n_reads is not None:
        lines.append(_line("Reads in FASTA", f"{n_reads:,}"))
    lines += [
        border_mid,
        _section("BLAST"),
        _line("Aligner", str(getattr(args, "aligner", "blast"))),
        _line("Database", str(getattr(args, "db", "") or getattr(args, "last_db_prefix", ""))),
        _line("Threads", str(getattr(args, "threads", 16))),
        _line("Max target seqs", str(getattr(args, "max_target_seqs", 2000))),
        _line("E-value", str(getattr(args, "evalue", "1e-5"))),
        _line("Inspection fields", "yes (qseq/sseq — required for damage & alignment viewer)"),
        border_mid,
        _section("RELIC-LCA"),
        _line("Min bitscore", str(lca.get("min_bitscore", 40.0))),
        _line("Top bitscore fraction", str(lca.get("top_bitscore_fraction", 0.80))),
        _note("  └─ 0.80 keeps hits within 80% of best (MEGAN uses 0.95 → reference-bias FP)"),
        _line("Min alignment length", str(lca.get("min_aln_len", 24)) + " bp"),
        _line("Min % identity", str(lca.get("min_pident", 80.0))),
        _line("Min query coverage", str(lca.get("min_qcov", 0.0))),
        _line("EM iterations", str(getattr(args, "em_iterations", 1))),
        _line("Min node posterior", str(asgn.get("min_node_posterior", 0.65))),
        _line("Min posterior margin", str(asgn.get("min_posterior_margin", 0.10))),
        _line("Min posterior (family)", str(asgn.get("min_posterior_family", 0.66))),
        _line("Min posterior (genus)", str(asgn.get("min_posterior_genus", 0.72))),
        _line("Min posterior (species)", str(asgn.get("min_posterior_species", 0.78))),
    ]

    if n_reads is not None:
        with_inspection = True
        est = _estimate_runtime(
            n_reads,
            getattr(args, "threads", 16),
            with_inspection,
            getattr(args, "db", None) or getattr(args, "last_db_prefix", None),
        )
        lines += [border_mid, _section("Runtime estimate")]
        for row in est.splitlines():
            lines.append(_note(row.strip()))

    lines.append(border_bot)
    print("\n".join(lines), flush=True)


def cmd_build_db(args: argparse.Namespace) -> int:
    """Build a targeted single-volume LAST database for a set of taxids."""
    import sys
    import logging
    logging.basicConfig(level=logging.INFO, format="%(message)s", stream=sys.stderr)

    from .builddb import build_subset_db

    # Collect seed taxids
    taxids: list = []
    if args.taxids:
        taxids.extend(args.taxids)
    if getattr(args, "taxid_file", None):
        from pathlib import Path as _P
        for line in _P(args.taxid_file).read_text().splitlines():
            line = line.strip()
            if line and not line.startswith("#"):
                taxids.append(line)
    if not taxids:
        print("Error: provide --taxids or --taxid-file.", file=sys.stderr)
        return 1

    print(
        f"[fillet build-db] Seed taxids: {len(taxids)}  |  "
        f"expand descendants: {not args.no_expand}  |  "
        f"threads: {args.threads}",
        file=sys.stderr,
    )
    print(f"[fillet build-db] BLAST DB dir: {args.blast_db_dir}", file=sys.stderr)
    print(f"[fillet build-db] Output prefix: {args.out_prefix}", file=sys.stderr)

    result = build_subset_db(
        taxids=taxids,
        nodes_path=args.nodes,
        blast_db_dir=args.blast_db_dir,
        output_prefix=args.out_prefix,
        expand_descendants=not args.no_expand,
        threads=args.threads,
        verbose=args.verbose,
    )
    print(
        f"[fillet build-db] Done.  "
        f"taxids: {result['n_taxids']:,}  |  "
        f"sequences: {result['n_sequences']:,}  |  "
        f"LAST DB: {result['output_prefix']}",
        file=sys.stderr,
    )
    return 0


def _parse_iso_ts(s: str) -> Optional[float]:
    """Parse ISO 8601 timestamp to epoch seconds, or None."""
    if not s:
        return None
    try:
        from datetime import datetime
        return datetime.fromisoformat(s).timestamp()
    except Exception:
        return None


def _fmt_elapsed(started: str, ended: str) -> str:
    """Return human-readable elapsed time between two ISO timestamps."""
    t0 = _parse_iso_ts(started)
    t1 = _parse_iso_ts(ended)
    if t0 is None or t1 is None:
        return ""
    secs = max(0, t1 - t0)
    h, rem = divmod(int(secs), 3600)
    m, s = divmod(rem, 60)
    if h:
        return f"{h}h{m:02d}m"
    return f"{m}m{s:02d}s"


def cmd_batch_status(args: argparse.Namespace) -> int:
    """Show batch progress summary without re-running."""
    import json as _json
    outdir = Path(args.outdir)
    progress_file = outdir / "batch_progress.json"
    stats_file = outdir / "batch_stats.tsv"

    if not progress_file.exists():
        print(f"No batch_progress.json found in {outdir}", file=sys.stderr)
        return 1

    state = _json.loads(progress_file.read_text(encoding="utf-8"))
    counts: Dict[str, int] = {"pending": 0, "running": 0, "complete": 0, "failed": 0, "screened": 0}
    rows = []
    for lib_id, info in state.items():
        status = info.get("status", "pending")
        counts[status] = counts.get(status, 0) + 1
        started = info.get("started_at", "")
        ended = info.get("completed_at") or info.get("failed_at", "")
        elapsed = _fmt_elapsed(started, ended)
        error = (info.get("error", "") or "")[:60]
        sqlite = info.get("sqlite", "")
        if sqlite:
            sqlite = Path(sqlite).name  # show filename only
        rows.append((lib_id, status, started[:16], ended[:16], elapsed, error, sqlite))

    total = len(rows)
    n_complete = counts.get("complete", 0)
    n_failed = counts.get("failed", 0)
    n_running = counts.get("running", 0)
    n_pending = counts.get("pending", 0)

    print(f"\n  Batch: {outdir}")
    print(f"  Total: {total}  |  Complete: {n_complete}  |  Failed: {n_failed}"
          f"  |  Running: {n_running}  |  Pending: {n_pending}")
    print()

    w = max((len(r[0]) for r in rows), default=20)
    marks = {"complete": "✓", "failed": "✗", "running": "→", "pending": "·", "screened": "~"}
    print(f"  {'':1}  {'Library':<{w}}  {'Status':<10}  {'Started':<16}  {'Ended':<16}  {'Elapsed':>8}  Detail")
    print("  " + "─" * (w + 70))
    for lib_id, status, started, ended, elapsed, error, sqlite in sorted(rows, key=lambda r: r[0]):
        mark = marks.get(status, " ")
        detail = error if error else sqlite
        print(f"  {mark}  {lib_id:<{w}}  {status:<10}  {started:<16}  {ended:<16}  {elapsed:>8}  {detail}")

    print()

    # Show paths to error logs for failed samples
    failed_libs = [(lib_id, state[lib_id]) for lib_id, *_ in rows if state[lib_id].get("status") == "failed"]
    if failed_libs:
        print("  Failed samples — error logs:")
        for lib_id, info in failed_libs:
            sample_dir = info.get("outdir", "")
            if sample_dir:
                err_log = Path(sample_dir) / "error.log"
                print(f"    {lib_id}: {err_log}")
        print()

    if stats_file.exists():
        print(f"  Stats spreadsheet: {stats_file}")
    runtime_file = outdir / "pipeline_runtime.json"
    if runtime_file.exists():
        rt = _json.loads(runtime_file.read_text(encoding="utf-8"))
        total_s = rt.get("cumulative_elapsed") or rt.get("total_elapsed", 0)
        from .progress import hms
        print(f"  Total runtime: {hms(total_s)}")
    return 0


def cmd_batch(args: argparse.Namespace) -> int:
    """Run the full metagenomics pipeline on a sequencing batch."""
    # Handle --generate-template
    if getattr(args, "generate_template", None):
        from .batchmeta import write_metadata_template
        path = args.generate_template
        write_metadata_template(path)
        print(f"  Template written to: {path}")
        print(f"  Fill in your library details, then run:")
        print(f"    fillet batch --metadata {path} --outdir /results/mybatch ...")
        return 0

    if not getattr(args, "metadata", None):
        print("Error: --metadata is required", file=sys.stderr)
        return 1
    if not getattr(args, "outdir", None):
        print("Error: --outdir is required", file=sys.stderr)
        return 1
    if not getattr(args, "nodes", None):
        print("Error: --nodes (NCBI nodes.dmp) is required", file=sys.stderr)
        return 1

    if not getattr(args, "yes", False):
        from .batchmeta import load_batch_metadata
        try:
            samples = load_batch_metadata(args.metadata)
            n = len(samples)
        except Exception:
            n = "?"
        print()
        print("  Fillet Batch Metagenomics Pipeline")
        print("  " + "─" * 60)
        print(f"  Metadata:   {args.metadata}")
        print(f"  Libraries:  {n}")
        print(f"  Output:     {args.outdir}")
        print(f"  Screener:   {args.screener}  (DB: {getattr(args, 'kraken_db', None) or 'N/A'})")
        print(f"  BLAST dir:  {getattr(args, 'blast_db_dir', None) or 'N/A'}")
        print(f"  Threads:    {args.threads}")
        print(f"  Resume:     {args.resume}")
        print()
        try:
            ans = input("  Proceed? [Y/n]: ").strip().lower()
        except EOFError:
            ans = "y"
        if ans and ans not in ("y", "yes"):
            print("  Aborted.")
            return 1

    from .batch import run_batch_pipeline
    run_batch_pipeline(
        metadata_path=args.metadata,
        outdir=args.outdir,
        kraken_db=getattr(args, "kraken_db", None),
        screener=args.screener,
        min_kraken_reads=args.min_kraken_reads,
        min_kraken_unique_kmers=args.min_kraken_unique_kmers,
        nodes_path=args.nodes,
        blast_db_dir=getattr(args, "blast_db_dir", None),
        last_db_prefix_base=getattr(args, "last_db_prefix_base", None),
        skip_build_if_exists=bool(getattr(args, "skip_build_db", False)),
        custom_ref_fastas=getattr(args, "custom_ref_fasta", None) or [],
        max_seed_rank=getattr(args, "max_seed_rank", "family"),
        max_seqs_per_taxon=getattr(args, "max_seqs_per_taxon", None) or None,
        chunked_db=getattr(args, "chunked_db", True),
        chunked_db_rank=getattr(args, "chunked_db_rank", None),
        parallel_chunks=getattr(args, "parallel_chunks", 1),
        chunk_timeout_sec=getattr(args, "chunk_timeout_sec", 1800),
        keep_chunk_mafs=getattr(args, "keep_chunk_mafs", True),
        last_train_file=getattr(args, "last_train_file", None),
        evalue=args.evalue,
        max_target_seqs=args.max_target_seqs,
        chunk_size=args.chunk_size,
        threads=args.threads,
        taxonomy_tsv=getattr(args, "taxonomy_tsv", None),
        names_path=getattr(args, "names", None),
        root_taxid=getattr(args, "root_taxid", "1"),
        sample_sheet=getattr(args, "sample_sheet", None),
        regional_taxa=getattr(args, "regional_taxa", None),
        config=getattr(args, "config", None),
        em_iterations=getattr(args, "em_iterations", 1),
        max_viewer_reads=getattr(args, "max_viewer_reads", 5000),
        min_length=getattr(args, "min_length", 24),
        dedup=getattr(args, "dedup", True),
        adapter_check=getattr(args, "adapter_check", True),
        preprocess_tool=getattr(args, "preprocess_tool", "auto"),
        resume=getattr(args, "resume", True),
        merge_at_end=getattr(args, "merge_at_end", True),
        stop_on_error=getattr(args, "stop_on_error", False),
        use_shared_db=getattr(args, "use_shared_db", True),
        fossil_table=getattr(args, "fossil_table", None),
        verbose=True,
    )
    return 0


def cmd_metagenome(args: argparse.Namespace) -> int:
    """Stage-by-stage metagenomics pipeline with nice progress display."""
    from .metagenome import run_metagenome_pipeline
    from .batchmeta import is_two_color_platform

    # Validate input: need exactly one of --query, --fastq-r1, or --bam
    # (unless --from-stage is used, in which case --query may be auto-discovered)
    has_query = bool(getattr(args, "query", None))
    has_fastq = bool(getattr(args, "fastq_r1", None))
    has_bam = bool(getattr(args, "bam", None))
    has_from_stage = bool(getattr(args, "from_stage", None))
    if sum([has_query, has_fastq, has_bam]) == 0 and not has_from_stage:
        print("Error: provide one of --query, --fastq-r1, or --bam", file=sys.stderr)
        return 1
    if sum([has_query, has_fastq, has_bam]) > 1:
        print("Error: --query, --fastq-r1, and --bam are mutually exclusive", file=sys.stderr)
        return 1

    # Resolve poly-G trim: explicit flag OR auto-detected from platform
    poly_g = getattr(args, "poly_g_trim", False)
    platform = getattr(args, "platform", None)
    if platform and not poly_g:
        poly_g = is_two_color_platform(platform)

    if not getattr(args, "yes", False):
        print()
        print("  Fillet metagenome pipeline")
        print("  " + "─" * 60)
        if has_bam:
            print(f"  Input BAM:  {args.bam}")
        elif has_fastq:
            print(f"  Input R1:   {args.fastq_r1}")
            if getattr(args, "fastq_r2", None):
                print(f"  Input R2:   {args.fastq_r2}")
        else:
            print(f"  Query:      {args.query}")
        print(f"  Sample:     {args.sample_id}")
        print(f"  Output:     {args.outdir}")
        if platform:
            poly_g_label = " (2-colour → poly-G trim ON)" if poly_g else ""
            print(f"  Platform:   {platform}{poly_g_label}")
        elif poly_g:
            print(f"  Poly-G trim: enabled")
        if getattr(args, "from_stage", None):
            print(f"  From stage: {args.from_stage}  (earlier stages skipped)")
            if getattr(args, "maf", None):
                print(f"  MAF:        {args.maf}")
        print(f"  Screener:   {args.screener}  (DB: {getattr(args, 'kraken_db', None) or 'N/A'})")
        print(f"  BLAST dir:  {getattr(args, 'blast_db_dir', None) or 'N/A'}")
        _train = getattr(args, 'last_train_file', None) or "(bundled nt_adna.train)"
        print(f"  LAST train: {_train}")
        print(f"  Threads:    {args.threads}")
        print()
        try:
            ans = input("  Proceed? [Y/n]: ").strip().lower()
        except EOFError:
            ans = "y"
        if ans and ans not in ("y", "yes"):
            print("  Aborted.")
            return 1

    run_metagenome_pipeline(
        fastq_r1=getattr(args, "fastq_r1", None),
        fastq_r2=getattr(args, "fastq_r2", None),
        bam=getattr(args, "bam", None),
        preprocess_tool=getattr(args, "preprocess_tool", "auto"),
        min_length=getattr(args, "min_length", 24),
        dedup=getattr(args, "dedup", True),
        adapter_check=getattr(args, "adapter_check", True),
        poly_g_trim=poly_g,
        query=getattr(args, "query", None),
        sample_id=args.sample_id,
        outdir=args.outdir,
        kraken_db=getattr(args, "kraken_db", None),
        screener=args.screener,
        min_kraken_reads=args.min_kraken_reads,
        min_kraken_unique_kmers=args.min_kraken_unique_kmers,
        max_kraken_dup_rate=getattr(args, "max_kraken_dup_rate", None),
        taxid_file=getattr(args, "taxid_file", None),
        nodes_path=args.nodes,
        blast_db_dir=getattr(args, "blast_db_dir", None),
        custom_ref_fastas=getattr(args, "custom_ref_fasta", None) or [],
        last_db_prefix=getattr(args, "last_db_prefix", None),
        skip_build_if_exists=bool(getattr(args, "skip_build_db", False)),
        max_seed_rank=getattr(args, "max_seed_rank", "family"),
        last_train_file=getattr(args, "last_train_file", None),
        last_min_score=getattr(args, "last_min_score", None),
        last_m=getattr(args, "last_m", None),
        evalue=args.evalue,
        max_target_seqs=args.max_target_seqs,
        max_hsps=getattr(args, "max_hsps", 1),
        chunk_size=args.chunk_size,
        resume=bool(args.resume),
        threads=args.threads,
        last_threads=getattr(args, "last_threads", None),
        taxonomy_tsv=getattr(args, "taxonomy_tsv", None),
        names_path=getattr(args, "names", None),
        root_taxid=getattr(args, "root_taxid", "1"),
        sample_sheet=getattr(args, "sample_sheet", None),
        regional_taxa=getattr(args, "regional_taxa", None),
        config=getattr(args, "config", None),
        em_iterations=getattr(args, "em_iterations", 1),
        dirichlet_alpha=getattr(args, "dirichlet_alpha", 0.0),
        damage_weight=getattr(args, "damage_weight", 0.0),
        no_deconvolve=getattr(args, "no_deconvolve", False),
        use_mini_db=not getattr(args, "no_mini_db", False),
        max_clade_seqs_per_taxon=getattr(args, "max_clade_seqs_per_taxon", 5000),
        max_clade_db_mb=getattr(args, "max_clade_db_mb", 50.0),
        mini_db_timeout_min=getattr(args, "mini_db_timeout_min", 30.0),
        max_parallel_mini_dbs=getattr(args, "max_parallel_mini_dbs", 1),
        curate_refs=not getattr(args, 'no_curate_refs', False),
        min_ref_length=getattr(args, 'min_ref_length', 300),
        max_ref_n_fraction=getattr(args, 'max_ref_n_fraction', 0.5),
        ref_scope=getattr(args, 'ref_scope', 'all'),
        ref_mode=getattr(args, 'ref_mode', 'full'),
        nt_index_path=getattr(args, 'nt_index', None),
        custom_ref_taxid_map=getattr(args, 'custom_ref_taxid_map', None),
        dedup_similarity=getattr(args, 'dedup_similarity', 0.90),
        max_seqs_per_taxon=getattr(args, "max_seqs_per_taxon", None) or None,
        chunked_db=getattr(args, "chunked_db", True),
        chunked_db_rank=getattr(args, "chunked_db_rank", None),
        parallel_chunks=getattr(args, "parallel_chunks", 1),
        chunk_timeout_sec=getattr(args, "chunk_timeout_sec", 1800),
        keep_chunk_mafs=getattr(args, "keep_chunk_mafs", True),
        yes=getattr(args, "yes", False),
        max_viewer_reads=getattr(args, "max_viewer_reads", 5000),
        no_viewer=getattr(args, "no_viewer", False),
        sqlite_viewer=getattr(args, "sqlite_viewer", True),
        from_stage=getattr(args, "from_stage", None),
        maf_path=getattr(args, "maf", None),
        verbose=getattr(args, "verbose", True),
        lca_overrides={
            **{key: getattr(args, key, None) for _f, _s, key, _c in QUICK_LCA_FLAGS},
            "suppress_unclassified_nodes": getattr(args, "suppress_unclassified_nodes", False),
            "cross_clade_graduated_off": getattr(args, "cross_clade_graduated_off", False),
        },
        mode=getattr(args, "mode", "auto"),
        uniqueness_index=getattr(args, "uniqueness_index", None),
    )
    return 0


def cmd_workflow(args: argparse.Namespace) -> int:
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    prefix = args.prefix or args.sample_id

    # Count reads and show settings summary before doing anything expensive
    cfg = load_config(args.config)
    cfg = apply_mode_profile(cfg, getattr(args, "mode", "auto"))
    cfg = apply_config_overrides(cfg, args)
    n_reads: int | None = None
    try:
        from .align import count_fasta_records
        n_reads = count_fasta_records(args.query)
    except Exception:
        pass  # non-fatal; estimate skipped

    if not getattr(args, "yes", False):
        _display_workflow_settings(args, cfg, n_reads)
        try:
            ans = input("\nProceed with these settings? [Y/n]: ").strip().lower()
        except EOFError:
            ans = "y"
        if ans and ans not in ("y", "yes"):
            print("[Fillet] Aborted by user.")
            return 1

    last_fmt = getattr(args, 'last_output_format', 'tab')
    aligned_ext = '.maf' if (args.aligner == 'last' and last_fmt == 'maf') else '.outfmt6'
    aligned = outdir / f"{prefix}.{args.aligner}{aligned_ext}"
    workdir = outdir / f"{prefix}.align_work"
    run_chunked_alignment(
        aligner=args.aligner,
        query=args.query,
        output=aligned,
        workdir=workdir,
        blast_db=args.db,
        blast_task=args.task,
        evalue=args.evalue,
        max_target_seqs=args.max_target_seqs,
        max_hsps=args.max_hsps,
        threads=args.threads,
        chunk_size=args.chunk_size,
        resume=bool(args.resume),
        with_inspection_fields=True,
        last_db_prefix=args.last_db_prefix,
        last_train_file=args.last_train_file,
        last_threads=args.last_threads,
        last_output_format=last_fmt,
    )
    if args.aligner == 'last' and last_fmt == 'tab' and not args.acc2taxid:
        print('[Fillet workflow] WARNING: LAST blasttab usually lacks taxids. Provide --acc2taxid or few/no hits may classify.', file=sys.stderr)
    classify_args = argparse.Namespace(
        input=[str(aligned)],
        sample_id=[args.sample_id],
        columns=None,
        query_fasta=[args.query],
        acc2taxid=args.acc2taxid,
        nodes=args.nodes,
        names=args.names,
        taxonomy_tsv=args.taxonomy_tsv,
        root_taxid=args.root_taxid,
        sample_sheet=args.sample_sheet,
        regional_taxa=args.regional_taxa,
        config=args.config,
        outdir=str(outdir),
        prefix=prefix,
        max_viewer_reads=args.max_viewer_reads,
        no_viewer=args.no_viewer,
        sqlite_viewer=args.sqlite_viewer,
        verbose=True if args.verbose else False,
        em_iterations=getattr(args, "em_iterations", 1),
        uniqueness_index=getattr(args, "uniqueness_index", None),
    )
    # transfer quick LCA override attributes
    for _flag, _section, key, _caster in QUICK_LCA_FLAGS:
        setattr(classify_args, key, getattr(args, key, None))
    return cmd_classify(classify_args)


def cmd_plot(args: argparse.Namespace) -> int:
    """Generate a publication-quality plot from taxon_summary.tsv files."""
    import matplotlib
    matplotlib.use("Agg")  # headless backend

    from .plot import (
        load_taxon_summary_tsv, bar_chart, stacked_bar_chart, bubble_plot,
        heatmap, stratigraphic_plot, stratigraphic_with_climate,
        stratigraphic_with_paleoclimate,
        damage_distribution, load_climate_series, save_figure,
    )

    rows = []
    for path in args.input:
        rows.extend(load_taxon_summary_tsv(path))

    if not rows:
        print("Error: no data loaded from input files.", file=sys.stderr)
        return 1

    figsize = (args.width or 12, args.height or 6)
    title = args.title or ""
    kwargs: Dict[str, Any] = dict(
        metric=args.metric,
        rank_filter=args.rank,
        min_reads=args.min_reads,
        top_n=args.top_n,
        taxa=args.taxa,
        samples=args.samples,
        title=title,
        figsize=figsize,
        normalize=args.normalize,
    )

    plot_type = args.plot_type
    if plot_type == "bar":
        fig = bar_chart(rows, **kwargs)
    elif plot_type == "stacked":
        fig = stacked_bar_chart(rows, **kwargs)
    elif plot_type == "bubble":
        del kwargs["normalize"]
        fig = bubble_plot(rows, **kwargs)
    elif plot_type == "heatmap":
        del kwargs["normalize"]
        fig = heatmap(rows, log_scale=args.log_scale, **kwargs)
    elif plot_type == "strat":
        fig = stratigraphic_plot(rows, **kwargs)
    elif plot_type == "strat-climate":
        show_lr04 = getattr(args, "show_lr04", False)
        show_gisp2_b = getattr(args, "show_gisp2_b", False)
        show_gisp2_d = getattr(args, "show_gisp2_d", False)
        show_mis = getattr(args, "show_mis", False)
        use_age = getattr(args, "use_age", False)
        dates_14c_path = getattr(args, "dates_14c", None)
        radiocarbon_dates = None
        if dates_14c_path:
            from .paleoclimate import load_14c_dates
            radiocarbon_dates = load_14c_dates(dates_14c_path)
        if show_lr04 or show_gisp2_b or show_gisp2_d or show_mis or radiocarbon_dates:
            # New path: bundled overlays
            del kwargs["normalize"]  # stratigraphic_with_paleoclimate handles it
            fig = stratigraphic_with_paleoclimate(
                rows,
                show_lr04=show_lr04,
                show_gisp2_b=show_gisp2_b,
                show_gisp2_d=show_gisp2_d,
                show_mis=show_mis,
                use_age=use_age,
                normalize=args.normalize,
                radiocarbon_dates=radiocarbon_dates,
                **kwargs,
            )
        else:
            # Legacy path: user-supplied CSV
            climate_rows = None
            if getattr(args, "climate_file", None):
                climate_rows = load_climate_series(args.climate_file)
            fig = stratigraphic_with_climate(rows, climate_rows,
                                             climate_label=getattr(args, "climate_label", "Climate proxy"),
                                             **kwargs)
    elif plot_type == "damage":
        fig = damage_distribution(rows, top_n=args.top_n, taxa=args.taxa,
                                  title=title or "Damage score distribution",
                                  figsize=figsize)
    else:
        print(f"Unknown plot type: {plot_type}", file=sys.stderr)
        return 1

    save_figure(fig, args.out, dpi=args.dpi)
    print(f"Saved {args.plot_type} plot to {args.out}")
    return 0


def cmd_update_taxonomy(args: argparse.Namespace) -> int:
    """Download fresh NCBI taxdump and extract nodes.dmp / names.dmp."""
    import tarfile
    import tempfile
    import urllib.request
    from datetime import date

    outdir = Path(args.outdir)
    if not args.no_date_subdir:
        outdir = outdir / date.today().isoformat()
    outdir.mkdir(parents=True, exist_ok=True)

    url = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
    tarball = outdir / "taxdump.tar.gz"

    print(f"  Downloading {url}")
    print(f"  → {tarball}")

    def _progress(count: int, block_size: int, total_size: int) -> None:
        if total_size > 0:
            pct = min(100, count * block_size * 100 // total_size)
            mb = count * block_size / 1_048_576
            total_mb = total_size / 1_048_576
            sys.stdout.write(f"\r  {pct:3d}%  {mb:.0f} / {total_mb:.0f} MB")
            sys.stdout.flush()

    urllib.request.urlretrieve(url, tarball, reporthook=_progress)
    print()

    want = {"nodes.dmp", "names.dmp", "merged.dmp", "delnodes.dmp"}
    print(f"  Extracting to {outdir} ...")
    with tarfile.open(tarball, "r:gz") as tf:
        for member in tf.getmembers():
            if member.name in want:
                tf.extract(member, path=outdir)
                print(f"    {member.name}  ({member.size / 1_048_576:.1f} MB)")

    if not args.keep_tarball:
        tarball.unlink()

    nodes = outdir / "nodes.dmp"
    names = outdir / "names.dmp"
    print(f"\n  Done.  nodes: {nodes}\n         names: {names}")
    print(f"\n  Use with:  --nodes {nodes} --names {names}")
    return 0


def cmd_batch_metagenome(args: argparse.Namespace) -> int:
    """Run batch metagenomics pipeline: shared DB per sample group."""
    import csv
    from .metagenome import run_batch_metagenome_pipeline

    # Parse sample sheet TSV
    samples = []
    with open(args.samples, "r", encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            row = {k.strip().lower(): v.strip() for k, v in row.items()}
            if not row.get("sample_id"):
                continue
            s = {
                "sample_id": row["sample_id"],
                "outdir": row["outdir"],
                "group_id": row.get("group_id", "default") or "default",
                "query": row.get("query") or None,
                "fastq_r1": row.get("fastq_r1") or row.get("r1") or None,
                "fastq_r2": row.get("fastq_r2") or row.get("r2") or None,
            }
            # Collect any per-sample LCA override columns present in the TSV row
            lca_cols = {key for _f, _s, key, _c in QUICK_LCA_FLAGS}
            per_sample_lca = {k: v for k, v in row.items() if k in lca_cols and v}
            if per_sample_lca:
                s["_lca_overrides"] = per_sample_lca
            samples.append(s)

    if not samples:
        print("Error: no samples found in sample sheet", file=sys.stderr)
        return 1

    groups = {}
    for s in samples:
        groups.setdefault(s["group_id"], []).append(s["sample_id"])
    print(f"  Batch metagenome: {len(samples)} samples, {len(groups)} DB group(s)")
    for gid, sids in groups.items():
        print(f"    group '{gid}': {', '.join(sids)}")
    print()

    if not getattr(args, "yes", False):
        try:
            ans = input("  Proceed? [Y/n]: ").strip().lower()
        except EOFError:
            ans = "y"
        if ans and ans not in ("y", "yes"):
            print("  Aborted.")
            return 1

    global_lca_overrides = {
        **{key: getattr(args, key, None) for _f, _s, key, _c in QUICK_LCA_FLAGS},
        "suppress_unclassified_nodes": getattr(args, "suppress_unclassified_nodes", False),
        "cross_clade_graduated_off": getattr(args, "cross_clade_graduated_off", False),
    }
    run_batch_metagenome_pipeline(
        samples,
        shared_db_outdir=args.shared_db_dir,
        kraken_db=getattr(args, "kraken_db", None),
        screener=args.screener,
        min_kraken_reads=args.min_kraken_reads,
        min_kraken_unique_kmers=args.min_kraken_unique_kmers,
        max_kraken_dup_rate=getattr(args, "max_kraken_dup_rate", None),
        nodes_path=args.nodes,
        names_path=getattr(args, "names", None),
        blast_db_dir=getattr(args, "blast_db_dir", None),
        max_seed_rank=getattr(args, "max_seed_rank", "family"),
        max_seqs_per_taxon=getattr(args, "max_seqs_per_taxon", None) or None,
        chunked_db=getattr(args, "chunked_db", True),
        chunked_db_rank=getattr(args, "chunked_db_rank", None),
        parallel_chunks=getattr(args, "parallel_chunks", 1),
        last_train_file=getattr(args, "last_train_file", None),
        last_m=getattr(args, "last_m", None),
        last_min_score=getattr(args, "last_min_score", None),
        evalue=args.evalue,
        max_target_seqs=args.max_target_seqs,
        chunk_size=args.chunk_size,
        resume=bool(args.resume),
        threads=args.threads,
        config=getattr(args, "config", None),
        em_iterations=getattr(args, "em_iterations", 1),
        no_deconvolve=getattr(args, "no_deconvolve", False),
        use_mini_db=not getattr(args, "no_mini_db", False),
        max_clade_seqs_per_taxon=getattr(args, "max_clade_seqs_per_taxon", 5000),
        max_clade_db_mb=getattr(args, "max_clade_db_mb", 50.0),
        mini_db_timeout_min=getattr(args, "mini_db_timeout_min", 30.0),
        max_parallel_mini_dbs=getattr(args, "max_parallel_mini_dbs", 1),
        curate_refs=not getattr(args, 'no_curate_refs', False),
        min_ref_length=getattr(args, 'min_ref_length', 300),
        max_ref_n_fraction=getattr(args, 'max_ref_n_fraction', 0.5),
        ref_scope=getattr(args, 'ref_scope', 'all'),
        ref_mode=getattr(args, 'ref_mode', 'full'),
        nt_index_path=getattr(args, 'nt_index', None),
        custom_ref_taxid_map=getattr(args, 'custom_ref_taxid_map', None),
        dedup_similarity=getattr(args, 'dedup_similarity', 0.90),
        no_viewer=getattr(args, "no_viewer", False),
        sqlite_viewer=getattr(args, "sqlite_viewer", True),
        lca_overrides=global_lca_overrides,
        uniqueness_index=getattr(args, "uniqueness_index", None),
        yes=getattr(args, "yes", True),
        verbose=getattr(args, "verbose", True),
    )
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.cmd == "classify":
        return cmd_classify(args)
    if args.cmd == "view":
        return cmd_view(args)
    if args.cmd == "init":
        return cmd_init(args)
    if args.cmd == "run-align":
        return cmd_run_align(args)
    if args.cmd == "run-blast":
        return cmd_run_blast(args)
    if args.cmd == "workflow":
        return cmd_workflow(args)
    if args.cmd == "metagenome":
        return cmd_metagenome(args)
    if args.cmd == "batch-metagenome":
        return cmd_batch_metagenome(args)
    if args.cmd == "batch":
        return cmd_batch(args)
    if args.cmd == "serve-viewer":
        serve_viewer(args.db, host=args.host, port=args.port)
        return 0
    if args.cmd == "gui":
        from .gui.app import main as gui_main
        return gui_main(db_path=getattr(args, "db", None))
    if args.cmd == "build-db":
        return cmd_build_db(args)
    if args.cmd == "batch-status":
        return cmd_batch_status(args)
    if args.cmd == "plot":
        return cmd_plot(args)
    if args.cmd == "update-taxonomy":
        return cmd_update_taxonomy(args)
    if args.cmd == "build-nt-index":
        from .builddb import build_nt_accession_index
        out = build_nt_accession_index(
            blast_db_dir=args.blast_db_dir,
            output_path=args.output,
            verbose=True,
        )
        print(f"NT accession index written: {out}", file=__import__("sys").stderr)
        return 0
    if args.cmd == "merge":
        from .viewer_sqlite import merge_viewer_dbs
        if len(args.input) < 2:
            print("Error: --input requires at least two databases.", file=__import__("sys").stderr)
            return 1
        print(f"Merging {len(args.input)} databases → {args.output}", file=__import__("sys").stderr)
        result = merge_viewer_dbs(
            args.input, args.output,
            select_samples=getattr(args, "select_samples", None),
            exclude_samples=getattr(args, "exclude_samples", None),
        )
        print(
            f"Done. {result['n_sources']} sources · "
            f"{len(result['samples'])} samples: {', '.join(result['samples'])}",
            file=__import__("sys").stderr,
        )
        if result["conflicts"]:
            print(
                f"Warning: {len(result['conflicts'])} duplicate sample(s) skipped "
                f"(first occurrence kept): {', '.join(result['conflicts'])}",
                file=__import__("sys").stderr,
            )
        return 0
    parser.print_help()
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
