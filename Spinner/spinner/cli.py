"""Command-line interface for Spinner.

Subcommands
-----------
audit           Annotate records and write decisions TSV (no FASTA output).
filter          Audit and write keep / review / reject FASTAs.
screen-taxonomy Run audit with taxonomy BLAST enabled.
screen-vector   Run audit with vector/UniVec BLAST enabled.
screen-windowed Run audit with windowed BLAST chimerism screen enabled.
cluster         Run audit with vsearch clustering enabled.
explain         Explain why a specific accession was kept / reviewed / rejected.
report          Re-generate summary TSV and HTML from an existing decisions TSV.
init-config     Copy starter configs to a new directory.
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

from . import VERSION
from .pipeline import run_pipeline
from .reporting import report_from_decisions


# ---------------------------------------------------------------------------
# Shared argument helpers
# ---------------------------------------------------------------------------

def _add_common(sp: argparse.ArgumentParser) -> None:
    """Add arguments shared by all pipeline subcommands."""
    sp.add_argument(
        "--fasta", nargs="+", required=True, metavar="FASTA",
        help="Input FASTA file(s), plain or gzipped. Multiple files are merged before processing.",
    )
    sp.add_argument(
        "--outprefix", required=True, metavar="PREFIX",
        help=(
            "Output file prefix, including directory path if desired. "
            "All output files are written as PREFIX.decisions.tsv, PREFIX.keep.fasta, etc."
        ),
    )
    sp.add_argument(
        "--mode", default="", choices=["reference_db", "bait_panel"], metavar="MODE",
        help=(
            "Run profile shortcut. "
            "'reference_db': full-stack curation for read mapping / metabarcoding / MEGAN databases "
            "(taxonomy BLAST + windowed chimerism + vsearch uchime + clustering, moderate QC). "
            "'bait_panel': full-stack curation for capture enrichment bait synthesis via FlyForge "
            "(same pipeline, stricter QC thresholds, chimeric sequences hard-rejected not reviewed). "
            "Auto-selects the matching bundled config. "
            "Providing --config as well overrides --mode."
        ),
    )
    sp.add_argument(
        "--config", default="", metavar="YAML",
        help=(
            "Path to a Spinner YAML config file. Overrides --mode if both are provided. "
            "Only the keys you specify are changed — all other settings use built-in defaults. "
            "Use 'Spinner init-config --outdir DIR' to copy starter configs for editing."
        ),
    )
    sp.add_argument(
        "--species-kingdom", default="", metavar="TSV",
        help=(
            "Optional TSV mapping species to kingdoms, produced by FlyGuide as "
            "*_species_kingdom.tsv. Enables kingdom-level taxonomy cross-checking "
            "(taxonomy_cross_kingdom hard-reject) and the kingdom breakdown in the HTML report. "
            "Without it, all records get kingdom=Unknown and cross-kingdom checks are disabled."
        ),
    )
    sp.add_argument(
        "--regions-config", default="", metavar="TSV",
        help=(
            "TSV with regex rules for marker/region classification (Mito, Plastid, NucMark, Other). "
            "Controls per-class caps and class-level capping. "
            "Defaults to the bundled configs/regions_config_example.tsv when omitted."
        ),
    )
    sp.add_argument(
        "--adapters", default="", metavar="TSV",
        help=(
            "TSV of adapter sequences to screen for (columns: name, sequence, max_mismatch, action). "
            "Covers Illumina, 454/Roche, Ion Torrent, BGI/MGI, and PacBio by default. "
            "Defaults to the bundled configs/adapters_default.tsv when omitted."
        ),
    )
    sp.add_argument(
        "--bad-keywords", default="", metavar="TSV",
        help=(
            "TSV of FASTA header keywords that flag problematic records "
            "(columns: keyword, action, reason). "
            "Defaults to the bundled configs/bad_keywords_ancient.tsv when omitted "
            "(flags synthetic constructs and vectors; reviews UNVERIFIED, PREDICTED, clone, plasmid)."
        ),
    )
    sp.add_argument(
        "--keep-temp", action="store_true", default=False,
        help=(
            "Keep temporary files after the run (keyed FASTA, windows FASTA). "
            "Useful for resuming a failed BLAST step or debugging."
        ),
    )


# ---------------------------------------------------------------------------
# Subcommand implementations
# ---------------------------------------------------------------------------

def _cmd_audit(args: argparse.Namespace) -> None:
    run_pipeline(args, filter_mode=False)


def _cmd_filter(args: argparse.Namespace) -> None:
    run_pipeline(args, filter_mode=True)


def _cmd_explain(args: argparse.Namespace) -> None:
    import csv
    matches = []
    with open(args.decisions, encoding="utf-8", errors="replace") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row.get("accession") == args.accession or row.get("record_key") == args.accession:
                matches.append(row)
    if not matches:
        raise SystemExit(f"Accession {args.accession!r} not found in {args.decisions}")
    if len(matches) > 1:
        print(
            f"Note: {len(matches)} rows found for accession {args.accession!r}. "
            f"Use --accession with the full record_key (e.g. {args.accession}__dup2) to select a specific one.\n"
        )
    for row in matches:
        reasons = [r for r in (row.get("reasons") or "").split(";") if r]
        print(f"Record key : {row.get('record_key')}")
        print(f"Accession  : {row.get('accession')}")
        print(f"Decision   : {row.get('decision')}")
        print(f"Score      : {row.get('decision_score')}")
        print(f"Class      : {row.get('marker_class')} / {row.get('region_id')}")
        print(f"Kingdom    : {row.get('kingdom')}")
        print(f"Species    : {row.get('species_guess')}")
        print(f"Reasons    :")
        for r in reasons:
            print(f"  - {r}")
        print(f"Header     : {row.get('header')}")
        print()


def _cmd_report(args: argparse.Namespace) -> None:
    report_from_decisions(args.decisions, args.outprefix)
    print(f"[Spinner] Report written: {args.outprefix}.summary.tsv")
    print(f"[Spinner] Report written: {args.outprefix}.summary.html")


def _cmd_init_config(args: argparse.Namespace) -> None:
    dest = Path(args.outdir)
    dest.mkdir(parents=True, exist_ok=True)
    src = Path(__file__).resolve().parent.parent / "configs"
    if not src.is_dir():
        raise SystemExit(f"Configs directory not found: {src}")
    for p in src.iterdir():
        if p.is_file():
            shutil.copy2(p, dest / p.name)
            print(f"[Spinner] Copied {p.name} -> {dest}")


# ---------------------------------------------------------------------------
# Argument parser
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="Spinner",
        description=(
            "TackleBox Spinner — reference sequence curation and QC.\n"
            "Curates NCBI / custom reference FASTAs for bait design, read mapping,\n"
            "BLAST databases, or metagenomic classifiers (MEGAN, Kraken, Bracken).\n\n"
            "Every input record is annotated and assigned KEEP / REVIEW / REJECT with\n"
            "a full explanation — nothing is silently discarded."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Common usage:\n"
            "\n"
            "  # Reference database (mapping / metabarcoding / MEGAN):\n"
            "  Spinner filter --fasta refs.fasta --species-kingdom sk.tsv \\\n"
            "    --mode reference_db --outprefix out/run_name\n"
            "\n"
            "  # Bait panel (input for FlyForge capture enrichment design):\n"
            "  Spinner filter --fasta refs.fasta --species-kingdom sk.tsv \\\n"
            "    --mode bait_panel --outprefix out/run_name\n"
            "\n"
            "  # Quick audit without writing FASTAs (preview decisions first):\n"
            "  Spinner audit --fasta refs.fasta --outprefix out/preview\n"
            "\n"
            "  # Explain a single record:\n"
            "  Spinner explain --decisions out/run_name.decisions.tsv --accession MK123456.1\n"
            "\n"
            "  # Copy bundled configs for project-specific customisation:\n"
            "  Spinner init-config --outdir my_configs\n"
            "\n"
            "BLAST performance tips:\n"
            "  - Set open file limit before running: ulimit -n 65536\n"
            "  - Set num_threads in your config (taxonomy_blast.num_threads, windowed_blast.num_threads)\n"
            "  - Download taxdb for taxonomy name resolution:\n"
            "      cd /path/to/blast_db && wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz\n"
            "      tar -xzf taxdb.tar.gz\n"
        ),
    )
    p.add_argument("--version", action="version", version=f"Spinner {VERSION}")
    sub = p.add_subparsers(dest="cmd", required=True, title="subcommands")

    # audit
    a = sub.add_parser(
        "audit",
        help="Annotate and write decisions TSV — no FASTA output. "
             "Use to preview decisions and tune settings before a full filter run. "
             "Accepts --mode and all filter flags.",
    )
    _add_common(a)
    a.set_defaults(func=_cmd_audit)

    # filter
    f = sub.add_parser(
        "filter",
        help="Full run: annotate, screen, and write keep / review / reject FASTAs. "
             "Use --mode reference_db or --mode bait_panel for full-stack runs with BLAST.",
    )
    _add_common(f)
    f.set_defaults(func=_cmd_filter)

    # screen-taxonomy
    st = sub.add_parser(
        "screen-taxonomy",
        help="Run audit with taxonomy BLAST sanity check enabled.",
    )
    _add_common(st)
    st.add_argument(
        "--blast-db", dest="taxonomy_blast_db", required=True, metavar="DB",
        help="BLAST database for taxonomy checks (e.g. RefSeq, nt, or custom curated DB).",
    )
    st.set_defaults(func=_cmd_audit)

    # screen-vector
    sv = sub.add_parser(
        "screen-vector",
        help="Run audit with vector/UniVec BLAST screening enabled.",
    )
    _add_common(sv)
    sv.add_argument(
        "--blast-db", dest="vector_blast_db", required=True, metavar="DB",
        help="Vector/UniVec BLAST database path.",
    )
    sv.set_defaults(func=_cmd_audit)

    # screen-windowed
    sw = sub.add_parser(
        "screen-windowed",
        help="Run audit with windowed BLAST chimerism screen enabled.",
    )
    _add_common(sw)
    sw.add_argument(
        "--blast-db", dest="windowed_blast_db", required=True, metavar="DB",
        help="BLAST database for windowed chimerism checks.",
    )
    sw.set_defaults(func=_cmd_audit)

    # cluster
    cl = sub.add_parser(
        "cluster",
        help="Run audit with vsearch clustering enabled.",
    )
    _add_common(cl)
    cl.set_defaults(func=_cmd_audit, enable_cluster=True)

    # explain
    e = sub.add_parser(
        "explain",
        help="Explain why a specific accession was kept / reviewed / rejected.",
    )
    e.add_argument("--decisions", required=True, metavar="TSV",
                   help="Path to a decisions TSV file.")
    e.add_argument("--accession", required=True, metavar="ACC",
                   help="Accession or record_key to look up.")
    e.set_defaults(func=_cmd_explain)

    # report
    r = sub.add_parser(
        "report",
        help="Re-generate summary TSV and HTML from an existing decisions TSV.",
    )
    r.add_argument("--decisions", required=True, metavar="TSV",
                   help="Path to an existing decisions TSV file.")
    r.add_argument("--outprefix", required=True, metavar="PREFIX",
                   help="Output file prefix for the regenerated report.")
    r.set_defaults(func=_cmd_report)

    # init-config
    ic = sub.add_parser(
        "init-config",
        help="Copy starter configs to a new directory for project-specific customisation.",
    )
    ic.add_argument("--outdir", required=True, metavar="DIR",
                    help="Destination directory for config copies.")
    ic.set_defaults(func=_cmd_init_config)

    return p


def main(argv=None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    # Dispatch to subcommand function stored in args.func.
    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()
        sys.exit(1)
