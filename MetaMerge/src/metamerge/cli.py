"""Command-line interface for MetaMerge.

Commands
--------
metamerge check
    Validate inputs and show a pre-run summary.  Useful before committing to a
    full run on a large metaDMG CSV.

metamerge run
    Run the full merge/classification workflow and write output files.

metamerge report
    Render heatmaps from plot-input TSVs produced by a previous ``run``.
    Requires Rscript on PATH.

Example usage::

    # Validate inputs
    metamerge check \\
        --megan-counts megan_counts.tsv \\
        --holi metadmg.csv \\
        --linker library_linker.csv

    # Run full workflow
    metamerge run \\
        --megan-counts megan_counts.tsv \\
        --holi metadmg.csv \\
        --linker library_linker.csv \\
        --config my_project.yaml \\
        --outdir results/ \\
        --yes

    # Render graph figures from a previous run
    metamerge report \\
        --input-dir results/report_inputs \\
        --outdir results/reports
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import textwrap
import time
from datetime import datetime
from pathlib import Path

import pandas as pd

from . import __version__
from .classify import build_merge
from .common_names import load_common_name_overrides
from .config import load_config
from .defaults import DEFAULT_CONFIG
from .io import load_holi, load_megan_counts
from .metadata import load_metadata
from .report import write_plot_inputs
from .workbook import write_workbook
from . import linker as linker_mod


# ─── Path helpers ────────────────────────────────────────────────────────────

def _resolve_counts_paths(inputs: list[str]) -> list[str]:
    """Resolve --counts inputs to a list of file paths.

    A single directory is expanded to all .tsv and .csv files it contains.
    Individual file paths are returned as-is after existence checks.
    """
    paths = []
    for inp in inputs:
        p = Path(inp)
        if not p.exists():
            print(f"Error: path does not exist: {p}", file=sys.stderr)
            raise SystemExit(1)
        if p.is_dir():
            found = sorted(p.glob("*.tsv")) + sorted(p.glob("*.csv"))
            if not found:
                print(f"Warning: no .tsv/.csv files found in directory: {p}", file=sys.stderr)
            paths.extend(str(f) for f in found)
        else:
            paths.append(str(p))
    return paths


def _clone_actions(source_parser: argparse.ArgumentParser, target_parser: argparse.ArgumentParser) -> None:
    """Clone user-facing arguments from one parser into another."""
    for action in source_parser._actions:
        if isinstance(action, (argparse._HelpAction, argparse._VersionAction)):
            continue
        opts = action.option_strings or []
        kwargs = {
            "dest": action.dest,
            "default": action.default,
            "required": action.required,
            "help": action.help,
            "nargs": action.nargs,
            "choices": action.choices,
            "metavar": action.metavar,
            "type": action.type,
        }
        if isinstance(action, argparse._StoreTrueAction):
            target_parser.add_argument(*opts, dest=action.dest, default=action.default, help=action.help, action="store_true")
        elif isinstance(action, argparse._StoreFalseAction):
            target_parser.add_argument(*opts, dest=action.dest, default=action.default, help=action.help, action="store_false")
        else:
            # remove None values that argparse dislikes in some contexts
            kwargs = {k:v for k,v in kwargs.items() if v is not None}
            target_parser.add_argument(*opts, **kwargs)


def _get_linker_path(args) -> str:
    """Return the linker path from CLI args, accepting a deprecated alias."""
    linker = getattr(args, "linker", None)
    metadata_alias = getattr(args, "metadata", None)
    if linker:
        return linker
    if metadata_alias:
        print("Warning: --metadata is deprecated for 'check' and 'run'; use --linker instead.", file=sys.stderr)
        return metadata_alias
    raise SystemExit("Error: a linker file is required. Use --linker <library_linker.csv>.")


# ─── Banner helpers ───────────────────────────────────────────────────────────

_SEP  = "=" * 62
_SEP2 = "-" * 62


def _banner() -> None:
    """Print the MetaMerge startup banner."""
    print(_SEP)
    print(f"  MetaMerge v{__version__}  |  Ensemble aDNA classification")
    print(f"  {datetime.now().strftime('%Y-%m-%d  %H:%M:%S')}")
    print(_SEP)


def _section(title: str) -> None:
    print(f"\n{title}")
    print(_SEP2)


def _fmt_seconds(seconds: float) -> str:
    if seconds < 60:
        return f"{seconds:.1f} s"
    m, s = divmod(int(seconds), 60)
    return f"{m} min {s} s"


# ─── Input summary ────────────────────────────────────────────────────────────

def summarize_inputs(
    megan_df: pd.DataFrame,
    holi_df: pd.DataFrame,
    meta_df: pd.DataFrame,
) -> dict:
    """Build a compact pre-run input summary dict."""
    matched_megan = sum(col in megan_df.columns for col in meta_df["megan_library_name"])
    n_pos = int(meta_df["is_positive_control"].sum()) if "is_positive_control" in meta_df.columns else 0
    return {
        "megan_taxa_rows":               int(len(megan_df)),
        "holi_rows":                     int(len(holi_df)),
        "metadata_rows":                 int(len(meta_df)),
        "megan_libraries":               int(meta_df["megan_library_name"].nunique()),
        "holi_libraries":                int(holi_df["sample"].nunique()),
        "negative_controls":             int(meta_df["is_negative_control"].sum()),
        "positive_controls":             n_pos,
        "matched_megan_library_columns": int(matched_megan),
    }


def print_input_summary(summary: dict) -> None:
    """Print the pre-run input summary to stdout."""
    _section("Input summary")
    labels = {
        "megan_taxa_rows":               "MEGAN taxon rows",
        "holi_rows":                     "Holi/metaDMG rows",
        "metadata_rows":                 "Linker rows",
        "megan_libraries":               "MEGAN libraries (linker)",
        "holi_libraries":                "Holi samples (in file)",
        "negative_controls":             "Negative controls (blanks)",
        "positive_controls":             "Positive controls",
        "matched_megan_library_columns": "Linker libs matched in MEGAN",
    }
    for key, label in labels.items():
        val = summary.get(key, "—")
        print(f"  {label:<38} {val}")


def prompt_continue(summary: dict, assume_yes: bool) -> None:
    """Print the pre-run summary and optionally prompt for confirmation."""
    print_input_summary(summary)
    print()
    if assume_yes:
        return
    answer = input("Continue with merge/classification? [y/N] ").strip().lower()
    if answer not in {"y", "yes"}:
        print("Aborted.")
        raise SystemExit(1)


# ─── Post-run summary ─────────────────────────────────────────────────────────

def _build_run_report(
    input_summary: dict,
    run_summary: dict,
    merged_df: pd.DataFrame,
    outdir: Path,
    elapsed: float,
    args,
    counts_display: str = "",
) -> str:
    """Build the full text of the run-summary report."""
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    status_order = [
        "Very high confidence", "High confidence", "Supported",
        "Tentative", "Weak support", "Blank-associated",
    ]
    status_counts = run_summary.get("status_counts", {})

    lines = [
        _SEP,
        f"  MetaMerge v{__version__}  —  Run Summary Report",
        f"  Generated: {now}",
        _SEP,
        "",
        "INPUTS",
        _SEP2,
        f"  MEGAN count matrix : {counts_display}",
        f"  Holi/metaDMG CSV   : {args.holi}",
        f"  Library linker     : {_get_linker_path(args)}",
        f"  Config             : {getattr(args, 'config', None) or 'defaults'}",
        f"  Common names       : {'online (NCBI → GBIF → iNat)' if getattr(args, 'online_common_names', False) else 'offline only'}",
        f"  Output directory   : {outdir}",
        "",
        "INPUT SUMMARY",
        _SEP2,
    ]
    labels = {
        "megan_taxa_rows":               "MEGAN taxon rows",
        "holi_rows":                     "Holi/metaDMG rows",
        "metadata_rows":                 "Linker rows",
        "megan_libraries":               "MEGAN libraries",
        "holi_libraries":                "Holi samples",
        "negative_controls":             "Negative controls (blanks)",
        "positive_controls":             "Positive controls",
        "matched_megan_library_columns": "Matched MEGAN columns",
    }
    for key, label in labels.items():
        lines.append(f"  {label:<38} {input_summary.get(key, '—')}")

    lines += [
        "",
        "CLASSIFICATION RESULTS",
        _SEP2,
        f"  Total taxa classified: {run_summary.get('n_taxa', '—')}",
        "",
    ]
    for s in status_order:
        n = status_counts.get(s, 0)
        bar = "#" * min(n, 40)
        lines.append(f"  {s:<28} {n:>5}  {bar}")

    # Unmatched taxa
    unmatched = run_summary.get("unmatched_taxa_examples", [])
    if unmatched:
        lines += [
            "",
            "TAXA WITHOUT HOLI TAX_PATH (no lineage info — classification unaffected)",
            _SEP2,
        ]
        for t in unmatched[:20]:
            lines.append(f"  {t}")
        if len(unmatched) == 20:
            lines.append("  … (showing first 20)")

    # Top confident taxa
    for status in ["Very high confidence", "High confidence"]:
        top = merged_df[merged_df["aDNA_support_status"] == status].head(10)
        if not top.empty:
            lines += ["", f"TOP TAXA — {status.upper()}", _SEP2]
            for _, r in top.iterrows():
                cn = f"  [{r.get('common_name', '')}]" if r.get("common_name") else ""
                lines.append(
                    f"  {r.get('scientific_name', '?'):<35}{cn}  "
                    f"(rank: {r.get('tax_rank','?')}, "
                    f"max_count: {r.get('megan_max_count', 0):.0f})"
                )

    lines += [
        "",
        "OUTPUTS",
        _SEP2,
        f"  {outdir}/",
    ]
    lines += [
        f"    metamerge_merged_support.xlsx   — workbook with classification",
        f"    metamerge_merged_support.tsv    — tab-separated text",
        f"    run_summary.json                — machine-readable summary",
        f"    warnings.tsv                    — run warnings",
        f"    report_inputs/                  — long-format tables for heatmaps",
        f"    reports/run_report.txt          — this report",
        f"    common_name_cache.json          — cached GBIF lookups (if online)",
        "",
        "TIMING",
        _SEP2,
        f"  Total runtime: {_fmt_seconds(elapsed)}",
        "",
        _SEP,
    ]
    return "\n".join(lines)


def print_run_complete(
    input_summary: dict,
    run_summary: dict,
    merged_df: pd.DataFrame,
    outdir: Path,
    elapsed: float,
    args,
    counts_display: str = "",
) -> str:
    """Print the run-complete summary to stdout and return the report text."""
    report = _build_run_report(input_summary, run_summary, merged_df, outdir, elapsed, args,
                               counts_display=counts_display)
    print()
    print(report)
    return report


# ─── Parser ──────────────────────────────────────────────────────────────────


def _threshold_default(key: str):
    return DEFAULT_CONFIG["thresholds"][key]


def _threshold_help(key: str, default, meaning: str, when_to_change: str = "") -> str:
    msg = f"{meaning} Default: {default}."
    if when_to_change:
        msg += f" Adjust when: {when_to_change}"
    return msg


def _has_validation_issues(summary: dict, metadata: pd.DataFrame, holi: pd.DataFrame) -> tuple[bool, list[str]]:
    issues: list[str] = []
    matched = summary.get("matched_megan_library_columns", 0)
    total = summary.get("megan_libraries", 0)
    if matched != total:
        issues.append(
            f"Only {matched}/{total} linker libraries matched MEGAN columns. Check that linker megan_library_name values exactly match the MEGAN headers."
        )
    if metadata["holi_library_name"].astype(str).str.startswith("REVIEW_NEEDED", na=False).any():
        bad = metadata.loc[metadata["holi_library_name"].astype(str).str.startswith("REVIEW_NEEDED", na=False), "merged_library_name"].astype(str).tolist()
        issues.append(
            "Some linker rows still have REVIEW_NEEDED holi_library_name values: " + ", ".join(bad[:8]) + (" …" if len(bad) > 8 else "")
        )
    if metadata["holi_library_name"].astype(str).eq("UNKNOWN").any() or ("sample_id" in metadata.columns and metadata["sample_id"].astype(str).eq("UNKNOWN").any()):
        issues.append("Some linker rows still contain UNKNOWN sample identifiers. Fill these in before running.")
    return (len(issues) > 0, issues)


def _handle_run_validation(summary: dict, metadata: pd.DataFrame, holi: pd.DataFrame, assume_yes: bool) -> None:
    print_input_summary(summary)
    has_issues, issues = _has_validation_issues(summary, metadata, holi)
    if not has_issues:
        print("\n  Validation passed with no blocking issues. Continuing automatically…")
        return

    _section("Blocking validation issues")
    for issue in issues:
        print(f"  - {issue}")
    print(textwrap.dedent("""

        Fix these issues in the linker first, then rerun `metamerge check` or `metamerge run`.
        Use --yes only when you intentionally want to continue despite these warnings
        (for example during debugging or when you know some libraries legitimately
        have NO_HOLI_DATA).
    """))
    if assume_yes:
        print("  --yes supplied: continuing despite validation issues.\n")
        return
    raise SystemExit(1)


def build_parser() -> argparse.ArgumentParser:
    """Construct the MetaMerge CLI argument parser."""
    parser = argparse.ArgumentParser(
        prog="metamerge",
        description=(
            f"MetaMerge v{__version__}: merge conservative MEGAN count matrices with Holi/metaDMG outputs, "
            "assign DNA-support categories, and generate report-ready tables/graphs."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent(f"""\
            Version
            -------
            MetaMerge v{__version__}

            Help tip
            --------
            Use `metamerge <subcommand> -h` for subcommand-specific help.
            Example: `metamerge run -h`

            Recommended workflow
            --------------------
            Step 0 — build the library linker from raw project metadata:
              metamerge linker --metadata metadata.xlsx --megan-counts megan_dir/ --holi metadmg.csv --out linker.csv

            Step 1 — validate the three core run inputs:
              metamerge check --megan-counts megan_dir/ --holi metadmg.csv --linker linker.csv

            Step 2 — run the merge/classifier and optionally render graphs:
              metamerge run --megan-counts megan_dir/ --holi metadmg.csv --linker linker.csv --outdir results/ --render-graphs

            Step 3 — rerender graphs from existing report tables if needed:
              metamerge report --input-dir results/report_inputs --outdir results/reports
        """),
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    linker_cmd = sub.add_parser(
        "linker",
        help="Build a MetaMerge linker from raw metadata + MEGAN columns.",
        description=f"MetaMerge v{__version__} — build a library linker from raw metadata, MEGAN counts, and optional Holi/metaDMG validation.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    _clone_actions(linker_mod.build_arg_parser(), linker_cmd)
    linker_cmd.add_argument("--version", action="version", version=f"metamerge linker {__version__}")

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument(
        "--megan-counts", required=True, nargs="+",
        help=(
            "One or more MEGAN count-matrix files (.xlsx, .csv, or .tsv), or a single directory containing them. "
            "When libraries are spread across multiple MEGAN runs, list all files here and MetaMerge will merge them on tax_id before classification."
        ),
    )
    common.add_argument("--holi", required=True, help="Path to the Holi/metaDMG CSV file used for per-library damage support.")
    common.add_argument(
        "--linker",
        help=(
            "Path to the library linker file (.csv, .tsv, or .xlsx). This is the output of `metamerge linker` and is the file consumed by `check` and `run`."
        ),
    )
    common.add_argument("--metadata", dest="metadata", help=argparse.SUPPRESS)
    common.add_argument(
        "--config",
        help=(
            "Optional YAML config file overriding default thresholds and settings. Use this when you want to change many parameters at once. "
            "Individual CLI threshold flags always take precedence over this file."
        ),
    )

    thr = common.add_argument_group(
        "threshold overrides",
        "Override individual classification thresholds without editing a config file. Defaults are tuned to conservative ancient-DNA screening and can be relaxed or tightened depending on project goals.",
    )
    thr.add_argument("--damage-min", type=float, metavar="FLOAT", help=_threshold_help("damage_min", _threshold_default("damage_min"), "Minimum estimated terminal damage fraction required for exact damage-supported calls.", "raise for stricter authentication; lower slightly for very low-coverage libraries."))
    thr.add_argument("--significance-min", type=float, metavar="FLOAT", help=_threshold_help("significance_min", _threshold_default("significance_min"), "Minimum damage significance required for exact damage-supported calls.", "raise if you want fewer false positives; lower if ancient signal is weak but expected."))
    thr.add_argument("--high-confidence-n-reads-min", type=int, metavar="INT", help=_threshold_help("high_confidence_n_reads_min", _threshold_default("high_confidence_n_reads_min"), "Minimum Holi/metaDMG reads for the top-tier damage-supported label.", "lower for sparse datasets; raise when only abundant hits should receive the top label."))
    thr.add_argument("--strong-count-min-reads", type=int, metavar="INT", help=_threshold_help("strong_count_min_reads", _threshold_default("strong_count_min_reads"), "Minimum MEGAN reads in a library to count as strong count support.", "lower for sparse targeted datasets; raise for broad shotgun datasets."))
    thr.add_argument("--strong-count-min-libraries", type=int, metavar="INT", help=_threshold_help("strong_count_min_libraries", _threshold_default("strong_count_min_libraries"), "Minimum number of positive non-control libraries for strong count support.", "set to 1 when only one library per sample exists; raise for stricter reproducibility."))
    thr.add_argument("--blank-absolute-min", type=int, metavar="INT", help=_threshold_help("blank_absolute_min", _threshold_default("blank_absolute_min"), "Absolute blank/control count at which contamination concern starts to matter.", "raise if blanks are noisy but expected; lower for ultra-clean workflows."))
    thr.add_argument("--blank-relative-min", type=float, metavar="FLOAT", help=_threshold_help("blank_relative_min", _threshold_default("blank_relative_min"), "Blank-to-sample count ratio threshold used in Blank-associated calls.", "lower for stricter contamination control; raise if true taxa also appear weakly in controls."))
    thr.add_argument("--tentative-min-reads", type=int, metavar="INT", help=_threshold_help("tentative_min_reads", _threshold_default("tentative_min_reads"), "Minimum reads for a taxon to move above Weak support when stronger evidence is missing."))
    thr.add_argument("--weak-support-min-reads", type=int, metavar="INT", help=_threshold_help("weak_support_min_reads", _threshold_default("weak_support_min_reads"), "Minimum reads required to retain a taxon as Weak support rather than dropping it entirely."))
    thr.add_argument("--qc-alignments-per-read-clean-max", type=float, metavar="FLOAT", help=_threshold_help("qc_alignments_per_read_clean_max", _threshold_default("qc_alignments_per_read_clean_max"), "Maximum alignments-per-read ratio still considered clean mapping behavior.", "lower if multimapping is a major concern."))
    thr.add_argument("--qc-alignments-per-read-caution-max", type=float, metavar="FLOAT", help=_threshold_help("qc_alignments_per_read_caution_max", _threshold_default("qc_alignments_per_read_caution_max"), "Above this alignments-per-read ratio, support is flagged as strong caution."))
    thr.add_argument("--qc-abs-rho-clean-max", type=float, metavar="FLOAT", help=_threshold_help("qc_abs_rho_clean_max", _threshold_default("qc_abs_rho_clean_max"), "Maximum absolute rho_Ac still treated as clean fit quality."))
    thr.add_argument("--qc-abs-rho-caution-max", type=float, metavar="FLOAT", help=_threshold_help("qc_abs_rho_caution_max", _threshold_default("qc_abs_rho_caution_max"), "Above this absolute rho_Ac value, fit quality is flagged as strong caution."))

    common.add_argument("--common-names-file", help="Optional CSV/TSV/XLSX with custom common-name overrides. Columns: scientific_name, tax_id, common_name (tax_id match takes priority).")
    common.add_argument(
        "--online-common-names", action="store_true",
        help=(
            "Query online taxonomy databases for missing common names.  "
            "Tries three sources in order: (1) NCBI Taxonomy eutils (using NCBI tax_id), "
            "(2) GBIF Species API (English vernacular names only), "
            "(3) iNaturalist Taxon API (English preferred_common_name).  "
            "Results are cached in common_name_cache.json so repeated runs are fast.  "
            "Adds a few seconds per uncached taxon."
        ),
    )
    common.add_argument("--outdir", default="metamerge_output", help="Output directory for workbook, TSV, JSON summaries, report inputs, and optional graphs. Default: metamerge_output/")
    common.add_argument("--yes", action="store_true", help="Non-interactive / force mode. Skip manual confirmation and continue even if validation reports linker issues. Useful for batch jobs and debugging.")

    check_cmd = sub.add_parser("check", parents=[common], help="Validate inputs and print a pre-run summary.", description=f"MetaMerge v{__version__} — validate MEGAN counts, Holi/metaDMG CSV, and linker alignment before running the full merge.")
    check_cmd.add_argument("--version", action="version", version=f"metamerge check {__version__}")

    run_cmd = sub.add_parser("run", parents=[common], help="Run the full merge/classification workflow.", description=f"MetaMerge v{__version__} — validate inputs, merge MEGAN + Holi/metaDMG, classify taxa, write outputs, and optionally render graphs.")
    run_cmd.add_argument("--version", action="version", version=f"metamerge run {__version__}")
    run_cmd.add_argument("--render-graphs", action="store_true", help="Render heatmaps and stacked-bar PDFs after writing report_inputs. Requires Rscript on PATH.")
    run_cmd.add_argument("--render-heatmaps", dest="render_graphs", action="store_true", help=argparse.SUPPRESS)
    run_cmd.add_argument(
        "--max-rank",
        dest="max_rank_for_reports",
        default=None,
        metavar="RANK",
        help=(
            "Most inclusive (highest/broadest) taxonomic rank to include in heatmap/report "
            "tables.  All ranks at or below this level are shown.  "
            "Accepted values: species, genus, family (default), order, class, phylum.  "
            "Use 'order' to include order-level taxa such as Proboscidea alongside "
            "family-level rows.  Use 'genus' to restrict figures to genus and below."
        ),
    )

    report_cmd = sub.add_parser("report", help="Render graphs from plot-input TSVs written by a previous run.", description=f"MetaMerge v{__version__} — render heatmaps and stacked-bar graphs from an existing report_inputs directory.")
    report_cmd.add_argument("--version", action="version", version=f"metamerge report {__version__}")
    report_cmd.add_argument("--input-dir", required=True, help="Directory containing report_inputs/*.tsv from a previous `metamerge run`.")
    report_cmd.add_argument("--outdir", default="metamerge_reports", help="Directory where graph PDFs and run_report.txt should be written. Default: metamerge_reports/")
    return parser


# ─── Sub-command implementations ──────────────────────────────────────────────

def run_heatmap_script(input_dir: Path, outdir: Path) -> None:
    """Run the bundled R graph script if Rscript is available on PATH."""
    rscript = shutil.which("Rscript")
    if not rscript:
        print(
            "\nNote: Rscript not found on PATH — skipping graph rendering.\n"
            "  Install R and ensure Rscript is on PATH, then re-run with:\n"
            "    metamerge report --input-dir ...",
            file=sys.stderr,
        )
        return
    script_path = Path(__file__).resolve().parents[2] / "r" / "metamerge_heatmap_report.R"
    if not script_path.exists():
        script_path = Path(__file__).resolve().parents[2] / "r" / "aedna_merge_heatmap_report.R"
    cmd = [rscript, str(script_path), "--input-dir", str(input_dir), "--outdir", str(outdir)]
    print("Running graph renderer:", " ".join(cmd))
    subprocess.run(cmd, check=True)


def build_warnings_df(
    metadata: pd.DataFrame,
    megan_df: pd.DataFrame,
    run_summary: dict,
) -> pd.DataFrame:
    """Collect run warnings and informational messages into a DataFrame."""
    warnings = []
    for lib in metadata["megan_library_name"]:
        if lib not in megan_df.columns:
            warnings.append({
                "level": "error",
                "message": f"Linker library missing from MEGAN counts: {lib}",
            })
    for taxon in run_summary.get("unmatched_taxa_examples", []):
        warnings.append({
            "level": "info",
            "message": f"No canonical Holi tax_path found for taxon: {taxon}",
        })
    if not warnings:
        warnings.append({"level": "info", "message": "No major warnings."})
    return pd.DataFrame(warnings)


def _apply_threshold_overrides(args, config: dict) -> None:
    """Apply any per-threshold CLI flag overrides to the config dict in-place.

    Per-threshold flags take precedence over both the built-in defaults and any
    --config YAML file.  Only flags that were explicitly set (non-None) are applied.
    """
    override_map = {
        "damage_min":                         getattr(args, "damage_min", None),
        "significance_min":                   getattr(args, "significance_min", None),
        "high_confidence_n_reads_min":        getattr(args, "high_confidence_n_reads_min", None),
        "strong_count_min_reads":             getattr(args, "strong_count_min_reads", None),
        "strong_count_min_libraries":         getattr(args, "strong_count_min_libraries", None),
        "blank_absolute_min":                 getattr(args, "blank_absolute_min", None),
        "blank_relative_min":                 getattr(args, "blank_relative_min", None),
        "tentative_min_reads":                getattr(args, "tentative_min_reads", None),
        "weak_support_min_reads":             getattr(args, "weak_support_min_reads", None),
        "qc_alignments_per_read_clean_max":   getattr(args, "qc_alignments_per_read_clean_max", None),
        "qc_alignments_per_read_caution_max": getattr(args, "qc_alignments_per_read_caution_max", None),
        "qc_abs_rho_clean_max":               getattr(args, "qc_abs_rho_clean_max", None),
        "qc_abs_rho_caution_max":             getattr(args, "qc_abs_rho_caution_max", None),
    }
    for key, val in override_map.items():
        if val is not None:
            config["thresholds"][key] = val


def command_check(args) -> None:
    """Validate inputs and print a pre-run summary (check sub-command)."""
    _banner()

    counts_paths = _resolve_counts_paths(args.megan_counts)
    counts_display = (
        counts_paths[0] if len(counts_paths) == 1
        else f"{len(counts_paths)} files ({Path(counts_paths[0]).parent})"
    )
    print(textwrap.dedent(f"""\

        Running:  metamerge check
        Purpose:  Validate that your three input files align correctly before
                  committing to a full merge run (which may take several minutes
                  on large metaDMG files).

        Inputs
        ------
          MEGAN counts : {counts_display}
          Holi CSV     : {args.holi}
          Linker       : {_get_linker_path(args)}
    """))

    print("Loading inputs — this may take a moment for large files…")
    t0 = time.time()
    config = load_config(args.config)
    config["taxonomy"]["online_common_names"] = bool(args.online_common_names)
    _apply_threshold_overrides(args, config)

    metadata = load_metadata(_get_linker_path(args), config)
    megan    = load_megan_counts(counts_paths, metadata, config)
    holi     = load_holi(args.holi, config)
    print(f"  Loaded in {_fmt_seconds(time.time() - t0)}.")

    summary = summarize_inputs(megan, holi, metadata)
    print_input_summary(summary)

    # Warn about linker issues.
    review_rows = metadata[
        metadata["holi_library_name"].str.startswith("REVIEW_NEEDED", na=False)
        | metadata["holi_library_name"].str.startswith("NO_HOLI_DATA", na=False)
    ]
    unknown_rows = metadata[metadata.get("sample_id", pd.Series(dtype=str)).eq("UNKNOWN")]

    if not review_rows.empty:
        print(
            f"\n  WARNING: {len(review_rows)} linker row(s) have holi_library_name "
            f"starting with REVIEW_NEEDED or NO_HOLI_DATA.\n"
            f"  These libraries will have no Holi/metaDMG support in the merge.\n"
            f"  Rows: {', '.join(review_rows['merged_library_name'].tolist()[:5])}"
        )
    if not unknown_rows.empty:
        print(
            f"\n  WARNING: {len(unknown_rows)} linker row(s) have sample_id=UNKNOWN.\n"
            f"  Review and correct the linker file before running."
        )

    _section("Validation result")
    matched = summary["matched_megan_library_columns"]
    total   = summary["megan_libraries"]
    if matched == total:
        print(f"  All {total} linker libraries found in the MEGAN count matrix.")
    else:
        print(
            f"  WARNING: only {matched}/{total} linker libraries matched MEGAN columns.\n"
            f"  Check that megan_library_name values match MEGAN column headers exactly."
        )

    print(textwrap.dedent(f"""
        Next steps
        ----------
        1. If the summary above looks correct, run the full merge:

             metamerge run \\
                 --megan-counts {" ".join(args.megan_counts)} \\
                 --holi         {args.holi} \\
                 --linker       {_get_linker_path(args)} \\
                 --outdir       results/

        2. To enable online common-name lookups via GBIF (optional):

             metamerge run ... --online-common-names

        3. To supply a custom config file overriding default thresholds:

             metamerge run ... --config my_project.yaml

           Key thresholds (edit in config YAML):
             damage_min              minimum C→T damage fraction (default 0.02)
             significance_min        minimum damage significance (default 2.0)
             high_confidence_n_reads_min  reads for Very high confidence (default 100)
             blank_max_fraction      blank/real ratio for Blank-associated (default 0.2)

        4. To render heatmaps after the run:

             metamerge run ... --render-graphs   (requires R + Rscript on PATH)

        See README.md for linker format details.
    """))


def command_run(args) -> None:
    """Run the full MetaMerge workflow (run sub-command)."""
    _banner()

    outdir       = Path(args.outdir)
    counts_paths = _resolve_counts_paths(args.megan_counts)
    counts_display = (
        counts_paths[0] if len(counts_paths) == 1
        else f"{len(counts_paths)} files ({Path(counts_paths[0]).parent})"
    )
    outdir.mkdir(parents=True, exist_ok=True)

    print(textwrap.dedent(f"""\

        Running:  metamerge run
        Output:   {outdir}/

        Inputs
        ------
          MEGAN counts : {counts_display}
          Holi CSV     : {args.holi}
          Linker       : {_get_linker_path(args)}
          Config       : {args.config or 'built-in defaults'}
          Common names : {'online (NCBI → GBIF → iNat)' if args.online_common_names else 'offline only'}

        Linker checklist (review before first use)
        ------------------------------------------
          - megan_library_name must exactly match a column header in the MEGAN
            count matrix (character-for-character, including the full file-path
            suffix that MEGAN sometimes exports).
          - holi_library_name must match a value in the 'sample' column of the
            Holi/metaDMG CSV.  Use NO_HOLI_DATA for libraries with no Holi data
            (MEGAN counts still included, damage support will be absent).
          - Negative controls (blank libraries) are detected if ANY of:
              is_negative_control = true / 1 / yes  (explicit boolean column)
              sample_type contains "blank"  (e.g. "extraction blank", "air blank")
              sample_type contains "negative"
              sample_type contains "control" but NOT "positive"
            Blank libraries drive Blank-associated classification.
          - Positive controls are detected if ANY of:
              is_positive_control = true / 1 / yes  (explicit boolean column)
              sample_type contains "positive"  (e.g. "positive control")
            Positive controls are plotted separately and never drive Blank-associated.
          - Real samples need none of the above — any row not flagged as a
            control is treated as a genuine biological sample.
          - If any row shows REVIEW_NEEDED in holi_library_name, fix it first.

        See README.md and docs/metadata_linker.md for full linker guidance.
    """))

    _section("Step 1 / 4  —  Loading inputs")
    t_total = time.time()
    t0 = time.time()
    config = load_config(args.config)
    config["taxonomy"]["online_common_names"] = bool(args.online_common_names)
    _apply_threshold_overrides(args, config)
    if getattr(args, "max_rank_for_reports", None):
        config["report"]["max_rank_for_reports"] = args.max_rank_for_reports

    metadata = load_metadata(_get_linker_path(args), config)
    print(f"  Linker      : {len(metadata)} libraries loaded.")

    print(f"  MEGAN counts: loading {len(counts_paths)} file(s)…", end=" ", flush=True)
    megan = load_megan_counts(counts_paths, metadata, config)
    print(f"{len(megan)} taxa.")

    print("  Holi CSV    : loading…", end=" ", flush=True)
    holi = load_holi(args.holi, config)
    print(f"{len(holi):,} rows, {holi['sample'].nunique()} samples.")
    print(f"  Done in {_fmt_seconds(time.time() - t0)}.")

    input_summary = summarize_inputs(megan, holi, metadata)
    _handle_run_validation(input_summary, metadata, holi, assume_yes=args.yes)

    _section("Step 2 / 4  —  Common names")
    common_overrides = load_common_name_overrides(args.common_names_file)
    if args.online_common_names:
        print(
            "  Online lookups enabled (NCBI Taxonomy → GBIF → iNaturalist).\n"
            "  Results will be cached in:\n"
            f"    {outdir}/common_name_cache.json\n"
            "  Network errors for any source are handled gracefully."
        )
    else:
        print(
            "  Using built-in common-name map only (offline mode).\n"
            "  Add --online-common-names to query NCBI, GBIF, and iNat."
        )
    if common_overrides:
        print(f"  User overrides loaded: {len(common_overrides)} entries.")

    _section("Step 3 / 4  —  Merge and classify")
    print(f"  Classifying {len(megan)} taxa across {len(metadata)} libraries…")
    t0 = time.time()
    merged, run_summary = build_merge(
        metadata=metadata,
        megan_df=megan,
        holi_df=holi,
        config=config,
        common_name_overrides=common_overrides,
        cache_path=outdir / "common_name_cache.json",
    )
    print(f"  Done in {_fmt_seconds(time.time() - t0)}.")

    _section("Step 4 / 4  —  Writing outputs")
    plot_input_paths = write_plot_inputs(merged, metadata, outdir=outdir, config=config)
    warnings_df      = build_warnings_df(metadata, megan, run_summary)

    prefix        = config["io"]["output_prefix"]
    workbook_path = outdir / f"{prefix}_merged_support.xlsx"
    tsv_path      = outdir / f"{prefix}_merged_support.tsv"

    write_workbook(
        merged_df=merged,
        metadata_df=metadata,
        run_summary=run_summary,
        warnings_df=warnings_df,
        outpath=workbook_path,
        config=config,
        plot_inputs=plot_input_paths,
    )
    merged.to_csv(tsv_path, sep="\t", index=False)
    warnings_df.to_csv(outdir / "warnings.tsv", sep="\t", index=False)

    with (outdir / "run_summary.json").open("w", encoding="utf-8") as fh:
        json.dump(run_summary, fh, indent=2, sort_keys=True)

    elapsed = time.time() - t_total
    report_text = print_run_complete(input_summary, run_summary, merged, outdir, elapsed, args,
                                     counts_display=counts_display)
    reports_dir = outdir / "reports"
    reports_dir.mkdir(parents=True, exist_ok=True)
    (reports_dir / "run_report.txt").write_text(report_text, encoding="utf-8")
    print(f"\n  Full report saved to: {reports_dir}/run_report.txt")

    print(textwrap.dedent(f"""
        What to do next
        ---------------
        1. Open {workbook_path.name} to review classifications.
           Start with the 'Key_results' sheet for the most immediately useful
           summary; 'Full_results' has all columns.
           Check the 'Warnings' sheet for any data-quality notes.

        2. Examine Blank-associated taxa carefully — high blank overlap may
           indicate modern contamination or ubiquitous environmental taxa.

        3. MetaMerge assigns categories based on DNA evidence only.  Additional
           lines of evidence (ecological plausibility, biogeographic range,
           macrofossil data) can be applied as a separate interpretive layer.

        4. If you want graph figures:
             metamerge report --input-dir {outdir}/report_inputs \\
                              --outdir {outdir}/reports
           (or re-run with --render-graphs if R is available)

        5. To adjust thresholds for a re-run, pass individual flags such as
           --damage-min, --significance-min, etc., or use --config for a full
           YAML override.  The run_summary.json records the exact config used.
    """))

    if getattr(args, "render_graphs", False):
        try:
            run_heatmap_script(outdir / "report_inputs", outdir / "reports")
        except subprocess.CalledProcessError:
            print(textwrap.dedent(f"""

                Graph rendering did not complete successfully.
                ---------------------------------------------
                The merge/classification stage already finished and the core outputs were written:
                  - {workbook_path}
                  - {tsv_path}
                  - {outdir / 'run_summary.json'}
                  - {outdir / 'report_inputs'}

                You do NOT need to rerun the full merge after fixing the R/graph issue.
                Once the graph problem is resolved, rerun only the reporting step:

                  metamerge report --input-dir {outdir / 'report_inputs'} --outdir {outdir / 'reports'}

                """))
            raise


def command_report(args) -> None:
    """Render heatmaps from existing plot-input TSVs (report sub-command)."""
    _banner()
    input_dir = Path(args.input_dir)
    outdir    = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    run_heatmap_script(input_dir=input_dir, outdir=outdir)


def main() -> None:
    """CLI entry point for the ``metamerge`` command."""
    # Convenience: allow `metamerge -h run` as a synonym for `metamerge run -h`.
    if len(sys.argv) > 2 and sys.argv[1] in {"-h", "--help"} and sys.argv[2] in {"linker", "check", "run", "report"}:
        sys.argv = [sys.argv[0], sys.argv[2], "--help"]
    # Let the dedicated linker parser handle its own help/usage cleanly.
    if len(sys.argv) > 1 and sys.argv[1] == "linker":
        linker_mod.main(sys.argv[2:])
        return

    parser = build_parser()
    args   = parser.parse_args()
    if args.command == "linker":
        linker_mod.main(sys.argv[2:])
    elif args.command == "check":
        command_check(args)
    elif args.command == "run":
        command_run(args)
    elif args.command == "report":
        command_report(args)
    else:
        parser.error(f"Unknown command: {args.command}")


if __name__ == "__main__":
    main()
