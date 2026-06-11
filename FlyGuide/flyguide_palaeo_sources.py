#!/usr/bin/env python3
"""
flyguide_palaeo_sources.py  —  one-command multi-source palaeo taxon builder

Pulls from Neotoma, PBDB, and optionally NOW, then merges into a single
deduplicated GBIF-like CSV ready for FlyGuide / NCBI-NT_Downloader.

Usage
-----
  python3 flyguide_palaeo_sources.py \\
      --region northern_hemisphere \\
      --period quaternary \\
      --organisms animals \\
      --out-prefix NH_quat_palaeo

  # With a pre-exported NOW table:
  python3 flyguide_palaeo_sources.py \\
      --region eurasia \\
      --period pleistocene \\
      --organisms mammals \\
      --now-input NOW_export.csv \\
      --out-prefix Eurasia_pleis_mammals

  # Smoke test (small PBDB limit, offline Neotoma fixture):
  python3 flyguide_palaeo_sources.py \\
      --region northern_hemisphere \\
      --period quaternary \\
      --organisms animals \\
      --limit 500 \\
      --out-prefix test_run

  # Chain directly into FlyGuide NCBI download:
  python3 flyguide_palaeo_sources.py \\
      --region north_america \\
      --period quaternary \\
      --organisms animals \\
      --out-prefix NA_quat_palaeo \\
      --run-flyguide your@email.org YOUR_API_KEY

Notes
-----
- PBDB is queried live (cached in .pbdb_cache/ by default).
- Neotoma is queried live (cached in .neotoma_cache/ by default).
- NOW requires a manually exported table (--now-input); it has no stable REST API.
- Organism vocabulary follows PBDB (richer); Neotoma receives a coarser equivalent.
  Broad categories 'animals', 'plants', 'all' work across all sources.
  Fine categories like 'mammals' are passed as-is to PBDB and broadened to
  'animals' for Neotoma (which does not sub-filter by clade).
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import textwrap
from typing import List, Optional, Sequence, Tuple

VERSION = "1.3.0"
HERE = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------------------
# Organism translation
# ---------------------------------------------------------------------------
# PBDB supports fine-grained clades; Neotoma only supports broad groups.
# Map PBDB vocabulary → Neotoma vocabulary.
_NEOTOMA_ORG_MAP = {
    "all": "all",
    "animals": "animals",
    "metazoa": "animals",
    "vertebrates": "animals",
    "chordates": "animals",
    "mammals": "animals",
    "birds": "animals",
    "fish": "animals",
    "reptiles": "animals",
    "reptiles_amphibians": "animals",
    "amphibians": "animals",
    "invertebrates": "animals",
    "arthropods": "animals",
    "molluscs": "animals",
    "plants": "plants",
    "fungi": "fungi",
    "diatoms": "diatoms",
    "protists": "protists",
}

# Organism group descriptions for --list-organisms
_ORGANISM_DESCRIPTIONS = {
    "all":                "All taxa (can be very large; use --limit for smoke tests)",
    "animals":            "All animals — Neotoma + PBDB Animalia",
    "vertebrates":        "All vertebrates [PBDB: Vertebrata; Neotoma: animals]",
    "chordates":          "All chordates [PBDB: Chordata; Neotoma: animals]",
    "mammals":            "Mammals [PBDB: Mammalia; Neotoma: animals; NOW: default]",
    "birds":              "Birds [PBDB: Aves; Neotoma: animals]",
    "fish":               "Fish [PBDB: Actinopterygii + more; Neotoma: animals]",
    "reptiles_amphibians":"Reptiles + amphibians [PBDB; Neotoma: animals]",
    "reptiles":           "Reptiles [PBDB: Reptilia; Neotoma: animals]",
    "amphibians":         "Amphibians [PBDB: Amphibia; Neotoma: animals]",
    "invertebrates":      "Invertebrates [PBDB: non-Chordata Animalia; Neotoma: animals]",
    "arthropods":         "Arthropods [PBDB: Arthropoda; Neotoma: animals]",
    "molluscs":           "Molluscs [PBDB: Mollusca; Neotoma: animals]",
    "plants":             "Plants [PBDB: Plantae; Neotoma: plants]",
    "fungi":              "Fungi [PBDB + Neotoma]",
    "diatoms":            "Diatoms [PBDB + Neotoma]",
    "protists":           "Protists [PBDB + Neotoma]",
}

# ---------------------------------------------------------------------------
# Period and region lists (delegated to child scripts for authoritative lists)
# ---------------------------------------------------------------------------

def _script(name: str) -> str:
    return os.path.join(HERE, name)


def _run(cmd: List[str], label: str, dry_run: bool = False, verbose: bool = True) -> int:
    if verbose:
        print(f"\n>>> {label}", flush=True)
        print("    " + " ".join(cmd), flush=True)
    if dry_run:
        print("    [dry-run: skipped]")
        return 0
    result = subprocess.run(cmd)
    return result.returncode


def _run_capture(cmd: List[str]) -> Tuple[int, str, str]:
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.returncode, result.stdout, result.stderr


def neotoma_organism(organism: str) -> str:
    """Translate PBDB-vocabulary organism to what Neotoma accepts."""
    return _NEOTOMA_ORG_MAP.get(organism.lower().replace("-", "_"), "animals")


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="flyguide_palaeo_sources.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent("""\
            One-command multi-source palaeo taxon builder.
            Queries Neotoma + PBDB (+ optional NOW), merges into a single
            deduplicated FlyGuide-compatible CSV.
        """),
        epilog=textwrap.dedent("""\
            Examples
            --------
              # Northern Hemisphere Quaternary animals (default):
              python3 flyguide_palaeo_sources.py --out-prefix NH_quat_palaeo

              # Eurasian Pleistocene mammals (with NOW):
              python3 flyguide_palaeo_sources.py \\
                  --region eurasia --period pleistocene --organisms mammals \\
                  --now-input NOW_export.csv --out-prefix Eurasia_pleis_mammals

              # Smoke test — small PBDB limit:
              python3 flyguide_palaeo_sources.py \\
                  --region tanzania --period quaternary --organisms mammals \\
                  --limit 200 --out-prefix test_tanzania

              # Chain into full NCBI download:
              python3 flyguide_palaeo_sources.py \\
                  --region north_america --period quaternary \\
                  --organisms animals --out-prefix NA_quat_palaeo \\
                  --run-flyguide your@email.org YOUR_NCBI_API_KEY
        """),
    )

    # Spatial / temporal filters
    filt = p.add_argument_group("Spatial / temporal filters")
    filt.add_argument("--region", default="northern_hemisphere",
                      help="Region preset (default: %(default)s)")
    filt.add_argument("--period", default="quaternary",
                      help="Period preset (default: %(default)s)")
    filt.add_argument("--organisms", default="animals",
                      help="Organism group (default: %(default)s). "
                           "Use --list-organisms to see options.")
    filt.add_argument("--status", choices=["all", "extinct", "extant"], default="all",
                      help="Neotoma taxon status filter (default: %(default)s)")

    # Source control
    src = p.add_argument_group("Source control")
    src.add_argument("--now-input", metavar="FILE",
                     help="Pre-exported NOW table (CSV/TSV/XLSX/ZIP). "
                          "NOW is skipped if not provided.")
    src.add_argument("--skip-neotoma", action="store_true",
                     help="Skip Neotoma query")
    src.add_argument("--skip-pbdb", action="store_true",
                     help="Skip PBDB query")
    src.add_argument("--limit", default=None, metavar="N",
                     help="Cap PBDB records per query (useful for smoke tests)")
    src.add_argument("--force-refresh", action="store_true",
                     help="Ignore cached API responses and re-query")

    # Output
    out = p.add_argument_group("Output")
    out.add_argument("--out-prefix", metavar="PREFIX", default=None,
                     help="Prefix for all output files. "
                          "Defaults to REGION_PERIOD_ORGANISMS_palaeo.")
    out.add_argument("--out", metavar="FILE", default=None,
                     help="Merged output CSV (default: PREFIX_merged_gbif.csv)")
    out.add_argument("--collapse", choices=["binomial", "trinomial", "genus"],
                     default="binomial",
                     help="Name deduplication level (default: %(default)s)")
    out.add_argument("--ncbi-name-mode",
                     choices=["binomial", "trinomial", "as-is"], default="binomial",
                     help="Name normalisation for NCBI queries (default: %(default)s)")

    # Downstream
    ds = p.add_argument_group("Downstream (optional)")
    ds.add_argument("--run-flyguide", nargs="+", metavar=("EMAIL", "API_KEY"),
                    help="Chain into flyguide.sh: provide EMAIL and optionally API_KEY")

    # Verbosity
    p.add_argument("--verbose", action="store_true", default=True)
    p.add_argument("--quiet", dest="verbose", action="store_false")
    p.add_argument("--dry-run", action="store_true",
                   help="Print commands without running them")

    # List modes
    p.add_argument("--list-regions", action="store_true",
                   help="List region presets (from PBDB module) and exit")
    p.add_argument("--list-periods", action="store_true",
                   help="List period presets (from PBDB module) and exit")
    p.add_argument("--list-organisms", action="store_true",
                   help="List organism presets and exit")
    p.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")

    return p


# ---------------------------------------------------------------------------
# Main logic
# ---------------------------------------------------------------------------

def main(argv: Optional[Sequence[str]] = None) -> int:
    p = build_parser()
    args = p.parse_args(argv)

    # ---- list modes --------------------------------------------------------
    if args.list_regions:
        rc, out, _ = _run_capture([sys.executable, _script("pbdb_to_gbif.py"), "--list-regions"])
        print(out)
        return rc

    if args.list_periods:
        rc, out, _ = _run_capture([sys.executable, _script("pbdb_to_gbif.py"), "--list-periods"])
        print(out)
        return rc

    if args.list_organisms:
        print(f"{'organism':<25} {'description'}")
        print("-" * 80)
        for k, desc in sorted(_ORGANISM_DESCRIPTIONS.items()):
            print(f"  {k:<23} {desc}")
        return 0

    # ---- derive output prefix / path ---------------------------------------
    prefix = args.out_prefix or f"{args.region}_{args.period}_{args.organisms}_palaeo"
    merged_csv = args.out or f"{prefix}_merged_gbif.csv"

    # ---- determine which sources to run -----------------------------------
    run_neotoma = not args.skip_neotoma
    run_pbdb = not args.skip_pbdb
    run_now = bool(args.now_input)

    if not run_neotoma and not run_pbdb and not run_now:
        print("ERROR: All sources skipped — nothing to do.", file=sys.stderr)
        return 1

    sources_running = (
        (["Neotoma"] if run_neotoma else []) +
        (["PBDB"] if run_pbdb else []) +
        (["NOW"] if run_now else [])
    )

    if args.verbose:
        print("=" * 70)
        print("  FlyGuide Palaeo Sources — multi-source taxon builder")
        print("=" * 70)
        print(f"  Region   : {args.region}")
        print(f"  Period   : {args.period}")
        print(f"  Organisms: {args.organisms}")
        print(f"  Status   : {args.status} (Neotoma only)")
        print(f"  Sources  : {', '.join(sources_running)}")
        print(f"  Output   : {merged_csv}")
        print("=" * 70)

    input_csvs: List[str] = []
    errors: List[str] = []

    # ---- Neotoma -----------------------------------------------------------
    if run_neotoma:
        neo_out = f"{prefix}_neotoma_gbif.csv"
        neo_org = neotoma_organism(args.organisms)
        if args.verbose and neo_org != args.organisms:
            print(f"\n[Neotoma] Note: '{args.organisms}' mapped to '{neo_org}' "
                  f"(Neotoma uses broad group categories)")
        cmd = [
            sys.executable, _script("neotoma_extinct_to_gbif.py"),
            "--region", args.region,
            "--period", args.period,
            "--organisms", neo_org,
            "--status", args.status,
            "--ncbi-name-mode", args.ncbi_name_mode,
            "--out", neo_out,
            "--write-flyguide-files",
            "--out-prefix", f"{prefix}_neotoma",
        ]
        if args.force_refresh:
            cmd.append("--force-refresh")
        if not args.verbose:
            cmd.append("--quiet")
        rc = _run(cmd, "Step 1/3 — Neotoma", dry_run=args.dry_run, verbose=args.verbose)
        if args.dry_run or rc == 0:
            input_csvs.append(neo_out)
        else:
            errors.append(f"Neotoma exited with code {rc}")
            if args.verbose:
                print(f"  [WARNING] Neotoma failed (exit {rc}) — will merge without it")

    # ---- PBDB --------------------------------------------------------------
    if run_pbdb:
        pbdb_out = f"{prefix}_pbdb_gbif.csv"
        cmd = [
            sys.executable, _script("pbdb_to_gbif.py"),
            "--region", args.region,
            "--period", args.period,
            "--organisms", args.organisms,
            "--ncbi-name-mode", args.ncbi_name_mode,
            "--out", pbdb_out,
            "--write-flyguide-files",
            "--out-prefix", f"{prefix}_pbdb",
        ]
        if args.limit:
            cmd += ["--limit", str(args.limit)]
        if args.force_refresh:
            cmd.append("--force-refresh")
        if not args.verbose:
            cmd += ["--quiet"]
        step = "Step 2/3" if run_neotoma else "Step 1"
        rc = _run(cmd, f"{step} — PBDB", dry_run=args.dry_run, verbose=args.verbose)
        if args.dry_run or rc == 0:
            input_csvs.append(pbdb_out)
        else:
            errors.append(f"PBDB exited with code {rc}")
            if args.verbose:
                print(f"  [WARNING] PBDB failed (exit {rc}) — will merge without it")

    # ---- NOW ---------------------------------------------------------------
    if run_now:
        now_out = f"{prefix}_now_gbif.csv"
        cmd = [
            sys.executable, _script("now_to_gbif.py"),
            "--input", args.now_input,
            "--region", args.region,
            "--period", args.period,
            "--ncbi-name-mode", args.ncbi_name_mode,
            "--out", now_out,
            "--write-flyguide-files",
            "--out-prefix", f"{prefix}_now",
        ]
        if not args.verbose:
            cmd += ["--quiet"]
        step = "Step 3/3" if (run_neotoma and run_pbdb) else "Step 2" if (run_neotoma or run_pbdb) else "Step 1"
        rc = _run(cmd, f"{step} — NOW", dry_run=args.dry_run, verbose=args.verbose)
        if args.dry_run or rc == 0:
            input_csvs.append(now_out)
        else:
            errors.append(f"NOW exited with code {rc}")
            if args.verbose:
                print(f"  [WARNING] NOW failed (exit {rc}) — will merge without it")

    # ---- Merge -------------------------------------------------------------
    if not args.dry_run and not input_csvs and not errors:
        print("\nERROR: No source produced output — cannot merge.", file=sys.stderr)
        return 1

    merge_step = len(sources_running) + 1
    cmd = [
        sys.executable, _script("flyguide_merge_palaeo_sources.py"),
        "--inputs", *input_csvs,
        "--out", merged_csv,
        "--collapse", args.collapse,
        "--out-prefix", prefix,
        "--write-flyguide-files",
    ]
    rc = _run(cmd, f"Step {merge_step}/{merge_step} — Merge",
              dry_run=args.dry_run, verbose=args.verbose)
    if rc != 0:
        errors.append(f"Merge exited with code {rc}")

    # ---- Optionally chain into flyguide.sh --------------------------------
    if args.run_flyguide and not args.dry_run and rc == 0:
        email = args.run_flyguide[0]
        api_key = args.run_flyguide[1] if len(args.run_flyguide) > 1 else None
        cmd = [os.path.join(HERE, "flyguide.sh"), merged_csv, prefix, email]
        if api_key:
            cmd.append(api_key)
        rc_fg = _run(cmd, "FlyGuide NCBI download", dry_run=False, verbose=args.verbose)
        if rc_fg != 0:
            errors.append(f"flyguide.sh exited with {rc_fg}")

    # ---- Summary ----------------------------------------------------------
    if args.verbose:
        print("\n" + "=" * 70)
        print("  Summary")
        print("=" * 70)
        if args.dry_run:
            print("  [dry-run mode: no files written]")
        else:
            for path in input_csvs:
                tag = (
                    "Neotoma" if "_neotoma_" in path else
                    "PBDB"    if "_pbdb_"    in path else
                    "NOW"     if "_now_"     in path else path
                )
                size = os.path.getsize(path) if os.path.exists(path) else 0
                print(f"  {tag:<10} {path}  ({size:,} bytes)")
            if os.path.exists(merged_csv):
                _count_rows(merged_csv)
            sp = f"{prefix}_species_search.txt"
            kd = f"{prefix}_species_kingdom.tsv"
            if os.path.exists(sp):
                with open(sp) as fh:
                    n = sum(1 for _ in fh)
                print(f"\n  Merged search names  : {n}  → {sp}")
            if os.path.exists(kd):
                print(f"  Kingdom map          : {kd}")
        if errors:
            print("\n  Warnings / errors:")
            for e in errors:
                print(f"    - {e}")
        print()

    return 0 if not errors else 1


def _count_rows(path: str) -> None:
    try:
        with open(path, "r", encoding="utf-8") as fh:
            n = sum(1 for _ in fh) - 1  # subtract header
        print(f"\n  Merged CSV           : {n} unique taxa → {path}")
    except OSError:
        pass


if __name__ == "__main__":
    raise SystemExit(main())
