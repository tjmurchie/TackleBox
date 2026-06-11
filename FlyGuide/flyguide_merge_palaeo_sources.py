#!/usr/bin/env python3
"""
Merge multiple FlyGuide-compatible palaeo-source CSVs.

Inputs can be outputs from neotoma_extinct_to_gbif.py, pbdb_to_gbif.py,
now_to_gbif.py, GBIF fossil exports that already have species/genus/kingdom,
or any similar CSV/TSV. The merged output remains GBIF-like and can be fed into
FlyGuide/gbif_prep_from_csv.py or written directly as FlyGuide species files.
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from collections import defaultdict
from typing import Dict, List, Optional, Sequence, Tuple

try:
    from _palaeo_tui import PalaeoProgress, print_header, print_done
except ImportError:
    PalaeoProgress = None  # type: ignore[assignment,misc]
    print_header = print_done = None  # type: ignore[assignment]

VERSION = "1.3.0"


def sniff(path: str) -> str:
    with open(path, "r", encoding="utf-8-sig", errors="replace") as fh:
        first = fh.read(4096)
    return "\t" if first.count("\t") > first.count(",") else ","


def lower_map(fields: Sequence[str]) -> Dict[str, str]:
    return {f.lower(): f for f in fields if f}


def read_rows(path: str) -> List[Dict[str, str]]:
    delim = sniff(path)
    with open(path, "r", encoding="utf-8-sig", errors="replace", newline="") as fh:
        reader = csv.DictReader(fh, delimiter=delim)
        rows = []
        for row in reader:
            clean = {str(k): ("" if v is None else str(v).strip()) for k, v in row.items() if k is not None}
            clean["__input_file"] = path
            rows.append(clean)
    return rows


def get(row: Dict[str, str], *names: str) -> str:
    lmap = {k.lower(): k for k in row}
    for name in names:
        key = lmap.get(name.lower())
        if key and row.get(key):
            return row[key].strip()
    return ""


def canonical_key(name: str, mode: str) -> str:
    tokens = re.sub(r"[^A-Za-z -]+", " ", name or "").split()
    if not tokens:
        return ""
    if mode == "genus":
        return tokens[0]
    if mode == "trinomial" and len(tokens) >= 3:
        return " ".join(tokens[:3])
    if len(tokens) >= 2:
        return " ".join(tokens[:2])
    return tokens[0]


def merge(inputs: Sequence[str], collapse: str, source_priority: Sequence[str],
          verbose: bool = True) -> List[Dict[str, str]]:
    grouped: Dict[str, Dict[str, str]] = {}
    source_rank = {s.lower(): i for i, s in enumerate(source_priority)}
    all_rows = [(path, row) for path in inputs for row in read_rows(path)]
    dash = PalaeoProgress(len(all_rows), label="input rows", verbose=verbose) if PalaeoProgress else None
    _last_key = ""
    for _i, (path, row) in enumerate(all_rows):
            name = get(row, "ncbi_search_name", "species", "scientificName")
            key = canonical_key(name, collapse)
            if not key:
                continue
            source = get(row, "source", "source_database") or os.path.basename(path)
            kingdom = get(row, "kingdom")
            phylum = get(row, "phylum")
            genus = get(row, "genus") or key.split()[0]
            existing = grouped.get(key)
            current_rank = source_rank.get(source.lower(), 999)
            existing_rank = source_rank.get((existing or {}).get("preferred_source", "").lower(), 999)
            if existing is None:
                grouped[key] = {
                    "species": key,
                    "genus": genus,
                    "kingdom": kingdom,
                    "phylum": phylum,
                    "scientificName": key,
                    "taxonRank": "species" if len(key.split()) >= 2 else "genus_or_higher",
                    "taxonKey": get(row, "taxonKey"),
                    "source": "merged_palaeo_sources",
                    "source_database": "merged_palaeo_sources",
                    "ncbi_search_name": key,
                    "canonical_binomial": canonical_key(key, "binomial"),
                    "name_cleaning_status": f"merged_from_{source}",
                    "preferred_source": source,
                    "source_files": path,
                    "source_databases": source,
                    "source_taxon_keys": get(row, "taxonKey"),
                    "regions": get(row, "region"),
                    "periods": get(row, "period"),
                }
            else:
                e = existing
                if current_rank < existing_rank:
                    e["preferred_source"] = source
                    e["taxonKey"] = get(row, "taxonKey") or e.get("taxonKey", "")
                if not e.get("kingdom") and kingdom:
                    e["kingdom"] = kingdom
                if not e.get("phylum") and phylum:
                    e["phylum"] = phylum
                def append_unique(field: str, val: str) -> None:
                    if not val:
                        return
                    items = [x for x in e.get(field, "").split(";") if x]
                    if val not in items:
                        items.append(val)
                    e[field] = ";".join(items)
                append_unique("source_files", path)
                append_unique("source_databases", source)
                append_unique("source_taxon_keys", get(row, "taxonKey"))
                append_unique("regions", get(row, "region"))
                append_unique("periods", get(row, "period"))
                e["name_cleaning_status"] = "merged_duplicate_sources"
            if dash and key and key != _last_key:
                _last_key = key
                dash.update(key, f"[{len(grouped):>6,} unique]")
    if dash:
        dash.finish()
    return [grouped[k] for k in sorted(grouped)]


def write_csv(path: str, rows: List[Dict[str, str]]) -> None:
    fields = [
        "species", "genus", "kingdom", "phylum", "scientificName", "taxonRank", "taxonKey",
        "source", "source_database", "ncbi_search_name", "canonical_binomial", "name_cleaning_status",
        "preferred_source", "source_databases", "source_taxon_keys", "regions", "periods", "source_files",
    ]
    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_flyguide_files(out_csv: str, rows: List[Dict[str, str]], out_prefix: Optional[str]) -> Tuple[str, str]:
    prefix = out_prefix or os.path.splitext(os.path.basename(out_csv))[0]
    species_path = f"{prefix}_species_search.txt"
    kingdom_path = f"{prefix}_species_kingdom.tsv"
    seen = {}
    for r in rows:
        name = r.get("ncbi_search_name") or r.get("species")
        if not name:
            continue
        seen[name] = (r.get("kingdom", ""), r.get("phylum", ""))
    with open(species_path, "w", encoding="utf-8") as fh:
        for name in sorted(seen):
            fh.write(name + "\n")
    with open(kingdom_path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        for name in sorted(seen):
            writer.writerow([name, seen[name][0], seen[name][1]])
    return species_path, kingdom_path


def main(argv: Optional[Sequence[str]] = None) -> int:
    p = argparse.ArgumentParser(
        description=(
            "Merge Neotoma/PBDB/NOW/GBIF-like palaeo taxon CSVs into a single "
            "deduplicated FlyGuide-compatible CSV. Each input must have at least "
            "a 'species' or 'ncbi_search_name' column. Source provenance is preserved "
            "in metadata columns; duplicates across sources are collapsed at the "
            "requested name level (binomial by default)."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples\n"
            "--------\n"
            "  # Merge Neotoma + PBDB + NOW outputs:\n"
            "  python3 flyguide_merge_palaeo_sources.py \\\n"
            "      --inputs neotoma_gbif.csv pbdb_gbif.csv now_gbif.csv \\\n"
            "      --out merged_palaeo_gbif.csv \\\n"
            "      --write-flyguide-files\n\n"
            "  # Trinomial collapse (keep subspecies where present):\n"
            "  python3 flyguide_merge_palaeo_sources.py \\\n"
            "      --inputs pbdb_gbif.csv now_gbif.csv \\\n"
            "      --collapse trinomial \\\n"
            "      --out merged_trinomial.csv\n\n"
            "  # Custom source priority (affects which source's metadata is preferred\n"
            "  # when the same name appears in multiple inputs):\n"
            "  python3 flyguide_merge_palaeo_sources.py \\\n"
            "      --inputs pbdb.csv neotoma.csv \\\n"
            "      --source-priority PBDB,Neotoma,NOW,GBIF \\\n"
            "      --out merged.csv --write-flyguide-files\n"
        ),
    )
    p.add_argument(
        "--inputs", nargs="+", required=True, metavar="FILE",
        help="Input CSV/TSV files to merge. Accepts any FlyGuide-compatible GBIF-like "
             "CSV: outputs of neotoma_extinct_to_gbif.py, pbdb_to_gbif.py, now_to_gbif.py, "
             "or standard GBIF downloads. Auto-detects comma/tab delimiter.",
    )
    p.add_argument(
        "--out", required=True, metavar="FILE",
        help="Output merged GBIF-like CSV. Contains at least: species, genus, kingdom, "
             "phylum, ncbi_search_name, source_databases, and name_cleaning_status.",
    )
    p.add_argument(
        "--collapse", choices=["binomial", "trinomial", "genus"], default="binomial",
        help="Name level for deduplication. 'binomial' (default) collapses to Genus species "
             "— best for NCBI retrievability since many extinct subspecies are not in NCBI. "
             "'trinomial' keeps Genus species subspecies when three tokens are present. "
             "'genus' collapses everything to the genus — useful for very coarse buckets.",
    )
    p.add_argument(
        "--source-priority", default="Neotoma,PBDB,NOW,GBIF", metavar="LIST",
        help="Comma-separated source names in preferred order. When the same taxon name "
             "appears in multiple inputs, the highest-priority source's metadata is used "
             "for kingdom/phylum. All sources are still listed in source_databases. "
             "Default: Neotoma,PBDB,NOW,GBIF",
    )
    p.add_argument(
        "--out-prefix", metavar="PREFIX",
        help="Output prefix for FlyGuide native files written by --write-flyguide-files. "
             "Defaults to the base name of --out (without extension).",
    )
    p.add_argument(
        "--write-flyguide-files", action="store_true",
        help="Also write PREFIX_species_search.txt (one name per line, sorted) and "
             "PREFIX_species_kingdom.tsv (name TAB kingdom TAB phylum). These are the "
             "files consumed directly by NCBI-NT_Downloader.pl.",
    )
    p.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    p.add_argument("--quiet", action="store_true", help="Suppress progress output")
    args = p.parse_args(argv)
    verbose = not args.quiet
    priority = [x.strip() for x in args.source_priority.split(",") if x.strip()]
    if print_header and verbose:
        print_header(
            f"FlyGuide palaeo merge  v{VERSION}",
            [
                ("Inputs", ", ".join(args.inputs)),
                ("Collapse", args.collapse),
                ("Priority", args.source_priority),
                ("Output", args.out),
            ],
        )
    rows = merge(args.inputs, args.collapse, priority, verbose=verbose)
    write_csv(args.out, rows)
    sp = kd = None
    if args.write_flyguide_files:
        sp, kd = write_flyguide_files(args.out, rows, args.out_prefix)
    if print_done and verbose:
        counts = [("Unique names written:", f"{len(rows):,}  →  {args.out}")]
        if sp:
            counts += [("FlyGuide species search:", sp), ("FlyGuide kingdom map:", kd)]
        print_done(f"Merged {len(args.inputs)} file(s)", counts)
    elif verbose:
        print(f"Merged {len(args.inputs)} files into {len(rows)} unique names -> {args.out}", file=sys.stderr)
        if sp:
            print(f"FlyGuide species search: {sp}", file=sys.stderr)
            print(f"FlyGuide kingdom map: {kd}", file=sys.stderr)
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
