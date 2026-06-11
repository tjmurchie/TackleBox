#!/usr/bin/env python3
"""
TackleBox: FlyGuide - GBIF prep helper

Given a GBIF CSV/TSV with species, genus, kingdom (and optionally phylum)
columns, produce:
  - OUTPREFIX_species_search.txt        one search name per line
  - OUTPREFIX_species_kingdom.tsv       species <TAB> kingdom [<TAB> phylum]
"""

import csv
import sys
import os

def main():
    if len(sys.argv) != 3:
        sys.stderr.write(
            "Usage: gbif_prep_from_csv.py GBIF_download.csv OUTPREFIX\n"
            "  GBIF_download.csv : GBIF occurrence or checklist file (CSV or TSV;\n"
            "                      must contain columns 'species', 'genus', 'kingdom';\n"
            "                      'phylum' is optional but used when present)\n"
            "  OUTPREFIX         : prefix for generated files (written in current directory)\n"
        )
        sys.exit(1)

    in_path    = sys.argv[1]
    out_prefix = sys.argv[2]

    if not os.path.exists(in_path):
        sys.stderr.write(f"ERROR: Input file not found: {in_path}\n")
        sys.exit(1)

    species_search_path  = f"{out_prefix}_species_search.txt"
    species_kingdom_path = f"{out_prefix}_species_kingdom.tsv"

    species_for_search    = set()
    species_kingdom_pairs = {}   # sp -> (kingdom, phylum)

    with open(in_path, "r", encoding="utf-8-sig", errors="replace") as f_in:
        first_line = f_in.readline()
        if not first_line:
            sys.stderr.write(f"ERROR: Input file appears to be empty: {in_path}\n")
            sys.exit(1)

        if "\t" in first_line and "," not in first_line:
            delimiter = "\t"
        elif "," in first_line and "\t" not in first_line:
            delimiter = ","
        else:
            delimiter = "\t" if first_line.count("\t") >= first_line.count(",") else ","

        f_in.seek(0)
        reader     = csv.DictReader(f_in, delimiter=delimiter)
        fieldnames = reader.fieldnames or []
        lower_map  = {name.lower(): name for name in fieldnames}

        required = ["species", "genus", "kingdom"]
        missing  = [c for c in required if c not in lower_map]
        if missing:
            sys.stderr.write(
                "ERROR: GBIF file must contain columns named "
                "'species', 'genus', and 'kingdom' (case-insensitive).\n"
                f"  Missing: {missing}\n"
                f"  Found: {fieldnames}\n"
            )
            sys.exit(1)

        col_species = lower_map["species"]
        col_genus   = lower_map["genus"]
        col_kingdom = lower_map["kingdom"]
        col_phylum  = lower_map.get("phylum")   # optional

        has_phylum = col_phylum is not None

        rows = 0
        for row in reader:
            rows += 1
            sp      = (row.get(col_species, "") or "").strip()
            genus   = (row.get(col_genus,   "") or "").strip()
            kingdom = (row.get(col_kingdom, "") or "").strip()
            phylum  = (row.get(col_phylum,  "") or "").strip() if has_phylum else ""

            if sp:
                species_for_search.add(sp)
                if kingdom:
                    # Later rows win for phylum if it differs, but kingdom is stable
                    if sp not in species_kingdom_pairs:
                        species_kingdom_pairs[sp] = (kingdom, phylum)
                    else:
                        # Fill in phylum if we previously had an empty one
                        existing_kd, existing_ph = species_kingdom_pairs[sp]
                        if not existing_ph and phylum:
                            species_kingdom_pairs[sp] = (existing_kd, phylum)
            else:
                if genus:
                    species_for_search.add(genus)

    # Write species search list (sorted)
    with open(species_search_path, "w", encoding="utf-8") as f_out:
        for name in sorted(species_for_search):
            f_out.write(name + "\n")

    # Write species <-> kingdom [<-> phylum] map
    with open(species_kingdom_path, "w", encoding="utf-8", newline="") as f_out:
        writer = csv.writer(f_out, delimiter="\t")
        for sp in sorted(species_kingdom_pairs):
            kingdom, phylum = species_kingdom_pairs[sp]
            writer.writerow([sp, kingdom, phylum])

    phylum_note = " (with phylum)" if has_phylum else " (no phylum column found)"
    sys.stderr.write(
        f"GBIF prep complete.\n"
        f"  Rows read from GBIF file  : {rows}\n"
        f"  Unique names for search   : {len(species_for_search)} -> {species_search_path}\n"
        f"  Unique species->kingdom{phylum_note}: {len(species_kingdom_pairs)} -> {species_kingdom_path}\n"
    )

if __name__ == "__main__":
    main()
