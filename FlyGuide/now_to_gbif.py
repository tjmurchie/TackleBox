#!/usr/bin/env python3
"""
FlyGuide NOW importer
======================

Create a FlyGuide/GBIF-like taxon CSV from NOW (New and Old Worlds Database
of Fossil Mammals) exported tables.

NOW has a browser UI and export support, but no documented stable public API
like PBDB, so this script is a robust IMPORTER for exported CSV/TSV/XLSX/ZIP
tables rather than a live API client. It accepts three input shapes:

  1. A single merged locality-species table   (--input)
  2. Separate species/localities/link tables  (--species-file/--localities-file/
                                                --locality-species-file)
  3. A directory or zip of exported files, auto-classified (--source-dir)

The output has the columns FlyGuide needs (`species`, `genus`, `kingdom`,
`phylum`) plus NOW metadata and optional native FlyGuide files:

  OUTPREFIX_species_search.txt
  OUTPREFIX_species_kingdom.tsv

This script intentionally uses only the Python standard library.

Examples
--------

# Single merged export:
python3 now_to_gbif.py \
  --input NOW_export.csv \
  --region eurasia \
  --period quaternary \
  --out Eurasia_quaternary_mammals_now_gbif.csv \
  --write-flyguide-files

# Three-table export:
python3 now_to_gbif.py \
  --species-file now_species.csv \
  --localities-file now_localities.csv \
  --locality-species-file now_locality_species.csv \
  --region holarctic \
  --period pleistocene \
  --out Holarctic_pleistocene_mammals_now_gbif.csv \
  --write-flyguide-files

# Directory/zip of exported files, auto-classified:
python3 now_to_gbif.py \
  --source-dir NOW_exports/ \
  --region northern_hemisphere \
  --period quaternary \
  --out NH_quaternary_mammals_now_gbif.csv \
  --write-flyguide-files

Notes
-----
NOW is a Cenozoic fossil MAMMAL database — every row this script emits has
kingdom=Animalia. A row in the output means: "this cleaned taxon/search
bucket appeared in a NOW export matching the requested filters." As with the
PBDB exporter, the goal is a large, auditable bucket of plausible NCBI search
names, not a definitive checklist.
"""

from __future__ import annotations

import argparse
import csv
import dataclasses
import json
import os
import re
import sys
import tempfile
import zipfile
from typing import Any, Dict, List, Optional, Sequence, Tuple

try:
    from _palaeo_tui import PalaeoProgress, print_header, print_step, print_done
except ImportError:
    PalaeoProgress = None  # type: ignore[assignment,misc]
    print_header = print_step = print_done = None  # type: ignore[assignment]

VERSION = "0.1.0-now"

NOW_EXPORT_INSTRUCTIONS = """\
How to export data from NOW (New and Old Worlds Database of Fossil Mammals)
----------------------------------------------------------------------------
NOW (https://nowdatabase.org, or the legacy helsinki.fi mirror) has no stable
public REST API, so exports must be produced from its browser search/results
UI:

1. Run a locality/species search for the region and time interval you want.
2. Use the results page's export/download option to save a table. NOW's
   export usually gives you one of:
     - A single merged locality-species table (species name + locality +
       country + coordinates + age columns in one file) -- use --input.
     - Separate "species", "localities", and "locality-species" (a.k.a.
       faunal list / link) tables -- use --species-file/--localities-file/
       --locality-species-file together.
   Some browser export flows bundle these into a single .zip -- point
   --source-dir at the extracted folder or the .zip itself and this script
   will classify the files by name and column signature.
3. NOW ages are usually reported in Ma (millions of years). If your export
   uses ka (thousands of years) or BP (years before present), pass
   --age-units ka or --age-units bp so ages are converted to Ma correctly.
4. Run this script against whatever you exported. Use --region/--period to
   filter to the area/interval you actually want (NOW exports are often
   broader than the search that produced them).
"""

# Keep these identical to pbdb_to_gbif.py's presets -- FlyGuide has no shared
# module between its standalone scripts, so the values are duplicated
# deliberately rather than imported, but they must stay in sync so a region
# or period preset means the same thing across every exporter.
PERIODS: Dict[str, Tuple[Optional[float], Optional[float], str]] = {
    "all": (None, None, "No age filter"),
    "quaternary": (0.0, 2.58, "Quaternary, ~0-2.58 Ma"),
    "pleistocene": (0.0117, 2.58, "Pleistocene, ~0.0117-2.58 Ma"),
    "holocene": (0.0, 0.0117, "Holocene, ~0-0.0117 Ma"),
    "late_pleistocene": (0.0117, 0.129, "Late Pleistocene, ~11.7-129 ka"),
    "middle_pleistocene": (0.129, 0.774, "Middle Pleistocene, ~129-774 ka"),
    "early_pleistocene": (0.774, 2.58, "Early Pleistocene, ~0.774-2.58 Ma"),
    "late_quaternary": (0.0, 0.129, "Late Quaternary, ~0-129 ka"),
    "lgm": (0.019, 0.0265, "Last Glacial Maximum, broad ~19-26.5 ka"),
    "last_glacial": (0.0117, 0.115, "Last glacial interval, broad ~11.7-115 ka"),
    "pliocene": (2.58, 5.333, "Pliocene, ~2.58-5.333 Ma"),
    "neogene": (2.58, 23.03, "Neogene excluding Quaternary-ish bucket, ~2.58-23.03 Ma"),
    "cenozoic": (0.0, 66.0, "Cenozoic, ~0-66 Ma"),
}

REGION_BBOXES: Dict[str, Tuple[float, float, float, float, str]] = {
    "global": (-180.0, -90.0, 180.0, 90.0, "global bbox"),
    "world": (-180.0, -90.0, 180.0, 90.0, "global bbox"),
    "northern_hemisphere": (-180.0, 0.0, 180.0, 90.0, "Northern Hemisphere bbox"),
    "north_hemisphere": (-180.0, 0.0, 180.0, 90.0, "Northern Hemisphere bbox"),
    "southern_hemisphere": (-180.0, -90.0, 180.0, 0.0, "Southern Hemisphere bbox"),
    "holarctic": (-180.0, 20.0, 180.0, 90.0, "coarse Holarctic-ish bbox, north of ~20 N"),
    "nearctic": (-170.0, 15.0, -50.0, 85.0, "coarse Nearctic/North America bbox"),
    "palearctic": (-25.0, 20.0, 180.0, 85.0, "coarse Palearctic/Eurasia+N Africa bbox"),
    "north_america": (-170.0, 5.0, -50.0, 85.0, "North America broad bbox"),
    "canada": (-142.0, 41.0, -52.0, 84.0, "Canada broad bbox"),
    "usa": (-170.0, 18.0, -65.0, 72.0, "United States incl. Alaska/Hawaii broad bbox"),
    "united_states": (-170.0, 18.0, -65.0, 72.0, "United States incl. Alaska/Hawaii broad bbox"),
    "eurasia": (-25.0, 0.0, 180.0, 85.0, "Eurasia broad bbox"),
    "europe": (-25.0, 34.0, 45.0, 72.0, "Europe broad bbox"),
    "asia": (25.0, -10.0, 180.0, 85.0, "Asia broad bbox"),
    "africa": (-20.0, -35.0, 55.0, 38.0, "Africa broad bbox"),
    "south_america": (-82.0, -56.0, -34.0, 13.0, "South America broad bbox"),
    "australia": (110.0, -45.0, 155.0, -10.0, "Australia broad bbox"),
    "tanzania": (29.0, -12.5, 41.5, -0.5, "Tanzania broad bbox"),
    "east_africa": (21.0, -12.5, 52.0, 18.5, "East Africa broad bbox"),
    "beringia": (-180.0, 50.0, -120.0, 75.0, "Beringia broad bbox, Alaska/Yukon side only"),
}
REGION_ALIASES = {k.replace("_", " "): k for k in REGION_BBOXES}
REGION_ALIASES.update({k.replace("_", "-"): k for k in REGION_BBOXES})

# Same name-cleaning vocabulary as pbdb_to_gbif.py, duplicated for the same
# reason as the presets above.
QUALIFIER_RE = re.compile(r"\b(cf|aff|nr|near|ex gr|gr|group|complex)\.?\b", re.IGNORECASE)
SP_RE = re.compile(r"\b(spp|sp|species|indet|indeterminate|undiff|undifferentiated)\.?\b", re.IGNORECASE)
TYPE_RE = re.compile(r"(?:-type\b|\btype\b|\bmorphotype\b)", re.IGNORECASE)
BAD_CHARS_RE = re.compile(r"[\[\]{}<>≈~;:=+*#@!]")
HYBRID_RE = re.compile(r"\b[x×]\b")
BAD_SINGLE_WORDS = {"unknown", "indeterminate", "undetermined", "unidentified", "other", "fossil", "sample"}
HIGHER_TAXON_SUFFIXES = (
    "idae", "inae", "ini", "aceae", "ales", "iformes", "phyta", "mycota",
    "oda", "ata", "opsida", "oidea", "ina", "inae", "inae/",
)


@dataclasses.dataclass
class CleanedName:
    search_name: str
    genus: str
    rank: str
    canonical_binomial: str
    cleaning_status: str


@dataclasses.dataclass
class RejectedName:
    source_id: str
    source_name: str
    reason: str
    context: str = ""


class NowImportError(RuntimeError):
    pass


def normalize_key(s: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", (s or "").lower())


def get_first(row: Dict[str, Any], names: Sequence[str]) -> str:
    for name in names:
        val = row.get(name)
        if val not in (None, ""):
            return str(val).strip()
    nmap = {normalize_key(k): k for k in row}
    for name in names:
        k = nmap.get(normalize_key(name))
        if k and row.get(k) not in (None, ""):
            return str(row[k]).strip()
    return ""


def to_float(value: Any) -> Optional[float]:
    if value in (None, "", "NA", "N/A", "null"):
        return None
    try:
        return float(str(value).strip())
    except ValueError:
        return None


def age_to_ma(value: Optional[float], units: str) -> Optional[float]:
    if value is None:
        return None
    if units == "ma":
        return value
    if units == "ka":
        return value / 1000.0
    if units == "bp":
        return value / 1_000_000.0
    raise NowImportError(f"Unknown --age-units {units!r}")


def parse_bbox(raw: str) -> Tuple[float, float, float, float, str]:
    parts = [p.strip() for p in raw.split(",")]
    if len(parts) != 4:
        raise argparse.ArgumentTypeError("bbox must be minLon,minLat,maxLon,maxLat")
    try:
        min_lon, min_lat, max_lon, max_lat = [float(p) for p in parts]
    except ValueError as exc:
        raise argparse.ArgumentTypeError("bbox values must be numeric") from exc
    return min_lon, min_lat, max_lon, max_lat, "custom bbox"


def resolve_region(region: str, bbox: Optional[str]) -> Tuple[str, Tuple[float, float, float, float, str]]:
    if bbox:
        return "custom_bbox", parse_bbox(bbox)
    key = (region or "global").strip().lower().replace("-", "_").replace(" ", "_")
    key = REGION_ALIASES.get(key, key)
    if key not in REGION_BBOXES:
        valid = ", ".join(sorted(REGION_BBOXES))
        raise SystemExit(f"Unknown region '{region}'. Valid presets: {valid}; or use --bbox.")
    return key, REGION_BBOXES[key]


def resolve_period(period: str, age_young_ma: Optional[float], age_old_ma: Optional[float]) -> Tuple[str, Optional[float], Optional[float], str]:
    if age_young_ma is not None or age_old_ma is not None:
        if age_young_ma is None or age_old_ma is None:
            raise SystemExit("Use both --age-young-ma and --age-old-ma for custom ages.")
        if age_young_ma < 0 or age_old_ma < 0 or age_young_ma > age_old_ma:
            raise SystemExit("Custom ages must satisfy 0 <= young <= old, in Ma.")
        return "custom", age_young_ma, age_old_ma, f"custom {age_young_ma}-{age_old_ma} Ma"
    key = (period or "quaternary").strip().lower().replace(" ", "_").replace("-", "_")
    if key not in PERIODS:
        valid = ", ".join(sorted(PERIODS))
        raise SystemExit(f"Unknown period '{period}'. Valid presets: {valid}; or use --age-young-ma/--age-old-ma.")
    young, old, desc = PERIODS[key]
    return key, young, old, desc


def clean_taxon_name(raw: str, mode: str = "binomial", keep_genus_only: bool = True,
                      keep_higher_taxa: bool = False, drop_morphotypes: bool = False) -> Optional[CleanedName]:
    """Identical cleaning logic to pbdb_to_gbif.py's clean_taxon_name (minus
    the slash-name splitting option, which NOW exports rarely need) -- kept
    in sync deliberately so the same source name collapses the same way
    regardless of which FlyGuide exporter produced it."""
    name = (raw or "").strip()
    if not name:
        return None
    original = name
    name = name.replace("†", "")
    name = name.replace("_", " ")
    name = BAD_CHARS_RE.sub(" ", name)
    name = re.sub(r"\s+", " ", name).strip()
    if not name:
        return None
    lowered = name.lower()
    if lowered in BAD_SINGLE_WORDS:
        return None
    if HYBRID_RE.search(name):
        name = HYBRID_RE.split(name)[0].strip()
    if drop_morphotypes and TYPE_RE.search(name):
        return None
    morph = bool(TYPE_RE.search(name))
    name = QUALIFIER_RE.sub(" ", name)
    had_sp_marker = bool(SP_RE.search(name))
    name = SP_RE.sub(" ", name)
    name = re.sub(r"\([^)]*\)", " ", name)
    name = re.sub(r"\b[a-z]\.?\s+(?=[A-Z][a-z])", " ", name)
    name = re.sub(r"\s+", " ", name).strip()
    if not name:
        return None
    tokens = [t for t in re.split(r"\s+", name) if t]
    if not tokens:
        return None
    genus = tokens[0]
    if not re.match(r"^[A-Z][A-Za-z-]+$", genus):
        if not keep_higher_taxa:
            return None
    if len(tokens) == 1:
        one = tokens[0]
        rank = "genus_or_higher"
        if one.lower().endswith(HIGHER_TAXON_SUFFIXES) and not keep_higher_taxa:
            return None
        if not keep_genus_only and not keep_higher_taxa:
            return None
        status = "genus_only" if had_sp_marker else "single_name"
        return CleanedName(search_name=one, genus=one, rank=rank, canonical_binomial="", cleaning_status=status)
    genus = tokens[0]
    epithet = tokens[1]
    if not re.match(r"^[a-z][a-z-]+$", epithet):
        if keep_higher_taxa:
            search = " ".join(tokens[:2]) if mode != "genus" else genus
            return CleanedName(search, genus, "higher_or_uncertain", "", "higher_or_uncertain")
        return None
    binomial = f"{genus} {epithet}"
    trinomial = binomial
    rank = "species"
    if len(tokens) >= 3 and re.match(r"^[a-z][a-z-]+$", tokens[2]):
        trinomial = f"{binomial} {tokens[2]}"
        rank = "subspecies_or_trinomial"
    if mode == "as-is":
        search = trinomial if rank == "subspecies_or_trinomial" else binomial
    elif mode == "trinomial":
        search = trinomial
    elif mode == "genus":
        search = genus
    else:
        search = binomial
    status_bits = ["cleaned"]
    if original != search:
        status_bits.append("collapsed" if search == binomial and trinomial != binomial else "normalized")
    if morph:
        status_bits.append("morphotype_kept")
    return CleanedName(search_name=search, genus=genus, rank=rank, canonical_binomial=binomial, cleaning_status=";".join(status_bits))


# --- Reading exported tables -------------------------------------------------

def sniff_delimiter(sample: str) -> str:
    try:
        return csv.Sniffer().sniff(sample, delimiters=",\t;").delimiter
    except csv.Error:
        return "\t" if sample.count("\t") > sample.count(",") else ","


def read_delimited(path: str) -> List[Dict[str, str]]:
    with open(path, "r", encoding="utf-8-sig", errors="replace", newline="") as fh:
        sample = fh.read(4096)
        fh.seek(0)
        delim = sniff_delimiter(sample)
        return list(csv.DictReader(fh, delimiter=delim))


def read_xlsx(path: str) -> List[Dict[str, str]]:
    """Minimal stdlib-only XLSX reader for a single flat sheet -- enough for
    a NOW browser export, not a general spreadsheet parser."""
    import xml.etree.ElementTree as ET

    ns = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    with zipfile.ZipFile(path) as zf:
        shared: List[str] = []
        if "xl/sharedStrings.xml" in zf.namelist():
            root = ET.fromstring(zf.read("xl/sharedStrings.xml"))
            for si in root.findall("m:si", ns):
                texts = si.findall(".//m:t", ns)
                shared.append("".join(t.text or "" for t in texts))
        sheet_names = [n for n in zf.namelist() if n.startswith("xl/worksheets/sheet") and n.endswith(".xml")]
        if not sheet_names:
            raise NowImportError(f"No worksheet found in {path}")
        sheet_xml = zf.read(sorted(sheet_names)[0])
        root = ET.fromstring(sheet_xml)
        rows_out: List[List[str]] = []
        for row_el in root.findall(".//m:sheetData/m:row", ns):
            cells: Dict[int, str] = {}
            for c_el in row_el.findall("m:c", ns):
                ref = c_el.get("r", "")
                col_letters = re.match(r"[A-Z]+", ref)
                col_idx = 0
                if col_letters:
                    for ch in col_letters.group():
                        col_idx = col_idx * 26 + (ord(ch) - ord("A") + 1)
                    col_idx -= 1
                v_el = c_el.find("m:v", ns)
                text = v_el.text if v_el is not None else ""
                if c_el.get("t") == "s" and text is not None:
                    try:
                        text = shared[int(text)]
                    except (ValueError, IndexError):
                        pass
                cells[col_idx] = text or ""
            if cells:
                width = max(cells) + 1
                rows_out.append([cells.get(i, "") for i in range(width)])
    if not rows_out:
        return []
    header = rows_out[0]
    return [dict(zip(header, r)) for r in rows_out[1:]]


def read_table(path: str) -> List[Dict[str, str]]:
    lower = path.lower()
    if lower.endswith(".xlsx"):
        return read_xlsx(path)
    return read_delimited(path)


def classify_export_file(path: str, rows: List[Dict[str, str]]) -> str:
    """Best-effort classification of one exported file as species/localities/
    link/merged, by filename hint first, falling back to column signature."""
    base = os.path.basename(path).lower()
    if "local" in base:
        return "localities"
    if "link" in base or "faunal" in base or "occ" in base:
        return "link"
    if "spec" in base or "taxa" in base:
        return "species"
    if not rows:
        return "merged"
    cols = {normalize_key(c) for c in rows[0].keys()}
    has_species = bool(cols & {normalize_key(x) for x in ("species", "speciesname", "taxon", "name")})
    has_locality = bool(cols & {normalize_key(x) for x in ("locality", "localityname", "site")})
    has_coords = bool(cols & {normalize_key(x) for x in ("lat", "latitude")})
    has_age = bool(cols & {normalize_key(x) for x in ("minma", "maxma", "age", "minage", "maxage")})
    if has_species and (has_locality or has_coords or has_age):
        return "merged"
    if has_locality and has_coords:
        return "localities"
    if has_species and not has_locality and not has_coords:
        return "species"
    return "link"


def gather_source_dir(source_dir: str, verbose: bool = True) -> List[Dict[str, str]]:
    """Load --source-dir (a directory, or a zip that is extracted first) and
    merge whatever species/localities/link/merged tables it contains into a
    single list of merged-style rows."""
    work_dir = source_dir
    tmp_ctx = None
    if os.path.isfile(source_dir) and source_dir.lower().endswith(".zip"):
        tmp_ctx = tempfile.TemporaryDirectory(prefix="now_to_gbif_")
        with zipfile.ZipFile(source_dir) as zf:
            zf.extractall(tmp_ctx.name)
        work_dir = tmp_ctx.name
    try:
        files = []
        for root, _dirs, names in os.walk(work_dir):
            for name in names:
                if name.lower().endswith((".csv", ".tsv", ".txt", ".xlsx")):
                    files.append(os.path.join(root, name))
        if not files:
            raise NowImportError(f"No CSV/TSV/XLSX files found under {source_dir}")
        buckets: Dict[str, List[Dict[str, str]]] = {"species": [], "localities": [], "link": [], "merged": []}
        for f in files:
            try:
                rows = read_table(f)
            except Exception as exc:  # noqa: BLE001 - a single bad file shouldn't abort the whole import
                if verbose:
                    print(f"  [WARNING] Could not read {f}: {exc}", file=sys.stderr)
                continue
            kind = classify_export_file(f, rows)
            buckets[kind].extend(rows)
            if verbose:
                print(f"  Classified {os.path.basename(f)} as '{kind}' ({len(rows)} rows)", file=sys.stderr)
        if buckets["merged"]:
            return buckets["merged"]
        return join_species_localities_link(buckets["species"], buckets["localities"], buckets["link"])
    finally:
        if tmp_ctx is not None:
            tmp_ctx.cleanup()


def join_species_localities_link(species_rows: List[Dict[str, str]], locality_rows: List[Dict[str, str]],
                                  link_rows: List[Dict[str, str]]) -> List[Dict[str, str]]:
    """Join separate species/localities/link tables into merged-style rows.
    Joins on the first shared-looking ID column found in each pair; falls
    back to species/locality name matching if no ID column is present."""
    def index_by(rows: List[Dict[str, str]], id_aliases: Sequence[str], name_aliases: Sequence[str]) -> Dict[str, Dict[str, str]]:
        idx: Dict[str, Dict[str, str]] = {}
        for row in rows:
            key = get_first(row, id_aliases) or get_first(row, name_aliases)
            if key:
                idx[key] = row
        return idx

    species_idx = index_by(
        species_rows,
        ["SpeciesID", "species_id", "SPECIES_ID", "id", "ID"],
        ["Species", "species", "SPECIES", "taxon", "name"],
    )
    locality_idx = index_by(
        locality_rows,
        ["LocalityID", "locality_id", "LOCALITY_ID", "id", "ID"],
        ["Locality", "locality", "LOCALITY", "site"],
    )

    if not link_rows:
        # No explicit link table -- assume species rows already carry
        # locality info directly (a "merged" table misclassified as species).
        return species_rows

    merged: List[Dict[str, str]] = []
    for link in link_rows:
        sp_key = get_first(link, ["SpeciesID", "species_id", "SPECIES_ID", "Species", "species", "SPECIES"])
        loc_key = get_first(link, ["LocalityID", "locality_id", "LOCALITY_ID", "Locality", "locality", "LOCALITY"])
        sp = species_idx.get(sp_key, {})
        loc = locality_idx.get(loc_key, {})
        row = dict(link)
        for k, v in sp.items():
            row.setdefault(k, v)
        for k, v in loc.items():
            row.setdefault(k, v)
        merged.append(row)
    return merged


def load_input_rows(args: argparse.Namespace) -> List[Dict[str, str]]:
    if args.fixture:
        with open(args.fixture, "r", encoding="utf-8") as fh:
            data = json.load(fh)
        return data if isinstance(data, list) else (data.get("records") or [])
    if args.input:
        return read_table(args.input)
    if args.source_dir:
        return gather_source_dir(args.source_dir, verbose=args.verbose)
    if args.species_file or args.localities_file or args.locality_species_file:
        species_rows = read_table(args.species_file) if args.species_file else []
        locality_rows = read_table(args.localities_file) if args.localities_file else []
        link_rows = read_table(args.locality_species_file) if args.locality_species_file else []
        return join_species_localities_link(species_rows, locality_rows, link_rows)
    raise SystemExit("Provide one of --input, --source-dir, or --species-file/--localities-file/--locality-species-file.")


# --- Aggregation --------------------------------------------------------------

def bbox_contains(lat: Optional[float], lon: Optional[float], bbox: Tuple[float, float, float, float, str]) -> bool:
    if lat is None or lon is None:
        return True  # don't silently drop records lacking coordinates here; --require-coords handles that explicitly
    min_lon, min_lat, max_lon, max_lat, _desc = bbox
    return min_lat <= lat <= max_lat and min_lon <= lon <= max_lon


def status_is_extinct_or_unknown(row: Dict[str, str]) -> bool:
    raw = get_first(row, ["Status", "status", "STATUS", "extant_extinct"])
    return True if not raw else raw.strip().lower() != "extant"


def aggregate_rows(rows: List[Dict[str, str]], args: argparse.Namespace, region_key: str, bbox: Tuple[float, float, float, float, str],
                    period_key: str, young_ma: Optional[float], old_ma: Optional[float],
                    region_desc: str, period_desc: str) -> Tuple[List[Dict[str, str]], List[RejectedName]]:
    grouped: Dict[str, Dict[str, Any]] = {}
    rejected: List[RejectedName] = []
    total = len(rows)
    dash = PalaeoProgress(total, label="records", verbose=args.verbose) if PalaeoProgress else None
    if dash and total:
        dash.set_current(f"Processing {total:,} NOW export rows", "starting")
    report_every = max(1, total // 200)

    for i, row in enumerate(rows):
        source_id = get_first(row, ["SpeciesID", "species_id", "LIDNUM", "id", "ID"])
        raw_name = get_first(row, ["Species", "species", "SPECIES", "taxon", "SciName", "name"])
        if not raw_name:
            rejected.append(RejectedName(source_id=source_id, source_name="", reason="no_species_name", context=json.dumps(row, ensure_ascii=False)[:500]))
            continue

        lat = to_float(get_first(row, ["Latitude", "latitude", "LAT", "lat"]))
        lon = to_float(get_first(row, ["Longitude", "longitude", "LONG", "LNG", "lon", "lng"]))
        if args.require_coords and (lat is None or lon is None):
            rejected.append(RejectedName(source_id=source_id, source_name=raw_name, reason="missing_coordinates"))
            continue
        if not bbox_contains(lat, lon, bbox):
            continue

        min_ma = age_to_ma(to_float(get_first(row, ["MinMa", "min_ma", "MIN_AGE", "AGE_MIN", "min_age"])), args.age_units)
        max_ma = age_to_ma(to_float(get_first(row, ["MaxMa", "max_ma", "MAX_AGE", "AGE_MAX", "max_age"])), args.age_units)
        if args.require_age and min_ma is None and max_ma is None:
            rejected.append(RejectedName(source_id=source_id, source_name=raw_name, reason="missing_age"))
            continue
        if young_ma is not None and old_ma is not None and min_ma is not None and max_ma is not None:
            # Record's own [min_ma, max_ma] age range must overlap the requested [young_ma, old_ma] window.
            if max_ma < young_ma or min_ma > old_ma:
                continue

        if args.status != "all":
            is_extinct = status_is_extinct_or_unknown(row)
            if args.status == "extant" and is_extinct:
                continue
            if args.status == "extinct" and not is_extinct:
                continue

        cleaned = clean_taxon_name(
            raw_name,
            mode=args.ncbi_name_mode,
            keep_genus_only=not args.drop_genus_only,
            keep_higher_taxa=args.keep_higher_taxa,
            drop_morphotypes=args.drop_morphotypes,
        )
        if not cleaned:
            rejected.append(RejectedName(source_id=source_id, source_name=raw_name, reason="could_not_clean_name"))
            continue

        key = cleaned.search_name
        country = get_first(row, ["Country", "country", "COUNTRY"])
        locality = get_first(row, ["Locality", "locality", "LOCALITY"])
        if key not in grouped:
            grouped[key] = {
                "species": cleaned.search_name,
                "genus": cleaned.genus,
                "kingdom": "Animalia",  # NOW is a fossil mammal database -- always Animalia
                "phylum": "Chordata",
                "scientificName": cleaned.search_name,
                "taxonRank": cleaned.rank,
                "taxonKey": source_id,
                "source": "NOW",
                "source_database": "New and Old Worlds Database of Fossil Mammals",
                "now_source_name": raw_name,
                "now_localities": set(),
                "now_countries": set(),
                "now_min_ma_youngest": min_ma,
                "now_max_ma_oldest": max_ma,
                "now_record_count": 0,
                "region": region_key,
                "region_description": region_desc,
                "period": period_key,
                "period_description": period_desc,
                "ncbi_search_name": cleaned.search_name,
                "canonical_binomial": cleaned.canonical_binomial,
                "name_cleaning_status": cleaned.cleaning_status,
            }
        g = grouped[key]
        g["now_record_count"] += 1
        if locality:
            g["now_localities"].add(locality)
        if country:
            g["now_countries"].add(country)
        if min_ma is not None:
            g["now_min_ma_youngest"] = min_ma if g["now_min_ma_youngest"] is None else min(g["now_min_ma_youngest"], min_ma)
        if max_ma is not None:
            g["now_max_ma_oldest"] = max_ma if g["now_max_ma_oldest"] is None else max(g["now_max_ma_oldest"], max_ma)
        if dash and (i % report_every == 0 or i == total - 1):
            dash.set_current(f"{key}", f"{len(grouped):,} unique taxa", completed=i + 1)
    if dash:
        dash.finish()

    out_rows: List[Dict[str, str]] = []
    for key in sorted(grouped):
        g = grouped[key]
        for set_col in ("now_localities", "now_countries"):
            g[set_col] = ";".join(sorted(g[set_col]))
        for col in ("now_min_ma_youngest", "now_max_ma_oldest"):
            g[col] = "" if g[col] is None else f"{g[col]:.6g}"
        for k, v in list(g.items()):
            if not isinstance(v, str):
                g[k] = str(v)
        out_rows.append(g)
    return out_rows, rejected


# --- Output -------------------------------------------------------------------

FIELDNAMES = [
    "species", "genus", "kingdom", "phylum", "scientificName", "taxonRank", "taxonKey",
    "source", "source_database", "now_source_name", "now_localities", "now_countries",
    "now_min_ma_youngest", "now_max_ma_oldest", "now_record_count",
    "region", "region_description", "period", "period_description",
    "ncbi_search_name", "canonical_binomial", "name_cleaning_status",
]


def write_csv(path: str, rows: List[Dict[str, str]]) -> None:
    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDNAMES, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_rejected(path: str, rejected: List[RejectedName]) -> None:
    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["source_id", "source_name", "reason", "context"])
        for r in rejected:
            writer.writerow([r.source_id, r.source_name, r.reason, r.context])


def write_flyguide_files(out_csv: str, rows: List[Dict[str, str]], out_prefix: Optional[str] = None) -> Tuple[str, str]:
    prefix = out_prefix or os.path.splitext(os.path.basename(out_csv))[0]
    species_path = f"{prefix}_species_search.txt"
    kingdom_path = f"{prefix}_species_kingdom.tsv"
    seen: Dict[str, Tuple[str, str]] = {}
    for row in rows:
        name = row.get("ncbi_search_name") or row.get("species")
        if not name:
            continue
        if name not in seen:
            seen[name] = (row.get("kingdom", ""), row.get("phylum", ""))
    with open(species_path, "w", encoding="utf-8") as fh:
        for name in sorted(seen):
            fh.write(name + "\n")
    with open(kingdom_path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        for name in sorted(seen):
            kingdom, phylum = seen[name]
            writer.writerow([name, kingdom, phylum])
    return species_path, kingdom_path


def run(args: argparse.Namespace) -> Dict[str, Any]:
    region_key, bbox = resolve_region(args.region, args.bbox)
    period_key, young_ma, old_ma, period_desc = resolve_period(args.period, args.age_young_ma, args.age_old_ma)
    region_desc = bbox[4]

    if print_header and args.verbose:
        period_label = f"{period_desc} ({young_ma}–{old_ma} Ma)" if young_ma is not None else period_desc
        print_header(
            f"FlyGuide NOW importer  v{VERSION}",
            [("Region", region_desc), ("Period", period_label), ("Output", args.out)],
        )

    input_rows = load_input_rows(args)
    rows, rejected = aggregate_rows(input_rows, args, region_key, bbox, period_key, young_ma, old_ma, region_desc, period_desc)
    if args.min_records > 1:
        rows = [r for r in rows if int(r.get("now_record_count", "0") or 0) >= args.min_records]

    write_csv(args.out, rows)
    rejected_path = os.path.splitext(args.out)[0] + ".rejected.tsv"
    write_rejected(rejected_path, rejected)
    outputs: Dict[str, str] = {"csv": args.out, "rejected": rejected_path}
    if args.write_flyguide_files:
        sp, kd = write_flyguide_files(args.out, rows, args.out_prefix)
        outputs["species_search"] = sp
        outputs["species_kingdom"] = kd

    summary = {
        "version": VERSION,
        "source": "NOW",
        "region": region_key,
        "region_bbox": list(bbox[:4]),
        "period": period_key,
        "period_bounds_ma": [young_ma, old_ma],
        "records_read": len(input_rows),
        "taxa_written": len(rows),
        "names_rejected": len(rejected),
        "outputs": outputs,
    }
    summary_path = os.path.splitext(args.out)[0] + ".summary.json"
    with open(summary_path, "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)
    outputs["summary"] = summary_path

    if print_done and args.verbose:
        counts = [
            ("Records read:", f"{len(input_rows):,}"),
            ("Unique search taxa written:", f"{len(rows):,}  →  {args.out}"),
            ("Rejected/noisy names:", f"{len(rejected):,}  →  {rejected_path}"),
        ]
        print_done("NOW import complete", counts)
    elif args.verbose:
        print("NOW import complete.", file=sys.stderr)
        print(f"  Records read: {len(input_rows)}", file=sys.stderr)
        print(f"  Unique search taxa written: {len(rows)} -> {args.out}", file=sys.stderr)
        print(f"  Rejected/noisy names: {len(rejected)} -> {rejected_path}", file=sys.stderr)
    return summary


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Import NOW fossil-mammal exported tables into a FlyGuide/GBIF-like CSV.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--input", help="Single merged locality-species export (CSV/TSV/XLSX)")
    p.add_argument("--species-file", help="NOW 'species' export table")
    p.add_argument("--localities-file", help="NOW 'localities' export table")
    p.add_argument("--locality-species-file", help="NOW locality-species/faunal-list link table")
    p.add_argument("--source-dir", help="Directory or .zip of NOW export files; auto-classified")
    p.add_argument("--explain-now-export", action="store_true", help="Print NOW export instructions and exit")

    p.add_argument("--region", default="northern_hemisphere", help="Built-in broad region preset")
    p.add_argument("--bbox", help="Custom minLon,minLat,maxLon,maxLat; overrides --region")
    p.add_argument("--period", default="quaternary", help="Broad time period preset")
    p.add_argument("--age-young-ma", type=float, help="Custom youngest age in Ma")
    p.add_argument("--age-old-ma", type=float, help="Custom oldest age in Ma")
    p.add_argument("--age-units", choices=["ma", "ka", "bp"], default="ma", help="Units of age columns in the export (converted to Ma internally)")
    p.add_argument("--require-coords", action="store_true", help="Drop records lacking latitude/longitude instead of keeping them")
    p.add_argument("--require-age", action="store_true", help="Drop records lacking any age value instead of keeping them")
    p.add_argument("--status", choices=["all", "extinct", "extant"], default="all", help="Best-effort extant/extinct filter using the export's status column")

    p.add_argument("--ncbi-name-mode", choices=["binomial", "trinomial", "as-is", "genus"], default="binomial", help="How to collapse names for NCBI searching")
    p.add_argument("--drop-genus-only", action="store_true", help="Drop genus-only names such as Canis sp.")
    p.add_argument("--keep-higher-taxa", action="store_true", help="Keep higher-taxon buckets; usually noisy for NCBI")
    p.add_argument("--drop-morphotypes", action="store_true", help="Drop -type/morphotype names")
    p.add_argument("--min-records", type=int, default=1, help="Require at least this many NOW records per cleaned taxon")

    p.add_argument("--out", help="Output GBIF-like CSV")
    p.add_argument("--out-prefix", help="Prefix for FlyGuide files; default derives from --out")
    p.add_argument("--write-flyguide-files", action="store_true", help="Also write *_species_search.txt and *_species_kingdom.tsv")

    p.add_argument("--fixture", help="Offline test fixture JSON with NOW-like records (list of dicts)")
    p.add_argument("--quiet", dest="verbose", action="store_false", help="Reduce progress output")
    p.add_argument("--verbose", dest="verbose", action="store_true")
    p.set_defaults(verbose=True)
    p.add_argument("--list-regions", action="store_true", help="List region presets and exit")
    p.add_argument("--list-periods", action="store_true", help="List period presets and exit")
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    p = build_parser()
    args = p.parse_args(argv)
    if args.explain_now_export:
        print(NOW_EXPORT_INSTRUCTIONS)
        return 0
    if args.list_regions:
        for k, v in sorted(REGION_BBOXES.items()):
            print(f"{k}\t{v[0]},{v[1]},{v[2]},{v[3]}\t{v[4]}")
        return 0
    if args.list_periods:
        for k, (young, old, desc) in sorted(PERIODS.items()):
            print(f"{k}\t{young}\t{old}\t{desc}")
        return 0
    if not args.out:
        p.error("--out is required unless --explain-now-export/--list-regions/--list-periods")
    if not any([args.input, args.source_dir, args.species_file, args.localities_file, args.locality_species_file, args.fixture]):
        p.error("Provide --input, --source-dir, --species-file/--localities-file/--locality-species-file, or --fixture")
    run(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
