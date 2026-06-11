#!/usr/bin/env python3
"""
FlyGuide PBDB exporter
======================

Create a FlyGuide/GBIF-like taxon CSV from Paleobiology Database (PBDB)
fossil occurrence data.

The output has the columns FlyGuide needs (`species`, `genus`, `kingdom`,
`phylum`) plus PBDB metadata and optional native FlyGuide files:

  OUTPREFIX_species_search.txt
  OUTPREFIX_species_kingdom.tsv

This script intentionally uses only the Python standard library.

Examples
--------

# Quaternary animal taxa from the Northern Hemisphere:
python3 pbdb_to_gbif.py \
  --region northern_hemisphere \
  --period quaternary \
  --organisms animals \
  --out NH_quaternary_animals_pbdb_gbif.csv \
  --write-flyguide-files

# Eurasian Pleistocene mammals, with binomial NCBI names:
python3 pbdb_to_gbif.py \
  --region eurasia \
  --period pleistocene \
  --organisms mammals \
  --ncbi-name-mode binomial \
  --out Eurasia_pleistocene_mammals_pbdb_gbif.csv \
  --write-flyguide-files

# A custom taxon/clade search:
python3 pbdb_to_gbif.py \
  --region africa \
  --period quaternary \
  --base-name Bovidae \
  --out Africa_quaternary_bovidae_pbdb_gbif.csv

Notes
-----
PBDB is a fossil occurrence database, not an NCBI-ready checklist. A row in the
output means: "this cleaned taxon/search bucket appeared in PBDB fossil
occurrences matching the requested filters." For ancient DNA reference hunting,
this is often exactly what you want: a large, auditable bucket of plausible names.
"""

from __future__ import annotations

import argparse
import csv
import dataclasses
import hashlib
import json
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import defaultdict
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

try:
    from _palaeo_tui import PalaeoProgress, print_header, print_step, print_done
except ImportError:
    PalaeoProgress = None  # type: ignore[assignment,misc]
    print_header = print_step = print_done = None  # type: ignore[assignment]

VERSION = "0.1.0-pbdb"
DEFAULT_BASE_URL = "https://paleobiodb.org/data1.2"

# PBDB ages use Ma.  Periods are deliberately broad buckets useful for aDNA/protein
# reference list generation rather than fine stratigraphic analysis.
PERIODS: Dict[str, Tuple[Optional[float], Optional[float], str, Optional[str]]] = {
    # key: (young_ma, old_ma, description, PBDB interval name when safe)
    "all": (None, None, "No age filter", None),
    "quaternary": (0.0, 2.58, "Quaternary, ~0-2.58 Ma", "Quaternary"),
    "pleistocene": (0.0117, 2.58, "Pleistocene, ~0.0117-2.58 Ma", "Pleistocene"),
    "holocene": (0.0, 0.0117, "Holocene, ~0-0.0117 Ma", "Holocene"),
    "late_pleistocene": (0.0117, 0.129, "Late Pleistocene, ~11.7-129 ka", "Late Pleistocene"),
    "middle_pleistocene": (0.129, 0.774, "Middle Pleistocene, ~129-774 ka", "Middle Pleistocene"),
    "early_pleistocene": (0.774, 2.58, "Early Pleistocene, ~0.774-2.58 Ma", "Early Pleistocene"),
    "late_quaternary": (0.0, 0.129, "Late Quaternary, ~0-129 ka", None),
    "lgm": (0.019, 0.0265, "Last Glacial Maximum, broad ~19-26.5 ka", None),
    "last_glacial": (0.0117, 0.115, "Last glacial interval, broad ~11.7-115 ka", None),
    "pliocene": (2.58, 5.333, "Pliocene, ~2.58-5.333 Ma", "Pliocene"),
    "neogene": (2.58, 23.03, "Neogene excluding Quaternary-ish bucket, ~2.58-23.03 Ma", "Neogene"),
    "cenozoic": (0.0, 66.0, "Cenozoic, ~0-66 Ma", "Cenozoic"),
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

# PBDB base_name presets. These are deliberately broad because the goal is
# generating search buckets, not making final biodiversity estimates.
ORGANISM_BASE_NAMES: Dict[str, List[str]] = {
    "all": [],
    "animals": ["Animalia"],
    "metazoa": ["Animalia"],
    "vertebrates": ["Vertebrata"],
    "chordates": ["Chordata"],
    "mammals": ["Mammalia"],
    "birds": ["Aves"],
    "fish": ["Actinopterygii", "Chondrichthyes", "Sarcopterygii", "Agnatha"],
    "reptiles_amphibians": ["Reptilia", "Amphibia"],
    "reptiles": ["Reptilia"],
    "amphibians": ["Amphibia"],
    "invertebrates": ["Animalia"],  # filtered locally to remove Chordata
    "arthropods": ["Arthropoda"],
    "molluscs": ["Mollusca"],
    "plants": ["Plantae"],
    "fungi": ["Fungi"],
    "diatoms": ["Bacillariophyta"],
    "protists": ["Foraminifera", "Radiolaria", "Diatomea"],
}

ANIMAL_PHYLA = {
    "Chordata", "Arthropoda", "Mollusca", "Annelida", "Brachiopoda", "Bryozoa",
    "Echinodermata", "Cnidaria", "Porifera", "Hemichordata", "Nematoda",
    "Platyhelminthes", "Priapulida", "Onychophora", "Tardigrada", "Rotifera",
}
PLANT_PHYLA_HINTS = {
    "Tracheophyta", "Streptophyta", "Bryophyta", "Anthophyta", "Magnoliophyta",
    "Pinophyta", "Lycopodiophyta", "Pteridophyta", "Charophyta", "Chlorophyta",
}
FUNGI_PHYLA_HINTS = {"Ascomycota", "Basidiomycota", "Glomeromycota", "Chytridiomycota"}

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

class PBDBError(RuntimeError):
    pass

class PBDBClient:
    def __init__(self, base_url: str = DEFAULT_BASE_URL, cache_dir: Optional[str] = None,
                 force_refresh: bool = False, retries: int = 3, sleep: float = 0.25,
                 timeout: int = 90, verbose: bool = True):
        self.base_url = base_url.rstrip("/")
        self.cache_dir = cache_dir
        self.force_refresh = force_refresh
        self.retries = retries
        self.sleep = sleep
        self.timeout = timeout
        self.verbose = verbose
        if cache_dir:
            os.makedirs(cache_dir, exist_ok=True)

    def _cache_path(self, endpoint: str, params: Dict[str, Any]) -> Optional[str]:
        if not self.cache_dir:
            return None
        key = json.dumps({"endpoint": endpoint, "params": params}, sort_keys=True, ensure_ascii=False)
        digest = hashlib.sha256(key.encode("utf-8")).hexdigest()[:24]
        safe = endpoint.strip("/").replace("/", "_") or "root"
        return os.path.join(self.cache_dir, f"{safe}_{digest}.json")

    def get(self, endpoint: str, params: Dict[str, Any]) -> Dict[str, Any]:
        clean_params = {k: v for k, v in params.items() if v is not None and v != "" and v != []}
        cache_path = self._cache_path(endpoint, clean_params)
        if cache_path and os.path.exists(cache_path) and not self.force_refresh:
            with open(cache_path, "r", encoding="utf-8") as fh:
                return json.load(fh)
        url = f"{self.base_url}/{endpoint.lstrip('/')}"
        query = urllib.parse.urlencode(clean_params, doseq=True)
        if query:
            url += "?" + query
        last_error: Optional[BaseException] = None
        for attempt in range(1, self.retries + 1):
            try:
                req = urllib.request.Request(url, headers={"User-Agent": f"FlyGuide-PBDB/{VERSION}"})
                with urllib.request.urlopen(req, timeout=self.timeout) as resp:
                    raw = resp.read().decode("utf-8", errors="replace")
                data = json.loads(raw)
                if cache_path:
                    with open(cache_path, "w", encoding="utf-8") as fh:
                        json.dump(data, fh, indent=2, ensure_ascii=False)
                time.sleep(self.sleep)
                return data
            except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, json.JSONDecodeError) as exc:
                last_error = exc
                if self.verbose:
                    print(f"PBDB request failed attempt {attempt}/{self.retries}: {url}\n  {exc}", file=sys.stderr)
                time.sleep(self.sleep * attempt * 2)
        raise PBDBError(f"PBDB request failed after {self.retries} attempts: {url}\n{last_error}")

    def occurrences(self, params: Dict[str, Any]) -> List[Dict[str, Any]]:
        data = self.get("occs/list.json", params)
        # PBDB usually returns {"records": [...]}.
        records = data.get("records") or data.get("data") or []
        if not isinstance(records, list):
            raise PBDBError(f"Unexpected PBDB response shape: keys={list(data.keys())}")
        return records


def normalize_key(s: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", (s or "").lower())


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


def resolve_period(period: str, age_young_ma: Optional[float], age_old_ma: Optional[float]) -> Tuple[str, Optional[float], Optional[float], str, Optional[str]]:
    if age_young_ma is not None or age_old_ma is not None:
        if age_young_ma is None or age_old_ma is None:
            raise SystemExit("Use both --age-young-ma and --age-old-ma for custom ages.")
        if age_young_ma < 0 or age_old_ma < 0 or age_young_ma > age_old_ma:
            raise SystemExit("Custom ages must satisfy 0 <= young <= old, in Ma.")
        return "custom", age_young_ma, age_old_ma, f"custom {age_young_ma}-{age_old_ma} Ma", None
    key = (period or "quaternary").strip().lower().replace(" ", "_").replace("-", "_")
    if key not in PERIODS:
        valid = ", ".join(sorted(PERIODS))
        raise SystemExit(f"Unknown period '{period}'. Valid presets: {valid}; or use --age-young-ma/--age-old-ma.")
    young, old, desc, interval = PERIODS[key]
    return key, young, old, desc, interval


def organism_base_names(organisms: str, custom_base_names: Sequence[str]) -> List[str]:
    bases: List[str] = []
    if custom_base_names:
        for item in custom_base_names:
            for part in item.split(","):
                part = part.strip()
                if part:
                    bases.append(part)
        return bases
    key = (organisms or "animals").strip().lower().replace("-", "_").replace(" ", "_")
    if key not in ORGANISM_BASE_NAMES:
        valid = ", ".join(sorted(ORGANISM_BASE_NAMES))
        raise SystemExit(f"Unknown organism preset '{organisms}'. Valid: {valid}; or use --base-name.")
    return ORGANISM_BASE_NAMES[key]


def clean_taxon_name(raw: str, mode: str = "binomial", keep_genus_only: bool = True,
                     keep_higher_taxa: bool = False, drop_morphotypes: bool = False,
                     split_slash_names: bool = False) -> Optional[CleanedName]:
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
    if split_slash_names and "/" in name:
        # Keep the first plausible side of simple Canis/Vulpes style names.
        name = name.split("/")[0].strip()
    name = QUALIFIER_RE.sub(" ", name)
    had_sp_marker = bool(SP_RE.search(name))
    name = SP_RE.sub(" ", name)
    name = re.sub(r"\([^)]*\)", " ", name)
    name = re.sub(r"\b[a-z]\.?\s+(?=[A-Z][a-z])", " ", name)  # stray abbreviated authors
    name = re.sub(r"\s+", " ", name).strip()
    if not name:
        return None
    tokens = [t for t in re.split(r"\s+", name) if t]
    if not tokens:
        return None
    genus = tokens[0]
    if not re.match(r"^[A-Z][A-Za-z-]+$", genus):
        # PBDB often returns capitalized higher taxa; reject numeric/noisy names.
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
    # Remove author/year-ish debris after first three words.
    genus = tokens[0]
    epithet = tokens[1]
    if not re.match(r"^[a-z][a-z-]+$", epithet):
        # This may be a higher taxon phrase, not a species.
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


def get_first(row: Dict[str, Any], names: Sequence[str]) -> str:
    for name in names:
        val = row.get(name)
        if val not in (None, ""):
            return str(val).strip()
    # fallback case-insensitive/normalized
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


def infer_kingdom_phylum(row: Dict[str, Any], organism_hint: str) -> Tuple[str, str]:
    kingdom = get_first(row, ["kingdom"])
    phylum = get_first(row, ["phylum"])
    klass = get_first(row, ["class"])
    if kingdom:
        return kingdom, phylum
    if phylum in ANIMAL_PHYLA or klass in {"Mammalia", "Aves", "Reptilia", "Amphibia", "Actinopterygii", "Chondrichthyes"}:
        return "Animalia", phylum or "Chordata"
    if phylum in PLANT_PHYLA_HINTS:
        return "Plantae", phylum
    if phylum in FUNGI_PHYLA_HINTS:
        return "Fungi", phylum
    hint = (organism_hint or "").lower()
    if hint in {"animals", "metazoa", "vertebrates", "mammals", "birds", "fish", "reptiles_amphibians", "reptiles", "amphibians", "invertebrates", "arthropods", "molluscs"}:
        if hint in {"mammals", "birds", "fish", "reptiles", "amphibians", "vertebrates", "reptiles_amphibians"}:
            return "Animalia", phylum or "Chordata"
        return "Animalia", phylum
    if hint == "plants":
        return "Plantae", phylum
    if hint == "fungi":
        return "Fungi", phylum
    return kingdom or "Unknown", phylum


def status_passes(row: Dict[str, Any], requested: str) -> bool:
    if requested == "all":
        return True
    raw = get_first(row, ["extant", "is_extant", "is extant", "taxon_extant"])
    if not raw:
        return True  # PBDB occs usually do not include this; do not silently discard.
    val = raw.lower()
    is_extant = val in {"1", "true", "yes", "y", "extant"}
    return (requested == "extant" and is_extant) or (requested == "extinct" and not is_extant)


def local_organism_filter(row: Dict[str, Any], organisms: str) -> bool:
    key = (organisms or "").lower().replace("-", "_").replace(" ", "_")
    if key != "invertebrates":
        return True
    phylum = get_first(row, ["phylum"])
    klass = get_first(row, ["class"])
    return phylum != "Chordata" and klass not in {"Mammalia", "Aves", "Reptilia", "Amphibia", "Actinopterygii", "Chondrichthyes"}


def build_pbdb_params(args: argparse.Namespace, base_name: Optional[str], region_info: Tuple[float, float, float, float, str],
                      young_ma: Optional[float], old_ma: Optional[float], interval: Optional[str]) -> Dict[str, Any]:
    min_lon, min_lat, max_lon, max_lat, _desc = region_info
    show_fields = args.show or "coords,classext,ident,coll,loc"
    params: Dict[str, Any] = {
        "vocab": "pbdb",
        "show": show_fields,
        "limit": args.limit,
    }
    if base_name:
        params["base_name"] = base_name
    if args.taxon_name:
        params["taxon_name"] = ",".join(args.taxon_name)
    if args.country_or_continent:
        params["cc"] = args.country_or_continent
    else:
        params.update({"lngmin": min_lon, "lngmax": max_lon, "latmin": min_lat, "latmax": max_lat})
    if args.prefer_interval and interval:
        params["interval"] = interval
    else:
        if young_ma is not None:
            params["min_ma"] = young_ma
        if old_ma is not None:
            params["max_ma"] = old_ma
    return params


def aggregate_records(records: Iterable[Dict[str, Any]], args: argparse.Namespace, region_key: str,
                      period_key: str, region_desc: str, period_desc: str,
                      progress: "Optional[Any]" = None) -> Tuple[List[Dict[str, str]], List[RejectedName]]:
    grouped: Dict[str, Dict[str, Any]] = {}
    rejected: List[RejectedName] = []
    _rec_list = list(records)
    _total = len(_rec_list)
    _dash = progress or (PalaeoProgress(_total, label="records", verbose=args.verbose) if PalaeoProgress else None)
    if _dash and _total:
        _dash.set_current(f"Processing {_total:,} occurrence records", "starting")
    _REPORT_EVERY = max(1, _total // 200)
    for _i, row in enumerate(_rec_list):
        if not status_passes(row, args.status):
            continue
        if not local_organism_filter(row, args.organisms):
            continue
        source_id = get_first(row, ["occurrence_no", "occurrence_id", "record_no"])
        accepted_name = get_first(row, ["accepted_name", "accepted_name_orig", "taxon_name"])
        identified_name = get_first(row, ["identified_name", "idn", "identified_name_orig"])
        raw_name = accepted_name or identified_name
        if not raw_name:
            genus_field = get_first(row, ["genus"])
            species_field = get_first(row, ["species", "specific_epithet"])
            raw_name = f"{genus_field} {species_field}".strip()
        cleaned = clean_taxon_name(
            raw_name,
            mode=args.ncbi_name_mode,
            keep_genus_only=not args.drop_genus_only,
            keep_higher_taxa=args.keep_higher_taxa,
            drop_morphotypes=args.drop_morphotypes,
            split_slash_names=args.split_slash_names,
        )
        if not cleaned:
            rejected.append(RejectedName(source_id=source_id, source_name=raw_name, reason="could_not_clean_name", context=json.dumps(row, ensure_ascii=False)[:500]))
            continue
        key = cleaned.search_name
        kingdom, phylum = infer_kingdom_phylum(row, args.organisms)
        max_ma = to_float(get_first(row, ["max_ma", "early_age", "age_max", "older_age"]))
        min_ma = to_float(get_first(row, ["min_ma", "late_age", "age_min", "younger_age"]))
        coll_no = get_first(row, ["collection_no", "collection_id"])
        country = get_first(row, ["cc", "country", "country_name"])
        lng = get_first(row, ["lng", "longitude"])
        lat = get_first(row, ["lat", "latitude"])
        if key not in grouped:
            grouped[key] = {
                "species": cleaned.search_name,
                "genus": cleaned.genus,
                "kingdom": kingdom,
                "phylum": phylum,
                "scientificName": cleaned.search_name,
                "taxonRank": cleaned.rank,
                "taxonKey": get_first(row, ["accepted_no", "identified_no"]),
                "source": "PBDB",
                "source_database": "Paleobiology Database",
                "pbdb_accepted_name": accepted_name,
                "pbdb_identified_names": set(),
                "pbdb_accepted_no": get_first(row, ["accepted_no"]),
                "pbdb_identified_nos": set(),
                "pbdb_occurrence_count": 0,
                "pbdb_collection_nos": set(),
                "pbdb_collection_count": 0,
                "pbdb_min_ma_youngest": min_ma,
                "pbdb_max_ma_oldest": max_ma,
                "pbdb_early_intervals": set(),
                "pbdb_late_intervals": set(),
                "pbdb_countries": set(),
                "pbdb_lng_example": lng,
                "pbdb_lat_example": lat,
                "region": region_key,
                "region_description": region_desc,
                "period": period_key,
                "period_description": period_desc,
                "ncbi_search_name": cleaned.search_name,
                "canonical_binomial": cleaned.canonical_binomial,
                "name_cleaning_status": cleaned.cleaning_status,
            }
        g = grouped[key]
        g["pbdb_occurrence_count"] += 1
        if identified_name:
            g["pbdb_identified_names"].add(identified_name)
        ident_no = get_first(row, ["identified_no"])
        if ident_no:
            g["pbdb_identified_nos"].add(ident_no)
        if coll_no:
            g["pbdb_collection_nos"].add(coll_no)
        if max_ma is not None:
            old = g.get("pbdb_max_ma_oldest")
            g["pbdb_max_ma_oldest"] = max_ma if old is None else max(old, max_ma)
        if min_ma is not None:
            young = g.get("pbdb_min_ma_youngest")
            g["pbdb_min_ma_youngest"] = min_ma if young is None else min(young, min_ma)
        early = get_first(row, ["early_interval"])
        late = get_first(row, ["late_interval"])
        if early:
            g["pbdb_early_intervals"].add(early)
        if late:
            g["pbdb_late_intervals"].add(late)
        if country:
            g["pbdb_countries"].add(country)
        if not g.get("phylum") and phylum:
            g["phylum"] = phylum
        if g.get("kingdom") in {"", "Unknown"} and kingdom:
            g["kingdom"] = kingdom
        if _dash and (_i % _REPORT_EVERY == 0 or _i == _total - 1):
            unique = len(grouped)
            _dash.set_current(
                f"{key}",
                f"{unique:,} unique taxa",
                completed=_i + 1,
            )
    if _dash:
        _dash.finish()
    rows: List[Dict[str, str]] = []
    for key in sorted(grouped):
        g = grouped[key]
        g["pbdb_collection_count"] = len(g["pbdb_collection_nos"])
        for set_col in ["pbdb_identified_names", "pbdb_identified_nos", "pbdb_collection_nos", "pbdb_early_intervals", "pbdb_late_intervals", "pbdb_countries"]:
            g[set_col] = ";".join(sorted(str(x) for x in g[set_col] if x))
        for col in ["pbdb_min_ma_youngest", "pbdb_max_ma_oldest"]:
            if g.get(col) is None:
                g[col] = ""
            else:
                g[col] = f"{g[col]:.6g}"
        for k, v in list(g.items()):
            if not isinstance(v, str):
                g[k] = str(v)
        rows.append(g)
    return rows, rejected


def write_csv(path: str, rows: List[Dict[str, str]]) -> None:
    fieldnames = [
        "species", "genus", "kingdom", "phylum", "scientificName", "taxonRank", "taxonKey",
        "source", "source_database", "pbdb_accepted_name", "pbdb_identified_names",
        "pbdb_accepted_no", "pbdb_identified_nos", "pbdb_occurrence_count", "pbdb_collection_count",
        "pbdb_collection_nos", "pbdb_min_ma_youngest", "pbdb_max_ma_oldest", "pbdb_early_intervals",
        "pbdb_late_intervals", "pbdb_countries", "pbdb_lng_example", "pbdb_lat_example",
        "region", "region_description", "period", "period_description", "ncbi_search_name",
        "canonical_binomial", "name_cleaning_status",
    ]
    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
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
    if out_prefix:
        prefix = out_prefix
    else:
        base = os.path.basename(out_csv)
        prefix = os.path.splitext(base)[0]
    species_path = f"{prefix}_species_search.txt"
    kingdom_path = f"{prefix}_species_kingdom.tsv"
    seen = {}
    for row in rows:
        name = row.get("ncbi_search_name") or row.get("species")
        if not name:
            continue
        if name not in seen:
            seen[name] = (row.get("kingdom", ""), row.get("phylum", ""))
        else:
            old_k, old_p = seen[name]
            seen[name] = (old_k or row.get("kingdom", ""), old_p or row.get("phylum", ""))
    with open(species_path, "w", encoding="utf-8") as fh:
        for name in sorted(seen):
            fh.write(name + "\n")
    with open(kingdom_path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        for name in sorted(seen):
            kingdom, phylum = seen[name]
            writer.writerow([name, kingdom, phylum])
    return species_path, kingdom_path


def load_fixture(path: str) -> List[Dict[str, Any]]:
    with open(path, "r", encoding="utf-8") as fh:
        data = json.load(fh)
    if isinstance(data, list):
        return data
    return data.get("records") or data.get("data") or []


def run(args: argparse.Namespace) -> Dict[str, Any]:
    region_key, region_info = resolve_region(args.region, args.bbox)
    min_lon, min_lat, max_lon, max_lat, region_desc = region_info
    period_key, young_ma, old_ma, period_desc, interval = resolve_period(args.period, args.age_young_ma, args.age_old_ma)
    bases = organism_base_names(args.organisms, args.base_name)

    if print_header and args.verbose:
        period_label = f"{period_desc} ({young_ma}–{old_ma} Ma)" if young_ma is not None else period_desc
        print_header(
            f"FlyGuide PBDB exporter  v{VERSION}",
            [
                ("Region", region_desc),
                ("Period", period_label),
                ("Group", f"{args.organisms}" + (f" ({', '.join(bases)})" if bases else "")),
                ("Output", args.out),
            ],
        )

    all_records: List[Dict[str, Any]] = []
    query_summaries: List[Dict[str, Any]] = []
    if args.fixture:
        all_records = load_fixture(args.fixture)
        query_summaries.append({"fixture": args.fixture, "records": len(all_records)})
    else:
        client = PBDBClient(args.base_url, args.cache_dir, args.force_refresh, args.retries, args.sleep, args.timeout, args.verbose)
        query_bases: List[Optional[str]] = bases if bases else [None]
        for base in query_bases:
            params = build_pbdb_params(args, base, region_info, young_ma, old_ma, interval)
            if args.verbose:
                print(f"  Querying PBDB API (base_name={base or 'all'}, limit=all — may take a moment)...", file=sys.stderr)
            t0 = time.time()
            records = client.occurrences(params)
            elapsed = time.time() - t0
            if args.verbose:
                print(f"  ✓ Fetched {len(records):,} occurrence records in {elapsed:.1f}s", file=sys.stderr)
            all_records.extend(records)
            query_summaries.append({"base_name": base or "", "records": len(records), "params": params})

    rows, rejected = aggregate_records(all_records, args, region_key, period_key, region_desc, period_desc)
    if args.min_occurrences > 1:
        rows = [r for r in rows if int(r.get("pbdb_occurrence_count", "0") or 0) >= args.min_occurrences]
    write_csv(args.out, rows)
    rejected_path = os.path.splitext(args.out)[0] + ".rejected.tsv"
    write_rejected(rejected_path, rejected)
    outputs = {"csv": args.out, "rejected": rejected_path}
    if args.write_flyguide_files:
        sp, kd = write_flyguide_files(args.out, rows, args.out_prefix)
        outputs["species_search"] = sp
        outputs["species_kingdom"] = kd
    summary = {
        "version": VERSION,
        "source": "PBDB",
        "region": region_key,
        "region_bbox": [min_lon, min_lat, max_lon, max_lat],
        "period": period_key,
        "period_bounds_ma": [young_ma, old_ma],
        "organisms": args.organisms,
        "base_names": bases,
        "records_read": len(all_records),
        "taxa_written": len(rows),
        "names_rejected": len(rejected),
        "queries": query_summaries,
        "outputs": outputs,
    }
    summary_path = os.path.splitext(args.out)[0] + ".summary.json"
    with open(summary_path, "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)
    outputs["summary"] = summary_path
    if print_done and args.verbose:
        counts = [
            ("Occurrence records read:", f"{len(all_records):,}"),
            ("Unique search taxa written:", f"{len(rows):,}  →  {args.out}"),
            ("Rejected/noisy names:", f"{len(rejected):,}  →  {rejected_path}"),
        ]
        if args.write_flyguide_files:
            counts += [
                ("FlyGuide species search:", outputs['species_search']),
                ("FlyGuide kingdom map:", outputs['species_kingdom']),
            ]
        print_done("PBDB export complete", counts)
    elif args.verbose:
        print(f"PBDB export complete.", file=sys.stderr)
        print(f"  Occurrence records read: {len(all_records)}", file=sys.stderr)
        print(f"  Unique search taxa written: {len(rows)} -> {args.out}", file=sys.stderr)
        print(f"  Rejected/noisy names: {len(rejected)} -> {rejected_path}", file=sys.stderr)
        if args.write_flyguide_files:
            print(f"  FlyGuide species search: {outputs['species_search']}", file=sys.stderr)
            print(f"  FlyGuide kingdom map: {outputs['species_kingdom']}", file=sys.stderr)
    return summary


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Export PBDB fossil occurrence taxa to a FlyGuide/GBIF-like CSV.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--region", default="northern_hemisphere", help="Built-in broad region preset")
    p.add_argument("--bbox", help="Custom minLon,minLat,maxLon,maxLat; overrides --region")
    p.add_argument("--period", default="quaternary", help="Broad time period preset")
    p.add_argument("--age-young-ma", type=float, help="Custom youngest age in Ma")
    p.add_argument("--age-old-ma", type=float, help="Custom oldest age in Ma")
    p.add_argument("--prefer-interval", action="store_true", help="Use PBDB named interval when available instead of min_ma/max_ma")
    p.add_argument("--organisms", default="animals", help="Organism preset")
    p.add_argument("--base-name", action="append", default=[], help="PBDB base_name/clade; repeat or comma-separate. Overrides --organisms")
    p.add_argument("--taxon-name", action="append", default=[], help="PBDB taxon_name filter; repeatable")
    p.add_argument("--country-or-continent", help="PBDB cc filter, e.g. NAM,EUR,US,CA; overrides bbox")
    p.add_argument("--status", choices=["all", "extinct", "extant"], default="all", help="Best-effort extant/extinct filter if PBDB returns an extant flag")
    p.add_argument("--ncbi-name-mode", choices=["binomial", "trinomial", "as-is", "genus"], default="binomial", help="How to collapse names for NCBI searching")
    p.add_argument("--drop-genus-only", action="store_true", help="Drop genus-only names such as Canis sp.")
    p.add_argument("--keep-higher-taxa", action="store_true", help="Keep higher-taxon buckets; usually noisy for NCBI")
    p.add_argument("--drop-morphotypes", action="store_true", help="Drop -type/morphotype names")
    p.add_argument("--split-slash-names", action="store_true", help="For simple slash names, keep first side")
    p.add_argument("--min-occurrences", type=int, default=1, help="Require at least this many PBDB occurrences per cleaned taxon")
    p.add_argument("--show", default="coords,classext,ident,coll,loc", help="PBDB show parameter")
    p.add_argument("--limit", default="all", help="PBDB limit parameter; use a small number for smoke tests")
    p.add_argument("--out", help="Output GBIF-like CSV")
    p.add_argument("--out-prefix", help="Prefix for FlyGuide files; default derives from --out")
    p.add_argument("--write-flyguide-files", action="store_true", help="Also write *_species_search.txt and *_species_kingdom.tsv")
    p.add_argument("--cache-dir", default=".pbdb_cache", help="Cache PBDB JSON responses")
    p.add_argument("--force-refresh", action="store_true", help="Ignore cached PBDB responses")
    p.add_argument("--base-url", default=DEFAULT_BASE_URL, help="PBDB API base URL")
    p.add_argument("--retries", type=int, default=3, help="HTTP retries")
    p.add_argument("--sleep", type=float, default=0.25, help="Sleep between HTTP requests")
    p.add_argument("--timeout", type=int, default=90, help="HTTP timeout seconds")
    p.add_argument("--fixture", help="Offline test fixture JSON with PBDB-like records")
    p.add_argument("--quiet", dest="verbose", action="store_false", help="Reduce progress output")
    p.add_argument("--verbose", dest="verbose", action="store_true")
    p.set_defaults(verbose=True)
    p.add_argument("--list-regions", action="store_true", help="List region presets and exit")
    p.add_argument("--list-periods", action="store_true", help="List period presets and exit")
    p.add_argument("--list-organisms", action="store_true", help="List organism presets and exit")
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    p = build_parser()
    args = p.parse_args(argv)
    if args.list_regions:
        for k, v in sorted(REGION_BBOXES.items()):
            print(f"{k}\t{v[0]},{v[1]},{v[2]},{v[3]}\t{v[4]}")
        return 0
    if args.list_periods:
        for k, (young, old, desc, interval) in sorted(PERIODS.items()):
            print(f"{k}\t{young}\t{old}\t{desc}\tinterval={interval or ''}")
        return 0
    if args.list_organisms:
        for k, bases in sorted(ORGANISM_BASE_NAMES.items()):
            print(f"{k}\t{','.join(bases) if bases else '(no base_name/all PBDB)'}")
        return 0
    if not args.out:
        p.error("--out is required unless listing presets")
    run(args)
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
