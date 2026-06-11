#!/usr/bin/env python3
"""
FlyGuide GBIF live query tool
==============================

Query the GBIF occurrence API for species lists directly from the command line.
No GBIF account or API key required for read-only searches.

Outputs a FlyGuide-compatible GBIF-like CSV (species/genus/kingdom/phylum) that
feeds directly into gbif_prep_from_csv.py and NCBI-NT_Downloader.pl.

Examples
--------

# Terrestrial northern hemisphere mammals, human observations only:
python3 gbif_query.py --taxon mammals --region northern_hemisphere --basis human --out NH_mammals_gbif.csv

# All vascular plants from Canada:
python3 gbif_query.py --taxon plants --country CA --out canada_plants_gbif.csv

# Birds from multiple continents:
python3 gbif_query.py --taxon birds --continent north_america,europe,asia --out holarctic_birds_gbif.csv

# Sharks and rays (Chondrichthyes) in the North Atlantic:
python3 gbif_query.py --taxon chondrichthyes --region north_atlantic --out north_atlantic_sharks_gbif.csv

# Whales and dolphins globally (includes open ocean — no country filter):
python3 gbif_query.py --taxon cetacea --region global --out global_cetacea_gbif.csv

# Custom scientific name or GBIF taxon key:
python3 gbif_query.py --taxon Bovidae --region africa --out africa_bovidae_gbif.csv
python3 gbif_query.py --taxon-key 359 --country US --out us_mammals_gbif.csv

# List all presets:
python3 gbif_query.py --list-taxon-groups
python3 gbif_query.py --list-regions
python3 gbif_query.py --list-basis
"""
from __future__ import annotations

import argparse
import concurrent.futures
import csv
import hashlib
import json
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from typing import Any, Dict, List, Optional, Sequence, Tuple

try:
    from _palaeo_tui import PalaeoProgress, print_header, print_done
except ImportError:
    PalaeoProgress = None  # type: ignore[assignment,misc]
    print_header = print_done = None  # type: ignore[assignment]

VERSION = "0.1.0-gbif"
GBIF_API = "https://api.gbif.org/v1"

# ── Taxon presets ──────────────────────────────────────────────────────────────
# Maps common names to GBIF backbone taxon keys.
# Users can also supply any scientific name or taxon key directly.
TAXON_PRESETS: Dict[str, Tuple[int, str, str, str]] = {
    # key: (gbif_taxon_key, display_name, kingdom, phylum)
    "animals":        (1,        "Animalia",         "Animalia",  ""),
    "animalia":       (1,        "Animalia",         "Animalia",  ""),
    "plants":         (6,        "Plantae",          "Plantae",   ""),
    "plantae":        (6,        "Plantae",          "Plantae",   ""),
    "fungi":          (5,        "Fungi",            "Fungi",     ""),
    "mammals":        (359,      "Mammalia",         "Animalia",  "Chordata"),
    "mammalia":       (359,      "Mammalia",         "Animalia",  "Chordata"),
    "birds":          (212,      "Aves",             "Animalia",  "Chordata"),
    "aves":           (212,      "Aves",             "Animalia",  "Chordata"),
    "reptiles":       (358,      "Reptilia",         "Animalia",  "Chordata"),
    "reptilia":       (358,      "Reptilia",         "Animalia",  "Chordata"),
    "amphibians":     (131,      "Amphibia",         "Animalia",  "Chordata"),
    "amphibia":       (131,      "Amphibia",         "Animalia",  "Chordata"),
    "fish":           (204,      "Actinopterygii",   "Animalia",  "Chordata"),
    "actinopterygii": (204,      "Actinopterygii",   "Animalia",  "Chordata"),
    "sharks":         (10,       "Chondrichthyes",   "Animalia",  "Chordata"),
    "chondrichthyes": (10,       "Chondrichthyes",   "Animalia",  "Chordata"),
    "insects":        (216,      "Insecta",          "Animalia",  "Arthropoda"),
    "insecta":        (216,      "Insecta",          "Animalia",  "Arthropoda"),
    "arachnids":      (367,      "Arachnida",        "Animalia",  "Arthropoda"),
    "arachnida":      (367,      "Arachnida",        "Animalia",  "Arthropoda"),
    "molluscs":       (52,       "Mollusca",         "Animalia",  "Mollusca"),
    "mollusca":       (52,       "Mollusca",         "Animalia",  "Mollusca"),
    "crustaceans":    (229,      "Malacostraca",     "Animalia",  "Arthropoda"),
    "whales":         (729,      "Cetacea",          "Animalia",  "Chordata"),
    "cetacea":        (729,      "Cetacea",          "Animalia",  "Chordata"),
    "vascular_plants":(7707728,  "Tracheophyta",     "Plantae",   "Tracheophyta"),
    "tracheophyta":   (7707728,  "Tracheophyta",     "Plantae",   "Tracheophyta"),
    "vertebrates":    (11418114, "Vertebrata",       "Animalia",  "Chordata"),
    "vertebrata":     (11418114, "Vertebrata",       "Animalia",  "Chordata"),
    "algae":          (3,        "Chromista",        "Chromista", ""),
    "chromista":      (3,        "Chromista",        "Chromista", ""),
    "bacteria":       (3,        "Bacteria",         "Bacteria",  ""),
}

# ── Geographic presets ─────────────────────────────────────────────────────────
# (min_lon, min_lat, max_lon, max_lat, description)
# Ocean basin presets cover approximate open-ocean areas for marine queries.
# Records in international waters carry no country code in GBIF, so no country
# filter is used — only the bounding box locates them.

REGION_BBOXES: Dict[str, Tuple[float, float, float, float, str]] = {
    # Land/continental regions
    "global":              (-180, -90,  180,  90,  "Global"),
    "northern_hemisphere": (-180,   0,  180,  90,  "Northern Hemisphere"),
    "southern_hemisphere": (-180, -90,  180,   0,  "Southern Hemisphere"),
    "holarctic":           (-180,  23.5, 180,  90, "Holarctic (>23.5°N)"),
    "north_america":       (-170,  14,  -50,  84,  "North America"),
    "central_america":     (-93,    7,  -77,  22,  "Central America"),
    "south_america":       (-82,  -56,  -34,  13,  "South America"),
    "americas":            (-170, -56,  -34,  84,  "The Americas"),
    "europe":              (-30,   34,   45,  72,  "Europe"),
    "eurasia":             (-30,   10,  180,  80,  "Eurasia"),
    "asia":                ( 25,  -10,  180,  80,  "Asia"),
    "east_asia":           ( 95,   20,  145,  55,  "East Asia"),
    "south_asia":          ( 60,    5,  100,  40,  "South Asia"),
    "southeast_asia":      ( 92,  -10,  140,  28,  "Southeast Asia"),
    "central_asia":        ( 50,   35,  100,  55,  "Central Asia"),
    "middle_east":         ( 25,   12,   65,  43,  "Middle East"),
    "africa":              (-20,  -35,   55,  38,  "Africa"),
    "east_africa":         ( 25,  -15,   51,  18,  "East Africa"),
    "sub_saharan_africa":  (-20,  -35,   55,  22,  "Sub-Saharan Africa"),
    "beringia":            (150,   55, -140,  75,  "Beringia"),
    "oceania":             (110,  -50,  180,  10,  "Oceania"),
    "australia":           (113,  -44,  154,  -10, "Australia"),
    "new_zealand":         (165,  -47,  179,  -34, "New Zealand"),
    "greenland":           (-73,   60,  -12,  84,  "Greenland"),
    "antarctica":          (-180, -90,  180,  -60, "Antarctica"),
    # Ocean basins — for marine taxa outside national boundaries
    "north_atlantic":      (-80,    0,   20,  70,  "North Atlantic Ocean"),
    "south_atlantic":      (-60,  -60,   20,   0,  "South Atlantic Ocean"),
    "atlantic":            (-80,  -60,   20,  70,  "Atlantic Ocean"),
    "north_pacific":       (120,    0,  180,  65,  "North Pacific Ocean (western)"),
    "north_pacific_east":  (-180,   0, -100,  65,  "North Pacific Ocean (eastern)"),
    "south_pacific":       (140,  -60,  180,   0,  "South Pacific Ocean (western)"),
    "south_pacific_east":  (-180, -60,  -70,   0,  "South Pacific Ocean (eastern)"),
    "pacific":             (120,  -60,  180,  65,  "Pacific Ocean (western half)"),
    "indian_ocean":        ( 20,  -60,  120,  30,  "Indian Ocean"),
    "southern_ocean":      (-180, -90,  180,  -60, "Southern Ocean"),
    "arctic_ocean":        (-180,  70,  180,  90,  "Arctic Ocean"),
    "mediterranean":       (-6,    30,   42,  47,  "Mediterranean Sea"),
    "caribbean":           (-90,    8,  -58,  28,  "Caribbean Sea"),
    "coral_triangle":      (110,  -12,  160,  20,  "Coral Triangle"),
}

REGION_ALIASES: Dict[str, str] = {
    "nh":          "northern_hemisphere",
    "sh":          "southern_hemisphere",
    "na":          "north_america",
    "sa":          "south_america",
    "au":          "australia",
    "nz":          "new_zealand",
    "n_atlantic":  "north_atlantic",
    "s_atlantic":  "south_atlantic",
    "n_pacific":   "north_pacific",
    "s_pacific":   "south_pacific",
    "ind_ocean":   "indian_ocean",
}

# ── Continent presets (GBIF continent enum values) ─────────────────────────────
CONTINENT_MAP: Dict[str, str] = {
    "north_america":  "NORTH_AMERICA",
    "south_america":  "SOUTH_AMERICA",
    "europe":         "EUROPE",
    "asia":           "ASIA",
    "africa":         "AFRICA",
    "oceania":        "OCEANIA",
    "antarctica":     "ANTARCTICA",
    "na":             "NORTH_AMERICA",
    "sa":             "SOUTH_AMERICA",
    "eu":             "EUROPE",
}

# Convenience multi-continent expansions
CONTINENT_GROUPS: Dict[str, List[str]] = {
    "northern_hemisphere": ["NORTH_AMERICA", "EUROPE", "ASIA"],
    "holarctic":           ["NORTH_AMERICA", "EUROPE", "ASIA"],
    "americas":            ["NORTH_AMERICA", "SOUTH_AMERICA"],
    "eurasia":             ["EUROPE", "ASIA"],
    "old_world":           ["EUROPE", "ASIA", "AFRICA"],
}

# ── Basis of record presets ────────────────────────────────────────────────────
BASIS_PRESETS: Dict[str, List[str]] = {
    "human":    ["HUMAN_OBSERVATION", "PRESERVED_SPECIMEN", "OBSERVATION",
                 "LITERATURE", "MATERIAL_CITATION"],
    "observed": ["HUMAN_OBSERVATION", "OBSERVATION"],
    "specimen": ["PRESERVED_SPECIMEN"],
    "machine":  ["MACHINE_OBSERVATION"],
    "fossil":   ["FOSSIL_SPECIMEN"],
    "living":   ["LIVING_SPECIMEN"],
    "any":      [],  # no filter
    "all":      [],  # alias
    "extant":   ["HUMAN_OBSERVATION", "PRESERVED_SPECIMEN", "OBSERVATION",
                 "LITERATURE", "MATERIAL_CITATION", "MACHINE_OBSERVATION",
                 "LIVING_SPECIMEN"],
}

# ── GBIF API client ────────────────────────────────────────────────────────────

class GBIFError(RuntimeError):
    pass


class GBIFClient:
    def __init__(self, cache_dir: Optional[str] = ".gbif_cache",
                 force_refresh: bool = False,
                 retries: int = 3, sleep: float = 0.5,
                 timeout: int = 60, verbose: bool = True):
        self.cache_dir = cache_dir
        self.force_refresh = force_refresh
        self.retries = retries
        self.sleep = sleep
        self.timeout = timeout
        self.verbose = verbose
        if cache_dir:
            os.makedirs(cache_dir, exist_ok=True)

    def _cache_path(self, url: str, params: Dict) -> Optional[str]:
        if not self.cache_dir:
            return None
        key = json.dumps({"url": url, "params": params}, sort_keys=True)
        digest = hashlib.sha256(key.encode()).hexdigest()[:24]
        return os.path.join(self.cache_dir, f"gbif_{digest}.json")

    def get(self, endpoint: str, params: Dict[str, Any]) -> Dict[str, Any]:
        clean = {k: v for k, v in params.items() if v is not None and v != "" and v != []}
        url = f"{GBIF_API}/{endpoint.lstrip('/')}"
        cache = self._cache_path(url, clean)
        if cache and os.path.exists(cache) and not self.force_refresh:
            with open(cache, encoding="utf-8") as fh:
                return json.load(fh)
        # Build query string — lists become repeated params for GBIF's OR behaviour
        parts = []
        for k, v in clean.items():
            if isinstance(v, list):
                for item in v:
                    parts.append((k, str(item)))
            else:
                parts.append((k, str(v)))
        query = urllib.parse.urlencode(parts)
        full_url = f"{url}?{query}" if query else url
        last_err: Optional[BaseException] = None
        for attempt in range(1, self.retries + 1):
            try:
                req = urllib.request.Request(
                    full_url, headers={"User-Agent": f"FlyGuide-GBIF/{VERSION}"})
                with urllib.request.urlopen(req, timeout=self.timeout) as resp:
                    data = json.loads(resp.read().decode("utf-8", errors="replace"))
                if cache:
                    with open(cache, "w", encoding="utf-8") as fh:
                        json.dump(data, fh, indent=2, ensure_ascii=False)
                time.sleep(self.sleep)
                return data
            except (urllib.error.URLError, urllib.error.HTTPError,
                    TimeoutError, json.JSONDecodeError) as exc:
                last_err = exc
                if self.verbose:
                    print(f"  GBIF request failed (attempt {attempt}/{self.retries}): {exc}",
                          file=sys.stderr)
                time.sleep(self.sleep * attempt * 2)
        raise GBIFError(f"GBIF request failed after {self.retries} attempts: {full_url}\n{last_err}")


# ── Taxon resolution ───────────────────────────────────────────────────────────

def resolve_taxon(name_or_key: str, client: GBIFClient) -> Tuple[int, str, str, str, str]:
    """
    Returns (taxon_key, display_name, kingdom, phylum, rank).
    Checks presets first, then tries GBIF species suggest API.
    """
    key = name_or_key.strip().lower().replace("-", "_").replace(" ", "_")

    # Check presets
    if key in TAXON_PRESETS:
        tk, display, kingdom, phylum = TAXON_PRESETS[key]
        return tk, display, kingdom, phylum, "preset"

    # Try parsing as a numeric taxon key
    try:
        tk = int(name_or_key.strip())
        data = client.get(f"species/{tk}", {})
        display = data.get("canonicalName") or data.get("scientificName", str(tk))
        kingdom = data.get("kingdom", "")
        phylum = data.get("phylum", "")
        rank = data.get("rank", "unknown").lower()
        return tk, display, kingdom, phylum, rank
    except ValueError:
        pass

    # Try GBIF species suggest (name lookup)
    data = client.get("species/suggest", {"q": name_or_key.strip(), "limit": 5})
    if not isinstance(data, list) or not data:
        raise GBIFError(f"Could not resolve taxon: {name_or_key!r}")
    best = data[0]
    tk = best.get("key") or best.get("nubKey")
    if not tk:
        raise GBIFError(f"No GBIF taxon key found for: {name_or_key!r}")
    display = best.get("canonicalName") or best.get("scientificName", name_or_key)
    kingdom = best.get("kingdom", "")
    phylum = best.get("phylum", "")
    rank = best.get("rank", "unknown").lower()
    return int(tk), display, kingdom, phylum, rank


# ── Species list fetching ──────────────────────────────────────────────────────

FACET_LIMIT = 10000  # per chunk; raise for broader coverage, lower for speed


def _extract_facet(result: Dict, field: str) -> List[Dict]:
    for f in result.get("facets", []):
        if f.get("field", "").upper() == field.upper():
            return f.get("counts", [])
    return []


def fetch_species_keys(client: GBIFClient, query_params: Dict[str, Any],
                       verbose: bool) -> Dict[str, int]:
    """
    Returns {speciesKey: occurrence_count} for all species matching the query.
    Uses facet=speciesKey (the only working GBIF species facet).
    Automatically chunks by taxonomic order when the facet limit is hit.
    """
    params = {**query_params, "limit": 0, "facet": "speciesKey",
              "facetLimit": FACET_LIMIT}
    result = client.get("occurrence/search", params)
    counts = _extract_facet(result, "SPECIES_KEY")

    if len(counts) < FACET_LIMIT:
        if verbose:
            print(f"  ✓ {len(counts):,} unique species keys retrieved", file=sys.stderr)
        return {c["name"]: c["count"] for c in counts}

    # Hit the facet limit — chunk by taxonomic order
    if verbose:
        print(f"  Species count hit facet limit ({FACET_LIMIT:,}). "
              f"Chunking by taxonomic order...", file=sys.stderr)
    return _fetch_chunked_keys(client, query_params, verbose)


def _fetch_chunked_keys(client: GBIFClient, query_params: Dict[str, Any],
                        verbose: bool) -> Dict[str, int]:
    """Chunk a large query by order (or family), collect all speciesKeys."""
    order_params = {**query_params, "limit": 0, "facet": "orderKey",
                    "facetLimit": 2000}
    result = client.get("occurrence/search", order_params)
    chunks = _extract_facet(result, "ORDER_KEY")
    chunk_param = "orderKey"
    chunk_label = "order"

    if not chunks:
        fam_params = {**query_params, "limit": 0, "facet": "familyKey",
                      "facetLimit": 5000}
        result = client.get("occurrence/search", fam_params)
        chunks = _extract_facet(result, "FAMILY_KEY")
        chunk_param = "familyKey"
        chunk_label = "family"

    all_keys: Dict[str, int] = {}
    dash = (PalaeoProgress(len(chunks), label=f"{chunk_label} chunks", verbose=verbose)
            if PalaeoProgress else None)

    for chunk in chunks:
        chunk_key = chunk["name"]
        chunk_count = chunk["count"]
        if dash:
            dash.set_current(f"{chunk_label} {chunk_key}", f"{chunk_count:,} occs")

        sp_params = {**query_params, chunk_param: chunk_key,
                     "limit": 0, "facet": "speciesKey", "facetLimit": FACET_LIMIT}
        r = client.get("occurrence/search", sp_params)
        for c in _extract_facet(r, "SPECIES_KEY"):
            key = c["name"]
            cnt = c["count"]
            if key not in all_keys or cnt > all_keys[key]:
                all_keys[key] = cnt

        if dash:
            dash.update(f"{chunk_label} {chunk_key}",
                        f"[{len(all_keys):>6,} unique so far]")
        elif verbose:
            print(f"  {chunk_label} {chunk_key}: total unique keys = {len(all_keys):,}",
                  file=sys.stderr)

    if dash:
        dash.finish()
    return all_keys


def resolve_species_keys(client: GBIFClient, keys: Dict[str, int],
                         verbose: bool,
                         workers: int = 10) -> List[Tuple[str, str, str, str, int]]:
    """
    Resolve {speciesKey: count} → [(name, genus, kingdom, phylum, count)].
    Calls /v1/species/{key} in parallel (default 10 threads).
    Results are cached on disk so subsequent runs are instant.
    """
    items = list(keys.items())
    dash = (PalaeoProgress(len(items), label="species", verbose=verbose)
            if PalaeoProgress else None)
    if dash:
        dash.set_current("resolving species names (parallel, cached)", "starting")

    resolved: List[Tuple[str, str, str, str, int]] = []

    def _resolve_one(item: Tuple[str, int]) -> Optional[Tuple[str, str, str, str, int]]:
        sp_key, occ_count = item
        try:
            data = client.get(f"species/{sp_key}", {})
        except GBIFError:
            return None
        name = data.get("canonicalName") or data.get("scientificName", "")
        genus = data.get("genus", name.split()[0] if name else "")
        kingdom = data.get("kingdom", "")
        phylum = data.get("phylum", "")
        if not name:
            return None
        return (name, genus, kingdom, phylum, occ_count)

    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(_resolve_one, item): item for item in items}
        for fut in concurrent.futures.as_completed(futures):
            result = fut.result()
            if result:
                resolved.append(result)
                if dash:
                    dash.update(result[0], f"[{result[4]:>8,} occs]")
                elif verbose and len(resolved) % 50 == 0:
                    print(f"  Resolved {len(resolved)}/{len(items)} species...",
                          file=sys.stderr)

    if dash:
        dash.finish()
    resolved.sort(key=lambda r: r[0])
    return resolved


# ── Name cleaning ──────────────────────────────────────────────────────────────

def clean_species_name(raw: str) -> Optional[str]:
    """Return a clean binomial or genus-only string, or None to reject."""
    name = re.sub(r"\s+", " ", raw.strip())
    # Remove author/year parentheticals and trailing noise
    name = re.sub(r"\([^)]*\)", " ", name)
    name = re.sub(r"\s+(sensu|auct|agg|sl|ss)\b.*", "", name, flags=re.I)
    name = re.sub(r"[^A-Za-z \-]", " ", name)
    tokens = [t for t in name.split() if len(t) > 1 and not t[0].isdigit()]
    if not tokens:
        return None
    # Keep genus + species only (binomial)
    if len(tokens) >= 2:
        return f"{tokens[0].capitalize()} {tokens[1].lower()}"
    return tokens[0].capitalize()


# ── Output helpers ─────────────────────────────────────────────────────────────

OUTPUT_FIELDS = [
    "species", "genus", "kingdom", "phylum", "scientificName", "taxonRank",
    "source", "source_database", "ncbi_search_name", "canonical_binomial",
    "gbif_occurrence_count", "gbif_taxon_name", "region", "country",
    "continent", "basis_of_record", "name_cleaning_status",
]


def build_rows(resolved: List[Tuple[str, str, str, str, int]],
               taxon_kingdom: str, taxon_phylum: str,
               region: str, countries: List[str],
               continents: List[str], basis: List[str]) -> Tuple[List[Dict], List[Tuple[str, str]]]:
    """
    resolved: [(name, genus, kingdom, phylum, occ_count)]
    taxon_kingdom/phylum: fallback when per-species values are blank.
    """
    rows: List[Dict] = []
    rejected: List[Tuple[str, str]] = []
    country_str = ";".join(countries)
    continent_str = ";".join(continents)
    basis_str = ";".join(basis) if basis else "any"

    for raw_name, genus, rec_kingdom, rec_phylum, occ_count in resolved:
        cleaned = clean_species_name(raw_name)
        if not cleaned:
            rejected.append((raw_name, "could_not_clean"))
            continue
        tokens = cleaned.split()
        clean_genus = genus or tokens[0]
        rank = "species" if len(tokens) >= 2 else "genus_or_higher"
        canonical = " ".join(tokens[:2]) if len(tokens) >= 2 else tokens[0]
        kingdom = rec_kingdom or taxon_kingdom
        phylum = rec_phylum or taxon_phylum
        rows.append({
            "species": cleaned,
            "genus": clean_genus,
            "kingdom": kingdom,
            "phylum": phylum,
            "scientificName": cleaned,
            "taxonRank": rank,
            "source": "GBIF",
            "source_database": "GBIF",
            "ncbi_search_name": cleaned,
            "canonical_binomial": canonical,
            "gbif_occurrence_count": str(occ_count),
            "gbif_taxon_name": raw_name,
            "region": region,
            "country": country_str,
            "continent": continent_str,
            "basis_of_record": basis_str,
            "name_cleaning_status": "ok" if cleaned == raw_name else "cleaned",
        })
    rows.sort(key=lambda r: r["species"])
    return rows, rejected


def write_csv(path: str, rows: List[Dict]) -> None:
    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=OUTPUT_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_rejected(path: str, rejected: List[Tuple[str, str]]) -> None:
    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["source_name", "reason"])
        writer.writerows(rejected)


def write_flyguide_files(out_csv: str, rows: List[Dict],
                         prefix: Optional[str]) -> Tuple[str, str]:
    base = prefix or os.path.splitext(os.path.basename(out_csv))[0]
    sp_path = f"{base}_species_search.txt"
    kd_path = f"{base}_species_kingdom.tsv"
    with open(sp_path, "w", encoding="utf-8") as fh:
        for r in rows:
            fh.write(r["ncbi_search_name"] + "\n")
    with open(kd_path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        for r in rows:
            writer.writerow([r["ncbi_search_name"], r["kingdom"], r["phylum"]])
    return sp_path, kd_path


# ── Query parameter builder ────────────────────────────────────────────────────

def build_query_params(taxon_key: int, region_bbox: Optional[Tuple],
                       countries: List[str], continents: List[str],
                       basis: List[str], year_range: Optional[str],
                       has_coordinate: bool) -> Dict[str, Any]:
    params: Dict[str, Any] = {"taxonKey": taxon_key}

    if region_bbox:
        min_lon, min_lat, max_lon, max_lat, _ = region_bbox
        params["decimalLatitude"] = f"{min_lat},{max_lat}"
        params["decimalLongitude"] = f"{min_lon},{max_lon}"

    if countries:
        params["country"] = countries  # list → repeated params → OR

    if continents and not region_bbox:
        params["continent"] = continents  # list → repeated params → OR

    if basis:
        params["basisOfRecord"] = basis  # list → repeated params → OR

    if year_range:
        params["year"] = year_range

    if has_coordinate:
        params["hasCoordinate"] = "true"
        params["hasGeospatialIssue"] = "false"

    return params


# ── Listing helpers ────────────────────────────────────────────────────────────

def list_taxon_groups() -> None:
    print("\nBuilt-in taxon group presets:\n")
    seen: Dict[int, str] = {}
    for name, (key, display, kingdom, phylum) in sorted(TAXON_PRESETS.items()):
        if key not in seen:
            seen[key] = name
            phy = f" / {phylum}" if phylum else ""
            print(f"  {name:<22}  GBIF key {key:<10}  {kingdom}{phy}")
    print("\nAlso accepts any scientific name (e.g. Bovidae, Nymphaeales)")
    print("or a numeric GBIF taxon key (e.g. --taxon-key 359)")


def list_regions() -> None:
    print("\nBuilt-in region presets:\n")
    print("  Land/continental regions:")
    ocean_keys = {"north_atlantic", "south_atlantic", "atlantic", "north_pacific",
                  "north_pacific_east", "south_pacific", "south_pacific_east",
                  "pacific", "indian_ocean", "southern_ocean", "arctic_ocean",
                  "mediterranean", "caribbean", "coral_triangle"}
    for key, bbox in REGION_BBOXES.items():
        tag = "  [ocean]" if key in ocean_keys else ""
        print(f"  {key:<25}  {bbox[4]}{tag}")
    print("\nOcean basin regions (no country filter — covers international waters):")
    for key in ocean_keys:
        if key in REGION_BBOXES:
            bbox = REGION_BBOXES[key]
            print(f"  {key:<25}  {bbox[4]}")
    print("\nAlso accepts: --bbox minLon,minLat,maxLon,maxLat")
    print("              --country CA,US,MX  (ISO 2-letter codes, comma-separated)")
    print("              --continent north_america,europe,asia  (comma-separated)")


def list_basis() -> None:
    print("\nBasis-of-record presets:\n")
    descriptions = {
        "human":    "HUMAN_OBSERVATION, PRESERVED_SPECIMEN, OBSERVATION, LITERATURE  (recommended for most uses)",
        "observed": "HUMAN_OBSERVATION, OBSERVATION only",
        "specimen": "PRESERVED_SPECIMEN only  (museum collections)",
        "machine":  "MACHINE_OBSERVATION only  (camera traps, acoustic, etc.)",
        "fossil":   "FOSSIL_SPECIMEN only",
        "living":   "LIVING_SPECIMEN only  (zoos, botanical gardens)",
        "extant":   "All non-fossil records",
        "any":      "No filter — includes fossils and all record types",
    }
    for k, desc in descriptions.items():
        print(f"  {k:<12}  {desc}")


# ── Main ───────────────────────────────────────────────────────────────────────

def main(argv: Optional[Sequence[str]] = None) -> int:
    # Handle listing flags before the main parser so they don't require --taxon/--out.
    _pre = argparse.ArgumentParser(add_help=False)
    _pre.add_argument("--list-taxon-groups", action="store_true")
    _pre.add_argument("--list-regions", action="store_true")
    _pre.add_argument("--list-basis", action="store_true")
    _pre_args, _ = _pre.parse_known_args(argv)
    if _pre_args.list_taxon_groups:
        list_taxon_groups()
        return 0
    if _pre_args.list_regions:
        list_regions()
        return 0
    if _pre_args.list_basis:
        list_basis()
        return 0

    p = argparse.ArgumentParser(
        description=(
            "Query GBIF occurrence data for species lists from the command line. "
            "No GBIF account or API key required. Outputs a FlyGuide-compatible "
            "GBIF-like CSV with species/genus/kingdom/phylum columns."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples\n"
            "--------\n"
            "  # Northern hemisphere mammals, human observations:\n"
            "  python3 gbif_query.py --taxon mammals --region northern_hemisphere \\\n"
            "      --basis human --out NH_mammals_gbif.csv --write-flyguide-files\n\n"
            "  # Vascular plants from Canada:\n"
            "  python3 gbif_query.py --taxon vascular_plants --country CA \\\n"
            "      --out canada_plants_gbif.csv\n\n"
            "  # Birds across the holarctic (multi-continent query):\n"
            "  python3 gbif_query.py --taxon birds \\\n"
            "      --continent north_america,europe,asia --basis human \\\n"
            "      --out holarctic_birds_gbif.csv\n\n"
            "  # Sharks and rays in the North Atlantic (open ocean included):\n"
            "  python3 gbif_query.py --taxon chondrichthyes \\\n"
            "      --region north_atlantic --out NAtlantic_sharks_gbif.csv\n\n"
            "  # Cetacea globally (covers high seas — no country filter):\n"
            "  python3 gbif_query.py --taxon cetacea --region global \\\n"
            "      --out global_cetacea_gbif.csv --write-flyguide-files\n\n"
            "  # Custom taxon by scientific name:\n"
            "  python3 gbif_query.py --taxon Bovidae --region africa \\\n"
            "      --out africa_bovidae_gbif.csv\n\n"
            "  # With coordinate quality filter and year range:\n"
            "  python3 gbif_query.py --taxon mammals --country US \\\n"
            "      --require-coords --year 2000,2024 --out US_mammals_gbif.csv\n"
        ),
    )

    # Taxon
    taxon_grp = p.add_mutually_exclusive_group(required=True)
    taxon_grp.add_argument(
        "--taxon", metavar="GROUP",
        help="Taxon group name (mammals, birds, plants, insects, …) or scientific name "
             "(Bovidae, Nymphaeales, …). Run --list-taxon-groups for built-in presets.",
    )
    taxon_grp.add_argument(
        "--taxon-key", type=int, metavar="KEY",
        help="GBIF numeric taxon key (e.g. 359 for Mammalia). Overrides --taxon.",
    )

    # Geography — mutually exclusive strategies
    geo_grp = p.add_argument_group("geographic filters")
    geo_grp.add_argument(
        "--region", metavar="PRESET",
        help="Named geographic preset (northern_hemisphere, north_atlantic, …). "
             "Run --list-regions for full list.",
    )
    geo_grp.add_argument(
        "--bbox", metavar="minLon,minLat,maxLon,maxLat",
        help="Custom bounding box. Can be combined with --country/--continent.",
    )
    geo_grp.add_argument(
        "--country", metavar="CC[,CC,…]",
        help="ISO 2-letter country codes, comma-separated (CA, US, MX, …). "
             "Multiple codes are OR'd.",
    )
    geo_grp.add_argument(
        "--continent", metavar="NAME[,NAME,…]",
        help="Continent names, comma-separated (north_america, europe, asia, …). "
             "Multiple are OR'd. Ignored when --region bbox is used.",
    )

    # Record filters
    filt_grp = p.add_argument_group("record filters")
    filt_grp.add_argument(
        "--basis", metavar="PRESET", default="human",
        help="Basis-of-record filter preset: human (default), observed, specimen, "
             "machine, fossil, extant, any. Run --list-basis for details.",
    )
    filt_grp.add_argument(
        "--year", metavar="START,END",
        help="Year range filter, e.g. 2000,2024.",
    )
    filt_grp.add_argument(
        "--require-coords", action="store_true",
        help="Only include records with coordinates and no geospatial issues.",
    )
    filt_grp.add_argument(
        "--min-occurrences", type=int, default=1, metavar="N",
        help="Exclude species with fewer than N GBIF occurrences matching the query.",
    )

    # Output
    out_grp = p.add_argument_group("output")
    out_grp.add_argument("--out", required=False, metavar="FILE",
                         help="Output GBIF-like CSV.")
    out_grp.add_argument("--out-prefix", metavar="PREFIX",
                         help="Prefix for FlyGuide species files (default: base name of --out).")
    out_grp.add_argument("--write-flyguide-files", action="store_true",
                         help="Also write PREFIX_species_search.txt and PREFIX_species_kingdom.tsv.")

    # API / cache
    api_grp = p.add_argument_group("API and caching")
    api_grp.add_argument("--cache-dir", default=".gbif_cache", metavar="DIR",
                          help="Directory for caching API responses. Default: .gbif_cache")
    api_grp.add_argument("--no-cache", action="store_true",
                          help="Disable caching (always fetch from GBIF).")
    api_grp.add_argument("--force-refresh", action="store_true",
                          help="Ignore cached responses and re-fetch from GBIF.")
    api_grp.add_argument("--sleep", type=float, default=0.3, metavar="SEC",
                          help="Seconds to sleep between API calls. Default: 0.3")

    # Listing
    p.add_argument("--list-taxon-groups", action="store_true",
                   help="List built-in taxon group presets and exit.")
    p.add_argument("--list-regions", action="store_true",
                   help="List built-in geographic region presets and exit.")
    p.add_argument("--list-basis", action="store_true",
                   help="List basis-of-record presets and exit.")
    p.add_argument("--verbose", action="store_true", default=True,
                   help="Show progress (default on).")
    p.add_argument("--quiet", dest="verbose", action="store_false",
                   help="Suppress progress output.")
    p.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")

    args = p.parse_args(argv)

    if not args.out:
        p.error("--out is required")

    # ── Resolve taxon ──────────────────────────────────────────────────────────
    client = GBIFClient(
        cache_dir=None if args.no_cache else args.cache_dir,
        force_refresh=args.force_refresh,
        sleep=args.sleep,
        verbose=args.verbose,
    )

    if args.taxon_key:
        if args.verbose:
            print(f"  Resolving GBIF taxon key {args.taxon_key}...", file=sys.stderr)
        taxon_key, taxon_display, kingdom, phylum, _ = resolve_taxon(
            str(args.taxon_key), client)
    else:
        if args.verbose:
            print(f"  Resolving taxon: {args.taxon!r}...", file=sys.stderr)
        taxon_key, taxon_display, kingdom, phylum, _ = resolve_taxon(args.taxon, client)

    # ── Resolve geography ──────────────────────────────────────────────────────
    region_bbox: Optional[Tuple] = None
    region_label = "global"

    if args.bbox:
        try:
            parts = [float(x) for x in args.bbox.split(",")]
            if len(parts) != 4:
                raise ValueError
            region_bbox = (*parts, "custom bbox")
            region_label = "custom_bbox"
        except ValueError:
            p.error("--bbox must be minLon,minLat,maxLon,maxLat (four numbers)")

    elif args.region:
        key = args.region.strip().lower().replace("-", "_")
        key = REGION_ALIASES.get(key, key)
        if key not in REGION_BBOXES:
            p.error(f"Unknown region: {args.region!r}. Run --list-regions for options.")
        region_bbox = REGION_BBOXES[key]
        region_label = key

    countries: List[str] = []
    if args.country:
        countries = [c.strip().upper() for c in args.country.split(",") if c.strip()]

    continents: List[str] = []
    if args.continent:
        for c in args.continent.split(","):
            c = c.strip().lower().replace("-", "_").replace(" ", "_")
            if c in CONTINENT_GROUPS:
                continents.extend(CONTINENT_GROUPS[c])
            elif c in CONTINENT_MAP:
                continents.append(CONTINENT_MAP[c])
            else:
                p.error(f"Unknown continent: {c!r}. "
                        f"Valid: {', '.join(CONTINENT_MAP)}")

    # ── Resolve basis of record ────────────────────────────────────────────────
    basis_key = (args.basis or "human").lower()
    if basis_key not in BASIS_PRESETS:
        p.error(f"Unknown --basis: {args.basis!r}. Run --list-basis for options.")
    basis = BASIS_PRESETS[basis_key]

    # ── Print header ───────────────────────────────────────────────────────────
    if print_header and args.verbose:
        geo_parts = []
        if region_label != "global":
            geo_parts.append(region_label)
        if countries:
            geo_parts.append(f"countries: {', '.join(countries)}")
        if continents:
            geo_parts.append(f"continents: {', '.join(continents)}")
        geo_str = ", ".join(geo_parts) if geo_parts else "global (no filter)"

        print_header(
            f"FlyGuide GBIF query  v{VERSION}",
            [
                ("Taxon", f"{taxon_display}  (key: {taxon_key})"),
                ("Kingdom", f"{kingdom}" + (f" / {phylum}" if phylum else "")),
                ("Geography", geo_str),
                ("Basis", basis_key + (f"  ({', '.join(basis)})" if basis else "  (any)")),
                ("Output", args.out),
            ],
        )

    # ── Build query and fetch ──────────────────────────────────────────────────
    query_params = build_query_params(
        taxon_key=taxon_key,
        region_bbox=region_bbox,
        countries=countries,
        continents=continents,
        basis=basis,
        year_range=args.year,
        has_coordinate=args.require_coords,
    )

    if args.verbose:
        print(f"  Querying GBIF occurrence API...", file=sys.stderr)

    # Check total occurrence count
    check = client.get("occurrence/search", {**query_params, "limit": 0})
    total_occs = check.get("count", 0)
    if args.verbose:
        print(f"  Total matching occurrences: {total_occs:,}", file=sys.stderr)

    # Step 1: get species keys via facet query (fast, chunked if needed)
    species_keys = fetch_species_keys(client, query_params, args.verbose)

    # Apply min-occurrence filter before name resolution to skip pointless API calls
    if args.min_occurrences > 1:
        before = len(species_keys)
        species_keys = {k: v for k, v in species_keys.items()
                        if v >= args.min_occurrences}
        if args.verbose:
            print(f"  Dropped {before - len(species_keys)} species keys below "
                  f"min-occurrences={args.min_occurrences}", file=sys.stderr)

    if args.verbose:
        print(f"  Resolving {len(species_keys):,} species keys to names "
              f"(cached after first run)...", file=sys.stderr)

    # Step 2: resolve keys → names + taxonomy (cached on disk)
    resolved = resolve_species_keys(client, species_keys, args.verbose)

    # ── Build output rows ──────────────────────────────────────────────────────
    rows, rejected = build_rows(
        resolved, kingdom, phylum,
        region_label, countries, continents, basis,
    )

    # ── Write outputs ──────────────────────────────────────────────────────────
    write_csv(args.out, rows)
    rejected_path = os.path.splitext(args.out)[0] + ".rejected.tsv"
    write_rejected(rejected_path, rejected)

    sp = kd = None
    if args.write_flyguide_files:
        sp, kd = write_flyguide_files(args.out, rows, args.out_prefix)

    # ── Summary ────────────────────────────────────────────────────────────────
    if print_done and args.verbose:
        counts = [
            ("GBIF occurrences matched:", f"{total_occs:,}"),
            ("Unique species written:", f"{len(rows):,}  →  {args.out}"),
            ("Names rejected/noisy:", f"{len(rejected):,}  →  {rejected_path}"),
        ]
        if sp:
            counts += [
                ("FlyGuide species search:", sp),
                ("FlyGuide kingdom map:", kd),
            ]
        print_done("GBIF query complete", counts)
    elif args.verbose:
        print(f"\nGBIF query complete.", file=sys.stderr)
        print(f"  GBIF occurrences matched : {total_occs:,}", file=sys.stderr)
        print(f"  Unique species written   : {len(rows):,} -> {args.out}", file=sys.stderr)
        print(f"  Names rejected/noisy     : {len(rejected):,} -> {rejected_path}", file=sys.stderr)
        if sp:
            print(f"  FlyGuide species search  : {sp}", file=sys.stderr)
            print(f"  FlyGuide kingdom map     : {kd}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
