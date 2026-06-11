#!/usr/bin/env python3
"""
FlyGuide Neotoma exporter
=========================

Create a GBIF-like taxon CSV from Neotoma occurrence data so FlyGuide's
existing GBIF prep + NCBI downloader can be reused for extinct/palaeo taxa.

This script intentionally uses only the Python standard library.

Typical examples
----------------

# Extinct animal taxa from the Quaternary in the Northern Hemisphere:
python neotoma_extinct_to_gbif.py \
  --region northern_hemisphere \
  --period quaternary \
  --organisms animals \
  --status extinct \
  --out northern_quaternary_animals_neotoma_gbif.csv \
  --write-flyguide-files

# Extinct plants from North America in the Holocene:
python neotoma_extinct_to_gbif.py \
  --region north_america \
  --period holocene \
  --organisms plants \
  --status extinct \
  --out na_holocene_plants_neotoma_gbif.csv

Notes
-----
Neotoma is an occurrence/palaeoecological database, not a global taxonomic
checklist. A result means: "a Neotoma occurrence in the requested space/time
whose taxon is marked extinct/selected and can be cleaned into a search bucket".
"""

from __future__ import annotations

import argparse
import csv
import dataclasses
import hashlib
import json
import math
import os
import re
import sys
import time
import textwrap
import urllib.error
import urllib.parse
import urllib.request
from collections import defaultdict
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

try:
    from _palaeo_tui import PalaeoProgress, print_header, print_done
except ImportError:
    PalaeoProgress = None  # type: ignore[assignment,misc]
    print_header = print_done = None  # type: ignore[assignment]

VERSION = "0.1.0-neotoma"
DEFAULT_BASE_URL = "https://api.neotomadb.org/v2.0"

# Geological period bounds in calendar years BP.  Neotoma occurrence queries use
# ageyoung (youngest) and ageold (oldest).  These are deliberately broad buckets.
PERIODS: Dict[str, Tuple[Optional[int], Optional[int], str]] = {
    "all": (None, None, "No age filter"),
    "quaternary": (0, 2_580_000, "Quaternary, ~0-2.58 Ma BP"),
    "pleistocene": (11_700, 2_580_000, "Pleistocene, ~11.7 ka-2.58 Ma BP"),
    "holocene": (0, 11_700, "Holocene, ~0-11.7 ka BP"),
    "late_pleistocene": (11_700, 129_000, "Late Pleistocene, ~11.7-129 ka BP"),
    "middle_pleistocene": (129_000, 774_000, "Middle Pleistocene, ~129-774 ka BP"),
    "early_pleistocene": (774_000, 2_580_000, "Early Pleistocene, ~774 ka-2.58 Ma BP"),
    "late_quaternary": (0, 129_000, "Late Quaternary, ~0-129 ka BP"),
    "lgm": (19_000, 26_500, "Last Glacial Maximum, broad ~19-26.5 ka BP"),
    "last_glacial": (11_700, 115_000, "Last glacial interval, broad ~11.7-115 ka BP"),
    "pliocene": (2_580_000, 5_333_000, "Pliocene, ~2.58-5.333 Ma BP"),
}

# Region presets are intentionally simple.  They are WKT polygons/bounding boxes
# that work with the Neotoma occurrence endpoint.  Use --geojson or --bbox for
# more careful biogeographic boundaries.
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

# Dataset type hints.  These are used only when --api-prefilter is auto/on.
# Unknown/new Neotoma dataset types are still caught by local taxagroup/name filters
# if --api-prefilter off or if the API rejects a hint and we retry without it.
ORGANISM_DATASET_HINTS: Dict[str, List[str]] = {
    "animals": ["vertebrate fauna", "mollusks", "insects", "beetles", "ostracodes", "cladocera"],
    "plants": ["pollen", "plant macrofossils"],
    "fungi": ["fungal spores"],
    "diatoms": ["diatoms"],
}

# Loose TaxaGroupID mapping.  Neotoma taxagroups change/grow; this mapping is
# conservative but not exhaustive.  Datasettype and taxon-name fallbacks help.
TAXAGROUP_MAP: Dict[str, Tuple[str, str, str]] = {
    # vertebrates
    "MAM": ("Animalia", "Chordata", "animals"),
    "AVE": ("Animalia", "Chordata", "animals"),
    "BRD": ("Animalia", "Chordata", "animals"),
    "FSH": ("Animalia", "Chordata", "animals"),
    "FIS": ("Animalia", "Chordata", "animals"),
    "HERP": ("Animalia", "Chordata", "animals"),
    "REP": ("Animalia", "Chordata", "animals"),
    "AMP": ("Animalia", "Chordata", "animals"),
    # invertebrates / arthropods / mollusks
    "BTL": ("Animalia", "Arthropoda", "animals"),
    "INS": ("Animalia", "Arthropoda", "animals"),
    "CHI": ("Animalia", "Arthropoda", "animals"),
    "CLA": ("Animalia", "Arthropoda", "animals"),
    "OST": ("Animalia", "Arthropoda", "animals"),
    "COP": ("Animalia", "Arthropoda", "animals"),
    "MOL": ("Animalia", "Mollusca", "animals"),
    # plants
    "VPL": ("Plantae", "Tracheophyta", "plants"),
    "PLT": ("Plantae", "Tracheophyta", "plants"),
    "BRY": ("Plantae", "Bryophyta", "plants"),
    "PTE": ("Plantae", "Tracheophyta", "plants"),
    "ALG": ("Plantae", "", "plants"),
    # fungi / protists / algae-ish groups
    "FNG": ("Fungi", "", "fungi"),
    "FUN": ("Fungi", "", "fungi"),
    "DIA": ("Chromista", "Bacillariophyta", "diatoms"),
    "FOR": ("Rhizaria", "Foraminifera", "protists"),
}

HIGHER_TAXON_SUFFIXES = (
    "idae", "inae", "ini", "aceae", "ales", "iformes", "phyta", "mycota",
    "oda", "ata", "opsida", "idae/", "aceae/"
)
BAD_SINGLE_WORDS = {
    "unknown", "indeterminate", "undetermined", "unidentified", "other",
    "charcoal", "organic", "inorganic", "water", "sediment", "sample",
}
QUALIFIER_RE = re.compile(r"\b(cf|aff|nr|near|ex gr|gr|group|complex)\.?\b", re.IGNORECASE)
SP_RE = re.compile(r"\b(spp|sp|species|undiff|undifferentiated)\.?\b", re.IGNORECASE)
TYPE_RE = re.compile(r"(?:-type\b|\btype\b|\bmorphotype\b)", re.IGNORECASE)
BAD_CHARS_RE = re.compile(r"[\[\]{}<>≈~;:,=+*#@!]")


@dataclasses.dataclass
class CleanedName:
    search_name: str
    genus: str
    rank: str
    canonical_binomial: str
    cleaning_status: str


@dataclasses.dataclass
class RejectedName:
    taxonid: str
    taxonname: str
    reason: str
    context: str = ""


class NeotomaError(RuntimeError):
    pass


class NeotomaClient:
    def __init__(self, base_url: str = DEFAULT_BASE_URL, cache_dir: Optional[str] = None,
                 force_refresh: bool = False, sleep: float = 0.15, retries: int = 3,
                 timeout: int = 60, verbose: bool = True):
        self.base_url = base_url.rstrip("/")
        self.cache_dir = cache_dir
        self.force_refresh = force_refresh
        self.sleep = sleep
        self.retries = retries
        self.timeout = timeout
        self.verbose = verbose
        if cache_dir:
            os.makedirs(cache_dir, exist_ok=True)

    def _cache_path(self, endpoint: str, params: Dict[str, Any]) -> Optional[str]:
        if not self.cache_dir:
            return None
        key = json.dumps({"endpoint": endpoint, "params": params}, sort_keys=True, ensure_ascii=False)
        digest = hashlib.sha256(key.encode("utf-8")).hexdigest()[:24]
        safe_endpoint = endpoint.strip("/").replace("/", "_") or "root"
        return os.path.join(self.cache_dir, f"{safe_endpoint}_{digest}.json")

    def get(self, endpoint: str, params: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
        params = {k: v for k, v in (params or {}).items() if v is not None and v != ""}
        cache_path = self._cache_path(endpoint, params)
        if cache_path and os.path.exists(cache_path) and not self.force_refresh:
            with open(cache_path, "r", encoding="utf-8") as fh:
                return json.load(fh)

        url = self.base_url + "/" + endpoint.lstrip("/")
        if params:
            url += "?" + urllib.parse.urlencode(params, doseq=True)

        last_error: Optional[BaseException] = None
        for attempt in range(1, self.retries + 1):
            try:
                req = urllib.request.Request(url, headers={"User-Agent": f"FlyGuide-Neotoma/{VERSION}"})
                with urllib.request.urlopen(req, timeout=self.timeout) as resp:
                    text = resp.read().decode("utf-8")
                payload = json.loads(text)
                if cache_path:
                    tmp = cache_path + ".tmp"
                    with open(tmp, "w", encoding="utf-8") as fh:
                        json.dump(payload, fh, ensure_ascii=False)
                    os.replace(tmp, cache_path)
                if self.sleep:
                    time.sleep(self.sleep)
                return payload
            except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
                last_error = exc
                if attempt < self.retries:
                    time.sleep(min(8.0, 0.75 * (2 ** (attempt - 1))))
        raise NeotomaError(f"Neotoma request failed after {self.retries} attempts: {url}\n{last_error}")

    def get_data(self, endpoint: str, params: Optional[Dict[str, Any]] = None) -> List[Dict[str, Any]]:
        payload = self.get(endpoint, params)
        status = str(payload.get("status", "")).lower()
        if status not in ("", "success"):
            raise NeotomaError(f"Neotoma status={payload.get('status')} message={payload.get('message')}")
        data = payload.get("data", [])
        if data is None:
            return []
        if not isinstance(data, list):
            raise NeotomaError(f"Expected list in response data for {endpoint}, got {type(data)}")
        return data

    def paged_get_data(self, endpoint: str, base_params: Optional[Dict[str, Any]] = None,
                       page_size: int = 5000, max_records: Optional[int] = None,
                       progress_label: str = "records") -> List[Dict[str, Any]]:
        all_rows: List[Dict[str, Any]] = []
        offset = 0
        page_num = 0
        estimated_total = max_records or 0
        dash: Optional[Any] = None
        while True:
            params = dict(base_params or {})
            params["limit"] = page_size
            params["offset"] = offset
            page = self.get_data(endpoint, params)
            if not page:
                break
            if max_records is not None and len(all_rows) + len(page) > max_records:
                page = page[: max_records - len(all_rows)]
            all_rows.extend(page)
            page_num += 1
            if len(page) == page_size and estimated_total == 0:
                estimated_total = page_size * 20
            if self.verbose and PalaeoProgress:
                if dash is None:
                    est = max(estimated_total, len(all_rows))
                    dash = PalaeoProgress(est, label=progress_label, verbose=True)
                if estimated_total and len(all_rows) > estimated_total:
                    dash.total = len(all_rows) + page_size
                dash.set_current(
                    f"page {page_num}  ({len(all_rows):,} {progress_label} so far)",
                    "fetching",
                )
            elif self.verbose:
                print(f"  fetched {len(all_rows):,} {progress_label}...", file=sys.stderr)
            if len(page) < page_size:
                break
            if max_records is not None and len(all_rows) >= max_records:
                break
            offset += page_size
        if dash:
            dash.total = len(all_rows)
            dash.completed = len(all_rows)
            dash.finish(f"Fetched {len(all_rows):,} {progress_label} in {page_num} page(s)")
        return all_rows


def norm_key(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", value.strip().lower()).strip("_")


def bbox_to_wkt(min_lon: float, min_lat: float, max_lon: float, max_lat: float) -> str:
    # WKT polygon coordinates are lon lat and must close the ring.
    return (
        f"POLYGON (({min_lon} {min_lat}, {max_lon} {min_lat}, "
        f"{max_lon} {max_lat}, {min_lon} {max_lat}, {min_lon} {min_lat}))"
    )


def make_wkts(min_lon: float, min_lat: float, max_lon: float, max_lat: float) -> List[str]:
    """Return one or two WKT polygons.  Bboxes that span the full longitude
    range (antimeridian-crossing, e.g. holarctic or northern_hemisphere) cause
    the Neotoma spatial API to silently return 0 results.  Detect this case and
    split at meridian 0 into western + eastern sub-queries instead.

    Additionally, the Neotoma occurrences API silently returns 0 results when
    min_lon <= -179 or max_lon >= 180 (empirically determined).  Split
    boundaries are clamped to [-178, 178] to avoid this."""
    if min_lon <= -170.0 and max_lon >= 170.0:
        w_min = max(min_lon, -178.0)  # avoid Neotoma ≤-179 silent-zero bug
        e_max = min(max_lon,  178.0)  # avoid Neotoma ≥180 silent-zero bug
        return [
            bbox_to_wkt(w_min, min_lat, 0.0, max_lat),
            bbox_to_wkt(0.0, min_lat, e_max, max_lat),
        ]
    return [bbox_to_wkt(min_lon, min_lat, max_lon, max_lat)]


def geojson_to_wkt(path: str) -> str:
    with open(path, "r", encoding="utf-8") as fh:
        obj = json.load(fh)
    if obj.get("type") == "Feature":
        obj = obj.get("geometry", {})
    if obj.get("type") == "FeatureCollection":
        features = obj.get("features") or []
        if not features:
            raise ValueError("GeoJSON FeatureCollection has no features")
        obj = features[0].get("geometry", {})
    gtype = obj.get("type")
    coords = obj.get("coordinates")
    if gtype == "Polygon":
        return "POLYGON (" + ", ".join(_ring_to_wkt(ring) for ring in coords) + ")"
    if gtype == "MultiPolygon":
        parts = []
        for polygon in coords:
            parts.append("(" + ", ".join(_ring_to_wkt(ring) for ring in polygon) + ")")
        return "MULTIPOLYGON (" + ", ".join(parts) + ")"
    raise ValueError(f"Unsupported GeoJSON geometry type: {gtype}. Use Polygon or MultiPolygon.")


def _ring_to_wkt(ring: Sequence[Sequence[float]]) -> str:
    if not ring:
        raise ValueError("Empty polygon ring")
    pts = [(float(x), float(y)) for x, y, *rest in ring]
    if pts[0] != pts[-1]:
        pts.append(pts[0])
    return "(" + ", ".join(f"{x} {y}" for x, y in pts) + ")"


def parse_csv_list(values: Sequence[str]) -> List[str]:
    out: List[str] = []
    for value in values or []:
        for part in str(value).split(","):
            part = part.strip()
            if part:
                out.append(part)
    return out


def parse_organisms(values: Sequence[str]) -> List[str]:
    raw = parse_csv_list(values)
    if not raw:
        return ["all"]
    out = []
    aliases = {
        "animal": "animals", "animals": "animals", "fauna": "animals", "vertebrates": "animals", "vertebrate": "animals",
        "plant": "plants", "plants": "plants", "flora": "plants", "vascular_plants": "plants",
        "fungus": "fungi", "fungi": "fungi", "fungal": "fungi",
        "diatom": "diatoms", "diatoms": "diatoms",
        "protist": "protists", "protists": "protists",
        "all": "all", "any": "all", "everything": "all",
    }
    for item in raw:
        key = norm_key(item)
        if key not in aliases:
            raise ValueError(f"Unknown organism group: {item}. Try animals, plants, fungi, diatoms, protists, all.")
        val = aliases[key]
        if val not in out:
            out.append(val)
    if "all" in out:
        return ["all"]
    return out


def truthy(value: Any) -> Optional[bool]:
    if value is None:
        return None
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return bool(value)
    s = str(value).strip().lower()
    if s in {"1", "true", "t", "yes", "y", "extinct"}:
        return True
    if s in {"0", "false", "f", "no", "n", "extant"}:
        return False
    return None


def first_present(row: Dict[str, Any], *names: str) -> Any:
    for name in names:
        if name in row:
            return row[name]
        low = name.lower()
        if low in row:
            return row[low]
        cap = name[:1].upper() + name[1:]
        if cap in row:
            return row[cap]
    return None


def fetch_taxa_table(client: NeotomaClient, page_size: int = 100000) -> Dict[str, Dict[str, Any]]:
    print("Fetching Neotoma taxa table...", file=sys.stderr)
    rows = client.paged_get_data(
        "/data/dbtables/taxa",
        {"count": "false"},
        page_size=page_size,
        progress_label="taxa rows",
    )
    taxa: Dict[str, Dict[str, Any]] = {}
    for row in rows:
        tid = first_present(row, "taxonid", "TaxonID")
        if tid is not None:
            taxa[str(tid)] = row
    print(f"  taxa table loaded: {len(taxa):,} taxa", file=sys.stderr)
    return taxa


def infer_taxonomy(taxon_meta: Optional[Dict[str, Any]], datasettype: str = "", taxonname: str = "") -> Tuple[str, str, str, str]:
    """Return kingdom, phylum, organism_group, inference_source."""
    tg = ""
    if taxon_meta:
        tg = str(first_present(taxon_meta, "taxagroupid", "TaxaGroupID") or "").strip().upper()
        if tg in TAXAGROUP_MAP:
            kd, ph, group = TAXAGROUP_MAP[tg]
            return kd, ph, group, f"taxagroupid:{tg}"
    dt = (datasettype or "").lower()
    if any(x in dt for x in ["vertebrate", "fauna", "fish", "herp", "mammal", "bird"]):
        return "Animalia", "Chordata", "animals", "datasettype"
    if any(x in dt for x in ["mollusk", "mollusc"]):
        return "Animalia", "Mollusca", "animals", "datasettype"
    if any(x in dt for x in ["insect", "beetle", "ostracod", "cladocera", "chironomid"]):
        return "Animalia", "Arthropoda", "animals", "datasettype"
    if any(x in dt for x in ["pollen", "plant", "macrofossil"]):
        return "Plantae", "Tracheophyta", "plants", "datasettype"
    if "fung" in dt or "spore" in dt:
        return "Fungi", "", "fungi", "datasettype"
    if "diatom" in dt:
        return "Chromista", "Bacillariophyta", "diatoms", "datasettype"
    tn = (taxonname or "").strip()
    if re.search(r"aceae\b", tn, re.IGNORECASE):
        return "Plantae", "", "plants", "name_suffix"
    if re.search(r"idae\b", tn, re.IGNORECASE):
        return "Animalia", "", "animals", "name_suffix"
    return "Unknown", "", "unknown", "unknown"


def organism_allowed(group: str, requested: Sequence[str]) -> bool:
    if not requested or "all" in requested:
        return True
    return group in requested


def clean_base_name(name: str) -> str:
    s = str(name or "").strip()
    s = s.replace("×", "x")
    s = re.sub(r"\s+", " ", s)
    # Drop parenthetical size/context comments, but leave taxonomic subgenus-like words outside.
    s = re.sub(r"\([^)]*\)", " ", s)
    s = s.replace("?", " ")
    s = QUALIFIER_RE.sub(" ", s)
    s = BAD_CHARS_RE.sub(" ", s)
    s = re.sub(r"\s+", " ", s).strip(" ._-;,")
    return s


def is_higher_taxon_word(word: str) -> bool:
    w = word.lower().strip(" .")
    return len(w) > 3 and any(w.endswith(sfx) for sfx in HIGHER_TAXON_SUFFIXES)


def clean_taxon_name(name: str, *, ncbi_name_mode: str = "binomial", keep_genus_only: bool = True,
                     keep_higher_taxa: bool = False, include_morphotypes: bool = True) -> Tuple[List[CleanedName], Optional[str]]:
    original = str(name or "").strip()
    if not original:
        return [], "empty taxon name"
    if "/" in original:
        return [], "slash/alternative taxon name; rerun with --split-slash-names to split"

    lower = original.lower().strip()
    if lower in BAD_SINGLE_WORDS:
        return [], "uninformative taxon name"

    morphotype = bool(TYPE_RE.search(original))
    s = clean_base_name(original)
    if not s:
        return [], "empty after cleaning qualifiers"

    # For pollen/plant morphotypes such as Ambrosia-type, use genus bucket unless disabled.
    if morphotype:
        if not include_morphotypes:
            return [], "morphotype/type taxon excluded"
        s = re.sub(TYPE_RE, " ", s).strip(" -_.")

    tokens = re.findall(r"[A-Za-z][A-Za-z\-]*", s)
    if not tokens:
        return [], "no Latin-like tokens"

    genus = tokens[0]
    if genus.lower() in BAD_SINGLE_WORDS:
        return [], "uninformative first token"
    genus = genus[:1].upper() + genus[1:]

    genus_only_markers = bool(SP_RE.search(original)) or len(tokens) == 1
    if genus_only_markers:
        if is_higher_taxon_word(genus):
            if keep_higher_taxa:
                return [CleanedName(genus, genus, "higher", genus, "higher_taxon_kept")], None
            return [], "higher taxon excluded"
        if keep_genus_only:
            return [CleanedName(genus, genus, "genus", genus, "genus_bucket")], None
        return [], "genus-only taxon excluded"

    if is_higher_taxon_word(genus) and len(tokens) == 1:
        if keep_higher_taxa:
            return [CleanedName(genus, genus, "higher", genus, "higher_taxon_kept")], None
        return [], "higher taxon excluded"

    species_epithet = tokens[1].lower()
    if species_epithet in {"sp", "spp", "species", "undiff", "undifferentiated", "type"}:
        if keep_genus_only:
            return [CleanedName(genus, genus, "genus", genus, "genus_bucket")], None
        return [], "genus-only taxon excluded"
    if not re.match(r"^[a-z][a-z\-]+$", species_epithet):
        if keep_genus_only:
            return [CleanedName(genus, genus, "genus", genus, "fallback_genus_bucket_bad_species_epithet")], None
        return [], "bad species epithet"

    binomial = f"{genus} {species_epithet}"
    search_name = binomial
    rank = "species"
    status = "binomial"

    if ncbi_name_mode == "trinomial" and len(tokens) >= 3:
        subsp = tokens[2].lower()
        if re.match(r"^[a-z][a-z\-]+$", subsp) and subsp not in {"subsp", "ssp", "var", "forma", "cf", "aff"}:
            search_name = f"{binomial} {subsp}"
            rank = "subspecies"
            status = "trinomial"
    elif ncbi_name_mode == "as-is":
        search_name = s
        rank = "taxon"
        status = "cleaned_as_is"

    return [CleanedName(search_name, genus, rank, binomial, status)], None


def expand_slash_name(name: str) -> List[str]:
    """Very conservative slash expansion.

    Ostrya/Carpinus -> [Ostrya, Carpinus]
    Bison bison/latifrons -> [Bison bison, Bison latifrons]
    Mustelidae/Mephitidae -> [Mustelidae, Mephitidae]
    """
    s = str(name or "").strip()
    if "/" not in s:
        return [s]
    s_clean = clean_base_name(s)
    parts = [p.strip() for p in s_clean.split("/") if p.strip()]
    if len(parts) < 2:
        return [s]
    first_tokens = parts[0].split()
    out = [parts[0]]
    if len(first_tokens) == 2:
        genus = first_tokens[0]
        for part in parts[1:]:
            ptoks = part.split()
            if len(ptoks) == 1 and ptoks[0][:1].islower():
                out.append(f"{genus} {ptoks[0]}")
            else:
                out.append(part)
    else:
        out.extend(parts[1:])
    return out


def build_occurrence_params(wkt: Optional[str], ageyoung: Optional[int], ageold: Optional[int],
                            datasettype: Optional[str] = None) -> Dict[str, Any]:
    params: Dict[str, Any] = {}
    if wkt:
        params["loc"] = wkt
    if ageyoung is not None:
        params["ageyoung"] = int(ageyoung)
    if ageold is not None:
        params["ageold"] = int(ageold)
    if datasettype:
        params["datasettype"] = datasettype
    return params


def get_occurrences_region_first(client: NeotomaClient, params: Dict[str, Any], page_size: int,
                                 max_occurrences: Optional[int], datasettypes: Sequence[str],
                                 api_prefilter: str) -> List[Dict[str, Any]]:
    if datasettypes and api_prefilter in {"on", "auto"}:
        all_occ: List[Dict[str, Any]] = []
        seen_occids = set()
        failed = False
        for dt in datasettypes:
            p = dict(params)
            p["datasettype"] = dt
            print(f"Fetching occurrences with datasettype prefilter: {dt}", file=sys.stderr)
            try:
                rows = client.paged_get_data("/data/occurrences", p, page_size=page_size,
                                            max_records=max_occurrences, progress_label="occurrences")
                for row in rows:
                    occid = row.get("occid")
                    key = str(occid) if occid is not None else json.dumps(row, sort_keys=True)[:200]
                    if key not in seen_occids:
                        seen_occids.add(key)
                        all_occ.append(row)
                if max_occurrences is not None and len(all_occ) >= max_occurrences:
                    return all_occ[:max_occurrences]
            except NeotomaError as exc:
                print(f"WARNING: datasettype-prefilter query failed for {dt}: {exc}", file=sys.stderr)
                failed = True
                break
        if not failed:
            return all_occ
        if api_prefilter == "on":
            raise NeotomaError("API datasettype prefilter failed and --api-prefilter=on was requested")
        print("Retrying occurrence query without datasettype prefilter...", file=sys.stderr)

    print("Fetching occurrences...", file=sys.stderr)
    return client.paged_get_data("/data/occurrences", params, page_size=page_size,
                                 max_records=max_occurrences, progress_label="occurrences")


def occurrence_taxon_fields(occ: Dict[str, Any]) -> Tuple[str, str, str, str, Any, Any, Any, Any]:
    sample = occ.get("sample") or {}
    site = occ.get("site") or {}
    age = occ.get("age") or {}
    taxonid = str(sample.get("taxonid") or "")
    taxonname = str(sample.get("taxonname") or "")
    datasettype = str(site.get("datasettype") or "")
    database = str(site.get("database") or "")
    siteid = site.get("siteid")
    datasetid = site.get("datasetid")
    ageyounger = age.get("ageyounger")
    ageolder = age.get("ageolder")
    return taxonid, taxonname, datasettype, database, siteid, datasetid, ageyounger, ageolder


def status_allowed(taxon_meta: Optional[Dict[str, Any]], requested: str) -> bool:
    if requested == "all":
        return True
    ext = truthy(first_present(taxon_meta or {}, "extinct", "Extinct"))
    if requested == "extinct":
        return ext is True
    if requested == "extant":
        return ext is False
    return True


def safe_num(value: Any) -> Optional[float]:
    try:
        if value is None or value == "":
            return None
        v = float(value)
        if math.isnan(v):
            return None
        return v
    except Exception:
        return None


def aggregate_occurrences(occurrences: List[Dict[str, Any]], taxa: Dict[str, Dict[str, Any]],
                          requested_status: str, requested_organisms: Sequence[str],
                          ncbi_name_mode: str, keep_genus_only: bool, keep_higher_taxa: bool,
                          include_morphotypes: bool, split_slash_names: bool,
                          min_occurrences: int, min_sites: int) -> Tuple[List[Dict[str, Any]], List[RejectedName], Dict[str, Any]]:
    groups: Dict[str, Dict[str, Any]] = {}
    rejected: List[RejectedName] = []
    counters = defaultdict(int)

    for occ in occurrences:
        counters["occurrences_seen"] += 1
        taxonid, taxonname, datasettype, database, siteid, datasetid, ageyounger, ageolder = occurrence_taxon_fields(occ)
        meta = taxa.get(str(taxonid))
        if not status_allowed(meta, requested_status):
            counters["status_filtered"] += 1
            continue
        kd, ph, org_group, infer_src = infer_taxonomy(meta, datasettype=datasettype, taxonname=taxonname)
        if not organism_allowed(org_group, requested_organisms):
            counters["organism_filtered"] += 1
            continue

        candidate_names = expand_slash_name(taxonname) if split_slash_names else [taxonname]
        any_clean = False
        for cname in candidate_names:
            cleaned, reason = clean_taxon_name(
                cname,
                ncbi_name_mode=ncbi_name_mode,
                keep_genus_only=keep_genus_only,
                keep_higher_taxa=keep_higher_taxa,
                include_morphotypes=include_morphotypes,
            )
            if reason:
                rejected.append(RejectedName(taxonid, taxonname, reason, context=datasettype))
                counters["name_rejected"] += 1
                continue
            for cn in cleaned:
                any_clean = True
                key = cn.search_name.lower()
                if key not in groups:
                    ext = truthy(first_present(meta or {}, "extinct", "Extinct"))
                    groups[key] = {
                        "species": cn.search_name,  # GBIF-like but intentionally FlyGuide/NCBI search bucket
                        "genus": cn.genus,
                        "kingdom": kd,
                        "phylum": ph,
                        "scientificName": cn.search_name,
                        "canonicalBinomial": cn.canonical_binomial,
                        "taxonRank": cn.rank,
                        "taxonKey": taxonid,
                        "source": "Neotoma",
                        "neotoma_taxonids": set(),
                        "neotoma_taxonnames": set(),
                        "neotoma_taxagroupids": set(),
                        "neotoma_extinct_values": set(),
                        "neotoma_occurrence_count": 0,
                        "neotoma_siteids": set(),
                        "neotoma_datasetids": set(),
                        "neotoma_databases": set(),
                        "neotoma_datasettypes": set(),
                        "neotoma_min_age_young": None,
                        "neotoma_max_age_old": None,
                        "region_taxonomy_inference": infer_src,
                        "name_cleaning_status": cn.cleaning_status,
                        "ncbiSearchName": cn.search_name,
                    }
                g = groups[key]
                g["neotoma_taxonids"].add(taxonid)
                g["neotoma_taxonnames"].add(taxonname)
                tg = str(first_present(meta or {}, "taxagroupid", "TaxaGroupID") or "")
                if tg:
                    g["neotoma_taxagroupids"].add(tg)
                extv = first_present(meta or {}, "extinct", "Extinct")
                if extv is not None:
                    g["neotoma_extinct_values"].add(str(extv))
                g["neotoma_occurrence_count"] += 1
                if siteid is not None:
                    g["neotoma_siteids"].add(str(siteid))
                if datasetid is not None:
                    g["neotoma_datasetids"].add(str(datasetid))
                if database:
                    g["neotoma_databases"].add(database)
                if datasettype:
                    g["neotoma_datasettypes"].add(datasettype)
                ay = safe_num(ageyounger)
                ao = safe_num(ageolder)
                if ay is not None:
                    cur = g["neotoma_min_age_young"]
                    g["neotoma_min_age_young"] = ay if cur is None else min(cur, ay)
                if ao is not None:
                    cur = g["neotoma_max_age_old"]
                    g["neotoma_max_age_old"] = ao if cur is None else max(cur, ao)
        if any_clean:
            counters["occurrences_retained_after_name_clean"] += 1

    rows: List[Dict[str, Any]] = []
    for g in groups.values():
        if g["neotoma_occurrence_count"] < min_occurrences:
            counters["min_occurrence_filtered_taxa"] += 1
            continue
        if len(g["neotoma_siteids"]) < min_sites:
            counters["min_site_filtered_taxa"] += 1
            continue
        out = dict(g)
        out["neotoma_taxonids"] = ";".join(sorted(out["neotoma_taxonids"], key=str))
        out["neotoma_taxonnames"] = "; ".join(sorted(out["neotoma_taxonnames"]))
        out["neotoma_taxagroupids"] = ";".join(sorted(out["neotoma_taxagroupids"]))
        out["neotoma_extinct_values"] = ";".join(sorted(out["neotoma_extinct_values"]))
        out["neotoma_site_count"] = len(out["neotoma_siteids"])
        out["neotoma_dataset_count"] = len(out["neotoma_datasetids"])
        out["neotoma_siteids"] = ";".join(sorted(out["neotoma_siteids"], key=str))
        out["neotoma_datasetids"] = ";".join(sorted(out["neotoma_datasetids"], key=str))
        out["neotoma_databases"] = ";".join(sorted(out["neotoma_databases"]))
        out["neotoma_datasettypes"] = ";".join(sorted(out["neotoma_datasettypes"]))
        rows.append(out)
    rows.sort(key=lambda r: (r.get("kingdom", ""), r.get("phylum", ""), r.get("species", "")))
    counters["unique_output_taxa"] = len(rows)
    return rows, rejected, dict(counters)


def write_csv(path: str, rows: List[Dict[str, Any]]) -> None:
    fieldnames = [
        "species", "genus", "kingdom", "phylum", "scientificName", "canonicalBinomial",
        "taxonRank", "taxonKey", "source", "ncbiSearchName",
        "neotoma_taxonids", "neotoma_taxonnames", "neotoma_taxagroupids", "neotoma_extinct_values",
        "neotoma_occurrence_count", "neotoma_site_count", "neotoma_dataset_count",
        "neotoma_min_age_young", "neotoma_max_age_old", "neotoma_databases", "neotoma_datasettypes",
        "region", "period", "ageyoung", "ageold", "region_taxonomy_inference", "name_cleaning_status",
        "neotoma_siteids", "neotoma_datasetids",
    ]
    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_rejected(path: str, rejected: List[RejectedName]) -> None:
    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["taxonid", "taxonname", "reason", "context"], delimiter="\t")
        writer.writeheader()
        for r in rejected:
            writer.writerow(dataclasses.asdict(r))


def write_flyguide_files(out_prefix: str, rows: List[Dict[str, Any]]) -> Tuple[str, str]:
    search_path = f"{out_prefix}_species_search.txt"
    kingdom_path = f"{out_prefix}_species_kingdom.tsv"
    seen = set()
    with open(search_path, "w", encoding="utf-8") as fh:
        for row in rows:
            name = row.get("ncbiSearchName") or row.get("species")
            if name and name not in seen:
                seen.add(name)
                fh.write(str(name) + "\n")
    with open(kingdom_path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        for row in rows:
            name = row.get("ncbiSearchName") or row.get("species")
            if name:
                writer.writerow([name, row.get("kingdom", "Unknown"), row.get("phylum", "")])
    return search_path, kingdom_path


def list_regions() -> None:
    print("Built-in region presets:\n")
    for key in sorted(REGION_BBOXES):
        minlon, minlat, maxlon, maxlat, desc = REGION_BBOXES[key]
        print(f"  {key:22s} bbox=({minlon}, {minlat}, {maxlon}, {maxlat})  {desc}")
    print("\nUse --bbox minLon,minLat,maxLon,maxLat or --geojson polygon.geojson for custom regions.")


def list_periods() -> None:
    print("Built-in period presets:\n")
    for key in sorted(PERIODS):
        young, old, desc = PERIODS[key]
        if young is None:
            print(f"  {key:20s} no age filter  {desc}")
        else:
            print(f"  {key:20s} ageyoung={young:<9} ageold={old:<9} {desc}")


def make_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Create a FlyGuide/GBIF-like CSV from Neotoma occurrence data.",
        epilog=textwrap.dedent("""
        Examples:
          python neotoma_extinct_to_gbif.py --region northern_hemisphere --period quaternary --organisms animals --status extinct --out NH_quaternary_animals.csv --write-flyguide-files
          python neotoma_extinct_to_gbif.py --region holarctic --period pleistocene --organisms animals --ncbi-name-mode binomial --out holarctic_pleistocene_animals.csv
          python neotoma_extinct_to_gbif.py --bbox -140,40,-50,85 --period holocene --organisms plants --status all --out canada_like_holocene_plants.csv
        """),
    )
    p.add_argument("--version", action="version", version=VERSION)
    p.add_argument("--list-regions", action="store_true", help="Print built-in region names and exit")
    p.add_argument("--list-periods", action="store_true", help="Print built-in period names and exit")
    p.add_argument("--region", default="global", help="Built-in region name, e.g. northern_hemisphere, holarctic, north_america, eurasia, tanzania")
    p.add_argument("--bbox", help="Custom bbox as minLon,minLat,maxLon,maxLat; overrides --region")
    p.add_argument("--geojson", help="Custom GeoJSON Polygon/MultiPolygon file; overrides --bbox and --region")
    p.add_argument("--period", default="quaternary", help="Period preset: all, quaternary, pleistocene, holocene, late_pleistocene, pliocene, etc.")
    p.add_argument("--age-young", type=int, help="Override youngest age BP; passed to Neotoma as ageyoung")
    p.add_argument("--age-old", type=int, help="Override oldest age BP; passed to Neotoma as ageold")
    p.add_argument("--organisms", action="append", default=[], help="Organism groups: animals, plants, fungi, diatoms, protists, all. Comma-separated or repeated.")
    p.add_argument("--status", choices=["extinct", "extant", "all"], default="extinct", help="Use Neotoma taxa.extinct status filter. Default: extinct")
    p.add_argument("--datasettype", action="append", default=[], help="Explicit Neotoma datasettype(s), comma-separated or repeated. Used as API prefilter if possible and in output metadata.")
    p.add_argument("--api-prefilter", choices=["auto", "on", "off"], default="auto", help="Use datasettype hints at the API level for speed. auto retries without them if rejected. Default: auto")
    p.add_argument("--ncbi-name-mode", choices=["binomial", "trinomial", "as-is"], default="binomial", help="Name bucket to export for NCBI search. Default collapses subspecies/splitter categories to binomial/genus buckets.")
    p.add_argument("--drop-genus-only", action="store_true", help="Do not keep genus-only buckets from sp./spp./undiff. names")
    p.add_argument("--keep-higher-taxa", action="store_true", help="Allow family/order/higher buckets such as Felidae or Poaceae")
    p.add_argument("--drop-morphotypes", action="store_true", help="Drop pollen/morphotype/type names instead of converting to genus buckets")
    p.add_argument("--split-slash-names", action="store_true", help="Split simple slash names such as Ostrya/Carpinus into separate buckets")
    p.add_argument("--min-occurrences", type=int, default=1, help="Minimum Neotoma occurrence rows required per output bucket")
    p.add_argument("--min-sites", type=int, default=1, help="Minimum unique Neotoma sites required per output bucket")
    p.add_argument("--page-size", type=int, default=5000, help="Neotoma API page size. Default: 5000")
    p.add_argument("--max-occurrences", type=int, default=0, help="Debug/safety cap on occurrences to fetch. 0 means unlimited")
    p.add_argument("--taxa-page-size", type=int, default=100000, help="Page size for taxa table fetch. Default: 100000")
    p.add_argument("--base-url", default=DEFAULT_BASE_URL, help="Neotoma API base URL")
    p.add_argument("--cache-dir", default=".neotoma_cache", help="Cache directory for API JSON pages. Use '' to disable.")
    p.add_argument("--force-refresh", action="store_true", help="Ignore cached Neotoma responses")
    p.add_argument("--sleep", type=float, default=0.15, help="Seconds to sleep between API calls")
    p.add_argument("--timeout", type=int, default=60, help="HTTP timeout per request")
    p.add_argument("--quiet", action="store_true", help="Less progress output")
    p.add_argument("--offline-occurrences-json", help="Testing: read occurrence JSON from a local Neotoma response/list instead of API")
    p.add_argument("--offline-taxa-json", help="Testing: read taxa JSON from a local Neotoma response/list instead of API")
    p.add_argument("--out", required=False, help="Output GBIF-like CSV")
    p.add_argument("--write-flyguide-files", action="store_true", help="Also write OUTPREFIX_species_search.txt and OUTPREFIX_species_kingdom.tsv")
    p.add_argument("--out-prefix", help="Prefix for --write-flyguide-files. Defaults to --out without extension")
    return p


def resolve_region(args: argparse.Namespace) -> Tuple[List[Optional[str]], str, str]:
    if args.geojson:
        return [geojson_to_wkt(args.geojson)], f"geojson:{args.geojson}", "custom GeoJSON"
    if args.bbox:
        parts = [float(x.strip()) for x in args.bbox.split(",")]
        if len(parts) != 4:
            raise ValueError("--bbox must be minLon,minLat,maxLon,maxLat")
        return make_wkts(*parts), f"bbox:{args.bbox}", "custom bbox"
    key_raw = str(args.region or "global").strip()
    key = norm_key(key_raw)
    key = REGION_ALIASES.get(key_raw.lower(), key)
    if key not in REGION_BBOXES:
        raise ValueError(f"Unknown region '{args.region}'. Run --list-regions.")
    minlon, minlat, maxlon, maxlat, desc = REGION_BBOXES[key]
    return make_wkts(minlon, minlat, maxlon, maxlat), key, desc


def resolve_period(args: argparse.Namespace) -> Tuple[Optional[int], Optional[int], str, str]:
    key = norm_key(args.period or "quaternary")
    if key not in PERIODS:
        raise ValueError(f"Unknown period '{args.period}'. Run --list-periods.")
    young, old, desc = PERIODS[key]
    if args.age_young is not None:
        young = args.age_young
    if args.age_old is not None:
        old = args.age_old
    if young is not None and old is not None and young > old:
        raise ValueError(f"age-young must be <= age-old. Got ageyoung={young}, ageold={old}.")
    return young, old, key, desc


def load_json_data(path: str) -> List[Dict[str, Any]]:
    with open(path, "r", encoding="utf-8") as fh:
        obj = json.load(fh)
    if isinstance(obj, dict) and "data" in obj:
        obj = obj["data"]
    if not isinstance(obj, list):
        raise ValueError(f"Expected JSON list or Neotoma response with data list in {path}")
    return obj


def datasettype_hints(requested_organisms: Sequence[str], explicit_datasettypes: Sequence[str]) -> List[str]:
    explicit = parse_csv_list(explicit_datasettypes)
    if explicit:
        return explicit
    if not requested_organisms or "all" in requested_organisms:
        return []
    hints: List[str] = []
    for org in requested_organisms:
        for dt in ORGANISM_DATASET_HINTS.get(org, []):
            if dt not in hints:
                hints.append(dt)
    return hints


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = make_argparser()
    args = parser.parse_args(argv)

    if args.list_regions:
        list_regions()
        return 0
    if args.list_periods:
        list_periods()
        return 0
    if not args.out:
        parser.error("--out is required unless --list-regions or --list-periods is used")

    wkts, region_label, region_desc = resolve_region(args)
    ageyoung, ageold, period_label, period_desc = resolve_period(args)
    organisms = parse_organisms(args.organisms)
    dtypes = datasettype_hints(organisms, args.datasettype)

    max_occ = args.max_occurrences if args.max_occurrences and args.max_occurrences > 0 else None
    cache_dir = args.cache_dir if args.cache_dir else None

    if not args.quiet:
        if print_header:
            params_list = [
                ("Region", f"{region_label} ({region_desc})"),
                ("Period", f"{period_label}  (ageyoung={ageyoung}, ageold={ageold})"),
                ("Organisms", ", ".join(organisms)),
                ("Status", args.status),
                ("Name mode", args.ncbi_name_mode),
                ("Output", args.out),
            ]
            if len(wkts) > 1:
                params_list.insert(1, ("Note", f"wide bbox split into {len(wkts)} sub-queries (antimeridian workaround)"))
            if dtypes:
                params_list.append(("Dataset types", ", ".join(dtypes)))
            print_header(f"FlyGuide Neotoma exporter  v{VERSION}", params_list)
        else:
            print("FlyGuide Neotoma exporter", file=sys.stderr)
            print(f"  region   : {region_label} ({region_desc})", file=sys.stderr)
            if len(wkts) > 1:
                print(f"  note     : wide bbox split into {len(wkts)} sub-queries (antimeridian workaround)", file=sys.stderr)
            print(f"  period   : {period_label} ({period_desc}); ageyoung={ageyoung}, ageold={ageold}", file=sys.stderr)
            print(f"  organisms: {', '.join(organisms)}", file=sys.stderr)
            print(f"  status   : {args.status}", file=sys.stderr)
            print(f"  name mode: {args.ncbi_name_mode}", file=sys.stderr)
            if dtypes:
                print(f"  datasettype API hints: {', '.join(dtypes)}", file=sys.stderr)

    if args.offline_taxa_json:
        taxa_rows = load_json_data(args.offline_taxa_json)
        taxa = {str(first_present(r, "taxonid", "TaxonID")): r for r in taxa_rows if first_present(r, "taxonid", "TaxonID") is not None}
    else:
        client = NeotomaClient(
            base_url=args.base_url,
            cache_dir=cache_dir,
            force_refresh=args.force_refresh,
            sleep=args.sleep,
            retries=3,
            timeout=args.timeout,
            verbose=not args.quiet,
        )
        taxa = fetch_taxa_table(client, page_size=args.taxa_page_size)

    if args.offline_occurrences_json:
        occurrences = load_json_data(args.offline_occurrences_json)
    else:
        client = NeotomaClient(
            base_url=args.base_url,
            cache_dir=cache_dir,
            force_refresh=args.force_refresh,
            sleep=args.sleep,
            retries=3,
            timeout=args.timeout,
            verbose=not args.quiet,
        )
        seen_occ_keys: set = set()
        occurrences: List[Dict[str, Any]] = []
        for wkt in wkts:
            chunk = get_occurrences_region_first(
                client, build_occurrence_params(wkt, ageyoung, ageold),
                page_size=args.page_size, max_occurrences=max_occ,
                datasettypes=dtypes, api_prefilter=args.api_prefilter,
            )
            for row in chunk:
                occid = row.get("occid")
                occ_key = str(occid) if occid is not None else json.dumps(row, sort_keys=True)[:200]
                if occ_key not in seen_occ_keys:
                    seen_occ_keys.add(occ_key)
                    occurrences.append(row)

    rows, rejected, counters = aggregate_occurrences(
        occurrences,
        taxa,
        requested_status=args.status,
        requested_organisms=organisms,
        ncbi_name_mode=args.ncbi_name_mode,
        keep_genus_only=not args.drop_genus_only,
        keep_higher_taxa=args.keep_higher_taxa,
        include_morphotypes=not args.drop_morphotypes,
        split_slash_names=args.split_slash_names,
        min_occurrences=args.min_occurrences,
        min_sites=args.min_sites,
    )

    for row in rows:
        row["region"] = region_label
        row["period"] = period_label
        row["ageyoung"] = "" if ageyoung is None else ageyoung
        row["ageold"] = "" if ageold is None else ageold

    out_csv = args.out
    write_csv(out_csv, rows)
    out_base = args.out_prefix or os.path.splitext(out_csv)[0]
    rejected_path = out_base + ".rejected.tsv"
    summary_path = out_base + ".summary.json"
    write_rejected(rejected_path, rejected)

    summary = {
        "version": VERSION,
        "region": region_label,
        "region_description": region_desc,
        "period": period_label,
        "period_description": period_desc,
        "ageyoung": ageyoung,
        "ageold": ageold,
        "organisms": organisms,
        "status": args.status,
        "ncbi_name_mode": args.ncbi_name_mode,
        "occurrences_fetched": len(occurrences),
        "taxa_table_rows": len(taxa),
        "output_taxa": len(rows),
        "rejected_names": len(rejected),
        "counters": counters,
        "outputs": {"csv": out_csv, "rejected": rejected_path, "summary": summary_path},
    }

    if args.write_flyguide_files:
        search_path, kingdom_path = write_flyguide_files(out_base, rows)
        summary["outputs"]["species_search"] = search_path
        summary["outputs"]["species_kingdom"] = kingdom_path

    with open(summary_path, "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)

    if print_done:
        counts = [
            ("Occurrences fetched:", f"{len(occurrences):,}"),
            ("Output taxon buckets:", f"{len(rows):,}"),
            ("Rejected names:", f"{len(rejected):,}  →  {rejected_path}"),
            ("CSV:", out_csv),
            ("Summary:", summary_path),
        ]
        if args.write_flyguide_files:
            counts += [
                ("FlyGuide search:", summary["outputs"]["species_search"]),
                ("FlyGuide kingdom map:", summary["outputs"]["species_kingdom"]),
            ]
        print_done("Neotoma export complete", counts)
    else:
        print("\nNeotoma export complete", file=sys.stderr)
        print(f"  Occurrences fetched : {len(occurrences):,}", file=sys.stderr)
        print(f"  Output taxon buckets: {len(rows):,}", file=sys.stderr)
        print(f"  Rejected names      : {len(rejected):,}", file=sys.stderr)
        print(f"  CSV                 : {out_csv}", file=sys.stderr)
        print(f"  Rejected log        : {rejected_path}", file=sys.stderr)
        print(f"  Summary             : {summary_path}", file=sys.stderr)
        if args.write_flyguide_files:
            print(f"  FlyGuide search     : {summary['outputs']['species_search']}", file=sys.stderr)
            print(f"  FlyGuide kingdom map: {summary['outputs']['species_kingdom']}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
