"""Optional common-name enrichment helpers for MetaMerge.

Common-name resolution is a convenience feature and does **not** affect any
classification decision.  Resolution proceeds through five layers, stopping at
the first layer that returns a name:

  1. User-supplied override file (CSV with scientific_name/tax_id/common_name).
  2. Small built-in convenience map in the package defaults.
  3. Best-effort online lookup via NCBI Taxonomy eutils (when tax_id is known).
  4. Best-effort online lookup via the GBIF Species API (by scientific name).
  5. Best-effort online lookup via the iNaturalist API (by scientific name).

Layers 3–5 are only attempted when ``--online-common-names`` is passed.
This design keeps the package fully functional offline and reproducible even
when any individual API is unavailable or changes.

All online results use **English-only** filtering — non-English vernacular names
from GBIF are not used as fallbacks (they were previously, causing e.g. French
names to appear).  NCBI and iNaturalist both return English names natively.

Online lookup notes
---------------------
NCBI Taxonomy (``efetch``):
  ``GET /entrez/eutils/efetch.fcgi?db=taxonomy&id=<tax_id>&retmode=xml``
  Parses ``<GenbankCommonName>`` first, then ``<OtherNames><CommonName>``.
  Only attempted when a numeric ``tax_id`` is available.

GBIF Species API:
  1. ``GET /v1/species/match?name=<scientific_name>`` to resolve a GBIF usage key.
  2. ``GET /v1/species/<usageKey>/vernacularNames`` for vernacular names.
  Returns English names only (language == "eng" or "en").

iNaturalist Taxon API:
  ``GET /v1/taxa?q=<scientific_name>&locale=en&per_page=5``
  Uses ``preferred_common_name`` from the top-ranked result that matches the
  queried scientific name closely.

Results from all sources are cached in a local JSON file (written to the
output directory) so that repeated runs do not re-query APIs for already-
looked-up taxa.  The cache stores ``{"name_or_source": {"name": ..., "src": ...}}``.

Override file format
---------------------
A plain CSV with columns: ``scientific_name``, ``tax_id``, ``common_name``.
Both scientific_name and tax_id are optional in each row, but at least one
must be provided.  See ``examples/common_names_template.csv`` for an example.
"""

from __future__ import annotations

import csv
import json
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Optional

import requests

from .utils import normalize_name


GBIF_MATCH_URL      = "https://api.gbif.org/v1/species/match"
GBIF_VERNACULAR_URL = "https://api.gbif.org/v1/species/{usageKey}/vernacularNames"
NCBI_EFETCH_URL     = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
INAT_TAXA_URL       = "https://api.inaturalist.org/v1/taxa"

_ENGLISH_LANG_CODES = {"eng", "en", "english"}


def load_common_name_overrides(path: str | None) -> dict:
    """Load a user-supplied common-name override CSV.

    Args:
        path: Path to a CSV with columns scientific_name, tax_id, common_name.
            Pass ``None`` to skip (returns an empty dict).

    Returns:
        Dict mapping ``("name", normalized_sci_name)`` or ``("tax_id", tax_id_str)``
        to the common-name string.
    """
    if not path:
        return {}
    table = {}
    with Path(path).open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            sci    = normalize_name(row.get("scientific_name") or row.get("tax_name"))
            tid    = str(row.get("tax_id") or "").strip()
            common = (row.get("common_name") or "").strip()
            if common:
                if sci:
                    table[("name", sci)] = common
                if tid:
                    table[("tax_id", tid)] = common
    return table


def load_common_name_cache(path: Path) -> dict:
    """Load a JSON cache of previous online lookups.

    Args:
        path: Path to the JSON cache file.

    Returns:
        Cache dict, or an empty dict if the file does not yet exist.
    """
    if path.exists():
        with path.open("r", encoding="utf-8") as handle:
            return json.load(handle)
    return {}


def save_common_name_cache(path: Path, cache: dict) -> None:
    """Persist the online common-name cache to disk.

    Args:
        path: Destination JSON path.
        cache: The cache dict to write.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(cache, handle, indent=2, sort_keys=True)


def _cache_get(cache: dict, key: str) -> tuple[Optional[str], str] | None:
    """Return ``(name, source)`` from cache, or ``None`` if key not present."""
    if key not in cache:
        return None
    entry = cache[key]
    if entry is None or entry == "":
        return None, "cache-miss"
    if isinstance(entry, dict):
        return entry.get("name") or None, entry.get("src", "cache")
    # Legacy flat string entries
    return entry or None, "cache"


def _cache_set(cache: dict, key: str, name: Optional[str], src: str) -> None:
    """Store a name + source in the cache."""
    cache[key] = {"name": name or "", "src": src}


def lookup_common_name_ncbi(
    tax_id: object,
    cache: dict,
) -> tuple[Optional[str], str]:
    """Look up a common name from the NCBI Taxonomy database using a tax_id.

    Queries ``eutils/efetch`` for the GenBank common name field.  Only
    attempted when a numeric ``tax_id`` is available; falls back gracefully
    for any network or parse error.

    Args:
        tax_id: NCBI taxonomy ID (int or string with digits).
        cache: Mutable in-memory lookup cache (keyed by ``"ncbi:<tax_id>"``).

    Returns:
        Tuple of ``(common_name_or_None, source_string)`` where source_string
        is one of ``"cache"``, ``"ncbi"``, ``"ncbi-no-match"``, or ``"ncbi-error"``.
    """
    tid = str(tax_id).strip() if tax_id is not None else ""
    if not tid or not tid.isdigit():
        return None, "ncbi-no-tax-id"

    cache_key = f"ncbi:{tid}"
    cached = _cache_get(cache, cache_key)
    if cached is not None:
        return cached

    try:
        resp = requests.get(
            NCBI_EFETCH_URL,
            params={"db": "taxonomy", "id": tid, "retmode": "xml"},
            timeout=15,
        )
        resp.raise_for_status()
        root = ET.fromstring(resp.text)

        # GenBank common name is the most curated field (usually English).
        gcn = root.find(".//GenbankCommonName")
        if gcn is not None and gcn.text and gcn.text.strip():
            name = gcn.text.strip()
            _cache_set(cache, cache_key, name, "ncbi")
            return name, "ncbi"

        # OtherNames/CommonName is second choice.
        for cn in root.findall(".//OtherNames/CommonName"):
            if cn.text and cn.text.strip():
                name = cn.text.strip()
                _cache_set(cache, cache_key, name, "ncbi")
                return name, "ncbi"

        _cache_set(cache, cache_key, None, "ncbi-no-match")
        return None, "ncbi-no-match"

    except Exception:
        _cache_set(cache, cache_key, None, "ncbi-error")
        return None, "ncbi-error"


def lookup_common_name_gbif(
    scientific_name: str,
    cache: dict,
) -> tuple[Optional[str], str]:
    """Best-effort GBIF vernacular name lookup for a single taxon.

    Returns **English names only** (language == "eng" or "en").  Non-English
    names from GBIF are deliberately ignored — they caused spurious foreign
    names to appear in previous versions.

    Results are written to ``cache`` keyed by ``"gbif:<normalized_name>"``
    so repeated calls are served from memory and can be persisted between runs.

    Args:
        scientific_name: Scientific name to look up.
        cache: Mutable dict used as an in-memory lookup cache.

    Returns:
        Tuple of ``(common_name_or_None, source_string)`` where source_string
        is one of ``"cache"``, ``"gbif"``, ``"gbif-no-match"``,
        ``"gbif-no-english-name"``, or ``"gbif-error"``.
    """
    key_norm = normalize_name(scientific_name)
    if not key_norm:
        return None, "missing-scientific-name"

    cache_key = f"gbif:{key_norm}"
    cached = _cache_get(cache, cache_key)
    if cached is not None:
        return cached

    try:
        match    = requests.get(GBIF_MATCH_URL, params={"name": scientific_name}, timeout=15)
        match.raise_for_status()
        payload   = match.json()
        usage_key = payload.get("usageKey")
        if not usage_key:
            _cache_set(cache, cache_key, None, "gbif-no-match")
            return None, "gbif-no-match"

        names_resp = requests.get(
            GBIF_VERNACULAR_URL.format(usageKey=usage_key), timeout=15
        )
        names_resp.raise_for_status()
        data = names_resp.json().get("results", [])

        for item in data:
            vern = (item.get("vernacularName") or "").strip()
            lang = (item.get("language") or "").lower().strip()
            if vern and lang in _ENGLISH_LANG_CODES:
                _cache_set(cache, cache_key, vern, "gbif")
                return vern, "gbif"

        _cache_set(cache, cache_key, None, "gbif-no-english-name")
        return None, "gbif-no-english-name"

    except Exception:
        _cache_set(cache, cache_key, None, "gbif-error")
        return None, "gbif-error"


def lookup_common_name_inaturalist(
    scientific_name: str,
    cache: dict,
) -> tuple[Optional[str], str]:
    """Look up a common name from the iNaturalist Taxon API.

    iNaturalist typically returns English common names via ``preferred_common_name``
    and has broad coverage across plants, animals, and fungi.

    Results are cached keyed by ``"inat:<normalized_name>"``.

    Args:
        scientific_name: Scientific name to look up.
        cache: Mutable in-memory lookup cache.

    Returns:
        Tuple of ``(common_name_or_None, source_string)`` where source_string
        is one of ``"cache"``, ``"inaturalist"``, ``"inaturalist-no-match"``,
        or ``"inaturalist-error"``.
    """
    key_norm = normalize_name(scientific_name)
    if not key_norm:
        return None, "missing-scientific-name"

    cache_key = f"inat:{key_norm}"
    cached = _cache_get(cache, cache_key)
    if cached is not None:
        return cached

    try:
        resp = requests.get(
            INAT_TAXA_URL,
            params={"q": scientific_name, "locale": "en", "per_page": 5},
            timeout=15,
        )
        resp.raise_for_status()
        results = resp.json().get("results", [])

        query_norm = key_norm.lower()
        for result in results:
            result_name = normalize_name(result.get("name", "")).lower()
            if result_name != query_norm:
                continue
            common = (result.get("preferred_common_name") or "").strip()
            if common:
                _cache_set(cache, cache_key, common, "inaturalist")
                return common, "inaturalist"

        _cache_set(cache, cache_key, None, "inaturalist-no-match")
        return None, "inaturalist-no-match"

    except Exception:
        _cache_set(cache, cache_key, None, "inaturalist-error")
        return None, "inaturalist-error"


# Keep for backward compatibility with tests / external callers that imported
# the old name directly.
def lookup_common_name_online(
    scientific_name: str,
    cache: dict,
    language: str = "eng",
) -> tuple[Optional[str], str]:
    """Deprecated alias: try GBIF then iNaturalist by scientific name.

    Prefer ``lookup_common_name_gbif`` or ``lookup_common_name_inaturalist``
    directly.  This wrapper is retained for backward compatibility.

    Args:
        scientific_name: Scientific name to look up.
        cache: Mutable in-memory lookup cache.
        language: Ignored (English-only is now enforced).

    Returns:
        ``(common_name_or_None, source_string)`` from the first source that
        returns an English name, or ``(None, "online-no-match")``.
    """
    name, src = lookup_common_name_gbif(scientific_name, cache)
    if name:
        return name, src
    name, src = lookup_common_name_inaturalist(scientific_name, cache)
    if name:
        return name, src
    return None, "online-no-match"


def resolve_common_name(
    scientific_name: str,
    tax_id: object,
    builtin_map: dict,
    override_map: dict,
    online: bool,
    cache: dict,
    language: str = "eng",
) -> tuple[Optional[str], str]:
    """Resolve a common name using the five-layer resolution hierarchy.

    Resolution order:
      1. ``override_map`` keyed by tax_id (exact numeric match).
      2. ``override_map`` keyed by normalized scientific name.
      3. ``builtin_map`` keyed by normalized scientific name.
      4. Online lookup (only if ``online=True``), tried in order:
         a. NCBI Taxonomy eutils (by tax_id, when available).
         b. GBIF Species API (by scientific name, English only).
         c. iNaturalist Taxon API (by scientific name, English preferred_name).
      5. ``(None, "not-requested")`` if online=False and no match found.

    Args:
        scientific_name: Scientific name of the taxon.
        tax_id: NCBI taxonomy ID (int or string), or None.
        builtin_map: Small built-in convenience map from defaults.
        override_map: User-supplied override dict (from load_common_name_overrides).
        online: Whether to allow online API lookups.
        cache: Mutable online-lookup cache dict.
        language: Preferred language (kept for API compatibility; English is
            always preferred regardless of this setting).

    Returns:
        Tuple of ``(common_name_or_None, source_string)``.
    """
    sci = normalize_name(scientific_name)
    tid = str(tax_id).strip() if tax_id is not None else ""

    # 1. User override by tax_id
    if tid and ("tax_id", tid) in override_map:
        return override_map[("tax_id", tid)], "override-tax-id"

    # 2. User override by normalized name
    if sci and ("name", sci) in override_map:
        return override_map[("name", sci)], "override-name"

    # 3. Built-in convenience map
    if sci in builtin_map:
        return builtin_map[sci], "builtin"

    if not online:
        return None, "not-requested"

    # 4a. NCBI eutils (by tax_id — most reliable for NCBI-curated English names)
    if tid and tid.isdigit():
        name, src = lookup_common_name_ncbi(tid, cache)
        if name:
            return name, src

    # 4b. GBIF (by scientific name — English vernacular names only)
    name, src = lookup_common_name_gbif(scientific_name, cache)
    if name:
        return name, src

    # 4c. iNaturalist (by scientific name — broad coverage, English preferred)
    name, src = lookup_common_name_inaturalist(scientific_name, cache)
    if name:
        return name, src

    return None, "online-no-match"
