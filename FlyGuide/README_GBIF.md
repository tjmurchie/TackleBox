# FlyGuide: GBIF Live Query (`gbif_query.py`)

_Part of TackleBox FlyGuide — v0.1.0-gbif_

---

## Overview

`gbif_query.py` lets you build FlyGuide species lists directly from GBIF occurrence
data on the command line, without navigating the GBIF website.  You specify a taxon
group and a geographic area; the tool queries the GBIF REST API, resolves species
names, and writes a GBIF-like CSV that feeds straight into `gbif_prep_from_csv.py`
and the rest of the FlyGuide pipeline.

No GBIF account or API key is required.  The tool uses only the public read-only
GBIF occurrence API, caches responses to disk, and has no Python dependencies
beyond the standard library.

---

## How it works

GBIF's occurrence API does not return a direct species list.  Instead, `gbif_query.py`
uses a two-step approach:

1. **Facet query** — asks the API for a `speciesKey` facet (integer keys, up to
   10 000 per request).  This returns how many occurrences each species has within
   the query area, without downloading individual records.

2. **Name resolution** — resolves each `speciesKey` to a canonical name, genus,
   kingdom, and phylum via `/v1/species/{key}`.  This is done in parallel (10
   threads by default) and every response is cached to disk, so a second run over
   the same taxon group is near-instant.

When a query contains more than 10 000 distinct species (e.g. all animals globally),
the tool automatically **chunks by taxonomic order** (or by family if orders are
unavailable), runs a separate facet query per chunk, and merges the results.  This
is transparent — you do not need to do anything differently.

---

## Output files

| File | Description |
|---|---|
| `<out>.csv` | Main FlyGuide-compatible GBIF-like CSV |
| `<out>.rejected.tsv` | Names that could not be cleaned to a binomial |
| `<prefix>_species_search.txt` | One name per line; consumed by `NCBI-NT_Downloader.pl` |
| `<prefix>_species_kingdom.tsv` | Name TAB kingdom TAB phylum; used for FASTA splitting |

The `_species_search.txt` and `_species_kingdom.tsv` files are only written when
`--write-flyguide-files` is passed.

### CSV columns

| Column | Description |
|---|---|
| `species` | Cleaned binomial (or genus if no species epithet) |
| `genus` | Genus name |
| `kingdom` | Kingdom (Animalia, Plantae, Fungi, …) |
| `phylum` | Phylum |
| `scientificName` | Same as `species` |
| `taxonRank` | `species` or `genus_or_higher` |
| `source` / `source_database` | `GBIF` |
| `ncbi_search_name` | Name as submitted to NCBI |
| `canonical_binomial` | First two tokens of `species` |
| `gbif_occurrence_count` | Number of GBIF occurrences matching the query |
| `gbif_taxon_name` | Raw name from GBIF (before cleaning) |
| `region` | Region preset used (or `global`) |
| `country` | Country filter(s) applied |
| `continent` | Continent filter(s) applied |
| `basis_of_record` | Basis-of-record filter applied |
| `name_cleaning_status` | `ok` (unchanged) or `cleaned` |

---

## Quickstart

```bash
# List available taxon groups, regions, and basis-of-record options
python3 gbif_query.py --list-taxon-groups
python3 gbif_query.py --list-regions
python3 gbif_query.py --list-basis

# Terrestrial northern hemisphere mammals (human observations only)
python3 gbif_query.py \
  --taxon mammals \
  --region northern_hemisphere \
  --basis human \
  --out NH_mammals_gbif.csv \
  --write-flyguide-files

# Feed straight into FlyGuide
./flyguide.sh \
  NH_mammals_gbif.csv \
  NH_mammals_refs \
  your@email.org \
  YOUR_NCBI_API_KEY
```

---

## All flags

### Taxon (required, choose one)

| Flag | Description |
|---|---|
| `--taxon GROUP` | Preset name (`mammals`, `birds`, `plants`, …) or any scientific name (`Bovidae`, `Nymphaeales`, `Cetacea`, …) |
| `--taxon-key KEY` | Numeric GBIF taxon key (e.g. `359` for Mammalia) |

### Geographic filters

| Flag | Description |
|---|---|
| `--region PRESET` | Named geographic preset — bounding box applied to all records |
| `--bbox minLon,minLat,maxLon,maxLat` | Custom bounding box |
| `--country CC[,CC,…]` | ISO 2-letter country codes, comma-separated; OR'd |
| `--continent NAME[,NAME,…]` | Continent names, comma-separated; OR'd |

Geographic flags can be combined.  `--region` and `--bbox` set a bounding box.
`--country` and `--continent` add additional GBIF filters on top of that.
When `--region` is an ocean basin preset, no country filter is added — open-ocean
records carry no country code in GBIF and are captured by the bounding box alone.

### Record filters

| Flag | Default | Description |
|---|---|---|
| `--basis PRESET` | `human` | Basis-of-record filter — see presets below |
| `--year START,END` | — | Year range, e.g. `2000,2024` |
| `--require-coords` | off | Only records with coordinates and no geospatial issues |
| `--min-occurrences N` | `1` | Drop species with fewer than N occurrences in the query area |

### Output

| Flag | Description |
|---|---|
| `--out FILE` | Output GBIF-like CSV (required) |
| `--out-prefix PREFIX` | Prefix for FlyGuide native files (default: base name of `--out`) |
| `--write-flyguide-files` | Also write `_species_search.txt` and `_species_kingdom.tsv` |

### API and caching

| Flag | Default | Description |
|---|---|---|
| `--cache-dir DIR` | `.gbif_cache` | Where to cache API responses |
| `--no-cache` | off | Disable caching; always fetch live |
| `--force-refresh` | off | Ignore cached responses and re-fetch |
| `--sleep SEC` | `0.3` | Seconds between API calls |

### Other

| Flag | Description |
|---|---|
| `--list-taxon-groups` | Print taxon group presets and exit |
| `--list-regions` | Print geographic presets and exit |
| `--list-basis` | Print basis-of-record presets and exit |
| `--quiet` | Suppress progress output |
| `--version` | Print version and exit |

---

## Taxon group presets

These common names map to GBIF backbone taxon keys.  Pass the name to `--taxon`
or run `--list-taxon-groups` to see the full table with GBIF keys.

| Preset | Group | Kingdom |
|---|---|---|
| `animals` / `animalia` | Animalia | Animalia |
| `plants` / `plantae` | Plantae | Plantae |
| `fungi` | Fungi | Fungi |
| `mammals` / `mammalia` | Mammalia | Animalia / Chordata |
| `birds` / `aves` | Aves | Animalia / Chordata |
| `reptiles` / `reptilia` | Reptilia | Animalia / Chordata |
| `amphibians` / `amphibia` | Amphibia | Animalia / Chordata |
| `fish` / `actinopterygii` | Actinopterygii | Animalia / Chordata |
| `sharks` / `chondrichthyes` | Chondrichthyes | Animalia / Chordata |
| `insects` / `insecta` | Insecta | Animalia / Arthropoda |
| `arachnids` / `arachnida` | Arachnida | Animalia / Arthropoda |
| `molluscs` / `mollusca` | Mollusca | Animalia / Mollusca |
| `crustaceans` | Malacostraca | Animalia / Arthropoda |
| `whales` / `cetacea` | Cetacea | Animalia / Chordata |
| `vascular_plants` / `tracheophyta` | Tracheophyta | Plantae |
| `vertebrates` / `vertebrata` | Vertebrata | Animalia / Chordata |
| `algae` / `chromista` | Chromista | Chromista |

You can also pass any scientific name (`--taxon Bovidae`, `--taxon Nymphaeales`) or
a numeric GBIF key (`--taxon-key 359`).  The tool will resolve it via the GBIF
species API.

---

## Geographic presets

Run `--list-regions` for the full list.  A selection:

### Land and continental regions

| Preset | Area |
|---|---|
| `global` | Entire world |
| `northern_hemisphere` | Latitudes 0–90°N |
| `holarctic` | >23.5°N (north of the Tropic of Cancer) |
| `north_america` | North America |
| `south_america` | South America |
| `americas` | North + South America |
| `europe` | Europe |
| `eurasia` | Europe + Asia |
| `asia` | Asia |
| `east_asia` | East Asia |
| `southeast_asia` | Southeast Asia |
| `africa` | Africa |
| `sub_saharan_africa` | Africa south of the Sahara |
| `oceania` | Australia + Pacific islands |
| `australia` | Australia only |
| `greenland` | Greenland |
| `antarctica` | Antarctica |
| `beringia` | Beringia (straddles the antimeridian) |

### Ocean basin regions

Records for marine taxa in international waters carry **no country code** in GBIF.
Omitting the country filter and using a bounding box is the correct approach —
the ocean basin presets do exactly this.

| Preset | Area |
|---|---|
| `north_atlantic` | North Atlantic Ocean |
| `south_atlantic` | South Atlantic Ocean |
| `atlantic` | Atlantic Ocean (both hemispheres) |
| `north_pacific` | North Pacific (western) |
| `north_pacific_east` | North Pacific (eastern) |
| `south_pacific` | South Pacific (western) |
| `south_pacific_east` | South Pacific (eastern) |
| `indian_ocean` | Indian Ocean |
| `southern_ocean` | Southern Ocean (south of 60°S) |
| `arctic_ocean` | Arctic Ocean (north of 70°N) |
| `mediterranean` | Mediterranean Sea |
| `caribbean` | Caribbean Sea |
| `coral_triangle` | Coral Triangle |

### Short aliases

`nh` → `northern_hemisphere`, `na` → `north_america`, `sa` → `south_america`,
`au` → `australia`, `nz` → `new_zealand`, `n_atlantic` → `north_atlantic`,
`ind_ocean` → `indian_ocean`.

### Country and continent filters

```bash
# Single country
python3 gbif_query.py --taxon birds --country CA --out canada_birds.csv

# Multiple countries (OR'd)
python3 gbif_query.py --taxon amphibians --country MX,GT,BZ,HN --out mesoamerica_amphibians.csv

# Continent (OR'd if multiple)
python3 gbif_query.py --taxon birds --continent north_america,europe,asia --out holarctic_birds.csv

# Convenience group — expands automatically to the member continents
python3 gbif_query.py --taxon mammals --continent eurasia --out eurasia_mammals.csv
```

Continent groups: `northern_hemisphere` → NA+EU+AS; `holarctic` → NA+EU+AS;
`americas` → NA+SA; `eurasia` → EU+AS; `old_world` → EU+AS+AF.

---

## Basis-of-record presets

| Preset | GBIF record types included |
|---|---|
| `human` (default) | HUMAN_OBSERVATION, PRESERVED_SPECIMEN, OBSERVATION, LITERATURE, MATERIAL_CITATION |
| `observed` | HUMAN_OBSERVATION, OBSERVATION |
| `specimen` | PRESERVED_SPECIMEN (museum collections) |
| `machine` | MACHINE_OBSERVATION (camera traps, acoustics) |
| `fossil` | FOSSIL_SPECIMEN |
| `living` | LIVING_SPECIMEN (zoos, botanical gardens) |
| `extant` | All non-fossil record types |
| `any` | No filter — includes fossils and all record types |

For most reference database work, `human` (the default) is appropriate.  Use
`any` if you want to include all GBIF records regardless of type.

---

## Combining with palaeo sources

`gbif_query.py` outputs the same GBIF-like CSV format as the palaeo exporters
(`neotoma_extinct_to_gbif.py`, `pbdb_to_gbif.py`, `now_to_gbif.py`).  You can
merge them with `flyguide_merge_palaeo_sources.py`:

```bash
# Generate extant species list from GBIF
python3 gbif_query.py \
  --taxon mammals \
  --region northern_hemisphere \
  --basis human \
  --out NH_mammals_extant_gbif.csv

# Generate extinct species list from Neotoma
python3 neotoma_extinct_to_gbif.py \
  --region northern_hemisphere \
  --period quaternary \
  --organisms animals \
  --status extinct \
  --out NH_mammals_neotoma_gbif.csv

# Merge into one deduplicated list
python3 flyguide_merge_palaeo_sources.py \
  --inputs NH_mammals_extant_gbif.csv NH_mammals_neotoma_gbif.csv \
  --out NH_mammals_merged_gbif.csv \
  --write-flyguide-files

# Run FlyGuide
./flyguide.sh \
  NH_mammals_merged_gbif.csv \
  NH_mammals_refs \
  your@email.org \
  YOUR_NCBI_API_KEY
```

---

## Caching

All GBIF API responses are cached to `.gbif_cache/` by default.  The cache is
keyed by the full request URL + parameters, so re-running the exact same query
is instant.  Name resolution (one API call per species key) benefits most from
the cache — a 243-species query that takes ~4 minutes the first time takes
0.17 s on the second run.

The `.gbif_cache/` directory is git-ignored and safe to delete if you want a
fresh start.

```bash
# Disable caching
python3 gbif_query.py --taxon birds --region europe --out europe_birds.csv --no-cache

# Force a fresh fetch (ignore existing cache)
python3 gbif_query.py --taxon birds --region europe --out europe_birds.csv --force-refresh
```

---

## Tutorial

Work through these examples in order to get a feel for the tool.  Each builds
on the one before.

### Step 1 — Explore what's available

```bash
python3 gbif_query.py --list-taxon-groups
python3 gbif_query.py --list-regions
python3 gbif_query.py --list-basis
```

### Step 2 — A small test query (fast)

Scottish freshwater fish — a small, well-documented group:

```bash
python3 gbif_query.py \
  --taxon fish \
  --country GB \
  --basis human \
  --out scotland_fish_gbif.csv
```

Expected: 30–60 species, completes in under 30 seconds.  Look at the CSV:

```bash
column -t -s, scotland_fish_gbif.csv | head -20
```

### Step 3 — Write FlyGuide files and inspect them

```bash
python3 gbif_query.py \
  --taxon fish \
  --country GB \
  --basis human \
  --out scotland_fish_gbif.csv \
  --write-flyguide-files \
  --out-prefix scotland_fish

wc -l scotland_fish_species_search.txt
head scotland_fish_species_search.txt
head scotland_fish_species_kingdom.tsv
```

These two files are what `NCBI-NT_Downloader.pl` and `split_fasta_by_kingdom_organelle_simple.pl`
consume directly.

### Step 4 — A continental-scale query

Northern hemisphere birds (human observations only):

```bash
python3 gbif_query.py \
  --taxon birds \
  --region northern_hemisphere \
  --basis human \
  --out NH_birds_gbif.csv \
  --write-flyguide-files \
  --out-prefix NH_birds
```

This may return 5 000–8 000 species and take 2–5 minutes on first run
(dominated by name resolution).  Subsequent runs use the cache and finish in
seconds.

### Step 5 — Marine taxa outside national boundaries

Sharks and rays in the North Atlantic, including international waters:

```bash
python3 gbif_query.py \
  --taxon sharks \
  --region north_atlantic \
  --basis any \
  --out north_atlantic_sharks_gbif.csv
```

The `north_atlantic` region preset applies a bounding box with no country filter,
so records from the high seas are included.

### Step 6 — Multi-continent query

Holarctic mammals (North America + Europe + Asia) in one command:

```bash
python3 gbif_query.py \
  --taxon mammals \
  --continent northern_hemisphere \
  --basis human \
  --out holarctic_mammals_gbif.csv \
  --write-flyguide-files
```

`northern_hemisphere` expands automatically to
`NORTH_AMERICA,EUROPE,ASIA` — three continent parameters OR'd together.

### Step 7 — Custom taxon by scientific name

Bovids (cattle, bison, antelopes, sheep, goats) in Africa:

```bash
python3 gbif_query.py \
  --taxon Bovidae \
  --region africa \
  --basis human \
  --out africa_bovidae_gbif.csv
```

### Step 8 — Filter by occurrence count

Some GBIF species keys belong to taxa with only a handful of records —
possibly misidentifications or data artefacts.  Drop them:

```bash
python3 gbif_query.py \
  --taxon mammals \
  --region north_america \
  --basis human \
  --min-occurrences 5 \
  --out NA_mammals_gbif.csv
```

### Step 9 — Pipe into FlyGuide

```bash
python3 gbif_query.py \
  --taxon mammals \
  --region north_america \
  --basis human \
  --out NA_mammals_gbif.csv \
  --write-flyguide-files \
  --out-prefix NA_mammals

./flyguide.sh \
  NA_mammals_gbif.csv \
  NA_mammals_refs \
  your@email.org \
  YOUR_NCBI_API_KEY
```

### Step 10 — Merge with palaeo sources

Combine extant + extinct Quaternary North American mammals into one list:

```bash
# Extant species from GBIF
python3 gbif_query.py \
  --taxon mammals \
  --region north_america \
  --basis human \
  --out NA_mammals_extant_gbif.csv

# Extinct species from Neotoma
python3 neotoma_extinct_to_gbif.py \
  --region north_america \
  --period quaternary \
  --organisms animals \
  --status extinct \
  --out NA_mammals_extinct_neotoma_gbif.csv

# Merge
python3 flyguide_merge_palaeo_sources.py \
  --inputs NA_mammals_extant_gbif.csv NA_mammals_extinct_neotoma_gbif.csv \
  --out NA_mammals_complete_gbif.csv \
  --write-flyguide-files \
  --out-prefix NA_mammals_complete

# Download references
./flyguide.sh \
  NA_mammals_complete_gbif.csv \
  NA_mammals_complete_refs \
  your@email.org \
  YOUR_NCBI_API_KEY
```

---

## Troubleshooting

**"Could not resolve taxon"** — the name is not in the GBIF backbone.  Try
`--taxon-key` with a known GBIF key, or check spelling.

**Very large queries time out** — increase `--sleep` to reduce request frequency,
or run with `--force-refresh` if a cached partial result is suspected.

**Species count is lower than expected** — GBIF only returns species with at
least one occurrence matching all filters.  Loosening `--basis` to `any` or
removing `--require-coords` often recovers more records.

**Cache seems stale** — run with `--force-refresh` to re-fetch all responses.
Or delete `.gbif_cache/` entirely and start fresh.

---

## Dependencies

- Python 3.7+, standard library only
- Internet access to `api.gbif.org`

No GBIF account or API key required.
