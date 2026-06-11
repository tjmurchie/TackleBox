# FlyGuide Neotoma exporter

_Part of TackleBox FlyGuide v1.3.0_

`neotoma_extinct_to_gbif.py` queries the Neotoma Paleoecology Database API and
produces a GBIF-like CSV that feeds directly into FlyGuide:

```text
Neotoma occurrences + taxa table
        ↓
neotoma_extinct_to_gbif.py
        ↓
GBIF-like CSV with species/genus/kingdom/phylum columns
        ↓
gbif_prep_from_csv.py → NCBI-NT_Downloader.pl
```

FlyGuide does not actually require GBIF itself — it requires a GBIF-shaped table
with at least `species`, `genus`, `kingdom`, `phylum` columns. This exporter
creates that table from Neotoma and also writes optional FlyGuide-native files:

```text
OUTPREFIX_species_search.txt
OUTPREFIX_species_kingdom.tsv
```

## Quick start

The fastest way to run Neotoma together with PBDB is the unified wrapper:

```bash
python3 flyguide_palaeo_sources.py \
  --region northern_hemisphere \
  --period quaternary \
  --organisms animals \
  --out NH_quaternary_extinct_animals \
  --write-flyguide-files
```

To run the Neotoma exporter alone:

```bash
python3 neotoma_extinct_to_gbif.py \
  --region northern_hemisphere \
  --period quaternary \
  --organisms animals \
  --status extinct \
  --out NH_quaternary_extinct_animals_neotoma_gbif.csv \
  --write-flyguide-files
```

Then run FlyGuide directly:

```bash
./flyguide.sh \
  --ncbi-name-mode species \
  NH_quaternary_extinct_animals_neotoma_gbif.csv \
  NH_quaternary_extinct_animals_refs \
  you@example.org \
  NCBI_API_KEY
```

Or use the Neotoma+FlyGuide convenience wrapper:

```bash
./flyguide_neotoma.sh \
  --region northern_hemisphere \
  --period quaternary \
  --organisms animals \
  --status extinct \
  --ncbi-name-mode binomial \
  -- NH_quaternary_extinct_animals_refs you@example.org NCBI_API_KEY
```

## Files

| File | Purpose |
|---|---|
| `neotoma_extinct_to_gbif.py` | Neotoma exporter. Python standard library only. |
| `flyguide_neotoma.sh` | Convenience wrapper: Neotoma export followed by FlyGuide. |
| `flyguide_palaeo_sources.py` | Unified palaeo wrapper: Neotoma + PBDB + optional NOW in one command. |
| `flyguide_merge_palaeo_sources.py` | Standalone merger for combining multiple palaeo source CSVs. |
| `flyguide.sh` | FlyGuide main pipeline script. |
| `NCBI-NT_Downloader.pl` | NCBI nucleotide downloader (supports `--name-mode species|trinomial|as-is`). |
| `tests/` | Offline unit tests and fixtures. |

## Region examples

Built-in regions are coarse bounding-box presets designed for palaeo reference-panel
building, not exact political or biogeographic boundaries.

```bash
python3 neotoma_extinct_to_gbif.py --list-regions
```

Useful examples:

```bash
--region northern_hemisphere
--region holarctic
--region north_america
--region eurasia
--region tanzania
--region east_africa
--region beringia
```

Custom region options:

```bash
--bbox minLon,minLat,maxLon,maxLat
--geojson my_polygon.geojson
```

Examples:

```bash
# Tanzania-ish bbox
--bbox 29,-12.5,41.5,-0.5

# More exact polygon
--geojson holarctic_custom.geojson
```

## Period examples

```bash
python3 neotoma_extinct_to_gbif.py --list-periods
```

Built-in presets:

```bash
--period quaternary
--period pleistocene
--period holocene
--period late_pleistocene
--period late_quaternary
--period lgm
--period pliocene
--period all
```

Override ages directly (Neotoma uses years before present):

```bash
--age-young 0 --age-old 50000
```

## Organism filters

```bash
--organisms animals
--organisms plants
--organisms fungi
--organisms diatoms
--organisms all
```

Comma-separated values are allowed:

```bash
--organisms animals,plants,fungi
```

Neotoma uses coarse organism groups (`animals`, `plants`, `fungi`, `diatoms`,
`protists`). Fine-grained clade names such as `mammals` or `vertebrates` are not
supported by Neotoma directly — the unified wrapper `flyguide_palaeo_sources.py`
handles this by translating them to `animals` for Neotoma while passing them
unmodified to PBDB, which supports full taxonomic filtering.

## Extinct/extant/all

Default:

```bash
--status extinct
```

Other options:

```bash
--status extant
--status all
```

The filter uses Neotoma's `taxa.extinct` field. If a taxon has no matching
taxa-table row and `--status extinct` is used, it is not retained.

Note: the unified `flyguide_palaeo_sources.py` wrapper defaults to `--status all`
so that both Neotoma and PBDB apply the same filter consistently.

## Name handling and subspecies

For FlyGuide/NCBI the main goal is retrievable search buckets, not strict taxonomic
realism. The default collapses trinomial/subspecies names to binomials:

```bash
--ncbi-name-mode binomial
```

Example:

```text
Bison bison antiquus  →  Bison bison
```

This is usually best for old splitter categories because many will not exist in NCBI
as valid subspecies. The original Neotoma name is retained in metadata columns.

Other modes:

```bash
--ncbi-name-mode trinomial   # keep Genus species subspecies when present
--ncbi-name-mode as-is       # cleaned original string
```

### Name mode naming: Python vs shell

There is a naming difference to be aware of when mixing the exporter with the
shell script directly:

| Context | Mode name | Meaning |
|---|---|---|
| `neotoma_extinct_to_gbif.py` | `binomial` | Collapse to Genus species |
| `flyguide.sh` | `species` | Same behaviour (legacy name) |
| `flyguide_neotoma.sh` | `binomial` | Translates internally to `species` for `flyguide.sh` |

The unified `flyguide_palaeo_sources.py` wrapper uses `binomial`/`trinomial`/`as-is`
throughout and handles the translation automatically.

The exporter retains genus-level buckets by default for names like `Canis spp.` or
`Pinus undiff.` because genus-level NCBI searches are often useful for broad
reference panels. To drop them:

```bash
--drop-genus-only
```

Morphotype/type names such as `Ambrosia-type` are converted to genus buckets by
default. To drop them:

```bash
--drop-morphotypes
```

Slash names are rejected by default. To split simple cases:

```bash
--split-slash-names
```

## Large searches

Broad queries such as:

```bash
--region northern_hemisphere --period quaternary --organisms all --status all
```

may be very large. Recommended workflow:

1. Run with a cap first:

```bash
python3 neotoma_extinct_to_gbif.py \
  --region northern_hemisphere \
  --period quaternary \
  --organisms animals \
  --status extinct \
  --max-occurrences 10000 \
  --out test_neotoma.csv
```

2. Inspect `test_neotoma.csv`, `test_neotoma.rejected.tsv`, and
   `test_neotoma.summary.json`.
3. Remove `--max-occurrences` for the full run.

The exporter caches API pages in `.neotoma_cache/` by default:

```bash
--force-refresh
```

ignores cached responses.

## Output columns

FlyGuide compatibility columns:

```text
species
genus
kingdom
phylum
```

Neotoma metadata columns:

```text
neotoma_taxonids
neotoma_taxonnames
neotoma_taxagroupids
neotoma_occurrence_count
neotoma_site_count
neotoma_dataset_count
neotoma_min_age_young
neotoma_max_age_old
neotoma_databases
neotoma_datasettypes
name_cleaning_status
ncbiSearchName
```

## Output files

For `--out example.csv`, the exporter writes:

```text
example.csv
example.rejected.tsv
example.summary.json
```

If `--write-flyguide-files` is used:

```text
example_species_search.txt
example_species_kingdom.tsv
```

The rejected file documents which Neotoma names were too messy or ambiguous to
convert into NCBI search names — worth reviewing before running large downloads.

## Testing

Offline tests:

```bash
python3 -m pytest -q tests
```

Manual offline smoke test:

```bash
python3 neotoma_extinct_to_gbif.py \
  --offline-taxa-json tests/data/neotoma_taxa_fixture.json \
  --offline-occurrences-json tests/data/neotoma_occurrences_fixture.json \
  --region northern_hemisphere \
  --period quaternary \
  --organisms animals \
  --status extinct \
  --out smoke.csv \
  --write-flyguide-files
```

Live API smoke test (requires internet):

```bash
python3 neotoma_extinct_to_gbif.py \
  --region tanzania \
  --period late_quaternary \
  --organisms animals \
  --status extinct \
  --max-occurrences 2000 \
  --out tanzania_late_quaternary_extinct_animals.csv
```

## Caveats

- Region presets are coarse. Use `--geojson` for publication-grade boundaries.
- Neotoma names often preserve original identifications and palaeo conventions
  (`cf.`, `aff.`, `sp.`, `undiff.`, `-type`, slash taxa). The exporter is
  intentionally transparent: it preserves original names and writes rejected
  names to a TSV rather than silently discarding them.
- This is not a definitive extinct-species checklist. It is an occurrence-derived
  list of Neotoma taxa within the requested spatial/temporal filters.
- NCBI retrieval is still dependent on what NCBI has and how records are named.
  Many palaeo names will have no NCBI sequence data — that is expected.
