# FlyGuide Tutorial: Building Regional Reference Databases

_TackleBox FlyGuide v1.3.0_

This tutorial walks through practical workflows from species list to downloaded
NCBI sequences.  It covers three routes:

1. **Extant species via GBIF** — query the GBIF occurrence API from the command line.
2. **Palaeo/extinct species** — pull from Neotoma, PBDB, and NOW.
3. **Combined** — merge extant + extinct into a single NCBI-ready list.

All scripts are standard-library Python 3.7+.  No pip installs required.

---

## Part 1: GBIF live queries (`gbif_query.py`)

### What it does

`gbif_query.py` queries GBIF's occurrence REST API for a species list.  You
give it a taxon group and a geographic area; it returns a GBIF-like CSV that
feeds straight into the FlyGuide pipeline.  No GBIF account or API key needed.
Responses are cached to disk so re-runs are instant.

---

### 1.1 — Explore what's available

Before running any query, check the built-in presets:

```bash
python3 gbif_query.py --list-taxon-groups
python3 gbif_query.py --list-regions
python3 gbif_query.py --list-basis
```

This shows you every preset you can pass to `--taxon`, `--region`, and `--basis`
without having to look anything up.

---

### 1.2 — A small test query

Scottish freshwater fish — a small, well-documented group, good for a quick test:

```bash
python3 gbif_query.py \
  --taxon fish \
  --country GB \
  --basis human \
  --out scotland_fish_gbif.csv
```

Expected output: 30–60 species, runs in under 30 seconds.  The first thing to
check is the CSV:

```bash
column -t -s, scotland_fish_gbif.csv | head -20
```

And look at any rejected names (entries the name-cleaner couldn't reduce to a
valid binomial):

```bash
cat scotland_fish_gbif.csv.rejected.tsv
```

Typical rejects are hybrids (`Salmo trutta x Salmo salar`), cultivar names,
or entries with only a genus.  These are written to `.rejected.tsv` alongside
the main CSV so you can audit them.

---

### 1.3 — Write FlyGuide-ready files

Add `--write-flyguide-files` to produce the two files that
`NCBI-NT_Downloader.pl` and `split_fasta_by_kingdom_organelle_simple.pl` consume:

```bash
python3 gbif_query.py \
  --taxon fish \
  --country GB \
  --basis human \
  --out scotland_fish_gbif.csv \
  --write-flyguide-files \
  --out-prefix scotland_fish
```

This writes:

- `scotland_fish_species_search.txt` — one name per line (the NCBI query list)
- `scotland_fish_species_kingdom.tsv` — name TAB kingdom TAB phylum

Inspect them:

```bash
wc -l scotland_fish_species_search.txt
head scotland_fish_species_search.txt
head scotland_fish_species_kingdom.tsv
```

---

### 1.4 — Run the full FlyGuide pipeline

Once you have the FlyGuide files, pass them to `flyguide.sh` to download NCBI
sequences:

```bash
./flyguide.sh \
  scotland_fish_gbif.csv \
  scotland_fish_refs \
  you@example.org \
  YOUR_NCBI_API_KEY
```

This produces:
- `scotland_fish_refs.fasta` — all downloaded sequences
- `scotland_fish_refs.Animal-Mito.fasta`, `scotland_fish_refs.Animal-NucMark.fasta`, etc.

---

### 1.5 — Larger queries

#### Terrestrial northern hemisphere mammals

```bash
python3 gbif_query.py \
  --taxon mammals \
  --region northern_hemisphere \
  --basis human \
  --out NH_mammals_gbif.csv \
  --write-flyguide-files \
  --out-prefix NH_mammals
```

Expect 4 000–6 000 species.  First run takes 2–5 minutes (name resolution);
subsequent runs use the cache and finish in seconds.

#### All vascular plants from Canada

```bash
python3 gbif_query.py \
  --taxon vascular_plants \
  --country CA \
  --out canada_plants_gbif.csv \
  --write-flyguide-files \
  --out-prefix canada_plants
```

#### Holarctic birds (multi-continent, OR'd)

```bash
python3 gbif_query.py \
  --taxon birds \
  --continent northern_hemisphere \
  --basis human \
  --out holarctic_birds_gbif.csv \
  --write-flyguide-files \
  --out-prefix holarctic_birds
```

`northern_hemisphere` expands to `NORTH_AMERICA,EUROPE,ASIA` automatically.

#### Custom taxon by scientific name

Any GBIF taxon name works:

```bash
python3 gbif_query.py \
  --taxon Bovidae \
  --region africa \
  --out africa_bovidae_gbif.csv
```

---

### 1.6 — Marine and high-seas taxa

Records in international waters carry no country code in GBIF.  Ocean basin
region presets apply a bounding box without a country filter, which is the
correct approach:

```bash
# Sharks and rays in the North Atlantic (open ocean included)
python3 gbif_query.py \
  --taxon chondrichthyes \
  --region north_atlantic \
  --basis any \
  --out north_atlantic_sharks_gbif.csv

# Cetacea globally — covers high seas
python3 gbif_query.py \
  --taxon cetacea \
  --region global \
  --out global_cetacea_gbif.csv

# Indian Ocean bony fish
python3 gbif_query.py \
  --taxon fish \
  --region indian_ocean \
  --out indian_ocean_fish_gbif.csv
```

---

### 1.7 — Filtering and quality control

```bash
# Only records with geolocated, quality-checked coordinates
python3 gbif_query.py \
  --taxon amphibians \
  --region north_america \
  --require-coords \
  --out NA_amphibians_gbif.csv

# Drop species with fewer than 5 GBIF occurrences
python3 gbif_query.py \
  --taxon insects \
  --country CA \
  --min-occurrences 5 \
  --out canada_insects_gbif.csv

# Restrict to a specific year range
python3 gbif_query.py \
  --taxon mammals \
  --country US \
  --year 2000,2024 \
  --out US_mammals_recent_gbif.csv

# Museum specimens only (no field observations)
python3 gbif_query.py \
  --taxon reptiles \
  --continent africa \
  --basis specimen \
  --out africa_reptiles_specimens_gbif.csv
```

---

## Part 2: Palaeo sources

Three databases cover complementary niches for extinct/Quaternary taxa:

| Source | Best for | How to access |
|---|---|---|
| **Neotoma** | Quaternary palaeoecology; vertebrates, pollen; strong North America | Live REST API |
| **PBDB** | Global fossil occurrences across all geologic time | Live REST API |
| **NOW** | Cenozoic fossil mammals; strong Eurasia/Africa | Browser export → `now_to_gbif.py` |

---

### 2.1 — Neotoma (`neotoma_extinct_to_gbif.py`)

#### Explore presets

```bash
python3 neotoma_extinct_to_gbif.py --list-regions
python3 neotoma_extinct_to_gbif.py --list-periods
python3 neotoma_extinct_to_gbif.py --list-organisms
```

#### Quaternary extinct animals from North America

```bash
python3 neotoma_extinct_to_gbif.py \
  --region north_america \
  --period quaternary \
  --organisms animals \
  --status extinct \
  --out NA_quaternary_extinct_neotoma_gbif.csv \
  --write-flyguide-files \
  --out-prefix NA_quaternary_extinct
```

What `--status extinct` means: the exporter keeps only Neotoma taxa flagged as
extinct or that have no extant GBIF match.  Use `--status all` to include both
extant and extinct.

#### Holarctic Pleistocene animals (broader geographic query)

```bash
python3 neotoma_extinct_to_gbif.py \
  --region holarctic \
  --period pleistocene \
  --organisms animals \
  --status extinct \
  --out holarctic_pleistocene_extinct_neotoma_gbif.csv
```

#### North American Holocene plants

```bash
python3 neotoma_extinct_to_gbif.py \
  --region north_america \
  --period holocene \
  --organisms plants \
  --status all \
  --out NA_holocene_plants_neotoma_gbif.csv
```

#### Name mode

By default (`--name-mode binomial`) all names are collapsed to two tokens for
NCBI retrievability.  Many extinct subspecies are not in NCBI as trinomials:

```bash
# Keep subspecies where present (trinomial mode)
python3 neotoma_extinct_to_gbif.py \
  --region north_america \
  --period pleistocene \
  --organisms animals \
  --name-mode trinomial \
  --out NA_pleistocene_trinomial_neotoma_gbif.csv
```

---

### 2.2 — PBDB (`pbdb_to_gbif.py`)

#### Explore presets

```bash
python3 pbdb_to_gbif.py --list-regions
python3 pbdb_to_gbif.py --list-periods
python3 pbdb_to_gbif.py --list-organisms
```

#### Northern hemisphere Quaternary animals

```bash
python3 pbdb_to_gbif.py \
  --region northern_hemisphere \
  --period quaternary \
  --organisms animals \
  --out NH_quaternary_animals_pbdb_gbif.csv \
  --write-flyguide-files
```

PBDB queries the [Paleobiology Database REST API](https://paleobiodb.org/data1.2)
and caches to `.pbdb_cache/`.

#### Eurasian Pleistocene mammals

```bash
python3 pbdb_to_gbif.py \
  --region eurasia \
  --period pleistocene \
  --organisms mammals \
  --out eurasia_pleistocene_mammals_pbdb_gbif.csv \
  --write-flyguide-files
```

#### Small smoke test before a large query

```bash
python3 pbdb_to_gbif.py \
  --region africa \
  --period quaternary \
  --organisms mammals \
  --limit 200 \
  --out africa_test_pbdb_gbif.csv \
  --verbose
```

---

### 2.3 — NOW (`now_to_gbif.py`)

NOW does not expose a stable public REST API.  Export a table from the
[NOW browser interface](https://nowdatabase.luomus.fi/), then pass it to
`now_to_gbif.py`.

```bash
# Get export instructions
python3 now_to_gbif.py --explain-now-export

# Import and filter a NOW export
python3 now_to_gbif.py \
  --input NOW_export.csv \
  --region eurasia \
  --period quaternary \
  --out eurasia_quaternary_now_gbif.csv \
  --write-flyguide-files

# Tanzania mammals from a NOW export
python3 now_to_gbif.py \
  --input NOW_export.csv \
  --region africa \
  --period cenozoic \
  --out africa_cenozoic_mammals_now_gbif.csv
```

---

### 2.4 — One-command all-sources wrapper (`flyguide_palaeo_sources.py`)

For most workflows, this is the recommended entry point.  It runs Neotoma +
PBDB (and optionally NOW) and merges the results automatically:

```bash
# Explore presets
python3 flyguide_palaeo_sources.py --list-regions
python3 flyguide_palaeo_sources.py --list-periods
python3 flyguide_palaeo_sources.py --list-organisms

# Northern hemisphere Quaternary extinct animals (default settings)
python3 flyguide_palaeo_sources.py \
  --out-prefix NH_quat_palaeo

# Eurasian Pleistocene mammals
python3 flyguide_palaeo_sources.py \
  --region eurasia \
  --period pleistocene \
  --organisms mammals \
  --out-prefix eurasia_pleis_mammals

# Add a pre-exported NOW table to the above
python3 flyguide_palaeo_sources.py \
  --region eurasia \
  --period pleistocene \
  --organisms mammals \
  --now-input NOW_export.csv \
  --out-prefix eurasia_pleis_mammals_with_now

# Chain directly into NCBI download
python3 flyguide_palaeo_sources.py \
  --region north_america \
  --period quaternary \
  --organisms animals \
  --out-prefix NA_quat_palaeo \
  --run-flyguide you@example.org YOUR_NCBI_API_KEY

# Smoke test (limit PBDB records, fast)
python3 flyguide_palaeo_sources.py \
  --region africa \
  --period quaternary \
  --organisms mammals \
  --limit 200 \
  --out-prefix africa_quat_test
```

---

### 2.5 — Merging sources manually (`flyguide_merge_palaeo_sources.py`)

When you run the exporters separately and want to combine them:

```bash
python3 flyguide_merge_palaeo_sources.py \
  --inputs \
    NA_quaternary_extinct_neotoma_gbif.csv \
    NH_quaternary_animals_pbdb_gbif.csv \
    eurasia_quaternary_now_gbif.csv \
  --collapse binomial \
  --out NH_quaternary_palaeo_merged_gbif.csv \
  --write-flyguide-files \
  --out-prefix NH_quaternary_palaeo_merged
```

The `--collapse` option controls deduplication level:

| Mode | Effect | When to use |
|---|---|---|
| `binomial` (default) | Collapse to genus + species | Best NCBI retrievability |
| `trinomial` | Keep genus + species + subspecies | When subspecies are in NCBI |
| `genus` | Collapse to genus only | Very coarse; not recommended for downloads |

Source priority (which database's metadata wins on conflict):

```bash
python3 flyguide_merge_palaeo_sources.py \
  --inputs pbdb.csv neotoma.csv \
  --source-priority PBDB,Neotoma,NOW,GBIF \
  --out merged.csv
```

---

## Part 3: Combined workflows (extant + extinct)

This is the most powerful use of FlyGuide for ancient DNA and sedaDNA work: a
single NCBI-ready list that includes both extant species observed in the region
today and extinct/palaeo taxa from fossil databases.

### 3.1 — Northern hemisphere mammals (extant + Quaternary extinct)

```bash
# Step 1: extant mammals from GBIF
python3 gbif_query.py \
  --taxon mammals \
  --region northern_hemisphere \
  --basis human \
  --out NH_mammals_extant_gbif.csv

# Step 2: Quaternary extinct mammals from Neotoma + PBDB
python3 flyguide_palaeo_sources.py \
  --region northern_hemisphere \
  --period quaternary \
  --organisms mammals \
  --status extinct \
  --out-prefix NH_mammals_extinct_palaeo

# Step 3: merge
python3 flyguide_merge_palaeo_sources.py \
  --inputs \
    NH_mammals_extant_gbif.csv \
    NH_mammals_extinct_palaeo_merged_gbif.csv \
  --out NH_mammals_complete_gbif.csv \
  --write-flyguide-files \
  --out-prefix NH_mammals_complete

# Step 4: download NCBI sequences
./flyguide.sh \
  NH_mammals_complete_gbif.csv \
  NH_mammals_complete_refs \
  you@example.org \
  YOUR_NCBI_API_KEY
```

---

### 3.2 — Eurasian sedaDNA panel (plants + vertebrates + palaeo)

A practical panel for Quaternary sedaDNA:

```bash
# Extant plants (GBIF)
python3 gbif_query.py \
  --taxon vascular_plants \
  --continent eurasia \
  --basis human \
  --out eurasia_plants_gbif.csv

# Extant vertebrates (GBIF)
python3 gbif_query.py \
  --taxon vertebrates \
  --continent eurasia \
  --basis human \
  --out eurasia_vertebrates_gbif.csv

# Pleistocene/Quaternary palaeo fauna (Neotoma + PBDB + NOW)
python3 flyguide_palaeo_sources.py \
  --region eurasia \
  --period quaternary \
  --organisms animals \
  --status extinct \
  --now-input NOW_export.csv \
  --out-prefix eurasia_quat_extinct_palaeo

# Merge everything
python3 flyguide_merge_palaeo_sources.py \
  --inputs \
    eurasia_plants_gbif.csv \
    eurasia_vertebrates_gbif.csv \
    eurasia_quat_extinct_palaeo_merged_gbif.csv \
  --out eurasia_sedaDNA_panel_gbif.csv \
  --write-flyguide-files \
  --out-prefix eurasia_sedaDNA_panel

# Download
./flyguide.sh \
  eurasia_sedaDNA_panel_gbif.csv \
  eurasia_sedaDNA_panel_refs \
  you@example.org \
  YOUR_NCBI_API_KEY
```

---

### 3.3 — Arctic/subarctic panel (good for permafrost eDNA)

```bash
# Extant holarctic fauna
python3 gbif_query.py \
  --taxon vertebrates \
  --region holarctic \
  --basis human \
  --out holarctic_vertebrates_gbif.csv

# Pleistocene megafauna
python3 flyguide_palaeo_sources.py \
  --region holarctic \
  --period pleistocene \
  --organisms animals \
  --status extinct \
  --out-prefix holarctic_pleistocene_extinct

# Beringian species specifically
python3 gbif_query.py \
  --taxon mammals \
  --region beringia \
  --basis any \
  --out beringia_mammals_gbif.csv

# Merge
python3 flyguide_merge_palaeo_sources.py \
  --inputs \
    holarctic_vertebrates_gbif.csv \
    holarctic_pleistocene_extinct_merged_gbif.csv \
    beringia_mammals_gbif.csv \
  --out arctic_panel_gbif.csv \
  --write-flyguide-files \
  --out-prefix arctic_panel
```

---

## Part 4: Convenience wrappers

### Neotoma → FlyGuide in one step (`flyguide_neotoma.sh`)

```bash
# Generate Neotoma CSV and immediately feed it to FlyGuide
./flyguide_neotoma.sh \
  --region north_america \
  --period quaternary \
  --organisms animals \
  --status extinct \
  -- NA_quat_extinct you@example.org YOUR_NCBI_API_KEY

# With GuideCheck enabled
./flyguide_neotoma.sh \
  --region holarctic \
  --period pleistocene \
  --organisms animals \
  --status extinct \
  -- holarctic_pleis you@example.org YOUR_NCBI_API_KEY --guidecheck
```

The `--` separates Neotoma options from FlyGuide arguments.

---

## Part 5: Tips and common patterns

### Inspect before downloading

For any large query, check the CSV and species count before launching the NCBI
download — it may run for hours:

```bash
# See how many species you have
wc -l NA_mammals_complete_species_search.txt

# Quick look at the merged list
head -30 NA_mammals_complete_gbif.csv | column -t -s,

# Check kingdom distribution
awk -F, 'NR>1 {print $3}' NA_mammals_complete_gbif.csv | sort | uniq -c | sort -rn
```

### Source provenance

The merged CSV always records where each name came from:

```bash
# See which sources contributed
awk -F, 'NR>1 {print $NF}' NA_mammals_complete_gbif.csv | sort | uniq -c
```

The `source_databases` column shows e.g. `Neotoma;PBDB` when the same name
appeared in both.

### Cache management

```bash
# First run downloads and caches everything
python3 gbif_query.py --taxon birds --region europe --out europe_birds.csv

# Re-run is instant (uses cache)
python3 gbif_query.py --taxon birds --region europe --out europe_birds.csv

# Force a fresh fetch (e.g. after new GBIF data is released)
python3 gbif_query.py --taxon birds --region europe --out europe_birds.csv --force-refresh

# Delete all caches to start fresh
rm -rf .gbif_cache/ .pbdb_cache/ .neotoma_cache/
```

### GuideCheck before downloading

Run GuideCheck on your merged list to see which taxa have NCBI entries before
committing to a long download:

```bash
bash guidecheck.sh \
  -i NA_mammals_complete_species_search.txt \
  -o NA_mammals_ncbi_guidecheck.tsv \
  --api-key YOUR_NCBI_API_KEY

# See how many have no NCBI data
awk -F'\t' '$8=="NONE"' NA_mammals_ncbi_guidecheck.tsv | wc -l
```

---

## Reference

| Script | Purpose | Full docs |
|---|---|---|
| `gbif_query.py` | GBIF live API species list query | `README_GBIF.md` |
| `neotoma_extinct_to_gbif.py` | Neotoma palaeoecology API exporter | `README_NEOTOMA.md` |
| `pbdb_to_gbif.py` | Paleobiology Database API exporter | `README_PBDB_NOW.md` |
| `now_to_gbif.py` | NOW fossil mammal database importer | `README_PBDB_NOW.md` |
| `flyguide_palaeo_sources.py` | One-command all-palaeo-sources wrapper | `README_PBDB_NOW.md` |
| `flyguide_merge_palaeo_sources.py` | Merge any combination of source CSVs | `README_PBDB_NOW.md` |
| `flyguide_neotoma.sh` | Neotoma → FlyGuide convenience wrapper | `README_NEOTOMA.md` |
| `flyguide.sh` | Main FlyGuide pipeline wrapper | `README.md` |
| `NCBI-NT_Downloader.pl` | NCBI nucleotide downloader | `README.md` |
