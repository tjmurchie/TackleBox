# FlyGuide palaeo-source exporters: PBDB and NOW

_Part of TackleBox FlyGuide v1.3.0_

FlyGuide includes two additional palaeo-source exporters for building large fossil
taxon buckets for ancient DNA and sedaDNA reference panel construction:

- `pbdb_to_gbif.py` — live Paleobiology Database API exporter.
- `now_to_gbif.py` — robust importer for exported NOW fossil mammal database tables.
- `flyguide_merge_palaeo_sources.py` — merger for combining Neotoma/PBDB/NOW/GBIF-like outputs.
- `flyguide_palaeo_sources.py` — unified one-command wrapper that queries all sources together.

All outputs are FlyGuide-compatible GBIF-like CSVs containing at least:

```text
species, genus, kingdom, phylum
```

They can also directly write the FlyGuide-native input files:

```text
OUTPREFIX_species_search.txt
OUTPREFIX_species_kingdom.tsv
```

## Quick start: all sources in one command

The easiest way to build a merged palaeo reference bucket is the unified wrapper:

```bash
python3 flyguide_palaeo_sources.py \
  --region northern_hemisphere \
  --period quaternary \
  --organisms animals \
  --out NH_quaternary_palaeo \
  --write-flyguide-files
```

This queries Neotoma and PBDB in parallel, merges the results at binomial level,
and writes the merged CSV plus FlyGuide-native species files. Add `--now-input`
with a NOW export file to include the NOW database as well.

---

## Philosophy

The goal is not to create a definitive taxonomic checklist. The goal is to create a
large, auditable, NCBI-searchable bucket of fossil/palaeo names for reference building.
Many names will not have NCBI records. Many subspecies/splitter names will collapse to
binomials by default. The original source names and IDs are preserved so the collapse
can be audited.

Default name behaviour is intentionally NCBI-friendly:

```text
Bison bison antiquus -> Bison bison
Mammuthus primigenius -> Mammuthus primigenius
Canis sp. -> Canis
```

Use `--ncbi-name-mode trinomial` when you want to keep subspecies/trinomials where
possible.

---

## Source comparison

| Source | Coverage | Access | Strengths |
|---|---|---|---|
| Neotoma | Quaternary, global | Live API | Palaeoecological context; site/dataset metadata |
| PBDB | Phanerozoic, global | Live API | Broad fossil coverage; invertebrates; deep time |
| NOW | Cenozoic mammals, global | Browser export | Strong Eurasian/African Cenozoic mammal coverage |
| GBIF fossil | Recent–Quaternary | Download | Museum specimen data; modern+fossil together |

---

## PBDB exporter

PBDB is the best broad global fossil occurrence source to add after Neotoma. It
supports spatial, temporal, and taxonomic filters through a public REST API.

### Basic examples

Northern Hemisphere Quaternary animals:

```bash
python3 pbdb_to_gbif.py \
  --region northern_hemisphere \
  --period quaternary \
  --organisms animals \
  --out NH_quaternary_animals_pbdb_gbif.csv \
  --write-flyguide-files
```

Holarctic Pleistocene mammals:

```bash
python3 pbdb_to_gbif.py \
  --region holarctic \
  --period pleistocene \
  --organisms mammals \
  --out Holarctic_pleistocene_mammals_pbdb_gbif.csv \
  --write-flyguide-files
```

Eurasian Quaternary vertebrates:

```bash
python3 pbdb_to_gbif.py \
  --region eurasia \
  --period quaternary \
  --organisms vertebrates \
  --out Eurasia_quaternary_vertebrates_pbdb_gbif.csv \
  --write-flyguide-files
```

Tanzania Quaternary mammals, small smoke-test size:

```bash
python3 pbdb_to_gbif.py \
  --region tanzania \
  --period quaternary \
  --organisms mammals \
  --limit 100 \
  --out Tanzania_quaternary_mammals_pbdb_gbif.csv \
  --write-flyguide-files \
  --verbose
```

Custom clade:

```bash
python3 pbdb_to_gbif.py \
  --region africa \
  --period quaternary \
  --base-name Bovidae \
  --out Africa_quaternary_bovidae_pbdb_gbif.csv \
  --write-flyguide-files
```

Custom age range in Ma:

```bash
python3 pbdb_to_gbif.py \
  --region eurasia \
  --age-young-ma 0 \
  --age-old-ma 0.05 \
  --organisms mammals \
  --out Eurasia_0_50ka_mammals_pbdb_gbif.csv
```

### Region presets

```bash
python3 pbdb_to_gbif.py --list-regions
```

Common presets:

- `northern_hemisphere`
- `holarctic`
- `north_america`
- `eurasia`
- `europe`
- `asia`
- `africa`
- `tanzania`
- `global`

Custom bounding box:

```bash
--bbox minLon,minLat,maxLon,maxLat
```

### Period presets

```bash
python3 pbdb_to_gbif.py --list-periods
```

Common presets:

- `quaternary`
- `pleistocene`
- `holocene`
- `late_pleistocene`
- `late_quaternary`
- `lgm`
- `pliocene`
- `all`

### Organism presets

```bash
python3 pbdb_to_gbif.py --list-organisms
```

Useful presets:

- `animals`
- `vertebrates`
- `mammals`
- `birds`
- `fish`
- `reptiles_amphibians`
- `invertebrates`
- `arthropods`
- `molluscs`
- `plants`
- `fungi`
- `all`

Internally these map to PBDB `base_name` taxonomic searches. Override with
`--base-name` for any clade not in the presets.

### PBDB status/extinction filtering

`--status all` is the default and safest. PBDB occurrence records do not always
include a simple extant/extinct flag. The script supports:

```bash
--status extinct
--status extant
--status all
```

but treats this as a best-effort filter only when an extant/status field is actually
returned by the API. For ancient DNA reference-building, `--period quaternary` plus
fossil occurrence filtering is usually more useful than strict extinct-only filtering.

---

## NOW importer

NOW is a Cenozoic fossil mammal database. The public site has a browser UI and export
support, but does not have a clearly documented stable JSON API like PBDB. Therefore
`now_to_gbif.py` is designed as a robust importer for exported CSV/TSV/XLSX/ZIP tables.

### Get export instructions

```bash
python3 now_to_gbif.py --explain-now-export
```

### Single merged export

If the NOW UI gives you a single locality-species table, run:

```bash
python3 now_to_gbif.py \
  --input NOW_export.csv \
  --region eurasia \
  --period quaternary \
  --out Eurasia_quaternary_mammals_now_gbif.csv \
  --write-flyguide-files
```

### Three-table export

If the NOW UI gives you separate species, localities, and locality-species/faunal-list
tables:

```bash
python3 now_to_gbif.py \
  --species-file now_species.csv \
  --localities-file now_localities.csv \
  --locality-species-file now_locality_species.csv \
  --region holarctic \
  --period pleistocene \
  --out Holarctic_pleistocene_mammals_now_gbif.csv \
  --write-flyguide-files
```

### Directory/zip inference

```bash
python3 now_to_gbif.py \
  --source-dir NOW_exports/ \
  --region northern_hemisphere \
  --period quaternary \
  --out NH_quaternary_mammals_now_gbif.csv \
  --write-flyguide-files
```

The script attempts to infer whether files are species/localities/link tables or a
merged table.

### Strict filtering

By default, records lacking coordinates or numeric ages are retained rather than
silently dropped. For stricter geographic/time filtering:

```bash
--require-coords --require-age
```

NOW ages are usually in Ma. If your export uses ka or BP:

```bash
--age-units ka
--age-units bp
```

### Using NOW with the unified wrapper

Pass your NOW export file to `flyguide_palaeo_sources.py` with `--now-input`:

```bash
python3 flyguide_palaeo_sources.py \
  --region holarctic \
  --period pleistocene \
  --organisms mammals \
  --now-input NOW_export.csv \
  --out Holarctic_pleistocene_mammals_all_sources \
  --write-flyguide-files
```

This merges Neotoma + PBDB + NOW in a single step.

---

## Merging sources manually

After running exporters separately, merge with:

```bash
python3 flyguide_merge_palaeo_sources.py \
  --inputs \
    NH_quaternary_animals_neotoma_gbif.csv \
    NH_quaternary_animals_pbdb_gbif.csv \
    NH_quaternary_mammals_now_gbif.csv \
  --collapse binomial \
  --out NH_quaternary_palaeo_merged_gbif.csv \
  --write-flyguide-files
```

Custom source priority (controls which source's kingdom/phylum metadata is preferred
when the same taxon appears in multiple inputs):

```bash
python3 flyguide_merge_palaeo_sources.py \
  --inputs pbdb.csv neotoma.csv \
  --source-priority PBDB,Neotoma,NOW,GBIF \
  --out merged.csv --write-flyguide-files
```

Then run FlyGuide:

```bash
./flyguide.sh \
  NH_quaternary_palaeo_merged_gbif.csv \
  NH_quaternary_palaeo_refs \
  your_email@example.org \
  YOUR_NCBI_API_KEY
```

---

## Recommended source strategy

For Quaternary/aDNA-focused work:

1. **Neotoma**: Quaternary palaeoecological occurrence taxa.
2. **PBDB**: broad global fossil coverage, especially outside Neotoma's strengths
   (invertebrates, deep time, global tropics).
3. **NOW**: strong mammal supplement, especially Eurasia/Africa/Cenozoic mammals.
4. **GBIF fossil specimens**: optional museum-specimen supplement.

Merge everything by binomial by default, keep source metadata, and let the NCBI step
fail gracefully for names with no sequence data.

---

## Output files

For `example.csv`, each exporter writes:

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

The rejected file documents which source names were too messy or ambiguous to convert
into NCBI search names — worth reviewing before running large downloads.
