# TackleBox: FlyGuide

**FlyGuide** is a TackleBox module for building region-specific NCBI nucleotide
reference panels from GBIF downloads, with optional NCBI "coverage checks" via
GuideCheck.

It was designed for metagenomic and sedimentary ancient DNA (sedaDNA) work,
where you want regionally sensible references to pre-map against (or to use as
BLAST databases) without hauling the entire nt database around.

---

## Overview

Given a GBIF occurrence or checklist download (CSV/TSV with `species`,
`genus`, and `kingdom` columns), FlyGuide:

1. Parses the GBIF file and builds:
   - a clean species/genus search list for NCBI,
   - a species↔kingdom mapping table.
2. Optionally runs **GuideCheck** to summarize what NCBI already has for each
   name (nuccore/SRA/assembly counts).
3. Uses `NCBI-NT_Downloader.pl` to fetch curated mitochondrial, plastid, and
   selected nuclear marker sequences from NCBI Nucleotide.
4. Splits the combined FASTA by kingdom and molecule type
   (mitochondrial / plastid / nuclear markers / other).

Outputs include:

- `OUTPREFIX_species_search.txt`
- `OUTPREFIX_species_kingdom.tsv`
- `OUTPREFIX_ncbi_guidecheck.tsv` (if GuideCheck enabled)
- `OUTPREFIX.fasta` (combined NCBI records)
- `OUTPREFIX.<Kingdom>-<Type>.fasta` (per-kingdom, per-type FASTAs)

---

## Scripts in this module

- `flyguide.sh`  
  Main wrapper that drives the full FlyGuide pipeline.

- `gbif_prep_from_csv.py`  
  GBIF pre-processor: reads a GBIF CSV/TSV and writes the search list plus
  species↔kingdom table.

- `NCBI-NT_Downloader.pl`  
  Perl script that queries NCBI Nucleotide via E-utilities and downloads a
  filtered, deduplicated FASTA.

- `split_fasta_by_kingdom_organelle_simple.pl`  
  Splits the combined FASTA into `<Kingdom>-<Type>` FASTAs.

- `guidecheck.sh`  
  Standalone NCBI "coverage check" tool for taxon lists. Used internally by
  `flyguide.sh`, but can also be run on its own.

---

## Dependencies

- bash
- Python 3.7+ (for `gbif_prep_from_csv.py`)
- Perl 5.x with `Bio::DB::EUtilities` (for `NCBI-NT_Downloader.pl`)
- `curl`, `jq`, `python3` (for `guidecheck.sh`)
- Internet access to NCBI E-utilities

An NCBI API key is strongly recommended.

---

## Quickstart

1. Export a GBIF CSV/TSV for your region (see **Preparing GBIF downloads** below).
2. From the directory with your GBIF file, run:

```bash
/path/to/TackleBox/flyguide/flyguide.sh \
  RegionX_GBIF.csv \
  RegionX_refs \
  your.email@example.org \
  YOUR_NCBI_API_KEY
```

3. Use the `RegionX_refs.*.fasta` files as references, and
   `RegionX_refs_ncbi_guidecheck.tsv` to understand which taxa have good
   coverage in NCBI.

To skip GuideCheck:

```bash
/path/to/TackleBox/flyguide/flyguide.sh --no-guidecheck \
  RegionX_GBIF.csv RegionX_refs your.email@example.org YOUR_NCBI_API_KEY
```

This will:

- Build `RegionX_refs_species_search.txt`
- Build `RegionX_refs_species_kingdom.tsv`
- Optionally build `RegionX_refs_ncbi_guidecheck.tsv`
- Download NCBI nucleotide records into `RegionX_refs.fasta`
- Split that FASTA into per-kingdom/per-type FASTAs.

---

## Preparing GBIF downloads

FlyGuide expects a GBIF CSV/TSV that includes at least `species`, `genus`, and
`kingdom` columns.

The GBIF web interface can change over time, but a typical workflow is:

1. Go to https://www.gbif.org/ and log in.
2. Use the map to zoom to your study region and apply a **geographic filter**
   (e.g. polygon or bounding box).
  ![Example GBIF map extent](docs/images/gbif_extent-filters.png)
3. Add **taxonomic filters** if desired (e.g. only `Plantae` + `Animalia`).
4. Choose an **"Occurrence"** or **"Checklist"** download that includes
   taxonomic fields.
  ![Example GBIF downloading species list](docs/images/gbif_download-specieslist.png)
5. Start the download and wait for the ZIP to be ready.
6. Unzip the download; identify the main occurrence/checklist file
   (often a `.csv` or `.tsv`).
7. Confirm that it contains `species`, `genus`, and `kingdom` columns
   (case-insensitive).

---

## GuideCheck (standalone usage)

Although `guidecheck.sh` is called automatically by `flyguide.sh`, you can also
use it on its own with any taxon list:

```bash
# One taxon name per line in taxa.txt
bash guidecheck.sh -i taxa.txt -o taxa_ncbi_guidecheck.tsv --api-key YOUR_NCBI_API_KEY
```

The output TSV has columns:

- `query_name`
- `taxid`
- `matched_name`
- `rank`
- `nuccore`
- `sra`
- `assembly`
- `status`

The `status` field summarizes availability:

- `NO_TAXID`
- `ASSEMBLY`
- `HAS_NUCCORE_AND_SRA`
- `SRA_ONLY`
- `NUCCORE_ONLY`
- `NONE`

---

## Repository layout (within TackleBox)

If FlyGuide lives inside a larger TackleBox repository, a typical structure is:

```text
TackleBox/
  flyguide/
    flyguide.sh
    gbif_prep_from_csv.py
    NCBI-NT_Downloader.pl
    split_fasta_by_kingdom_organelle_simple.pl
    guidecheck.sh
    README.md
  spinner/         # future module for reference filtering
  flyforge/        # future bait design module for in-house baits
  fillet/   	   # future metagenomic pipeline
  README.md        # high-level TackleBox overview
  LICENSE
```

---

## Planned extensions

Planned FlyGuide extensions include:

- Support for sourcing taxon names from the **Neotoma Paleoecology Database**
  in addition to GBIF, and feeding those into the same NCBI workflow.
- Deeper integration with other TackleBox modules (e.g. Spinner, CatchToFillet)
  to streamline reference construction → mapping → classification.

If you use FlyGuide in a publication, please cite the associated paper and this
repository URL.
