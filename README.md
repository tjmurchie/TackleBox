# TackleBox

<p align="center">
  <img src="assets/TackleBox.png" alt="TackleBox header" width="700">
</p>

***TackleBox*** is a collection of modules for ancient DNA research developed at the Hakai Institute in Campbell River, BC. Targeted ancient DNA workflows share more than a few similarities with fishing (especially capture enrichment with baits/probes). The tools you use to target your DNA of interest are often a major determinant of downstream success. TackleBox is intended to fill practical gaps in the palaeogenomics software landscape and help researchers maximize their chances of catching the big one. May your lines be tight!

## Current modules

- `FlyGuide/` – Build regional NCBI nucleotide reference panels from GBIF exports or Neotoma paleoecological data, with an optional NCBI nucleotide “coverage check” (`GuideCheck`). Includes a Neotoma exporter (`neotoma_extinct_to_gbif.py`) for generating extinct/paleo taxon lists from broad spatial and temporal queries.
- `FlyForge/` – Design RNA bait panels with customizable parameters for in-house bait synthesis, including design summaries, FASTA outputs, and oligo-pool files for ordering. Also includes `FlyForgeAudit` for evaluating pre-designed panels and augmenting existing panels with new targets.
- `MetaMerge/` – Merge MEGAN taxonomic classifications with Holi/metaDMG classifications and damage profiles, along with other associated sample metadata (depth, age, site, sample type, etc.) to produce taxon-by-taxon aDNA damage support spreadsheets with associated metadata and auto-generated heatmaps and stacked bar plots.
- `Spinner/` – Curate NCBI and custom reference FASTAs before use in bait design, short-read mapping, or metagenomic classification. Optimised for ancient DNA shotgun metagenomics and capture enrichment across all kingdoms. Every input record is audited and assigned to KEEP / REVIEW / REJECT with a full reasons trail — nothing is silently discarded. Features include: per-sequence QC (length, N-fraction, complexity, duplicates), adapter and vector contamination scanning, BLAST/MMSeqs2 taxonomy sanity checking with optional NCBI taxdump lineage support, windowed chimerism detection for long sequences, vsearch uchime chimera screening, sole representative rescue for rare/extinct taxa, and score-based decisions with configurable thresholds. Defaults are calibrated for ancient DNA (high-N sequences are penalised, not rejected; capping disabled; UNVERIFIED records flagged for review rather than rejection).

## Planned modules

- `Fillet/` – BLAST/LAST/MEGAN-based metagenomic classification pipeline.

Each module has its own README with installation instructions, dependencies, usage, and examples.
