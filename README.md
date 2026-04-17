# TackleBox

<p align="center">
  <img src="assets/TackleBox.png" alt="TackleBox header" width="700">
</p>

***TackleBox*** is a collection of modules for ancient DNA research developed at the Hakai Institute in Campbell River, BC. Targeted ancient DNA workflows share more than a few similarities with fishing (especially capture enrichment with baits/probes). The tools you use to target your DNA of interest are often a major determinant of downstream success. TackleBox is intended to fill practical gaps in the palaeogenomics software landscape and help researchers maximize their chances of catching the big one. May your lines be tight!

## Current modules

- `FlyGuide/` – Build regional NCBI nucleotide reference panels from GBIF exports, with an optional NCBI nucleotide “coverage check” (`GuideCheck`).
- `FlyForge/` – Design RNA bait panels with customizable parameters for in-house bait synthesis, including design summaries, FASTA outputs, and oligo-pool files for ordering. Also includes `FlyForgeAudit` for evaluating pre-designed panels and augmenting existing panels with new targets.
- `MetaMerge/` – Merge MEGAN taxonomic classifications with Holi/metaDMG classifications and damage profiles, along with other associated sample metadata (depth, age, site, sample type, etc.) to produce taxon-by-taxon aDNA damage support spreadsheets with associated metadata and auto-generated heatmaps and stacked bar plots.

## Planned modules

- `Spinner/` – Filter regional reference FASTAs from larger databases to remove problematic public sequences.
- `Fillet/` – BLAST/LAST/MEGAN-based metagenomic classification pipeline.

Each module has its own README with installation instructions, dependencies, usage, and examples.
