# TackleBox: FlyForge

<p align="center">
  <img src="../assets/FlyForge.png" alt="TackleBox: FlyForge header" width="700">
</p>

**FlyForge** is a bait/probe designer for hybridization capture enrichment with integrated support for **in-house RNA bait synthesis from DNA oligo pools**. It is built for both focused projects (for example, single mitogenomes or locus panels) and larger multi-target panels.

FlyForge is part of the **TackleBox** suite.

## What FlyForge does

- tiles probes at a user-defined density
- preprocesses references (`U -> T`, short `N` runs -> `T`, optional repeat soft-masking)
- optionally removes complementary regions from the input references in the CARPDM style
- filters probes by ambiguity, masked fraction, Tm, internal LguI/BspQI sites, perfect reverse complements, BLAST off-targets, self-BLAST redundancy, and `cd-hit-est` clustering
- constructs **order-ready oligo-pool sequences** for in-house synthesis
- validates the final bare bait set back against the input references with per-target coverage plots and summary tables

Validation note: final validation assigns each bait **one primary on-target placement**. For FlyForge-generated bait IDs that encode source coordinates, those coordinates are used to avoid artificial coverage inflation in repetitive regions (for example mitochondrial control-region repeats).

## Oligo-pool design model

FlyForge follows the CARPDM oligo-pool strategy for in-house synthesis of RNA baits.

### Final oligo structure

```text
5'-GCTAATACGACTCACTATAGGG [probe] [reverse-complement of amplification primer]-3'
```

Where:

- `GCTAATACGACTCACTATAGGG` is the **22 nt T7 promoter-containing primer sequence** used in CARPDM
- the selected amplification primer is always designed to end in `GCTCTTCG`
- therefore its reverse complement begins with `CGAAGAGC`, producing the intended CARPDM-compatible BspQI/LguI-associated tail on the final oligo

### Wet-lab logic

The intended molecular workflow is:

1. PCR amplify the oligo pool
2. digest with **LguI / BspQI** to generate the correct transcription template boundary
3. perform **T7 in vitro transcription**
4. purify the RNA bait pool

## Requirements

### Python

- Python 3.8+
- `biopython`
- `primer3-py`
- `matplotlib`
- `seaborn`
- `pandas`
- `numpy`
- `tqdm`

### External tools

- **BLAST+** (`blastn`, `makeblastdb`)
- **CD-HIT** (`cd-hit-est`)

## Installation

```bash
conda create -n flyforge python=3.12
conda activate flyforge
conda install -c bioconda blast cd-hit
pip install biopython primer3-py matplotlib seaborn pandas numpy tqdm

chmod +x FlyForge.py FlyForge FlyForgeAudit.py FlyForgeAudit
```

## Testing

```bash
pip install pytest
cd FlyForge
pytest              # fast tests (tiling math, Tm calculation) -- runs in ~1 second
pytest -m slow       # also runs the real oligo-pool assembly test (requires blastn, ~1-4 minutes)
```

The fast suite covers the tiling window/step arithmetic (including circular-genome
wraparound), `compute_tm()`'s behavior on both normal and pathological input,
`--circular-ids` resolution, self-repeat masking (both its own pure-function
correctness and the `--skip-self-mask`/`--no-skip-self-mask` CLI default), and
`--redundancy-prereduce-identity`'s `>=0.80` validation boundary, and the
`--min-cluster-identity`/`--hard-max-baits` binary-search logic and CLI validation. The
`slow` suite runs `design_opool()` end to end against real BLAST+/primer3 and verifies
the literal oligo sequence structure that gets sent for synthesis
(`5'-[T7 promoter][probe][primer reverse-complement]-3'`), including that the probe
sequence itself survives unmodified in the final assembled oligo, plus real end-to-end
`--redundancy-prereduce-identity` and `--hard-max-baits` runs against actual
`cd-hit-est`/`vsearch`/`blastn`.

## Quick start

### Standard FlyForge design

```bash
FlyForge \
  -i reference.fasta \
  --prefix my_panel \
  --output-dir flyforge_output \
  --bait-length 80 \
  --tiling-density 4 \
  --min-tm 50 \
  --threads 8
```

### Standard design with a custom exclusion database

```bash
FlyForge \
  -i targets/*.fasta \
  --prefix env_panel \
  --output-dir flyforge_output \
  --bait-length 80 \
  --tiling-density 3 \
  --min-tm 50 \
  --remove-complements \
  --blast-db /path/to/exclusion_db \
  --blast-min-pident 80 \
  --blast-max-hits 5 \
  --threads 16
```

### Skip oligo-pool construction

```bash
FlyForge \
  -i reference.fasta \
  --prefix bare_probes_only \
  --no-opool
```

## Core outputs

A typical run produces:

- `PREFIX_final_baits.fa` — final bare bait sequences
- `PREFIX_probes.fna` — bare bait FASTA used for validation and QC
- `PREFIX_oligo_pool.fna` — full order-ready oligo sequences
- `PREFIX_amplification_primers.fna` — T7 primer and selected second primer
- `PREFIX_final_blast.xml` — validation BLAST output
- `PREFIX_target_info.csv` — per-target validation metrics
- `PREFIX_probe_info.csv` — per-probe QC metrics
- `PREFIX_target_probe_pairs.csv` — target:probe pairings from validation BLAST
- `PREFIX_per_ref_stats.tsv` — per-reference bait counts plus **final validated** coverage statistics
- `PREFIX_summary.tsv` — parameters and step-wise run summary
- `PREFIX_plots/` — violin plots and individual target coverage plots
- `PREFIX_progress.log` — log of the full run

## Important behavior in this version

### Validation BLAST

FlyForge validates the final bare probes against the input references with BLAST using `-dust no` and `-soft_masking false`. This avoids undercounting valid probe hits in low-complexity regions during the post-design validation step.

### O-pool failure behavior

FlyForge now **fails loudly** if it cannot identify a valid second amplification primer or if primer-selection BLAST output cannot be parsed. That is deliberate: a failed o-pool design should not silently produce something that could be mistaken for an order-ready synthesis file.

### Self-repeat masking default (v1.3.0+)

`--skip-self-mask` now defaults to **skipped** (previously masking was on by default).
Self-repeat masking counts k-mers across the *entire* input FASTA, not per-reference --
correct for a genuinely small/few-genome reference set (it finds real within-genome
repeats: tandem repeats, transposons, low-complexity runs), but at real multi-species
panel scale it stops detecting genuine repeats and instead flags ordinary,
cross-species-conserved coding sequence as "repetitive" -- discarding exactly the loci
that let one bait cross-capture DNA from related species. Pass `--no-skip-self-mask` to
restore the original masked-by-default behavior, which is still the right call for a
small/few-genome reference set (see the single-organism example below).

## Recommended use patterns

### Single-organism ancient DNA or degraded DNA capture

A small/few-genome reference set is exactly where self-repeat masking is still the right
call (see above) -- `--no-skip-self-mask` restores it explicitly, since the tool-wide
default is now skipped:

```bash
FlyForge \
  -i reference.fasta \
  --prefix my_species \
  --bait-length 80 \
  --tiling-density 4 \
  --min-tm 50 \
  --remove-complements \
  --no-skip-self-mask \
  --circular \
  --threads 8
```

### Multi-species or environmental panel design

Self-repeat masking is skipped by default here (the right call at this scale -- see
above), so no extra flag is needed. `--max-baits`/`--min-tiling-density` are worth
setting explicitly for a real affordability/pool-size target -- FlyForge auto-adjusts
tiling density to fit, rather than tiling at a fixed density and hoping the result lands
near budget:

```bash
FlyForge \
  -i targets/*.fasta \
  --prefix multi_taxon_panel \
  --bait-length 80 \
  --tiling-density 3 \
  --min-tm 50 \
  --remove-complements \
  --blast-db /path/to/exclusion_db \
  --blast-min-pident 80 \
  --blast-max-hits 5 \
  --cluster-identity 0.95 \
  --max-baits 20000 \
  --min-tiling-density 1.0 \
  --threads 16
```

### Very large multi-species panels (tens of thousands of references)

At real large-scale, masking off (the new default above) can leave `self_blast_filter()`
facing a candidate pool so large its own redundancy-removal step takes hours (confirmed:
2h39m on a real 726,110-candidate panel). `--redundancy-prereduce-identity` (v1.4.0+,
default disabled) adds a cheap `cd-hit-est` pre-pass immediately before that step,
cutting the same real case to under 90 seconds while landing closer to budget besides.
Must be `>=0.80` (cd-hit-est's own hard floor in nucleotide mode):

```bash
FlyForge \
  -i targets/*.fasta \
  --prefix multi_taxon_panel \
  --bait-length 80 \
  --tiling-density 3 \
  --min-tm 50 \
  --remove-complements \
  --max-baits 20000 \
  --min-tiling-density 1.0 \
  --redundancy-prereduce-identity 0.80 \
  --threads 16
```

### Hard bait-count cap (Twist-tier-aligned affordability targets)

When a panel needs to fit a specific, non-negotiable bait count (e.g. a Twist pricing
tier boundary), `--min-cluster-identity`/`--hard-max-baits` (v1.5.0+, both required
together, default disabled) replace the trio above entirely with a single
divergence-bounded clustering step. Give it the loosest divergence tolerance still
scientifically acceptable for the project (e.g. `0.75` = 25% divergence) as
`--min-cluster-identity`, and the real ceiling as `--hard-max-baits` -- the pipeline
binary-searches identity within that floor for the tightest cap-respecting result, and
reports plainly (`INFEASIBLE` in the progress log) if even the loosest tolerance can't
reach the target, rather than silently exceeding it or crashing:

```bash
FlyForge \
  -i targets/*.fasta \
  --prefix affordable_panel \
  --bait-length 80 \
  --tiling-density 3 \
  --min-tm 50 \
  --remove-complements \
  --min-cluster-identity 0.75 \
  --hard-max-baits 60000 \
  --threads 16
```

If this reports infeasible at your real project's scale, that's a real signal the target
needs less scope (fewer species/loci/markers), not a tighter identity search -- reducing
reference diversity to hit a bait-count budget is not something this flag will do for
you.

## Companion module: FlyForgeAudit

This repository now also includes **FlyForgeAudit**, a companion module for:

- auditing an existing bait set against one or more references
- screening a panel against an avoid database
- designing the **minimal spike-in bait set** needed to extend an existing panel to new organisms

See `FlyForgeAudit_README.md` for full documentation.

## Citation

If you use FlyForge, cite this software and the CARPDM publication that established the in-house oligo-pool synthesis strategy.

## Recent updates
- `--min-cluster-identity`/`--hard-max-baits` (v1.5.0) provide a real, guaranteed hard cap
  on final bait count via divergence-bounded vsearch clustering, replacing the old soft
  `--max-baits`/masking-based size controls for panels that must fit an exact budget.
- Circular targets can now be enabled with `--circular` (all references) or `--circular-ids ref1,ref2` (selected references). This is useful for mitochondrial, plastid, viral, plasmid, and other circular genomes so baits can wrap across the linearized ends.
- FlyForgeAudit now distinguishes **actionable flags** from **informational notes** in the terminal summary and prints the recommendation text directly at the end of the run.
- FlyForgeAudit now includes an `opool` mode that converts an existing bare-bait FASTA directly into an order-ready oligo pool and amplification-primer FASTA without redesigning the panel.
- FlyForgeAudit now uses a FlyForge-style progress dashboard so the current step, elapsed time, ETA, and completed stages are visible during audit and augment runs.
