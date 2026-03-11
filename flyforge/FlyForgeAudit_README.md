# TackleBox: FlyForgeAudit

**FlyForgeAudit** is the companion analysis and augmentation module for **FlyForge**.

It has two modes:

- `audit` — evaluate an **existing bare bait set** against one or more references and reproduce FlyForge-style validation outputs
- `augment` — start from an **existing bait set** plus **new target references**, then design the **minimal additional bait set** needed to bring those new targets up to the requested coverage threshold

## Why this exists

FlyForge designs a panel from scratch.

FlyForgeAudit answers the next practical questions:

- *How good is my current bait set against this reference set?*
- *Where are the coverage gaps?*
- *Are there obvious off-target hits against organisms I want to avoid?*
- *How many extra baits do I actually need to spike in for a new organism set?*

## Important input rule

FlyForgeAudit expects **bare bait sequences only** as the starting panel input.

Use:

- `*_final_baits.fa`
- or another FASTA containing just the bait sequences

Do **not** give it an oligo-pool FASTA that already includes the T7 and amplification-primer tails.

---

## Installation

FlyForgeAudit uses the same environment as FlyForge.

```bash
conda create -n flyforge python=3.12
conda activate flyforge
conda install -c bioconda blast cd-hit
pip install biopython primer3-py matplotlib seaborn pandas numpy tqdm

chmod +x FlyForge.py FlyForge FlyForgeAudit.py FlyForgeAudit
```

---

## Mode 1: audit

### What it does

- BLASTs the existing bait set against one or more references
- writes FlyForge-style validation outputs:
  - `target_info.csv`
  - `probe_info.csv`
  - `target_probe_pairs.csv`
  - `per_ref_stats.tsv`
  - coverage plots and violin plots
- optionally screens the panel against a curated avoid database
- writes recommendation tables and a summary file

### Example

```bash
FlyForgeAudit audit \
  --baits existing_panel.fa \
  --reference refs/*.fasta \
  --prefix panel_audit \
  --output-dir audit_out \
  --desired-coverage-depth 1 \
  --coverage-fraction-goal 0.95 \
  --avoid-db /path/to/exclude_db \
  --avoid-min-pident 80 \
  --threads 8
```

### Key audit outputs

- `PREFIX_final_baits.fa` — copy of the audited bait panel
- `PREFIX_final_blast.xml` — validation BLAST output
- `PREFIX_target_info.csv` — per-target coverage metrics
- `PREFIX_probe_info.csv` — per-probe QC metrics
- `PREFIX_target_probe_pairs.csv` — BLAST-derived target:probe pairings
- `PREFIX_per_ref_stats.tsv` — per-reference validated coverage summary
- `PREFIX_avoid_hits.tsv` — optional avoid-database hits
- `PREFIX_recommendations.tsv` and `PREFIX_recommendations.txt`
- `PREFIX_summary.tsv`
- `PREFIX_plots/`
- `PREFIX_progress.log`

---

## Mode 2: augment

### What it does

- audits how well the **existing panel** covers the **new target set**
- identifies bases that fall below the requested minimum coverage depth
- greedily proposes the smallest practical additional bait set to cover those deficits
- filters those new baits using FlyForge-style QC and optional BLAST filters
- screens proposed new baits against the **existing panel** to avoid duplicate / near-duplicate / complement conflicts
- outputs both:
  - the **extra bare bait FASTA**
  - the **extra order-ready oligo-pool FASTA** for synthesis
- writes a merged panel audit for the augmented panel against the new targets

### Example

```bash
FlyForgeAudit augment \
  --existing-baits existing_panel.fa \
  --new-targets new_refs/*.fasta \
  --prefix panel_spikein \
  --output-dir spikein_out \
  --min-existing-coverage 1 \
  --coverage-fraction-goal 0.95 \
  --min-tm 50 \
  --threads 8
```

### Extra outputs produced in augment mode

- `PREFIX_extra_final_baits.fa` — new spike-in baits only
- `PREFIX_extra_oligo_pool.fna` — order-ready extra oligo pool
- `PREFIX_extra_amplification_primers.fna` — primer pair for the extra oligo pool
- `PREFIX_existing_final_baits.fa` — copy of the original panel
- `PREFIX_final_baits.fa` — merged panel (existing + extra)
- `PREFIX_target_info.csv`, `PREFIX_probe_info.csv`, `PREFIX_per_ref_stats.tsv`, `PREFIX_plots/` — audit outputs for the merged panel against the new targets
- `PREFIX_extra_avoid_hits.tsv` — optional avoid-database hits for the extra bait set
- `PREFIX_avoid_hits.tsv` — optional avoid-database hits for the merged panel
- `PREFIX_recommendations.tsv`
- `PREFIX_summary.tsv`

---

## Coverage logic in augment mode

Augment mode uses the current bait panel to estimate coverage over the new targets, then calculates a **coverage deficit array** for positions below the requested threshold.

It then proposes new bait windows greedily by repeatedly choosing the window that covers the largest remaining deficit. After that, the proposed baits are filtered through the same kinds of checks used by FlyForge:

- ambiguous bases
- masked fraction
- Tm
- internal LguI/BspQI motifs
- perfect reverse-complement conflicts
- optional BLAST exclusion filtering
- optional CARPDM-style specificity filtering
- optional self-BLAST redundancy collapse
- optional `cd-hit-est` clustering

This process repeats for up to `--max-augment-iterations` rounds.

---

## Shared useful options

### Reference preparation

- `--remove-complements` — run FlyForge/CARPDM-style complementary-target cleanup before analysis
- `--skip-self-mask` — skip internal repeat masking
- `--repeat-k`, `--repeat-threshold` — repeat masker settings

### Off-target review

- `--avoid-db` — curated avoid-database BLAST screen
- `--avoid-min-pident` — percent identity threshold for flagging hits
- `--avoid-max-hits` — maximum reported hits per bait

### Augment-only practical options

- `--min-existing-coverage` — minimum target coverage depth to satisfy with the augmented panel
- `--max-augment-iterations` — maximum number of iterative spike-in rounds
- `--no-opool` — skip order-ready extra oligo-pool generation

### Proposed-bait filtering options in augment mode

- `--min-tm`
- `--ambiguous-cutoff`
- `--max-masked-frac`
- `--blast-db`
- `--blast-min-pident`
- `--blast-max-hits`
- `--specificity-db`
- `--cluster-identity`
- `--cluster-overlap`
- `--no-cluster`
- `--no-redundancy`

---

## Notes on interpretation

### `per_ref_stats.tsv` in audit mode

In audit mode, the `per_ref_stats.tsv` file is used as a **validated panel-performance summary**, not a full design-history table. It reports the current panel’s validated bait coverage/count metrics for each reference.

### Avoid-database hits

The avoid-database screen is meant for **review and triage**, not blind automatic deletion. Strong hits are a prompt to inspect whether those baits should be replaced, relaxed, or retained depending on the biological goals of the panel.

### Oligo-pool output in augment mode

The order-ready oligo-pool FASTA generated in augment mode applies to the **extra spike-in baits only**. That is the synthesis-ready file for the incremental add-on set.

---

## Typical workflows

### Audit an already ordered bait panel

```bash
FlyForgeAudit audit \
  --baits ordered_panel_final_baits.fa \
  --reference target_refs/*.fasta \
  --prefix ordered_panel_qc \
  --output-dir ordered_panel_qc_out \
  --threads 8
```

### Extend an existing panel to a new species set

```bash
FlyForgeAudit augment \
  --existing-baits existing_panel_final_baits.fa \
  --new-targets new_species/*.fasta \
  --prefix new_species_spikein \
  --output-dir new_species_spikein_out \
  --min-existing-coverage 1 \
  --min-tm 50 \
  --avoid-db /path/to/exclude_db \
  --avoid-min-pident 80 \
  --threads 8
```

---

## Relationship to FlyForge

- Use **FlyForge** when designing a panel from scratch.
- Use **FlyForgeAudit** when reviewing or extending an existing panel.
