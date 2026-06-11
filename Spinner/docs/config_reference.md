# Spinner Configuration Reference

This document describes every config key in the Spinner YAML config files.
Defaults are defined in `spinner/config.py` and can be overridden by any user
YAML — only keys you want to change need to appear in your file.

All defaults are calibrated for **ancient DNA metagenomics and capture
enrichment** as of v0.7.0.

---

## `run`

General runtime settings.

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `profile_name` | str | `"default"` | Informational label written to logs. |
| `progress_bars` | bool | `true` | Print progress bars to stderr. |
| `write_reject_fasta` | bool | `true` | Write `OUTPREFIX.reject.fasta`. |
| `write_review_fasta` | bool | `true` | Write `OUTPREFIX.review.fasta`. |
| `fail_on_missing_external_tool` | bool | `false` | Exit non-zero if blastn/vsearch is missing. If false, warn and skip. |
| `keep_temp_files` | bool | `false` | Retain temporary files (`*.tmp.keyed.fasta`, `*.windows.fasta`). Useful for debugging. |

---

## `inputs`

Default paths for config-specified input files. CLI flags (`--adapters`, etc.)
override these.

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `regions_config` | str | `""` | Path to regions_config TSV. |
| `species_kingdom` | str | `""` | Path to species_kingdom TSV. |
| `adapters_tsv` | str | `""` | Path to adapters TSV. |
| `bad_keywords_tsv` | str | `""` | Path to bad_keywords TSV. |

---

## `steps`

Toggle whole pipeline stages on or off.

| Key | Default | When to enable | Requires |
|-----|---------|----------------|----------|
| `basic_qc` | `true` | Always | — |
| `classify_regions` | `true` | Always | — |
| `bad_keyword_screen` | `true` | Almost always | `--bad-keywords` TSV |
| `adapter_screen` | `true` | Almost always | `--adapters` TSV |
| `vector_screen` | `false` | Publication / bait-design | UniVec BLAST DB |
| `taxonomy_blast` | `false` | Taxonomy sanity checking | BLAST DB or MMSeqs2 DB |
| `windowed_blast` | `false` | Chimerism detection in long seqs | BLAST DB |
| `chimera_screen` | `false` | De novo chimera detection | vsearch |
| `cluster` | `false` | Reducing haplotype redundancy | vsearch |
| `cap_references` | `false` | Controlling over-representation | — |
| `fcs_adaptor` | `false` | Assembly-grade adapter check | FCS-adaptor hook |
| `fcs_gx` | `false` | Assembly contamination | FCS-GX hook |
| `report` | `true` | Always useful | — |

---

## `basic_qc`

Per-sequence quality checks. Failures add reasons that influence the decision.

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `min_length` | int | `50` | Minimum sequence length (bp) → `length_below_min` → REJECT. |
| `max_length` | int | `400000` | Maximum sequence length (bp) → `length_above_max` → REJECT. Increase for full organelle genomes. |
| `min_length_by_class` | dict | `{}` | Per-marker-class length floors. E.g. `{NucMark: 200, Mito: 100}`. Adds reason `length_below_class_min` → REJECT. Applied in addition to `min_length`. |
| `max_n_fraction` | float | `0.20` | Maximum fraction of N bases (0.0–1.0). Above this → `n_fraction_high` (score penalty only; **not** a hard reject). |
| `max_non_iupac_fraction` | float | `0.01` | Maximum fraction of non-IUPAC characters → `non_iupac_fraction_high` → REJECT. |
| `min_shannon_entropy` | float | `1.2` | Minimum Shannon entropy of ACGT bases → `low_complexity` → REVIEW. 0 = mono-nucleotide run; 2.0 = perfectly even. |
| `max_homopolymer_run` | int | `60` | Maximum homopolymer run length (bp) → `homopolymer_long` → REVIEW. N runs are excluded. |
| `remove_exact_duplicate_sequences` | bool | `true` | Flag second+ records with identical sequences (by SHA-256) → `duplicate_sequence` → REJECT. |
| `remove_duplicate_accessions` | bool | `true` | Flag second+ records with the same accession → `duplicate_accession` → REJECT. |

---

## `classification`

Marker/region classification behaviour.

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `default_class_if_no_match` | str | `"Other"` | Class assigned when no region rule or builtin heuristic matches. |
| `coerce_unknown_plastid_to_plant` | bool | `true` | If kingdom is Unknown but class is Plastid, coerce kingdom to Plant. |

---

## `adapter_screen`

Controls adapter/primer contamination scanning.

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `max_mismatch` | int | `2` | Global maximum Hamming mismatches per adapter match. Per-adapter overrides from the TSV take precedence. |
| `scan_reverse_complement` | bool | `true` | Also scan the reverse complement of each adapter. |
| `reject_internal_hits` | bool | `true` | Internal adapter hits → `adapter_internal` → REJECT. |
| `review_terminal_hits` | bool | `true` | Terminal adapter hits → `adapter_terminal` → REVIEW. |
| `terminal_window_bp` | int | `25` | Positions within this distance of either end are terminal. |
| `min_adapter_match_length` | int | `12` | Ignore adapters shorter than this after loading from TSV. |
| `default_action` | str | `"reject"` | Fallback action if not specified per-adapter in TSV. |

---

## `bad_keywords`

Controls FASTA header keyword scanning.

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `case_sensitive` | bool | `false` | If false, matching is case-insensitive. |
| `default_action` | str | `"review"` | Fallback action if keyword has no explicit action in the TSV. |

The bundled `configs/bad_keywords.tsv` is aDNA-optimised.  `UNVERIFIED`,
`PREDICTED`, and `low quality` are review (not reject) — these labels appear
on many legitimate ancient DNA submissions.

---

## `vector_screen`

BLAST-based vector / UniVec contamination screening.

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `blast_db` | str | `""` | Path to a blastn-formatted UniVec or custom vector database. |
| `blast_task` | str | `"blastn-short"` | BLAST task. Use `blastn-short` for short adapter-like sequences. |
| `min_hit_length` | int | `20` | Ignore hits shorter than this (bp). |
| `min_pident` | float | `90.0` | Minimum % identity to call a vector hit. |
| `reject_internal_hits` | bool | `true` | Internal hits → `vector_internal` → REJECT. |
| `review_terminal_hits` | bool | `true` | Terminal hits → `vector_terminal` → REVIEW. |
| `terminal_window_bp` | int | `25` | Terminal zone width in bp. |
| `outfmt` | str | (see config) | BLAST output format. Do not change unless you update the parser. |

---

## `taxonomy_blast`

Taxonomy sanity checking via BLAST or MMSeqs2.

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `method` | str | `"blastn"` | Search tool: `"blastn"` or `"mmseqs2"`. Both produce identical output for the parser. |
| `blast_db` | str | `""` | Path to BLAST DB or MMSeqs2 taxonomy-indexed DB. |
| `blast_task` | str | `"megablast"` | BLAST task (ignored for MMSeqs2). |
| `num_threads` | int | `1` | CPU threads. Set to 8–32 on a cluster. |
| `max_target_seqs` | int | `10` | Maximum hits per query. |
| `max_hsps` | int | `5` | Maximum high-scoring pairs per subject (BLAST only). |
| `evalue` | str | `"1e-10"` | E-value cutoff. |
| `min_qcov` | float | `50.0` | Minimum query coverage % to consider a hit. |
| `min_pident` | float | `70.0` | Minimum % identity. Ancient sequences tolerate lower identity. |
| `expected_taxon_level` | str | `"genus"` | Comparison level: `"species"` or `"genus"`. |
| `reject_cross_kingdom` | bool | `true` | Reject if top hit is in a different kingdom (accurate only with taxdump). |
| `review_if_no_expected_match` | bool | `false` | Flag with REVIEW when no hit matches expected taxon. Off by default: rare taxa are often absent from NCBI databases. |
| `taxdump_dir` | str | `""` | Path to NCBI taxdump directory containing `nodes.dmp` and `names.dmp`. Enables lineage-aware kingdom comparison; handles NCBI 2024+ rank reorganisation. |
| `outfmt` | str | (see config) | BLAST outfmt 6 columns. Do not change unless you update the parser. |

### MMSeqs2 notes

For taxonomy columns (`taxid`, `taxname`) to work, the database must be
indexed with `mmseqs createtaxdb`.  Without taxonomy indexing, string-matching
mode is used automatically.  The `--format-output` flag used by Spinner
produces tab-separated columns in the same order as the BLAST parser expects.

### Taxonomy status values

| Status | Meaning |
|--------|---------|
| `PASS_SPECIES` | Top hit's scientific name contains the expected species. |
| `PASS_GENUS` | Top hit's scientific name contains the expected genus. |
| `NO_EXPECTED_MATCH` | No qualifying hit matched the expected taxon. |
| `NOT_CHECKED` | No search was run (step disabled or no hits above thresholds). |
| `REJECT_CROSS_KINGDOM` | Hit kingdom differs from expected kingdom (requires taxdump). |

---

## `windowed_blast`

Chimerism / mixed-origin detection for long sequences.

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `blast_db` | str | `""` | BLAST DB for window queries. Falls back to `taxonomy_blast.blast_db` if empty. |
| `blast_task` | str | `"megablast"` | BLAST task. |
| `enabled_for_min_length` | int | `1000` | Only apply windowed BLAST to sequences at least this long. |
| `window_size` | int | `500` | Sliding window size (bp). |
| `window_step` | int | `250` | Step between windows (bp). |
| `max_target_seqs` | int | `3` | BLAST hits per window. |
| `max_hsps` | int | `1` | Maximum HSPs per subject (BLAST only). |
| `evalue` | str | `"1e-10"` | E-value cutoff. |
| `min_pident` | float | `80.0` | Minimum % identity to use a window hit. |
| `num_threads` | int | `1` | CPU threads for windowed BLAST. |
| `min_conflicting_windows` | int | `2` | Minimum discordant windows to flag a record as conflicting. |
| `conflict_action` | str | `"review"` | Action for flagged records: `"review"` or `"reject"`. |
| `taxdump_comparison_rank` | str | `"family"` | Rank for lineage-aware conflict detection when taxdump is loaded. Options: `"genus"`, `"family"`, `"order"`, `"class"`. Falls back to genus-string matching if taxdb is unavailable. |
| `outfmt` | str | (see config) | BLAST outfmt 6 columns. Do not change unless you update the parser. |

---

## `chimera_screen`

vsearch uchime chimera detection.

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `method` | str | `"uchime_denovo"` | `"uchime_denovo"` (no reference DB needed) or `"uchime_ref"` (reference-based, higher sensitivity). |
| `reference_db` | str | `""` | Path to reference FASTA for `uchime_ref` mode. |
| `vsearch_path` | str | `"vsearch"` | Path or command name for vsearch. |
| `reject_chimeras` | bool | `true` | Y verdict → `chimera_detected` → REJECT. |
| `review_borderline` | bool | `true` | ? verdict → `chimera_borderline` → REVIEW. |
| `abskew` | float | `2.0` | vsearch `--abskew` parameter. Lower values increase sensitivity. |

---

## `cluster`

vsearch / cd-hit-est clustering to reduce haplotype redundancy.

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `method` | str | `"vsearch"` | Clustering tool: `"vsearch"` or `"cdhit"`. |
| `identity` | float | `0.99` | Sequence identity threshold (0.0–1.0). |
| `by` | list | `["species_guess", "marker_class"]` | Group records by these fields before clustering. |
| `vsearch_path` | str | `"vsearch"` | Path or command name for vsearch. |
| `cdhit_path` | str | `"cd-hit-est"` | Path or command name for cd-hit-est. |
| `max_reps_per_cluster` | int | `1` | Maximum centroid representatives per cluster. |
| `nonrepresentative_action` | str | `"review"` | Action for non-centroid members: `"review"` or `"reject"`. |

---

## `capping`

Species × marker representation capping.  Disabled by default for aDNA work
(all haplotype diversity is valuable for capture enrichment and mapping).

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `mode` | str | `"species_marker"` | `"species_marker"` (per species × class) or `"species"` (per species total). |
| `uncapped_classes` | list | `[]` | Marker classes exempt from capping. |
| `max_per_species_marker` | dict | `{Mito: 20, Plastid: 20, NucMark: 10, Other: 5}` | Maximum records per species per marker class when `cap_references: true`. |
| `max_per_species_total` | int | `50` | Maximum total records per species across all classes. |
| `cap_action` | str | `"review"` | Action for records exceeding cap: `"review"` or `"reject"`. |
| `rescue_sole_representatives` | bool | `true` | After all filtering, if a species has 0 KEEP records, promote its best REVIEW record to KEEP with reason `sole_representative`. Critical for rare/extinct taxa. |

---

## `decision_rules`

Controls which reasons force REJECT and which prevent KEEP.

| Key | Type | Description |
|-----|------|-------------|
| `hard_reject_reasons` | list | If any of these reasons is present, record is REJECT regardless of score. |
| `review_reasons` | list | If any of these reasons is present, record cannot be KEEP even with a high score. |
| `score_thresholds.keep_min` | int | Minimum score to be KEEP (default **65**). |
| `score_thresholds.review_min` | int | Minimum score to be REVIEW (default **30**). Below this → REJECT. |

### Default `hard_reject_reasons`

```
adapter_internal, vector_internal, duplicate_accession, duplicate_sequence,
length_below_min, length_above_max, length_below_class_min,
non_iupac_fraction_high, bad_keyword_reject, taxonomy_cross_kingdom,
chimera_detected, fcs_adaptor_hit, fcs_gx_contaminant
```

Note: `n_fraction_high` is intentionally **absent** from hard reject reasons.
High-N sequences from rare or extinct species receive a score penalty only.

### Default `review_reasons`

```
adapter_terminal, vector_terminal, low_complexity, homopolymer_long,
bad_keyword_review, taxonomy_no_expected_match, taxonomy_not_checked,
windowed_blast_conflict, chimera_borderline, cluster_nonrepresentative,
cap_exceeded, fcs_adaptor_review, fcs_gx_review
```

---

## `scoring`

Score adjustments applied per reason.  All values are integers added to the
starting score (`start: 100`).  Penalties are negative; bonuses are positive.

| Reason | Default | Notes |
|--------|---------|-------|
| `start` | `100` | Initial score for every record. |
| `adapter_internal` | `−100` | Hard reject |
| `adapter_terminal` | `−40` | Review reason |
| `vector_internal` | `−100` | Hard reject |
| `vector_terminal` | `−50` | Review reason |
| `duplicate_accession` | `−100` | Hard reject |
| `duplicate_sequence` | `−100` | Hard reject |
| `length_below_min` | `−100` | Hard reject |
| `length_above_max` | `−100` | Hard reject |
| `length_below_class_min` | `−100` | Hard reject |
| `n_fraction_high` | `−20` | **Soft penalty only** — not in hard_reject_reasons |
| `non_iupac_fraction_high` | `−40` | Hard reject |
| `low_complexity` | `−25` | Review reason |
| `homopolymer_long` | `−20` | Review reason |
| `bad_keyword_reject` | `−100` | Hard reject |
| `bad_keyword_review` | `−30` | Review reason |
| `taxonomy_cross_kingdom` | `−100` | Hard reject |
| `taxonomy_no_expected_match` | `−10` | Review reason; reduced penalty since rare taxa are often absent |
| `taxonomy_same_species` | `+20` | Bonus |
| `taxonomy_same_genus` | `+10` | Bonus |
| `refseq_preferred` | `+10` | Accession starts with NC_, NM_, NR_, NZ_, NG_, XM_, XR_, etc. |
| `voucher_keyword` | `+5` | Header contains "voucher", "holotype", "type strain", etc. |
| `complete_organelle` | `+5` | Header contains "complete genome", "complete cds", etc. |
| `chimera_detected` | `−100` | Hard reject |
| `chimera_borderline` | `−60` | Review reason |
| `sole_representative` | `+10` | Rescued as only KEEP for its species × class |
| `cluster_representative` | `+5` | vsearch cluster centroid |
| `cluster_nonrepresentative` | `−25` | Review reason |
| `cap_exceeded` | `−30` | Review reason (only relevant if `cap_references: true`) |
| `windowed_blast_conflict` | `−60` | Review reason |
| `fcs_adaptor_hit` | `−100` | Hard reject |
| `fcs_adaptor_review` | `−40` | Review reason |
| `fcs_gx_contaminant` | `−100` | Hard reject |
| `fcs_gx_review` | `−40` | Review reason |

---

## `fcs_adaptor` / `fcs_gx`

External hook configuration. Spinner does not assume a fixed output format.
Post-process tool output into a simple TSV with `accession` (or `qseqid`) and
`status` columns, then point the config at it.

| Key | Type | Description |
|-----|------|-------------|
| `command` | str | Shell command template. Placeholders: `{input_fasta}`, `{outprefix}`, `{label}`. |
| `results_tsv` | str | Path to a normalised TSV. If empty, defaults to `OUTPREFIX.{label}.tsv`. |
| `reject_on_any_hit` | bool | (`fcs_adaptor`) Reject on any hit rather than only internal hits. |
| `tax_id` | str | (`fcs_gx`) NCBI taxid for the expected organism. |
| `reject_divisions` | list | (`fcs_gx`) FCS-GX contamination divisions to hard-reject. |

---

## Example: minimal aDNA metagenomics run

```yaml
# project.yml
run:
  profile_name: "my_project"

basic_qc:
  min_length: 100
  min_length_by_class:
    NucMark: 150
    Mito: 100

steps:
  taxonomy_blast: true
  windowed_blast: false    # only useful for seqs >= 1000 bp

taxonomy_blast:
  method: blastn           # or mmseqs2 for faster runs
  blast_db: /data/databases/nt/nt
  taxdump_dir: /data/databases/taxdump
  num_threads: 16
  min_pident: 70.0
```

Run with:

```bash
./Spinner filter \
  --fasta my_refs.fasta \
  --species-kingdom my_species_kingdom.tsv \
  --adapters configs/adapters_default.tsv \
  --bad-keywords configs/bad_keywords.tsv \
  --config project.yml \
  --outprefix results/my_project
```
