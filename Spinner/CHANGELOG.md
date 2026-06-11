# Spinner Changelog

## v0.7.0 (2026-05-05)

Ancient DNA first — Spinner's default config is now calibrated for aDNA
metagenomics and capture enrichment.  All changes are backwards-compatible:
existing custom YAML configs still override any default.  Also adds MMSeqs2
as a fast BLAST alternative and lineage-aware windowed BLAST conflict detection.

### Changed (aDNA-first defaults)

- **`n_fraction_high` removed from `hard_reject_reasons`**: High-N sequences
  now receive a score penalty (−20) rather than an automatic REJECT.  A high-N
  reference for a rare or extinct species is better than no reference at all.
- **`max_n_fraction`**: 0.05 → 0.20.  Ancient assemblies routinely have 10–20%
  N in gap regions from low-coverage mapping.
- **`cap_references`**: True → False.  All haplotype diversity is retained by
  default; capping is counterproductive for capture-enrichment and mapping databases.
- **`rescue_sole_representatives`**: False → True.  If a species ends up with
  zero KEEP records, its best REVIEW record is promoted to KEEP automatically.
- **`score_thresholds.keep_min`**: 80 → 65.  Imperfect references are still
  valuable; preference is to REVIEW rather than REJECT.
- **`score_thresholds.review_min`**: 50 → 30.  Keeps the REJECT band tight.
- **`taxonomy_blast.min_pident`**: 80.0 → 70.0.  Divergent ancient sequences
  routinely have reduced identity to modern NCBI representatives.
- **`taxonomy_blast.review_if_no_expected_match`**: True → False.  Rare or
  extinct taxa are often absent from NCBI BLAST databases; not penalising them.
- **`scoring.n_fraction_high`**: −50 → −20 (soft penalty only).
- **`scoring.taxonomy_no_expected_match`**: −35 → −10 (rare taxa often absent).
- **`configs/bad_keywords.tsv`**: UNVERIFIED, PREDICTED, and low quality changed
  from reject → review.  Removed: uncultured, metagenome, environmental sample,
  whole genome shotgun, scaffold, contig, chromosome — all legitimate in aDNA
  reference databases.  Retained as reject: synthetic construct, vector.
- **`configs/spinner_ancient_metagenome.yml`**: Now documents the aDNA profile
  (most settings match the new defaults).

### Added

- **MMSeqs2 support** (`taxonomy_blast.method: mmseqs2`): `run_mmseqs()` in
  `external.py` calls `mmseqs easy-search` with `--format-output` matching the
  existing BLAST column layout — no changes to `parse_tax_blast()` needed.
  MMSeqs2 is 10–100× faster than BLAST on large nt databases.  Taxonomy columns
  (taxid, taxname) require a taxonomy-indexed MMSeqs2 database
  (`mmseqs createtaxdb`); without it, string-matching mode is used.
- **`taxonomy_blast.num_threads`** added to DEFAULT_CONFIG (was only in some
  YAML profiles).
- **`windowed_blast.num_threads`** added to DEFAULT_CONFIG.
- **Lineage-aware windowed BLAST conflict detection**
  (`windowed_blast.taxdump_comparison_rank: family`): when `taxdump_dir` is
  configured, windowed BLAST windows are now compared at the specified rank
  (default: `family`) using the loaded `TaxdumpDB` rather than genus-level
  string matching.  Falls back to string matching if taxdb is unavailable or
  a window's taxid is not found.  The taxdb loaded during the taxonomy BLAST
  stage is automatically passed to `parse_windowed_blast()` — no extra config
  needed.
- **4 new pytest tests** (138 → 142 total): windowed BLAST taxdb conflict,
  same-family no-conflict, missing taxid fallback, MMSeqs2 format compatibility.

---

## v0.6.0 (2026-05-05)

Major feature release for ancient DNA metagenomics and capture enrichment
pipelines.  All changes are backwards-compatible with existing configs.

### Added

- **Ancient DNA config profile** (`configs/spinner_ancient_metagenome.yml`):
  Calibrated for shotgun and capture-enrichment reference databases covering
  all kingdoms.  Key differences from default: capping disabled, N-fraction is
  a soft score penalty (not a hard reject), BLAST thresholds relaxed for
  divergent ancient sequences, sole-representative rescue enabled, lower keep/
  review score thresholds.
- **Ancient keyword list** (`configs/bad_keywords_ancient.tsv`):
  Removes `UNVERIFIED`, `PREDICTED`, `low quality`, `whole genome shotgun`,
  `scaffold`, `contig`, `uncultured`, `metagenome`, and `environmental sample`
  from the bad-keyword filter — all legitimate in aDNA reference databases.
  Only `synthetic construct` and `vector` remain as hard rejects.
- **Sole representative rescue** (`capping.rescue_sole_representatives`):
  After all filtering, any species × marker class group with zero KEEP records
  has its best-scoring REVIEW record promoted to KEEP with reason
  `sole_representative`.  Prevents losing the only reference for rare or
  extinct taxa.  Enabled by default in the aDNA config; off by default elsewhere.
- **Per-class minimum length** (`basic_qc.min_length_by_class`):
  Dict of marker-class-specific length floors, e.g. `{NucMark: 200, Mito: 100}`.
  Adds hard-reject reason `length_below_class_min`.  Global `min_length` still
  applies independently.
- **Gzipped FASTA input**: `parse_fasta()` transparently opens `.fasta.gz` files
  using `gzip.open` — no manual decompression needed.
- **Chimera detection via vsearch uchime** (`steps.chimera_screen`):
  New pipeline step that runs `vsearch --uchime_denovo` (no reference DB needed)
  or `vsearch --uchime_ref` (reference-based).  Chimeric sequences receive reason
  `chimera_detected` (hard reject by default); borderline sequences receive
  `chimera_borderline` (review).  Config section: `chimera_screen`.
- **Score bonuses now implemented** (were defined in config but never triggered):
  - `refseq_preferred` (+10): accession starts with NC_, NM_, NR_, NZ_, NG_,
    XM_, XR_, WP_, XP_, NP_, AC_, NW_, or NT_.
  - `voucher_keyword` (+5): header contains "voucher", "type strain", "holotype",
    "paratype", "type specimen", "neotype", "lectotype", or "syntype".
  - `complete_organelle` (+5): header contains "complete genome", "complete
    mitochondrial", "complete plastid", "complete sequence", "complete cds", or
    "complete coding sequence".
- **Per-species coverage audit in HTML report**: new table showing KEEP / REVIEW /
  REJECT counts per identified species, sorted so species with zero KEEP records
  appear first and are highlighted in amber.  Also written to `summary.tsv` under
  the `species_coverage` section.
- **23 new pytest tests** (115 → 138 total): per-class min length, gzipped FASTA,
  sole representative rescue (6 cases), score bonus reasons (8 cases), and
  uchime output parsing (6 cases).

### Changed

- `n_fraction_high` is now configurable as either a hard reject (default, backwards
  compatible) or a soft score penalty by removing it from `decision_rules.
  hard_reject_reasons` in a custom config.  The aDNA profile uses −20 penalty only.
- `chimera_detected` and `chimera_borderline` added to `hard_reject_reasons` and
  `review_reasons` respectively in DEFAULT_CONFIG.
- `length_below_class_min` added to `hard_reject_reasons` and scoring (−100).
- `sole_representative` (+10) and `chimera_detected` (−100) / `chimera_borderline`
  (−60) added to scoring table.
- Chimera screen step added to pipeline step banner output.
- `rescue_sole_representatives` log section always shown in pipeline output.

---

## v0.5.1 (2026-05-05)

### Fixed

- **Critical: TaxdumpDB cross-kingdom false positives with NCBI 2024+ taxdumps.**
  NCBI reorganised their taxonomy in 2023-2024: the former `superkingdom` rank
  (Bacteria / Eukaryota / Archaea) was renamed to `domain`, and a new `kingdom`
  rank was introduced within Bacteria (e.g. Pseudomonadati, Bacillati) and
  Eukaryota.  `TaxdumpDB.get_kingdom()` previously returned these intermediate
  kingdom names (e.g. "Pseudomonadati") instead of the domain-level name
  ("Bacteria"), causing all bacterial records to be falsely flagged as
  REJECT_CROSS_KINGDOM when compared against the expected kingdom "Bacteria".
  Fixed by adding `get_domain()` which prefers the `domain` rank then falls
  back to `superkingdom` for older taxdumps. `get_kingdom()` is now an alias
  for `get_domain()`.

### Added

- **Verbose pipeline output**: run header lists enabled/disabled steps; each
  BLAST stage prints the database path and thread count; post-annotation stats
  (duplicates, adapters, keywords, N-fraction); per-stage timing; and a
  detailed final summary with decisions by class and kingdom, top reasons,
  all output file paths with descriptions, and "What to do next" hints.
- **`num_threads` support for BLAST**: `taxonomy_blast.num_threads` and
  `windowed_blast.num_threads` config keys pass `-num_threads` to blastn.
  Defaults to 1 (safe); set to 8–32 on a cluster for much faster throughput.
- **New config profiles**:
  - `configs/spinner_with_nt_blast.yml` — full nt + taxdump (general use)
  - `configs/spinner_with_nt_blast_amplicons.yml` — nt + taxdump, windowed
    BLAST disabled (best for short amplicons: NucMark/ITS/18S/COI/cytb)
  - `configs/spinner_CO1_blast.yml` — BOLD+GenBank combined COI database
- **`utils.section()`**: lightweight sub-section header for grouping related
  log lines within a stage.

## v0.5.0 (2026-05-04)

Major refactor from single-file to modular Python package, plus full test suite.

### Added

- **Modular package structure** (`spinner/` directory with 18 focused modules):
  - `utils.py` — terminal helpers, progress bars
  - `seq_utils.py` — IUPAC, revcomp, GC, entropy, homopolymer
  - `config.py` — `DEFAULT_CONFIG`, `load_config`, `deep_update`
  - `fasta.py` — FASTA parsing and writing
  - `regions.py` — region/marker classification, species-kingdom loading
  - `adapters.py` — adapter scanning (Hamming + revcomp)
  - `keywords.py` — bad-keyword scanning
  - `annotation.py` — `Annotation` dataclass + main annotation loop
  - `external.py` — BLAST runner, generic hook parser
  - `vector_screen.py` — vector BLAST result parser
  - `taxonomy_blast.py` — taxonomy BLAST parser, `TaxdumpDB`, windowed BLAST
  - `clustering.py` — vsearch UC parser + run wrapper
  - `decisions.py` — scoring and decision assignment
  - `capping.py` — species × marker capping
  - `reporting.py` — TSV/HTML summary writing, `report_from_decisions()`
  - `pipeline.py` — main pipeline orchestration
  - `cli.py` — argument parser and subcommand dispatch
- **115 pytest tests** covering all core functionality (no external tools required).
- **`windowed_status` field** added to `Annotation` and `decisions.tsv`.
- **`TaxdumpDB` class** for NCBI taxdump lineage lookups (optional; requires `nodes.dmp` + `names.dmp`).
- **`report` subcommand**: regenerate summary TSV + HTML from an existing `decisions.tsv`.
- **`--keep-temp` flag**: prevent cleanup of temporary keyed and windowed FASTAs.
- **Reproducibility outputs per run**:
  - `OUTPREFIX.command.txt` — exact command line
  - `OUTPREFIX.run_config.resolved.yml` — full resolved config snapshot
- **Improved HTML summary report**: styled static HTML with decision stat boxes, per-class breakdown, adapter hit summary, bad keyword summary, taxonomy status table, windowed BLAST status table.
- **Improved example data** (`examples/minimal_input/example_refs.fasta`):
  - Added clean nuclear marker (18S rRNA)
  - Added duplicate-sequence record under different accession (SEQDUPB.1)
  - Added `UNVERIFIED:` bad-keyword record
  - All sequences are distinct (no accidental hash collisions)
- **`explain` handles duplicate accessions**: shows all matching rows with a note suggesting `record_key` for disambiguation.
- `run_config.keep_temp_files` config key to control temp file cleanup.

### Changed

- `Spinner.py` is now a thin launcher (`from spinner.cli import main`).
- `Annotation.as_dict()` now includes `windowed_status` column.
- `parse_tax_blast()` now accepts an optional `TaxdumpDB` for lineage-aware comparison.
- `parse_windowed_blast()` now sets `windowed_status` (was previously only adding reasons).
- HTML report is static and dependency-free (no external CSS/JS libraries).
- `run_pipeline()` creates output parent directories automatically.

### Fixed

- `score_decide()` called twice (before and after capping) to incorporate `cap_exceeded` reason correctly.
- Taxonomy BLAST parser now correctly filters on `min_qcov` threshold.
- `parse_vector_blast()` converts 1-based BLAST qstart/qend to 0-based position for terminal check.

---

## v0.4.0-alpha

This was a broad scaffold release intended for hardening.

- Single-file `Spinner.py` implementation (~654 lines).
- Core FASTA parsing, QC, adapter/keyword screening, decisions TSV, HTML summary.
- External tool hooks for vector BLAST, taxonomy BLAST, windowed BLAST, vsearch, FCS-adaptor, FCS-GX.
- Dedicated `screen-taxonomy`, `screen-vector`, `screen-windowed`, `cluster` subcommands.
- Keyed FASTA writing to handle duplicate accessions in BLAST/clustering.
- Expanded YAML config profiles: default, bait-design, mapping-db, high-confidence.

---

## v0.1.0 (initial release)

- TackleBox: Spinner proof-of-concept.
- Basic FASTA ingest, species/kingdom guessing, marker classification.
- Length / N / entropy / homopolymer / non-IUPAC checks.
- Duplicate accession and exact duplicate sequence detection.
- Adapter scanning (forward + reverse complement, Hamming mismatches).
- Bad keyword scanning.
- Score-based KEEP/REVIEW/REJECT decisions.
- Species × marker capping.
- `decisions.tsv`, `summary.tsv`, `summary.html` outputs.
- `keep.fasta`, `review.fasta`, `reject.fasta` outputs.
- `explain` and `init-config` subcommands.
