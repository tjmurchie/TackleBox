# Spinner Changelog

## 2026-08-30 (follow-up: report/console output now matches the accept/reject-only default)

Same-day follow-up to the decision-model change below. `reporting.py`
(`write_summary_tsv`/`write_summary_html`/`write_split_fastas`) and `pipeline.py`'s
console output (`_print_final_summary`, the run-mode banner, the numbered "what to do
next" hints) now take/derive `three_state_mode` and stop showing REVIEW-related
output entirely in the real default: no more permanently-empty REVIEW stat box/table
column/summary line, and `.review.fasta` is no longer written at all (previously always
written, always empty). `three_state_mode: true` restores every REVIEW-related output
exactly as it was before. `Spinner report` (`report_from_decisions`, which regenerates a
report from an existing decisions.tsv with no original cfg available) auto-detects which
mode to render by checking whether any real REVIEW decision is present in the file.

183 passed (was 182), ruff clean.

## 2026-08-30 (real behavior change: accept/reject only is now the default -- REVIEW retired)

Spinner's decision model changes from three states (KEEP/REVIEW/REJECT) to two
(KEEP/REJECT) by default. Motivation: Fillet's `bait-eval` (a separate, comprehensive
identity/capture-competition check run against the actual finished baits, against the
full NT database) is now the authoritative final check in this pipeline -- Spinner's own
job is cheap, fast triage, not the last word, so it no longer needs a REVIEW tier forcing
manual attention before trusting its own calls.

**What changed** (`decisions.py::score_decide()`): `hard_reject_reasons` still always
REJECT, unchanged. Every other reason in `review_reasons` still applies its existing
score penalty, but no longer blocks KEEP on its own -- a single `score_thresholds.keep_min`
cutoff decides KEEP vs REJECT. Set the new `decision_rules.three_state_mode: true` to
restore the exact old behavior (`score_thresholds.review_min` only has an effect in that
mode).

**Real, deliberate side effect worth knowing about**: `cap_action: "review"` (one of
capping's two modes) is now genuinely soft in the default mode -- a cap-exceeded record
whose score otherwise clears `keep_min` will still reach KEEP (previously, hitting a
review_reason always blocked KEEP regardless of score). Project configs that want capping
to be a hard, enforced limit should use `cap_action: "reject"` instead, which is
unaffected by this change.

**`rescue_sole_representatives()` retargeted, not removed**: in the old REVIEW-based
world, a species with zero KEEP records could have its best REVIEW candidate promoted.
With REVIEW gone by default, this function now promotes the best REJECT candidate that
was NOT hard-rejected (i.e. rejected purely on score, not for a structural reason like
contamination/chimera/duplicate/adapter-vector hit) -- preserving the same "don't lose
the only reference for a rare species" protection this function has always existed for.
`three_state_mode: true` restores the exact old REVIEW-based rescue candidate pool.

182 passed (was 174), ruff clean (4 pre-existing findings, unchanged). Still to come in
this same round: `reporting.py`/`pipeline.py`'s console/file output still reference
REVIEW literally (harmless in the new default -- just always empty/zero -- but not yet
cleaned up for clarity); the planned NR-protein escalation extension and Fillet-backed
reject-audit are separately scoped, not yet built.

## 2026-08-29 (complete organelle genomes were permanently blocked from auto-KEEP)

Found via forensic analysis of a real complete-mitogenome cluster centroid (accession
EU153454.1, 16,480 bp) from a live PalaeoSCOPE run: it scored 110 (well above
`keep_min: 65`), was correctly selected as its 166-member cluster's centroid, and was
clean on every other QC axis, yet stayed stuck at REVIEW instead of reaching KEEP
automatically. Root cause: `taxonomy_blast.max_query_length` (10,000 bp in real project
configs) excludes full organelle genomes from the taxonomy verification search entirely,
per that setting's own long-standing comment ("skip sequences >10 kb — full organelles
don't need kingdom verification"). But length-exempted records used to be tagged
`taxonomy_not_checked` -- the EXACT SAME reason a record genuinely SUBMITTED to the search
but given no hit receives -- and since `taxonomy_not_checked` is unconditionally in
`decision_rules.review_reasons`, every length-exempted complete genome was permanently
blocked from auto-KEEP regardless of score. This directly undermined the goal of an
end-to-end pipeline: the very records most confidently correct (complete organelle
genomes) required manual review every time, while genuinely ambiguous short/partial
records got the identical treatment.

Fixed: length-exempted records with `marker_class` of `Mito` or `Plastid` (the two classes
the exemption's own rationale actually applies to -- a random contaminant essentially
never assembles into a complete, correctly-headered, correctly-sized organelle genome) now
get a new, distinct `taxonomy_exempt_length` reason instead, scored neutrally (`0`) in
`config.py` and deliberately **not** added to `review_reasons`/`hard_reject_reasons`, so a
genuinely clean, high-scoring complete genome can reach KEEP on its own merit. Records of
other marker classes (`NucMark`, `Other`) that happen to also exceed `max_query_length`
(e.g. a large nuclear scaffold) are deliberately left untagged and keep the normal
`taxonomy_not_checked` review-forcing safety net -- they have no equivalent
"can't-plausibly-be-a-contaminant" guarantee. New `taxonomy_blast.py::mark_length_exempt_records()`
implements the marker-class-scoped tagging (called from `pipeline.py`'s taxonomy_blast
stage before `parse_tax_blast()` runs); `parse_tax_blast()`'s per-record loop now skips
records already tagged this way so they are never double-tagged with the blocking reason.

174 passed (was 166), ruff clean (18 pre-existing findings in touched files, unchanged).
Adversarially verified by an independent fresh agent before commit given this changes
live default scoring behavior for every Spinner user, not just PalaeoSCOPE.

Verified against real data: re-ran the exact same `Spinner filter` command against the
real beringia_organelle_test_v2 project's already-fetched FlyGuide output (reusing the
cached `taxonomy_blast.tsv`, so this isolates the fix's effect from the search itself).
Exactly one record's decision changed: EU153454.1, REVIEW -> KEEP, decision_score
unchanged at 110, reasons correctly `complete_organelle;taxonomy_exempt_length;
cluster_representative` (was `...;taxonomy_not_checked;...`). All 700 other REJECT/REVIEW/
KEEP calls are byte-identical to the pre-fix run -- zero unintended promotions.

## 2026-08-28 (a missing/broken taxonomy tool was silently MORE permissive than a working one)

Found via a live PalaeoSCOPE run at real Holarctic-beetle scale (11 species, taxonomically
much more diverse than earlier single-genus tests) whose `mmseqs` binary was missing from
`PATH` in that specific run's environment: `taxonomy_blast.py::parse_tax_blast()` used to
`return` immediately when its input file did not exist (e.g. because the search subprocess
itself failed to run at all, caught and logged as a warning by `pipeline.py`) -- skipping
its own per-record loop entirely, so **no record ever got the `"taxonomy_not_checked"`
reason added**, even though every one of them was genuinely never checked. Since only that
reason string (not the `taxonomy_status` field alone) drives review/reject scoring, a
missing/broken taxonomy tool silently produced a **more permissive** panel than either a
working search or an explicitly-disabled one -- the same review-forcing safety net a
record deliberately exempted by `max_query_length` already correctly gets. Confirmed via a
real accession: a genuine 16kb+ complete-genome record correctly got `taxonomy_not_checked`
(and was correctly barred from auto-KEEP by it) when the search ran and exempted it by
length, but was MISSING that reason entirely -- and reached KEEP -- when the whole search
silently failed due to the missing binary.

Fixed: `parse_tax_blast()` now reads the file only when it exists, leaving `hits` empty
otherwise, so its own per-record loop (which already correctly treats "not present in
`hits`" as `taxonomy_not_checked` for the length-exemption case) now runs unconditionally
and applies the identical reason for the tool-failure case too. `pipeline.py`'s calling
site now always calls `parse_tax_blast()` rather than gating it on the output file's
existence, with a clearer log message ("No taxonomy_blast output available (search failed
to run) -- marking all candidates taxonomy_not_checked") for the failure case.

Verified against real data with a deliberately mmseqs-less `PATH` and a fresh output
prefix (to rule out Spinner's own resume-cache reusing a prior successful run's output):
all 151 real records in the Holarctic beetle test now correctly get `taxonomy_not_checked`
when the tool is missing, and the resulting KEEP/REVIEW/REJECT tally (11/43/97) exactly
matches the tally produced when `mmseqs` is genuinely present and running -- confirming the
fix makes a broken environment behave conservatively (same outcome as a real, working
check finding nothing conclusive) rather than silently more permissively (18 KEEP under the
old, broken behavior).

166 passed (was 164), ruff clean (24 pre-existing findings, unchanged, confirmed via git
stash diff).

## 2026-08-27 (real-data curation bugfixes, found via PalaeoSCOPE's first live e2e run)

Four real, confirmed bugs found via forensic analysis of real `curated_refs.decisions.tsv`
output from an actual `bait_panel`-mode run against live NCBI data (not synthetic test
fixtures) — the first real user of a couple of these code paths surfaced genuine, previously
latent defects. 164 passed (was 142), ruff clean, zero regressions to any existing config's
behavior.

- **`guess_species()` (`regions.py`) fabricated a bogus pseudo-species from legacy
  abbreviated-genus headers** (e.g. `"M.primigenius"`, short for `"Mammuthus primigenius"`,
  smashed into one token) — it mistook the abbreviation for a genus and the next ordinary
  description word (e.g. `"mitochondrial"`) for a species epithet, fabricating
  `"M.primigenius mitochondrial"` as its own singleton species group. Since clustering/
  capping group by `(species_guess, marker_class)`, this let a 93bp unverified fragment
  bypass all real competition and get auto-promoted to KEEP via
  `rescue_sole_representatives()`, ahead of 38+ genuine complete mitogenomes for the same
  real species. Fixed via a new `build_genus_abbrev_index()` pre-pass (built from the same
  input batch's own correctly-parsed full-genus headers) that `guess_species()` can now use
  to resolve an abbreviated header back to its real species; degrades safely (no
  fabrication) when no matching full-genus sibling exists in the batch.
- **`capping.cap_action: "reject"` was dead code** (`capping.py`, `decisions.py`,
  `config.py`) — `cap_refs()` set `a.decision = "REJECT"` for an over-cap record under
  `cap_action: "reject"`, but the mandatory second `score_decide()` call immediately after
  it unconditionally recomputes `a.decision` from `a.reasons` alone, discarding whatever
  `cap_refs()` just set. Since `"cap_exceeded"` was only in `review_reasons`, not
  `hard_reject_reasons`, `cap_action: "reject"` had *zero* real effect versus the default
  `"review"` in every case — confirmed via two real accessions that should have been
  impossible to KEEP under a `max_per_species_marker: 0` + `cap_action: reject` config, but
  reached KEEP anyway via rescue. Fixed by having `cap_refs()` add a new, distinct
  `"cap_exceeded_reject"` reason (kept alongside the existing generic `"cap_exceeded"`)
  specifically under `cap_action: "REJECT"`, now listed in `decision_rules.hard_reject_reasons`
  — every existing config uses `cap_action: "review"` and is completely unaffected
  (regression-tested for byte-identical behavior).
- **Legacy non-coding marker sequences (rRNA/tRNA/D-loop/introns/intergenic spacers) never
  got taxonomically verified** — the only taxonomy check was a 6-frame-translated protein
  search with no ORF to find in non-coding sequence, leaving the large majority of a
  real legacy-marker-heavy dataset at `taxonomy_not_checked` by default. Fixed by wiring in
  Spinner's own existing (but previously unused) nucleotide-mode search machinery: new
  opt-in `taxonomy_blast.nt_fallback_db` config key (default `""`, fully inert unless
  explicitly set) re-runs records left `taxonomy_not_checked` by the primary search through
  a real nucleotide megablast. The intentional `max_query_length` skip of full organelle
  genomes from taxonomy checking (documented in bundled configs as "full organelles don't
  need kingdom verification") is untouched by this change.
- **`"complete_organelle"` reason false-triggered on partial records** (`annotation.py`) —
  `_COMPLETE_TERMS` included the bare phrases `"complete sequence"`/`"complete cds"`/
  `"complete coding sequence"`, which GenBank headers routinely use to describe ONE small
  named sub-feature (a single gene or spacer) within a record that is explicitly partial
  everywhere else. Tightened to 4 unambiguous whole-genome-scale phrases (`"complete
  genome"`, `"complete mitochondrial genome"`, `"complete plastid genome"`, `"complete
  chloroplast genome"`) — calibrated against a real 718-record dataset, where 181 of 229
  old-term matches were confirmed real false positives from sub-feature phrasing.

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
