# Changelog

## Unreleased
- `compute_tm()` Tm-calculation failures (empty, all-ambiguous, or invalid-
  IUPAC-code sequence, e.g. a window that's entirely `R`/`Y`/`W`/etc. -- `N`
  runs of 1-10 are already converted to `T` during preprocessing, so this
  specifically affects longer ambiguity runs or other IUPAC codes) are now
  counted across a run and reported as an explicit warning once tiling
  finishes. Filtering behavior is unchanged (these baits are still excluded,
  same as before) -- only visibility changed: previously a Tm-calculation
  failure was indistinguishable from a genuinely low-Tm bait anywhere in the
  output. `compute_tm()`/`tile_sequence()` gained an optional `on_failure`/
  `on_tm_failure` callback parameter (backward compatible, defaults to `None`).
- added FlyForge's first automated test suite (`tests/`, `pytest.ini`): tiling
  window/step arithmetic including circular-genome wraparound, `compute_tm()`'s
  behavior on normal and pathological input plus the new failure-visibility
  callback, and an end-to-end (real BLAST+/primer3, marked
  `slow`) check that `design_opool()`'s final assembled oligo sequence has the
  correct structure and that the probe sequence itself survives unmodified

## FlyForge v1.2.0 / FlyForgeAudit v1.2.0
- added a third `opool` mode to FlyForgeAudit so an existing bare-bait FASTA can be converted directly into an order-ready oligo pool without redesigning the panel
- `opool` writes a copied bait FASTA, the generated oligo-pool FASTA, amplification primers, probe QC table, recommendations, summary, and progress log
- `opool` summaries distinguish actionable flags from informational notes and report compatibility warnings such as internal BspQI/LguI motifs in the input bait set
- added actionable-flags vs informational-notes reporting in FlyForgeAudit terminal summaries and summary TSVs
- printed recommendation text directly in the end-of-run FlyForgeAudit console summary
- added FlyForge-style progress dashboards to FlyForgeAudit audit and augment modes
- added optional circular-reference handling to FlyForge and FlyForgeAudit (`--circular`, `--circular-ids`) so wraparound bait design/analysis can cover linearized sequence ends
- extended circular-reference BLAST validation subjects by `bait_length - 1` to support wraparound mapping during validation/audit

## v1.1.1 / FlyForgeAudit v1.0.1 - 2026-03-12

### Validation / audit math
- fixed inflated coverage depth in repetitive regions by assigning **one primary on-target placement per bait** instead of summing every valid BLAST HSP
- when bait identifiers encode positional provenance (for example `ref|b123|pos281-360` or trusted coordinate-style suffixes), validation now uses that metadata for coverage placement
- restored expected circular-panel behavior for terminal wraparound baits during audit/validation coverage calculations
- `target_probe_pairs.csv`, per-target probe counts, and `num_targets` are now based on unique primary assignments rather than raw HSP counts

### Augment
- corrected the coverage-deficit math used by `augment` so new-bait proposals are driven by primary panel coverage instead of repetitive-region HSP pileups

### Notes
- this specifically fixes the false 8x/175x-style inflation seen when auditing panels against low-complexity mitochondrial control regions

## v1.1.0 - 2026-03-11

### FlyForge
- corrected T7 promoter documentation to 22 nt (`GCTAATACGACTCACTATAGGG`)
- clarified CARPDM wet-lab workflow language (PCR -> LguI/BspQI digest -> T7 transcription)
- made oligo-pool primer selection fail loudly when no valid second primer is found or primer-selection BLAST parsing fails
- switched validation BLAST to `-dust no -soft_masking false` to avoid undercounting valid low-complexity probe hits
- stopped splitting target identifiers on `|` during validation parsing
- refreshed `per_ref_stats.tsv` from final validation coverage rather than pre-filter tiled coverage
- added `contains_lgui_site` to `probe_info.csv`
- fixed skipped self-mask progress handling

### FlyForgeAudit
- added new companion module with `audit` and `augment` modes
- `audit` reproduces FlyForge-style panel validation outputs for an existing bait set
- `augment` designs a minimal spike-in bait set for new targets and emits extra order-ready oligo-pool FASTAs
- added avoid-database screening, recommendation outputs, and merged-panel summaries

### Documentation
- refreshed main README
- added `FlyForgeAudit_README.md`
- updated release copies and software metadata
