# Changelog

## FlyForge v1.3.0 / FlyForgeAudit v1.3.0 - 2026-09-01 (real default change: --skip-self-mask now defaults to True)

**Real finding, from a live ~18,500-reference multi-species PalaeoSCOPE Phase B panel**:
`self_repeat_softmask()` counts k-mers ACROSS THE ENTIRE input FASTA (every reference
from every sequence combined), not per-reference. This is exactly right for the
small/few-reference case this was originally built for (detecting genuine within-genome
structural repeats -- tandem repeats, transposons, low-complexity runs, before tiling
across them). But at real multi-species panel scale it stops detecting genuine repeats
and instead flags ordinary, biologically CONSERVED single-copy coding sequence -- shared
across species precisely because it's a good marker, not because it's repetitive -- as
"repetitive". Confirmed on a real 18,541-reference panel: 93.7% of all bases masked
overall, and a real 651bp COI barcode gene (the single most standard animal
DNA-barcoding marker in existence) masked at 100%. Hybridization capture actively WANTS
conserved loci -- that's what lets one bait cross-capture DNA from related species --
so masking them out works directly against the goal of a small panel with broad
cross-capture reach.

**`--skip-self-mask` now defaults to `True`** (masking skipped) instead of `False`.
`--no-skip-self-mask` restores the exact original masked-by-default behavior for anyone
who still wants it -- a genuinely small/few-genome reference set is exactly where
within-genome repeat detection remains the right call. Implemented via
`argparse.BooleanOptionalAction`, so `--skip-self-mask` passed explicitly (as existing
callers/scripts already do) behaves identically to before -- only the *bare* invocation
(neither flag given) changes.

**Real before/after comparison** on a genuine 200-record subset of a real curated
multi-species ANIMAL reference panel (32 of which are real COI barcode genes): masked
(old default) produced 1,844 final baits, 73.5% of references got zero final baits, and
29 of 32 COI records got zero baits (20 total COI-derived baits). Skipped (new default)
produced 11,491 final baits (6.2x more), 24.5% zero-bait references, and only 8 of 32
COI records got zero baits (584 total COI-derived baits, ~29x more). Full real ANIMAL
panel numbers and reasoning: see PalaeoSCOPE's own CHANGELOG.md, 2026-09-01 entries.

**FlyForgeAudit got the identical fix**, since it reuses FlyForge's own
`self_repeat_softmask()` but had its own independent `--skip-self-mask` flag/default (in
`add_shared_filtering_args()`, shared by its `audit` and `augment` subcommands) -- fixed
for consistency, same `BooleanOptionalAction` approach, same new default.

`FlyForge.py`'s `main()` argument-parser construction extracted into a new
`build_arg_parser()` function so real CLI defaults/combinations are directly testable,
not just reachable through a full pipeline run (`add_shared_filtering_args()` in
FlyForgeAudit.py was already separately testable, no refactor needed there). New
`tests/test_self_repeat_mask.py` (13 tests total, both tools' first-ever test coverage
of this function/flag): the new CLI default for both FlyForge and FlyForgeAudit,
backward-compat of the explicit `--skip-self-mask` flag, `--no-skip-self-mask` opt-out,
and `self_repeat_softmask()`'s own pure-function correctness (threshold behavior,
N-handling, case-insensitivity, and a direct demonstration of the real
cross-species-conservation finding this change is about).

Also fixed in this release: `pyproject.toml` said `1.2.1` but FlyForge.py's in-file
`__version__` constant was still `1.2.0`, and FlyForgeAudit.py's own `__version__` was
still `1.2.0` too (both drifted during/since the prior release) -- all three now read
`1.3.0`.

46 passed (was 33), 9 deselected (slow, unaffected -- also re-run directly and confirmed
green: `tests/test_design_opool.py`, real BLAST+/primer3, 9 passed in ~5 min). Ruff:
same 2 pre-existing findings, unchanged.

## v1.2.1 - 2026-08-31 (fix: --circular-ids crash when --remove-complements fragments a target)

Real bug, found via a live ~10,000-species PalaeoSCOPE Phase B run: `--remove-complements`
can chop a single input record into 2+ linear sub-fragments (each renamed with a
`_[start:end]` suffix) when it hits its own self-BLAST minus-strand region. If that same
accession was also requested via `--circular-ids`, `parse_circular_id_set()` resolving the
bare accession against the now-fragmented ID set became genuinely ambiguous (it matches
multiple fragments) and this used to raise `RuntimeError`, aborting the ENTIRE FlyForge
run and losing every other reference's already-completed tiling work along with it (real
proof case: `NC_011818.1`, an Accipiter gentilis mitochondrion).

Once a circular molecule has been cut into 2+ linear sub-fragments at arbitrary
internal-repeat boundaries, none of the resulting pieces are topologically circular
anymore anyway (the original end-wraps-to-start join point isn't preserved by any single
fragment) -- so the scientifically correct outcome is to drop that hint from the circular
set entirely, not guess which fragment to keep treating as circular. `parse_circular_id_set()`
now logs a warning and skips an unresolvable/ambiguous hint (that reference is simply
tiled linearly instead) rather than raising; gained an optional `log_fn` parameter (the
main pipeline now passes its own `log`).

New test suite: `tests/test_circular_ids.py`. 33 passed (was 25), 9 deselected (slow).

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
