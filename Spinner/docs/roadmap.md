# Spinner Roadmap

This document describes known limitations, planned improvements, and ideas for
future development.

---

## Current limitations — what Spinner does NOT check

Understanding what is and is not currently detected helps you decide whether
additional manual curation is needed for your project.

### 1. Biological plausibility of marker × kingdom combinations

Spinner classifies each record into a marker class (`Mito`, `Plastid`,
`NucMark`, `Other`) based on header keywords, and assigns a kingdom from the
species-kingdom TSV or heuristic guessing.  It does not currently cross-check
whether these two assignments are biologically consistent:

| Combination | Why it is suspicious |
|-------------|----------------------|
| Bacteria + `Mito` | Bacteria have no mitochondria |
| Bacteria + `Plastid` | Bacteria have no plastids |
| Archaea + `Mito` or `Plastid` | Same; these organelles are eukaryotic |
| Animal + `Plastid` | Animals have no chloroplasts |
| Fungal + `Plastid` | Fungi have no plastids |

These combinations can arise from mis-annotated NCBI submissions, copy-paste
errors in the submission header, or errors in the species-kingdom lookup.
Currently Spinner would assign a RefSeq bonus or pass taxonomy BLAST for such a
record if the sequence itself happens to BLAST to something in the same kingdom.
A `marker_kingdom_inconsistent` reason with a REVIEW action would catch this
class of error without any external tools.

**Planned**: add `basic_qc.check_marker_kingdom_consistency: true` with a configurable
table of forbidden combinations.

---

### 2. GC content extremes

Spinner reports `gc_fraction` in `decisions.tsv` but does not filter on it.
For **bait probe design**, sequences with GC < 20% or GC > 80% produce poor
hybridisation baits: low-GC probes have low melting temperatures; high-GC
probes form stable secondary structures and hairpins.  For **mapping references**,
extreme-GC regions suffer from PCR amplification bias, so references in those
regions may be systematically under-represented in your reads regardless of
their taxonomic correctness.

**Planned**: `basic_qc.min_gc_fraction` and `basic_qc.max_gc_fraction` (defaults
off; suggested range 0.20–0.80 for bait design).  Adds reason
`gc_fraction_extreme` → REVIEW.

---

### 3. Terminal vs internal N distribution

The current `n_fraction_high` check treats all N bases equally.  For bait
design, the location of Ns matters:

- **Terminal Ns** (first/last ~50 bp): common assembly trimming artefact; the
  usable internal sequence may still make a good bait.
- **Internal Ns**: represent true assembly gaps; an internal N in the middle of
  a bait probe disrupts hybridisation and should be flagged more strongly.
- **N-runs**: a single 30-N gap is less problematic than 30 scattered Ns (which
  suggests many independent low-coverage positions).

**Planned**: separate tracking of `n_terminal_count` vs `n_internal_count`, with
configurable thresholds.  An internal gap run longer than `max_internal_n_run`
would add reason `internal_n_gap` → REVIEW.

---

### 4. IUPAC ambiguity fraction (beyond non-IUPAC)

The `non_iupac_fraction_high` check rejects sequences with characters outside
the IUPAC DNA alphabet.  However, IUPAC ambiguity codes within the accepted set
(R, Y, S, W, K, M, B, D, H, V) are not currently tracked separately.

For **short-read mapping**, aligners (BWA, Bowtie2) treat IUPAC codes as
wildcard positions, which can create false-positive alignments at ambiguous
sites.  A reference with 5% IUPAC ambiguity codes distributed throughout the
sequence will produce systematically noisier pileups than an N-heavy sequence
(which aligners typically soft-mask).

**Planned**: `basic_qc.max_iupac_ambiguity_fraction` (separate from N-fraction),
adding reason `iupac_ambiguity_high` → REVIEW.

---

### 5. Sequence length plausibility by marker

`min_length_by_class` sets a floor for each marker class.  But the upper bound
of plausibility is not enforced semantically.  A 200 bp sequence labelled
"complete mitochondrial genome" in the header is impossible (animal mtDNA is
~16 kb), yet it would pass all current filters and receive the
`complete_organelle` bonus if the header contains "complete".

**Planned**: a `basic_qc.expected_length_by_class` dict with `(min, max)` ranges
representing biologically plausible lengths for each marker class.  Sequences
outside the range receive `length_implausible_for_class` → REVIEW rather than
REJECT (the sequence may still be useful, just suspicious).

Example defaults:
```yaml
expected_length_by_class:
  Mito:    [1000, 200000]   # complete mtDNA to small fragment
  Plastid: [1000, 200000]   # complete plastome to fragment
  NucMark: [100, 10000]     # ITS, 18S, 28S fragments
```

---

### 6. Nuclear mitochondrial DNA (NUMTs)

NUMTs are segments of mitochondrial DNA that have been incorporated into the
nuclear genome.  In some taxa (particularly birds and some insects) they are
numerous and highly similar to genuine mtDNA.  A NUMT submitted to GenBank as
"mitochondrial" would pass all Spinner checks including taxonomy BLAST (because
the sequence is genuinely from the expected taxon), yet it produces incorrect
mapping pileups and false metagenomic classifications.

Detection requires either:
- Comparison to a high-quality nuclear genome reference (FCS-GX can help with
  contaminant detection generally)
- Phylogenetic placement within a curated mtDNA reference tree
- Long-read sequencing context (NUMTs are nuclear so they are flanked by nuclear
  sequence)

This is a known hard problem.  **No simple Spinner check can reliably detect
NUMTs without a well-annotated reference genome for the taxon in question.**
Noted here so users are aware.

---

### 7. Repetitive sequences (for bait design)

Tandem repeats, microsatellites, and transposable element-derived sequences
produce non-specific bait probes that capture off-target DNA.  The current
`low_complexity` check (Shannon entropy) catches extreme mono-nucleotide runs
but does not detect tandem repeats (e.g., `ATCATCATCATC...`) or short inverted
repeats that produce hairpin structures.

**Planned (Tier 5)**: optional integration with a simple built-in tandem repeat
scanner (detecting runs of a short period ≤ 10 bp) that adds reason
`tandem_repeat_detected` → REVIEW.  Full RepeatMasker integration is out of
scope.

---

### 8. Sequences from mixed-source extracts

Ancient DNA metagenomics databases sometimes include sequences assembled from
environmental or mixed-extract samples where the taxonomic assignment is
probabilistic rather than definitive.  The `metagenome` and `environmental
sample` keywords were deliberately removed from the default bad keyword list
(v0.7.0) because they are common on legitimate ancient DNA reference sequences.

The trade-off: a small number of genuinely ambiguous sequences may be retained.
If your project requires strict single-source provenance, re-add these keywords
to your `bad_keywords.tsv`:

```
metagenome           review    Metagenomic origin; provenance uncertain
environmental sample review    Environmental source; taxonomy may be approximate
```

---

## Tier 5 — Planned features

These are well-defined improvements that require implementation work but have no
fundamental design barriers.

### Infrastructure

**pip-installable package** (`pyproject.toml`): Remove the dependency on running
Spinner from a fixed directory.  Add an entry point so `spinner` becomes a
system command after `pip install .`.  This is the single most important change
for making TackleBox usable outside its home directory.

**GitHub Actions CI**: A `.github/workflows/test.yml` that runs the pytest suite
on every push and pull request.  Catches regressions immediately.

### New QC checks

**GC content filter** (`basic_qc.min_gc_fraction` / `max_gc_fraction`):
Filter or flag sequences with extreme GC content.  Most useful for bait design
(default suggested: 0.20–0.80).  GC fraction is already computed and stored in
`Annotation.gc_fraction` — only the threshold check and reason need adding.

**Kingdom × marker class consistency** (`basic_qc.check_marker_kingdom_consistency`):
Flag biologically impossible combinations (Bacteria + Mito, Animal + Plastid,
etc.) as `marker_kingdom_inconsistent` → REVIEW.  No external tools needed.

**Internal N-gap tracking**: Separate terminal vs internal N fractions, and flag
long internal N runs (`max_internal_n_run`) that would disrupt bait hybridisation
or create alignment artefacts.

**IUPAC ambiguity fraction**: Track R/Y/S/W/K/M/B/D/H/V separately from N,
add `max_iupac_ambiguity_fraction` threshold.

**Sequence length plausibility** (`expected_length_by_class`): Flag sequences
that are implausibly short or long for their classified marker type.

### New pipeline steps

**Diamond support** (`taxonomy_blast.method: diamond`): Diamond is a fast
translated-search aligner used for aDNA metagenomic classification.  Adding it
as a third method option (alongside `blastn` and `mmseqs2`) would be a natural
extension.  The output parser would need a new column mapping but the pipeline
dispatch pattern is already in place.

**Kraken2 / Bracken output hook**: A post-processing step that reads a Kraken2
report and flags records whose Kraken2 classification conflicts with the expected
taxon.  This would give a second independent taxonomy check without requiring
a BLAST database.

**Tandem repeat scanning**: A built-in lightweight scanner for short-period
tandem repeats (period ≤ 10 bp, run length ≥ 50 bp).  Adds reason
`tandem_repeat_detected` → REVIEW.  Relevant for bait design.

### Output improvements

**Per-class FASTA output**: Write one `keep.CLASS.fasta` per marker class
automatically (e.g., `keep.Mito.fasta`, `keep.NucMark.fasta`).  Useful when
downstream pipelines differ by marker type.

**Richer HTML report**: Add per-marker bar charts, N-fraction histogram, length
distribution, GC distribution, and interactive species coverage table with
sorting and filtering.  Currently the report is tables only.

**`--dry-run` flag**: Validate config and print what *would* happen without
writing any files.  Useful before long BLAST runs.

---

## Design principles — what Spinner intentionally does NOT do

These are out-of-scope by design, not by omission:

- **Phylogenetic placement**: Spinner checks taxonomy by header text and BLAST
  top hit.  It does not build trees or place sequences in a reference phylogeny.
  Use Epa-ng or IQTREE for placement.

- **NCBI metadata correction**: Spinner cannot correct errors in NCBI submission
  records.  It flags them.  Correcting NCBI records requires submitting an
  update to GenBank directly.

- **De novo assembly QC**: Spinner curates pre-assembled reference sequences.
  It is not a genome assembler QC tool (use QUAST or CheckM for that).

- **Read-level QC**: Spinner operates on assembled or curated FASTA references.
  It is not a raw read quality trimmer (use fastp, Trim Galore!, or AdapterRemoval
  for raw reads).
