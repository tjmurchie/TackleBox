# Spinner Algorithm Notes

Technical notes on how Spinner's internal algorithms work.

---

## Duplicate detection

### Accession duplicates

Records are assigned a unique `record_key` during parsing.  The first record
with accession X gets `record_key = X`; the second gets `X__dup2`; the third
`X__dup3`, and so on.  All records appear in `decisions.tsv` with their own
row.  The second and later copies receive `duplicate_accession = True` and
reason `duplicate_accession` → REJECT.

### Exact sequence duplicates

After IUPAC normalisation (upper-case), a SHA-256 hash is computed for each
sequence.  If a hash is seen more than once, the second and later copies
receive `duplicate_sequence = True` and reason `duplicate_sequence` → REJECT.
This is independent of the accession: two records with different accessions
but identical sequences both get flagged.

---

## Adapter scanning

1. For each adapter in the TSV:
   a. Attempt an exact `str.find()` of the adapter sequence on the forward strand.
   b. If `scan_reverse_complement = true`, repeat for the reverse complement.
   c. If no exact hit: slide a window of `len(adapter)` across the query and
      apply Hamming distance ≤ `max_mismatch`.  Stop at first hit.
2. Record the 0-based position and strand of the first hit.
3. Classify as terminal if `pos ≤ terminal_window_bp` or
   `pos + len(adapter) ≥ len(seq) - terminal_window_bp`.
4. Only the first adapter hit per record is reported in `decisions.tsv`.
   If multiple adapters match, the one that appears first in the TSV wins.

**Terminal vs internal**: terminal hits are treated as REVIEW by default
because they could be insert-size artefacts where the sequencer read into
the adapter at the end of a short fragment.  Internal hits (adapter buried
in the middle of the sequence) are treated as REJECT because they imply
the record is composed partly of adapter sequence, which is dangerous for
bait design.

---

## Adapter Hamming search complexity

The sliding-window Hamming search is O(len(seq) × len(adapter)) per adapter,
per strand.  For typical NCBI records (<10 kb) and a small adapter list
(<20 adapters), this is fast.  For very large FASTAs or many adapters,
consider pre-filtering with an exact match step or reducing `max_mismatch` to 0.

---

## Shannon entropy

Only ACGT bases are counted (N and IUPAC ambiguities are excluded).

    H = -sum(p_i * log2(p_i)) for i in {A, C, G, T}

Values:
- 0.0 = all same base (e.g. polyA tail)
- 2.0 = perfectly even ACGT distribution
- A typical coding sequence ≈ 1.8–2.0
- A telomeric repeat (e.g. TTAGGG) ≈ 1.5

The default `min_shannon_entropy = 1.2` flags extreme low-complexity sequences.
Adjust for your marker type — ITS sequences can have lower complexity than
coding genes.

---

## Homopolymer detection

A simple linear scan counts consecutive identical bases (ACGT only, after
upper-casing).  The maximum run length is reported.  Runs of N are not counted
against the homopolymer threshold (they are counted against `max_n_fraction`).

---

## Scoring and decisions

Each record starts at score 100.  Each reason in `decision_rules.scoring`
adds (positive or negative) to the score.

Decision assignment (after capping):

1. If any reason is in `hard_reject_reasons` → **REJECT** (regardless of score).
2. Else if score ≥ `keep_min` (default 65) AND no reason is in `review_reasons` → **KEEP**.
3. Else if score ≥ `review_min` (default 30) → **REVIEW**.
4. Else → **REJECT**.

`score_decide()` is called twice: once before capping to set initial decisions
(used by `cap_refs()` to skip already-rejected records), and once after capping
to incorporate `cap_exceeded` reasons.

### n_fraction_high treatment

`n_fraction_high` is **not** in `hard_reject_reasons`.  High-N sequences from
rare or extinct taxa receive a score penalty of −20 only.  A high-N reference
from a woolly mammoth is preferable to no reference at all.  To restore hard-
reject behaviour, add `n_fraction_high` to `decision_rules.hard_reject_reasons`
in your config.

---

## Capping

1. Records are grouped by `(species_guess or genus_guess, marker_class)`.
2. Within each group, records are sorted by:
   - KEEP before REVIEW
   - Higher `decision_score` first
   - Longer sequence first
   - Lower `n_fraction` first
   - Accession alphabetically (tie-break)
3. The top N records (where N = `max_per_species_marker[marker_class]`) are kept.
4. The rest receive `cap_exceeded` reason → REVIEW (or REJECT if `cap_action = reject`).
5. Already-REJECT records are excluded from group counting.

---

## Taxonomy BLAST sanity check

### String-matching mode (default)

The top BLAST hit's scientific name is compared to the expected species/genus:

1. If `species_guess` appears (case-insensitive) in the scientific name → `PASS_SPECIES`.
2. Else if `genus_guess` appears → `PASS_GENUS`.
3. Else → `NO_EXPECTED_MATCH`.

This is crude but fast and requires no external data.

### Taxdump lineage mode (when `taxdump_dir` is set)

1. The `staxids` column of the BLAST output is parsed.
2. The first taxid is resolved to its full lineage using `TaxdumpDB`.
3. The lineage kingdom is compared to the expected kingdom from `species_kingdom.tsv`.
4. If kingdoms differ → `REJECT_CROSS_KINGDOM`.
5. Otherwise, genus-level comparison using lineage data.

Loading NCBI taxdump takes ~30 seconds and ~2 GB RAM for the full database.
For most projects, string-matching mode is sufficient.

---

## RefSeq vs nt vs custom DB

**RefSeq** is a curated, non-redundant database with verified annotations.  It
is an excellent conservative sanity check: if a record has no RefSeq match for
its expected genus, something may be wrong.  The limitation is coverage: many
regional or non-model taxa are absent or underrepresented.  Do not
automatically reject records lacking a RefSeq match — use REVIEW instead.

**nt** (NCBI nucleotide) has much broader coverage but includes all public
submissions, including the problematic sequences you may be trying to detect.
Use nt when you need sensitivity for sparse taxa.

**Custom curated DB**: build using `makeblastdb -dbtype nucl -in curated.fasta`.
Best for targeted regional pipelines where you know the expected taxon set.

---

## Windowed BLAST / chimerism screen

1. Only records longer than `enabled_for_min_length` are processed.
2. The sequence is divided into overlapping windows of `window_size` bp with
   `window_step` bp between starts.  Window IDs follow the pattern
   `RECORD_ID|winN|start-end`.
3. Windows are BLASTed against the configured database.
4. The best hit per window is recorded.
5. If ≥ `min_conflicting_windows` windows map to different taxa →
   `WINDOWED_CONFLICT` + reason `windowed_blast_conflict`.

### Conflict detection modes

**String-matching mode** (no taxdump): the first two words (genus + species)
of each window's top-hit scientific name are compared.

**Lineage mode** (taxdump loaded): conflict detection uses
`taxdump_comparison_rank` (default: `family`).  The taxid from the BLAST output
is resolved to its lineage in `TaxdumpDB`, and the name at the configured rank
is compared across windows.  This correctly distinguishes, e.g., *Rangifer* and
*Bos* (both Cervidae/Bovidae, both valid hits in the same region) from a true
chimera mixing animal and plant sequences.  Falls back to genus-string matching
if a window's taxid is absent from the taxdump.

The loaded `TaxdumpDB` from the taxonomy_blast stage is passed automatically to
`parse_windowed_blast()` — no extra configuration is needed.

---

## vsearch UC file format

Spinner's `parse_uc()` reads vsearch UC files with columns (tab-separated):

```
type   cluster_num   length/pident   ...   query_label   target_label
```

- Type `S` (seed): query_label is the centroid.
- Type `H` (hit): query_label is a non-centroid member; target_label is the centroid.
- Type `C` (cluster summary): ignored by Spinner.

Centroids receive `cluster_role = "centroid"` and reason `cluster_representative`.
Members receive `cluster_role = "member"` and reason `cluster_nonrepresentative`.

The keyed FASTA (written as `OUTPREFIX.tmp.keyed.fasta`) uses Spinner's internal
`record_key` as the FASTA ID to avoid ambiguity when duplicate accessions exist.

---

## FCS hook generic TSV parser

The `parse_generic_hook_table()` function reads any TSV with an identifier
column (`accession`, `qseqid`, `query`, or `record_key`) and one or more
status/action columns.  If any cell contains the strings `reject`,
`contaminant`, `adapter`, `vector`, or `foreign`, the hard reason is added;
otherwise the review reason is added.

This flexible design lets Spinner interoperate with any contamination tool
whose output can be post-processed into a simple two-column TSV.

---

## What Spinner checks — and what it does not

### Checks currently implemented

| Category | Check | Reason(s) added |
|----------|-------|-----------------|
| Length | Too short / too long / below class minimum | `length_below_min`, `length_above_max`, `length_below_class_min` |
| Sequence quality | High N-fraction (soft penalty) | `n_fraction_high` |
| Sequence quality | Non-IUPAC characters | `non_iupac_fraction_high` |
| Sequence quality | Low Shannon entropy (mono-nucleotide bias) | `low_complexity` |
| Sequence quality | Long homopolymer run | `homopolymer_long` |
| Redundancy | Exact duplicate sequences (SHA-256) | `duplicate_sequence` |
| Redundancy | Duplicate accession numbers | `duplicate_accession` |
| Contamination | Adapter contamination (Hamming scan, internal/terminal) | `adapter_internal`, `adapter_terminal` |
| Contamination | Vector contamination (UniVec BLAST) | `vector_internal`, `vector_terminal` |
| Contamination | Chimeric sequences (vsearch uchime) | `chimera_detected`, `chimera_borderline` |
| Contamination | Windowed BLAST conflict (long sequences) | `windowed_blast_conflict` |
| Contamination | FCS-adaptor / FCS-GX (hook-based) | `fcs_adaptor_hit`, `fcs_gx_contaminant` |
| Taxonomy | BLAST/MMSeqs2 top hit vs expected taxon | `taxonomy_same_species`, `taxonomy_same_genus`, `taxonomy_no_expected_match`, `taxonomy_cross_kingdom` |
| Header | Bad keyword scan | `bad_keyword_reject`, `bad_keyword_review` |
| Scoring | RefSeq accession prefix bonus | `refseq_preferred` |
| Scoring | Voucher / type specimen keyword bonus | `voucher_keyword` |
| Scoring | Complete organelle keyword bonus | `complete_organelle` |
| Rescue | Sole representative promotion | `sole_representative` |

### Known gaps — categories of bad references that can pass current filters

**Biological plausibility of marker × kingdom**: Spinner does not cross-check
whether the assigned marker class is biologically possible for the assigned
kingdom.  Bacteria classified as `Mito` (bacteria have no mitochondria), or
animals classified as `Plastid` (animals have no chloroplasts), are currently
not flagged.  These combinations arise from header mis-annotation.  A
`marker_kingdom_inconsistent` check is planned (see `docs/roadmap.md`).

**GC content extremes**: For bait probe design, sequences with GC < 20% or
> 80% produce poor hybridisation kinetics.  Spinner reports `gc_fraction` in
`decisions.tsv` but does not filter on it.  A `min_gc_fraction` /
`max_gc_fraction` threshold is planned.

**Terminal vs internal N distribution**: High N-fraction sequences receive a
soft score penalty regardless of whether the Ns are at the ends (trimming
artefact, often still usable for baits) or distributed internally (assembly
gaps, disrupts bait hybridisation).  Planned: separate `n_terminal_fraction`
vs `n_internal_fraction` tracking.

**IUPAC ambiguity codes (R, Y, S, W, K, M, B, D, H, V)**: Spinner rejects
characters outside the full IUPAC DNA alphabet, but does not separately track
the fraction of ambiguous-but-valid IUPAC codes.  These can cause ambiguous
alignments in BWA/Bowtie2.  Planned: `max_iupac_ambiguity_fraction` threshold.

**Length implausibility for marker class**: A 200 bp record labelled "complete
mitochondrial genome" passes all current filters.  Planned: `expected_length_by_class`
ranges to flag implausible combinations.

**NUMTs (nuclear mitochondrial DNA inserts)**: Segments of mtDNA incorporated
into the nuclear genome that are accidentally submitted as mitochondrial
references.  No simple sequence-only check can reliably detect these without a
well-annotated reference nuclear genome for the taxon.

**Tandem repeats**: Short-period tandem repeats (microsatellites, satellite
DNA) that reduce bait specificity are not specifically detected.  Extreme
cases are caught by the Shannon entropy check, but moderate repeats pass
through.

For full details on each gap and planned mitigations, see `docs/roadmap.md`.
