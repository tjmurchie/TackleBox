# TackleBox: Spinner

<p align="center">
  <img src="../assets/Spinner.png" alt="TackleBox: Spinner header" width="700">
</p>

# Spinner — Reference Sequence Curation for Ancient DNA

**Spinner** curates NCBI / custom reference FASTA databases before they are used
for bait design, short-read mapping, BLAST databases, or metagenomic classifiers
(MEGAN, Kraken, Bracken).  It is optimised for **ancient DNA shotgun metagenomics
and capture enrichment** across all kingdoms.

**Core rule: never silently discard a reference.** Every input record receives an
audit row in `decisions.tsv` explaining exactly why it was kept, reviewed, or
rejected.  Rejected records are written to `reject.fasta`, not deleted.

Current version: **v0.7.0**

---

## Where Spinner fits

```
GBIF / regional taxon list
      ↓
FlyGuide  (NCBI nucleotide downloader — downloads reference FASTAs per taxon)
      ↓
raw reference FASTAs  (*_refs.fasta, *_species_kingdom.tsv)
      ↓
Spinner   ← you are here
      ↓                              ↓
curated keep.fasta            review.fasta + reject.fasta
+ decisions.tsv               (inspect before discarding)
      ↓
FlyForge bait design  /  mapping reference  /  MEGAN taxonomy database
```

Spinner can also be used standalone on any multi-FASTA file without FlyGuide
outputs.

---

## Installation

```bash
git clone https://github.com/your-org/TackleBox.git
cd TackleBox/Spinner
chmod +x Spinner Spinner.py

pip install pyyaml          # only required dependency
```

**Python ≥ 3.8** required.

### Optional external tools

Only needed when the corresponding pipeline steps are enabled:

| Tool | Step | Notes |
|------|------|-------|
| `blastn` (NCBI BLAST+) | taxonomy sanity check, vector screen, windowed chimerism | Must be in `$PATH` |
| `mmseqs` (MMSeqs2) | taxonomy sanity check (fast alternative to BLAST) | 10–100× faster than BLAST; needs taxonomy-indexed DB |
| `vsearch` | clustering / redundancy reduction, chimera detection | |
| `makeblastdb` | building local BLAST databases | |
| UniVec BLAST DB | vector/adapter contamination screen | Download from NCBI |
| NCBI FCS tools | assembly-grade adapter/contamination check | Hook-based integration |

If an external tool is absent, Spinner warns and skips the step rather than
crashing (controlled by `run.fail_on_missing_external_tool`).

---

## Quickstart

```bash
cd TackleBox/Spinner

./Spinner filter \
  --fasta examples/minimal_input/example_refs.fasta \
  --species-kingdom examples/minimal_input/example_species_kingdom.tsv \
  --outprefix /tmp/example_spinner
```

`--regions-config`, `--adapters`, and `--bad-keywords` are optional — when omitted,
Spinner auto-resolves them to the bundled defaults (`regions_config_example.tsv`,
`adapters_default.tsv`, `bad_keywords_ancient.tsv`).

For full-stack runs with BLAST, chimera detection, and clustering, use `--mode`:

```bash
# Curated reference database (mapping / metabarcoding / MEGAN)
./Spinner filter \
  --fasta refs.fasta \
  --species-kingdom species_kingdom.tsv \
  --mode reference_db \
  --outprefix out/run_name

# Bait panel (input for FlyForge capture enrichment design)
./Spinner filter \
  --fasta refs.fasta \
  --species-kingdom species_kingdom.tsv \
  --mode bait_panel \
  --outprefix out/run_name
```

`--mode` auto-selects the matching bundled config. Providing `--config` explicitly
overrides `--mode`.

This writes:

```
/tmp/example_spinner.decisions.tsv           ← one row per input record
/tmp/example_spinner.summary.tsv             ← decision and reason counts
/tmp/example_spinner.summary.html            ← human-readable report
/tmp/example_spinner.keep.fasta              ← high-confidence references
/tmp/example_spinner.review.fasta            ← records needing manual check
/tmp/example_spinner.reject.fasta            ← flagged records (auditable)
/tmp/example_spinner.command.txt             ← exact command for reproducibility
/tmp/example_spinner.run_config.resolved.yml ← full resolved config snapshot
```

### Explain a single record

```bash
./Spinner explain \
  --decisions /tmp/example_spinner.decisions.tsv \
  --accession NC_000001.1
```

### Re-generate the HTML/TSV report from an existing decisions table

```bash
./Spinner report \
  --decisions /tmp/example_spinner.decisions.tsv \
  --outprefix /tmp/example_spinner_report
```

### Copy starter configs to customise for your project

```bash
./Spinner init-config --outdir my_project_configs
```

---

## Subcommands

| Subcommand | Description |
|------------|-------------|
| `audit` | Annotate and write decisions TSV — no FASTA output. Use to preview decisions before filtering. |
| `filter` | Audit + write keep / review / reject FASTAs. |
| `screen-taxonomy` | Audit with taxonomy BLAST/MMSeqs2 sanity check enabled (requires `--blast-db`). |
| `screen-vector` | Audit with vector/UniVec BLAST screening enabled (requires `--blast-db`). |
| `screen-windowed` | Audit with windowed chimerism BLAST screen enabled (requires `--blast-db`). |
| `cluster` | Audit with vsearch clustering enabled. |
| `explain` | Explain why a specific accession was kept / reviewed / rejected. |
| `report` | Re-generate summary TSV and HTML from an existing decisions TSV. |
| `init-config` | Copy starter YAML and TSV configs to a directory for project-specific customisation. |

Both `audit` and `filter` accept `--mode reference_db` or `--mode bait_panel` to
auto-select the appropriate full-stack config (see [Run modes](#run-modes-reference_db-vs-bait_panel)).

---

## Input files

### FASTA (`--fasta`)

One or more FASTA files, plain or gzipped (`.fasta.gz`).  Multi-line sequences
are handled.  Duplicate accessions are preserved as separate rows in
`decisions.tsv` with unique record keys (`ACC__dup2`, `ACC__dup3`, …).

### Species-kingdom TSV (`--species-kingdom`)

Tab-separated file mapping species to kingdoms. Produced by FlyGuide as
`*_species_kingdom.tsv`. Columns:

```
species    genus    kingdom
```

Kingdom values are normalised to: `Animal`, `Plant`, `Fungi`, `Bacteria`,
`Archaea`, `Protist`, `Unknown`.

### Regions config TSV (`--regions-config`)

Controls marker/region classification. Each row is a regex rule applied to FASTA
headers. Example:

```
region_id    class      enabled_default    regex
MITO         Mito       1                  mitochondri[on]|mitochondrial
PLASTID      Plastid    1                  chloroplast|plastid|plastome
NUCRDNA_18S  NucMark    1                  \b18S\b|small subunit ribosomal
COI          Mito       0                  \bCOI\b|\bcox1\b
```

If no rule matches, built-in heuristics apply.

### Adapters TSV (`--adapters`)

Columns: `name`, `sequence`, `max_mismatch`, `action`. Comment lines beginning
with `#` are ignored.

```
name                         sequence                                   max_mismatch   action
Illumina_TruSeq_Universal    AGATCGGAAGAGCACACGTCTGAACTCCAGTCA          2              reject
Nextera_ME                   CTGTCTCTTATACACATCT                         1              reject
```

The bundled `configs/adapters_default.tsv` covers sequencing platforms whose
adapter-contaminated sequences appear in public reference databases:

| Platform | Sequences included |
|---|---|
| Illumina (TruSeq / Nextera) | TruSeq Universal + Indexed, Nextera ME, P5/P7 partials, short common prefix |
| 454 / Roche (FLX Titanium) | Adapter A + B — common in GenBank submissions from 2005–2014 |
| Ion Torrent (PGM / Proton) | P1 adapter + A adapter |
| BGI / MGI (DNBSEQ / MGISEQ) | Read 1 + Read 2 + short common prefix |
| PacBio | SMRTbell universal hairpin adapter |

### Bad keywords TSV (`--bad-keywords`)

Columns: `keyword`, `action` (`reject`|`review`), `reason`. Matched
case-insensitively by default against FASTA headers.

The bundled `configs/bad_keywords.tsv` is optimised for ancient DNA:
`UNVERIFIED`, `PREDICTED`, and `low quality` are flagged for **review** (not
rejected), since these labels are common on legitimate ancient DNA submissions.
Only `synthetic construct` and `vector` are hard rejects.

---

## Output files

### `OUTPREFIX.decisions.tsv`

The primary audit table — one row per input record.  Key columns:

| Column | Description |
|--------|-------------|
| `accession` | NCBI accession from FASTA header |
| `record_key` | Internal unique ID (`ACC__dup2` for second copy of same accession) |
| `header` | Full original FASTA header |
| `length` | Sequence length in bp |
| `species_guess` | Species guessed from header |
| `kingdom` | Kingdom from species-kingdom TSV (or heuristic) |
| `marker_class` | `Mito` / `Plastid` / `NucMark` / `Other` |
| `n_fraction` | Fraction of N bases |
| `shannon_entropy` | Sequence complexity (0 = mono-nucleotide; 2.0 = even ACGT) |
| `adapter_hit` | True if an adapter was detected |
| `taxonomy_status` | `PASS_SPECIES` / `PASS_GENUS` / `NO_EXPECTED_MATCH` / `NOT_CHECKED` |
| `windowed_status` | `WINDOWED_OK` / `WINDOWED_CONFLICT` / `NOT_CHECKED` |
| `decision_score` | Composite score (starts at 100, adjusted by each reason) |
| `decision` | **KEEP** / **REVIEW** / **REJECT** |
| `reasons` | Semicolon-separated list of reasons that influenced the decision |

### Decision values

- **KEEP**: record passed all configured filters and its score is ≥ `keep_min` (default 65).
- **REVIEW**: record has a soft flag (low complexity, terminal adapter, borderline score).  Inspect before including.
- **REJECT**: record failed a hard rule (duplicate, length, non-IUPAC, bad keyword reject, etc.).

### `OUTPREFIX.keep.fasta` / `.review.fasta` / `.reject.fasta`

Original headers and sequences are preserved exactly.

---

## Configuration

Spinner is controlled by a YAML config (optional — sensible aDNA defaults are
built in) and TSV files.

### All defaults are aDNA-first

As of v0.7.0, `DEFAULT_CONFIG` is calibrated for ancient DNA metagenomics:

| Setting | Default | Reason |
|---------|---------|--------|
| `max_n_fraction` | `0.20` | Ancient assemblies routinely have 10–20% N from low-coverage gap filling |
| `n_fraction_high` in `hard_reject_reasons` | **excluded** | High-N sequences from rare/extinct taxa are better than no reference |
| `cap_references` | `false` | All haplotype diversity is retained; capping reduces capture breadth |
| `rescue_sole_representatives` | `true` | If a species has 0 KEEP records, its best REVIEW is promoted |
| `keep_min` | `65` | Imperfect references are still valuable |
| `review_min` | `30` | Prefer REVIEW over REJECT for borderline records |
| `taxonomy_blast.min_pident` | `70.0` | Divergent ancient sequences tolerate lower identity |
| `review_if_no_expected_match` | `false` | Rare taxa are often absent from NCBI databases |

To override any default, create a minimal YAML with only the keys you want to
change:

```yaml
# my_project.yml — only override what you need
basic_qc:
  min_length: 100
taxonomy_blast:
  blast_db: /data/databases/nt/nt
  taxdump_dir: /data/databases/taxdump
  num_threads: 16
```

### Config profiles (bundled)

| Profile | `--mode` shortcut | Use case |
|---------|---|----------|
| `spinner_reference_db.yml` | `reference_db` | Full-stack curation for read mapping / metabarcoding / MEGAN databases |
| `spinner_bait_panel.yml` | `bait_panel` | Full-stack curation for capture enrichment bait synthesis via FlyForge |
| `spinner_ancient_metagenome.yml` | — | Named aDNA profile (matches defaults; use for record-keeping) |
| `spinner_with_nt_blast.yml` | — | NCBI nt + taxdump + windowed BLAST for long sequences |
| `spinner_with_nt_blast_amplicons.yml` | — | nt + taxdump, windowed BLAST disabled (short amplicons) |
| `spinner_CO1_blast.yml` | — | COI / BOLD+GenBank combined database |
| `spinner_bait_design.yml` | — | Strict QC thresholds, no BLAST (legacy bait profile) |
| `spinner_mapping_db.yml` | — | Inclusive: preserve maximum taxon diversity |
| `spinner_high_confidence.yml` | — | Conservative: publication-grade curation |

> **Note:** `spinner_reference_db.yml`, `spinner_bait_panel.yml`, and configs
> containing local database paths are excluded from version control via `.gitignore`.
> They are kept locally for reproducibility.

### Toggling pipeline steps

```yaml
steps:
  basic_qc: true            # length, N-fraction, complexity, duplicates
  classify_regions: true    # marker class (Mito/Plastid/NucMark/Other)
  adapter_screen: true      # Illumina/Nextera adapter scanning
  bad_keyword_screen: true  # header keyword filtering
  vector_screen: false      # UniVec BLAST — enable for bait-design pipelines
  taxonomy_blast: false     # BLAST/MMSeqs2 sanity check — requires a database
  windowed_blast: false     # chimerism detection for long sequences (≥1000 bp)
  chimera_screen: false     # vsearch uchime chimera detection
  cluster: false            # vsearch clustering to reduce haplotype redundancy
  cap_references: false     # per-species cap (disabled by default for aDNA)
  fcs_adaptor: false        # NCBI FCS-adaptor hook
  fcs_gx: false             # NCBI FCS-GX contamination hook
  report: true
```

---

## Run modes: `reference_db` vs `bait_panel`

Both modes run the full pipeline — taxonomy BLAST, windowed chimerism detection,
vsearch uchime chimera screen, and clustering — on top of the standard adapter,
keyword, and basic QC screens. They differ in how strict each threshold is and
what happens to chimeric sequences.

### Core difference in purpose

**`reference_db`** prioritises coverage and diversity. Every species should have
representation, rare taxa are valuable, and borderline sequences are flagged for
review rather than discarded. A slightly imperfect sequence from a rare plant
species is more useful than no reference.

**`bait_panel`** prioritises physical bait quality. Every output sequence will be
synthesized as a real reagent. A bad bait does not just get ignored — it
physically captures the wrong target, contaminates your enrichment, and corrupts
downstream results. FlyForge performs its own design-level filtering (tiling, GC
optimization, off-target specificity), so Spinner's job for bait input is
specifically to eliminate sequences that would produce bad baits.

### Setting-by-setting comparison

| Setting | `reference_db` | `bait_panel` | Why it differs |
|---|---|---|---|
| `min_length` | 50 bp | **80 bp** | 80 bp is the practical minimum for a stable hybridizable bait |
| `max_n_fraction` | 5% | **1%** | N-masked positions in a bait cannot hybridize — wasted synthesis budget |
| `min_shannon_entropy` | 1.2 | **1.5** | Low-complexity sequences hybridize non-specifically across the genome |
| `max_homopolymer_run` | 60 | **30** | Long homopolymer runs cause synthesis failures |
| `windowed_blast conflict_action` | `review` | **`reject`** | A chimeric reference is auditable; a chimeric bait actively captures wrong targets |
| `cluster identity` | 0.99 | **0.98** | FlyForge handles tiling density; Spinner just eliminates near-duplicates |
| `max_per_species_marker` (Mito/Plastid) | 20 / 20 | **8 / 8** | FlyForge optimises per-target bait count; Spinner feeding it 20 near-identical sequences wastes work |
| `rescue_sole_representatives` | true | true | Both modes retain at least one record per species for coverage |

### Chimerism detection — two complementary screens

Both modes run both chimera screens; they catch different things:

| Screen | Works on | Detects |
|---|---|---|
| `windowed_blast` | Sequences ≥ 1000 bp only | Mixed-origin sequences: discordant taxonomy across windows of the same sequence (e.g. plastid genome with an accidental mitochondrial insert) |
| `chimera_screen` (vsearch uchime) | Any length, especially short amplicons | PCR chimeras: artificial sequences produced when two templates recombine during amplification and the product is submitted to GenBank |

The two screens rarely overlap — windowed BLAST catches assembly-level problems
in long sequences; uchime catches amplification-level problems in short ones.

---

## Taxonomy BLAST / MMSeqs2

### Using BLAST

```yaml
steps:
  taxonomy_blast: true
taxonomy_blast:
  method: blastn
  blast_db: /data/databases/nt/nt
  taxdump_dir: /data/databases/taxdump   # optional: enables lineage-aware kingdom check
  num_threads: 8
  min_pident: 70.0
  min_qcov: 50.0
```

### Using MMSeqs2 (10–100× faster)

MMSeqs2 requires a database created with `mmseqs createdb` and indexed with
`mmseqs createtaxdb` for taxonomy columns to work.

```yaml
steps:
  taxonomy_blast: true
taxonomy_blast:
  method: mmseqs2
  blast_db: /data/databases/mmseqs2/nt_tax
  num_threads: 16
  min_pident: 70.0
```

The output format is parsed identically to BLAST — no other changes needed.

### Database choice

| Database | Pros | Cons | Best for |
|----------|------|------|----------|
| RefSeq | Curated, non-redundant | Incomplete for sparse/regional/ancient taxa | Conservative sanity checks |
| nt | Broadest coverage | Noisier; includes bad records you may be detecting | Sparse taxa, aDNA |
| Custom curated | Tailored to project | Requires curation effort | Regional pipelines |

---

## Windowed BLAST (chimerism detection)

For sequences ≥ 1000 bp, Spinner tiles overlapping 500 bp windows and BLASTs
each independently.  If ≥ 2 windows map to different taxa, the record is flagged
as `WINDOWED_CONFLICT`.

When a taxdump is loaded, conflict detection uses lineage comparison at the rank
configured by `windowed_blast.taxdump_comparison_rank` (default: `family`):

```yaml
steps:
  windowed_blast: true
windowed_blast:
  blast_db: /data/databases/nt/nt   # falls back to taxonomy_blast.blast_db
  taxdump_comparison_rank: family   # genus / family / order / class
  num_threads: 8
```

---

## Chimera detection

Uses vsearch uchime (de novo or reference-based):

```yaml
steps:
  chimera_screen: true
chimera_screen:
  method: uchime_denovo      # no reference DB needed
  # method: uchime_ref       # use reference_db for higher sensitivity
  # reference_db: /path/to/reference.fasta
  reject_chimeras: true
  review_borderline: true
```

---

## Sole representative rescue

After all filtering, if a species would have zero KEEP records, the best
available REVIEW record is promoted to KEEP.  This prevents database holes for
rare, extinct, or poorly-sequenced taxa.

Enabled by default (`capping.rescue_sole_representatives: true`).  To disable:

```yaml
capping:
  rescue_sole_representatives: false
```

---

## Decision scoring

Every record starts with score `100`.  Each reason adds a positive or negative
modifier:

| Reason | Effect |
|--------|--------|
| `adapter_internal` | −100 → REJECT |
| `duplicate_sequence` | −100 → REJECT |
| `bad_keyword_reject` | −100 → REJECT |
| `n_fraction_high` | −20 (soft penalty only) |
| `bad_keyword_review` | −30 |
| `adapter_terminal` | −40 |
| `taxonomy_no_expected_match` | −10 |
| `taxonomy_same_species` | +20 |
| `taxonomy_same_genus` | +10 |
| `refseq_preferred` | +10 (accession starts with NC_, NM_, NR_, etc.) |
| `voucher_keyword` | +5 (header contains "voucher", "holotype", etc.) |
| `complete_organelle` | +5 (header contains "complete genome", etc.) |
| `sole_representative` | +10 (rescued as only reference for species) |

Score ≥ 65 and no review reason → **KEEP**
Score ≥ 30 → **REVIEW**
Score < 30 or any hard-reject reason → **REJECT**

---

## Reproducibility

Every run writes:

- `OUTPREFIX.command.txt` — exact command line
- `OUTPREFIX.run_config.resolved.yml` — full resolved config including defaults

These files allow an exact re-run from any decisions table.  Include them as
supplementary material in publications.

---

## Testing

```bash
pip install pytest pyyaml
cd TackleBox/Spinner
python -m pytest tests/ -v
```

142 tests covering: FASTA parsing, gzipped FASTA, duplicate detection, adapter
internal/terminal logic, bad keyword scanning, region classification, per-class
minimum length, score bonuses (RefSeq, voucher, complete organelle), decision
scoring, capping, sole representative rescue, taxonomy BLAST parsers, windowed
BLAST (with and without taxdb), chimera detection (uchime), and end-to-end
pipeline runs.

---

## Troubleshooting

**`pyyaml` not found**: `pip install pyyaml`

**`blastn` not found**: Install NCBI BLAST+ and add to `$PATH`.  Spinner warns
and skips BLAST steps by default (`run.fail_on_missing_external_tool: false`).

**All records getting `taxonomy_not_checked`**: BLAST step is enabled but no
hits passed the `min_pident`/`min_qcov` thresholds.  Check your database path
and try lowering `min_pident` to `60.0` for a diagnostic run.

**Many `NO_EXPECTED_MATCH`**: Normal for rare/ancient taxa.  With default
settings this does not force REVIEW — check `review_if_no_expected_match: false`.

**Duplicate accessions**: Expected — both records appear in `decisions.tsv`.
The second copy gets `duplicate_accession` → REJECT.  Use `record_key` (not
`accession`) to disambiguate in `explain`.

**`cap_rank: 0`**: Records already REJECT before capping are excluded from cap
counting and retain `cap_rank: 0`.

---

## Roadmap

See [`docs/roadmap.md`](docs/roadmap.md) for the full list of planned
improvements and known limitations.  Highlights:

**Planned QC additions (Tier 5)**
- GC content range filter (`min_gc_fraction` / `max_gc_fraction`) — catches
  bait-design failures due to extreme GC without any external tools
- Kingdom × marker class biological consistency check — flags records where
  the assigned marker is impossible for the assigned kingdom (bacteria + Mito,
  animal + Plastid, etc.)
- Internal vs terminal N distribution — distinguish assembly trimming artefacts
  from internal assembly gaps
- IUPAC ambiguity fraction threshold — R/Y/S/W/K/M codes tracked separately
  from N; relevant for mapping reference quality
- Length plausibility by marker class — flag "complete mitochondrial genome"
  records that are only 200 bp

**Planned pipeline additions**
- Diamond as a third taxonomy search method (`method: diamond`)
- Kraken2 classification cross-check hook
- Per-class FASTA output (`keep.Mito.fasta`, `keep.NucMark.fasta`, etc.)
- Richer HTML report with length/GC/N-fraction distributions
- `pip install .` support (`pyproject.toml`)
- GitHub Actions CI

**Known hard problems** (no simple fix without external reference data)
- NUMT detection (nuclear mitochondrial DNA inserts)
- Tandem repeat / microsatellite content for bait specificity

---

## Citation

Spinner is part of TackleBox.  If you use it in a publication, please:

1. Cite the TackleBox repository and Spinner version (`./Spinner --version`).
2. Include the `OUTPREFIX.run_config.resolved.yml` as supplementary material.
3. Note the BLAST/MMSeqs2 database name and download date.

---

## Development

See `docs/algorithm_notes.md` for implementation details and
`docs/config_reference.md` for a full description of every config key.

Contributions and bug reports are welcome via GitHub Issues.
