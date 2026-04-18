# MetaMerge v1.1.0

**MetaMerge** is a command-line tool for ensemble ancient-DNA classification.
It merges a conservative **BLASTn → MEGAN7** taxon count matrix with
**Holi / metaDMG** damage assessments to assign DNA-evidence confidence
categories to each detected taxon.

The core idea: MEGAN7 with strict LCA parameters is deliberately conservative —
it will miss real taxa, but what it does call is trustworthy.  Holi/metaDMG
characterises DNA damage (C→T substitutions) for a broader set of alignments,
providing evidence that a taxon's signal is genuinely ancient rather than modern
contamination.  MetaMerge combines both lines of evidence into a single,
auditable output with explicit confidence levels.

---

## How it works

```
MEGAN7 count matrix  ──┐
                       ├──▶  MetaMerge  ──▶  merged workbook + TSV + reports
Holi/metaDMG CSV     ──┤
Library linker CSV   ──┘
```

The **library linker** is a small crosswalk CSV that maps each MEGAN column
name to its corresponding Holi sample name and carries project metadata
(site, sample type, is_negative_control, …).  MetaMerge uses the linker to
join the two workflows in a library-aware way without requiring identical
column naming conventions.

---

## DNA support categories

| Category | Criteria |
|---|---|
| **Very high confidence** | Exact Holi damage support (≥ 100 reads, strong sig., clean QC, strong MEGAN counts) |
| **High confidence** | Exact Holi damage support (below the ≥ 100 read threshold or with mild QC caution) |
| **Supported** | Exact Holi support or strong MEGAN counts or low-rank lineage-consistent damage support |
| **Tentative** | Weak / incomplete DNA evidence; above the tentative-floor read count |
| **Weak support** | Minimal non-control evidence; below the tentative floor |
| **Blank-associated** | Substantial blank overlap relative to real-library counts |

"Blank-associated" always takes priority regardless of damage evidence.

---

## Installation

```bash
# Clone or copy the repository, then install in editable mode:
pip install -e .

# Verify:
metamerge --help
```

Requires Python ≥ 3.10 and the packages listed in `pyproject.toml`
(`pandas`, `openpyxl`, `PyYAML`, `tqdm`, `requests`, `numpy`).

---

## Workflow

### Step 0 — Generate the library linker (first run only)

If you are starting from a project metadata spreadsheet and MEGAN TSV files,
use the bundled helper script to generate the linker automatically:

```bash
python scripts/make_metamerge_linker.py \
    --meta   project_metadata.xlsx \
    --megan  megan_counts_dir/ \
    --metadmg metaDMG_output.csv \
    --out    library_linker.csv
```

The script reads the MEGAN column headers, strips the long file-path suffixes
to recover library stems, matches them to your metadata by library ID, and
derives Holi sample names by applying a hierarchy of normalizations
(extraction-suffix stripping, underscore→hyphen conversion, leading-zero
removal).  Any rows it cannot resolve automatically are flagged as
`REVIEW_NEEDED` for manual correction.

See `scripts/make_metamerge_linker.py --help` and
[docs/metadata_linker.md](docs/metadata_linker.md) for details.

### Step 1 — Validate inputs

```bash
metamerge check \
    --counts megan_counts.tsv \
    --holi   metaDMG_output.csv \
    --meta   library_linker.csv
```

Prints a pre-run summary: how many libraries match, how many taxa are present
in MEGAN vs Holi, and any mismatches that would cause silent data loss.

### Step 2 — Run the merge

```bash
metamerge run \
    --counts megan_counts.tsv \
    --holi   metaDMG_output.csv \
    --meta   library_linker.csv \
    --config config/defaults.yaml \
    --outdir results/
```

Optional flags:
- `--online-common-names` — query NCBI → GBIF → iNaturalist for English common names
- `--max-rank family` — most inclusive (broadest) rank in heatmap tables (`species`,
  `genus`, `family` *(default)*, `order`, `class`, …)
- `--render-heatmaps` — render PDFs immediately after the merge

Outputs:
- `results/<prefix>_merged_support.xlsx` — publication-ready Excel workbook
- `results/<prefix>_merged_support.tsv` — tab-separated text version
- `results/<prefix>_summary.json` — machine-readable run summary
- `results/report_inputs/` — long-format plotting tables for R heatmaps

### Step 3 — (Optional) Render heatmaps

```bash
metamerge run ... --render-graphs
# or separately:
metamerge report \
    --input-dir results/report_inputs \
    --outdir    results/reports
```

Produces grouped heatmaps, damage heatmaps, and stacked-bar plots (animals,
plants, fungi, microbes) using the bundled R script
`r/metamerge_heatmap_report.R`.

After the run completes, a `sample_order.csv` file is written to
`results/report_inputs/sample_order.csv`.  Edit this file and pass it to
`metamerge report --sample-order` to customise the layout of every plot
(see [Customising plot layout](#customising-plot-layout) below).

### Step 4 — (Optional) Customise plot layout

```bash
# 1. Edit results/report_inputs/sample_order.csv in any spreadsheet tool
# 2. Re-render with the edited file:
metamerge report \
    --input-dir results/report_inputs \
    --outdir    results/reports_custom \
    --sample-order results/report_inputs/sample_order.csv
```

---

## Customising plot layout

Every `metamerge run` writes a wide-format CSV to
`<outdir>/report_inputs/sample_order.csv`.  Each **column** is a plot group
(e.g. `A51_sediments`, `NHB_bones`, `negative_controls`) and each **row**
within a column is a `plot_sample_id` in the order you want it to appear on
the x-axis of every heatmap, damage heatmap, and stacked-bar plot.

### What you can do by editing the file

| Edit | Effect |
|---|---|
| Reorder rows within a column | Changes the left-to-right sample order within that plot group |
| Move a sample ID to a different column | Re-assigns that sample to the new plot group |
| Place the same sample ID in **multiple** columns | Sample appears in each of those group plots simultaneously |
| Delete a row entirely | Excludes that sample from all plots |
| Add a row | Adds a sample to a group (it will appear as an empty column if it has no detections) |

> **Zero-detection samples** are always included in `sample_order.csv` and
> will appear as empty columns on the plots — an absent signal is still
> scientifically meaningful.

### Replicate count labels

When a biological sample is represented by more than one technical library
(e.g. two sequencing chemistries), the x-axis label is automatically
suffixed with `(x2)`, `(x3)`, etc.  The count reflects **all** technical
libraries for that sample, not just those that had detections in a
particular broad group.

### Re-rendering after edits

```bash
metamerge report \
    --input-dir <outdir>/report_inputs \
    --outdir    <outdir>/reports_custom \
    --sample-order <outdir>/report_inputs/sample_order.csv
```

The `--sample-order` flag is optional.  Without it, the default group
assignments and alphabetical ordering from the original run are used.

---

## Library linker format

The linker is a CSV with one row per sequencing library.  Required columns:

| Column | Description |
|---|---|
| `megan_library_name` | **Exact** MEGAN column header (character-for-character match) |
| `holi_library_name` | Corresponding `sample` value in the metaDMG / Holi CSV |
| `merged_library_name` | Human-readable unique name used in the merged workbook |
| `is_negative_control` | `true` or `false` |

Optional but recommended:

| Column | Description |
|---|---|
| `sample_id` | Biological sample identifier |
| `biological_sample_id` | Collapsed ID used for plotting (multiple libraries from one sample) |
| `sample_type` | Free-text type label — drives control detection (see below) |
| `is_positive_control` | `true` or `false` — explicit positive-control flag |
| `site` | Site or location code |
| `context` | Archaeological context |
| `group` | Broad group for report ordering (e.g. `SITE1_sediments`, `SITE2_bones`) |
| `chemistry` | Sequencing chemistry (e.g. `Shg`, `E-PNmp`) |
| `notes` | Free-text notes |

### Control detection keywords

MetaMerge detects library type automatically from the linker columns.  The
rules are evaluated in order; a row that matches an earlier rule is not
re-evaluated by later rules.

**Negative controls** (blank libraries) — any row where:
- `is_negative_control` = `true`, `1`, or `yes`  *(explicit flag — highest priority)*
- `sample_type` contains `"blank"` (e.g. `"extraction blank"`, `"air blank"`)
- `sample_type` contains `"negative"` (e.g. `"negative control"`)
- `sample_type` contains `"control"` but **not** `"positive"`

Negative-control libraries drive **Blank-associated** classification when their
read counts are high relative to real-library counts.  Their `holi_library_name`
is automatically set to `blanks`; damage assessment is not applied.

**Positive controls** — any row where:
- `is_positive_control` = `true`, `1`, or `yes`  *(explicit flag)*
- `sample_type` contains `"positive"` (e.g. `"positive control"`)

Positive controls are plotted in a separate `positive_controls` group and never
contribute to Blank-associated classification.

**Real samples** — all other rows (not matched by any rule above).

> **Tip:** The safest practice is to set `is_negative_control` / `is_positive_control`
> explicitly in the linker rather than relying on `sample_type` keyword detection.
> Use `sample_type` keywords as a convenient cross-check, not as a substitute.

See `examples/coprolite_demo_metadata.csv` for a worked example.

---

## Configuration

MetaMerge ships with sensible defaults in `config/defaults.yaml`.  All
thresholds can be overridden in a project-level YAML:

```yaml
thresholds:
  damage_min: 0.02          # Minimum damage to count as damage-supported
  significance_min: 2.0     # Minimum significance score
  high_confidence_n_reads_min: 100  # Read threshold for Very high confidence
  tentative_min_reads: 5    # Minimum reads for Tentative (vs Weak support)
  blank_max_fraction: 0.2   # Max blank-to-real ratio before blank-associated

lineage:
  rank_levels: [species, genus, family, order, class]
  lineage_max_steps: 2
  lineage_max_rank_level: family
  forbid_broad_names: [Mammalia, Bilateria, Boreoeutheria]
```

Pass a custom config with `--config my_config.yaml`.  Any key you specify
overrides the corresponding default; unspecified keys keep their defaults.

---

## Common names

Common-name resolution is optional and separate from classification.  When
`--online-common-names` is passed, MetaMerge queries three databases in order
and stops at the first English match:

1. **User-supplied override CSV** (`--common-names-file`) — matched by `tax_id`
   then by scientific name.
2. **Built-in convenience map** — covers frequently encountered aDNA taxa.
3. **NCBI Taxonomy** (`eutils/efetch`) — queries the `GenbankCommonName` field
   using the NCBI `tax_id`.  Most reliable for curated vertebrate / major taxon
   names.
4. **GBIF Species API** — queries vernacular names, English-language entries only
   (`language == "eng"` or `"en"`).  Non-English GBIF names are discarded.
5. **iNaturalist Taxon API** — searches by scientific name, uses
   `preferred_common_name` (English locale).  Good coverage for plants, animals,
   and fungi.
6. Blank if nothing is found.

All online results are cached to `common_name_cache.json` in the output
directory.  Subsequent runs reuse the cache and skip re-querying APIs.

The workflow is fully usable offline — the first two layers (override + built-in)
never require internet access.

---

## Holi / metaDMG columns used

MetaMerge reads the following columns from the Holi/metaDMG CSV:

| Column | Purpose |
|---|---|
| `sample` | Library identifier (matched via linker) |
| `tax_id` | Primary taxon identifier for exact matching |
| `tax_name` | Taxon name (fallback if tax_id absent) |
| `tax_rank` | Rank string |
| `tax_path` | Semicolon-delimited lineage path |
| `N_reads` | Read count |
| `N_alignments` | Alignment count |
| `damage` | Estimated C→T damage fraction |
| `significance` | Damage significance (SDs from zero) |
| `rho_Ac` | Fit-correlation metric |
| `MAP_valid` | Whether the MAP fit is valid |

Mapping/fit quality is summarised internally as `clean`, `caution`, or
`strong caution` using configurable thresholds on `N_alignments/N_reads`,
`|rho_Ac|`, and `MAP_valid`.

---

## MEGAN export format

MetaMerge expects a **wide count matrix** — one row per taxon, one column per
library — in TSV or CSV format.  The taxon-ID column may be named `tax_id`,
`taxon_id`, `TaxID`, `NCBI_taxid`, or `#Datasets` (MEGAN7's default).

Having a `tax_id` column is strongly preferred because it allows unambiguous
matching to Holi records even when taxon names differ between databases.
If only taxon names are present, matching falls back to name + rank, then
name alone.

---

## Repository layout

```text
MetaMerge/
├── README.md
├── pyproject.toml
├── config/
│   └── defaults.yaml           # Default thresholds and settings
├── docs/
│   ├── workflow.md             # Detailed workflow description
│   └── metadata_linker.md     # Linker format and generation guide
├── examples/
│   ├── coprolite_demo_metadata.csv   # Example library linker
│   ├── common_names_template.csv     # Template for custom common names
│   └── example_config.yaml           # Example project config
├── r/
│   └── metamerge_heatmap_report.R    # R heatmap rendering script
├── scripts/
│   └── make_metamerge_linker.py      # Linker generator (project setup)
├── src/
│   └── metamerge/
│       ├── __init__.py
│       ├── classify.py         # Core classification engine
│       ├── cli.py              # Command-line interface
│       ├── common_names.py     # Common-name resolution
│       ├── config.py           # Config loading and merging
│       ├── defaults.py         # Default column aliases and settings
│       ├── holi.py             # Holi/metaDMG matching and lineage support
│       ├── io.py               # MEGAN and Holi file loading
│       ├── metadata.py         # Library linker loading and validation
│       ├── report.py           # Long-format plotting table generation
│       ├── utils.py            # Utility helpers
│       └── workbook.py         # Excel workbook writer
└── tests/
    ├── test_classification.py  # classify_status() unit tests
    ├── test_matching.py        # Lineage-support logic unit tests
    └── test_smoke.py           # End-to-end smoke tests
```

---

## Running the tests

```bash
pip install -e ".[dev]"   # or just: pip install pytest
pytest
```

All 14 tests should pass.

---

## Limitations

- Direct BIOM 1 parsing is not implemented.
- Ecology, biogeography, and macrofossil corroboration are intentionally
  excluded from the core classifier.
- Lineage support requires `tax_path` to be present in the Holi/metaDMG output.
- Online common-name lookups query NCBI, GBIF, and iNaturalist sequentially;
  any of these services may be rate-limited, unavailable, or change their API.

---

## Suggested workflow for new projects

1. Export a wide MEGAN count matrix (TSV preferred, with `#Datasets` tax_id column).
2. Export the full Holi/metaDMG CSV.
3. Run `make_metamerge_linker.py` to generate the library linker from your
   project metadata.
4. Manually correct any `REVIEW_NEEDED` or `UNKNOWN` rows in the linker.
5. Run `metamerge check` to validate inputs.
6. Run `metamerge run --render-graphs` to produce the merged workbook and
   initial heatmaps.
7. Review the output, warnings, and blank-association flags.
8. Edit `report_inputs/sample_order.csv` to customise sample ordering and
   group assignments, then re-render with `metamerge report --sample-order`.
9. Apply ecology/macrofossil interpretation as a separate scientific layer.

---

## Citation

If you use MetaMerge in published work, please cite this repository and the
underlying tools (MEGAN7, metaDMG/Holi) that produce the inputs.

---

## TackleBox

MetaMerge is part of **TackleBox** — a collection of scriptable modules for
targeted palaeogenomic research developed at the Hakai Institute.  Other
modules in the collection include:

- **FlyGuide** — constructs regional NCBI nucleotide reference panels from
  GBIF exports
- **FlyForge** — designs RNA bait panels for targeted aDNA capture

MetaMerge fits naturally alongside the **Fillet** module (metagenomic
classification via BLAST/MEGAN), providing the post-classification ensemble
confidence layer on top of raw MEGAN output.

---

## Author and acknowledgments

**Author:** Tyler Murchie (Hakai Institute and McMaster University)

MetaMerge was developed with the assistance of Claude (Anthropic) and ChatGPT
(OpenAI) for code generation, documentation, and debugging support.

---

## License

MIT — see `LICENSE`.
