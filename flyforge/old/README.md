# TackleBox: FlyForge

**FlyForge** is a comprehensive bait/probe designer for hybridization capture enrichment with integrated support for in-house probe synthesis via oligo pools. It designs custom bait panels for single-species targets (e.g., mitogenomes, nuclear loci) or large multi-target projects with thousands of reference sequences (e.g., environmental DNA panels, gene family surveys).

FlyForge is part of the **TackleBox** suite of tools for targeted enrichment workflows.

---

## Features

- **Flexible tiling** with configurable bait length and density
- **Iterative density adjustment** to approach a target bait count (`--max-baits`)
- **Preprocessing**: uppercase, U→T, N-run replacement, self-repeat softmasking
- **Length-aware target handling**: drop short references, pad near-length references
- **Input complementary region removal** via iterative self-BLAST (recommended for multi-gene inputs)
- **Multi-stage quality filtering**: ambiguous bases, repeat-masked fraction, melting temperature
- **LguI/BspQI restriction site filtering** for oligo pool compatibility
- **Perfect complement removal** from the probe set
- **Dual BLAST filtering modes**:
  - Percent-identity filter (`--blast-db`) — remove probes matching off-target databases above a configurable %ID threshold
  - nident-based specificity filter (`--specificity-db`) — CARPDM-style taxonomy-aware filtering against NCBI nt
- **Self-BLAST redundancy collapsing** with iterative greedy algorithm
- **cd-hit-est clustering** to collapse near-identical baits
- **O-pool design**: T7 transcription promoter + optimized amplification primer + BspQI cut site, following Hackenberger et al. (2024)
- **Primer design**: exhaustive candidate generation with Tm matching, hairpin/homo/heterodimer checking, and BLAST-based selection
- **BLAST validation**: probes vs. input sequences with position-level coverage tracking
- **Comprehensive visualization**: violin plots (GC, Tm, coverage), per-target coverage line plots (PNG + SVG)
- **Both bare probes AND oligo pool outputs** in a single run
- **Progress dashboard** with elapsed time, ETA, step checklist, and progress bar
- **End-of-run summary** printed to screen with key statistics

---

## Requirements

- **Python** 3.8+
- **Python packages**: biopython, primer3-py, matplotlib, seaborn, pandas, numpy, tqdm
- **External tools**:
  - [BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) (`blastn`, `makeblastdb` in PATH)
  - [CD-HIT](http://weizhongli-lab.org/cd-hit/) (`cd-hit-est` in PATH)

### Installation (conda)

```bash
conda create -n flyforge python=3.12
conda activate flyforge

# External tools
conda install -c bioconda blast cd-hit

# Python packages
pip install biopython primer3-py matplotlib seaborn pandas numpy tqdm
```

Make executable:

```bash
chmod +x FlyForge.py FlyForge
```

---

## Quick Start

```bash
# Basic bait design (80 bp probes, 4x tiling, with o-pool synthesis)
FlyForge -i target.fasta --prefix my_baits --tiling-density 4 --threads 4

# Skip o-pool design (ordering pre-synthesized RNA probes directly)
FlyForge -i target.fasta --prefix my_baits --no-opool

# With BLAST off-target exclusion (e.g., against a bacteria database)
FlyForge -i target.fasta --prefix my_baits \
    --blast-db /path/to/exclusion_db --blast-min-pident 80 --threads 8

# With CARPDM-style specificity filter against NCBI nt
FlyForge -i target.fasta --prefix my_baits \
    --specificity-db /path/to/nt --threads 8
```

---

## Recommended Settings

Below are recommended configurations for common use cases in targeted enrichment research. All examples assume 80 bp probes (the default) and oligo pool synthesis. Adjust `--threads` to match your system.

### Scenario 1: Single-species target enrichment from ancient DNA

**Use case**: Enriching a single target organism (e.g., mammal mitogenome, nuclear loci) from a degraded, low-DNA sample dominated by environmental bacteria.

**Key considerations**:
- High tiling density ensures good coverage even with degraded, fragmented DNA
- Quality filters remove low-complexity and repetitive probes that could bind non-specifically
- No off-target BLAST filter is needed for most single-species cases — probes designed from a mammal/bird/plant reference are inherently too divergent from bacteria to cross-hybridize (typically <50% sequence identity)
- If your target includes highly conserved loci (e.g., rRNA, EF1α), consider adding `--blast-db` with a bacterial database to remove probes with unexpected cross-homology
- Permissive enough to cross-capture closely related species (e.g., congenerics), which is often desirable for aDNA work

```bash
FlyForge -i reference.fasta --prefix my_species \
    --bait-length 80 \
    --tiling-density 4 \
    --min-tm 50 \
    --remove-complements \
    --threads 8
```

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| `--tiling-density 4` | 20 bp step size | Dense coverage compensates for DNA fragmentation and loss |
| `--min-tm 50` | 50°C minimum Tm | Removes probes with poor hybridization potential |
| `--remove-complements` | on | Removes self-complementary regions that could cause probe-probe binding |

**Optional additions**:
- `--max-baits N` — constrain pool size to budget (auto-adjusts tiling density)
- `--blast-db /path/to/bacteria_db --blast-min-pident 80 --blast-max-hits 5` — remove probes with ≥80% identity to bacteria (only needed if targets include conserved loci)
- `--no-opool` — skip oligo pool design if ordering pre-synthesized probes

### Scenario 2: Multi-species target enrichment from environmental samples

**Use case**: Enriching multiple target organisms (e.g., plants, mammals, birds, fish) from environmental DNA (eDNA) or ancient environmental samples (sedaDNA), while avoiding enrichment of bacteria, fungi, nematodes, and other non-target organisms.

**Key considerations**:
- `--remove-complements` is critical for multi-gene/multi-species inputs to prevent cross-complementary probe interference
- The `--blast-db` percent-identity filter is the recommended way to exclude non-target organisms — point it at a BLAST database of organisms you want to *avoid* enriching
- Use a moderate `--blast-min-pident` threshold (75–80%) to remove probes with meaningful cross-homology while preserving probes that may cross-capture related target species
- Self-BLAST redundancy and cd-hit-est clustering reduce probe set size and remove redundant probes across targets
- `--max-baits` helps stay within oligo pool ordering limits (e.g., Twist Bioscience 200 nt max oligo, various pool size tiers)

```bash
FlyForge -i targets/*.fasta --prefix env_panel \
    --bait-length 80 \
    --tiling-density 3 \
    --min-tm 50 \
    --remove-complements \
    --blast-db /path/to/exclusion_db \
    --blast-min-pident 80 \
    --blast-max-hits 5 \
    --cluster-identity 0.95 \
    --threads 16
```

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| `--tiling-density 3` | ~27 bp step size | Good coverage while managing probe count across many targets |
| `--min-tm 50` | 50°C minimum Tm | Ensures probes with good hybridization thermodynamics |
| `--remove-complements` | on | Prevents probe self-annealing artifacts from multi-gene inputs |
| `--blast-db` | exclusion DB | Database of organisms to avoid enriching |
| `--blast-min-pident 80` | 80% identity cutoff | Removes probes with strong homology to non-targets while allowing cross-capture of closely related species (congenerics, confamiliars) |
| `--blast-max-hits 5` | 5 hits per query | Checks more potential off-target matches per probe |
| `--cluster-identity 0.95` | 95% clustering | Collapses near-identical probes across overlapping targets |

**Optional additions**:
- `--max-baits N` — hard cap on final probe count to fit pool size/budget constraints
- `--probe-num-cutoff N` — target count for self-BLAST redundancy collapsing
- `--specificity-db /path/to/nt --threads N` — CARPDM-style nident filter against NCBI nt (see note below)

### Building an exclusion database

For Scenarios 1 and 2, the `--blast-db` filter works with any standard BLAST nucleotide database. To build a custom exclusion database from organisms you want to avoid enriching:

```bash
# Download or collect FASTA sequences of non-target organisms
# (bacteria, fungi, nematodes, etc.)
cat bacteria.fasta fungi.fasta nematodes.fasta > exclusion_seqs.fasta

# Build BLAST database
makeblastdb -in exclusion_seqs.fasta -dbtype nucl -out exclusion_db
```

Alternatively, you can use pre-built databases. For example, to use NCBI nt as the exclusion database, set `--blast-db /path/to/nt`. Any probe with ≥`--blast-min-pident`% identity to *anything* in the database will be removed, so choose your database and threshold carefully:

- **Narrow exclusion DB** (e.g., bacteria + fungi only): You can use a lower threshold (~75%) since every hit is an organism you want to exclude
- **Broad exclusion DB** (e.g., full NCBI nt): You need a higher threshold (~90%) since the database also contains your target organisms and relatives — only remove probes matching clearly unrelated sequences

### A note on the `--specificity-db` filter

The `--specificity-db` filter uses the nident-based approach from CARPDM (Hackenberger et al. 2024). When pointed at a database named `nt`, it enables taxonomy-aware filtering that:
- Retrieves the `sskingdom` field for each BLAST hit
- Removes probes matching **Archaea, Eukaryota, or Viruses** within an intermediate identity range (62.5%–97.5% of probe length)
- **Preserves probes matching Bacteria** regardless of identity

This was designed for **antimicrobial resistance gene enrichment** where bacterial targets are desired and eukaryotic/viral cross-hybridization is not. It is **not appropriate** for enrichment of eukaryotic targets (mammals, birds, plants, etc.) where you want to *keep* eukaryotic matches and *exclude* bacteria.

For eukaryotic target enrichment, use the `--blast-db` percent-identity filter with a custom exclusion database as described above.

---

## Pipeline Overview

1. **Preprocessing**
   - Uppercase, U→T conversion, N-run (1–10 bp) → T-run replacement
   - Count modified bases per reference

2. **Complementary Region Removal** (optional, `--remove-complements`)
   - Iterative self-BLAST to find and remove complementary regions between input targets
   - Prevents self-hybridization in probe sets from multi-gene inputs

3. **Self-Repeat Softmasking** (optional, skip with `--skip-self-mask`)
   - k-mer counting across all references
   - Positions covered by repetitive k-mers are lowercase-masked

4. **Tiling**
   - Sliding window with configurable bait length and step size
   - Iterative density adjustment if `--max-baits` is set
   - Short sequences padded or dropped as appropriate

5. **Ambiguous Base Filter**
   - Remove baits with ≥ `--ambiguous-cutoff` ambiguous bases

6. **Repeat-Masked Fraction Filter**
   - Remove baits with masked fraction > `--max-masked-frac`

7. **Melting Temperature Filter** (optional, `--min-tm`)
   - Remove baits with Tm below threshold using RNA/DNA nearest-neighbor model

8. **LguI/BspQI Site Filter** (when building o-pool)
   - Remove baits containing GAAGAGC or GCTCTTC to prevent unwanted cutting during synthesis

9. **Complementary Bait Removal**
   - Remove one of each pair of perfect reverse-complement baits

10. **BLAST Percent-Identity Filter** (optional, `--blast-db`)
    - Remove baits with ≥ `--blast-min-pident` identity to off-target database

11. **BLAST Specificity Filter** (optional, `--specificity-db`)
    - CARPDM-style nident-based filtering
    - Taxonomy-aware mode for NCBI nt (removes non-bacterial hits within identity range)
    - See "A note on the `--specificity-db` filter" above

12. **Self-BLAST Redundancy Filter** (optional, skip with `--no-redundancy`)
    - Self-BLAST to identify complementary probes (nident ≥ 30, minus strand)
    - Iterative identity-based redundancy collapsing (greedy algorithm)
    - Configurable probe count target (`--probe-num-cutoff`)

13. **cd-hit-est Clustering** (optional, skip with `--no-cluster`)
    - Collapse highly similar baits by identity and overlap

14. **O-Pool Design** (optional, skip with `--no-opool`)
    - Append T7 promoter (5′-GCTAATACGACTCACTATAGGG-3′) to 5′ end
    - Generate and screen candidate amplification primers (18 bp):
      - Tm matched to T7 primer (±0.1°C)
      - No homopolymer runs (AAA, TTT, GGG, CCC)
      - No LguI/BspQI sites
      - Low hairpin, homodimer, heterodimer Tm (< 20°C)
    - BLAST primers against probe set; select primer with fewest identities
    - Append reverse complement of selected primer to 3′ end

    Final oligo structure:
    ```
    5′-GCTAATACGACTCACTATAGGG [80 bp probe] RC_primer -3′
       T7 promoter             core probe   amplification primer
    ```

15. **BLAST Validation**
    - BLAST final probes against original input targets
    - Per-position coverage tracking
    - Per-target and per-probe summary statistics

16. **Visualization**
    - Individual target coverage plots (PNG + SVG)
    - Violin plots: probe GC, Tm, target count, target length, coverage proportion/depth/stdev

---

## Output Files

```
output_dir/
├── {prefix}_final_baits.fa               # Bare probe sequences (for QC/comparison)
├── {prefix}_probes.fna                    # Bare probe sequences (same content)
├── {prefix}_oligo_pool.fna                # Full oligo sequences with T7 + primer (for ordering)
├── {prefix}_amplification_primers.fna     # T7 and amplification primer sequences
├── {prefix}_summary.tsv                   # Parameters + step-wise statistics
├── {prefix}_per_ref_stats.tsv             # Per-reference coverage and bait counts
├── {prefix}_progress.log                  # Step timings and ETA log
├── {prefix}_target_info.csv               # Per-target summary (length, GC, coverage)
├── {prefix}_probe_info.csv                # Per-probe summary (GC, Tm, target count)
├── {prefix}_target_probe_pairs.csv        # Target–probe pairings
├── {prefix}_final_blast.xml               # Raw BLAST validation results
└── {prefix}_plots/
    ├── individual_target_coverages/        # Per-target coverage line plots (PNG + SVG)
    ├── probe_gc.{png,svg}
    ├── probe_tm.{png,svg}
    ├── probe_num_targets.{png,svg}
    ├── target_len.{png,svg}
    ├── target_gc.{png,svg}
    ├── target_probe_count.{png,svg}
    ├── target_coverage_prop.{png,svg}
    ├── target_coverage_depth.{png,svg}
    └── target_coverage_stdev.{png,svg}
```

---

## Command-Line Reference

### Core I/O

| Flag | Default | Description |
|------|---------|-------------|
| `-i, --input` | Required | Input FASTA file(s) |
| `--prefix` | Required | Prefix for all output files |
| `--output-dir` | `flyforge_output` | Output directory |
| `--progress-log` | Auto | Progress log file path |

### Design Parameters

| Flag | Default | Description |
|------|---------|-------------|
| `--bait-length` | 80 | Bait/probe length in bp |
| `--tiling-density` | 3.0 | Tiling density (step = bait_length / density) |
| `--omit-short-leq` | 70 | Drop references with length ≤ this |
| `--pad-min` | 71 | Pad references ≥ this and < bait_length to full bait length |
| `--min-tm` | 0 (off) | Minimum melting temperature filter (recommended: 50) |

### Filtering

| Flag | Default | Description |
|------|---------|-------------|
| `--ambiguous-cutoff` | 10 | Remove baits with ≥ N ambiguous bases |
| `--max-masked-frac` | 0.25 | Max repeat-masked fraction |
| `--remove-complements` | off | Remove complementary regions in input sequences |

### Self-Repeat Masking

| Flag | Default | Description |
|------|---------|-------------|
| `--repeat-k` | 15 | k-mer size for masking |
| `--repeat-threshold` | 3 | Min count to flag k-mer as repeat |
| `--skip-self-mask` | off | Skip masking (assume pre-masked input) |

### BLAST Percent-Identity Filter

| Flag | Default | Description |
|------|---------|-------------|
| `--blast-db` | None | BLAST database for percent-identity off-target filtering |
| `--blast-evalue` | 1e-5 | BLAST e-value cutoff |
| `--blast-min-pident` | 90.0 | Min percent identity to remove a bait |
| `--blast-max-hits` | 1 | max_target_seqs per query |

### BLAST Specificity Filter (CARPDM-style)

| Flag | Default | Description |
|------|---------|-------------|
| `--specificity-db` | None | BLAST database for nident-based specificity filtering. If database basename is `nt`, enables taxonomy-aware mode (designed for bacterial AMR enrichment) |

### Redundancy & Clustering

| Flag | Default | Description |
|------|---------|-------------|
| `--no-redundancy` | off | Skip self-BLAST redundancy filter |
| `--probe-num-cutoff` | 100000 | Target probe count for redundancy collapsing |
| `--no-cluster` | off | Skip cd-hit-est clustering |
| `--cluster-identity` | 0.95 | cd-hit-est identity threshold |
| `--cluster-overlap` | 0.83 | cd-hit-est overlap threshold |

### O-Pool Synthesis

| Flag | Default | Description |
|------|---------|-------------|
| `--no-opool` | off | Skip oligo pool design (T7 + primer) |

### Performance & Constraints

| Flag | Default | Description |
|------|---------|-------------|
| `--max-baits` | None | Target max bait count (auto-adjusts tiling density) |
| `--min-tiling-density` | 1.0 | Floor for density auto-adjustment |
| `--max-total-bp` | 0 (off) | Abort if total reference bp exceeds this |
| `--no-coverage` | off | Skip per-base coverage stats |
| `--threads` | 1 | Threads for BLAST and cd-hit-est |

---

## Oligo Pool Synthesis Background

FlyForge designs oligo pools for in-house RNA probe synthesis following the approach of Hackenberger et al. (2024). The key components:

1. **T7 Promoter** (`5′-GCTAATACGACTCACTATAGGG-3′`): Appended to the 5′ end of each probe. Includes a GC clamp and three terminal guanines that are incorporated as the first nucleotides of the RNA transcript during in vitro transcription.

2. **BspQI Cut Site** (`GCTCTTC`): Embedded in the amplification primer. After PCR amplification of the oligo pool, BspQI digestion cleaves downstream of the recognition site, producing properly-sized probe templates with defined 3′ ends.

3. **Amplification Primer**: An 18 bp primer is computationally designed and selected to minimize cross-hybridization with the probe set. It is appended (reverse complemented) to the 3′ end of each oligo. Primer selection criteria:
   - Tm matched to T7 primer (±0.1°C) for uniform amplification
   - No homopolymer runs (≥3 identical bases)
   - No internal LguI/BspQI restriction sites
   - Low hairpin, homodimer, and heterodimer Tm (< 20°C)
   - BLAST-based selection of the candidate with fewest total identities to the probe set

The final oligo pool can be ordered from a provider like Twist Bioscience at the base price (up to 200 nt per oligo at 80 bp probe length, giving a 120 nt total oligo: 21 nt T7 + 80 nt probe + 18 nt RC primer + 1 nt BspQI overhang).

### Workflow after ordering

1. Receive oligo pool from Twist Bioscience (or similar provider)
2. PCR amplify the pool using the T7 and amplification primers (output in `{prefix}_amplification_primers.fna`)
3. Perform in vitro transcription with T7 RNA polymerase to produce biotinylated RNA probes
4. Digest with BspQI to remove the primer-derived 3′ tail
5. Use RNA probes for hybridization capture enrichment

---

## Citation

If you use FlyForge in your research, please cite this software and the CARPDM publication for the oligo pool synthesis methodology:

> Hackenberger D, Imtiaz H, Raphenya AR, Smith KW, Alcock BP, Poinar HN, Wright GD, McArthur AG. 2024. CARPDM: cost-effective antibiotic resistome profiling of metagenomic samples using targeted enrichment. *bioRxiv* 2024.03.27.587061. doi:[10.1101/2024.03.27.587061](https://doi.org/10.1101/2024.03.27.587061)

See `CITATION.cff` for machine-readable citation metadata.

---

## License

MIT License. See [LICENSE](LICENSE) for details.
