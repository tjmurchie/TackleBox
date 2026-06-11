# Spinner Tutorial: Running and Testing All Features

This tutorial walks through every major Spinner feature, starting with the
bundled toy data (no external tools needed) and working up to full BLAST and
database runs on real data.  All commands run from `TackleBox/Spinner/`.

---

## Setup

```bash
cd TackleBox/Spinner
chmod +x Spinner Spinner.py
pip install pyyaml          # only required dependency

./Spinner --version         # Spinner 0.7.0
./Spinner --help            # shows all subcommands
```

Run the test suite to confirm a working install:

```bash
python -m pytest tests/ -v
# Expected: 142 passed in ~1 s
```

---

## The example data

`examples/minimal_input/example_refs.fasta` has 8 deliberately crafted records.
With v0.7.0 aDNA-first defaults the expected outcomes are:

| Record | What it represents | Expected outcome |
|--------|-------------------|------------------|
| `NC_000001.1` | Clean mitochondrial reference | KEEP |
| `BADADAPT.1` | Illumina TruSeq adapter near start | REVIEW (terminal adapter) |
| `NUC18S.1` | Clean 18S nuclear marker | KEEP |
| `DUP.1` (first) | Clean plastid reference | KEEP |
| `DUP.1` (second) | Same accession used twice | REJECT (duplicate accession) |
| `SEQDUPB.1` | Different accession, identical sequence to NC_000001.1 | REJECT (duplicate sequence) |
| `LOWN.1` | 100% N-bases (80 Ns) | REVIEW (n_fraction_high −20 + homopolymer_long; score kept above review_min=30) |
| `UNVERIFIED.1` | Header contains "UNVERIFIED:" | REVIEW (bad_keyword_review; UNVERIFIED → review by default) |

> **v0.7.0 note:** `LOWN.1` and `UNVERIFIED.1` are now **REVIEW** (not REJECT).
> High-N sequences receive a score penalty; UNVERIFIED records are flagged for
> inspection rather than discarded.  Both decisions are intentional for aDNA work.

---

## Part 1: Core QC — no external tools

### Step 1: Minimal filter run

```bash
mkdir -p /tmp/spinner_test

./Spinner filter \
  --fasta examples/minimal_input/example_refs.fasta \
  --species-kingdom examples/minimal_input/example_species_kingdom.tsv \
  --regions-config configs/regions_config_example.tsv \
  --adapters configs/adapters_default.tsv \
  --bad-keywords configs/bad_keywords.tsv \
  --outprefix /tmp/spinner_test/step1
```

**Verify** the run produced all expected files:

```bash
ls /tmp/spinner_test/step1.*
```

Should see:
```
step1.command.txt
step1.decisions.tsv
step1.keep.fasta
step1.reject.fasta
step1.review.fasta
step1.run_config.resolved.yml
step1.summary.html
step1.summary.tsv
```

### Step 2: Inspect the decisions table

```bash
# Quick view — decisions column and reasons
python3 -c "
import csv, sys
with open('/tmp/spinner_test/step1.decisions.tsv') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        print(row['record_key'].ljust(20), row['decision'].ljust(8), row['reasons'])
"
```

Expected output:
```
NC_000001.1          KEEP     refseq_preferred;complete_organelle
BADADAPT.1           REVIEW   adapter_terminal
NUC18S.1             KEEP     (no reasons)
DUP.1                KEEP     (no reasons)
DUP.1__dup2          REJECT   duplicate_accession
SEQDUPB.1            REJECT   duplicate_sequence
LOWN.1               REVIEW   n_fraction_high;low_complexity;homopolymer_long
UNVERIFIED.1         REVIEW   bad_keyword_review
```

### Step 3: Explain a single record

```bash
./Spinner explain \
  --decisions /tmp/spinner_test/step1.decisions.tsv \
  --accession LOWN.1
```

### Step 4: Open the HTML report

```bash
firefox /tmp/spinner_test/step1.summary.html &
```

Look for: decision counts, marker class breakdown, per-species coverage table
(0-KEEP species highlighted in amber), adapter summary, top reasons table.

### Step 5: Verify the resolved config

```bash
cat /tmp/spinner_test/step1.run_config.resolved.yml | grep -A5 "decision_rules:"
```

Confirm `n_fraction_high` is NOT in `hard_reject_reasons`.

---

## Part 2: Audit mode (no FASTA output)

```bash
./Spinner audit \
  --fasta examples/minimal_input/example_refs.fasta \
  --adapters configs/adapters_default.tsv \
  --bad-keywords configs/bad_keywords.tsv \
  --outprefix /tmp/spinner_test/step2_audit
```

`audit` writes `decisions.tsv` but no FASTAs.  Useful for testing config
changes before a full filter run.

---

## Part 3: Custom config overrides

Create a small override YAML to test config changes:

```bash
cat > /tmp/test_config.yml << 'EOF'
# Restore strict n_fraction behaviour for testing
decision_rules:
  hard_reject_reasons:
    - adapter_internal
    - vector_internal
    - duplicate_accession
    - duplicate_sequence
    - length_below_min
    - length_above_max
    - length_below_class_min
    - n_fraction_high           # re-add: hard reject high-N records
    - non_iupac_fraction_high
    - bad_keyword_reject
    - taxonomy_cross_kingdom
    - chimera_detected
    - fcs_adaptor_hit
    - fcs_gx_contaminant
EOF

./Spinner audit \
  --fasta examples/minimal_input/example_refs.fasta \
  --bad-keywords configs/bad_keywords.tsv \
  --config /tmp/test_config.yml \
  --outprefix /tmp/spinner_test/step3_strict
```

**Verify** `LOWN.1` is now REJECT:

```bash
python3 -c "
import csv
with open('/tmp/spinner_test/step3_strict.decisions.tsv') as f:
    rows = {r['record_key']: r['decision'] for r in csv.DictReader(f, delimiter='\t')}
print('LOWN.1:', rows['LOWN.1'])   # should be REJECT
"
```

---

## Part 4: Sole representative rescue

Test that a species with only REVIEW records gets rescued:

```bash
# Create a test FASTA with one species, no clean sequences
cat > /tmp/test_sole_rep.fasta << 'EOF'
>HIGHN_SP.1 Mammuthus primigenius mitochondrion complete genome
NNNNNNNNNNNNNNNNNNNNACGTACGTACGTNNNNNNNNNNNNNNNNNNNNACGTACGT
EOF

./Spinner audit \
  --fasta /tmp/test_sole_rep.fasta \
  --outprefix /tmp/spinner_test/step4_rescue

python3 -c "
import csv
with open('/tmp/spinner_test/step4_rescue.decisions.tsv') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        print(row['record_key'], row['decision'], row['reasons'])
"
```

If the record's score sits between review_min (30) and keep_min (65), and it's
the only record for that species, `sole_representative` should appear in the
reasons and decision should be KEEP.

---

## Part 5: Gzipped FASTA input

```bash
gzip -k examples/minimal_input/example_refs.fasta

./Spinner audit \
  --fasta examples/minimal_input/example_refs.fasta.gz \
  --outprefix /tmp/spinner_test/step5_gz

wc -l /tmp/spinner_test/step5_gz.decisions.tsv
# Should be 9 lines (1 header + 8 records)
```

---

## Part 6: Per-class minimum length

```bash
cat > /tmp/test_classmin.yml << 'EOF'
basic_qc:
  min_length_by_class:
    NucMark: 500    # reject short nuclear markers
    Mito: 100
EOF

./Spinner audit \
  --fasta examples/minimal_input/example_refs.fasta \
  --regions-config configs/regions_config_example.tsv \
  --config /tmp/test_classmin.yml \
  --outprefix /tmp/spinner_test/step6_classmin

# NUC18S.1 (NucMark, ~120bp) should now be REJECT with length_below_class_min
python3 -c "
import csv
with open('/tmp/spinner_test/step6_classmin.decisions.tsv') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        if row['marker_class'] == 'NucMark':
            print(row['record_key'], row['decision'], row['reasons'])
"
```

---

## Part 7: Explain and report subcommands

```bash
# explain by accession
./Spinner explain \
  --decisions /tmp/spinner_test/step1.decisions.tsv \
  --accession NC_000001.1

# regenerate report (e.g. after manually editing decisions.tsv)
./Spinner report \
  --decisions /tmp/spinner_test/step1.decisions.tsv \
  --outprefix /tmp/spinner_test/step7_report

ls /tmp/spinner_test/step7_report.*
```

---

## Part 8: init-config

```bash
./Spinner init-config --outdir /tmp/spinner_configs

ls /tmp/spinner_configs/
# Should contain all bundled .yml and .tsv files
```

---

## Part 9: Taxonomy BLAST sanity check

Requires `blastn` and a BLAST database. Test with a small local DB first.

### Build a minimal test BLAST DB

```bash
# Grab a small FASTA of known sequences for testing
# (Replace with your actual database path for real runs)
BLAST_DB=/path/to/your/blast/database/nt

./Spinner screen-taxonomy \
  --fasta examples/minimal_input/example_refs.fasta \
  --species-kingdom examples/minimal_input/example_species_kingdom.tsv \
  --blast-db "$BLAST_DB" \
  --outprefix /tmp/spinner_test/step9_taxonomy
```

Or via config:

```bash
cat > /tmp/test_blast.yml << 'EOF'
steps:
  taxonomy_blast: true
taxonomy_blast:
  blast_db: /path/to/blast/nt
  taxdump_dir: /path/to/taxdump
  num_threads: 8
  min_pident: 70.0
EOF

./Spinner audit \
  --fasta my_refs.fasta \
  --species-kingdom my_species_kingdom.tsv \
  --config /tmp/test_blast.yml \
  --outprefix /tmp/spinner_test/step9_blast_cfg
```

**Verify results:**

```bash
python3 -c "
import csv
from collections import Counter
with open('/tmp/spinner_test/step9_blast_cfg.decisions.tsv') as f:
    statuses = Counter(r['taxonomy_status'] for r in csv.DictReader(f, delimiter='\t'))
for k, v in statuses.most_common():
    print(f'{v:>6}  {k}')
"
```

Expected for clean data: mostly `PASS_SPECIES` / `PASS_GENUS` / `NOT_CHECKED`
(records shorter than BLAST thresholds).

---

## Part 10: MMSeqs2 (fast taxonomy search)

MMSeqs2 must be installed and a taxonomy-indexed database must exist.

```bash
cat > /tmp/test_mmseqs.yml << 'EOF'
steps:
  taxonomy_blast: true
taxonomy_blast:
  method: mmseqs2
  blast_db: /path/to/mmseqs2/nt_tax_db
  num_threads: 16
  min_pident: 70.0
EOF

./Spinner audit \
  --fasta my_refs.fasta \
  --species-kingdom my_species_kingdom.tsv \
  --config /tmp/test_mmseqs.yml \
  --outprefix /tmp/spinner_test/step10_mmseqs
```

The output format is identical to BLAST — compare `step9_blast_cfg.decisions.tsv`
and `step10_mmseqs.decisions.tsv` and results should match closely.

---

## Part 11: Windowed BLAST (chimerism detection)

Only runs on sequences ≥ 1000 bp by default.

```bash
cat > /tmp/test_windowed.yml << 'EOF'
steps:
  taxonomy_blast: true
  windowed_blast: true
taxonomy_blast:
  blast_db: /path/to/blast/nt
  taxdump_dir: /path/to/taxdump
  num_threads: 8
windowed_blast:
  taxdump_comparison_rank: family
  num_threads: 8
EOF

./Spinner audit \
  --fasta my_long_refs.fasta \
  --config /tmp/test_windowed.yml \
  --outprefix /tmp/spinner_test/step11_windowed

# Check windowed BLAST results
python3 -c "
import csv
from collections import Counter
with open('/tmp/spinner_test/step11_windowed.decisions.tsv') as f:
    statuses = Counter(r['windowed_status'] for r in csv.DictReader(f, delimiter='\t'))
for k, v in statuses.most_common():
    print(f'{v:>6}  {k}')
"
```

---

## Part 12: Chimera screen (vsearch uchime)

Requires vsearch to be installed.

```bash
cat > /tmp/test_chimera.yml << 'EOF'
steps:
  chimera_screen: true
chimera_screen:
  method: uchime_denovo
  reject_chimeras: true
  review_borderline: true
EOF

./Spinner audit \
  --fasta my_refs.fasta \
  --config /tmp/test_chimera.yml \
  --outprefix /tmp/spinner_test/step12_chimera

python3 -c "
import csv
from collections import Counter
with open('/tmp/spinner_test/step12_chimera.decisions.tsv') as f:
    reasons = Counter(r for row in csv.DictReader(f, delimiter='\t')
                      for r in row['reasons'].split(';') if r)
print('chimera_detected:', reasons.get('chimera_detected', 0))
print('chimera_borderline:', reasons.get('chimera_borderline', 0))
"
```

---

## Part 13: Clustering (vsearch)

Requires vsearch to be installed.

```bash
./Spinner cluster \
  --fasta examples/minimal_input/example_refs.fasta \
  --species-kingdom examples/minimal_input/example_species_kingdom.tsv \
  --outprefix /tmp/spinner_test/step13_cluster

python3 -c "
import csv
from collections import Counter
with open('/tmp/spinner_test/step13_cluster.decisions.tsv') as f:
    roles = Counter(r['cluster_role'] for r in csv.DictReader(f, delimiter='\t'))
print(dict(roles))
"
```

---

## Iterative testing strategy

Use this order when setting up a new project:

### Stage 1 — Toy data validation (15 minutes)
Run Steps 1–4 above on the example data.  If KEEP/REVIEW/REJECT counts match
expectations, your install is working correctly.

### Stage 2 — Small real subset (~100 records)
Pick one well-characterised taxon group from your actual FASTA (e.g. all
mammalian mitochondrial sequences).  Run `audit` (no FASTA output) and inspect
`decisions.tsv`:
- What fraction are KEEP?  (Target: ≥ 70% for curated databases)
- What are the top 5 reasons?  (Expect: `duplicate_sequence`, `n_fraction_high`,
  `duplicate_accession`)
- Are any clean records getting flagged?  Use `explain` to investigate.

### Stage 3 — Config tuning
Adjust thresholds based on Stage 2 observations:
- Too many high-N records from aDNA targets → confirm `n_fraction_high` is NOT
  a hard reject (check resolved config)
- Unexpected REJECTs for rare taxa → check `rescue_sole_representatives: true`
- Too many `NO_EXPECTED_MATCH` → confirm `review_if_no_expected_match: false`

### Stage 4 — Full database run
Run `filter` on the complete FASTA.  Check the HTML report:
- Per-species coverage table: 0-KEEP species should be few; amber rows indicate
  where sole representative rescue fired
- Top reasons: if > 30% of records share a single reason, investigate
- Kingdom breakdown: cross-kingdom counts should be near zero (check if
  `taxonomy_cross_kingdom` appears unexpectedly)

### Stage 5 — BLAST validation (optional but recommended)
Enable `taxonomy_blast` on the `keep.fasta` output from Stage 4.  A high
`PASS_GENUS` rate confirms the database is taxonomically consistent.

### Validation checks after each stage

```bash
# Count decisions
python3 -c "
import csv; from collections import Counter
with open('OUTPREFIX.decisions.tsv') as f:
    c = Counter(r['decision'] for r in csv.DictReader(f, delimiter='\t'))
total = sum(c.values())
for d in ('KEEP','REVIEW','REJECT'):
    print(f'{d}: {c[d]:,} ({100*c[d]/total:.1f}%)')
"

# Top 10 reasons
python3 -c "
import csv; from collections import Counter
with open('OUTPREFIX.decisions.tsv') as f:
    c = Counter(r for row in csv.DictReader(f, delimiter='\t')
                for r in row['reasons'].split(';') if r)
for r, n in c.most_common(10):
    print(f'{n:>6}  {r}')
"

# All REJECT records
python3 -c "
import csv
with open('OUTPREFIX.decisions.tsv') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        if row['decision'] == 'REJECT':
            print(row['record_key'].ljust(25), row['reasons'])
" | sort | uniq -c | sort -rn | head -30
```

---

## Common issues and what to check

| Observation | Likely cause | What to do |
|-------------|--------------|------------|
| Many `taxonomy_not_checked` | BLAST step disabled | Enable `steps.taxonomy_blast` |
| Many `NO_EXPECTED_MATCH` | Rare taxa absent from DB | Normal; `review_if_no_expected_match: false` (default) |
| Clean records REJECT with `n_fraction_high` | Old config with n_fraction in hard_reject | Check resolved config; remove from `hard_reject_reasons` |
| All records REVIEW with `bad_keyword_review` | `UNVERIFIED` in old bad_keywords.tsv | Use new `configs/bad_keywords.tsv` (v0.7.0) |
| Unexpected `duplicate_sequence` counts | Large 16S / ITS databases with NCBI redundancy | Expected; ~90% of 16S sequences in nt are duplicates |
| 0 records rescued by sole representative | All species already have ≥1 KEEP | Good — no rescue needed |
| BLAST hangs | Large nt database on NFS mount, I/O bound | Add `--keep-temp` and use `screen-taxonomy` subcommand to resume |
