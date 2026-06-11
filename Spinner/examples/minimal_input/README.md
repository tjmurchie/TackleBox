# Minimal Spinner example

From the `Spinner/` directory:

```bash
./Spinner filter \
  --fasta examples/minimal_input/example_refs.fasta \
  --species-kingdom examples/minimal_input/example_species_kingdom.tsv \
  --regions-config configs/regions_config_example.tsv \
  --adapters configs/adapters_default.tsv \
  --bad-keywords configs/bad_keywords.tsv \
  --config configs/spinner_default.yml \
  --outprefix examples/expected_outputs/example_spinner
```

Expected outputs include:

- `example_spinner.decisions.tsv`
- `example_spinner.summary.tsv`
- `example_spinner.summary.html`
- `example_spinner.keep.fasta`
- `example_spinner.review.fasta`
- `example_spinner.reject.fasta`
