"""Configuration loading and defaults.

Call ``load_config(path)`` to get a fully-resolved config dict merging
DEFAULT_CONFIG with any user YAML overrides.
"""
from __future__ import annotations

from typing import Any, Dict

DEFAULT_CONFIG: Dict[str, Any] = {
    "run": {
        "profile_name": "default",
        "progress_bars": True,
        "write_reject_fasta": True,
        "write_review_fasta": True,
        "fail_on_missing_external_tool": False,
        "keep_temp_files": False,
    },
    "inputs": {
        "regions_config": "",
        "species_kingdom": "",
        "adapters_tsv": "",
        "bad_keywords_tsv": "",
    },
    "steps": {
        "basic_qc": True,
        "classify_regions": True,
        "bad_keyword_screen": True,
        "adapter_screen": True,
        "vector_screen": False,
        "taxonomy_blast": False,
        "windowed_blast": False,
        "chimera_screen": False,        # vsearch uchime chimera detection
        "cluster": False,
        "cap_references": False,
        "fcs_adaptor": False,
        "fcs_gx": False,
        "report": True,
    },
    "basic_qc": {
        "min_length": 50,
        "max_length": 400000,
        "min_length_by_class": {},      # e.g. {"NucMark": 200, "Mito": 100} — overrides min_length per class
        "max_n_fraction": 0.20,
        "max_non_iupac_fraction": 0.01,
        "min_shannon_entropy": 1.2,
        "max_homopolymer_run": 60,
        "remove_exact_duplicate_sequences": True,
        "remove_duplicate_accessions": True,
        "rescue_cross_species_duplicates": True,  # promote one rep per species when all records are cross-species dups
    },
    "classification": {
        "default_class_if_no_match": "Other",
        "coerce_unknown_plastid_to_plant": True,
    },
    "adapter_screen": {
        "max_mismatch": 2,
        "scan_reverse_complement": True,
        "reject_internal_hits": True,
        "review_terminal_hits": True,
        "terminal_window_bp": 25,
        "min_adapter_match_length": 12,
        "default_action": "reject",
        "trim_and_rescreen": False,     # trim at adapter position and re-check; opt-in
    },
    "bad_keywords": {
        "case_sensitive": False,
        "default_action": "review",
    },
    "vector_screen": {
        "method": "blastn",
        "blast_db": "",
        "blast_task": "blastn-short",
        "min_hit_length": 20,
        "min_pident": 90.0,
        "reject_internal_hits": True,
        "review_terminal_hits": True,
        "terminal_window_bp": 25,
        "outfmt": "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore",
    },
    "taxonomy_blast": {
        "method": "blastn",             # "blastn" or "mmseqs2"
        "blast_db": "",
        "blast_task": "megablast",
        "search_type": 1,               # MMSeqs2 only: 1=amino acid, 2=translated nuc→protein, 3=nuc-nuc
        "max_query_length": 0,          # skip sequences longer than this for taxonomy search (0 = no limit)
        "tmp_dir": "",                  # MMSeqs2 temp dir — defaults to system /tmp; set to local path if /tmp is small
        "max_target_seqs": 10,
        "max_hsps": 5,
        "evalue": "1e-10",
        "min_qcov": 50.0,
        "min_pident": 70.0,
        "num_threads": 1,
        "expected_taxon_level": "genus",
        "reject_cross_kingdom": True,
        "review_if_no_expected_match": False,
        "outfmt": "6 qseqid saccver pident length qlen qstart qend evalue bitscore staxids sscinames",
        "taxdump_dir": "",
        # Multi-database cross-kingdom escalation: sequences rejected by the primary DB
        # are re-checked against NR protein then NT nucleotide before final rejection.
        "escalate_cross_kingdom": False,
        "nr_protein_db": "",            # MMSeqs2 NR protein DB path for level-1 escalation
        "nr_mmseqs_binary": "",         # override mmseqs binary for NR search (if different env)
        "nt_blast_db": "",              # BLAST NT nucleotide DB path for level-2 escalation
        "escalation_nr_min_pident": 0.0,   # 0 = use min_pident; override for NR-specific threshold
        "escalation_nt_min_pident": 80.0,  # megablast-appropriate threshold for NT nucleotide check
    },
    "windowed_blast": {
        "blast_db": "",
        "blast_task": "megablast",
        "enabled_for_min_length": 1000,
        "window_size": 500,
        "window_step": 250,
        "max_target_seqs": 3,
        "max_hsps": 1,
        "evalue": "1e-10",
        "min_pident": 80.0,
        "num_threads": 1,
        "min_conflicting_windows": 2,
        "conflict_action": "review",
        "taxdump_comparison_rank": "family",  # rank for lineage-aware conflict detection
        "outfmt": "6 qseqid saccver pident length qlen qstart qend evalue bitscore staxids sscinames",
    },
    "chimera_screen": {
        "method": "uchime_denovo",      # "uchime_denovo" (no DB needed) or "uchime_ref" (needs reference_db)
        "reference_db": "",             # path to reference FASTA for uchime_ref mode
        "vsearch_path": "vsearch",
        "reject_chimeras": True,        # Y verdict → chimera_detected (hard reject by default)
        "review_borderline": True,      # ? verdict → chimera_borderline (review)
        "abskew": 2.0,                  # vsearch --abskew; lower = more sensitive
    },
    "cluster": {
        "method": "vsearch",
        "identity": 0.99,
        "by": ["species_guess", "marker_class"],
        "vsearch_path": "vsearch",
        "cdhit_path": "cd-hit-est",
        "max_reps_per_cluster": 1,
        "nonrepresentative_action": "review",
    },
    "capping": {
        "mode": "species_marker",
        "uncapped_classes": [],
        "max_per_species_marker": {"Mito": 20, "Plastid": 20, "NucMark": 10, "Other": 5},
        "max_per_species_total": 50,
        "cap_action": "review",
        "rescue_sole_representatives": True,   # if a species has 0 KEEP after all filtering,
                                               # promote the best REVIEW record to KEEP
    },
    "decision_rules": {
        "hard_reject_reasons": [
            "adapter_internal",
            "vector_internal",
            "duplicate_accession",
            "duplicate_sequence",
            "length_below_min",
            "length_above_max",
            "length_below_class_min",
            # n_fraction_high is NOT a hard reject: ancient assemblies can have high N content.
            # A high-N sequence from a rare or extinct species is still better than no reference.
            "non_iupac_fraction_high",
            "bad_keyword_reject",
            "taxonomy_cross_kingdom",
            "chimera_detected",
            "fcs_adaptor_hit",
            "fcs_gx_contaminant",
        ],
        "review_reasons": [
            "adapter_terminal",
            "vector_terminal",
            "low_complexity",
            "homopolymer_long",
            "bad_keyword_review",
            "taxonomy_no_expected_match",
            "taxonomy_not_checked",
            "windowed_blast_conflict",
            "chimera_borderline",
            "cluster_nonrepresentative",
            "cap_exceeded",
            "fcs_adaptor_review",
            "fcs_gx_review",
        ],
        "score_thresholds": {"keep_min": 65, "review_min": 30},
    },
    "scoring": {
        "start": 100,
        "adapter_internal": -100,
        "adapter_terminal": -40,
        "vector_internal": -100,
        "vector_terminal": -50,
        "duplicate_accession": -100,
        "duplicate_sequence": -100,
        "length_below_min": -100,
        "length_above_max": -100,
        "length_below_class_min": -100,
        "n_fraction_high": -20,         # soft penalty only (not in hard_reject_reasons)
        "non_iupac_fraction_high": -40,
        "low_complexity": -25,
        "homopolymer_long": -20,
        "bad_keyword_reject": -100,
        "bad_keyword_review": -30,
        "taxonomy_cross_kingdom": -100,
        "taxonomy_no_expected_match": -10,  # reduced: rare/ancient taxa often absent from BLAST DB
        "taxonomy_same_species": 20,
        "taxonomy_same_genus": 10,
        "refseq_preferred": 10,
        "voucher_keyword": 5,
        "complete_organelle": 5,
        "chimera_detected": -100,
        "chimera_borderline": -60,
        "rescued_duplicate": 0,         # neutral — let other criteria decide KEEP/REVIEW
        "adapter_trimmed_rescued": 5,   # small positive: sequence survived trimming and re-screen
        "taxonomy_rescued_nr_protein": 0,   # cleared cross-kingdom by NR protein
        "taxonomy_rescued_nt_blast": 0,     # cleared cross-kingdom by NT nucleotide
        "sole_representative": 10,
        "cluster_representative": 5,
        "cluster_nonrepresentative": -25,
        "cap_exceeded": -30,
        "windowed_blast_conflict": -60,
        "fcs_adaptor_hit": -100,
        "fcs_adaptor_review": -40,
        "fcs_gx_contaminant": -100,
        "fcs_gx_review": -40,
    },
    "fcs_adaptor": {
        "command": "",
        "results_tsv": "",
        "reject_on_any_hit": True,
    },
    "fcs_gx": {
        "command": "",
        "results_tsv": "",
        "tax_id": "",
        "reject_divisions": [],
    },
}


def deep_update(base: dict, override: dict) -> dict:
    """Recursively merge *override* into *base*, returning a new dict."""
    out = dict(base)
    for k, v in (override or {}).items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = deep_update(out[k], v)
        else:
            out[k] = v
    return out


def load_yaml(path: str) -> dict:
    if not path:
        return {}
    try:
        import yaml
    except ImportError:
        raise SystemExit("PyYAML is required. Install with: pip install pyyaml")
    with open(path, encoding="utf-8") as f:
        return yaml.safe_load(f) or {}


def load_config(path: str) -> dict:
    """Return DEFAULT_CONFIG merged with the YAML at *path* (if any)."""
    if not path:
        return dict(DEFAULT_CONFIG)
    return deep_update(DEFAULT_CONFIG, load_yaml(path))
