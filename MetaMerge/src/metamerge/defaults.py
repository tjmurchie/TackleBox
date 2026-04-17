"""Default configuration for the MetaMerge package.

All thresholds and settings listed here can be overridden by a user-supplied
YAML config file passed to the CLI with ``--config``.  The keys and structure
must be preserved in any override file; only the values need to change.

Threshold rationale
-------------------
damage_min (0.01)
    Minimum metaDMG/Holi *damage* estimate to call a taxon exact-damage-supported.
    metaDMG reports damage as the C→T substitution rate at position 1 of reads
    mapping to that taxon.  Even 1% is a meaningful signal when combined with a
    high significance score.

significance_min (2.0)
    Minimum metaDMG/Holi *significance* (z-score of the damage estimate).
    Combined with damage_min this filters out taxa whose damage estimate is
    non-zero only because of noise.

high_confidence_n_reads_min (100)
    Minimum metaDMG N_reads for the "Very high confidence" tier.
    Requires at least 100 reads to have contributed to the damage estimate,
    ensuring the estimate is statistically robust.

strong_count_min_reads (50)
    Minimum MEGAN count in at least one non-control library for the "strong
    count support" flag.

strong_count_min_libraries (2)
    Minimum number of non-control libraries that must be positive for the
    "strong count support" flag.  Requiring presence in ≥2 libraries filters
    single-library noise.

blank_absolute_min (10)
    If any blank library has ≥10 MEGAN reads for a taxon, blank concern is
    flagged regardless of the real-library count.

blank_relative_min (0.10)
    If the maximum blank count / maximum real count ≥ 10%, blank concern is
    flagged.  Together with blank_absolute_min this allows a proportional
    assessment.

qc_alignments_per_read_clean_max (10.0)
    N_alignments / N_reads below this ratio is considered "clean" from a
    multimapping perspective.

qc_alignments_per_read_caution_max (20.0)
    Above this ratio the mapping is flagged as "strong caution" because
    excessive multimapping suggests the reads may not be genuinely placed.

qc_abs_rho_clean_max (0.30)
    |rho_Ac| below this value is considered a clean fit correlation.

qc_abs_rho_caution_max (0.50)
    |rho_Ac| above this value triggers "strong caution" — the damage model fit
    is poorly correlated with the data.

lineage_max_steps (2)
    Maximum number of rank steps between a focal taxon and a lineage-supporting
    candidate taxon.  Prevents long-range lineage support (e.g., a class-level
    hit supporting a species).

lineage_max_rank_level ("family")
    Broadest rank allowed for lineage support.  Taxa at order level or above
    are too broad to provide meaningful support for specific taxa.

tentative_min_reads (5)
    Minimum max real count for the "Tentative" category.

weak_support_min_reads (1)
    Minimum max real count to be reported at all (above 0).
"""

DEFAULT_CONFIG = {
    "io": {
        # Sheet name for Excel MEGAN count matrices; None = first sheet.
        "counts_sheet": None,
        # Sheet name for Excel metadata files; None = first sheet.
        "metadata_sheet": None,
        # Prefix for all output files.
        "output_prefix": "metamerge",
    },

    # Column-name aliases allow the package to work with different export
    # formats without requiring users to rename their columns.
    # Each list is searched in order; the first match wins.
    "column_aliases": {
        # MEGAN typically exports tax_id as '#Datasets' in TSV count matrices.
        "tax_id": [
            "tax_id", "taxon_id", "TaxID", "Taxon ID", "NCBI_taxid", "#Datasets",
        ],
        "tax_name": [
            "tax_name", "taxon_name", "Taxon-node", "Taxon-node ", "Taxon",
            "Name", "Scientific_name",
        ],
        "tax_rank": [
            "tax_rank", "rank", "Taxon_rank", "Taxonomic_rank",
        ],
        "sample": ["sample", "library", "Sample"],
        # Library-linker columns (must appear in the metadata/linker CSV).
        "megan_library_name": ["megan_library_name"],
        "holi_library_name":  ["holi_library_name"],
        "merged_library_name": ["merged_library_name"],
        # Optional metadata columns carried through to the workbook.
        "sample_id":          ["sample_id"],
        "sample_type":        ["sample_type", "type"],
        "is_negative_control": ["is_negative_control", "negative_control"],
        "is_positive_control": ["is_positive_control", "positive_control"],
        "site":   ["site"],
        "age":    ["age"],
        "depth":  ["depth"],
        "group":  ["group"],
        "notes":  ["notes"],
    },

    # Columns that MetaMerge loads from the Holi/metaDMG CSV.
    # All must be present; the loader will raise a clear error if any are missing.
    "holi_required_columns": [
        "sample",
        "tax_id",
        "tax_name",
        "tax_rank",
        "N_reads",
        "N_alignments",
        "damage",
        "significance",
        "rho_Ac",
        "MAP_valid",
        "tax_path",
    ],

    "thresholds": {
        "damage_min":                       0.01,
        "significance_min":                 2.0,
        "high_confidence_n_reads_min":      100,
        "strong_count_min_reads":           50,
        "strong_count_min_libraries":       2,
        "blank_absolute_min":               10,
        "blank_relative_min":               0.10,
        "qc_alignments_per_read_clean_max": 10.0,
        "qc_alignments_per_read_caution_max": 20.0,
        "qc_abs_rho_clean_max":             0.30,
        "qc_abs_rho_caution_max":           0.50,
        "lineage_max_steps":                2,
        "lineage_max_rank_level":           "family",
        "tentative_min_reads":              5,
        "weak_support_min_reads":           1,
    },

    "lineage": {
        # Set to False to disable the lineage-support check entirely.
        "enabled": True,
        # Taxon names that are too broad to provide meaningful lineage support.
        # Add project-specific over-broad names here.
        "forbid_broad_names": [
            "root", "cellular organisms", "Eukaryota", "Bacteria", "Archaea",
            "Bilateria", "Deuterostomia", "Protostomia", "Amniota", "Mammalia",
            "Theria", "Eutheria", "Placentalia", "Boreoeutheria", "Euarchontoglires",
            "Laurasiatheria", "Chordata", "Metazoa", "Animalia",
        ],
        # Rank hierarchy from most specific to most general.
        # This list controls the level arithmetic used by lineage support.
        "rank_levels": [
            "no rank", "clade", "subspecies", "species", "species group",
            "subgenus", "genus", "subtribe", "tribe", "subfamily", "family",
            "superfamily", "infraorder", "suborder", "order", "superorder",
            "class", "phylum", "kingdom", "domain", "superkingdom",
        ],
    },

    "taxonomy": {
        # Set to True to allow best-effort online GBIF common-name lookups.
        # The workflow is fully functional without internet access.
        "online_common_names": False,
        # IETF language tag for preferred common-name language.
        "common_name_language": "eng",
        # Small built-in convenience map; expand as needed.
        # Keys are scientific names (exact string match after normalize_name).
        "builtin_common_name_map": {
            "Homo":              "human",
            "Bison":             "bison",
            "Equus":             "horse",
            "Puma":              "puma",
            "Mammuthus":         "mammoth",
            "Urocitellus":       "ground squirrel",
            "Lepus":             "hare",
            "Rangifer tarandus": "caribou",
            "Betula":            "birch",
            "Picea":             "spruce",
            "Vaccinium":         "blueberry",
            "Salix":             "willow",
            "Achillea":          "yarrow",
            "Megalonyx":         "ground sloth",
            "Ondatra":           "muskrat",
            "Castor":            "beaver",
            "Oncorhynchus":      "Pacific salmon",
            "Ceratophyllum":     "hornwort",
            "Typha":             "cattail",
            "Phragmites":        "common reed",
        },
    },

    "report": {
        # Whether to auto-render heatmaps after a run (requires Rscript on PATH).
        "generate_heatmaps": False,
        # Maximum number of taxa per broad group in heatmap tables.
        "top_n": 30,
        # Most inclusive (broadest) taxonomic rank included in heatmap tables.
        "max_rank_for_reports": "family",
        # Minimum status included in heatmap tables.
        "min_status_for_reports": "Supported",
        # Whether to include negative-control libraries as columns in plots.
        "include_negative_controls": True,
        # Whether to include positive-control libraries as columns in plots.
        "include_positive_controls": True,
    },
}
