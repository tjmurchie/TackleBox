from __future__ import annotations

from dataclasses import dataclass, asdict, replace as dc_replace
from typing import Dict, Iterable, List, Optional, Tuple, Set
import json
import math

from .io import Hit
from .taxonomy import Taxonomy

# Rank levels for cross-clade detection.  Mirrors builddb._RANK_LEVEL (subset).
# Used to determine whether the LCA of a candidate hit with the best hit is above
# the order rank (level 16), which triggers the cross_clade_hit_penalty.
_RANK_LEVEL: Dict[str, int] = {
    "species": 3, "genus": 7, "family": 11,
    "order": 16, "superorder": 17,
    "class": 20, "superclass": 21,
    "phylum": 23, "kingdom": 26, "superkingdom": 27,
}

# NCBI superkingdom / kingdom anchor taxids for per-taxon filter dispatch.
# These match the same constants in builddb.py — kept separate to avoid a
# circular import between the two modules.
_SK_BACTERIA     = "2"
_SK_ARCHAEA      = "2157"
_SK_EUKARYOTA    = "2759"
_SK_METAZOA      = "33208"   # Animals
_SK_VIRIDIPLANTAE = "33090"  # Plants
_SK_FUNGI        = "4751"
_SK_VIRUSES      = "10239"

# Kingdoms for which stk/SES stacking filters are safe.  Eukaryotes using
# organelle/barcode references have high stk by construction — not a filter signal.
_STACKING_FILTER_KINGDOMS = frozenset({"bacteria", "archaea", "viruses", "unknown"})
# Kingdoms where short-reference filters (min_reads_fungi, damage-weight down-scaling)
# apply: ITS/18S references are too short for stk/SES to discriminate.
_SHORT_REF_KINGDOMS = frozenset({"fungi", "eukaryota"})

# Detection order: most specific first.
_KINGDOM_DETECT_ORDER: List[Tuple[str, str]] = [
    ("animals",   _SK_METAZOA),
    ("plants",    _SK_VIRIDIPLANTAE),
    ("fungi",     _SK_FUNGI),
    ("bacteria",  _SK_BACTERIA),
    ("archaea",   _SK_ARCHAEA),
    ("eukaryota", _SK_EUKARYOTA),
    ("viruses",   _SK_VIRUSES),
]


def _taxon_kingdom(taxid: str, taxonomy: "Taxonomy") -> str:
    """Return the broad kingdom for *taxid* using the Taxonomy object.

    Returns one of: ``"bacteria"``, ``"archaea"``, ``"animals"``,
    ``"plants"``, ``"fungi"``, ``"eukaryota"`` (other/protists),
    ``"viruses"``, or ``"unknown"``.
    """
    ancestor_set = set(taxonomy.ancestors(taxid))
    for kingdom, anchor in _KINGDOM_DETECT_ORDER:
        if anchor in ancestor_set:
            return kingdom
    return "unknown"


@dataclass
class ReadAssignment:
    read_id: str
    sample_id: str
    assigned_taxid: str
    assigned_name: str
    assigned_rank: str
    posterior: float
    entropy: float
    status: str
    reason: str
    n_hits_raw: int
    n_hits_used: int
    best_hit_taxid: str
    best_hit_name: str
    best_bitscore: float
    best_pident: float
    best_qcov: Optional[float]
    read_len: Optional[int] = None
    read_sequence: str = ""
    best_subject_id: str = ""
    best_subject_name: str = ""
    best_aln_len: int = 0
    best_qstart: Optional[int] = None
    best_qend: Optional[int] = None
    best_sstart: Optional[int] = None
    best_send: Optional[int] = None
    best_slen: Optional[int] = None
    second_path_taxid: str = ""
    second_path_name: str = ""
    second_path_posterior: float = 0.0
    lca_taxid: str = ""
    lca_name: str = ""
    alignment_evidence_json: str = "[]"
    # aDNA damage metrics (computed from qseq/sseq when available)
    damage_ct_5p: float = 0.0   # C→T rate at 5' terminus
    damage_ga_3p: float = 0.0   # G→A rate at 3' terminus
    damage_score: float = 0.0   # combined damage score (mean of above)
    damage_computable: bool = False  # True when qseq/sseq was present — needed for unbiased mean
    # EM refinement flag
    em_refined: int = 0         # 1 if this assignment was updated by EM pass
    # Experimental Bayesian fields (populated by BayesianRelicLCAClassifier only)
    damage_log_lr: float = 0.0  # log LR: P(damage|ancient_taxon) / P(damage|modern)
    # Weighted mean uniqueness of hits that determined this assignment (1.0 = fully unique).
    # Near 0 → reads align only to conserved loci shared by many taxa (FP risk).
    # Only populated when a uniqueness index is loaded (else stays 1.0).
    mean_hit_uniqueness: float = 1.0


@dataclass
class TaxonSummary:
    taxid: str
    name: str
    rank: str
    sample_id: str
    direct_hard_reads: int = 0
    direct_weighted_reads: float = 0.0
    cumulative_hard_reads: int = 0
    cumulative_weighted_reads: float = 0.0
    mean_posterior: float = 0.0
    conflicted_reads: int = 0
    negative_weighted_reads: float = 0.0
    blank_fraction: float = 0.0
    evidence_reads: int = 0
    unique_references: int = 0
    best_reference_breadth: float = 0.0
    mean_reference_breadth: float = 0.0
    best_reference_span: int = 0   # max(send)-min(sstart) on the most-covered reference; 0 when no coords
    max_stack_depth: int = 0
    top_locus_fraction: float = 0.0
    # stack_concentration: peak stack depth divided by total evidence reads.
    # Values near 1 → nearly all reads pile on one locus (suspicious).
    # Values near 0 → reads are spread across the reference (credible).
    # Robust to multi-reference scenarios; works even when slen is absent.
    stack_concentration: float = 0.0
    # Phase 1 stack-quality metrics (reference-length-aware)
    # coverage_gini: Gini coefficient of 50 bp-bin depth profile across all references.
    # 0 = reads perfectly spread; ~1 = all reads at one locus (FP fingerprint).
    # Includes zero bins for uncovered regions when reference length is known.
    coverage_gini: float = 0.0
    # stack_excess_score (SES): how many times more concentrated is the peak read density
    # than expected under uniform coverage across the full reference length?
    # SES = max_50bp_bin * total_ref_length / (evidence_reads * 50)
    # SES ≈ 1 → uniform; SES >> 1 → stacking.
    stack_excess_score: float = 0.0
    # n_covered_windows: count of distinct (accession, 50bp-bin) pairs with ≥1 read.
    # Directly measures "multiple distinct genomic regions". Stacking FPs concentrate all
    # reads in 1-2 windows on one accession; genuine TPs spread across 10-1000s of windows.
    n_covered_windows: int = 0
    # mean_windows_per_ref: n_covered_windows / unique_references. Normalises by accession
    # count so multi-accession taxa are not penalised for adding new references.
    mean_windows_per_ref: float = 0.0
    # max_per_ref_ses: maximum per-accession SES across all reference accessions.
    # Global SES can be diluted when a stacked accession is outnumbered by clean ones.
    # Per-ref SES = (max_50bp_bin_for_ref * slen) / (ref_n_reads * 50).
    max_per_ref_ses: float = 0.0
    # Total length (bp) of all reference sequences for this taxon (sum of accession slens).
    effective_ref_length: int = 0
    # Composite authenticity score (0-1, higher = more likely genuine TP).
    # Phase 1 weights; will be calibrated in Phase 2 against benchmark ground truth.
    composite_authenticity: float = 0.0
    # aDNA damage summary
    mean_damage_score: float = 0.0
    max_damage_score: float = 0.0
    flags: str = ""
    # Normalized abundance (Item 17)
    reads_per_million: float = 0.0     # RPM (hard reads / total input reads × 1e6)
    weighted_per_million: float = 0.0  # weighted reads / total input reads × 1e6
    relative_abundance: float = 0.0   # cumulative_weighted / total_assigned_weighted
    # Multi-proxy authenticity support (Item 3)
    eco_support: bool = False   # taxon in regional reference DB
    pal_support: bool = False   # matches palynology table
    fos_support: bool = False   # matches fossil occurrence table
    fos_evidence_text: str = ""  # semicolon-joined evidence strings from fossil records
    n_support_lines: int = 0    # total independent lines of positive evidence
    authenticity_tier: int = 0  # 1 (strongest) … 5 (present only); 0 = rejected/flagged
    authenticity_badge: str = ""  # ★ ● ◆ ▲ ○ ✕ or ""  (T1–T5, T0)
    authenticity_tier_pch: int = 0  # R/ggplot2 pch value: T0→4, T1→8, T2→16, T3→18, T4→17, T5→1
    # Evidence quality tier: HIGH / MEDIUM / LOW (computed in summarize_assignments).
    # HIGH: reads>=10, stk<0.05, n_covered_windows>=5 (or no alignment data), no warning flags.
    # MEDIUM: reads>=3, stk<0.15, n_covered_windows>=2 (or no alignment data).
    # LOW: everything else that survived hard filters.
    confidence_tier: str = ""
    # Stack-exclusion robustness metrics (populated when stack_exclusion_enabled=True).
    # non_stack_reads: reads with ≥1 alignment interval bin outside the identified stack zone.
    # stack_excluded_breadth: best reference breadth after removing stacked 50bp bins.
    # best_ref_length: length (bp) of best reference accession; 0 when slen unavailable.
    # stacking_status: "" = not assessed; "stacking_present_but_supported" = stacking
    #   detected but assignment holds without stacked reads; "stacking_fp" = collapsed
    #   (row removed from output); "short_reference_unassessable" = ref < min length.
    non_stack_reads: int = 0
    stack_excluded_breadth: float = 0.0
    best_ref_length: int = 0
    stacking_status: str = ""
    # Weighted mean uniqueness of hits classified to this taxon (posterior-weighted across reads).
    # Near 0 → assignments dominated by conserved-locus hits (FP risk, likely plant cp gene FP).
    # Near 1 → reads aligning to unique genomic regions (genuine TP signature).
    # Only populated when a uniqueness index is loaded (else stays 1.0).
    mean_hit_uniqueness: float = 1.0
    # Experimental Bayesian fields (populated by BayesianRelicLCAClassifier only)
    proportion_posterior_mean: float = 0.0  # Dirichlet posterior mean proportion
    proportion_ci_low: float = 0.0          # 95% credible interval lower bound
    proportion_ci_high: float = 0.0         # 95% credible interval upper bound
    estimated_damage_rate: float = 0.0      # estimated C→T rate for this taxon


def clamp(x: float, lo: float = 0.0, hi: float = 1.0) -> float:
    return max(lo, min(hi, x))


def entropy(probs: Iterable[float]) -> float:
    vals = [p for p in probs if p > 0]
    return -sum(p * math.log(p, 2) for p in vals)


def _safe_float(x: object, default: float = 0.0) -> float:
    try:
        if x is None or x == "":
            return default
        return float(x)
    except Exception:
        return default


def _safe_int(x: object, default: int = 0) -> int:
    try:
        if x is None or x == "":
            return default
        return int(float(x))
    except Exception:
        return default


def compute_damage_pattern(
    alignments: List[Tuple[str, str]],
    n_pos: int = 15,
) -> Dict[str, List[float]]:
    """Compute per-position C→T (5') and G→A (3') rates across many reads.

    Used for mapDamage-style visualisation.  Returns a dict:
        {
          "ct_5p": [rate_pos1, rate_pos2, ..., rate_posN],   # C→T at 5' terminus
          "ga_3p": [rate_pos1, rate_pos2, ..., rate_posN],   # G→A at 3' terminus
          "n_reads": int,                                    # reads contributing
        }

    Each element at index i gives the misincorporation frequency at position i+1
    from the respective terminus.  Positions with no reference C (5') or G (3') in
    any read have rate 0.0.

    Parameters
    ----------
    alignments : list of (qseq, sseq) pairs (BLAST alignment strings, no gaps handled)
    n_pos      : number of positions from each terminus to report (default 15)
    """
    ct_counts  = [0] * n_pos   # C→T observations at each 5' position
    ct_denoms  = [0] * n_pos   # C observations at each 5' position
    ga_counts  = [0] * n_pos   # G→A observations at each 3' position
    ga_denoms  = [0] * n_pos   # G observations at each 3' position
    n_reads    = 0

    for qseq, sseq in alignments:
        if not qseq or not sseq or len(qseq) != len(sseq):
            continue
        q = qseq.upper()
        s = sseq.upper()
        pairs = [(qi, si) for qi, si in zip(q, s) if qi != '-' and si != '-']
        if not pairs:
            continue
        n_reads += 1
        # 5' C→T
        for pos_idx, (qi, si) in enumerate(pairs[:n_pos]):
            if si == 'C':
                ct_denoms[pos_idx] += 1
                if qi == 'T':
                    ct_counts[pos_idx] += 1
        # 3' G→A (count from the 3' end inward)
        for pos_idx, (qi, si) in enumerate(reversed(pairs[-n_pos:])):
            if si == 'G':
                ga_denoms[pos_idx] += 1
                if qi == 'A':
                    ga_counts[pos_idx] += 1

    ct_rates = [c / d if d > 0 else 0.0 for c, d in zip(ct_counts, ct_denoms)]
    ga_rates = [c / d if d > 0 else 0.0 for c, d in zip(ga_counts, ga_denoms)]
    return {"ct_5p": ct_rates, "ga_3p": ga_rates, "n_reads": n_reads}


def compute_read_damage(
    qseq: str,
    sseq: str,
    strand: str = "plus",
    n_term: int = 5,
) -> Tuple[float, float, float]:
    """Compute classic aDNA deamination damage from BLAST qseq/sseq alignment strings.

    Returns (ct_5prime_rate, ga_3prime_rate, combined_score).

    Ancient DNA shows characteristic deamination at read termini:
    - 5' end: cytosine deamination → C→T substitutions (reference C, query T)
    - 3' end: G→A on the minus strand appears as G→A in plus-strand hits

    Gaps ('-') in qseq or sseq are skipped. Both strings must be equal length
    (as output by BLAST -outfmt 6 with qseq/sseq fields).

    Parameters
    ----------
    qseq:   Query sequence string from BLAST output
    sseq:   Subject (reference) sequence string from BLAST output
    strand: 'plus' or 'minus' — BLAST already reverse-complements minus-strand
            query sequences, so the 5'/3' orientation is preserved in both cases
    n_term: Number of terminal positions to examine (default 5)
    """
    if not qseq or not sseq or len(qseq) != len(sseq):
        return 0.0, 0.0, 0.0

    q = qseq.upper()
    s = sseq.upper()

    # Build list of aligned (query, subject) base pairs, skipping gaps
    pairs = [(qi, si) for qi, si in zip(q, s) if qi != '-' and si != '-']
    if not pairs:
        return 0.0, 0.0, 0.0

    # 5' end: reference C → query T
    five_pairs = pairs[:n_term]
    ct_denom = sum(1 for _q, _s in five_pairs if _s == 'C')
    ct_count = sum(1 for _q, _s in five_pairs if _s == 'C' and _q == 'T')
    ct_rate = ct_count / ct_denom if ct_denom > 0 else 0.0

    # 3' end: reference G → query A
    three_pairs = pairs[-n_term:]
    ga_denom = sum(1 for _q, _s in three_pairs if _s == 'G')
    ga_count = sum(1 for _q, _s in three_pairs if _s == 'G' and _q == 'A')
    ga_rate = ga_count / ga_denom if ga_denom > 0 else 0.0

    combined = (ct_rate + ga_rate) / 2.0
    return ct_rate, ga_rate, combined


def _best_hit_damage(hits: List[Hit], n_term: int = 5) -> Tuple[float, float, float, bool]:
    """Compute damage from the best hit that has qseq/sseq available.

    Returns (ct_5p, ga_3p, combined, computable) where computable=True means
    qseq was present and the rates are trustworthy (even if all zero).
    """
    for h in hits:
        if h.qseq and h.sseq:
            ct, ga, dmg = compute_read_damage(h.qseq, h.sseq, strand=h.sstrand or "plus", n_term=n_term)
            return ct, ga, dmg, True
    return 0.0, 0.0, 0.0, False


class RelicLCAClassifier:
    """Probabilistic/weighted LCA resolver for ancient metagenomic reads.

    RELIC-LCA (Read Evidence-weighted Likelihoods for Inferred Classification).
    Fillet does not try to make a brittle best-hit choice. It scores each candidate
    hit, aggregates evidence up the taxonomy, and walks down the tree only while
    one child clade remains stable enough to justify going deeper.

    v0.6 additions:
    - aDNA damage detection from qseq/sseq alignment strings
    - Sample-level EM refinement via classify_groups_with_em()
    """

    def __init__(
        self,
        taxonomy: Taxonomy,
        config: Dict,
        regional_taxa: Optional[Dict[str, Dict[str, str]]] = None,
        reference_metadata: Optional[Dict[str, Dict[str, str]]] = None,
        uniqueness_index=None,
    ):
        self.taxonomy = taxonomy
        self.cfg = config
        self.regional_taxa = regional_taxa or {}
        self.reference_metadata = reference_metadata or {}
        self._uniqueness_index = uniqueness_index
        self._cache_config()

    def _cache_config(self) -> None:
        """Cache config scalars as instance attrs to avoid repeated dict lookups in hot loops."""
        f = self.cfg.get("hit_filters", {})
        self._min_bitscore   = float(f.get("min_bitscore", 0))
        self._max_evalue     = float(f.get("max_evalue", 1e99))
        self._min_aln_len    = int(f.get("min_aln_len", 0))
        self._min_pident     = float(f.get("min_pident", 0))
        self._min_read_len   = int(f.get("min_read_len", 0) or 0)
        self._max_read_len   = int(f.get("max_read_len", 0) or 0)
        self._top_bitscore_frac    = float(f.get("top_bitscore_fraction", 0.80))
        self._max_hits_per_read    = int(f.get("max_hits_per_read", 5000))
        # Default 1: keep only the best hit per reference taxon before LCA scoring.
        # Multiple hits from the same taxon (e.g. 14 human loci for one conserved read)
        # are not independent evidence — summing them creates NT-abundance bias.
        self._max_hits_per_taxon   = int(f.get("max_hits_per_taxon", 1))
        self._min_qcov             = float(f.get("min_qcov", 0.0))
        w = self.cfg.get("weights", {})
        self._w_bitscore   = float(w.get("bitscore", 1.0))
        self._w_qcov       = float(w.get("query_coverage", 0.85))
        self._w_identity   = float(w.get("identity", 0.65))
        self._w_aln_len    = float(w.get("alignment_length", 0.30))
        self._w_read_len   = float(w.get("read_length", 0.20))
        self._w_regional   = float(w.get("regional_prior", 0.0))
        self._w_database   = float(w.get("database_weight", 0.0))
        self._w_em         = float(w.get("em_support", 0.30))
        self._w_coherence  = float(w.get("taxonomic_coherence", 0.60))
        self._w_isolated   = float(w.get("isolated_hit_penalty", 0.45))
        self._cross_clade_penalty = float(w.get("cross_clade_hit_penalty", 0.25))
        # rank level above which a hit's LCA with the plurality-winner is considered cross-clade.
        # With graduated penalties enabled this becomes the lowest level that receives any penalty.
        # Default 16 = order: cross-order and above get graduated suppression.
        self._cross_clade_rank_threshold = int(w.get("cross_clade_rank_threshold", 16))
        # Graduated cross-clade penalty: severity scales with phylogenetic distance.
        # True (default): use the factor table below instead of a single flat factor.
        # False: fall back to the original binary threshold + flat factor.
        self._cross_clade_graduated = bool(w.get("cross_clade_graduated", True))
        # Cap hits at the genus level: at most this many hits from any single genus (i.e.,
        # all species-level hits under the same genus share a single slot).  This prevents
        # NT over-representation of a genus (e.g. 14 Epinephelus species) from swamping
        # a correctly-matching taxon with fewer NT sequences.  0 = disabled.
        self._max_hits_per_genus = int(f.get("max_hits_per_genus", 1))
        # Soft genus cap: when True, hits that score within top_bitscore_frac% of the best
        # hit in the same genus bypass the hard per-genus cap.  This lets near-identical
        # sister species (e.g. B. amyloliquefaciens vs B. velezensis) both enter LCA scoring
        # so the read resolves to their common genus rather than arbitrarily picking one.
        # For reads at species-discriminating loci the sister falls below the tbf threshold
        # and is excluded normally — correct species assignment is preserved.
        self._genus_soft_cap_tbf = bool(f.get("genus_soft_cap_tbf", True))
        # Per-instance cache for genus-ancestor lookups (populated lazily during classification).
        if not hasattr(self, "_genus_cache"):
            self._genus_cache: Dict[str, str] = {}
        a = self.cfg.get("assignment", {})
        self._min_posterior_margin      = float(a.get("min_posterior_margin", 0.10))
        self._min_node_posterior        = float(a.get("min_node_posterior", 0.65))
        # Cache per-rank posterior thresholds into a dict for fast _walk_tree lookups.
        # Defaults are intentionally conservative at high taxonomic ranks to prevent
        # NT-database-coverage FPs: reads from conserved loci that BLAST primarily to
        # a non-target superorder/order (e.g. Laurasiatheria) because the target species
        # is underrepresented in NT for those loci.  Genuine detections have posteriors
        # well above these thresholds (typically 0.94+); NT-bias artifacts cluster at
        # 0.80-0.89.  Overridable per-project via config `assignment.min_posterior_<rank>`.
        self._rank_thresholds: Dict[str, float] = {
            "superorder": 0.90,
            "order": 0.85,
        }
        for k, v in a.items():
            if k.startswith("min_posterior_") and k not in ("min_posterior_margin",):
                rank_key = k[len("min_posterior_"):]
                self._rank_thresholds[rank_key] = float(v)
        self._em_reclass_threshold      = float(a.get("em_reclass_posterior_threshold", 0.80))
        self._max_entropy_species       = float(a.get("max_entropy_for_species", 1.80))
        self._max_entropy_genus         = float(a.get("max_entropy_for_genus", 2.20))
        self._max_entropy_family        = float(a.get("max_entropy_for_family", 3.00))
        self._max_entropy_order         = float(a.get("max_entropy_for_order", 3.50))
        # Sample-coherence lift: after EM, reads assigned to taxa with very low sample
        # support relative to the dominant clade are lifted to their common ancestor with
        # the dominant taxon (or to unclassified if the LCA is at/above order level).
        # min_em_coherence_fraction=0 disables the lift.
        self._em_coherence_lift_fraction = float(a.get("min_em_coherence_fraction", 0.0))
        # If the LCA of the assigned taxon with the dominant taxon is at or above this
        # rank level, set the read to unclassified rather than assigning to a very broad
        # ancestor.  Default: class (20) → reads lifted to class+ go to unclassified.
        # Set to 0 to always assign to LCA (never go unclassified via this path).
        self._em_coherence_max_lca_level = int(a.get("em_coherence_max_lca_level", 20))
        # Hard rank cap: after all EM and coherence-lift processing, any assignment whose
        # rank level (walking up through clade/no-rank nodes) exceeds this value is
        # demoted to unclassified.  This eliminates broad-ancestor FPs (e.g. Boreoeutheria,
        # Rodentia) for samples expected to contain only a single target taxon.
        # 0 = disabled (default).  11 = family-level cap (allows genus, tribe, family).
        self._max_assignment_rank_level = int(a.get("max_assignment_rank_level", 0))
        # Hit-diversity boost: raise walk-tree thresholds when candidates span many classes.
        self._hit_diversity_alpha       = float(w.get("hit_diversity_alpha", 0.10))
        self._hit_diversity_min_classes = int(w.get("hit_diversity_min_classes", 1))
        v = self.cfg.get("viewer", {})
        self._max_detail_hits           = int(v.get("max_alignment_hits_per_read", 80))
        d = self.cfg.get("damage", {})
        self._damage_n_term             = int(d.get("n_terminal_bases", 5))
        # Damage-aware scoring: "auto", "ancient", or "modern".
        self._damage_mode               = str(d.get("damage_mode", "auto"))
        self._ancient_damage_threshold  = float(d.get("ancient_damage_threshold", 0.05))
        # Bits to add back per C→T (5') or G→A (3') damage mismatch when computing
        # damage-adjusted bitscores for candidate filtering.  Prevents correct-taxon
        # alignments from being excluded by the top_bitscore_fraction threshold solely
        # because terminal damage mismatches lower their raw bitscore relative to a
        # competing taxon that naturally carries T/A at those positions.
        # Default 8.0 ≈ LAST transition mismatch penalty + lost match score for DNA.
        self._damage_bitscore_correction = float(d.get("bitscore_correction_per_mismatch", 8.0))
        # Per-sample antiquity state (populated at runtime by classify_groups_with_em).
        if not hasattr(self, "_ancient_samples"):
            self._ancient_samples: Set[str] = set()

    # ── Graduated cross-clade factor table ────────────────────────────────────
    # Keys are the *inclusive upper bound* of lca_level for that band.
    # Reads: "lca_level ≤ key → return this factor".
    _GRADUATED_FACTORS: List[Tuple[int, float]] = [
        (11, 1.00),   # ≤ family: same order or closer — no penalty
        (16, 0.60),   # ≤ order: cross-order within same class (e.g. Rodentia vs Lagomorpha)
        (20, 0.25),   # ≤ class: cross-class within phylum (e.g. mammal vs reptile)
        (23, 0.08),   # ≤ phylum: cross-phylum within kingdom (mammal vs fish in Chordata — KEY CASE)
        (99, 0.02),   # above: cross-kingdom / superkingdom (animal vs plant)
    ]

    def _graduated_cross_clade_factor(self, lca_level: int) -> float:
        """Return a graduated suppression factor for a hit whose LCA with the
        plurality-winner falls at *lca_level*.  Falls back to the original flat
        factor if graduated mode is disabled."""
        if not self._cross_clade_graduated:
            return (
                self._cross_clade_penalty
                if lca_level > self._cross_clade_rank_threshold
                else 1.0
            )
        if lca_level < 0:
            return 1.0
        for cap, factor in self._GRADUATED_FACTORS:
            if lca_level <= cap:
                return factor
        return 0.02

    def _genus_ancestor(self, taxid: str) -> str:
        """Return the genus-level ancestor of *taxid*, or *taxid* itself when it
        is already at or above genus level.  Results are cached per-instance."""
        cached = self._genus_cache.get(taxid)
        if cached is not None:
            return cached
        _genus_level = _RANK_LEVEL["genus"]  # 7
        ancs = self.taxonomy.ancestors(taxid, include_self=True)  # root → leaf order
        result = taxid
        # Walk leaf→root: first ancestor whose rank level ≥ genus is the genus
        # (or the nearest coarser group if genus is missing from the lineage).
        for anc in reversed(ancs):
            lvl = _RANK_LEVEL.get(self.taxonomy.rank(anc), -1)
            if lvl >= _genus_level:
                result = anc
                break
        self._genus_cache[taxid] = result
        return result

    def _hit_passes(self, h: Hit) -> bool:
        if h.bitscore < self._min_bitscore:
            return False
        if h.evalue > self._max_evalue:
            return False
        if h.aln_len < self._min_aln_len:
            return False
        if h.pident < self._min_pident:
            return False
        if self._min_qcov > 0 and h.query_coverage is not None and h.query_coverage < self._min_qcov:
            return False
        if self._min_read_len and h.qlen is not None and h.qlen < self._min_read_len:
            return False
        if self._max_read_len and h.qlen is not None and h.qlen > self._max_read_len:
            return False
        return self.taxonomy.exists(h.taxid)

    def _regional_prior(self, h: Hit) -> float:
        if not self.regional_taxa:
            return 1.0
        keys = set(self.taxonomy.ancestors(h.taxid, include_self=True))
        keys.add(self.taxonomy.name(h.taxid))
        for k in keys:
            if k in self.regional_taxa:
                row = self.regional_taxa[k]
                raw = row.get("weight") or row.get("prior") or row.get("regional_weight") or "1.25"
                try:
                    return max(0.05, float(raw))
                except Exception:
                    return 1.25
        return 1.0

    def _database_prior(self, h: Hit) -> float:
        if not self.reference_metadata:
            return 1.0
        row = self.reference_metadata.get(h.subject_id) or self.reference_metadata.get(h.subject_id.split(".", 1)[0])
        if not row:
            return 1.0
        raw = row.get("weight") or row.get("reference_weight") or row.get("quality_weight")
        if raw is None:
            return 1.0
        try:
            return max(0.01, float(raw))
        except Exception:
            return 1.0

    def _damage_corrected_pident(self, h: Hit) -> float:
        """Return pident with terminal C→T (5') and G→A (3') ignored as authentic deamination.

        Deamination mismatches counted at the first/last n_terminal_bases positions are
        added back to the match count, so ancient reads are not unfairly penalised for
        authentic damage signatures.  Falls back to raw pident when qseq/sseq unavailable.
        """
        if not h.qseq or not h.sseq or len(h.qseq) != len(h.sseq):
            return h.pident
        q = h.qseq.upper()
        s = h.sseq.upper()
        pairs = [(qi, si) for qi, si in zip(q, s) if qi != '-' and si != '-']
        if not pairs:
            return h.pident
        n_term = self._damage_n_term
        damage_mismatches = 0
        for qi, si in pairs[:n_term]:
            if si == 'C' and qi == 'T':
                damage_mismatches += 1
        for qi, si in pairs[-n_term:]:
            if si == 'G' and qi == 'A':
                damage_mismatches += 1
        aln_len = max(h.aln_len, 1)
        original_matches = (h.pident / 100.0) * aln_len
        corrected_matches = min(float(aln_len), original_matches + damage_mismatches)
        return corrected_matches / aln_len * 100.0

    def _apply_uniqueness_weights(self, hits: List[Hit]) -> None:
        """Set h.uniqueness_weight for each hit using the loaded uniqueness index.

        No-op when no index is loaded (weights remain 1.0).  Uses the per-50bp-window
        mean uniqueness for the [sstart, send] interval of each hit.
        """
        if self._uniqueness_index is None:
            return
        from fillet.builddb import lookup_uniqueness
        for h in hits:
            if h.sstart is not None and h.send is not None:
                h.uniqueness_weight = lookup_uniqueness(
                    self._uniqueness_index, h.subject_id, h.sstart, h.send
                )

    def _damage_adjusted_bitscore(self, h: Hit) -> float:
        """Return bitscore boosted by estimated terminal-damage mismatch correction.

        C→T (5') and G→A (3') mismatches at read termini are expected aDNA artefacts.
        These mismatches lower the LAST/BLAST bitscore of the correct-taxon alignment,
        potentially excluding it from the candidate set when a competing taxon naturally
        carries T (or A) at those positions and therefore scores higher.  Estimating the
        bitscore the hit would have without those penalties allows the candidate-filtering
        threshold to re-include alignments that are only below the bar because of damage.
        Falls back to raw bitscore when alignment strings are unavailable.
        """
        if not h.qseq or not h.sseq or len(h.qseq) != len(h.sseq):
            return h.bitscore
        q = h.qseq.upper()
        s = h.sseq.upper()
        pairs = [(qi, si) for qi, si in zip(q, s) if qi != '-' and si != '-']
        if not pairs:
            return h.bitscore
        n_term = self._damage_n_term
        n_damage = 0
        for qi, si in pairs[:n_term]:
            if si == 'C' and qi == 'T':
                n_damage += 1
        for qi, si in pairs[-n_term:]:
            if si == 'G' and qi == 'A':
                n_damage += 1
        if n_damage == 0:
            return h.bitscore
        return h.bitscore + n_damage * self._damage_bitscore_correction

    def _compute_diversity_boost(self, anc_map: Dict[str, List[str]]) -> float:
        """Return a threshold boost proportional to the phylogenetic spread of candidate hits.

        Conserved-region reads produce hits across many distantly related organisms (fungi,
        fish, mammals).  Boosting the walk-tree thresholds in proportion to the number of
        distinct classes prevents over-specific assignment when evidence is phylogenetically
        ambiguous.  Boost = alpha × max(0, n_distinct_classes − min_classes).
        """
        if self._hit_diversity_alpha <= 0 or not anc_map:
            return 0.0
        classes: Set[str] = set()
        for taxid in anc_map:
            cls = self.taxonomy.ranked_ancestor(taxid, ["class"])
            if cls:
                classes.add(cls)
        n = len(classes)
        return self._hit_diversity_alpha * max(0, n - self._hit_diversity_min_classes)

    def _base_score(self, h: Hit, max_bits: float, em_support: Optional[Dict[Tuple[str, str], float]] = None, sample_id: str = "", anc_cache: Optional[Dict[str, List[str]]] = None, is_ancient: bool = False) -> float:
        bits_factor = clamp((h.bitscore * h.uniqueness_weight) / max(max_bits, 1.0))
        qcov = h.query_coverage
        qcov_factor = 1.0 if qcov is None else clamp(qcov)
        pident = self._damage_corrected_pident(h) if is_ancient else h.pident
        id_factor = clamp(pident / 100.0)
        len_factor = clamp(math.sqrt(max(h.aln_len, 0) / 80.0))
        read_len_factor = 1.0
        if h.qlen:
            # Long enough ancient reads are useful, but do not make longer modern-ish reads dominate.
            read_len_factor = clamp(math.sqrt(max(h.qlen, 0) / 70.0), 0.15, 1.35)
        eps = 1e-9
        score = (bits_factor + eps) ** self._w_bitscore
        score *= (qcov_factor + eps) ** self._w_qcov
        score *= (id_factor + eps) ** self._w_identity
        score *= (len_factor + eps) ** self._w_aln_len
        score *= (read_len_factor + eps) ** self._w_read_len
        if self._w_regional:
            score *= self._regional_prior(h) ** self._w_regional
        if self._w_database:
            score *= self._database_prior(h) ** self._w_database
        # EM prior: multiply by sample-level support for this taxon (normalized)
        if em_support and sample_id and self._w_em > 0:
            ancs = anc_cache[h.taxid] if anc_cache else self.taxonomy.ancestors(h.taxid, include_self=True)
            lineage_support = max(
                (em_support.get((sample_id, anc), 0.0) for anc in ancs),
                default=0.0,
            )
            # Soft prior: scale from 1.0 (no support) to 2.0 (strong support).
            # We use log(1 + support) to avoid runaway feedback on already dominant taxa.
            em_factor = 1.0 + math.log1p(max(lineage_support, 0.0)) * 0.25
            score *= em_factor ** self._w_em
        return score

    def _candidate_hits(self, hits: List[Hit], is_ancient: bool = False) -> List[Hit]:
        passed = [h for h in hits if self._hit_passes(h)]
        if not passed:
            return []
        passed.sort(key=lambda x: x.bitscore, reverse=True)
        if is_ancient:
            # Use damage-adjusted bitscores for threshold computation so that hits
            # whose raw bitscore is suppressed only by terminal C→T/G→A damage
            # mismatches are not excluded from the candidate set.  The adjusted score
            # is used for filtering only; downstream scoring still uses raw bitscore.
            adj_bits = [self._damage_adjusted_bitscore(h) for h in passed]
            max_bits = max(adj_bits)
            threshold = max_bits * self._top_bitscore_frac
            candidates = [h for h, ab in zip(passed, adj_bits) if ab >= threshold]
        else:
            max_bits = passed[0].bitscore
            threshold = max_bits * self._top_bitscore_frac
            candidates = [h for h in passed if h.bitscore >= threshold]
        if self._max_hits_per_taxon > 0:
            per_taxon: Dict[str, int] = {}
            deduped = []
            for h in candidates:
                n = per_taxon.get(h.taxid, 0)
                if n < self._max_hits_per_taxon:
                    deduped.append(h)
                    per_taxon[h.taxid] = n + 1
            candidates = deduped
        # Genus-level cap: collapse species-level hits from the same genus.
        # Hard cap: at most max_hits_per_genus representatives per genus (prevents
        # NT-abundance bias from genera with 100s of sequences, e.g. 14 Epinephelus spp.).
        # Soft tie-break (genus_soft_cap_tbf=True, default): hits that score within
        # top_bitscore_frac% of the best hit in the same genus bypass the hard cap.
        # This resolves near-identical sister pairs (e.g. B. amyloliquefaciens vs
        # B. velezensis) to their common genus via LCA rather than arbitrarily picking one.
        # At species-discriminating loci the sister falls below the tbf threshold and is
        # excluded normally, so genuine species-level assignments are preserved.
        if self._max_hits_per_genus > 0:
            # First pass: find best bitscore within each genus
            genus_best_bits: Dict[str, float] = {}
            for h in candidates:
                genus = self._genus_ancestor(h.taxid)
                if genus not in genus_best_bits or h.bitscore > genus_best_bits[genus]:
                    genus_best_bits[genus] = h.bitscore
            per_genus: Dict[str, int] = {}
            genus_deduped = []
            for h in candidates:
                genus = self._genus_ancestor(h.taxid)
                n = per_genus.get(genus, 0)
                within_hard_cap = n < self._max_hits_per_genus
                # Soft tie-break: admit hits tied with the genus best even past the hard cap
                within_genus_tbf = (
                    self._genus_soft_cap_tbf
                    and h.bitscore >= genus_best_bits[genus] * self._top_bitscore_frac
                )
                if within_hard_cap or within_genus_tbf:
                    genus_deduped.append(h)
                    per_genus[genus] = n + 1
            candidates = genus_deduped
        return candidates[:self._max_hits_per_read]

    def _score_hits(self, hits: List[Hit], em_support: Optional[Dict[Tuple[str, str], float]] = None, sample_id: str = "", is_ancient: bool = False) -> Tuple[List[float], Dict[str, List[str]]]:
        if not hits:
            return [], {}
        max_bits = max(h.bitscore * h.uniqueness_weight for h in hits)
        # Precompute ancestors once per unique taxid to avoid redundant cache lookups
        anc_map: Dict[str, List[str]] = {}
        for h in hits:
            if h.taxid not in anc_map:
                anc_map[h.taxid] = self.taxonomy.ancestors(h.taxid, include_self=True)
        # Use a list instead of dict — keys are 0..n-1 so list indexing is faster.
        base: List[float] = [self._base_score(h, max_bits, em_support=em_support, sample_id=sample_id, anc_cache=anc_map, is_ancient=is_ancient) for h in hits]
        total = sum(base) or 1.0
        # Group base scores by taxid before propagating to ancestors.
        # Same numeric result as the per-hit loop, but O(unique_taxids × lineage)
        # instead of O(all_hits × lineage) — a large win when many hits share a taxid.
        base_by_taxid: Dict[str, float] = {}
        for i, h in enumerate(hits):
            base_by_taxid[h.taxid] = base_by_taxid.get(h.taxid, 0.0) + base[i]
        lineage_support: Dict[str, float] = {}
        for taxid, bscore in base_by_taxid.items():
            for anc in anc_map[taxid]:
                lineage_support[anc] = lineage_support.get(anc, 0.0) + bscore
        coherence_w = self._w_coherence
        isolated_penalty = self._w_isolated
        # Precompute normalised support to avoid per-hit divisions inside max().
        support_norm: Dict[str, float] = {anc: v / total for anc, v in lineage_support.items()}
        # Coherence and parent_share are taxid-level properties — compute once per unique
        # taxid rather than once per hit.  With cap40 data and many duplicate-taxid hits,
        # this can reduce the coherence loop by 10-20x.
        coh_by_taxid: Dict[str, float] = {}
        parent_share_by_taxid: Dict[str, float] = {}
        for taxid, ancs in anc_map.items():
            informative = ancs[2:] if len(ancs) > 2 else ancs
            coh_by_taxid[taxid] = max([support_norm.get(a, 0.0) for a in informative]) if informative else 1.0
            parent = self.taxonomy.parent(taxid)
            parent_share_by_taxid[taxid] = support_norm.get(parent, 0.0) if parent else 1.0
        multi_hit = len(hits) > 1
        penalty_factor = max(0.05, 1.0 - isolated_penalty)

        # Cross-clade suppression: hits whose LCA with the plurality-winner taxid is at or
        # above the cross_clade_rank_threshold receive a graduated suppression factor.
        # Using the plurality winner (highest total base score) rather than the single
        # best-bitscore hit avoids penalising e.g. bovids when the single best hit is a
        # dolphin at marginally higher bitscore.
        # Graduated mode: factor scales continuously with phylogenetic distance:
        #   ≤ family (11)  → 1.00 (no penalty)
        #   ≤ order  (16)  → 0.60 (mild: cross-order within same class)
        #   ≤ class  (20)  → 0.25 (moderate: e.g. mammal vs reptile)
        #   ≤ phylum (23)  → 0.08 (strong: e.g. mammal vs fish in Chordata)
        #   above          → 0.02 (near-zero: cross-kingdom)
        ref_taxid = max(base_by_taxid, key=lambda t: base_by_taxid[t])
        ref_anc_set: frozenset = frozenset(anc_map[ref_taxid])
        cross_clade_weight_by_taxid: Dict[str, float] = {}
        for taxid, ancs in anc_map.items():
            if taxid == ref_taxid:
                cross_clade_weight_by_taxid[taxid] = 1.0
                continue
            # LCA = deepest ancestor shared with the plurality-winner (ancs is root-first, leaf-last)
            lca_txid = next((t for t in reversed(ancs) if t in ref_anc_set), None)
            if lca_txid is None:
                cross_clade_weight_by_taxid[taxid] = self._graduated_cross_clade_factor(99)
                continue
            # Walk up from lca_txid until we find a named Linnean rank.
            cur = lca_txid
            lca_level = -1
            for _ in range(20):
                lvl = _RANK_LEVEL.get(self.taxonomy.rank(cur), -1)
                if lvl >= 0:
                    lca_level = lvl
                    break
                p = self.taxonomy.parent(cur)
                if not p or p == cur:
                    lca_level = 99
                    break
                cur = p
            cross_clade_weight_by_taxid[taxid] = self._graduated_cross_clade_factor(lca_level)

        scores: List[float] = []
        for i, h in enumerate(hits):
            coh = coh_by_taxid[h.taxid]
            score = base[i] * ((0.5 + 0.5 * coh) ** coherence_w)
            if parent_share_by_taxid[h.taxid] < 0.25 and multi_hit:
                score *= penalty_factor
            score *= cross_clade_weight_by_taxid.get(h.taxid, 1.0)
            scores.append(score)
        return scores, anc_map

    def _node_support(self, hits: List[Hit], scores: List[float], anc_map: Optional[Dict[str, List[str]]] = None) -> Tuple[Dict[str, float], float]:
        node: Dict[str, float] = {}
        total = sum(scores) or 0.0
        if anc_map is None:
            anc_map = {}
            for h in hits:
                if h.taxid not in anc_map:
                    anc_map[h.taxid] = self.taxonomy.ancestors(h.taxid, include_self=True)
        # Group scores by taxid first, then propagate once per unique taxid.
        score_by_taxid: Dict[str, float] = {}
        for i, h in enumerate(hits):
            score_by_taxid[h.taxid] = score_by_taxid.get(h.taxid, 0.0) + scores[i]
        for taxid, s in score_by_taxid.items():
            for anc in anc_map[taxid]:
                node[anc] = node.get(anc, 0.0) + s
        return node, total

    def _threshold_for_rank(self, rank: str) -> float:
        return self._rank_thresholds.get(rank.replace(" ", "_"), self._min_node_posterior)

    def _walk_tree(self, node_support: Dict[str, float], total: float, diversity_boost: float = 0.0) -> Tuple[str, float, str, float, str]:
        margin_min = self._min_posterior_margin
        root = self.taxonomy.root_taxid if self.taxonomy.root_taxid in self.taxonomy.taxa else next(iter(self.taxonomy.taxa))
        cur = root
        reason_bits: List[str] = []
        while True:
            children = [c for c in self.taxonomy.taxa.get(cur).children if node_support.get(c, 0.0) > 0] if cur in self.taxonomy.taxa else []
            if not children:
                break
            ranked = sorted(children, key=lambda c: node_support.get(c, 0.0), reverse=True)
            best = ranked[0]
            best_post = node_support.get(best, 0.0) / total if total else 0.0
            second = ranked[1] if len(ranked) > 1 else ""
            second_post = node_support.get(second, 0.0) / total if second and total else 0.0
            margin = best_post - second_post
            thresh = self._threshold_for_rank(self.taxonomy.rank(best)) + diversity_boost
            if best_post >= thresh and margin >= margin_min:
                cur = best
                reason_bits.append(f"descended to {self.taxonomy.name(cur)} ({best_post:.3f}; margin {margin:.3f})")
            else:
                reason_bits.append(
                    f"stopped above {self.taxonomy.name(best)}: posterior {best_post:.3f}, margin {margin:.3f}, threshold {thresh:.3f}"
                )
                break
        post = node_support.get(cur, 0.0) / total if total else 0.0
        parent = self.taxonomy.parent(cur)
        siblings = [c for c in self.taxonomy.taxa.get(parent, self.taxonomy.taxa.get(cur)).children if c != cur and node_support.get(c, 0.0) > 0] if parent in self.taxonomy.taxa else []
        if siblings:
            sec = max(siblings, key=lambda c: node_support.get(c, 0.0))
            sec_post = node_support.get(sec, 0.0) / total if total else 0.0
        else:
            sec, sec_post = "", 0.0
        return cur, post, sec, sec_post, "; ".join(reason_bits)

    # Ranks at or above family are considered "discordant" for chimera detection.
    _DISCORDANT_RANKS: Set[str] = {
        "superkingdom", "kingdom", "phylum", "class", "order", "family",
    }

    def _detect_chimeric_segments(
        self,
        candidates: List[Hit],
        min_segment_len: int = 15,
        max_overlap_fraction: float = 0.30,
        bitscore_frac: float = 0.80,
    ) -> Optional[str]:
        """Return a status tag if top hits cover discordant, non-overlapping query segments.

        Chimeric reads arise from ligation artefacts: the 5' portion of a read maps
        confidently to taxon A while the 3' portion maps to an unrelated taxon B.
        BLAST reports these as separate hits with different qstart/qend ranges.

        Detection criteria (all must hold):
        - At least two top-scoring hits (bitscore ≥ bitscore_frac * max) have qstart/qend.
        - Their query intervals overlap by ≤ max_overlap_fraction of the shorter interval.
        - Each segment is ≥ min_segment_len bp.
        - Their LCA in the taxonomy tree is at family rank or above (truly discordant).

        Returns a string like "chimeric_segments:1-30(Homo sapiens)/35-65(Ursus arctos)@family"
        or None if no chimeric pattern is detected.
        """
        if not candidates:
            return None

        max_bits = candidates[0].bitscore  # already sorted by bitscore desc
        threshold = max_bits * bitscore_frac
        # Deduplicate by taxid: keep the hit with the widest query span per taxid.
        # Hits from the same taxon are never discordant, so only one representative
        # per taxid is needed.  This reduces O(n_hits²) to O(n_taxa²).
        best_per_taxid: Dict[str, Hit] = {}
        for h in candidates:
            if h.bitscore < threshold or h.qstart is None or h.qend is None:
                continue
            span = abs(h.qend - h.qstart)
            prev = best_per_taxid.get(h.taxid)
            if prev is None or span > abs(prev.qend - prev.qstart):  # type: ignore[operator]
                best_per_taxid[h.taxid] = h
        top = list(best_per_taxid.values())[:50]  # hard cap — one rep/taxon is enough
        if len(top) < 2:
            return None

        for i, ha in enumerate(top):
            a0 = min(ha.qstart, ha.qend)
            a1 = max(ha.qstart, ha.qend)
            len_a = a1 - a0 + 1
            if len_a < min_segment_len:
                continue
            for hb in top[i + 1:]:
                if hb.taxid == ha.taxid:
                    continue
                b0 = min(hb.qstart, hb.qend)
                b1 = max(hb.qstart, hb.qend)
                len_b = b1 - b0 + 1
                if len_b < min_segment_len:
                    continue
                overlap = max(0, min(a1, b1) - max(a0, b0) + 1)
                if overlap / min(len_a, len_b) > max_overlap_fraction:
                    continue  # segments overlap too much — not separate regions
                lca = self.taxonomy.lca([ha.taxid, hb.taxid]) or "0"
                lca_rank = self.taxonomy.rank(lca)
                if lca_rank not in self._DISCORDANT_RANKS:
                    continue  # same genus/species level disagreement, not a chimera
                name_a = self.taxonomy.name(ha.taxid)
                name_b = self.taxonomy.name(hb.taxid)
                return f"chimeric_segments:{a0}-{a1}({name_a})/{b0}-{b1}({name_b})@{lca_rank}"
        return None

    def _hit_details(self, hits: List[Hit], scores: List[float], max_hits: int) -> str:
        total = sum(scores) or 1.0
        rows = []
        for i, h in enumerate(hits[:max_hits]):
            interval = h.subject_interval
            # Per-hit damage if qseq/sseq available
            ct, ga, dmg = (0.0, 0.0, 0.0)
            if h.qseq and h.sseq:
                ct, ga, dmg = compute_read_damage(h.qseq, h.sseq, strand=h.sstrand or "plus")
            rows.append({
                "subject_id": h.subject_id,
                "subject_name": h.subject_name or "",
                "taxid": h.taxid,
                "taxon_name": self.taxonomy.name(h.taxid),
                "rank": self.taxonomy.rank(h.taxid),
                "bitscore": h.bitscore,
                "pident": h.pident,
                "evalue": h.evalue,
                "aln_len": h.aln_len,
                "mismatch": h.mismatch,
                "gapopen": h.gapopen,
                "qstart": h.qstart,
                "qend": h.qend,
                "qlen": h.qlen,
                "qcov": h.query_coverage,
                "sstart": h.sstart,
                "send": h.send,
                "slen": h.slen,
                "sstrand": h.sstrand or "",
                "qseq": h.qseq or "",
                "sseq": h.sseq or "",
                "btop": h.btop or "",
                "subject_interval": list(interval) if interval else None,
                "hit_weight": scores[i],
                "hit_posterior": scores[i] / total,
                "damage_ct_5p": round(ct, 4),
                "damage_ga_3p": round(ga, 4),
                "damage_score": round(dmg, 4),
            })
        return json.dumps(rows, ensure_ascii=False, separators=(",", ":"))

    def classify_read(
        self,
        read_id: str,
        hits: List[Hit],
        sample_id: str = "sample",
        sequence: str = "",
        em_support: Optional[Dict[Tuple[str, str], float]] = None,
        em_refined: bool = False,
    ) -> ReadAssignment:
        n_raw = len(hits)
        # Apply uniqueness weights before candidate filtering (weights don't affect
        # candidate selection — that still uses raw bitscore — only scoring).
        self._apply_uniqueness_weights(hits)
        # Determine antiquity before candidate filtering so that damage-adjusted
        # bitscores can be used to widen the candidate set for ancient reads.
        is_ancient = (
            self._damage_mode == "ancient"
            or (self._damage_mode == "auto" and sample_id in self._ancient_samples)
        )
        candidates = self._candidate_hits(hits, is_ancient=is_ancient)
        max_detail_hits = self._max_detail_hits
        if not candidates:
            read_len = len(sequence) if sequence else (hits[0].qlen if hits else None)
            return ReadAssignment(
                read_id=read_id,
                sample_id=sample_id,
                assigned_taxid="0",
                assigned_name="unclassified",
                assigned_rank="no rank",
                posterior=0.0,
                entropy=0.0,
                status="unclassified:no_candidate_hits",
                reason="No hits passed minimal junk filters or mapped to known taxonomy.",
                n_hits_raw=n_raw,
                n_hits_used=0,
                best_hit_taxid="",
                best_hit_name="",
                best_bitscore=0.0,
                best_pident=0.0,
                best_qcov=None,
                read_len=read_len,
                read_sequence=sequence,
                lca_taxid="",
                lca_name="",
            )
        candidates.sort(key=lambda h: h.bitscore, reverse=True)
        scores, anc_map = self._score_hits(candidates, em_support=em_support, sample_id=sample_id, is_ancient=is_ancient)
        node_support, total = self._node_support(candidates, scores, anc_map=anc_map)
        diversity_boost = self._compute_diversity_boost(anc_map)
        assigned, post, sec, sec_post, reason = self._walk_tree(node_support, total, diversity_boost=diversity_boost)
        leaf_probs: Dict[str, float] = {}
        for i, h in enumerate(candidates):
            leaf_probs[h.taxid] = leaf_probs.get(h.taxid, 0.0) + scores[i]
        probs = [v / total for v in leaf_probs.values()] if total else []
        ent = entropy(probs)
        unique_taxids = list(anc_map)  # anc_map keys are already unique taxids
        lca = self.taxonomy.lca(unique_taxids) or "0"
        status = "assigned"
        if len({self.taxonomy.ranked_ancestor(t, ["class", "order", "family"]) or t for t in unique_taxids}) > 1:
            status += ";conflicted_cross_clade"
        rank = self.taxonomy.rank(assigned)
        if rank == "species" and ent > self._max_entropy_species:
            status += ";high_entropy_species"
        if rank == "genus" and ent > self._max_entropy_genus:
            status += ";high_entropy_genus"
        if rank == "family" and ent > self._max_entropy_family:
            status += ";high_entropy_family"
        if rank == "order" and ent > self._max_entropy_order:
            status += ";high_entropy_order"
        chimera_tag = self._detect_chimeric_segments(candidates)
        if chimera_tag:
            status += f";{chimera_tag}"
        if em_refined:
            status += ";em_refined"
        best = candidates[0]
        read_len = len(sequence) if sequence else best.qlen
        # Compute aDNA damage from best-scoring hit that has qseq/sseq
        damage_n_term = self._damage_n_term
        ct_5p, ga_3p, dmg, dmg_computable = _best_hit_damage(candidates, n_term=damage_n_term)
        # Posterior-weighted mean uniqueness across candidate hits.
        # Reflects the average uniqueness of the alignments that drove this assignment.
        _uq_num = sum(scores[i] * candidates[i].uniqueness_weight for i in range(len(candidates)))
        mean_uq = _uq_num / total if total > 0 else 1.0
        return ReadAssignment(
            read_id=read_id,
            sample_id=sample_id,
            assigned_taxid=assigned,
            assigned_name=self.taxonomy.name(assigned),
            assigned_rank=rank,
            posterior=post,
            entropy=ent,
            status=status,
            reason=reason,
            n_hits_raw=n_raw,
            n_hits_used=len(candidates),
            best_hit_taxid=best.taxid,
            best_hit_name=self.taxonomy.name(best.taxid),
            best_bitscore=best.bitscore,
            best_pident=best.pident,
            best_qcov=best.query_coverage,
            read_len=read_len,
            read_sequence=sequence,
            best_subject_id=best.subject_id,
            best_subject_name=best.subject_name or "",
            best_aln_len=best.aln_len,
            best_qstart=best.qstart,
            best_qend=best.qend,
            best_sstart=best.sstart,
            best_send=best.send,
            best_slen=best.slen,
            second_path_taxid=sec,
            second_path_name=self.taxonomy.name(sec) if sec else "",
            second_path_posterior=sec_post,
            lca_taxid=lca,
            lca_name=self.taxonomy.name(lca),
            alignment_evidence_json=self._hit_details(candidates, scores, max_detail_hits),
            damage_ct_5p=ct_5p,
            damage_ga_3p=ga_3p,
            damage_score=dmg,
            damage_computable=dmg_computable,
            em_refined=int(em_refined),
            mean_hit_uniqueness=mean_uq,
        )

    def classify_groups(
        self,
        groups: Iterable[Tuple[str, List[Hit]]],
        sample_id: str,
        sequences: Optional[Dict[str, str]] = None,
    ) -> List[ReadAssignment]:
        sequences = sequences or {}
        return [
            self.classify_read(read_id, hits, sample_id=sample_id, sequence=sequences.get(read_id, ""))
            for read_id, hits in groups
        ]

    def classify_groups_with_em(
        self,
        groups: Iterable[Tuple[str, List[Hit]]],
        sample_id: str,
        sequences: Optional[Dict[str, str]] = None,
        em_iterations: int = 1,
    ) -> List[ReadAssignment]:
        """Classify reads with optional sample-level EM support refinement.

        EM pass logic (1-2 iterations recommended):
        1. First-pass classification (streaming — one read at a time).
        2. Build a taxon support map from all first-pass posteriors.
        3. Re-classify only ambiguous reads (posterior < threshold or conflicted)
           using the support map as an additional prior.

        Memory model: only hits for ambiguous reads are retained between passes.
        High-confidence reads (posterior >= threshold and not conflicted) have
        their hits discarded immediately after first-pass classification, so
        memory usage scales with the fraction of ambiguous reads rather than
        total input size.  This makes EM safe for very large samples.
        """
        sequences = sequences or {}
        em_reclass_threshold = self._em_reclass_threshold

        if em_iterations <= 0:
            # Pure streaming — no EM, lowest possible memory footprint.
            assignments = [
                self.classify_read(rid, hits, sample_id=sample_id, sequence=sequences.get(rid, ""))
                for rid, hits in groups
            ]
            return self._apply_rank_cap(assignments)

        # Streaming first pass: classify every read, retain hits only for reads
        # that might benefit from EM re-classification (ambiguous / conflicted).
        assignments: List[ReadAssignment] = []
        ambiguous_hits: Dict[str, List[Hit]] = {}

        for rid, hits in groups:
            a = self.classify_read(rid, hits, sample_id=sample_id, sequence=sequences.get(rid, ""))
            assignments.append(a)
            if (
                a.assigned_taxid != "0"
                and (a.posterior < em_reclass_threshold or "conflicted" in a.status)
            ):
                ambiguous_hits[rid] = hits
            # hits for high-confidence reads are discarded immediately

        # Phase 4: damage-aware antiquity detection for "auto" mode.
        # Compute mean damage score across first-pass assigned reads; if the sample
        # is ancient (mean_damage > threshold), mark it so the EM pass uses
        # damage-corrected pident for re-classification.
        if self._damage_mode == "auto" and sample_id not in self._ancient_samples:
            damage_vals = [a.damage_score for a in assignments if a.assigned_taxid != "0" and a.damage_computable]
            if damage_vals and (sum(damage_vals) / len(damage_vals)) > self._ancient_damage_threshold:
                self._ancient_samples.add(sample_id)

        for _iter in range(em_iterations):
            # Build sample support from all first-pass assignments.
            em_support: Dict[Tuple[str, str], float] = {}
            for a in assignments:
                if a.assigned_taxid == "0":
                    continue
                key: Tuple[str, str] = (sample_id, a.assigned_taxid)
                em_support[key] = em_support.get(key, 0.0) + a.posterior
                for anc in self.taxonomy.ancestors(a.assigned_taxid, include_self=False):
                    akey: Tuple[str, str] = (sample_id, anc)
                    em_support[akey] = em_support.get(akey, 0.0) + a.posterior * 0.5

            # Re-classify only the buffered ambiguous reads.
            new_assignments: List[ReadAssignment] = []
            for a in assignments:
                hits = ambiguous_hits.get(a.read_id)
                if hits is not None:
                    refined = self.classify_read(
                        a.read_id,
                        hits,
                        sample_id=sample_id,
                        sequence=sequences.get(a.read_id, a.read_sequence or ""),
                        em_support=em_support,
                        em_refined=False,
                    )
                    actually_changed = (
                        refined.assigned_taxid != a.assigned_taxid
                        or abs(refined.posterior - a.posterior) > 0.02
                    )
                    if actually_changed:
                        clean_status = refined.status.replace(";em_refined", "")
                        refined = dc_replace(
                            refined,
                            em_refined=1,
                            status=clean_status + ";em_refined",
                        )
                    new_assignments.append(refined)
                else:
                    new_assignments.append(a)
            assignments = new_assignments

            # Shrink the ambiguous buffer for subsequent iterations.
            if _iter < em_iterations - 1:
                ambiguous_hits = {
                    a.read_id: ambiguous_hits[a.read_id]
                    for a in assignments
                    if a.read_id in ambiguous_hits
                    and a.assigned_taxid != "0"
                    and (a.posterior < em_reclass_threshold or "conflicted" in a.status)
                }

        # Sample-coherence lift (optional, controlled by min_em_coherence_fraction).
        # After EM, reads assigned to taxa with very low direct assignment support
        # relative to the dominant taxon are lifted to their LCA with the dominant taxon.
        # Uses DIRECT posterior sums (not the EM-inflated ancestor em_support) to avoid
        # false positives from shared ancestors (e.g., squirrel support at Glires does
        # not protect Ochotona assignments from being lifted).
        if em_iterations > 0 and self._em_coherence_lift_fraction > 0:
            # Build a direct assignment support map: taxid → sum of posteriors.
            direct_support_map_lift: Dict[str, float] = {}
            for a in assignments:
                if a.assigned_taxid != "0":
                    direct_support_map_lift[a.assigned_taxid] = (
                        direct_support_map_lift.get(a.assigned_taxid, 0.0) + a.posterior
                    )

            if direct_support_map_lift:
                dom_lift_txid = max(direct_support_map_lift, key=lambda t: direct_support_map_lift[t])
                dom_lift_support = direct_support_map_lift[dom_lift_txid]
                lift_threshold = dom_lift_support * self._em_coherence_lift_fraction
                dom_lift_ancs: Set[str] = set(
                    self.taxonomy.ancestors(dom_lift_txid, include_self=True)
                )

                corrected: List[ReadAssignment] = []
                for a in assignments:
                    if a.assigned_taxid == "0":
                        corrected.append(a)
                        continue
                    # Compare DIRECT assignment support for this taxon vs threshold.
                    a_direct = direct_support_map_lift.get(a.assigned_taxid, 0.0)
                    if a_direct >= lift_threshold:
                        corrected.append(a)
                        continue
                    # Low support — find LCA of this taxon with the dominant taxon.
                    taxon_ancs_lift = self.taxonomy.ancestors(a.assigned_taxid, include_self=True)
                    # taxon_ancs_lift is root-first, leaf-last; walk leaf→root.
                    lca_lift = next(
                        (t for t in reversed(taxon_ancs_lift) if t in dom_lift_ancs), None
                    )
                    if lca_lift is None:
                        # No shared ancestor — set unclassified.
                        corrected.append(dc_replace(
                            a,
                            assigned_taxid="0",
                            assigned_name="unclassified",
                            assigned_rank="no rank",
                            status=a.status + ";em_coherence_lift:unclassified",
                        ))
                        continue
                    # Walk up from lca_lift to a named Linnean rank if needed.
                    lca_lift_rank = self.taxonomy.rank(lca_lift)
                    if _RANK_LEVEL.get(lca_lift_rank, -1) < 0:
                        cur = self.taxonomy.parent(lca_lift)
                        for _ in range(10):
                            if not cur or cur == lca_lift:
                                break
                            lvl = _RANK_LEVEL.get(self.taxonomy.rank(cur), -1)
                            if lvl >= 0:
                                lca_lift = cur
                                lca_lift_rank = self.taxonomy.rank(cur)
                                break
                            lca_lift = cur
                            cur = self.taxonomy.parent(cur)
                    # If the LCA is the same as the assigned taxon (assignment IS a shared
                    # broad ancestor), keep as-is to avoid lifting already-broad calls.
                    if lca_lift == a.assigned_taxid:
                        corrected.append(a)
                        continue
                    # If LCA is at or above the max-level threshold, set unclassified
                    # rather than polluting the report with very broad ancestor calls.
                    max_lca = self._em_coherence_max_lca_level
                    lca_level_val = _RANK_LEVEL.get(lca_lift_rank, -1)
                    if max_lca > 0 and lca_level_val >= max_lca:
                        corrected.append(dc_replace(
                            a,
                            assigned_taxid="0",
                            assigned_name="unclassified",
                            assigned_rank="no rank",
                            status=a.status + f";em_coherence_lift:unclassified(lca_too_broad:{self.taxonomy.name(lca_lift)})",
                        ))
                    else:
                        corrected.append(dc_replace(
                            a,
                            assigned_taxid=lca_lift,
                            assigned_name=self.taxonomy.name(lca_lift),
                            assigned_rank=lca_lift_rank,
                            status=a.status + f";em_coherence_lift:{a.assigned_name}->{self.taxonomy.name(lca_lift)}",
                        ))
                assignments = corrected

        return self._apply_rank_cap(assignments)

    def _apply_rank_cap(self, assignments: List[ReadAssignment]) -> List[ReadAssignment]:
        """Demote assignments above max_assignment_rank_level to unclassified.

        Clade / no-rank nodes are resolved by walking up to the nearest named
        Linnean ancestor before comparing against the cap level.
        """
        max_rank = self._max_assignment_rank_level
        if max_rank <= 0:
            return assignments
        capped: List[ReadAssignment] = []
        for a in assignments:
            if a.assigned_taxid == "0":
                capped.append(a)
                continue
            rank = self.taxonomy.rank(a.assigned_taxid)
            level = _RANK_LEVEL.get(rank, -1)
            if level < 0:
                cur = self.taxonomy.parent(a.assigned_taxid)
                for _ in range(25):
                    if not cur or cur == a.assigned_taxid:
                        level = 99
                        break
                    level = _RANK_LEVEL.get(self.taxonomy.rank(cur), -1)
                    if level >= 0:
                        break
                    cur = self.taxonomy.parent(cur)
            if level > max_rank:
                capped.append(dc_replace(
                    a,
                    assigned_taxid="0",
                    assigned_name="unclassified",
                    assigned_rank="no rank",
                    status=a.status + f";rank_cap:unclassified({a.assigned_name}@{rank})",
                ))
            else:
                capped.append(a)
        return capped


def _merge_intervals(intervals: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    if not intervals:
        return []
    intervals = sorted((min(a, b), max(a, b)) for a, b in intervals)
    merged = [intervals[0]]
    for a, b in intervals[1:]:
        la, lb = merged[-1]
        if a <= lb + 1:
            merged[-1] = (la, max(lb, b))
        else:
            merged.append((a, b))
    return merged


def _max_overlap(intervals: List[Tuple[int, int]]) -> int:
    events: List[Tuple[int, int]] = []
    for a, b in intervals:
        events.append((a, 1))
        events.append((b + 1, -1))
    depth = max_depth = 0
    for _pos, delta in sorted(events):
        depth += delta
        max_depth = max(max_depth, depth)
    return max_depth


def _gini(values: List[float]) -> float:
    """Gini coefficient (0=uniform, ~1=all concentrated).

    For coverage depth profiles: 0 means reads are evenly spread,
    ~1 means nearly all depth is at one position (stacking FP pattern).
    """
    n = len(values)
    if n == 0:
        return 0.0
    total = sum(values)
    if total == 0.0:
        return 0.0
    sorted_v = sorted(values)
    # Standard sorted-values formula (0-indexed i):
    # G = (2 * sum((i+1)*v_i)) / (n * total) - (n+1)/n
    weighted = sum((i + 1) * v for i, v in enumerate(sorted_v))
    return (2.0 * weighted) / (n * total) - (n + 1.0) / n


def _gini_sparse(nonzero_vals: List[float], n_total: int) -> float:
    """Gini from a sparse representation: nonzero_vals are the occupied-bin counts,
    n_total is the total number of bins (including zeros).

    Mathematically equivalent to _gini() on a zero-padded list but avoids
    materializing O(n_total) entries when n_total >> len(nonzero_vals).
    Zeros sort first and contribute 0 to the weighted sum, so only the
    non-zero values (sorted ascending) need to be enumerated, offset by n_zero.
    """
    if n_total == 0:
        return 0.0
    total = sum(nonzero_vals)
    if total == 0.0:
        return 0.0
    n_nonzero = len(nonzero_vals)
    n_zero = n_total - n_nonzero
    sorted_nz = sorted(nonzero_vals)
    # In the full sorted array, zeros occupy positions 0..(n_zero-1) (contributing 0),
    # then non-zero values occupy positions n_zero..(n_total-1).
    weighted = sum((n_zero + j + 1) * v for j, v in enumerate(sorted_nz))
    return (2.0 * weighted) / (n_total * total) - (n_total + 1.0) / n_total


def _composite_authenticity(row: "TaxonSummary", damage_weight: float = 1.0) -> float:
    """Composite authenticity score (0-1, higher = more likely a genuine TP).

    Phase 1 weights — will be calibrated against garg_up_bac benchmark in Phase 2.
    Combines: posterior (LCA fit), coverage_gini (uniformity), stack_excess_score
    (reference-length-normalised stacking), breadth, damage, read count,
    mean_windows_per_ref (distinct genomic regions), and max_per_ref_ses.

    ``damage_weight`` scales the damage component (0.0–1.0). For short-reference
    kingdoms (fungi ITS, protist 18S) the damage signal is unreliable with <10 reads
    on a ~1 kb reference; set damage_weight≈0.3 and the freed weight shifts to breadth.
    """
    import math

    # Anti-stacking: SES=1 → uniform (good), SES>>1 → stacked (bad).
    ses_ok = clamp(1.0 / max(row.stack_excess_score, 1.0))

    # Gini: 0=uniform (good), 1=stacked (bad).
    gini_ok = 1.0 - clamp(row.coverage_gini)

    # Breadth (0-1 already).
    breadth_ok = clamp(row.best_reference_breadth)

    # Posterior (0-1).
    post_ok = clamp(row.mean_posterior)

    # Damage — aDNA authenticity marker.
    dmg_ok = clamp(row.mean_damage_score)

    # Statistical support: saturates at ~10 000 reads.
    reads_ok = clamp(math.log10(max(row.evidence_reads, 1)) / 4.0)

    # Distinct genomic regions per reference. Saturates at 100 windows/ref (log scale).
    # Stacking FPs land ~1 window/ref; genuine TPs cover 10-1000s of windows.
    windows_ok = clamp(math.log10(max(row.mean_windows_per_ref, 1.0)) / 2.0)

    # Per-reference SES: catches stacking diluted in global SES by many clean accessions.
    per_ref_ses_ok = clamp(1.0 / max(row.max_per_ref_ses, 1.0))

    # Damage weight scaling: freed weight (1 − damage_weight) × 0.07 shifts to breadth.
    _dmg_w = clamp(damage_weight)
    _breadth_w = 0.10 + 0.07 * (1.0 - _dmg_w)

    # Weighted sum (weights sum to 1.0 — Phase 1 placeholder; calibrate in Phase 2).
    return clamp(
        0.25 * post_ok
        + 0.20 * windows_ok
        + 0.15 * gini_ok
        + 0.15 * ses_ok
        + _breadth_w * breadth_ok
        + 0.07 * _dmg_w * dmg_ok
        + 0.05 * per_ref_ses_ok
        + 0.03 * reads_ok
    )


def _interval_stats(intervals_by_ref: Dict[str, List[Tuple[int, int, Optional[int]]]]) -> Tuple[int, float, float, int, int, float, int, float, int, int, int, float]:
    """Return n_refs, best_breadth, mean_breadth, max_stack, best_span, top_locus_fraction,
    evidence_reads, coverage_gini, max_50bp_bin, total_ref_length,
    n_covered_windows, max_per_ref_ses."""
    if not intervals_by_ref:
        return 0, 0.0, 0.0, 0, 0, 0.0, 0, 0.0, 0, 0, 0, 0.0
    breadths: List[float] = []
    max_stack = 0
    top_locus_fraction = 0.0
    total_intervals = sum(len(v) for v in intervals_by_ref.values())
    # Track span for the reference with the most reads (most likely to be diagnostic).
    best_span = 0
    best_span_ref_count = 0
    # Phase 1 coverage-quality metrics — sparse representation for Gini.
    # Instead of materializing a zero-padded list (O(n_bins) per ref, up to 100K/ref),
    # accumulate only the non-zero bin counts and total bin count. _gini_sparse()
    # computes the identical Gini coefficient in O(k log k) where k = occupied bins.
    gini_nonzero: List[float] = []  # non-zero 50bp bin counts across all refs
    gini_n_total: int = 0           # total bin count (including uncovered bins)
    max_50bp_bin: int = 0
    total_ref_length: int = 0
    # Track which ref accessions have already contributed their slen to total_ref_length.
    _seen_ref_slens: Dict[str, int] = {}
    # Reference-relative stacking metrics.
    n_covered_windows: int = 0       # distinct (accession, 50bp-bin) pairs with ≥1 read
    per_ref_ses_vals: List[float] = []  # per-accession SES; max taken after loop
    for _ref, vals in intervals_by_ref.items():
        intervals = [(a, b) for a, b, _slen in vals]
        slens = [s for _a, _b, s in vals if s]
        slen = max(slens) if slens else None
        merged = _merge_intervals(intervals)
        covered = sum(b - a + 1 for a, b in merged)
        if slen:
            breadths.append(clamp(covered / max(1, slen)))
            # Accumulate reference length (once per accession).
            if _ref not in _seen_ref_slens:
                _seen_ref_slens[_ref] = slen
                total_ref_length += slen
        # If slen is unavailable, do not fake genome/reference breadth. The intervals
        # still contribute to stacking metrics, but breadth remains 0 and is flagged.
        max_stack = max(max_stack, _max_overlap(intervals))
        # Stacking proxy: many reads starting in the same 10 bp bin on the same reference.
        bins: Dict[int, int] = {}
        for a, _b in intervals:
            bins[a // 10] = bins.get(a // 10, 0) + 1
        if vals:
            top_locus_fraction = max(top_locus_fraction, max(bins.values()) / len(vals))
        # Span on this reference (requires no slen).
        if len(vals) >= best_span_ref_count:
            span = max(b for a, b in intervals) - min(a for a, b in intervals)
            if len(vals) > best_span_ref_count or span > best_span:
                best_span = span
                best_span_ref_count = len(vals)
        # 50bp-bin counts for Gini and SES computation.
        bins_50: Dict[int, int] = {}
        for a, _b in intervals:
            bins_50[a // 50] = bins_50.get(a // 50, 0) + 1
        ref_max_50 = max(bins_50.values()) if bins_50 else 0
        max_50bp_bin = max(max_50bp_bin, ref_max_50)
        if slen:
            n_bins = max(1, (slen + 49) // 50)
        else:
            # No slen — treat occupied bins as the full set (no zero-padding possible).
            n_bins = len(bins_50)
        gini_nonzero.extend(float(c) for c in bins_50.values())
        gini_n_total += n_bins
        # Reference-relative stacking metrics.
        n_covered_windows += len(bins_50)
        ref_n_reads = len(vals)
        if slen and ref_n_reads > 0:
            per_ref_ses_vals.append((ref_max_50 * slen) / (ref_n_reads * 50))
    best_b = max(breadths) if breadths else 0.0
    mean_b = sum(breadths) / len(breadths) if breadths else 0.0
    coverage_gini = _gini_sparse(gini_nonzero, gini_n_total)
    max_per_ref_ses = max(per_ref_ses_vals) if per_ref_ses_vals else 0.0
    return len(intervals_by_ref), best_b, mean_b, max_stack, best_span, top_locus_fraction, total_intervals, coverage_gini, max_50bp_bin, total_ref_length, n_covered_windows, max_per_ref_ses


def _stack_excluded_stats(
    intervals_by_ref: Dict[str, List[Tuple[int, int, Optional[int]]]],
    depth_multiplier: float = 5.0,
    min_abs_stack_depth: int = 3,
) -> Tuple[int, float, int]:
    """Return (non_stack_reads, best_stack_excluded_breadth, best_ref_length).

    For each reference accession that has a known length (slen):
    1. Build 50bp bin depths.
    2. Identify stacked bins: depth > depth_multiplier × mean_depth AND depth >= min_abs_stack_depth.
       The absolute floor prevents sparse coverage (1 read per bin on a long reference) from
       being classified as stacked relative to a near-zero mean.
    3. Count reads with ≥1 interval bin falling outside the stacked zone.
    4. Compute reference breadth after zeroing stacked positions.

    Returns metrics for the reference with the highest stack-excluded breadth.
    Returns (0, 0.0, 0) when no reference length is available.
    """
    if not intervals_by_ref:
        return 0, 0.0, 0

    total_non_stack: int = 0
    best_excl_breadth: float = 0.0
    best_ref_len: int = 0

    for _ref, vals in intervals_by_ref.items():
        intervals = [(a, b) for a, b, _slen in vals]
        slens = [s for _a, _b, s in vals if s]
        slen = max(slens) if slens else None
        if not intervals or not slen:
            continue

        # 50bp bin depth map for this reference.
        bins_50: Dict[int, int] = {}
        for a, b in intervals:
            for bk in range(a // 50, b // 50 + 1):
                bins_50[bk] = bins_50.get(bk, 0) + 1

        # Mean depth zero-padded to reference length.
        n_bins = max(1, (slen + 49) // 50)
        mean_depth = sum(bins_50.values()) / n_bins
        rel_threshold = depth_multiplier * max(mean_depth, 1e-9)
        stack_set: Set[int] = {
            bk for bk, cnt in bins_50.items()
            if cnt > rel_threshold and cnt >= min_abs_stack_depth
        }

        # Reads with ≥1 interval bin outside the stack zone.
        ref_non_stack = 0
        for a, b in intervals:
            if set(range(a // 50, b // 50 + 1)) - stack_set:
                ref_non_stack += 1
        total_non_stack += ref_non_stack

        # Stack-excluded breadth: walk each interval, keep only segments in non-stacked bins.
        clean: List[Tuple[int, int]] = []
        for a, b in intervals:
            cur = a
            while cur <= b:
                bk = cur // 50
                seg_end = min(b, (bk + 1) * 50 - 1)
                if bk not in stack_set:
                    clean.append((cur, seg_end))
                cur = seg_end + 1
        covered = sum(e - s + 1 for s, e in _merge_intervals(clean)) if clean else 0
        excl_breadth = clamp(covered / slen)
        if excl_breadth > best_excl_breadth or (
            excl_breadth == best_excl_breadth and slen > best_ref_len
        ):
            best_excl_breadth = excl_breadth
            best_ref_len = slen

    return total_non_stack, best_excl_breadth, best_ref_len


def summarize_assignments(
    taxonomy: Taxonomy,
    assignments: List[ReadAssignment],
    sample_roles: Optional[Dict[str, str]] = None,
    config: Optional[Dict] = None,
    contaminants: Optional["set[str]"] = None,
    regional_taxa: Optional[Dict[str, Dict[str, str]]] = None,
    palynology_taxa: Optional["set[str]"] = None,
    fossil_taxa: Optional["set[str]"] = None,
    sample_depths: Optional[Dict[str, int]] = None,
    fossil_taxa_by_sample: Optional[Dict[str, "set[str]"]] = None,
    fos_evidence_by_sample: Optional[Dict[str, "list[str]"]] = None,
) -> List[TaxonSummary]:
    sample_roles = sample_roles or {}
    config = config or {}
    contaminants = contaminants or set()
    regional_taxa = regional_taxa or {}
    palynology_taxa = palynology_taxa or set()
    fossil_taxa = fossil_taxa or set()
    sample_depths = sample_depths or {}
    fossil_taxa_by_sample = fossil_taxa_by_sample or {}
    fos_evidence_by_sample = fos_evidence_by_sample or {}
    # Tally total reads per sample (all assignments, assigned + unassigned) for normalization.
    _total_reads_per_sample: Dict[str, int] = {}
    _assigned_wt_per_sample: Dict[str, float] = {}
    for a in assignments:
        _total_reads_per_sample[a.sample_id] = _total_reads_per_sample.get(a.sample_id, 0) + 1
        if a.assigned_taxid != "0":
            _assigned_wt_per_sample[a.sample_id] = _assigned_wt_per_sample.get(a.sample_id, 0.0) + a.posterior
    # direct_* track counts for reads assigned exactly to that taxid (not descendants).
    direct_hard: Dict[Tuple[str, str], int] = {}
    direct_wt: Dict[Tuple[str, str], float] = {}
    # cum_* are propagated up to all ancestors in one pass — avoids descendants() traversal.
    cum_hard: Dict[Tuple[str, str], int] = {}
    cum_wt: Dict[Tuple[str, str], float] = {}
    # Running sum accumulators instead of lists — same mean/max result, O(1) memory per taxon.
    cum_post_sum: Dict[Tuple[str, str], float] = {}
    cum_conf: Dict[Tuple[str, str], int] = {}
    cum_dmg_sum: Dict[Tuple[str, str], float] = {}
    cum_dmg_cnt: Dict[Tuple[str, str], int] = {}
    cum_dmg_max: Dict[Tuple[str, str], float] = {}
    # Posterior-weighted sum of per-read mean_hit_uniqueness, for computing taxon-level mean.
    cum_uq_sum: Dict[Tuple[str, str], float] = {}
    coverage: Dict[Tuple[str, str], Dict[str, List[Tuple[int, int, Optional[int]]]]] = {}

    for a in assignments:
        if a.assigned_taxid == "0":
            continue
        dkey = (a.sample_id, a.assigned_taxid)
        direct_hard[dkey] = direct_hard.get(dkey, 0) + 1
        direct_wt[dkey] = direct_wt.get(dkey, 0.0) + a.posterior
        has_coords = a.best_subject_id and a.best_sstart is not None and a.best_send is not None
        if has_coords:
            interval: Tuple[int, int, Optional[int]] = (
                min(int(a.best_sstart), int(a.best_send)),
                max(int(a.best_sstart), int(a.best_send)),
                a.best_slen,
            )
        is_conflicted = "conflicted" in a.status
        # Include ALL reads where qseq was present (damage_computable=True), even when
        # damage_score=0.0. Excluding zero-damage reads creates selection bias: reads
        # where no C happened to fall at the 5' terminus get excluded, leaving only
        # those that had a C there — which inflates the mean on non-damaged samples.
        dmg_available = a.damage_computable
        # Coverage intervals tracked at the direct-assignment level only (not propagated to
        # ancestors), keeping memory O(reads) instead of O(reads × ancestor_depth).
        # Ancestor nodes derive coverage stats solely from directly-assigned reads; their
        # propagated read counts still appear in cum_hard/cum_wt as normal.
        if has_coords:
            coverage.setdefault(dkey, {}).setdefault(a.best_subject_id, []).append(interval)
        # Propagate to all ancestors (including self) in one ancestors() call.
        for anc in taxonomy.ancestors(a.assigned_taxid, include_self=True):
            akey = (a.sample_id, anc)
            cum_hard[akey] = cum_hard.get(akey, 0) + 1
            cum_wt[akey] = cum_wt.get(akey, 0.0) + a.posterior
            cum_post_sum[akey] = cum_post_sum.get(akey, 0.0) + a.posterior
            if is_conflicted:
                cum_conf[akey] = cum_conf.get(akey, 0) + 1
            if dmg_available:
                cum_dmg_sum[akey] = cum_dmg_sum.get(akey, 0.0) + a.damage_score
                cum_dmg_cnt[akey] = cum_dmg_cnt.get(akey, 0) + 1
                d = a.damage_score
                if d > cum_dmg_max.get(akey, 0.0):
                    cum_dmg_max[akey] = d
            cum_uq_sum[akey] = cum_uq_sum.get(akey, 0.0) + a.posterior * a.mean_hit_uniqueness

    rows: Dict[Tuple[str, str], TaxonSummary] = {}
    for (sample, taxid) in cum_hard:
        key = (sample, taxid)
        n_hard = cum_hard[key]
        post_sum = cum_post_sum.get(key, 0.0)
        dmg_cnt = cum_dmg_cnt.get(key, 0)
        n_refs, best_breadth, mean_breadth, max_stack, best_span, top_locus_frac, ev_reads, \
            cov_gini, max_50bp_bin, total_ref_len, n_cov_win, max_pref_ses = \
            _interval_stats(coverage.get(key, {}))
        mean_dmg = cum_dmg_sum.get(key, 0.0) / dmg_cnt if dmg_cnt > 0 else 0.0
        max_dmg = cum_dmg_max.get(key, 0.0)
        _W = 50  # bin width for SES calculation
        _ses = (max_50bp_bin * total_ref_len) / (ev_reads * _W) if ev_reads > 0 and total_ref_len > 0 else 0.0
        rows[key] = TaxonSummary(
            taxid=taxid,
            name=taxonomy.name(taxid),
            rank=taxonomy.rank(taxid),
            sample_id=sample,
            direct_hard_reads=direct_hard.get(key, 0),
            direct_weighted_reads=direct_wt.get(key, 0.0),
            cumulative_hard_reads=cum_hard[key],
            cumulative_weighted_reads=cum_wt.get(key, 0.0),
            mean_posterior=post_sum / n_hard if n_hard > 0 else 0.0,
            conflicted_reads=cum_conf.get(key, 0),
            evidence_reads=ev_reads,
            unique_references=n_refs,
            best_reference_breadth=best_breadth,
            mean_reference_breadth=mean_breadth,
            best_reference_span=best_span,
            max_stack_depth=max_stack,
            top_locus_fraction=top_locus_frac,
            stack_concentration=max_stack / ev_reads if ev_reads > 0 else 0.0,
            coverage_gini=cov_gini,
            stack_excess_score=_ses,
            n_covered_windows=n_cov_win,
            mean_windows_per_ref=n_cov_win / max(1, n_refs),
            max_per_ref_ses=max_pref_ses,
            effective_ref_length=total_ref_len,
            mean_damage_score=mean_dmg,
            max_damage_score=max_dmg,
            mean_hit_uniqueness=cum_uq_sum.get(key, post_sum) / post_sum if post_sum > 0 else 1.0,
        )

    summary_cfg = config.get("summary", {}) if config else {}
    # Pre-compute kingdom for each unique taxid in rows.  Used for:
    # (1) composite_authenticity damage weight scaling per kingdom, and
    # (2) stacking/SES filter gating (filters only safe for bacteria/archaea).
    _row_kingdoms: Dict[str, str] = {}
    for _, _tid in rows:
        if _tid not in _row_kingdoms:
            _row_kingdoms[_tid] = _taxon_kingdom(_tid, taxonomy)

    # Compute composite_authenticity for all rows (requires all other fields populated).
    # For short-reference kingdoms (fungi/ITS, protist/18S) the damage signal is unreliable
    # with few reads on a ~1 kb reference — down-weight damage and up-weight breadth instead.
    _composite_dmg_w_fungi = float(summary_cfg.get("composite_damage_weight_fungi", 1.0))
    for (_, _tid), row in rows.items():
        _kdom = _row_kingdoms.get(_tid, "unknown")
        _dmg_w = _composite_dmg_w_fungi if _kdom in _SHORT_REF_KINGDOMS else 1.0
        row.composite_authenticity = _composite_authenticity(row, damage_weight=_dmg_w)

    negative_samples = {s for s, role in sample_roles.items() if str(role).lower() in {"negative", "neg", "blank", "control_negative", "extraction_blank", "library_blank"}}
    neg_by_tax: Dict[str, float] = {}
    for (sample, taxid), row in rows.items():
        if sample in negative_samples:
            neg_by_tax[taxid] = neg_by_tax.get(taxid, 0.0) + row.cumulative_weighted_reads
    min_reads_flag = int(summary_cfg.get("min_cumulative_reads_for_confident_call", 3))
    min_breadth_flag = float(summary_cfg.get("min_best_reference_breadth_for_confident_call", 0.0))
    max_top_locus = float(summary_cfg.get("max_top_locus_fraction", 0.80))
    # Hard filter: remove taxa where stack_concentration >= threshold (conserved-locus FPs).
    # stack_concentration = max_stack_depth / evidence_reads. FP taxa from conserved mtDNA
    # positions have all reads at one site (stk ≈ 0.2–1.0), while genuine taxa with reads
    # spread across a genome have stk ≈ 0.004–0.02. Default 1.01 = disabled.
    max_stk_hard = float(summary_cfg.get("max_stack_concentration_filter", 1.01))
    if max_stk_hard <= 1.0:
        stacked_fp_keys = [
            k for k, row in rows.items()
            if row.stack_concentration >= max_stk_hard
            and row.cumulative_hard_reads >= min_reads_flag
            and _row_kingdoms.get(row.taxid, "unknown") in _STACKING_FILTER_KINGDOMS
        ]
        for k in stacked_fp_keys:
            del rows[k]
    # Hard filter: reference-relative stacking (default 0.0 = disabled for both).
    # min_mean_windows_per_ref: remove taxa whose reads don't cover multiple distinct
    #   50bp genomic windows per accession. Stacking FPs sit in 1-2 windows; genuine
    #   TPs with even a partial genome spread across many windows.
    # max_per_ref_ses_filter: per-accession SES cap. Global SES can be diluted when
    #   one stacked accession is outnumbered by many clean ones; per-ref SES catches it.
    min_win_per_ref = float(summary_cfg.get("min_mean_windows_per_ref", 0.0))
    max_pref_ses_hard = float(summary_cfg.get("max_per_ref_ses_filter", 0.0))
    if min_win_per_ref > 0.0 or max_pref_ses_hard > 0.0:
        ref_stacking_fp_keys = [
            k for k, row in rows.items()
            if row.cumulative_hard_reads >= min_reads_flag
            and row.evidence_reads > 0  # skip ancestor-only rows (dkey-only coverage gives ev=0)
            and _row_kingdoms.get(row.taxid, "unknown") in _STACKING_FILTER_KINGDOMS
            and (
                (min_win_per_ref > 0.0 and row.mean_windows_per_ref < min_win_per_ref)
                or (max_pref_ses_hard > 0.0 and row.max_per_ref_ses > max_pref_ses_hard)
            )
        ]
        for k in ref_stacking_fp_keys:
            del rows[k]
    # Global SES hard filter (standalone): remove taxa with stack_excess_score above threshold.
    # Use --ses-filter alone only when your TP organisms are well-represented across the whole
    # genome in NT; some genuine taxa (e.g. organisms with unusual NT coverage) can have high SES.
    # Suggested value: 2000 for most bacterial scenarios (TP max observed: ~1244, S. maltophilia).
    max_ses_hard = float(summary_cfg.get("max_stack_excess_score_filter", 0.0))
    if max_ses_hard > 0.0:
        ses_fp_keys = [
            k for k, row in rows.items()
            if row.cumulative_hard_reads >= min_reads_flag
            and row.stack_excess_score > max_ses_hard
            and _row_kingdoms.get(row.taxid, "unknown") in _STACKING_FILTER_KINGDOMS
        ]
        for k in ses_fp_keys:
            del rows[k]
    # Combined SES + low-windows filter (safer): remove taxa where reads are extremely stacked
    # (high SES) AND cover very few distinct genomic windows (low mean_windows_per_ref).
    # This combination is the fingerprint of conserved-locus cross-clade FPs (Laurasiatheria,
    # Actinopterygii, plant superclades).
    #
    # WARNING: UNSAFE when simulated/expected taxa include large-genome eukaryotes (animals,
    # plants). TP-ancestors of animals/plants also show SES > 40K AND windows < 15 because
    # reads from many species are assigned to high-level LCA nodes, each reference getting only
    # a few windows covered. E.g. Eukaryota (SES=1.2M, wins=14.4), Metazoa (SES=813K, wins=7.6),
    # Boreoeutheria (SES=917K, wins=6.8), Canis genus (SES=116K, wins=12.8) — all would be
    # wrongly removed. Only apply this filter for bacteria-only datasets where the windows/ref
    # cutoff reliably separates large-genome cross-clade FPs from genuine bacterial TP-ancestors.
    # Calibration for bac-only: --ses-combined-filter 40000 --ses-combined-max-windows 15.
    ses_combined_threshold = float(summary_cfg.get("ses_combined_filter", 0.0))
    ses_combined_max_windows = float(summary_cfg.get("ses_combined_max_windows_per_ref", 5.0))
    if ses_combined_threshold > 0.0:
        combined_ses_fp_keys = [
            k for k, row in rows.items()
            if row.cumulative_hard_reads >= min_reads_flag
            and row.stack_excess_score > ses_combined_threshold
            and row.mean_windows_per_ref < ses_combined_max_windows
            and _row_kingdoms.get(row.taxid, "unknown") in _STACKING_FILTER_KINGDOMS
        ]
        for k in combined_ses_fp_keys:
            del rows[k]
    # Stack-exclusion robustness check: for each taxon showing meaningful stacking
    # (stack_concentration > 0.05), remove reads from stacked 50bp bins and test
    # whether the assignment still holds.  Kingdom-agnostic — safe for eukaryotes
    # because short references (organelle/barcode) are flagged rather than removed.
    # Disabled by default; enable with stack_exclusion_enabled=True in [summary].
    stack_excl_enabled = bool(summary_cfg.get("stack_exclusion_enabled", False))
    if stack_excl_enabled:
        stack_excl_min_ref = int(summary_cfg.get("stack_exclusion_min_ref_length", 1000))
        stack_excl_k = float(summary_cfg.get("stack_exclusion_depth_multiplier", 5.0))
        stack_excl_min_abs = int(summary_cfg.get("stack_exclusion_min_abs_stack_depth", 3))
        stack_excl_min_reads = int(summary_cfg.get("stack_exclusion_min_non_stack_reads", min_reads_flag))
        stack_excl_min_breadth = float(summary_cfg.get("stack_exclusion_min_non_stack_breadth", 0.02))
        stack_fp_keys = []
        for k, row in rows.items():
            if row.cumulative_hard_reads < min_reads_flag or row.evidence_reads == 0:
                continue
            cov = coverage.get(k, {})
            if not cov:
                continue
            ns_reads, ns_breadth, ref_len = _stack_excluded_stats(
                cov, depth_multiplier=stack_excl_k, min_abs_stack_depth=stack_excl_min_abs
            )
            row.non_stack_reads = ns_reads
            row.stack_excluded_breadth = ns_breadth
            row.best_ref_length = ref_len
            if ref_len < stack_excl_min_ref:
                # Reference too short to reliably detect within-locus stacking.
                row.stacking_status = "short_reference_unassessable"
                continue
            # Only engage when meaningful stacking is present: stk > 10% AND at least
            # min_abs_stack_depth reads concentrated at a single position.
            approx_max_stack = round(row.stack_concentration * max(row.evidence_reads, 1))
            if row.stack_concentration <= 0.10 or approx_max_stack < stack_excl_min_abs:
                continue
            # Robustness test: does sufficient non-stack evidence remain?
            if ns_reads >= stack_excl_min_reads and ns_breadth >= stack_excl_min_breadth:
                row.stacking_status = "stacking_present_but_supported"
            else:
                row.stacking_status = "stacking_fp"
                stack_fp_keys.append(k)
        for k in stack_fp_keys:
            del rows[k]
    # Per-kingdom minimum-reads filter for short-reference kingdoms (fungi ITS, protists).
    # Stacking/SES filters are disabled for these kingdoms, so a low-read hit has no
    # secondary quality gate.  min_reads_fungi defaults to the global min_reads when not set,
    # keeping the filter transparent unless explicitly enabled.
    min_reads_fungi = int(summary_cfg.get("min_reads_fungi", min_reads_flag))
    if min_reads_fungi > min_reads_flag:
        fungi_low_reads_keys = [
            k for k, row in rows.items()
            if _row_kingdoms.get(row.taxid, "unknown") in _SHORT_REF_KINGDOMS
            and row.cumulative_hard_reads < min_reads_fungi
        ]
        for k in fungi_low_reads_keys:
            del rows[k]
    # Composite-authenticity hard filter: remove taxa below the minimum composite score.
    # composite_authenticity is a weighted combination of breadth, posterior, damage, and SES
    # (range 0–1). Low values indicate reads concentrated on conserved loci shared across many
    # taxa (e.g. plant cp gene families), which is the hallmark of cross-family FPs.
    # Only applied to taxa with cumulative_hard_reads >= min_reads_flag to protect low-read
    # but genuine detections (e.g. trace ancient DNA from rare species). Default 0.0 = disabled.
    min_comp_hard = float(summary_cfg.get("min_composite_authenticity", 0.0))
    if min_comp_hard > 0.0:
        comp_low_keys = [
            k for k, row in rows.items()
            if row.composite_authenticity < min_comp_hard
            and row.cumulative_hard_reads >= min_reads_flag
        ]
        for k in comp_low_keys:
            del rows[k]
    # Name-based node suppression: remove NCBI taxonomy "bucket" nodes that are always noise.
    # Three categories:
    #   (1) "X (in: Y)" — NCBI's leftover-grouping notation for sequences not classified to a
    #       named species within clade Y; these are never the target taxon.
    #   (2) "environmental samples" — a meta-taxon for bulk environmental sequences; not informative.
    #   (3) "uncultured X" — uncharacterised isolates; uninformative at species level.
    # These nodes do not appear as ancestors of named species in NCBI taxonomy, so suppressing
    # them cannot create false negatives. Default True when enabled via mode profile.
    if bool(summary_cfg.get("suppress_unclassified_nodes", False)):
        _NCBI_BUCKET_NODES = frozenset({
            "unclassified sequences", "unclassified entries",
            "other sequences", "environmental samples",
            "artificial sequences", "synthetic constructs",
        })
        def _is_noise_node(name: str) -> bool:
            lname = name.lower()
            return (
                " (in: " in name
                or lname.startswith("environmental samples")
                or lname.startswith("uncultured ")
                or lname.startswith("unclassified sequences")
                or lname.startswith("unclassified entries")
                or lname in _NCBI_BUCKET_NODES
            )
        noise_keys = [k for k, row in rows.items() if _is_noise_node(row.name)]
        for k in noise_keys:
            del rows[k]
    phylo_factor = int(summary_cfg.get("incongruent_with_dominant_dominant_factor", 10))
    phylo_min = phylo_factor * min_reads_flag
    # A family must have >= this fraction of sample direct reads to be "dominant".
    # Prevents low-count off-target families from self-qualifying as dominant.
    incongruent_dominant_frac = float(summary_cfg.get("incongruent_dominant_fraction", 0.01))

    # Build per-sample set of "dominant" family-level ancestors (taxa with many reads).
    # Used to flag low-count taxa that are phylogenetically remote from the main signal.
    _FAM_RANKS = ["family", "superfamily"]
    _ORD_RANKS = ["order", "suborder", "infraorder", "superorder"]

    def _fam_or_ord(taxid: str) -> Optional[str]:
        fam = taxonomy.ranked_ancestor(taxid, _FAM_RANKS)
        return fam if fam else taxonomy.ranked_ancestor(taxid, _ORD_RANKS)

    # Per-sample total direct reads (for fraction-based dominant threshold).
    sample_direct_totals: Dict[str, int] = {}
    for (sample, _), row in rows.items():
        sample_direct_totals[sample] = sample_direct_totals.get(sample, 0) + row.direct_hard_reads

    dominant_clades: Dict[str, Set[str]] = {}  # sample_id → set of family/order taxids
    for (sample, taxid), row in rows.items():
        sample_total = sample_direct_totals.get(sample, 0)
        dom_threshold = max(float(phylo_min), incongruent_dominant_frac * sample_total)
        if (sample not in negative_samples
                and row.direct_hard_reads > 0
                and row.cumulative_hard_reads >= dom_threshold):
            clade = _fam_or_ord(taxid)
            if clade:
                dominant_clades.setdefault(sample, set()).add(clade)

    for row in rows.values():
        neg = neg_by_tax.get(row.taxid, 0.0)
        row.negative_weighted_reads = neg
        denom = row.cumulative_weighted_reads + neg
        row.blank_fraction = neg / denom if denom > 0 else 0.0
        flags = []
        if neg > 0 and row.sample_id not in negative_samples:
            flags.append("negative_control_overlap")
        if row.conflicted_reads > 0:
            flags.append("conflicted_reads_present")
        if row.cumulative_hard_reads < min_reads_flag:
            flags.append("low_read_support")
        if row.evidence_reads > 0 and row.unique_references > 0 and row.best_reference_breadth == 0:
            flags.append("no_reference_length_for_breadth")
        if min_breadth_flag and row.best_reference_breadth < min_breadth_flag and row.cumulative_hard_reads >= min_reads_flag:
            flags.append("low_reference_breadth")
        # stacked_locus_warning: fires when most reads pile on a single locus.
        # Uses stack_concentration (peak depth / total reads) rather than top_locus_fraction
        # because top_locus_fraction is per-reference and mis-fires on multi-reference taxa.
        if row.cumulative_hard_reads >= min_reads_flag and row.stack_concentration >= max_top_locus:
            flags.append("stacked_locus_warning")
        # high_ses_warning: fires when stack_excess_score is high relative to a typical TP.
        # SES > 2000 is the fingerprint of conserved-locus cross-clade FPs; TPs rarely exceed
        # ~1250 even for unevenly-covered organisms. Also fires for mean_windows_per_ref < 3
        # (fewer than 3 distinct genomic regions covered — consistent with conserved-locus hits).
        if row.cumulative_hard_reads >= min_reads_flag and (
            row.stack_excess_score > 2000.0
            or (row.mean_windows_per_ref > 0.0 and row.mean_windows_per_ref < 3.0)
        ):
            flags.append("high_ses_warning")
        # Flag taxa whose family is absent from the sample's dominant clades.
        if (row.sample_id not in negative_samples
                and row.direct_hard_reads > 0):
            dom = dominant_clades.get(row.sample_id, set())
            if dom:
                clade = _fam_or_ord(row.taxid)
                if clade and clade not in dom:
                    flags.append("incongruent_with_dominant")
        if contaminants and (row.taxid in contaminants or row.name in contaminants):
            flags.append("user_contaminant")
        if row.stacking_status == "stacking_present_but_supported":
            flags.append("stacking_present_but_supported")
        elif row.stacking_status == "short_reference_unassessable":
            flags.append("short_reference_unassessable")
        row.flags = ";".join(flags)

    # ── Normalization (Item 17) ──────────────────────────────────────────────
    for row in rows.values():
        # Use user-supplied depth if provided, else total reads seen in assignments.
        depth = sample_depths.get(row.sample_id) or _total_reads_per_sample.get(row.sample_id, 0)
        if depth > 0:
            row.reads_per_million = row.cumulative_hard_reads / depth * 1_000_000
            row.weighted_per_million = row.cumulative_weighted_reads / depth * 1_000_000
        total_wt = _assigned_wt_per_sample.get(row.sample_id, 0.0)
        if total_wt > 0:
            row.relative_abundance = row.cumulative_weighted_reads / total_wt

    # ── Multi-proxy authenticity tier (Item 3) ───────────────────────────────
    _auth_cfg = config.get("authenticity", {})
    _dmg_tier1 = float(_auth_cfg.get("damage_tier1", 0.05))
    _dmg_tier2 = float(_auth_cfg.get("damage_tier2", 0.02))
    _min_breadth = float(_auth_cfg.get("min_breadth", 0.05))
    _min_posterior = float(_auth_cfg.get("min_posterior", 0.75))
    _max_blank = float(_auth_cfg.get("max_blank_fraction", 0.10))
    _badges = {1: "★", 2: "●", 3: "◆", 4: "▲", 5: "○", 0: "✕"}
    _pch    = {1: 8,   2: 16,  3: 18,  4: 17,  5: 1,  0: 4}

    for row in rows.values():
        if row.sample_id in negative_samples:
            continue
        # Check eco/pal/fos support by taxid, name, or ancestors
        anc_and_self = set(taxonomy.ancestors(row.taxid, include_self=True)) | {row.name}
        row.eco_support = regional_taxa and bool(anc_and_self & set(regional_taxa.keys()))
        row.pal_support = palynology_taxa and bool(anc_and_self & palynology_taxa)
        # Per-sample (site/age-resolved) fossil support takes priority over the batch-wide set
        _sample_fos = fossil_taxa_by_sample.get(row.sample_id) if fossil_taxa_by_sample else None
        if _sample_fos is not None:
            row.fos_support = bool(anc_and_self & _sample_fos)
            if row.fos_support:
                row.fos_evidence_text = "; ".join(fos_evidence_by_sample.get(row.sample_id, []))
        else:
            row.fos_support = fossil_taxa and bool(anc_and_self & fossil_taxa)
        # Count independent positive evidence lines
        ev = [
            row.mean_damage_score >= _dmg_tier2,
            row.best_reference_breadth >= _min_breadth,
            row.mean_posterior >= _min_posterior,
            # Blank fraction only counts as positive evidence when negatives exist in the batch.
            bool(negative_samples) and row.blank_fraction <= _max_blank and row.sample_id not in negative_samples,
            bool(row.eco_support),
            bool(row.pal_support),
            bool(row.fos_support),
        ]
        n = sum(ev)
        row.n_support_lines = n
        flagged = "user_contaminant" in row.flags or row.blank_fraction > 0.5
        if flagged:
            t = 0
        elif n >= 4 and row.mean_damage_score >= _dmg_tier1:
            t = 1
        elif n >= 3 and row.mean_damage_score >= _dmg_tier2:
            t = 2
        elif n >= 2:
            t = 3
        elif n == 1:
            t = 4
        else:
            t = 5
        row.authenticity_tier      = t
        row.authenticity_badge     = _badges[t]
        row.authenticity_tier_pch  = _pch[t]

    # ── Confidence tier ──────────────────────────────────────────────────────────
    # Summarises evidence quality in a 3-level label visible to end users.
    # HIGH: well-supported detection with genome-wide spread.
    # MEDIUM: moderate support, passes basic quality checks.
    # LOW: marginal support; present in output but treat with caution.
    for row in rows.values():
        reads = row.cumulative_hard_reads
        stk = row.stack_concentration
        mwin = row.mean_windows_per_ref
        warn_flags = {"stacked_locus_warning", "high_ses_warning", "low_read_support"}
        has_warning = bool(warn_flags & set(row.flags.split(";")) if row.flags else set())
        no_aln = row.evidence_reads == 0  # ancestor-only rows with no direct alignments
        if (reads >= 10
                and stk < 0.05
                and (no_aln or mwin >= 5)
                and not has_warning):
            row.confidence_tier = "HIGH"
        elif (reads >= 3
              and stk < 0.15
              and (no_aln or mwin >= 2)):
            row.confidence_tier = "MEDIUM"
        else:
            row.confidence_tier = "LOW"

    return sorted(rows.values(), key=lambda r: (r.sample_id, -r.cumulative_weighted_reads, r.name))


def assignments_to_rows(assignments: List[ReadAssignment]) -> List[Dict[str, object]]:
    return [asdict(a) for a in assignments]


def summaries_to_rows(summaries: List[TaxonSummary]) -> List[Dict[str, object]]:
    return [asdict(s) for s in summaries]


