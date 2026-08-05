"""Generalized, source-agnostic confidence classification for MetaMerge.

``classify.classify_status()`` was written around exactly two evidence sources
(MEGAN counts + Holi/metaDMG damage) as a single hardcoded cascade. This module
adds a source-agnostic layer on top of it, supporting any non-empty subset of
{MEGAN, Holi, Fillet} evidence for a given taxon, without touching
``classify_status()`` itself:

  - When both MEGAN and Holi evidence are present for a taxon, ``classify_status_v2``
    calls the original ``classify_status()`` directly and unmodified -- so the
    2-source case is byte-for-byte identical to pre-Fillet-integration behavior
    by construction, not merely by testing. This is the load-bearing backward-
    compatibility guarantee described in FILLET_METAMERGE_INTEGRATION_PLAN.md.
  - When Fillet evidence is also present alongside a fully-corroborated MEGAN+Holi
    call, a new upgraded status distinguishes 3-source corroboration from a
    pairwise-only one.
  - When one or more of the three sources is entirely absent for a taxon (e.g. a
    real Holi+Fillet-only run with no MEGAN at all), a dedicated classification
    path handles that specific combination -- there are 7 non-empty subsets of a
    3-element set; each has its own handler below.
  - A new continuous ``ensemble_support_score`` (0-1) is always computed
    alongside the discrete status, reflecting how many independent sources
    corroborate the call and how strongly -- mirroring Fillet's own established
    design philosophy of never collapsing continuous and discrete evidence into
    one number, applied here at the ensemble level. A genuine cross-classifier
    taxonomic disagreement (Fillet's own call anchors to a different taxon than
    MEGAN/Holi's, beyond the existing lineage-support rank-cap) is not a
    separate dead-end status -- it subtracts a penalty from this score while the
    full underlying per-source evidence remains visible in the output record.

Evidence-source signal vocabulary
-----------------------------------
Each source contributes a ``SourceSignals`` bundle. All signal fields default to
``None``, meaning "this source did not supply this signal" -- distinct from
``False`` ("this source supplied this signal and it evaluated negative"). This
distinction matters: a source that was never run for a taxon must never
silently be treated as "checked this and it failed."
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

from .classify import classify_status


@dataclass
class SourceSignals:
    """Evidence signals contributed by ONE source (MEGAN, Holi, or Fillet)
    for one taxon anchor. ``present=False`` means this source has no data
    at all for this taxon (as opposed to data that evaluated negatively)."""

    present: bool = False
    exact_damage_support: Optional[bool] = None
    exact_damage_support_ge100: Optional[bool] = None
    strong_count_support: Optional[bool] = None
    some_count_support: Optional[bool] = None
    max_count: float = 0.0
    qc_label: Optional[str] = None
    blank_associated: Optional[bool] = None
    # Fillet-specific signals (None for MEGAN/Holi sources).
    authenticated: Optional[bool] = None
    eco_support: Optional[bool] = None
    pal_support: Optional[bool] = None
    fos_support: Optional[bool] = None


@dataclass
class TaxonEvidence:
    """All evidence for one taxon anchor across up to 3 sources, plus the
    taxonomy-reconciliation outcome between them."""

    megan: SourceSignals = field(default_factory=SourceSignals)
    holi: SourceSignals = field(default_factory=SourceSignals)
    fillet: SourceSignals = field(default_factory=SourceSignals)
    lineage_support: bool = False
    max_real_count: float = 0.0
    # True if Fillet's own taxonomic call anchors to a different taxon than
    # MEGAN/Holi's, beyond the existing conservative lineage-support rank-cap
    # (holi.is_meaningful_low_rank_lineage_support's own thresholds) -- computed
    # by the caller (build_merge/the taxonomy-reconciliation layer), not by this
    # module, since it requires the full lineage-matching machinery in holi.py.
    discordant: bool = False


def _classify_megan_holi(ev: TaxonEvidence, thresholds: dict) -> tuple[str, str]:
    """The original, untouched 2-source path -- calls classify_status() directly."""
    return classify_status(
        exact_damage_support=bool(ev.holi.exact_damage_support),
        exact_damage_support_ge100=bool(ev.holi.exact_damage_support_ge100),
        strong_count_support=bool(ev.megan.strong_count_support),
        some_count_support=bool(ev.megan.some_count_support),
        lineage_support=ev.lineage_support,
        blank_associated=bool(ev.megan.blank_associated),
        qc_label=ev.holi.qc_label or "not-applicable",
        max_real_count=ev.max_real_count,
        thresholds=thresholds,
    )


def _classify_holi_fillet(ev: TaxonEvidence, thresholds: dict) -> tuple[str, str]:
    """Holi+Fillet, no MEGAN. Tyler's expected primary real-world 2-source case
    going forward. There is no MEGAN count matrix, so "count support" comes
    from Fillet's own direct_hard_reads instead of MEGAN's per-library counts;
    Fillet's own authentication (composite_authenticity + authenticity_tier)
    substitutes for MEGAN's "strong_count_support" role in the original
    cascade's structure, since it is itself already a read-count-aware signal."""
    fillet_count_support = bool(ev.fillet.strong_count_support)
    if ev.holi.exact_damage_support_ge100 and (fillet_count_support or ev.fillet.authenticated):
        if ev.fillet.authenticated and (ev.fillet.eco_support or ev.fillet.pal_support or ev.fillet.fos_support):
            return (
                "Very high confidence",
                "Holi exact damage-supported (>=100 reads); Fillet-authenticated with independent eco/pal/fos support",
            )
        return (
            "Very high confidence",
            "Holi exact damage-supported (>=100 reads); Fillet count/authenticity support",
        )
    if ev.holi.exact_damage_support and (fillet_count_support or ev.fillet.authenticated):
        return "High confidence", "Holi exact damage-supported; Fillet count/authenticity support"
    if ev.holi.exact_damage_support or fillet_count_support or ev.fillet.authenticated or ev.lineage_support:
        basis = []
        if ev.holi.exact_damage_support:
            basis.append("Holi damage-supported (exact)")
        if ev.fillet.authenticated:
            basis.append("Fillet-authenticated")
        if fillet_count_support:
            basis.append("Fillet strong counts")
        if ev.lineage_support:
            basis.append("lineage-supported")
        return "Supported", "; ".join(basis)
    if ev.max_real_count >= thresholds["tentative_min_reads"]:
        return "Tentative", "weak/incomplete DNA support"
    if ev.max_real_count >= thresholds["weak_support_min_reads"]:
        return "Weak support", "very weak DNA support"
    return "Weak support", "no non-control support"


def _classify_megan_fillet(ev: TaxonEvidence, thresholds: dict) -> tuple[str, str]:
    """MEGAN+Fillet, no Holi -- Fillet's own damage assessment substitutes for
    Holi's exact-damage-support role."""
    fillet_damage_ok = bool(ev.fillet.exact_damage_support)
    if fillet_damage_ok and ev.megan.strong_count_support:
        basis = "Fillet damage-supported; strong MEGAN counts"
        if ev.fillet.authenticated:
            basis += "; Fillet-authenticated"
        return "High confidence", basis
    if fillet_damage_ok or ev.megan.strong_count_support or ev.fillet.authenticated or (ev.lineage_support and ev.megan.some_count_support):
        basis = []
        if fillet_damage_ok:
            basis.append("Fillet damage-supported")
        if ev.megan.strong_count_support:
            basis.append("strong MEGAN counts")
        elif ev.megan.some_count_support:
            basis.append("some MEGAN counts")
        if ev.fillet.authenticated:
            basis.append("Fillet-authenticated")
        if ev.lineage_support:
            basis.append("lineage-supported")
        return "Supported", "; ".join(basis)
    if ev.max_real_count >= thresholds["tentative_min_reads"]:
        return "Tentative", "weak/incomplete DNA support"
    if ev.max_real_count >= thresholds["weak_support_min_reads"]:
        return "Weak support", "very weak DNA support"
    return "Weak support", "no non-control support"


def _classify_single_source(source: SourceSignals, source_name: str, max_real_count: float, thresholds: dict) -> tuple[str, str]:
    """Only one of the three sources has any data for this taxon at all. No
    cross-classifier corroboration is possible, so this path is deliberately
    conservative -- the ceiling is "Supported", never "High"/"Very high
    confidence", since those tiers are meant to reflect convergent independent
    evidence, which a single source can never provide by definition."""
    if source_name == "fillet":
        if source.authenticated and (source.eco_support or source.pal_support or source.fos_support):
            return "Supported", "Fillet-authenticated with independent eco/pal/fos support (single-source: Fillet only)"
        if source.authenticated:
            return "Supported", "Fillet-authenticated (single-source: Fillet only)"
    elif source_name == "holi" and source.exact_damage_support:
        return "Supported", "Holi exact damage-supported (single-source: Holi only)"
    elif source_name == "megan" and source.strong_count_support:
        return "Supported", "strong MEGAN counts (single-source: MEGAN only)"

    if max_real_count >= thresholds["tentative_min_reads"]:
        return "Tentative", f"weak/incomplete DNA support (single-source: {source_name} only)"
    if max_real_count >= thresholds["weak_support_min_reads"]:
        return "Weak support", f"very weak DNA support (single-source: {source_name} only)"
    return "Weak support", f"no non-control support (single-source: {source_name} only)"


_DEFAULT_ENSEMBLE_WEIGHTS = {
    "megan_weight": 0.30,
    "holi_weight": 0.35,
    "fillet_weight": 0.35,
    "discordance_penalty": 0.40,
    "fillet_support_line_bonus": 0.05,
}


def classify_status_v2(
    ev: TaxonEvidence, config: dict,
) -> tuple[str, str, float]:
    """Source-agnostic classification: works for any non-empty subset of
    {MEGAN, Holi, Fillet} evidence present for a taxon.

    Args:
        ev: The taxon's evidence bundle.
        config: Full MetaMerge config dict (uses config["thresholds"] and
            config.get("ensemble_score", {})).

    Returns (status, basis_summary, ensemble_support_score). status uses the
    same vocabulary as classify_status() (Very high confidence / High
    confidence / Supported / Tentative / Weak support / Blank-associated) plus
    one new value ("Very high confidence (3-source corroborated)") reachable
    only when all three sources are present and mutually agree.
    """
    thresholds = config["thresholds"]

    if ev.megan.blank_associated or ev.holi.blank_associated:
        return "Blank-associated", "blank-associated", 0.0

    have_megan, have_holi, have_fillet = ev.megan.present, ev.holi.present, ev.fillet.present

    if have_megan and have_holi:
        status, basis = _classify_megan_holi(ev, thresholds)
    elif have_holi and have_fillet:
        status, basis = _classify_holi_fillet(ev, thresholds)
    elif have_megan and have_fillet:
        status, basis = _classify_megan_fillet(ev, thresholds)
    elif have_megan:
        status, basis = _classify_single_source(ev.megan, "megan", ev.max_real_count, thresholds)
    elif have_holi:
        status, basis = _classify_single_source(ev.holi, "holi", ev.max_real_count, thresholds)
    elif have_fillet:
        status, basis = _classify_single_source(ev.fillet, "fillet", ev.max_real_count, thresholds)
    else:
        return "Weak support", "no evidence source data for this taxon", 0.0

    # 3-source corroboration upgrade: only reachable when all three sources
    # independently support the taxon AND their taxonomic calls agree.
    if (
        have_megan and have_holi and have_fillet
        and status in ("Very high confidence", "High confidence")
        and ev.fillet.authenticated
        and (ev.fillet.eco_support or ev.fillet.pal_support or ev.fillet.fos_support)
        and not ev.discordant
    ):
        status = "Very high confidence (3-source corroborated)"
        basis += "; Fillet-authenticated with independent eco/pal/fos support, all 3 sources agree"

    if ev.discordant:
        basis += "; NOTE: Fillet's taxonomic call disagrees with MEGAN/Holi beyond the lineage rank-cap"

    score = compute_ensemble_support_score(ev, config.get("ensemble_score", {}))
    return status, basis, score


def compute_ensemble_support_score(ev: TaxonEvidence, ensemble_weights: dict | None = None) -> float:
    """Continuous 0-1 multi-method-agreement score, computed alongside (not
    instead of) the discrete status above. Each present source contributes its
    configured weight only when ITS OWN evidence for this taxon clears that
    source's own support bar; a genuine cross-classifier taxonomic
    disagreement subtracts a penalty rather than merely omitting a source's
    weight, since discordance is positive evidence of a problem, not simply an
    absence of corroboration. Fillet's independent eco/pal/fos support lines
    each add a small bonus on top, since they are evidence dimensions
    MEGAN/Holi have no equivalent of at all.

    Args:
        ev: The taxon's evidence bundle.
        ensemble_weights: config["ensemble_score"] sub-dict; falls back to
            _DEFAULT_ENSEMBLE_WEIGHTS for any key not present (or if None).
    """
    weights = {**_DEFAULT_ENSEMBLE_WEIGHTS, **(ensemble_weights or {})}
    megan_w  = weights["megan_weight"]
    holi_w   = weights["holi_weight"]
    fillet_w = weights["fillet_weight"]
    discordance_penalty = weights["discordance_penalty"]
    support_line_bonus  = weights["fillet_support_line_bonus"]

    score = 0.0
    if ev.megan.present and ev.megan.strong_count_support:
        score += megan_w
    elif ev.megan.present and ev.megan.some_count_support:
        score += megan_w * 0.5

    if ev.holi.present and ev.holi.exact_damage_support_ge100:
        score += holi_w
    elif ev.holi.present and ev.holi.exact_damage_support:
        score += holi_w * 0.6

    if ev.fillet.present and ev.fillet.authenticated:
        score += fillet_w
        for flag in (ev.fillet.eco_support, ev.fillet.pal_support, ev.fillet.fos_support):
            if flag:
                score += support_line_bonus

    if ev.discordant:
        score -= discordance_penalty

    return max(0.0, min(1.0, score))
