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


def classify_status(
    exact_damage_support: bool,
    exact_damage_support_ge100: bool,
    strong_count_support: bool,
    some_count_support: bool,
    lineage_support: bool,
    blank_associated: bool,
    qc_label: str,
    max_real_count: float,
    thresholds: dict,
) -> tuple[str, str]:
    """Assign a DNA-support category and short basis summary (original,
    untouched 2-source MEGAN+Holi cascade).

    Applies the classification rules in priority order.  Blank-associated always
    wins; Very high confidence requires the strictest combination of signals.

    Args:
        exact_damage_support: True if any real library has exact Holi damage+sig.
        exact_damage_support_ge100: True if any supporting row has N_reads ≥ threshold.
        strong_count_support: True if count evidence meets the strong-count thresholds.
        some_count_support: True if at least one real library is positive.
        lineage_support: True if any real library has lineage-consistent support.
        blank_associated: True if blank evidence dominates.
        qc_label: QC label from compute_qc_label ("clean", "caution",
            "strong caution").
        max_real_count: Maximum MEGAN count across all real libraries.
        thresholds: The "thresholds" sub-dict from the MetaMerge config.

    Returns:
        Tuple of ``(status_string, basis_summary_string)``.
    """
    if blank_associated:
        return "Blank-associated", "blank-associated"

    if exact_damage_support_ge100 and strong_count_support and qc_label != "strong caution":
        return (
            "Very high confidence",
            "exact damage-supported; strong counts; >=100 Holi reads",
        )

    if exact_damage_support and strong_count_support:
        if qc_label == "clean":
            return "High confidence", "exact damage-supported; strong counts"
        return "High confidence", f"exact damage-supported; strong counts; {qc_label}"

    if (
        exact_damage_support
        or strong_count_support
        or (lineage_support and some_count_support)
    ):
        basis = []
        if exact_damage_support:
            basis.append("damage-supported (exact)")
        if lineage_support:
            basis.append("lineage-supported")
        if strong_count_support:
            basis.append("strong counts")
        elif some_count_support:
            basis.append("some counts")
        if qc_label != "clean" and exact_damage_support:
            basis.append(qc_label)
        return "Supported", "; ".join(basis) if basis else "supported"

    if max_real_count >= thresholds["tentative_min_reads"]:
        return "Tentative", "weak/incomplete DNA support"

    if max_real_count >= thresholds["weak_support_min_reads"]:
        return "Weak support", "very weak DNA support"

    return "Weak support", "no non-control support"


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


def _classify_primary_additional(
    primary_authenticated: bool,
    primary_strong: bool,
    primary_extra_lines: bool,
    primary_name: str,
    additional_strong: bool,
    additional_some: bool,
    additional_name: str,
    lineage_support: bool,
    max_real_count: float,
    thresholds: dict,
) -> tuple[str, str]:
    """Generic 2-source cascade for a "primary" source that supplies BOTH
    counts and its own authentication decision, plus an "additional" source
    that can only ever upgrade the resulting tier, never gate it.

    Used for both Holi+Fillet and MEGAN+Fillet (Tyler's design, 2026-08-05):
    by default Fillet is the primary source in both, since its own
    composite_authenticity/authenticity_tier + eco/pal/fos lines already
    combine multiple evidence dimensions on their own (unlike a bare MEGAN
    count or a bare Holi damage/significance pair) -- Fillet alone can reach
    any tier here on its own strength. The role assignment is swappable via
    config["source_roles"], in which case the caller passes the OTHER
    source's own fields into the primary_* parameters instead.

    Args:
        primary_authenticated: Primary source's own pass/fail decision.
        primary_strong: Primary source's own "strong count support" signal.
        primary_extra_lines: An independent corroborating line the primary
            source itself supplies (e.g. Fillet's eco/pal/fos support) --
            False for sources with no such extra line (e.g. Holi, MEGAN).
        primary_name: Display name for basis text (e.g. "Fillet").
        additional_strong: Additional source's strongest support signal.
        additional_some: Additional source's weaker/baseline support signal.
        additional_name: Display name for basis text (e.g. "Holi").
        lineage_support: Whether lineage-consistent support exists.
        max_real_count: Best available real-count signal across sources.
        thresholds: The "thresholds" sub-dict of the MetaMerge config.
    """
    # "Very high confidence" needs the primary source authenticated PLUS at
    # least 2 of its 3 possible corroborating signals (its own strong counts,
    # its own extra eco/pal/fos-style line, the additional source's strongest
    # signal) -- deliberately not a hard requirement on primary_strong
    # specifically, so the additional source's own strong corroboration can
    # substitute for it (e.g. Fillet authenticated + eco support + Holi's
    # >=100-read exact damage support reaches the top tier even without
    # Fillet's own read count independently clearing its strong-count bar).
    strength_signal_count = sum([primary_strong, primary_extra_lines, additional_strong])
    if primary_authenticated and strength_signal_count >= 2:
        basis = [f"{primary_name}-authenticated"]
        if primary_strong:
            basis.append(f"strong {primary_name} counts")
        if primary_extra_lines:
            basis.append("independent eco/pal/fos support")
        if additional_strong:
            basis.append(f"{additional_name} strongly corroborates")
        return "Very high confidence", "; ".join(basis)

    if primary_authenticated and (primary_strong or primary_extra_lines or additional_some):
        basis = [f"{primary_name}-authenticated"]
        if primary_strong:
            basis.append(f"strong {primary_name} counts")
        if primary_extra_lines:
            basis.append("independent eco/pal/fos support")
        if additional_some:
            basis.append(f"{additional_name} corroborates")
        return "High confidence", "; ".join(basis)

    if primary_authenticated or primary_strong or additional_some or lineage_support:
        basis = []
        if primary_authenticated:
            basis.append(f"{primary_name}-authenticated")
        if primary_strong:
            basis.append(f"strong {primary_name} counts")
        if additional_some:
            basis.append(f"{additional_name} support")
        if lineage_support:
            basis.append("lineage-supported")
        return "Supported", "; ".join(basis) if basis else "supported"

    if max_real_count >= thresholds["tentative_min_reads"]:
        return "Tentative", "weak/incomplete DNA support"
    if max_real_count >= thresholds["weak_support_min_reads"]:
        return "Weak support", "very weak DNA support"
    return "Weak support", "no non-control support"


def _classify_holi_fillet(ev: TaxonEvidence, thresholds: dict, primary: str = "fillet") -> tuple[str, str]:
    """Holi+Fillet, no MEGAN. Tyler's expected primary real-world 2-source case
    going forward. Default role assignment (primary="fillet"): Fillet is the
    counts+support source, Holi is additional corroborating support -- override
    with config["source_roles"]["holi_fillet_primary"]="holi" to swap them."""
    fillet_extra_lines = bool(ev.fillet.eco_support or ev.fillet.pal_support or ev.fillet.fos_support)
    if primary == "holi":
        return _classify_primary_additional(
            primary_authenticated=bool(ev.holi.exact_damage_support),
            primary_strong=bool(ev.holi.exact_damage_support_ge100),
            primary_extra_lines=False,
            primary_name="Holi",
            additional_strong=bool(ev.fillet.authenticated and fillet_extra_lines),
            additional_some=bool(ev.fillet.authenticated or ev.fillet.strong_count_support),
            additional_name="Fillet",
            lineage_support=ev.lineage_support,
            max_real_count=ev.max_real_count,
            thresholds=thresholds,
        )
    return _classify_primary_additional(
        primary_authenticated=bool(ev.fillet.authenticated),
        primary_strong=bool(ev.fillet.strong_count_support),
        primary_extra_lines=fillet_extra_lines,
        primary_name="Fillet",
        additional_strong=bool(ev.holi.exact_damage_support_ge100),
        additional_some=bool(ev.holi.exact_damage_support),
        additional_name="Holi",
        lineage_support=ev.lineage_support,
        max_real_count=ev.max_real_count,
        thresholds=thresholds,
    )


def _classify_megan_fillet(ev: TaxonEvidence, thresholds: dict, primary: str = "fillet") -> tuple[str, str]:
    """MEGAN+Fillet, no Holi. Default role assignment (primary="fillet"):
    Fillet is the counts+support source, MEGAN's own counts are additional
    corroborating support -- override with
    config["source_roles"]["megan_fillet_primary"]="megan" to swap them."""
    fillet_extra_lines = bool(ev.fillet.eco_support or ev.fillet.pal_support or ev.fillet.fos_support)
    if primary == "megan":
        return _classify_primary_additional(
            primary_authenticated=bool(ev.megan.strong_count_support),
            primary_strong=bool(ev.megan.strong_count_support),
            primary_extra_lines=False,
            primary_name="MEGAN",
            additional_strong=bool(ev.fillet.authenticated and fillet_extra_lines),
            additional_some=bool(ev.fillet.authenticated or ev.fillet.exact_damage_support),
            additional_name="Fillet",
            lineage_support=bool(ev.lineage_support and ev.megan.some_count_support),
            max_real_count=ev.max_real_count,
            thresholds=thresholds,
        )
    return _classify_primary_additional(
        primary_authenticated=bool(ev.fillet.authenticated),
        primary_strong=bool(ev.fillet.strong_count_support),
        primary_extra_lines=fillet_extra_lines,
        primary_name="Fillet",
        additional_strong=bool(ev.megan.strong_count_support),
        additional_some=bool(ev.megan.some_count_support),
        additional_name="MEGAN",
        lineage_support=ev.lineage_support,
        max_real_count=ev.max_real_count,
        thresholds=thresholds,
    )


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

# Minimum agreeing/present sources for the "N-source corroborated" upgrade
# tier below. Kept as a module constant rather than a config threshold since
# it's a structural property of the status vocabulary (what "corroborated"
# means), not a per-project tuning knob like the count/damage thresholds.
_MIN_SOURCES_FOR_CORROBORATION_BONUS = 3


def _source_supports(signals: SourceSignals, name: str) -> bool:
    """Whether one source's OWN evidence, on its own terms, counts as
    supporting a taxon -- the equal-weight building block for
    ``compute_methods_agreement``, deliberately independent of
    ``compute_ensemble_support_score``'s trust-weighted formula (Tyler,
    2026-08-05: "two separate metrics", not one score trying to do both)."""
    if not signals.present:
        return False
    if name == "fillet":
        return bool(signals.authenticated)
    if name == "holi":
        return bool(signals.exact_damage_support)
    if name == "megan":
        return bool(signals.strong_count_support)
    return False


def compute_methods_agreement(ev: TaxonEvidence) -> tuple[int, int]:
    """Equal-weight "how many methods agree" tally, computed alongside (not
    instead of) ``compute_ensemble_support_score``. Every present source that
    clears its OWN authentication bar counts exactly 1, regardless of which
    tool it is or how it's weighted in the trust-weighted score -- this is
    the metric that genuinely generalizes to N classifiers and matches
    Tyler's own framing of the ensemble score as fundamentally "how many
    methods agree."

    Returns:
        Tuple of ``(agreement_count, sources_present_count)``.
    """
    sources = [(ev.megan, "megan"), (ev.holi, "holi"), (ev.fillet, "fillet")]
    sources_present = sum(1 for s, _ in sources if s.present)
    agreement_count = sum(1 for s, name in sources if s.present and _source_supports(s, name))
    return agreement_count, sources_present


def classify_status_v2(
    ev: TaxonEvidence, config: dict,
) -> tuple[str, str, float, int, float]:
    """Source-agnostic classification: works for any non-empty subset of
    {MEGAN, Holi, Fillet} evidence present for a taxon.

    Args:
        ev: The taxon's evidence bundle.
        config: Full MetaMerge config dict (uses config["thresholds"],
            config.get("ensemble_score", {}), and
            config.get("source_roles", {})).

    Returns:
        Tuple of ``(status, basis_summary, ensemble_support_score,
        methods_agreement_count, methods_agreement_fraction)``. status uses
        the same vocabulary as classify_status() (Very high confidence / High
        confidence / Supported / Tentative / Weak support / Blank-associated)
        plus one new value ("Very high confidence (N-source corroborated)")
        reachable when at least 3 present sources independently support AND
        agree on the call, for whatever N that turns out to be (today capped
        at 3 since only MEGAN/Holi/Fillet exist; generalizes automatically if
        more sources are added later).
    """
    thresholds   = config["thresholds"]
    source_roles = config.get("source_roles", {})

    if ev.megan.blank_associated or ev.holi.blank_associated:
        return "Blank-associated", "blank-associated", 0.0, 0, 0.0

    have_megan, have_holi, have_fillet = ev.megan.present, ev.holi.present, ev.fillet.present
    agreement_count, sources_present = compute_methods_agreement(ev)
    agreement_fraction = (agreement_count / sources_present) if sources_present else 0.0

    if have_megan and have_holi:
        status, basis = _classify_megan_holi(ev, thresholds)
    elif have_holi and have_fillet:
        status, basis = _classify_holi_fillet(
            ev, thresholds, primary=source_roles.get("holi_fillet_primary", "fillet")
        )
    elif have_megan and have_fillet:
        status, basis = _classify_megan_fillet(
            ev, thresholds, primary=source_roles.get("megan_fillet_primary", "fillet")
        )
    elif have_megan:
        status, basis = _classify_single_source(ev.megan, "megan", ev.max_real_count, thresholds)
    elif have_holi:
        status, basis = _classify_single_source(ev.holi, "holi", ev.max_real_count, thresholds)
    elif have_fillet:
        status, basis = _classify_single_source(ev.fillet, "fillet", ev.max_real_count, thresholds)
    else:
        return "Weak support", "no evidence source data for this taxon", 0.0, 0, 0.0

    # N-source corroboration upgrade: reachable when at least
    # _MIN_SOURCES_FOR_CORROBORATION_BONUS present sources independently
    # support the taxon. Generalizes the original hardcoded "all 3 of MEGAN/
    # Holi/Fillet present and agree" check to the equal-weight methods_
    # agreement tally, so it scales automatically if more sources are added.
    if (
        sources_present >= _MIN_SOURCES_FOR_CORROBORATION_BONUS
        and agreement_count >= _MIN_SOURCES_FOR_CORROBORATION_BONUS
        and status in ("Very high confidence", "High confidence")
        and not ev.discordant
    ):
        status = f"Very high confidence ({agreement_count}-source corroborated)"
        basis += f"; {agreement_count} of {sources_present} independent methods agree"

    if ev.discordant:
        basis += "; NOTE: Fillet's taxonomic call disagrees with MEGAN/Holi beyond the lineage rank-cap"

    score = compute_ensemble_support_score(ev, config.get("ensemble_score", {}))
    return status, basis, score, agreement_count, agreement_fraction


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
