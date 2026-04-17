"""Tests for MetaMerge lineage-support logic (is_meaningful_low_rank_lineage_support).

Verifies that:
  - Direct ancestor/descendant relationships at allowed ranks return True.
  - Broad/forbidden clade names are rejected.
  - Same-name taxa are not treated as self-supporting.
  - Rank-step limits are enforced.
"""

from metamerge.holi import is_meaningful_low_rank_lineage_support


_CFG_BASE = {
    "lineage": {
        "forbid_broad_names": ["Mammalia", "Bilateria", "Boreoeutheria"],
        "rank_levels": ["species", "genus", "family", "order", "class"],
    },
    "thresholds": {
        "lineage_max_steps": 2,
        "lineage_max_rank_level": "family",
    },
}


def test_family_to_genus_accepted():
    """A genus within a family should provide valid lineage support."""
    assert is_meaningful_low_rank_lineage_support(
        focal_tax_name="Bovidae",
        focal_tax_rank="family",
        focal_path="Animalia;Chordata;Mammalia;Artiodactyla;Bovidae",
        candidate_tax_name="Bison",
        candidate_tax_rank="genus",
        candidate_path="Animalia;Chordata;Mammalia;Artiodactyla;Bovidae;Bison",
        cfg=_CFG_BASE,
    )


def test_broad_lineage_support_rejected():
    """Broad clade names (Boreoeutheria) in the forbid list must not count."""
    assert not is_meaningful_low_rank_lineage_support(
        focal_tax_name="Puma",
        focal_tax_rank="genus",
        focal_path="Animalia;Chordata;Mammalia;Carnivora;Felidae;Puma",
        candidate_tax_name="Boreoeutheria",
        candidate_tax_rank="class",
        candidate_path="Animalia;Chordata;Mammalia;Boreoeutheria",
        cfg=_CFG_BASE,
    )


def test_same_name_rejected():
    """A taxon cannot provide lineage support for itself."""
    assert not is_meaningful_low_rank_lineage_support(
        focal_tax_name="Bison",
        focal_tax_rank="genus",
        focal_path="Animalia;Chordata;Mammalia;Artiodactyla;Bovidae;Bison",
        candidate_tax_name="Bison",
        candidate_tax_rank="genus",
        candidate_path="Animalia;Chordata;Mammalia;Artiodactyla;Bovidae;Bison",
        cfg=_CFG_BASE,
    )


def test_too_many_rank_steps_rejected():
    """A lineage relationship spanning more than lineage_max_steps is rejected.

    With rank_levels = [species, genus, family, order, class] and max_steps=2:
      species (index 0) → order (index 3) = 3 steps, which exceeds the limit.
    Note: species → family is exactly 2 steps (allowed); species → order is 3 (rejected).
    """
    assert not is_meaningful_low_rank_lineage_support(
        focal_tax_name="Bison bison",
        focal_tax_rank="species",
        focal_path="Animalia;Chordata;Mammalia;Artiodactyla;Bovidae;Bison;Bison bison",
        candidate_tax_name="Artiodactyla",
        candidate_tax_rank="order",
        candidate_path="Animalia;Chordata;Mammalia;Artiodactyla",
        cfg=_CFG_BASE,
    )


def test_no_lineage_relationship_rejected():
    """Taxa that share no ancestor/descendant relationship are rejected."""
    assert not is_meaningful_low_rank_lineage_support(
        focal_tax_name="Puma",
        focal_tax_rank="genus",
        focal_path="Animalia;Chordata;Mammalia;Carnivora;Felidae;Puma",
        candidate_tax_name="Bison",
        candidate_tax_rank="genus",
        candidate_path="Animalia;Chordata;Mammalia;Artiodactyla;Bovidae;Bison",
        cfg=_CFG_BASE,
    )
