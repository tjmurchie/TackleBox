"""Marker / region classification and species-kingdom lookup.

Regions are classified primarily using user-supplied regex rules from a TSV
(``regions_config.tsv``).  If no rule matches, built-in heuristics based on
common NCBI header keywords are used as a fallback.
"""
from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from typing import Dict, List, Tuple

from .utils import warn


@dataclass
class RegionRule:
    region_id: str
    klass: str
    enabled: bool
    regex: str
    ncbi_title_clause: str = ""


def load_regions(path: str) -> List[RegionRule]:
    rules: List[RegionRule] = []
    if not path:
        return rules
    try:
        with open(path, encoding="utf-8", errors="replace") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                rid = (row.get("region_id") or "").strip()
                klass = (row.get("class") or row.get("klass") or "Other").strip()
                en = (row.get("enabled_default") or row.get("enabled") or "1").strip()
                regex = (row.get("regex") or "").strip()
                title = (row.get("ncbi_title_clause") or "").strip()
                if rid:
                    rules.append(
                        RegionRule(
                            rid,
                            klass,
                            en.lower() not in ("0", "false", "no"),
                            regex,
                            title,
                        )
                    )
    except FileNotFoundError:
        warn(f"Regions config not found: {path}")
    return rules


def classify(header: str, rules: List[RegionRule], cfg: dict) -> Tuple[str, str]:
    """Return (marker_class, region_id) for *header*.

    1. Try each enabled rule from *rules* in order.
    2. Fall back to built-in heuristics.
    """
    for rule in rules:
        if not rule.enabled or not rule.regex:
            continue
        try:
            if re.search(rule.regex, header, re.I):
                return rule.klass or "Other", rule.region_id
        except re.error:
            warn(f"Bad regex in regions config for {rule.region_id!r}: {rule.regex!r}")

    h = header.lower()
    if "mitochondr" in h or "[location=mitochondrion]" in h:
        return "Mito", "MITO_BUILTIN"
    if any(x in h for x in ("chloroplast", "plastid", "plastome")) or "[location=chloroplast]" in h:
        return "Plastid", "PLASTID_BUILTIN"
    if re.search(
        r"\b18s\b|\b28s\b|small subunit ribosomal|ssu rrna|ssu rna"
        r"|large subunit ribosomal|lsu rrna|lsu rna|\bits1\b|\bits2\b"
        r"|internal transcribed spacer|5\.8s|histone h3|\bh3\b",
        h,
    ):
        return "NucMark", "NUCMARK_BUILTIN"

    default_class = cfg.get("classification", {}).get("default_class_if_no_match", "Other")
    return default_class, "OTHER_BUILTIN"


def guess_species(header: str) -> Tuple[str, str]:
    """Return (species, genus) guessed from an NCBI-style FASTA header.

    Recognises two common NCBI formats:
    - ``>ACCESSION Genus species ... [organism=...]``
    - ``>ACCESSION ... [Genus species ...]`` (taxonomy bracket at end)
    """
    h = header[1:] if header.startswith(">") else header
    fields = h.split()
    # Format: >ACC Genus species ...
    if (
        len(fields) >= 3
        and re.match(r"^[A-Z][\w.-]+$", fields[1])
        and re.match(r"^[a-z][\w.-]+", fields[2].strip(",;"))
    ):
        return f"{fields[1]} {fields[2].strip(',;')}", fields[1]
    # Format with taxonomy bracket at end: [Genus species]
    m = re.search(r"\[([A-Z][\w.-]+\s+[a-z][\w.-]+)[^\]]*\]\s*$", h)
    if m:
        sp = m.group(1)
        return sp, sp.split()[0]
    return "", ""


def norm_kingdom(k: str) -> str:
    """Normalise a kingdom/phylum string to one of the canonical values."""
    s = (k or "").lower()
    if "animal" in s or s == "metazoa":
        return "Animal"
    if "plant" in s or "plantae" in s or "viridiplantae" in s:
        return "Plant"
    if "fung" in s:
        return "Fungi"
    if "bacteria" in s:
        return "Bacteria"
    if "archaea" in s:
        return "Archaea"
    if "protozoa" in s:
        return "Protozoa"
    if "protist" in s or "chromista" in s:
        return "Protist"
    return "Unknown"


def load_species_kingdom(path: str) -> Tuple[Dict[str, str], Dict[str, str]]:
    """Return (species_to_kingdom, genus_to_kingdom) dicts (keys lower-cased).

    Accepts TSV files with or without a header row.  When a header is present
    the columns ``species``, ``genus``, and ``kingdom`` are located by name;
    otherwise columns are assumed to be in order species, [genus,] kingdom.
    """
    sp2k: Dict[str, str] = {}
    g2k: Dict[str, str] = {}
    if not path:
        return sp2k, g2k
    try:
        with open(path, encoding="utf-8", errors="replace") as f:
            sample = f.readline()
            f.seek(0)
            if sample.lower().startswith("species"):
                reader = csv.DictReader(f, delimiter="\t")
                fields = {x.lower(): x for x in (reader.fieldnames or [])}
                for row in reader:
                    sp = (row.get(fields.get("species", ""), "") or "").strip()
                    genus = (row.get(fields.get("genus", ""), "") or "").strip()
                    kd = norm_kingdom(
                        (row.get(fields.get("kingdom", ""), "") or "").strip()
                    )
                    if sp:
                        sp2k[sp.lower()] = kd
                        genus = genus or sp.split()[0]
                    if genus:
                        g2k.setdefault(genus.lower(), kd)
            else:
                for line in f:
                    if not line.strip():
                        continue
                    parts = line.rstrip("\n\r").split("\t")
                    if len(parts) >= 3:
                        sp, genus, kd = parts[0].strip(), parts[1].strip(), norm_kingdom(parts[2])
                    elif len(parts) >= 2:
                        sp, kd = parts[0].strip(), norm_kingdom(parts[1])
                        genus = sp.split()[0] if sp else ""
                    else:
                        continue
                    if sp:
                        sp2k[sp.lower()] = kd
                    if genus:
                        g2k.setdefault(genus.lower(), kd)
    except FileNotFoundError:
        warn(f"Species-kingdom file not found: {path}")
    return sp2k, g2k
