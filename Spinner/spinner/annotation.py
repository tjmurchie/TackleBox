"""Annotation dataclass and the main per-record annotation loop.

``Annotation`` holds every field that ends up as a column in decisions.tsv.
``annotate()`` is the core per-record loop that populates Annotation objects
from FastaRecords using the loaded QC and classification helpers.
"""
from __future__ import annotations

import dataclasses
import re
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional

from .adapters import find_adapter, is_terminal, load_adapters
from .keywords import load_keywords
from .regions import classify, guess_species, load_regions, load_species_kingdom, norm_kingdom
from .seq_utils import IUPAC_DNA, entropy, gc_frac, max_hpoly, seq_hash
from .utils import progress

# RefSeq accession prefixes — these sequences are NCBI-curated.
_REFSEQ_PREFIXES = (
    "NC_", "NM_", "NR_", "NZ_", "NG_", "XM_", "XR_", "WP_",
    "XP_", "NP_", "AC_", "NW_", "NT_",
)

# Header terms that indicate a voucher or type specimen — high-value aDNA sources.
_VOUCHER_TERMS = (
    "voucher", "type strain", "holotype", "paratype",
    "type specimen", "neotype", "lectotype", "syntype",
)

# Header terms that indicate a complete organelle / genome assembly.
_COMPLETE_TERMS = (
    "complete genome", "complete mitochondrial", "complete plastid",
    "complete sequence", "complete cds", "complete coding sequence",
)


@dataclass
class Annotation:
    # Identity
    accession: str
    record_key: str
    source_file: str
    header: str
    length: int
    seq_sha256: str

    # Taxonomy
    species_guess: str = ""
    genus_guess: str = ""
    kingdom: str = "Unknown"

    # Marker classification
    marker_class: str = "Other"
    region_id: str = ""

    # Basic QC metrics
    n_count: int = 0
    n_fraction: float = 0.0
    non_iupac_count: int = 0
    non_iupac_fraction: float = 0.0
    gc_fraction: float = 0.0
    shannon_entropy: float = 0.0
    max_homopolymer: int = 0

    # Duplicate flags
    duplicate_accession: bool = False
    duplicate_sequence: bool = False

    # Adapter screen
    adapter_hit: bool = False
    adapter_name: str = ""
    adapter_position: str = ""
    adapter_internal: bool = False
    adapter_terminal: bool = False

    # Vector screen (BLAST-based)
    vector_hit: bool = False
    vector_internal: bool = False
    vector_terminal: bool = False

    # Bad keywords
    bad_keyword_hit: bool = False
    bad_keyword: str = ""
    bad_keyword_action: str = ""

    # Taxonomy BLAST
    taxonomy_status: str = "NOT_CHECKED"
    taxonomy_top_hit: str = ""
    taxonomy_top_name: str = ""
    taxonomy_top_pident: str = ""
    taxonomy_top_length: str = ""
    taxonomy_top_staxids: str = ""

    # Windowed BLAST / chimerism
    windowed_status: str = "NOT_CHECKED"

    # Clustering
    cluster_id: str = ""
    cluster_role: str = ""

    # Capping
    cap_rank: int = 0

    # Decision
    decision_score: int = 100
    decision: str = "KEEP"
    reasons: List[str] = field(default_factory=list)

    def add_reason(self, r: str) -> None:
        if r and r not in self.reasons:
            self.reasons.append(r)

    def as_dict(self) -> dict:
        d = dataclasses.asdict(self)
        d["reasons"] = ";".join(self.reasons)
        return d


def annotate(
    records,
    cfg: dict,
    species_kingdom: Optional[str] = None,
    regions_config: Optional[str] = None,
    adapters_tsv: Optional[str] = None,
    bad_keywords_tsv: Optional[str] = None,
) -> Dict[str, Annotation]:
    """Annotate all records and return a dict keyed by record_key.

    Loads species-kingdom mapping, region rules, adapters, and keywords from
    their respective files, then runs per-record QC checks in a single pass.
    """
    sp2k, g2k = load_species_kingdom(species_kingdom or "")
    regions = load_regions(regions_config or "")
    adapters = load_adapters(adapters_tsv or "")
    keywords = load_keywords(bad_keywords_tsv or "")

    ann: Dict[str, Annotation] = {}
    seen_acc: set = set()
    seen_hash: dict = {}

    show_progress = cfg.get("run", {}).get("progress_bars", True)

    for rec in progress(records, "Annotating", show_progress):
        s = rec.seq_upper
        acc = rec.accession

        # --- taxonomy / kingdom ---
        sp, gen = guess_species(rec.header)
        kd = sp2k.get(sp.lower(), g2k.get(gen.lower(), "Unknown"))

        # --- marker classification ---
        klass, rid = classify(rec.header, regions, cfg)
        if kd == "Unknown" and klass == "Plastid" and cfg.get("classification", {}).get(
            "coerce_unknown_plastid_to_plant", True
        ):
            kd = "Plant"

        h = seq_hash(s)
        a = Annotation(
            accession=acc,
            record_key=rec.id,
            source_file=rec.source_file,
            header=rec.header,
            length=len(s),
            seq_sha256=h,
            species_guess=sp,
            genus_guess=gen,
            kingdom=kd,
            marker_class=klass,
            region_id=rid,
        )

        # --- basic sequence metrics ---
        a.n_count = s.count("N")
        a.n_fraction = a.n_count / len(s) if s else 1.0
        a.non_iupac_count = sum(1 for b in s if b not in IUPAC_DNA)
        a.non_iupac_fraction = a.non_iupac_count / len(s) if s else 1.0
        a.gc_fraction = gc_frac(s)
        a.shannon_entropy = entropy(s)
        a.max_homopolymer = max_hpoly(s)

        # --- duplicate detection ---
        bqc = cfg.get("basic_qc", {})
        if bqc.get("remove_duplicate_accessions", True):
            if acc in seen_acc:
                a.duplicate_accession = True
                a.add_reason("duplicate_accession")
            else:
                seen_acc.add(acc)
        if bqc.get("remove_exact_duplicate_sequences", True):
            if h in seen_hash:
                a.duplicate_sequence = True
                a.add_reason("duplicate_sequence")
            else:
                seen_hash[h] = acc

        # --- basic QC filters ---
        if cfg.get("steps", {}).get("basic_qc", True):
            if a.length < int(bqc.get("min_length", 0)):
                a.add_reason("length_below_min")
            if a.length > int(bqc.get("max_length", 10**12)):
                a.add_reason("length_above_max")
            # Per-class minimum length (overrides global min_length for specific classes).
            class_min = int(bqc.get("min_length_by_class", {}).get(klass, 0))
            if class_min and a.length < class_min:
                a.add_reason("length_below_class_min")
            if a.n_fraction > float(bqc.get("max_n_fraction", 1.0)):
                a.add_reason("n_fraction_high")
            if a.non_iupac_fraction > float(bqc.get("max_non_iupac_fraction", 1.0)):
                a.add_reason("non_iupac_fraction_high")
            if a.shannon_entropy < float(bqc.get("min_shannon_entropy", 0.0)):
                a.add_reason("low_complexity")
            if a.max_homopolymer > int(bqc.get("max_homopolymer_run", 10**9)):
                a.add_reason("homopolymer_long")

        # --- scoring bonuses (applied regardless of basic_qc step toggle) ---
        # RefSeq accessions are NCBI-curated — higher confidence than GenBank.
        if acc.startswith(_REFSEQ_PREFIXES):
            a.add_reason("refseq_preferred")
        # Voucher / type specimens are high-value reference sources.
        header_lower = rec.header.lower()
        if any(t in header_lower for t in _VOUCHER_TERMS):
            a.add_reason("voucher_keyword")
        # Complete organelle / genome assemblies are better mapping references.
        if any(t in header_lower for t in _COMPLETE_TERMS):
            a.add_reason("complete_organelle")

        # --- adapter screen ---
        if cfg.get("steps", {}).get("adapter_screen", True) and adapters:
            win = int(cfg.get("adapter_screen", {}).get("terminal_window_bp", 25))
            for ad in adapters:
                hit = find_adapter(s, ad, cfg)
                if hit is not None:
                    pos, strand = hit
                    a.adapter_hit = True
                    a.adapter_name = ad.name
                    a.adapter_position = f"{pos + 1}-{pos + len(ad.sequence)}({strand})"
                    if is_terminal(pos, len(ad.sequence), len(s), win):
                        a.adapter_terminal = True
                        a.add_reason("adapter_terminal")
                    else:
                        a.adapter_internal = True
                        a.add_reason("adapter_internal")
                    break  # report first hit only

        # --- bad keyword screen ---
        if cfg.get("steps", {}).get("bad_keyword_screen", True) and keywords:
            case_sensitive = cfg.get("bad_keywords", {}).get("case_sensitive", False)
            htxt = rec.header if case_sensitive else rec.header.lower()
            for kw in keywords:
                needle = kw.keyword if case_sensitive else kw.keyword.lower()
                if needle in htxt:
                    a.bad_keyword_hit = True
                    a.bad_keyword = kw.keyword
                    a.bad_keyword_action = kw.action
                    a.add_reason(
                        "bad_keyword_reject" if kw.action == "reject" else "bad_keyword_review"
                    )
                    break  # report first hit only

        ann[rec.id] = a

    return ann


# ---------------------------------------------------------------------------
# Post-annotation rescue helpers
# ---------------------------------------------------------------------------

def _trim_at_adapter(seq: str, adapter_pos_str: str, min_length: int) -> Optional[str]:
    """Split *seq* at an adapter hit, returning the longest valid flanking piece.

    *adapter_pos_str* is the 1-based "start-end(strand)" string stored in
    Annotation.adapter_position.  Returns None when neither flank meets
    *min_length*.
    """
    m = re.match(r"(\d+)-(\d+)\([+-]\)", adapter_pos_str)
    if not m:
        return None
    start_0 = int(m.group(1)) - 1  # 0-based start of adapter
    end_0 = int(m.group(2))        # 0-based exclusive end (= 1-based inclusive end)
    five_prime = seq[:start_0]
    three_prime = seq[end_0:]
    valid = [s for s in (five_prime, three_prime) if len(s) >= min_length]
    return max(valid, key=len) if valid else None


def rescue_cross_species_duplicates(ann: Dict[str, "Annotation"], cfg: dict) -> int:
    """Promote one representative for each species with all-duplicate-sequence records.

    Global deduplication marks every second-and-later occurrence of a sequence
    as duplicate_sequence.  When Species B's sequences are all exact copies of
    Species A's, Species B ends up with zero KEEP candidates even though the
    sequence is legitimate — it is simply attributed to two taxa at that locus.

    This function un-flags the best-quality record for each such species so it
    can proceed through external screens and downstream filtering.

    Selection priority: RefSeq accession > longest sequence > fewest Ns.
    Returns the count of rescued records.
    """
    if not cfg.get("basic_qc", {}).get("rescue_cross_species_duplicates", True):
        return 0

    by_species: Dict[str, List[str]] = defaultdict(list)
    for key, a in ann.items():
        sp = a.species_guess.strip()
        if sp:
            by_species[sp].append(key)

    rescued = 0
    for sp, keys in by_species.items():
        non_dup = [k for k in keys
                   if not ann[k].duplicate_sequence and not ann[k].duplicate_accession]
        if non_dup:
            continue  # species already has real unique records

        candidates = [k for k in keys
                      if ann[k].duplicate_sequence
                      and not ann[k].duplicate_accession
                      and "rescued_duplicate" not in ann[k].reasons]
        if not candidates:
            continue  # nothing to rescue (already rescued or all accession-dups)

        def _rank(key: str) -> tuple:
            a = ann[key]
            is_rs = int(a.accession.startswith(_REFSEQ_PREFIXES))
            return (is_rs, a.length, -a.n_count)

        best = max(candidates, key=_rank)
        a = ann[best]
        a.duplicate_sequence = False
        if "duplicate_sequence" in a.reasons:
            a.reasons.remove("duplicate_sequence")
        a.add_reason("rescued_duplicate")
        rescued += 1

    return rescued


def trim_and_rescreen_adapters(
    records: list,
    ann: Dict[str, "Annotation"],
    cfg: dict,
    adapters_tsv: str,
) -> int:
    """Trim internal adapter hits and re-screen the trimmed sequences.

    For each record flagged ``adapter_internal``, the sequence is split at the
    adapter position (keeping the longer clean flank) and re-checked for
    adapter contamination.  Records whose trimmed sequence is adapter-free and
    meets the minimum length are rescued: their annotation is updated and the
    record's sequence is replaced in-place so that external tools receive the
    clean trimmed version.

    Records where trimming would still leave internal adapter contamination, or
    where both flanks are too short, remain as-is.

    Enabled only when ``adapter_screen.trim_and_rescreen: true`` in config.
    Returns the count of rescued records.
    """
    if not cfg.get("adapter_screen", {}).get("trim_and_rescreen", False):
        return 0

    adapters = load_adapters(adapters_tsv or "")
    if not adapters:
        return 0

    min_length = int(cfg.get("basic_qc", {}).get("min_length", 50))
    terminal_win = int(cfg.get("adapter_screen", {}).get("terminal_window_bp", 25))
    record_by_key = {r.id: r for r in records}
    rescued = 0

    for key, a in ann.items():
        if not a.adapter_internal or "adapter_internal" not in a.reasons:
            continue
        r = record_by_key.get(key)
        if r is None:
            continue

        s = r.seq_upper
        trimmed = _trim_at_adapter(s, a.adapter_position, min_length)
        if trimmed is None:
            continue

        # Re-screen the trimmed sequence for any remaining adapter contamination.
        new_hit = False
        new_name = ""
        new_pos_str = ""
        new_internal = False
        new_terminal = False
        for ad in adapters:
            hit = find_adapter(trimmed, ad, cfg)
            if hit is not None:
                pos, strand = hit
                new_hit = True
                new_name = ad.name
                new_pos_str = f"{pos + 1}-{pos + len(ad.sequence)}({strand})"
                if is_terminal(pos, len(ad.sequence), len(trimmed), terminal_win):
                    new_terminal = True
                else:
                    new_internal = True
                break

        if new_internal:
            continue  # still internally contaminated — cannot rescue

        # Rescue: mutate the record's sequence so downstream tools see the clean version.
        r.seq = trimmed

        # Update annotation flags.
        a.adapter_hit = new_hit
        a.adapter_name = new_name
        a.adapter_position = new_pos_str
        a.adapter_internal = False
        a.adapter_terminal = new_terminal
        if "adapter_internal" in a.reasons:
            a.reasons.remove("adapter_internal")
        if new_terminal and "adapter_terminal" not in a.reasons:
            a.add_reason("adapter_terminal")
        a.add_reason("adapter_trimmed_rescued")

        # Recompute sequence metrics on the trimmed sequence.
        a.length = len(trimmed)
        a.n_count = trimmed.count("N")
        a.n_fraction = a.n_count / len(trimmed) if trimmed else 1.0
        a.non_iupac_count = sum(1 for b in trimmed if b not in IUPAC_DNA)
        a.non_iupac_fraction = a.non_iupac_count / len(trimmed) if trimmed else 1.0
        a.gc_fraction = gc_frac(trimmed)
        a.shannon_entropy = entropy(trimmed)
        a.max_homopolymer = max_hpoly(trimmed)

        # Re-check length constraint after trimming.
        if a.length < int(cfg.get("basic_qc", {}).get("min_length", 0)):
            a.add_reason("length_below_min")

        rescued += 1

    return rescued
