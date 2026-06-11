"""Taxonomy BLAST sanity checking with optional NCBI taxdump lineage support.

Two levels of taxonomy comparison are supported:

1. **String matching** (always available): the top-hit scientific name is
   compared against the expected genus/species from the header or
   species_kingdom TSV.  This is fast but crude.

2. **Taxdump lineage comparison** (when ``taxonomy_blast.taxdump_dir`` is set):
   NCBI ``nodes.dmp`` and ``names.dmp`` are loaded into a TaxdumpDB.  Hit
   taxids are resolved to full lineages and compared at the genus, family, and
   kingdom levels.  This enables proper cross-kingdom rejection.

Expected BLAST outfmt:
  6 qseqid saccver pident length qlen qstart qend evalue bitscore staxids sscinames

Windowed BLAST (chimerism screen):
  For long sequences, overlapping windows are BLASTed.  Windows whose top-hit
  names diverge are flagged as windowed_blast_conflict.
"""
from __future__ import annotations

import os
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

from .annotation import Annotation
from .utils import info, warn


# ---------------------------------------------------------------------------
# Taxdump database
# ---------------------------------------------------------------------------

class TaxdumpDB:
    """In-memory lookup table built from NCBI taxdump nodes.dmp / names.dmp.

    Only the fields needed for kingdom-level and rank-level comparison are
    stored.  Loading ~2 GB taxdump files takes ~30 s and ~2 GB RAM; this is
    optional and only activated when ``taxdump_dir`` is configured.
    """

    def __init__(self, taxdump_dir: str) -> None:
        self.parent: Dict[int, int] = {}
        self.rank: Dict[int, str] = {}
        self.name: Dict[int, str] = {}
        self._load(taxdump_dir)

    def _load(self, taxdump_dir: str) -> None:
        nodes_path = os.path.join(taxdump_dir, "nodes.dmp")
        names_path = os.path.join(taxdump_dir, "names.dmp")

        if not os.path.exists(nodes_path) or not os.path.exists(names_path):
            raise FileNotFoundError(
                f"taxdump files not found in {taxdump_dir!r}. "
                "Expected nodes.dmp and names.dmp from NCBI taxdump."
            )

        info(f"Loading taxdump nodes from {nodes_path}")
        with open(nodes_path, encoding="utf-8", errors="replace") as f:
            for line in f:
                parts = line.split("|")
                if len(parts) < 3:
                    continue
                try:
                    taxid = int(parts[0].strip())
                    parent = int(parts[1].strip())
                    rank = parts[2].strip()
                    self.parent[taxid] = parent
                    self.rank[taxid] = rank
                except ValueError:
                    continue

        info(f"Loading taxdump names from {names_path}")
        with open(names_path, encoding="utf-8", errors="replace") as f:
            for line in f:
                parts = line.split("|")
                if len(parts) < 4:
                    continue
                try:
                    taxid = int(parts[0].strip())
                    name = parts[1].strip()
                    name_class = parts[3].strip()
                    if name_class == "scientific name":
                        self.name[taxid] = name
                except ValueError:
                    continue

        info(f"Taxdump loaded: {len(self.parent):,} nodes, {len(self.name):,} named taxa")

    def get_lineage(self, taxid: int) -> List[Tuple[str, str]]:
        """Return list of (rank, scientific_name) from species up to root."""
        lineage: List[Tuple[str, str]] = []
        seen: set = set()
        current = taxid
        while current not in seen and current in self.parent:
            seen.add(current)
            rank = self.rank.get(current, "no rank")
            name = self.name.get(current, "")
            lineage.append((rank, name))
            if current == 1:  # root
                break
            current = self.parent[current]
        return lineage

    def get_rank_name(self, taxid: int, target_rank: str) -> str:
        """Return the scientific name at *target_rank* in the lineage of *taxid*."""
        for rank, name in self.get_lineage(taxid):
            if rank == target_rank:
                return name
        return ""

    def get_domain(self, taxid: int) -> str:
        """Return the top-level domain for *taxid*.

        Prefers the "domain" rank (NCBI 2024+ taxdump), then falls back to
        "superkingdom" (older taxdumps).  Returns one of: Bacteria, Eukaryota,
        Archaea, Viruses — or "" if not found.

        Note: NCBI reorganised their taxonomy in 2023-2024.  The former
        "superkingdom" rank (Bacteria / Eukaryota / Archaea) was renamed to
        "domain", and a new "kingdom" rank was introduced within Bacteria
        (e.g. Pseudomonadati, Bacillati) and within Eukaryota.  Callers
        should use get_domain() for cross-kingdom comparisons and NOT rely on
        "kingdom" or "superkingdom" alone.
        """
        lineage = self.get_lineage(taxid)
        # First pass: prefer the "domain" rank (2024+ taxdumps)
        for rank, name in lineage:
            if rank == "domain":
                return name
        # Fallback: "superkingdom" (pre-2024 taxdumps)
        for rank, name in lineage:
            if rank == "superkingdom":
                return name
        return ""

    def get_kingdom(self, taxid: int) -> str:
        """Alias for get_domain() — kept for backwards compatibility."""
        return self.get_domain(taxid)

    def is_same_rank(self, taxid1: int, taxid2: int, rank: str) -> bool:
        """Return True if both taxids share the same name at *rank*."""
        n1 = self.get_rank_name(taxid1, rank)
        n2 = self.get_rank_name(taxid2, rank)
        return bool(n1 and n2 and n1.lower() == n2.lower())


# ---------------------------------------------------------------------------
# Taxonomy BLAST parser
# ---------------------------------------------------------------------------

def parse_tax_blast(
    path: str,
    ann: Dict[str, Annotation],
    cfg: dict,
    taxdb: Optional[TaxdumpDB] = None,
) -> None:
    """Read taxonomy BLAST TSV and annotate *ann* in-place.

    Expected columns (outfmt 6):
      qseqid saccver pident length qlen qstart qend evalue bitscore staxids sscinames

    Only the top hit per query (first line) is used for the taxonomy decision.
    Subsequent hits are ignored.
    """
    tb = cfg.get("taxonomy_blast", {})
    min_qcov = float(tb.get("min_qcov", 50.0))
    min_pid = float(tb.get("min_pident", 80.0))
    reject_cross_kingdom = tb.get("reject_cross_kingdom", True)

    seen_top: set = set()  # record first hit per query
    hits: Dict[str, dict] = {}  # qid -> parsed hit dict

    if not os.path.exists(path):
        return

    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            parts = line.rstrip().split("\t")
            if len(parts) < 11:
                continue
            qid, sacc, pid_s, length_s, qlen_s = parts[:5]
            evalue_s, bits_s, stax_s, sciname = parts[7], parts[8], parts[9], parts[10]

            if qid in seen_top:
                continue
            try:
                pid = float(pid_s)
                length = int(length_s)
                qlen = int(qlen_s)
            except ValueError:
                continue

            if pid < min_pid:
                continue
            # MMSeqs2 translated search (search_type=2, nuc→protein) reports alnlen
            # in amino-acid units while qlen is in nucleotides.  Multiply by 3 to
            # convert so that min_qcov is applied on a consistent nt-equivalent scale.
            search_type = int(tb.get("search_type", 1))
            eff_len = length * 3 if search_type == 2 else length
            qcov = eff_len / qlen * 100 if qlen else 0.0
            if qcov < min_qcov:
                continue

            seen_top.add(qid)
            hits[qid] = {
                "sacc": sacc,
                "pid": pid,
                "length": length,
                "qlen": qlen,
                "evalue": evalue_s,
                "bits": bits_s,
                "staxids": stax_s,
                "sciname": sciname.strip(),
            }

    for qid, a in ann.items():
        if qid not in hits:
            # Query was not in BLAST output (no hits above thresholds).
            a.add_reason("taxonomy_not_checked")
            continue

        h = hits[qid]
        a.taxonomy_top_hit = h["sacc"]
        a.taxonomy_top_name = h["sciname"]
        a.taxonomy_top_pident = str(h["pid"])
        a.taxonomy_top_length = str(h["length"])
        a.taxonomy_top_staxids = h["staxids"]

        hs = h["sciname"].lower()
        sp = a.species_guess.lower()
        gen = a.genus_guess.lower()

        # --- taxdump-based comparison (preferred when available) ---
        if taxdb and h["staxids"] and h["staxids"] != "N/A":
            try:
                staxid = int(h["staxids"].split(";")[0])
                hit_kingdom = taxdb.get_kingdom(staxid)
                hit_genus = taxdb.get_rank_name(staxid, "genus")

                # Cross-kingdom check.
                if reject_cross_kingdom and hit_kingdom and a.kingdom != "Unknown":
                    expected_kd = a.kingdom
                    # Map Spinner kingdom names to taxdump superkingdom names.
                    kd_map = {
                        "Animal": "Eukaryota",
                        "Plant": "Eukaryota",
                        "Fungi": "Eukaryota",
                        "Protist": "Eukaryota",
                        "Bacteria": "Bacteria",
                        "Archaea": "Archaea",
                    }
                    expected_super = kd_map.get(expected_kd, "")
                    # A cross-kingdom hit is one that does not share the same
                    # superkingdom as the expected kingdom.
                    if expected_super and hit_kingdom not in ("", expected_super):
                        a.taxonomy_status = "REJECT_CROSS_KINGDOM"
                        a.add_reason("taxonomy_cross_kingdom")
                        continue

                if gen and hit_genus.lower() == gen:
                    a.taxonomy_status = "PASS_GENUS"
                    a.add_reason("taxonomy_same_genus")
                elif sp and sp in hs:
                    a.taxonomy_status = "PASS_SPECIES"
                    a.add_reason("taxonomy_same_species")
                else:
                    a.taxonomy_status = "NO_EXPECTED_MATCH"
                    a.add_reason("taxonomy_no_expected_match")
                continue
            except (ValueError, AttributeError):
                pass  # fall through to string matching

        # --- string-matching fallback ---
        if sp and sp in hs:
            a.taxonomy_status = "PASS_SPECIES"
            a.add_reason("taxonomy_same_species")
        elif gen and gen in hs:
            a.taxonomy_status = "PASS_GENUS"
            a.add_reason("taxonomy_same_genus")
        else:
            a.taxonomy_status = "NO_EXPECTED_MATCH"
            a.add_reason("taxonomy_no_expected_match")


# ---------------------------------------------------------------------------
# Windowed BLAST (chimerism screen)
# ---------------------------------------------------------------------------

def parse_tax_blast_escalation(
    path: str,
    ann: Dict[str, Annotation],
    cfg: dict,
    taxdb: Optional[TaxdumpDB],
    source_label: str,
    nt_mode: bool = False,
) -> int:
    """Re-evaluate cross-kingdom rejections against a fallback database.

    For each record previously flagged ``taxonomy_cross_kingdom``, checks
    whether the escalation search found an expected-kingdom hit.  If so, the
    cross-kingdom reason is removed and the record is re-classified as
    ``ESCALATION_RESCUED`` with a ``taxonomy_rescued_{source_label}`` reason.

    *nt_mode* selects nucleotide-to-nucleotide comparison (NT BLAST) rather
    than translated nuc→protein (MMSeqs2 Swiss-Prot / NR).  The default
    pident threshold for NT mode is 80 % (megablast territory); for protein
    mode the value from ``taxonomy_blast.min_pident`` is used (default 30 %).

    Returns the count of rescued records.
    """
    tb = cfg.get("taxonomy_blast", {})
    min_qcov = float(tb.get("min_qcov", 50.0))
    default_pid = 80.0 if nt_mode else float(tb.get("min_pident", 30.0))
    min_pid = float(tb.get(
        f"escalation_{'nt' if nt_mode else 'nr'}_min_pident", default_pid
    ))
    search_type = 1 if nt_mode else int(tb.get("search_type", 1))

    if not os.path.exists(path):
        return 0

    seen_top: set = set()
    hits: Dict[str, dict] = {}

    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            parts = line.rstrip().split("\t")
            if len(parts) < 11:
                continue
            qid, sacc, pid_s, length_s, qlen_s = parts[:5]
            stax_s, sciname = parts[9], parts[10]
            if qid in seen_top:
                continue
            try:
                pid = float(pid_s)
                length = int(length_s)
                qlen = int(qlen_s)
            except ValueError:
                continue
            if pid < min_pid:
                continue
            eff_len = length * 3 if search_type == 2 else length
            qcov = eff_len / qlen * 100 if qlen else 0.0
            if qcov < min_qcov:
                continue
            seen_top.add(qid)
            hits[qid] = {
                "sacc": sacc, "pid": pid,
                "staxids": stax_s, "sciname": sciname.strip(),
            }

    kd_map = {
        "Animal": "Eukaryota", "Plant": "Eukaryota",
        "Fungi": "Eukaryota", "Protist": "Eukaryota",
        "Bacteria": "Bacteria", "Archaea": "Archaea",
    }

    rescued = 0
    for qid, a in ann.items():
        if "taxonomy_cross_kingdom" not in a.reasons:
            continue
        if qid not in hits:
            continue

        h = hits[qid]
        found_expected = False

        if taxdb and h["staxids"] and h["staxids"] not in ("N/A", "0", ""):
            try:
                staxid = int(h["staxids"].split(";")[0])
                hit_kingdom = taxdb.get_kingdom(staxid)
                expected_super = kd_map.get(a.kingdom, "")
                if expected_super and hit_kingdom == expected_super:
                    found_expected = True
            except (ValueError, AttributeError):
                hs = h["sciname"].lower()
                gen = a.genus_guess.lower()
                sp = a.species_guess.lower()
                if (gen and gen in hs) or (sp and sp in hs):
                    found_expected = True
        else:
            hs = h["sciname"].lower()
            gen = a.genus_guess.lower()
            sp = a.species_guess.lower()
            if (gen and gen in hs) or (sp and sp in hs):
                found_expected = True

        if found_expected:
            if "taxonomy_cross_kingdom" in a.reasons:
                a.reasons.remove("taxonomy_cross_kingdom")
            a.taxonomy_status = "ESCALATION_RESCUED"
            a.add_reason(f"taxonomy_rescued_{source_label}")
            a.taxonomy_top_hit = h["sacc"]
            a.taxonomy_top_name = h["sciname"]
            a.taxonomy_top_pident = str(h["pid"])
            a.taxonomy_top_staxids = h["staxids"]
            rescued += 1

    return rescued


def make_windowed_fasta(records, path: str, cfg: dict) -> int:
    """Write a sliding-window FASTA from long *records* for chimerism checks.

    Returns the number of windows written.
    """
    wb = cfg.get("windowed_blast", {})
    min_len = int(wb.get("enabled_for_min_length", 1000))
    win = int(wb.get("window_size", 500))
    step = int(wb.get("window_step", 250))
    n = 0

    with open(path, "w", encoding="utf-8") as out:
        for r in records:
            seq = r.seq_upper
            if len(seq) < min_len:
                continue
            last_start = max(0, len(seq) - win)
            starts = list(range(0, last_start + 1, step))
            if not starts or starts[-1] != last_start:
                starts.append(last_start)
            for wi, start in enumerate(starts, 1):
                sub = seq[start : start + win]
                # Skip very short tail windows.
                if len(sub) < min(100, max(1, win // 2)):
                    continue
                n += 1
                out.write(f">{r.id}|win{wi}|{start + 1}-{start + len(sub)}\n{sub}\n")

    return n


def parse_windowed_blast(
    path: str,
    ann: Dict[str, Annotation],
    cfg: dict,
    taxdb: Optional["TaxdumpDB"] = None,
) -> None:
    """Flag records whose windows have incompatible top-hit taxa.

    When *taxdb* is provided, windows are compared at the rank specified by
    ``windowed_blast.taxdump_comparison_rank`` (default: ``family``).  This
    gives accurate conflict detection across deeply divergent taxa without
    relying on string matching of scientific names.

    Without *taxdb*, comparison falls back to matching the first two words
    (genus + species) of each window's top-hit scientific name.
    """
    wb = cfg.get("windowed_blast", {})
    min_conf = int(wb.get("min_conflicting_windows", 2))
    action = wb.get("conflict_action", "review")
    min_pid = float(wb.get("min_pident", 80.0))
    comparison_rank = wb.get("taxdump_comparison_rank", "family")

    if not os.path.exists(path):
        return

    # Store (sciname, taxid_str) per window — take only the first (best) hit.
    first_by_window: Dict[str, Tuple[str, str]] = {}
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            parts = line.rstrip().split("\t")
            if len(parts) < 11:
                continue
            qid = parts[0]
            try:
                pid = float(parts[2])
            except ValueError:
                continue
            if pid < min_pid:
                continue
            sciname = parts[10].strip()
            taxid_str = parts[9].strip() if len(parts) > 9 else ""
            if qid not in first_by_window:
                first_by_window[qid] = (sciname, taxid_str)

    by_parent: Dict[str, List[Tuple[str, str]]] = defaultdict(list)
    for qid, hit in first_by_window.items():
        parent = qid.split("|", 1)[0]
        by_parent[parent].append(hit)

    for parent, hits in by_parent.items():
        if parent not in ann:
            continue
        a = ann[parent]
        if len(hits) < min_conf:
            a.windowed_status = "WINDOWED_OK"
            continue

        # Determine comparison labels (one per window).
        if taxdb is not None:
            labels: set = set()
            for sciname, taxid_str in hits:
                label = ""
                if taxid_str and taxid_str not in ("", "N/A", "0"):
                    try:
                        tid = int(taxid_str.split(";")[0])
                        label = taxdb.get_rank_name(tid, comparison_rank)
                    except (ValueError, AttributeError):
                        pass
                if not label:
                    # Fall back to genus-level string for this window.
                    label = " ".join((sciname or "").split()[:1]).lower()
                labels.add(label.lower())
            labels.discard("")
            is_conflict = len(labels) >= 2
        else:
            # String-matching fallback: compare first two words (genus species).
            taxa = {" ".join((nm or "").split()[:2]).lower() for nm, _ in hits if nm}
            taxa.discard("")
            is_conflict = len(taxa) >= 2

        if is_conflict:
            a.windowed_status = "WINDOWED_CONFLICT"
            a.add_reason("windowed_blast_conflict")
            if action == "reject":
                a.add_reason("taxonomy_cross_kingdom")
        else:
            a.windowed_status = "WINDOWED_OK"

    # Mark records that were long enough but had no windows hit.
    wb_min_len = int(wb.get("enabled_for_min_length", 1000))
    for key, a in ann.items():
        if a.windowed_status == "NOT_CHECKED" and a.length >= wb_min_len:
            a.windowed_status = "WINDOWED_OK"  # long but no BLAST hits


def load_taxdb(cfg: dict) -> Optional[TaxdumpDB]:
    """Load TaxdumpDB if ``taxonomy_blast.taxdump_dir`` is configured, else None."""
    td = cfg.get("taxonomy_blast", {}).get("taxdump_dir", "")
    if not td:
        return None
    try:
        return TaxdumpDB(td)
    except FileNotFoundError as e:
        warn(str(e))
        return None
