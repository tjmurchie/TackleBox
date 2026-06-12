from __future__ import annotations

import argparse
import datetime
import json
import sqlite3
import threading
import time
import urllib.parse
from dataclasses import asdict
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set

from .relic_lca import ReadAssignment, TaxonSummary, compute_damage_pattern
from .taxonomy import Taxonomy
from .viewer import _blank_metrics, build_viewer_payload

SCHEMA_VERSION = 2


def _j(x: object) -> str:
    return json.dumps(x, ensure_ascii=False, separators=(",", ":"))


def _nodes_dmp_mtime_iso(nodes_dmp_path: Optional[str | Path]) -> Optional[str]:
    """Return ISO-8601 mtime of nodes.dmp, or None if unavailable."""
    if not nodes_dmp_path:
        return None
    try:
        mtime = Path(nodes_dmp_path).stat().st_mtime
        return datetime.datetime.fromtimestamp(mtime, tz=datetime.timezone.utc).strftime("%Y-%m-%d")
    except OSError:
        return None


def write_viewer_sqlite(
    db_path: str | Path,
    taxonomy: Taxonomy,
    assignments: List[ReadAssignment],
    summaries: List[TaxonSummary],
    *,
    max_embed_reads: int = 0,
    nodes_dmp_path: Optional[str | Path] = None,
    deconvolution_rows: Optional[List[Dict]] = None,
    dark_taxids_tsv: Optional[Path] = None,
    run_config: Optional[Dict] = None,
) -> None:
    """Write an SQLite-backed Fillet viewer database.

    The DB stores the full taxonomy skeleton required for supported nodes, per-sample
    taxon metrics, and all read/alignment evidence. The HTML viewer served by
    `fillet serve-viewer` loads the tree first and fetches supporting reads only
    when a node is selected.
    """
    db_path = Path(db_path)
    db_path.parent.mkdir(parents=True, exist_ok=True)
    if db_path.exists():
        db_path.unlink()
    con = sqlite3.connect(str(db_path))
    cur = con.cursor()
    cur.executescript(
        """
        PRAGMA journal_mode=DELETE;
        CREATE TABLE meta(key TEXT PRIMARY KEY, value TEXT NOT NULL);
        CREATE TABLE taxon(
            taxid TEXT PRIMARY KEY,
            parent_taxid TEXT,
            name TEXT,
            rank TEXT
        );
        CREATE INDEX idx_taxon_parent ON taxon(parent_taxid);
        CREATE TABLE summary(
            taxid TEXT,
            sample_id TEXT,
            metrics_json TEXT NOT NULL,
            flags TEXT DEFAULT '',
            PRIMARY KEY(taxid, sample_id)
        );
        CREATE INDEX idx_summary_taxid ON summary(taxid);
        CREATE INDEX idx_summary_sample ON summary(sample_id);
        CREATE TABLE read_assignment(
            read_key TEXT PRIMARY KEY,
            read_id TEXT,
            sample_id TEXT,
            assigned_taxid TEXT,
            assigned_name TEXT,
            assigned_rank TEXT,
            posterior REAL,
            entropy REAL,
            status TEXT,
            reason TEXT,
            read_len INTEGER,
            read_sequence TEXT,
            payload_json TEXT NOT NULL
        );
        CREATE INDEX idx_read_taxid ON read_assignment(assigned_taxid);
        CREATE INDEX idx_read_sample ON read_assignment(sample_id);
        CREATE INDEX idx_read_taxid_post ON read_assignment(assigned_taxid, posterior DESC);
        CREATE TABLE deconvolution(
            clade_taxid TEXT,
            sample_id TEXT,
            species_taxid TEXT,
            species_name TEXT,
            proportion_mean REAL,
            proportion_ci_low REAL,
            proportion_ci_high REAL,
            bayes_factor REAL,
            supporting_reads INTEGER,
            prior_evidence TEXT,
            PRIMARY KEY(clade_taxid, sample_id, species_taxid)
        );
        CREATE INDEX idx_deconv_clade ON deconvolution(clade_taxid);
        CREATE TABLE IF NOT EXISTS dark_taxon(
            taxid TEXT PRIMARY KEY,
            reason TEXT NOT NULL,
            chunk_id TEXT
        );
        """
    )
    samples = sorted({a.sample_id for a in assignments} | {s.sample_id for s in summaries})
    meta_rows = [
        ("schema_version", str(SCHEMA_VERSION)),
        ("created_at", time.strftime("%Y-%m-%d %H:%M:%S")),
        ("root_taxid", taxonomy.root_taxid),
        ("samples", _j(samples)),
        ("viewer_note", "SQLite-backed Fillet viewer: reads are fetched on demand."),
    ]
    nodes_date = _nodes_dmp_mtime_iso(nodes_dmp_path)
    if nodes_date:
        meta_rows.append(("nodes_dmp_date", nodes_date))
    if run_config:
        meta_rows.append(("run_config", _j(run_config)))
    cur.executemany("INSERT INTO meta(key,value) VALUES(?,?)", meta_rows)

    supported: Set[str] = set()
    for s in summaries:
        if s.cumulative_hard_reads > 0 or s.cumulative_weighted_reads > 0:
            supported.update(taxonomy.ancestors(s.taxid, include_self=True))
    if taxonomy.root_taxid in taxonomy.taxa:
        supported.add(taxonomy.root_taxid)
    for a in assignments:
        if a.assigned_taxid != "0":
            supported.update(taxonomy.ancestors(a.assigned_taxid, include_self=True))

    tax_rows = []
    for t in supported:
        if not taxonomy.exists(t):
            continue
        tax_rows.append((t, taxonomy.parent(t), taxonomy.name(t), taxonomy.rank(t)))
    cur.executemany("INSERT OR REPLACE INTO taxon VALUES(?,?,?,?)", tax_rows)

    for s in summaries:
        metrics = {
            "hard_reads": s.cumulative_hard_reads,
            "direct_hard_reads": s.direct_hard_reads,
            "weighted_reads": s.cumulative_weighted_reads,
            "mean_posterior": s.mean_posterior,
            "conflicted_reads": s.conflicted_reads,
            "negative_weighted_reads": s.negative_weighted_reads,
            "blank_fraction": s.blank_fraction,
            "evidence_reads": s.evidence_reads,
            "unique_references": s.unique_references,
            "best_reference_breadth": s.best_reference_breadth,
            "mean_reference_breadth": s.mean_reference_breadth,
            "best_reference_span": s.best_reference_span,
            "max_stack_depth": s.max_stack_depth,
            "top_locus_fraction": s.top_locus_fraction,
            "stack_concentration": s.stack_concentration,
            "mean_damage_score": s.mean_damage_score,
            "max_damage_score": s.max_damage_score,
            "direct_weighted_reads": s.direct_weighted_reads,
            "mean_windows_per_ref": getattr(s, "mean_windows_per_ref", 0.0),
            "n_covered_windows": getattr(s, "n_covered_windows", 0),
            "composite_authenticity": getattr(s, "composite_authenticity", 0.0),
            "authenticity_tier": getattr(s, "authenticity_tier", None),
            "authenticity_badge": getattr(s, "authenticity_badge", None),
            "authenticity_tier_pch": getattr(s, "authenticity_tier_pch", 0),
            "reads_per_million": getattr(s, "reads_per_million", None),
            "weighted_per_million": getattr(s, "weighted_per_million", None),
            "relative_abundance": getattr(s, "relative_abundance", None),
            # Multi-proxy support flags (force bool — fields can be set/dict when
            # the upstream and-expression short-circuits on an empty collection)
            "eco_support": bool(getattr(s, "eco_support", False)),
            "pal_support": bool(getattr(s, "pal_support", False)),
            "fos_support": bool(getattr(s, "fos_support", False)),
            "fos_evidence_text": str(getattr(s, "fos_evidence_text", "") or ""),
            "n_support_lines": int(getattr(s, "n_support_lines", 0) or 0),
            # Experimental Bayesian fields
            "proportion_posterior_mean": getattr(s, "proportion_posterior_mean", 0.0),
            "proportion_ci_low": getattr(s, "proportion_ci_low", 0.0),
            "proportion_ci_high": getattr(s, "proportion_ci_high", 0.0),
            "estimated_damage_rate": getattr(s, "estimated_damage_rate", 0.0),
        }
        cur.execute(
            "INSERT OR REPLACE INTO summary(taxid,sample_id,metrics_json,flags) VALUES(?,?,?,?)",
            (s.taxid, s.sample_id, _j(metrics), s.flags),
        )

    for i, a in enumerate(assignments, start=1):
        try:
            hits = json.loads(a.alignment_evidence_json or "[]")
        except Exception:
            hits = []
        payload = {
            "read_id": a.read_id,
            "sample_id": a.sample_id,
            "assigned_taxid": a.assigned_taxid,
            "assigned_name": a.assigned_name,
            "assigned_rank": a.assigned_rank,
            "posterior": a.posterior,
            "entropy": a.entropy,
            "status": a.status,
            "reason": a.reason,
            "read_len": a.read_len,
            "read_sequence": a.read_sequence,
            "best_subject_id": a.best_subject_id,
            "best_hit_name": a.best_hit_name,
            "best_bitscore": a.best_bitscore,
            "best_pident": a.best_pident,
            "best_qcov": a.best_qcov,
            "best_sstart": a.best_sstart,
            "best_send": a.best_send,
            "best_slen": a.best_slen,
            "lca_name": a.lca_name,
            "lca_taxid": a.lca_taxid,
            "damage_ct_5p": a.damage_ct_5p,
            "damage_ga_3p": a.damage_ga_3p,
            "damage_score": a.damage_score,
            "em_refined": a.em_refined,
            "hits": hits,
        }
        key = f"{a.sample_id}::{a.read_id}::{i}"
        cur.execute(
            """INSERT INTO read_assignment
            (read_key,read_id,sample_id,assigned_taxid,assigned_name,assigned_rank,posterior,entropy,status,reason,read_len,read_sequence,payload_json)
            VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?)""",
            (key, a.read_id, a.sample_id, a.assigned_taxid, a.assigned_name, a.assigned_rank,
             a.posterior, a.entropy, a.status, a.reason, a.read_len, a.read_sequence, _j(payload)),
        )
    if deconvolution_rows:
        cur.executemany(
            """INSERT OR REPLACE INTO deconvolution
            (clade_taxid, sample_id, species_taxid, species_name,
             proportion_mean, proportion_ci_low, proportion_ci_high,
             bayes_factor, supporting_reads, prior_evidence)
            VALUES(?,?,?,?,?,?,?,?,?,?)""",
            [
                (
                    r.get("clade_taxid", ""),
                    r.get("sample_id", ""),
                    r.get("species_taxid", ""),
                    r.get("species_name", ""),
                    float(r.get("proportion_mean", 0.0)),
                    float(r.get("proportion_ci_low", 0.0)),
                    float(r.get("proportion_ci_high", 0.0)),
                    float(r.get("bayes_factor", 1.0)),
                    int(r.get("supporting_reads", 0)),
                    str(r.get("prior_evidence", "")),
                )
                for r in deconvolution_rows
            ],
        )

    if dark_taxids_tsv is not None:
        _dark_path = Path(dark_taxids_tsv)
        if _dark_path.exists():
            import csv
            with open(_dark_path, encoding="utf-8", newline="") as _fh:
                _reader = csv.DictReader(_fh, delimiter="\t")
                _dark_rows = [
                    (row["taxid"], row["reason"], row.get("chunk_id", ""))
                    for row in _reader
                ]
            if _dark_rows:
                cur.executemany(
                    "INSERT OR REPLACE INTO dark_taxon(taxid, reason, chunk_id) VALUES(?,?,?)",
                    _dark_rows,
                )

    # Cache virtual branch counts from Python (no SQL scan needed at load time).
    # Count from the assignments list we already have in memory.
    _vbc: Dict[str, int] = {"no_hits": 0, "em_demoted": 0, "lca_broad": 0, "rank_capped": 0}
    for _a in assignments:
        if _a.assigned_taxid != "0":
            continue
        _st = _a.status or ""
        if _st.startswith("unclassified:no_candidate_hits"):
            _vbc["no_hits"] += 1
        elif "em_coherence_lift:unclassified" in _st and "lca_too_broad" not in _st:
            _vbc["em_demoted"] += 1
        elif "lca_too_broad" in _st:
            _vbc["lca_broad"] += 1
        elif "rank_cap:unclassified" in _st:
            _vbc["rank_capped"] += 1
    cur.execute("INSERT OR REPLACE INTO meta(key,value) VALUES(?,?)",
                ("virtual_branch_counts", _j(_vbc)))

    con.commit()
    con.close()


def _dict_rows(cur: sqlite3.Cursor) -> List[Dict[str, object]]:
    cols = [d[0] for d in cur.description]
    return [dict(zip(cols, row)) for row in cur.fetchall()]


# ── Virtual branch (unclassified / excluded reads) ────────────────────────────

_VIRTUAL_TOOLTIPS: Dict[str, str] = {
    "no_hits": (
        "Reads that produced no database alignments. These reads were presented to the "
        "classifier but no BLAST/LAST hit was found against the reference database. "
        "Possible causes: (1) organism absent from the reference DB, (2) excessive "
        "sequence degradation beyond recognisable similarity, (3) non-biological sequence "
        "(adapters, PhiX spike-in), (4) novel or highly divergent organism. High counts "
        "here may indicate reference database gaps or library preparation issues."
    ),
    "em_demoted": (
        "Reads that initially had database hits but were removed by the EM coherence "
        "filter. The EM algorithm determined these reads were statistically inconsistent "
        "with the main assignment signal for their candidate taxon — likely stochastic "
        "cross-hits rather than genuine sequence similarity. This is the primary "
        "anti-contamination mechanism in RELIC-LCA. A high fraction of EM-demoted reads "
        "for a taxon suggests that taxon's signal is noise rather than genuine presence."
    ),
    "lca_broad": (
        "Reads whose EM-refined LCA (Lowest Common Ancestor) resolved to a taxonomic "
        "level too broad to be informative — for example 'root', 'cellular organisms', "
        "or a superkingdom. The assignment was suppressed to avoid artifactual inflation "
        "of high-level clades. Common causes: reads mapping to highly conserved genes "
        "shared across distant lineages (e.g. ribosomal RNA, core metabolic enzymes), "
        "or multi-kingdom cross-mapping from a mixed sample."
    ),
    "rank_capped": (
        "Reads whose best LCA was above the minimum reportable taxonomic rank. The "
        "classifier requires a minimum taxonomic resolution for an assignment to be "
        "reported. These reads had valid alignments, but the LCA resolved at a rank "
        "above the configured minimum (e.g., class or phylum level). Adjusting "
        "--lca-rank-cap or using a more targeted reference DB may recover these "
        "assignments at finer resolution."
    ),
    "contaminants": (
        "Taxa flagged as laboratory or environmental contaminants. These may have been "
        "identified by the classifier's internal contamination flags (user_contaminant), "
        "or manually marked via right-click in this viewer. Contaminant taxa remain "
        "visible in the main tree for traceability and are aggregated here for review.\n\n"
        "To mark a taxon as a contaminant: right-click it in the tree and choose "
        "Mark as → Contaminant."
    ),
    "root": (
        "A summary branch for reads that were not assigned to any specific taxon. "
        "Click each sub-category to browse the individual reads and their alignment "
        "details in the Inspector panel.\n\n"
        "• No alignments — no database hit found\n"
        "• EM demoted — coherence-filtered post-EM\n"
        "• LCA too broad — mapped to a non-informative clade\n"
        "• Rank-capped — LCA above minimum rank threshold\n"
        "• Contaminants — user-flagged or auto-detected contaminants"
    ),
}


def _build_virtual_branch(
    con: sqlite3.Connection,
    meta: Optional[Dict[str, str]] = None,
) -> Optional[Dict[str, object]]:
    """Build the synthetic __unclassified__ branch from the read_assignment table.

    Uses counts cached in ``meta`` when available (written at classification time).
    Falls back to SQL COUNT queries for old files, but skips entirely on large
    databases (>500 MB) where unindexed LIKE scans would hang over NFS.

    Returns None if the table is absent, the status column is missing (old schema),
    or there are no unclassified/contaminant reads at all.
    """
    # ── Fast path: use counts cached in meta at write time ────────────────────
    if meta:
        _cached_json = meta.get("virtual_branch_counts")
        if _cached_json:
            try:
                _c = json.loads(_cached_json)
                n_no_hits     = int(_c.get("no_hits",     0))
                n_em_demoted  = int(_c.get("em_demoted",  0))
                n_lca_broad   = int(_c.get("lca_broad",   0))
                n_rank_capped = int(_c.get("rank_capped", 0))
                # Contaminants are user-flagged at runtime; query is cheap (summary only)
                n_contaminants = 0
                try:
                    contam_rows = con.execute(
                        "SELECT COUNT(*) FROM summary WHERE flags LIKE '%user_contaminant%'"
                    ).fetchone()
                    n_contaminants = int(contam_rows[0]) if contam_rows else 0
                except Exception:
                    pass
                total_unclassified = n_no_hits + n_em_demoted + n_lca_broad + n_rank_capped
                if total_unclassified == 0 and n_contaminants == 0:
                    return None
                # Skip building contaminant detail nodes — just use the total
                return _build_virtual_nodes(
                    n_no_hits, n_em_demoted, n_lca_broad, n_rank_capped, n_contaminants
                )
            except Exception:
                pass  # fall through to SQL path

    # ── Slow path: SQL COUNT queries (old files without cached counts) ─────────
    # Skip if the database is too large to scan quickly over NFS.
    try:
        _page_count = con.execute("PRAGMA page_count").fetchone()[0]
        _page_size  = con.execute("PRAGMA page_size").fetchone()[0]
        if _page_count * _page_size > 500 * 1024 * 1024:
            return None  # >500 MB — skip to avoid hanging on NFS
    except Exception:
        pass

    try:
        cols = {r[1] for r in con.execute("PRAGMA table_info(read_assignment)").fetchall()}
        if "status" not in cols:
            return None
    except Exception:
        return None

    def _count(where: str) -> int:
        try:
            row = con.execute(f"SELECT COUNT(*) FROM read_assignment WHERE {where}").fetchone()
            return int(row[0]) if row else 0
        except Exception:
            return 0

    n_no_hits     = _count("assigned_taxid='0' AND status LIKE 'unclassified:no_candidate_hits%'")
    n_em_demoted  = _count(
        "assigned_taxid='0' AND status LIKE '%em_coherence_lift:unclassified'"
        " AND status NOT LIKE '%lca_too_broad%'"
    )
    n_lca_broad   = _count("assigned_taxid='0' AND status LIKE '%lca_too_broad%'")
    n_rank_capped = _count("assigned_taxid='0' AND status LIKE '%rank_cap:unclassified%'")

    # Contaminant taxa: user_contaminant flag in the summary flags column
    try:
        contam_rows = con.execute(
            "SELECT s.taxid, COUNT(ra.read_key) as n "
            "FROM summary s "
            "LEFT JOIN read_assignment ra ON ra.assigned_taxid = s.taxid "
            "WHERE s.flags LIKE '%user_contaminant%' "
            "GROUP BY s.taxid"
        ).fetchall()
        n_contaminants = sum(int(r["n"] or 0) for r in contam_rows)
    except Exception:
        n_contaminants = 0

    total_unclassified = n_no_hits + n_em_demoted + n_lca_broad + n_rank_capped
    if total_unclassified == 0 and n_contaminants == 0:
        return None

    return _build_virtual_nodes(n_no_hits, n_em_demoted, n_lca_broad, n_rank_capped, n_contaminants)


def _build_virtual_nodes(
    n_no_hits: int,
    n_em_demoted: int,
    n_lca_broad: int,
    n_rank_capped: int,
    n_contaminants: int,
) -> Optional[Dict[str, object]]:
    """Assemble the virtual __unclassified__ tree node from pre-computed counts."""
    total_unclassified = n_no_hits + n_em_demoted + n_lca_broad + n_rank_capped
    if total_unclassified == 0 and n_contaminants == 0:
        return None

    def _vnode(taxid: str, name: str, category: str, hard_reads: int,
               children: Optional[List] = None) -> Dict[str, object]:
        return {
            "taxid": taxid,
            "name": name,
            "rank": "virtual",
            "is_virtual_node": True,
            "virtual_category": category,
            "virtual_tooltip": _VIRTUAL_TOOLTIPS.get(category, ""),
            "metrics": {"hard_reads": hard_reads, "weighted_reads": 0.0,
                        "direct_hard_reads": hard_reads},
            "flags": [],
            "children": children or [],
        }

    children: List[Dict[str, object]] = [
        _vnode("__no_hits__",      "No alignments", "no_hits",     n_no_hits),
        _vnode("__em_demoted__",   "EM demoted",    "em_demoted",  n_em_demoted),
        _vnode("__lca_broad__",    "LCA too broad", "lca_broad",   n_lca_broad),
        _vnode("__rank_capped__",  "Rank-capped",   "rank_capped", n_rank_capped),
        _vnode("__contaminants__", "Contaminants",  "contaminants", n_contaminants),
    ]
    return _vnode(
        "__unclassified__",
        "Unclassified & Excluded",
        "root",
        total_unclassified,
        children=children,
    )


class ViewerDB:
    def __init__(self, path: str | Path):
        self.path = str(path)
        import urllib.parse as _up
        _uri = "file:" + _up.quote(str(path), safe="/:@") + "?mode=ro"
        try:
            self.con = sqlite3.connect(_uri, uri=True, check_same_thread=False)
        except Exception:
            # Fallback: some SQLite builds / NFS mounts reject URI mode
            self.con = sqlite3.connect(str(path), check_same_thread=False)
        self.con.row_factory = sqlite3.Row
        # Read-path tuning: large page cache + temp tables in memory
        try:
            self.con.execute("PRAGMA cache_size=-131072")   # 128 MB
            self.con.execute("PRAGMA temp_store=2")          # MEMORY
            self.con.execute("PRAGMA query_only=1")
        except Exception:
            pass
        self._children_cache: Dict[str, List[str]] | None = None
        self._lock = threading.Lock()  # serialise con.execute() across threads

    def meta(self) -> Dict[str, str]:
        rows = self.con.execute("SELECT key,value FROM meta").fetchall()
        return {r["key"]: r["value"] for r in rows}

    def taxonomy_age_warning(self, warn_after_days: int = 180) -> Optional[str]:
        """Return a warning string if the stored taxonomy is older than warn_after_days, else None."""
        nodes_date = self.meta().get("nodes_dmp_date")
        if not nodes_date:
            return None
        try:
            dt = datetime.datetime.strptime(nodes_date, "%Y-%m-%d").date()
        except ValueError:
            return None
        age_days = (datetime.date.today() - dt).days
        if age_days >= warn_after_days:
            months = age_days // 30
            return (
                f"Taxonomy database is {months} months old (nodes.dmp from {nodes_date}). "
                "Consider running `fillet update-taxonomy` to download a fresh NCBI taxonomy."
            )
        return None

    def children_map(self) -> Dict[str, List[str]]:
        if self._children_cache is not None:
            return self._children_cache
        with self._lock:
            if self._children_cache is None:   # double-check after acquiring
                m: Dict[str, List[str]] = {}
                rows = self.con.execute("SELECT taxid,parent_taxid FROM taxon").fetchall()
                for r in rows:
                    taxid = str(r["taxid"])
                    parent = str(r["parent_taxid"])
                    if parent == taxid:
                        continue
                    m.setdefault(parent, []).append(taxid)
                self._children_cache = m
        return self._children_cache

    def descendants(self, taxid: str, include_self: bool = True) -> Set[str]:
        out: Set[str] = set([taxid]) if include_self else set()
        stack = list(self.children_map().get(taxid, []))
        while stack:
            t = stack.pop()
            if t in out:
                continue
            out.add(t)
            stack.extend(self.children_map().get(t, []))
        return out

    def payload(self) -> Dict[str, object]:
        meta = self.meta()
        root = meta.get("root_taxid", "1")
        samples = json.loads(meta.get("samples", "[]"))
        taxa = {str(r["taxid"]): dict(r) for r in self.con.execute("SELECT * FROM taxon")}
        summaries = self.con.execute("SELECT * FROM summary").fetchall()
        by_tax = {tid: _blank_metrics() for tid in taxa}
        flags: Dict[str, Set[str]] = {tid: set() for tid in taxa}
        for r in summaries:
            tid = str(r["taxid"]); sid = str(r["sample_id"])
            metrics = json.loads(r["metrics_json"])
            m = by_tax.setdefault(tid, _blank_metrics())
            for k in ["hard_reads", "direct_hard_reads", "weighted_reads", "direct_weighted_reads",
                      "conflicted_reads", "negative_weighted_reads", "evidence_reads"]:
                m[k] = (m.get(k, 0) or 0) + metrics.get(k, 0)
            for k in ["unique_references", "max_stack_depth", "n_covered_windows"]:
                m[k] = max(m.get(k, 0) or 0, metrics.get(k, 0) or 0)
            for k in ["mean_posterior", "blank_fraction", "best_reference_breadth",
                      "mean_reference_breadth", "top_locus_fraction",
                      "mean_damage_score", "max_damage_score", "stack_concentration",
                      "estimated_damage_rate", "mean_windows_per_ref", "composite_authenticity"]:
                m[k] = max(float(m.get(k) or 0), float(metrics.get(k) or 0))
            # Best (lowest numbered) authenticity tier across samples
            tier = metrics.get("authenticity_tier")
            if tier is not None:
                cur_tier = m.get("authenticity_tier")
                if cur_tier is None or int(tier) < int(cur_tier):
                    m["authenticity_tier"] = int(tier)
                    m["authenticity_badge"] = metrics.get("authenticity_badge") or ""
            # Normalised read counts — sum across samples
            for k in ("reads_per_million", "weighted_per_million", "relative_abundance"):
                v = metrics.get(k)
                if v is not None:
                    m[k] = (m.get(k) or 0) + float(v)
            m.setdefault("by_sample", {})[sid] = metrics
            for f in str(r["flags"] or "").split(";"):
                if f:
                    flags.setdefault(tid, set()).add(f)
        children = self.children_map()
        def node(tid: str) -> Dict[str, object]:
            tr = taxa[tid]
            kids = [c for c in children.get(tid, []) if c in taxa]
            kids.sort(key=lambda c: float(by_tax.get(c, {}).get("weighted_reads", 0) or 0), reverse=True)
            return {
                "taxid": tid,
                "name": tr.get("name") or tid,
                "rank": tr.get("rank") or "no rank",
                "metrics": by_tax.get(tid, _blank_metrics()),
                "flags": sorted(flags.get(tid, set())),
                "children": [node(c) for c in kids],
            }
        if root not in taxa:
            root = next(iter(taxa)) if taxa else "1"
        tree = node(root)
        virtual = _build_virtual_branch(self.con, meta=meta)
        if virtual is not None:
            tree["children"].append(virtual)
        return {"root_taxid": root, "samples": samples, "tree": tree, "reads": [], "viewer_note": meta.get("viewer_note", "")}

    def reads(
        self,
        taxids: Iterable[str],
        sample: str = "__ALL__",
        limit: int = 250,
        descendants: bool = True,
        ordered: bool = True,
        max_desc: int = 0,
    ) -> List[Dict[str, object]]:
        """Fetch reads assigned under *taxids*.

        Parameters
        ----------
        ordered   : sort results by posterior DESC then read_len DESC.
                    Pass False for fast unordered retrieval (no full-scan sort).
                    When False, the in-memory final sort is also skipped.
        max_desc  : cap the number of descendant taxids consulted (0 = no cap).
                    Useful when the caller only needs a representative sample and
                    the selected node has thousands of descendants (avoids huge
                    multi-chunk queries that are slow to sort on large DBs).
        """
        # Route virtual-branch taxids (prefixed with __) to their own SQL logic
        taxid_list = list(taxids)
        if any(str(t).startswith("__") for t in taxid_list):
            result: List[Dict[str, object]] = []
            for vtid in taxid_list:
                result.extend(self._reads_virtual(str(vtid), sample, limit))
            if ordered:
                result.sort(key=lambda r: (-(r.get("posterior") or 0), -(r.get("read_len") or 0)))
            return result[:limit]

        wanted: Set[str] = set()
        for t in taxid_list:
            if descendants:
                desc = self.descendants(str(t), include_self=True)
                if max_desc > 0 and len(desc) > max_desc:
                    desc = set(list(desc)[:max_desc])
                wanted.update(desc)
            else:
                wanted.add(str(t))
        if not wanted:
            return []
        wanted_list = list(wanted)
        # SQLite SQLITE_LIMIT_VARIABLE_NUMBER defaults to 999; chunk IN clause to stay safe.
        _CHUNK = 990
        all_rows: list = []
        for i in range(0, len(wanted_list), _CHUNK):
            chunk = wanted_list[i : i + _CHUNK]
            ph = ",".join("?" for _ in chunk)
            params: List[object] = list(chunk)
            sql = f"SELECT payload_json FROM read_assignment WHERE assigned_taxid IN ({ph})"
            if sample and sample != "__ALL__":
                sql += " AND sample_id=?"
                params.append(sample)
            if ordered:
                sql += " ORDER BY posterior DESC, read_len DESC"
            sql += " LIMIT ?"
            params.append(int(limit))
            with self._lock:
                rows = self.con.execute(sql, params).fetchall()
            all_rows.extend(rows)
            if len(all_rows) >= limit * 2:
                break
        result = [json.loads(r["payload_json"]) for r in all_rows]
        if ordered:
            result.sort(key=lambda r: (-(r.get("posterior") or 0), -(r.get("read_len") or 0)))
        return result[:limit]

    def _reads_virtual(self, vtid: str, sample: str, limit: int) -> List[Dict[str, object]]:
        """Fetch reads for a virtual-branch category node (taxid starts with __)."""
        _STATUS_CONDS = {
            "__no_hits__":      "assigned_taxid='0' AND status LIKE 'unclassified:no_candidate_hits%'",
            "__em_demoted__":   "assigned_taxid='0' AND status LIKE '%em_coherence_lift:unclassified'"
                                " AND status NOT LIKE '%lca_too_broad%'",
            "__lca_broad__":    "assigned_taxid='0' AND status LIKE '%lca_too_broad%'",
            "__rank_capped__":  "assigned_taxid='0' AND status LIKE '%rank_cap:unclassified%'",
            "__unclassified__": "assigned_taxid='0'",
        }
        if vtid == "__contaminants__":
            try:
                ctids = [
                    str(r["taxid"])
                    for r in self.con.execute(
                        "SELECT DISTINCT taxid FROM summary WHERE flags LIKE '%user_contaminant%'"
                    ).fetchall()
                ]
            except Exception:
                return []
            if not ctids:
                return []
            ph = ",".join("?" for _ in ctids)
            sql = f"SELECT payload_json FROM read_assignment WHERE assigned_taxid IN ({ph})"
            params: List[object] = list(ctids)
            if sample and sample != "__ALL__":
                sql += " AND sample_id=?"
                params.append(sample)
            sql += " ORDER BY posterior DESC LIMIT ?"
            params.append(limit)
            with self._lock:
                rows = self.con.execute(sql, params).fetchall()
            return [json.loads(r["payload_json"]) for r in rows]

        cond = _STATUS_CONDS.get(vtid)
        if not cond:
            return []
        params_plain: List[object] = []
        sql = f"SELECT payload_json FROM read_assignment WHERE {cond}"
        if sample and sample != "__ALL__":
            sql += " AND sample_id=?"
            params_plain.append(sample)
        sql += " ORDER BY posterior DESC LIMIT ?"
        params_plain.append(limit)
        with self._lock:
            try:
                rows = self.con.execute(sql, params_plain).fetchall()
            except Exception:
                return []
        return [json.loads(r["payload_json"]) for r in rows]

    def coverage(self, taxids: Iterable[str], sample: str = "__ALL__", limit: int = 10000) -> Dict[str, object]:
        """Return per-reference coverage intervals for reads assigned under taxids."""
        reads = self.reads(taxids, sample=sample, limit=limit, descendants=True)
        subjects: Dict[str, Dict[str, object]] = {}
        for r in reads:
            sid = r.get("best_subject_id") or ""
            if not sid or not r.get("best_sstart") or not r.get("best_send"):
                continue
            if sid not in subjects:
                subjects[sid] = {
                    "subject_id": sid,
                    "name": r.get("best_hit_name") or r.get("best_subject_name") or sid,
                    "slen": r.get("best_slen") or 0,
                    "intervals": [],
                }
            subjects[sid]["intervals"].append({
                "a": int(min(r["best_sstart"], r["best_send"])),
                "b": int(max(r["best_sstart"], r["best_send"])),
                "posterior": float(r.get("posterior") or 0),
                "read_id": r.get("read_id") or "",
            })
        # Sort by interval count descending, return top 10
        sorted_subs = sorted(subjects.values(), key=lambda x: len(x["intervals"]), reverse=True)[:10]
        return {"subjects": sorted_subs, "total_reads": len(reads)}

    def damage_pattern(self, taxids: Iterable[str], sample: str = "__ALL__", n_pos: int = 15, limit: int = 10000) -> Dict[str, object]:
        """Compute per-position C→T (5') and G→A (3') damage rates across reads.

        Reads the qseq/sseq fields stored in each read's hits (alignment_evidence_json)
        and accumulates per-position substitution rates for mapDamage-style visualisation.
        Returns 0-rate arrays when no qseq/sseq data is available (e.g. cap40 BLAST).
        """
        reads = self.reads(taxids, sample=sample, limit=limit, descendants=True)
        alignments = []
        for r in reads:
            # Use damage_ct_5p/ga_3p from the top-level read fields first;
            # extract qseq/sseq from stored hits for per-position computation
            for hit in r.get("hits") or []:
                qseq = hit.get("qseq") or ""
                sseq = hit.get("sseq") or ""
                if qseq and sseq:
                    alignments.append((qseq, sseq))
                    break  # use best hit only (hits ordered by bitscore)
        pattern = compute_damage_pattern(alignments, n_pos=n_pos)
        pattern["total_reads"] = len(reads)
        return pattern

    def extract(self, taxids: Iterable[str], sample: str = "__ALL__", fmt: str = "fasta", limit: int = 100000) -> tuple[str, str]:
        reads = self.reads(taxids, sample=sample, limit=limit, descendants=True)
        if fmt == "tsv":
            lines = ["sample_id\tread_id\tassigned_taxid\tassigned_name\tposterior\tread_sequence"]
            for r in reads:
                seq = str(r.get("read_sequence") or "")
                lines.append("\t".join([str(r.get("sample_id","")), str(r.get("read_id","")), str(r.get("assigned_taxid","")), str(r.get("assigned_name","")), str(r.get("posterior",0)), seq]))
            return "text/tab-separated-values", "\n".join(lines) + "\n"
        lines = []
        for r in reads:
            seq = str(r.get("read_sequence") or "")
            if not seq:
                continue
            header = f">{r.get('sample_id','')}|{r.get('read_id','')}|taxid={r.get('assigned_taxid','')}|posterior={float(r.get('posterior') or 0):.4f}"
            lines.append(header)
            lines.append(seq)
        return "text/plain", "\n".join(lines) + ("\n" if lines else "")


def merge_viewer_dbs(
    input_paths: List[str | Path],
    output_path: str | Path,
    select_samples: List[str] = None,
    exclude_samples: List[str] = None,
) -> Dict[str, object]:
    """Merge two or more Fillet viewer SQLite databases into a single output database.

    Each source database contributes its own samples. If the same sample_id appears
    in more than one source the first occurrence wins (later ones are skipped and
    reported in the returned ``conflicts`` list).

    Returns a dict: ``{"samples": [...], "conflicts": [...], "n_sources": N}``
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.exists():
        output_path.unlink()

    out_con = sqlite3.connect(str(output_path))
    out_cur = out_con.cursor()
    out_cur.executescript(
        """
        PRAGMA journal_mode=DELETE;
        CREATE TABLE meta(key TEXT PRIMARY KEY, value TEXT NOT NULL);
        CREATE TABLE taxon(
            taxid TEXT PRIMARY KEY,
            parent_taxid TEXT,
            name TEXT,
            rank TEXT
        );
        CREATE INDEX idx_taxon_parent ON taxon(parent_taxid);
        CREATE TABLE summary(
            taxid TEXT,
            sample_id TEXT,
            metrics_json TEXT NOT NULL,
            flags TEXT DEFAULT '',
            PRIMARY KEY(taxid, sample_id)
        );
        CREATE INDEX idx_summary_taxid ON summary(taxid);
        CREATE INDEX idx_summary_sample ON summary(sample_id);
        CREATE TABLE read_assignment(
            read_key TEXT PRIMARY KEY,
            read_id TEXT,
            sample_id TEXT,
            assigned_taxid TEXT,
            assigned_name TEXT,
            assigned_rank TEXT,
            posterior REAL,
            entropy REAL,
            status TEXT,
            reason TEXT,
            read_len INTEGER,
            read_sequence TEXT,
            payload_json TEXT NOT NULL
        );
        CREATE INDEX idx_read_taxid ON read_assignment(assigned_taxid);
        CREATE INDEX idx_read_sample ON read_assignment(sample_id);
        CREATE INDEX idx_read_taxid_post ON read_assignment(assigned_taxid, posterior DESC);
        CREATE TABLE deconvolution(
            clade_taxid TEXT,
            sample_id TEXT,
            species_taxid TEXT,
            species_name TEXT,
            proportion_mean REAL,
            proportion_ci_low REAL,
            proportion_ci_high REAL,
            bayes_factor REAL,
            supporting_reads INTEGER,
            prior_evidence TEXT,
            PRIMARY KEY(clade_taxid, sample_id, species_taxid)
        );
        CREATE INDEX idx_deconv_clade ON deconvolution(clade_taxid);
        CREATE TABLE dark_taxon(
            taxid TEXT PRIMARY KEY,
            reason TEXT NOT NULL,
            chunk_id TEXT
        );
        """
    )

    all_samples: List[str] = []
    seen_samples: Set[str] = set()
    conflicts: List[str] = []
    root_taxid = "1"

    for i, src_path in enumerate(input_paths):
        src_con = sqlite3.connect(str(Path(src_path)))
        meta = dict(src_con.execute("SELECT key, value FROM meta").fetchall())
        src_samples: List[str] = json.loads(meta.get("samples", "[]"))

        if i == 0:
            root_taxid = meta.get("root_taxid", "1")

        # Apply sample filters before registering samples from this source
        if select_samples:
            src_samples = [s for s in src_samples if s in select_samples]
        if exclude_samples:
            src_samples = [s for s in src_samples if s not in exclude_samples]

        for sid in src_samples:
            if sid not in seen_samples:
                all_samples.append(sid)
                seen_samples.add(sid)
            else:
                conflicts.append(sid)

        # The set of sample IDs that passed the filter for this source DB
        active_samples: Set[str] = set(src_samples)

        # taxon — INSERT OR IGNORE so the first DB's taxonomy wins for conflicts
        out_cur.executemany(
            "INSERT OR IGNORE INTO taxon(taxid, parent_taxid, name, rank) VALUES(?,?,?,?)",
            src_con.execute("SELECT taxid, parent_taxid, name, rank FROM taxon").fetchall(),
        )

        # summary — keyed by (taxid, sample_id); only copy rows for active samples
        rows = [
            r for r in src_con.execute(
                "SELECT taxid, sample_id, metrics_json, flags FROM summary"
            ).fetchall()
            if r[1] in active_samples
        ]
        out_cur.executemany(
            "INSERT OR IGNORE INTO summary(taxid, sample_id, metrics_json, flags) VALUES(?,?,?,?)",
            rows,
        )

        # read_assignment — keyed by read_key = sample_id:read_id; only copy active samples
        rows = [
            r for r in src_con.execute(
                """SELECT read_key, read_id, sample_id, assigned_taxid, assigned_name,
                   assigned_rank, posterior, entropy, status, reason,
                   read_len, read_sequence, payload_json
                   FROM read_assignment"""
            ).fetchall()
            if r[2] in active_samples  # r[2] = sample_id column
        ]
        out_cur.executemany(
            """INSERT OR IGNORE INTO read_assignment(
                read_key, read_id, sample_id, assigned_taxid, assigned_name,
                assigned_rank, posterior, entropy, status, reason,
                read_len, read_sequence, payload_json
            ) VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?)""",
            rows,
        )

        # deconvolution — only copy rows for active samples (if table exists in source)
        src_tables = {r[0] for r in src_con.execute(
            "SELECT name FROM sqlite_master WHERE type='table'"
        ).fetchall()}
        if "deconvolution" in src_tables:
            deconv_rows = [
                r for r in src_con.execute(
                    """SELECT clade_taxid, sample_id, species_taxid, species_name,
                       proportion_mean, proportion_ci_low, proportion_ci_high,
                       bayes_factor, supporting_reads, prior_evidence
                       FROM deconvolution"""
                ).fetchall()
                if r[1] in active_samples
            ]
            if deconv_rows:
                out_cur.executemany(
                    """INSERT OR IGNORE INTO deconvolution(
                        clade_taxid, sample_id, species_taxid, species_name,
                        proportion_mean, proportion_ci_low, proportion_ci_high,
                        bayes_factor, supporting_reads, prior_evidence
                    ) VALUES(?,?,?,?,?,?,?,?,?,?)""",
                    deconv_rows,
                )

        # dark_taxon — INSERT OR IGNORE (first source wins for duplicate taxids)
        if "dark_taxon" in src_tables:
            dark_rows = src_con.execute(
                "SELECT taxid, reason, chunk_id FROM dark_taxon"
            ).fetchall()
            if dark_rows:
                out_cur.executemany(
                    "INSERT OR IGNORE INTO dark_taxon(taxid, reason, chunk_id) VALUES(?,?,?)",
                    dark_rows,
                )

        src_con.close()

    out_cur.executemany(
        "INSERT INTO meta(key, value) VALUES(?,?)",
        [
            ("schema_version", str(SCHEMA_VERSION)),
            ("created_at", time.strftime("%Y-%m-%d %H:%M:%S")),
            ("root_taxid", root_taxid),
            ("samples", json.dumps(all_samples)),
            ("viewer_note", f"Merged from {len(input_paths)} databases."),
            ("merged_sources", json.dumps([str(p) for p in input_paths])),
        ],
    )
    out_con.commit()
    out_con.close()

    return {"samples": all_samples, "conflicts": conflicts, "n_sources": len(input_paths)}


def serve_viewer(db_path: str | Path, host: str = "127.0.0.1", port: int = 8765) -> None:
    from importlib import resources
    db = ViewerDB(db_path)
    template = resources.files("fillet.data").joinpath("viewer_template.html").read_text(encoding="utf-8")
    html = template.replace("__FILLET_DATA_PAYLOAD__", "null")
    html = html.replace("const DATA_EMBEDDED_MODE = true;", "const DATA_EMBEDDED_MODE = false;")

    class Handler(BaseHTTPRequestHandler):
        def _send(self, code: int, body: str | bytes, ctype: str = "application/json") -> None:
            if isinstance(body, str):
                body = body.encode("utf-8")
            self.send_response(code)
            self.send_header("Content-Type", ctype + "; charset=utf-8")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)

        def log_message(self, fmt: str, *args: object) -> None:
            print("[Fillet viewer] " + fmt % args)

        def do_GET(self) -> None:  # noqa: N802
            url = urllib.parse.urlparse(self.path)
            qs = urllib.parse.parse_qs(url.query)
            try:
                if url.path in {"/", "/viewer"}:
                    self._send(200, html, "text/html")
                elif url.path == "/api/payload":
                    self._send(200, _j(db.payload()))
                elif url.path == "/api/reads":
                    taxids = qs.get("taxid", []) + qs.get("taxids", [])
                    if len(taxids) == 1 and "," in taxids[0]:
                        taxids = [x for x in taxids[0].split(",") if x]
                    sample = qs.get("sample", ["__ALL__"])[0]
                    limit = int(qs.get("limit", ["250"])[0])
                    self._send(200, _j({"reads": db.reads(taxids, sample=sample, limit=limit)}))
                elif url.path == "/api/coverage":
                    taxids = qs.get("taxid", []) + qs.get("taxids", [])
                    if len(taxids) == 1 and "," in taxids[0]:
                        taxids = [x for x in taxids[0].split(",") if x]
                    sample = qs.get("sample", ["__ALL__"])[0]
                    limit = int(qs.get("limit", ["5000"])[0])
                    self._send(200, _j(db.coverage(taxids, sample=sample, limit=limit)))
                elif url.path == "/api/damage_pattern":
                    taxids = qs.get("taxid", []) + qs.get("taxids", [])
                    if len(taxids) == 1 and "," in taxids[0]:
                        taxids = [x for x in taxids[0].split(",") if x]
                    sample = qs.get("sample", ["__ALL__"])[0]
                    n_pos = int(qs.get("n_pos", ["15"])[0])
                    limit = int(qs.get("limit", ["10000"])[0])
                    self._send(200, _j(db.damage_pattern(taxids, sample=sample, n_pos=n_pos, limit=limit)))
                elif url.path == "/api/extract":
                    taxids = qs.get("taxid", []) + qs.get("taxids", [])
                    if len(taxids) == 1 and "," in taxids[0]:
                        taxids = [x for x in taxids[0].split(",") if x]
                    sample = qs.get("sample", ["__ALL__"])[0]
                    fmt = qs.get("format", ["fasta"])[0]
                    ctype, body = db.extract(taxids, sample=sample, fmt=fmt)
                    self._send(200, body, ctype)
                else:
                    self._send(404, _j({"error": "not found"}))
            except Exception as e:
                self._send(500, _j({"error": str(e)}))

    httpd = ThreadingHTTPServer((host, port), Handler)
    print(f"[Fillet viewer] Serving {db_path}")
    print(f"[Fillet viewer] Open http://{host}:{port}/")
    httpd.serve_forever()
