"""Clustering support: vsearch UC file parsing, cluster run, and uchime chimera detection."""
from __future__ import annotations

import shutil
from typing import Dict

from .annotation import Annotation
from .external import run
from .utils import BlastTicker, info, warn


def parse_uc(path: str, ann: Dict[str, Annotation]) -> None:
    """Parse a vsearch UC file and annotate cluster membership in *ann*.

    UC format columns (tab-separated):
      0: type  (S=seed/centroid, H=hit/member, C=cluster summary)
      1: cluster number
      2: length or identity
      8: query label
      9: target label (centroid id, or * for seeds)

    Centroids receive cluster_role="centroid" and the cluster_representative
    reason; non-centroid members receive cluster_role="member" and
    cluster_nonrepresentative.
    """
    if not path:
        return
    cent: Dict[str, str] = {}   # cluster_number -> centroid_id
    q2c: Dict[str, str] = {}    # query_id -> cluster_number

    try:
        with open(path, encoding="utf-8", errors="replace") as f:
            for line in f:
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.rstrip().split("\t")
                if len(parts) < 9:
                    continue
                typ, cid, qid = parts[0], parts[1], parts[8]
                q2c[qid] = cid
                if typ == "S":
                    cent[cid] = qid
    except FileNotFoundError:
        warn(f"Clustering UC file not found: {path}")
        return

    for key, a in ann.items():
        if key not in q2c:
            continue
        cid = q2c[key]
        a.cluster_id = cid
        if cent.get(cid) == key:
            a.cluster_role = "centroid"
            a.add_reason("cluster_representative")
        else:
            a.cluster_role = "member"
            a.add_reason("cluster_nonrepresentative")


def run_vsearch(
    input_fasta: str,
    outprefix: str,
    cfg: dict,
    ann: Dict[str, Annotation],
    total_queries: int = 0,
) -> None:
    """Run vsearch clustering on *input_fasta* and annotate *ann*."""
    cl = cfg.get("cluster", {})
    vsearch_path = cl.get("vsearch_path", "vsearch")

    if not shutil.which(vsearch_path):
        msg = f"vsearch not found at {vsearch_path!r}; clustering skipped."
        if cfg.get("run", {}).get("fail_on_missing_external_tool", False):
            raise RuntimeError(msg)
        warn(msg)
        return

    uc_path = outprefix + ".clusters.uc"
    cen_path = outprefix + ".cluster_centroids.fasta"
    ident = str(cl.get("identity", 0.99))
    ticker_label = f"vsearch clustering  (identity={ident})"

    try:
        with BlastTicker(ticker_label, output_file=uc_path,
                         total_queries=total_queries, avg_hsp=1.0, burn_in=30.0):
            run(
                [vsearch_path, "--cluster_fast", input_fasta, "--id", ident,
                 "--centroids", cen_path, "--uc", uc_path],
                outprefix + ".vsearch.log",
            )
        parse_uc(uc_path, ann)
    except RuntimeError as e:
        warn(f"Clustering failed: {e}")


# ---------------------------------------------------------------------------
# Chimera detection via vsearch uchime
# ---------------------------------------------------------------------------

def parse_uchimeout(
    path: str,
    ann: Dict[str, Annotation],
    reject_chimeras: bool = True,
    review_borderline: bool = True,
) -> tuple:
    """Parse a vsearch --uchimeout file and annotate chimera status in *ann*.

    vsearch uchimeout columns (tab-separated, 18 fields):
      0: score  1: query_label  2: parentA  3: parentB  4: top_parent
      5-16: identity/overlap stats  17: verdict (Y=chimera, N=clean, ?=borderline)

    Returns (chimera_count, borderline_count).
    """
    chimeras = 0
    borderline = 0
    try:
        with open(path, encoding="utf-8", errors="replace") as f:
            for line in f:
                parts = line.rstrip().split("\t")
                if len(parts) < 18:
                    continue
                qid = parts[1]
                verdict = parts[17]
                if qid not in ann:
                    continue
                if verdict == "Y":
                    ann[qid].add_reason("chimera_detected")
                    chimeras += 1
                elif verdict == "?" and review_borderline:
                    ann[qid].add_reason("chimera_borderline")
                    borderline += 1
    except FileNotFoundError:
        warn(f"Uchime output file not found: {path}")
    return chimeras, borderline


def run_uchime(
    input_fasta: str,
    outprefix: str,
    cfg: dict,
    ann: Dict[str, Annotation],
    total_queries: int = 0,
) -> tuple:
    """Run vsearch uchime chimera detection and annotate *ann*.

    Supports two modes controlled by ``chimera_screen.method``:
    - ``uchime_denovo``: de novo detection within the input set (no reference DB needed).
    - ``uchime_ref``: reference-based detection using ``chimera_screen.reference_db``.

    Returns (chimera_count, borderline_count).
    """
    cs = cfg.get("chimera_screen", {})
    vsearch_path = cs.get("vsearch_path", "vsearch")
    method = cs.get("method", "uchime_denovo")
    ref_db = cs.get("reference_db", "")
    reject_chimeras = cs.get("reject_chimeras", True)
    review_borderline = cs.get("review_borderline", True)
    abskew = str(cs.get("abskew", 2.0))

    if not shutil.which(vsearch_path):
        msg = f"vsearch not found at {vsearch_path!r}; chimera screen skipped."
        if cfg.get("run", {}).get("fail_on_missing_external_tool", False):
            raise RuntimeError(msg)
        warn(msg)
        return 0, 0

    uchimeout_path = outprefix + ".uchimeout.tsv"
    chimeras_path = outprefix + ".chimeras.fasta"
    borderline_path = outprefix + ".chimeras_borderline.fasta"
    nonchimeras_path = outprefix + ".nonchimeras.fasta"

    if method == "uchime_ref" and ref_db:
        info(f"  Method: uchime_ref  |  Reference DB: {ref_db}")
        cmd = [
            vsearch_path, "--uchime_ref", input_fasta,
            "--db", ref_db,
            "--abskew", abskew,
            "--uchimeout", uchimeout_path,
            "--chimeras", chimeras_path,
            "--borderline", borderline_path,
            "--nonchimeras", nonchimeras_path,
        ]
    else:
        if method == "uchime_ref" and not ref_db:
            warn("chimera_screen.method=uchime_ref but reference_db not set — falling back to uchime_denovo")
        info("  Method: uchime_denovo")
        cmd = [
            vsearch_path, "--uchime_denovo", input_fasta,
            "--abskew", abskew,
            "--uchimeout", uchimeout_path,
            "--chimeras", chimeras_path,
            "--borderline", borderline_path,
            "--nonchimeras", nonchimeras_path,
        ]

    ticker_label = f"vsearch uchime  ({method})"
    try:
        with BlastTicker(ticker_label, output_file=uchimeout_path,
                         total_queries=total_queries, avg_hsp=1.0, burn_in=30.0):
            run(cmd, outprefix + ".uchime.log")
        return parse_uchimeout(uchimeout_path, ann, reject_chimeras, review_borderline)
    except RuntimeError as e:
        warn(f"Chimera screen failed: {e}")
        return 0, 0
