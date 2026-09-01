"""Clustering support: vsearch UC file parsing, cluster run, and uchime chimera detection."""
from __future__ import annotations

import os
import shutil
from collections import defaultdict
from typing import Dict, List

from .annotation import Annotation
from .external import run
from .fasta import FastaRecord, parse_fasta, write_fasta
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


def _remap_uc_cluster_ids(path: str, prefix: str) -> List[str]:
    """Read a per-group vsearch UC file, returning its lines with the cluster-number
    column (field 1) prefixed by *prefix* so cluster IDs stay globally unique once
    merged with every other group's own independently-numbered (0, 1, 2, ...) clusters.
    Returns ``[]`` (not an error) if *path* doesn't exist -- a group whose vsearch call
    failed simply contributes no cluster annotations, matching this module's existing
    degrade-rather-than-abort convention for a single failed clustering run."""
    lines: List[str] = []
    try:
        with open(path, encoding="utf-8", errors="replace") as f:
            for line in f:
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9:
                    continue
                parts[1] = f"{prefix}{parts[1]}"
                lines.append("\t".join(parts) + "\n")
    except FileNotFoundError:
        pass
    return lines


def run_vsearch(
    input_fasta: str,
    outprefix: str,
    cfg: dict,
    ann: Dict[str, Annotation],
    total_queries: int = 0,
) -> None:
    """Run vsearch clustering on *input_fasta* and annotate *ann*.

    Real bug, found 2026-09-01 via a live ~10,000-species PalaeoSCOPE Phase B run:
    `cluster.by` (default `["species_guess", "marker_class"]` in `config.py` -- not
    something a caller opted out of, every real project inherits this default) was
    declared in the config schema but never actually read anywhere in this function --
    every run clustered the ENTIRE input FASTA in ONE global `vsearch --cluster_fast`
    call regardless. At smoke-test scale (30 species, a few hundred candidate records)
    this was invisible; at real production scale (367,014 candidate records for one
    panel) it meant vsearch had to greedily compare every record against an
    ever-growing pool of centroids accumulated across ALL ~10,000 unrelated species
    combined, instead of the intended small, independent per-(species, marker_class)
    comparison pools -- observed live: 26+ CPU-hours and still 0 bytes of output after
    3.5+ wall-clock hours, PalaeoSCOPE's own progress estimator projecting a ~100 HOUR
    ETA. Clustering at 98%+ identity across genuinely different species was never going
    to produce a meaningful cross-species match anyway (real biological sequence
    divergence between species almost always exceeds a 2% identity threshold at these
    markers) -- so the global scope bought nothing semantically, only a combinatorial
    blowup no test at smoke-test scale could have caught.

    Fixed by honoring `cluster.by`: when set (the real default), records are grouped by
    the named `Annotation` attributes and each group gets its OWN independent
    `vsearch --cluster_fast` invocation -- typically tens of records each rather than
    hundreds of thousands, and vsearch's own greedy comparison pool never spans
    unrelated species. A singleton group (exactly one record) is trivially its own
    centroid and skips invoking vsearch entirely (no comparison possible with only one
    sequence) -- a common case at real scale (many species have very few candidate
    sequences for a given marker class) and a meaningful additional saving on its own.
    Per-group UC files are merged into the SAME `<outprefix>.clusters.uc`/
    `.cluster_centroids.fasta` paths this function has always written, with each
    group's own independently-numbered clusters remapped to a globally-unique ID
    (`g<index>_<original>`) so merging never collides two different groups' "cluster 0".
    Per-group work is resumable: a group whose own UC file already exists on disk from
    a prior partial run is not re-clustered. `cluster.by: []`/unset falls back to the
    exact previous single-global-call behavior, unchanged, for any config that
    deliberately wants the old behavior (or for the (unlikely) case where a project
    intentionally wants cross-species clustering).

    Groups are currently processed sequentially, not in a subprocess pool, matching
    this codebase's existing `run_mmseqs()` per-batch precedent -- the architectural
    fix (bounding each vsearch call's own input size) is what eliminates the real
    blowup; further parallelizing dispatch across groups is a possible future
    optimization if real-world timing after this fix is still not fast enough, not
    something this fix assumes is necessary.
    """
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
    by_fields = cl.get("by") or []

    if not by_fields:
        # Previous behavior, unchanged: one global clustering call across every record.
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
        return

    records = parse_fasta(input_fasta)
    groups: Dict[tuple, List[FastaRecord]] = defaultdict(list)
    for r in records:
        a = ann.get(r.id)
        groups[tuple(getattr(a, f, "") if a else "" for f in by_fields)].append(r)

    n_singletons = sum(1 for recs in groups.values() if len(recs) == 1)
    n_clustered_groups = len(groups) - n_singletons
    info(f"  Clustering by {by_fields}: {len(groups):,} groups "
         f"({n_singletons:,} singletons skip vsearch entirely, "
         f"{n_clustered_groups:,} groups need real clustering)")

    group_dir = outprefix + ".vsearch_groups"
    os.makedirs(group_dir, exist_ok=True)
    log_path = outprefix + ".vsearch.log"

    # Written incrementally (not buffered in memory then written once at the end) for
    # two real reasons: BlastTicker estimates progress by watching uc_path grow, so a
    # write-once-at-the-end would show 0% for the entire run and then jump to 100%; and
    # holding every record's full sequence in memory at once (hundreds of thousands of
    # records, some full organelle genomes) is real, avoidable memory pressure.
    ticker_label = f"vsearch clustering  (identity={ident}, {n_clustered_groups:,} groups)"
    with BlastTicker(ticker_label, output_file=uc_path, total_queries=len(records),
                     avg_hsp=1.0, burn_in=30.0), \
         open(uc_path, "w", encoding="utf-8") as uc_out, \
         open(cen_path, "w", encoding="utf-8") as cen_out:
        for gi, recs in enumerate(groups.values()):
            prefix = f"g{gi}_"
            if len(recs) == 1:
                r = recs[0]
                uc_out.write(f"S\t{prefix}0\t*\t*\t*\t*\t*\t*\t{r.id}\t*\n")
                cen_out.write(r.header if r.header.startswith(">") else ">" + r.header)
                cen_out.write("\n")
                s = r.seq_upper
                for i in range(0, len(s), 80):
                    cen_out.write(s[i:i + 80] + "\n")
                continue

            group_uc = os.path.join(group_dir, f"group_{gi:06d}.uc")
            group_cen = os.path.join(group_dir, f"group_{gi:06d}.centroids.fasta")
            if not (os.path.exists(group_uc) and os.path.getsize(group_uc) > 0):
                group_fasta = os.path.join(group_dir, f"group_{gi:06d}.fasta")
                write_fasta(recs, group_fasta)
                try:
                    run(
                        [vsearch_path, "--cluster_fast", group_fasta, "--id", ident,
                         "--centroids", group_cen, "--uc", group_uc],
                        log_path, verbose=False,
                    )
                except RuntimeError as e:
                    warn(f"Clustering failed for group {gi} ({len(recs)} records): {e}")
                    continue
            uc_out.writelines(_remap_uc_cluster_ids(group_uc, prefix))
            uc_out.flush()
            if os.path.exists(group_cen):
                with open(group_cen, encoding="utf-8") as gf:
                    shutil.copyfileobj(gf, cen_out)
                cen_out.flush()

    parse_uc(uc_path, ann)


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
