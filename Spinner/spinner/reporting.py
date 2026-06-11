"""Summary reports: TSV, HTML, and split FASTA outputs.

The HTML report is fully static (no external dependencies) and includes:
  - Decision summary boxes
  - Decisions by marker class table
  - Top rejection/review reasons
  - Adapter hit summary
  - Bad keyword hit summary
  - Taxonomy status summary (if run)
  - Windowed BLAST summary (if run)
  - Output file list
"""
from __future__ import annotations

import csv
import datetime
import html as html_mod
from collections import Counter
from typing import Dict, List

from .annotation import Annotation
from .fasta import FastaRecord, write_fasta


# ---------------------------------------------------------------------------
# TSV summary
# ---------------------------------------------------------------------------

def _species_coverage(ann: Dict[str, Annotation]) -> List[tuple]:
    """Return per-species decision counts sorted by 0-KEEP species first, then total desc.

    Returns list of (species, kingdom, keep, review, reject, total, note) tuples.
    Only includes records with an identified species (non-empty, non-Unknown).
    """
    from collections import defaultdict
    sp_data: dict = defaultdict(lambda: {"kingdom": "Unknown", "KEEP": 0, "REVIEW": 0, "REJECT": 0})
    for a in ann.values():
        sp = a.species_guess or a.genus_guess or ""
        if not sp or sp.lower() == "unknown":
            continue
        sp_data[sp]["kingdom"] = a.kingdom if a.kingdom != "Unknown" else sp_data[sp]["kingdom"]
        sp_data[sp][a.decision] += 1

    rows = []
    for sp, d in sp_data.items():
        keep, review, reject = d["KEEP"], d["REVIEW"], d["REJECT"]
        total = keep + review + reject
        note = "NO KEEP" if keep == 0 else ""
        rows.append((sp, d["kingdom"], keep, review, reject, total, note))

    rows.sort(key=lambda r: (r[6] == "", -r[5]))  # 0-KEEP first, then by total desc
    return rows


def write_summary_tsv(ann: Dict[str, Annotation], outprefix: str) -> None:
    counts = Counter(a.decision for a in ann.values())
    reasons = Counter(r for a in ann.values() for r in a.reasons)
    by_class = Counter((a.decision, a.marker_class) for a in ann.values())
    adapter_hits = Counter(a.adapter_name for a in ann.values() if a.adapter_hit)
    kw_hits = Counter(a.bad_keyword for a in ann.values() if a.bad_keyword_hit)
    taxonomy_counts = Counter(a.taxonomy_status for a in ann.values())

    with open(outprefix + ".summary.tsv", "w", encoding="utf-8") as out:
        out.write("section\tkey\tvalue\n")
        for k in ("KEEP", "REVIEW", "REJECT"):
            out.write(f"decision\t{k}\t{counts.get(k, 0)}\n")
        for (d, c), n in sorted(by_class.items()):
            out.write(f"decision_by_class\t{d}:{c}\t{n}\n")
        for r, n in reasons.most_common():
            out.write(f"reason\t{r}\t{n}\n")
        for name, n in adapter_hits.most_common():
            out.write(f"adapter_hit\t{name}\t{n}\n")
        for kw, n in kw_hits.most_common():
            out.write(f"bad_keyword_hit\t{kw}\t{n}\n")
        for status, n in taxonomy_counts.most_common():
            out.write(f"taxonomy_status\t{status}\t{n}\n")
        # Per-species coverage
        for row in _species_coverage(ann):
            sp, kd, keep, review, reject, total, note = row
            out.write(f"species_coverage\t{sp}\t{keep}/{review}/{reject}"
                      f"{'  NO_KEEP' if note else ''}\n")


# ---------------------------------------------------------------------------
# HTML summary
# ---------------------------------------------------------------------------

_CSS = """
body{font-family:monospace,Courier;max-width:1100px;margin:0 auto;padding:24px;color:#222}
h1{border-bottom:3px solid #333;padding-bottom:6px}
h2{border-bottom:1px solid #ccc;color:#555;margin-top:2em}
table{border-collapse:collapse;width:100%;margin-bottom:1.2em}
th,td{border:1px solid #ddd;padding:7px 10px;text-align:left}
th{background:#f0f0f0;font-weight:bold}
tr:nth-child(even){background:#fafafa}
tr.no-keep{background:#fff3cd}
tr.no-keep td{color:#856404}
.stat-wrap{display:flex;gap:16px;flex-wrap:wrap;margin-bottom:1.5em}
.stat-box{flex:1;min-width:150px;border:1px solid #ddd;border-radius:6px;
          padding:18px;text-align:center;background:#fafafa}
.stat-n{font-size:2.4em;font-weight:bold;line-height:1}
.KEEP .stat-n{color:#2a7}
.REVIEW .stat-n{color:#b70}
.REJECT .stat-n{color:#c33}
.KEEP .stat-label{color:#2a7}
.REVIEW .stat-label{color:#b70}
.REJECT .stat-label{color:#c33}
.stat-label{font-size:1.1em;margin-top:6px}
.none{color:#999;font-style:italic}
"""


def _table(headers: list, rows: list, row_class_fn=None) -> str:
    if not rows:
        return '<p class="none">None</p>'
    th = "".join(f"<th>{html_mod.escape(str(h))}</th>" for h in headers)
    body = ""
    for row in rows:
        cls = f' class="{row_class_fn(row)}"' if row_class_fn and row_class_fn(row) else ""
        body += f"<tr{cls}>" + "".join(f"<td>{html_mod.escape(str(c))}</td>" for c in row) + "</tr>"
    return f"<table><tr>{th}</tr>{body}</table>"


def write_summary_html(ann: Dict[str, Annotation], outprefix: str, version: str = "") -> None:
    from spinner import VERSION
    ver = version or VERSION
    counts = Counter(a.decision for a in ann.values())
    reasons = Counter(r for a in ann.values() for r in a.reasons)
    by_class: Counter = Counter()
    for a in ann.values():
        by_class[(a.decision, a.marker_class)] += 1

    adapter_hits = Counter(a.adapter_name for a in ann.values() if a.adapter_hit)
    adapter_rows = [(name, n) for name, n in adapter_hits.most_common()]

    kw_hits: Counter = Counter()
    for a in ann.values():
        if a.bad_keyword_hit:
            kw_hits[(a.bad_keyword, a.bad_keyword_action)] += 1
    kw_rows = [(kw, act, n) for (kw, act), n in kw_hits.most_common()]

    taxonomy_counts = Counter(a.taxonomy_status for a in ann.values())
    tax_rows = [(st, n) for st, n in taxonomy_counts.most_common()]
    tax_run = any(st != "NOT_CHECKED" for st in taxonomy_counts)

    windowed_counts = Counter(a.windowed_status for a in ann.values())
    wind_run = any(st != "NOT_CHECKED" for st in windowed_counts)
    wind_rows = [(st, n) for st, n in windowed_counts.most_common()]

    # Per-species coverage table (species, kingdom, KEEP, REVIEW, REJECT, total, note)
    sp_rows_raw = _species_coverage(ann)
    sp_rows = [(sp, kd, k, rv, rj, tot, note) for sp, kd, k, rv, rj, tot, note in sp_rows_raw]
    no_keep_count = sum(1 for r in sp_rows if r[6] == "NO KEEP")

    n_total = len(ann)
    ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # Class × decision table.
    classes = ["Mito", "Plastid", "NucMark", "Other"]
    decisions = ["KEEP", "REVIEW", "REJECT"]
    class_rows = [
        [klass] + [by_class.get((d, klass), 0) for d in decisions]
        for klass in classes
        if any(by_class.get((d, klass), 0) for d in decisions)
    ]

    # Top reasons table.
    reason_rows = [(r, n) for r, n in reasons.most_common(30)]

    # Output file list.
    out_files = [
        outprefix + ".decisions.tsv",
        outprefix + ".summary.tsv",
        outprefix + ".keep.fasta",
        outprefix + ".review.fasta",
        outprefix + ".reject.fasta",
    ]

    def stat_box(decision: str) -> str:
        n = counts.get(decision, 0)
        return (
            f'<div class="stat-box {decision}">'
            f'<div class="stat-n">{n}</div>'
            f'<div class="stat-label">{decision}</div>'
            f"</div>"
        )

    lines: List[str] = [
        "<!DOCTYPE html>",
        "<html lang='en'>",
        f"<head><meta charset='utf-8'><title>Spinner {ver} Summary</title>",
        f"<style>{_CSS}</style></head>",
        "<body>",
        f"<h1>Spinner {ver} Summary Report</h1>",
        f"<p>Generated: {ts} &nbsp;|&nbsp; Total records: <strong>{n_total}</strong></p>",
        "<h2>Decision summary</h2>",
        '<div class="stat-wrap">',
        stat_box("KEEP"),
        stat_box("REVIEW"),
        stat_box("REJECT"),
        "</div>",
        "<h2>Decisions by marker class</h2>",
        _table(["Marker class"] + decisions, class_rows),
        "<h2>Top reasons (keep / review / reject)</h2>",
        _table(["Reason", "Count"], reason_rows),
        "<h2>Adapter hits</h2>",
        _table(["Adapter name", "Records hit"], adapter_rows) if adapter_hits else '<p class="none">No adapter hits detected.</p>',
        "<h2>Bad keyword hits</h2>",
        _table(["Keyword", "Action", "Count"], kw_rows) if kw_hits else '<p class="none">No bad keyword hits.</p>',
        "<h2>Taxonomy BLAST status</h2>",
        _table(["Status", "Count"], tax_rows) if tax_run else '<p class="none">Taxonomy BLAST was not run.</p>',
        "<h2>Windowed BLAST status</h2>",
        _table(["Status", "Count"], wind_rows) if wind_run else '<p class="none">Windowed BLAST was not run.</p>',
        "<h2>Per-species coverage</h2>",
        (f'<p><strong>{no_keep_count} species with zero KEEP records</strong> (highlighted in amber) '
         f'— check these for sole-representative rescue or manual review.</p>'
         if no_keep_count else ""),
        _table(
            ["Species", "Kingdom", "KEEP", "REVIEW", "REJECT", "Total", "Note"],
            sp_rows,
            row_class_fn=lambda r: "no-keep" if r[6] == "NO KEEP" else "",
        ) if sp_rows else '<p class="none">No identified species in this dataset.</p>',
        "<h2>Output files</h2>",
        "<ul>" + "".join(f"<li><code>{html_mod.escape(f)}</code></li>" for f in out_files) + "</ul>",
        "</body></html>",
    ]

    with open(outprefix + ".summary.html", "w", encoding="utf-8") as out:
        out.write("\n".join(lines))


# ---------------------------------------------------------------------------
# Split FASTA output
# ---------------------------------------------------------------------------

def write_split_fastas(
    records: List[FastaRecord],
    ann: Dict[str, Annotation],
    outprefix: str,
    cfg: dict,
) -> None:
    """Write keep / review / reject FASTA files."""
    by_key = {r.id: r for r in records}
    tasks = [
        ("KEEP", "keep", True),
        ("REVIEW", "review", cfg.get("run", {}).get("write_review_fasta", True)),
        ("REJECT", "reject", cfg.get("run", {}).get("write_reject_fasta", True)),
    ]
    for decision, suffix, enabled in tasks:
        if not enabled:
            continue
        recs = [
            by_key[a.record_key]
            for a in ann.values()
            if a.decision == decision and a.record_key in by_key
        ]
        write_fasta(recs, f"{outprefix}.{suffix}.fasta")


# ---------------------------------------------------------------------------
# Report-from-decisions (for `Spinner report` subcommand)
# ---------------------------------------------------------------------------

def report_from_decisions(decisions_path: str, outprefix: str) -> None:
    """Re-generate summary TSV and HTML from an existing decisions TSV.

    Useful for reformatting or updating a report without re-running the pipeline.
    """
    ann: Dict[str, Annotation] = {}
    with open(decisions_path, encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            # Reconstruct a minimal Annotation from the TSV row.
            a = Annotation(
                accession=row.get("accession", ""),
                record_key=row.get("record_key", row.get("accession", "")),
                source_file=row.get("source_file", ""),
                header=row.get("header", ""),
                length=int(row.get("length") or 0),
                seq_sha256=row.get("seq_sha256", ""),
                species_guess=row.get("species_guess", ""),
                genus_guess=row.get("genus_guess", ""),
                kingdom=row.get("kingdom", "Unknown"),
                marker_class=row.get("marker_class", "Other"),
                region_id=row.get("region_id", ""),
                n_count=int(row.get("n_count") or 0),
                n_fraction=float(row.get("n_fraction") or 0),
                non_iupac_count=int(row.get("non_iupac_count") or 0),
                non_iupac_fraction=float(row.get("non_iupac_fraction") or 0),
                gc_fraction=float(row.get("gc_fraction") or 0),
                shannon_entropy=float(row.get("shannon_entropy") or 0),
                max_homopolymer=int(row.get("max_homopolymer") or 0),
                duplicate_accession=row.get("duplicate_accession", "").lower() == "true",
                duplicate_sequence=row.get("duplicate_sequence", "").lower() == "true",
                adapter_hit=row.get("adapter_hit", "").lower() == "true",
                adapter_name=row.get("adapter_name", ""),
                adapter_position=row.get("adapter_position", ""),
                adapter_internal=row.get("adapter_internal", "").lower() == "true",
                adapter_terminal=row.get("adapter_terminal", "").lower() == "true",
                vector_hit=row.get("vector_hit", "").lower() == "true",
                vector_internal=row.get("vector_internal", "").lower() == "true",
                vector_terminal=row.get("vector_terminal", "").lower() == "true",
                bad_keyword_hit=row.get("bad_keyword_hit", "").lower() == "true",
                bad_keyword=row.get("bad_keyword", ""),
                bad_keyword_action=row.get("bad_keyword_action", ""),
                taxonomy_status=row.get("taxonomy_status", "NOT_CHECKED"),
                taxonomy_top_hit=row.get("taxonomy_top_hit", ""),
                taxonomy_top_name=row.get("taxonomy_top_name", ""),
                taxonomy_top_pident=row.get("taxonomy_top_pident", ""),
                taxonomy_top_length=row.get("taxonomy_top_length", ""),
                taxonomy_top_staxids=row.get("taxonomy_top_staxids", ""),
                windowed_status=row.get("windowed_status", "NOT_CHECKED"),
                cluster_id=row.get("cluster_id", ""),
                cluster_role=row.get("cluster_role", ""),
                cap_rank=int(row.get("cap_rank") or 0),
                decision_score=int(row.get("decision_score") or 100),
                decision=row.get("decision", "KEEP"),
                reasons=[r for r in (row.get("reasons") or "").split(";") if r],
            )
            ann[a.record_key] = a

    write_summary_tsv(ann, outprefix)
    write_summary_html(ann, outprefix)
