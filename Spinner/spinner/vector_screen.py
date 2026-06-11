"""Vector / UniVec BLAST screening.

parse_vector_blast() reads a blastn tabular output and annotates each record
with vector_hit / vector_internal / vector_terminal flags.

The expected outfmt is:
  6 qseqid sseqid pident length qstart qend sstart send evalue bitscore
"""
from __future__ import annotations

from typing import Dict

from .annotation import Annotation
from .adapters import is_terminal


def parse_vector_blast(path: str, ann: Dict[str, Annotation], cfg: dict) -> None:
    """Read the blastn vector screen TSV and annotate *ann* in-place."""
    vs = cfg.get("vector_screen", {})
    win = int(vs.get("terminal_window_bp", 25))
    min_len = int(vs.get("min_hit_length", 20))
    min_pid = float(vs.get("min_pident", 90.0))

    try:
        with open(path, encoding="utf-8", errors="replace") as f:
            for line in f:
                parts = line.rstrip().split("\t")
                if len(parts) < 6:
                    continue
                qid, _sseqid, pid_s, length_s, qstart_s, qend_s = parts[:6]
                # Record key may have been sanitised; find first match.
                if qid not in ann:
                    continue
                try:
                    pid = float(pid_s)
                    length = int(length_s)
                    qstart = int(qstart_s)
                    qend = int(qend_s)
                except ValueError:
                    continue
                if pid < min_pid or length < min_len:
                    continue
                # Convert 1-based BLAST coords to 0-based position.
                pos = min(qstart, qend) - 1
                a = ann[qid]
                a.vector_hit = True
                if is_terminal(pos, length, a.length, win):
                    a.vector_terminal = True
                    a.add_reason("vector_terminal")
                else:
                    a.vector_internal = True
                    a.add_reason("vector_internal")
    except FileNotFoundError:
        pass
