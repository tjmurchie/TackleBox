"""FASTA parsing and writing.

Duplicate accessions are preserved as separate records. The second and later
copies receive a record_key of ``ACCESSION__dup2``, ``ACCESSION__dup3``, etc.
so that every input record has a unique internal identifier in decisions.tsv.
"""
from __future__ import annotations

import gzip
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import List

from .seq_utils import IUPAC_DNA  # noqa: F401 – re-exported for callers


@dataclass
class FastaRecord:
    id: str          # internal record_key (unique within a run)
    header: str      # original FASTA header line including leading '>'
    seq: str         # raw sequence (may be mixed case, multi-line joined)
    source_file: str = ""

    @property
    def accession(self) -> str:
        h = self.header[1:] if self.header.startswith(">") else self.header
        return h.split()[0]

    @property
    def seq_upper(self) -> str:
        return re.sub(r"\s+", "", self.seq.upper())


def parse_fasta(paths) -> List[FastaRecord]:
    """Parse one or more FASTA files, returning one FastaRecord per sequence.

    Duplicate accessions are handled without silently collapsing rows: the
    second occurrence of accession X gets record_key ``X__dup2``, the third
    ``X__dup3``, and so on.
    """
    recs: List[FastaRecord] = []
    if isinstance(paths, (str, Path)):
        paths = [paths]
    for p in paths:
        header = None
        chunks: List[str] = []
        opener = gzip.open if str(p).endswith(".gz") else open
        with opener(p, "rt", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                line = line.rstrip("\r\n")
                if not line:
                    continue
                if line.startswith(">"):
                    if header is not None:
                        recs.append(
                            FastaRecord(
                                id=header[1:].split()[0],
                                header=header,
                                seq="".join(chunks),
                                source_file=str(p),
                            )
                        )
                    header = line
                    chunks = []
                else:
                    chunks.append(line.strip())
            if header is not None:
                recs.append(
                    FastaRecord(
                        id=header[1:].split()[0],
                        header=header,
                        seq="".join(chunks),
                        source_file=str(p),
                    )
                )

    # Assign stable unique record_keys for duplicate accessions.
    seen: Counter = Counter()
    for r in recs:
        acc = r.accession
        seen[acc] += 1
        r.id = acc if seen[acc] == 1 else f"{acc}__dup{seen[acc]}"
    return recs


def write_fasta(recs: List[FastaRecord], path: str, wrap: int = 80) -> None:
    """Write records to *path* using original headers and upper-cased sequences."""
    with open(path, "w", encoding="utf-8") as out:
        for r in recs:
            out.write(r.header if r.header.startswith(">") else ">" + r.header)
            out.write("\n")
            s = r.seq_upper
            for i in range(0, len(s), wrap):
                out.write(s[i : i + wrap] + "\n")


def write_keyed_fasta(recs: List[FastaRecord], path: str, wrap: int = 80) -> None:
    """Write FASTA using Spinner record IDs as sequence identifiers.

    This avoids ambiguity when duplicate accessions are present. The original
    header text is retained after a space so the file remains human-readable.
    Used as query input for BLAST and clustering so that result parsing maps
    correctly back to the internal record_key.
    """
    with open(path, "w", encoding="utf-8") as out:
        for r in recs:
            safe_id = re.sub(r"[^A-Za-z0-9_.:-]", "_", r.id)
            original = r.header[1:] if r.header.startswith(">") else r.header
            out.write(f">{safe_id} original={original}\n")
            s = r.seq_upper
            for i in range(0, len(s), wrap):
                out.write(s[i : i + wrap] + "\n")
