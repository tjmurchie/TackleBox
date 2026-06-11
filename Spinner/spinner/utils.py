"""Terminal helpers: logging, progress bars, timing."""
from __future__ import annotations

import os
import sys
import threading
import time
from typing import Optional


def eprint(*a, **kw):
    print(*a, file=sys.stderr, **kw)


def info(msg: str) -> None:
    eprint(f"[Spinner] {msg}")


def warn(msg: str) -> None:
    eprint(f"[Spinner:WARNING] {msg}")


def stage(msg: str) -> None:
    eprint("\n" + "=" * 78 + f"\n[Spinner] {msg}\n" + "=" * 78)


def section(msg: str) -> None:
    """Print a lightweight sub-section header (no double-rule)."""
    eprint(f"\n[Spinner] --- {msg} ---")


def fmt_seconds(s: float) -> str:
    s = int(max(0, s))
    h = s // 3600
    m = (s % 3600) // 60
    sec = s % 60
    return f"{h:02d}h:{m:02d}m:{sec:02d}s" if h else f"{m:02d}m:{sec:02d}s"


class BlastTicker:
    """Context manager that shows a progress bar on stderr while a blocking
    BLAST / MMSeqs2 subprocess runs.

    Because BLAST is a black box, progress is estimated by monitoring the
    output TSV file as BLAST writes results to it.  Each hit line written
    represents roughly one HSP; combined with ``total_queries`` and an
    assumed average HSP count (``avg_hsp``), an ETA is calculated after a
    short burn-in period.

    ETA is clearly marked as approximate — actual HSP counts per query vary
    with the database and sequence type.  The progress bar is a useful
    order-of-magnitude guide, not a precision clock.

    Parameters
    ----------
    label:
        Display label shown in the progress line.
    output_file:
        Path to the BLAST TSV output file.  Monitored for line count growth.
    total_queries:
        Total number of query sequences submitted to BLAST.
    avg_hsp:
        Assumed average output lines per query (including queries with no hits,
        which write 0 lines).  Default 5 is conservative for megablast vs NT
        with ``-max_target_seqs 10 -max_hsps 5``.  Use 2 for windowed BLAST
        (``-max_target_seqs 3 -max_hsps 1``).
    burn_in:
        Seconds before ETA is shown.  Used to accumulate enough data for a
        stable rate estimate.  Defaults to 60 s.
    interval:
        Refresh interval in seconds.  Defaults to 5 s.

    Usage::

        with BlastTicker("Taxonomy BLAST", output_file=out_t,
                         total_queries=73_589, avg_hsp=5):
            run_blast(query, db, out_t, cfg)
    """

    _FRAMES = ["|", "/", "-", "\\"]
    _BAR_WIDTH = 28

    def __init__(
        self,
        label: str,
        output_file: str = "",
        total_queries: int = 0,
        avg_hsp: float = 5.0,
        burn_in: float = 60.0,
        interval: float = 5.0,
    ) -> None:
        self.label = label
        self.output_file = output_file
        self.total_queries = total_queries
        self.avg_hsp = avg_hsp
        self.burn_in = burn_in
        self.interval = interval
        self.batch_info: Optional[list] = None  # [completed_batches, total_batches]
        self._stop: Optional[threading.Event] = None
        self._thread: Optional[threading.Thread] = None
        self._start: float = 0.0

    # ------------------------------------------------------------------
    def _count_lines(self) -> int:
        """Count newlines in the output file using fast binary reads."""
        if not self.output_file:
            return 0
        try:
            with open(self.output_file, "rb") as fh:
                return sum(
                    chunk.count(b"\n")
                    for chunk in iter(lambda: fh.read(1 << 20), b"")
                )
        except OSError:
            return 0

    def __enter__(self) -> "BlastTicker":
        self._start = time.time()
        self._stop = threading.Event()
        self._thread = threading.Thread(target=self._run, daemon=True)
        self._thread.start()
        return self

    def __exit__(self, *_) -> None:
        if self._stop:
            self._stop.set()
        if self._thread:
            self._thread.join()
        elapsed = time.time() - self._start
        bar = "#" * self._BAR_WIDTH
        sys.stderr.write(
            f"\r[Spinner]   \u2713  {self.label}\n"
            f"[Spinner]       [{bar}] 100.0%  [{fmt_seconds(elapsed)}]" + " " * 20 + "\n"
        )
        sys.stderr.flush()

    def _run(self) -> None:
        i = 0
        # burn-in snapshot: (elapsed, lines) recorded once burn-in completes
        burn_snap: Optional[tuple] = None

        while not self._stop.wait(self.interval):  # type: ignore[union-attr]
            elapsed = time.time() - self._start
            frame = self._FRAMES[i % len(self._FRAMES)]
            lines = self._count_lines()

            if self.total_queries and self.output_file:
                total_lines_est = self.total_queries * self.avg_hsp
                frac = min(1.0, lines / total_lines_est) if total_lines_est > 0 else 0.0
                pct = frac * 100.0
                done_blocks = int(self._BAR_WIDTH * frac)
                bar = "#" * done_blocks + "-" * (self._BAR_WIDTH - done_blocks)

                if self.batch_info and self.batch_info[1] > 0:
                    done, total = self.batch_info[0], self.batch_info[1]
                    batch_str = f"  batch {min(done + 1, total)}/{total}"
                else:
                    batch_str = ""

                if elapsed >= self.burn_in:
                    if burn_snap is None:
                        burn_snap = (elapsed, lines)
                    b_elapsed, b_lines = burn_snap
                    dt = elapsed - b_elapsed
                    if dt > 0 and lines > b_lines:
                        rate = (lines - b_lines) / dt        # lines/sec since burn-in
                        remaining = max(0.0, (total_lines_est - lines) / rate)
                        eta_str = fmt_seconds(remaining)
                    else:
                        eta_str = "calculating..."
                    status = (
                        f"{batch_str}"
                        f"  elapsed {fmt_seconds(elapsed)}"
                        f"  ETA ~{eta_str}  (estimate: ~{pct:.0f}% done)"
                    )
                else:
                    remaining_burn = self.burn_in - elapsed
                    status = (
                        f"{batch_str}"
                        f"  elapsed {fmt_seconds(elapsed)}"
                        f"  ({remaining_burn:.0f}s burn-in before ETA)"
                    )

                line = f"\r[Spinner]   {frame}  [{bar}] {pct:5.1f}%{status}" + " " * 5
            else:
                line = (
                    f"\r[Spinner]   {frame}  {self.label}"
                    f"  elapsed {fmt_seconds(elapsed)} ..."
                    + " " * 10
                )

            sys.stderr.write(line)
            sys.stderr.flush()
            i += 1


def progress(items, label: str = "Progress", enabled: bool = True):
    """Yield items from *items* with an optional terminal progress bar."""
    total = len(items) if hasattr(items, "__len__") else None
    start = time.time()
    last = 0.0
    for i, x in enumerate(items, 1):
        if enabled and total:
            now = time.time()
            if i == 1 or i == total or now - last > 0.25:
                frac = i / total
                width = 30
                done = int(width * frac)
                bar = "#" * done + "-" * (width - done)
                eta = (now - start) * (1 - frac) / frac if frac else 0
                sys.stderr.write(
                    f"\r{label} [{bar}] {i}/{total} ({frac * 100:5.1f}%)"
                    f" elapsed {fmt_seconds(now - start)} ETA {fmt_seconds(eta)}"
                )
                sys.stderr.flush()
                last = now
        yield x
    if enabled and total:
        elapsed = time.time() - start
        sys.stderr.write(f"\r{label} [{'#' * 30}] {total}/{total} (100.0%)"
                         f" elapsed {fmt_seconds(elapsed)} ETA 00m:00s\n")
        sys.stderr.flush()
