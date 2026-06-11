"""
Shared TUI helpers for FlyGuide palaeo exporters.

Matches the visual style of NCBI-NT_Downloader.pl:
  ─── separator lines
  ═══ section headers
  ✓   completed items
  ▶   current item
  █░  progress bar
  In-place redraw on TTY; one-line fallback for pipes/logs.
"""
from __future__ import annotations

import os
import sys
import time
from typing import List, Optional, Tuple

RECENT_MAX = 5
_BURNIN = 3  # updates before ETA is shown


def _term_width() -> int:
    try:
        return max(os.get_terminal_size().columns, 60)
    except OSError:
        return 80


def _fmt_time(seconds: float) -> str:
    s = int(seconds)
    h, rem = divmod(s, 3600)
    m, sec = divmod(rem, 60)
    if h > 0:
        return f"{h}h {m:02d}m {sec:02d}s"
    if m > 0:
        return f"{m}m {sec:02d}s"
    return f"{sec}s"


# ── Static header / footer helpers ────────────────────────────────────────────

def print_header(title: str, params: List[Tuple[str, str]], file=sys.stderr) -> None:
    w = _term_width()
    sep = "═" * w
    print(sep, file=file)
    print(f"  {title}", file=file)
    print(sep, file=file)
    if params:
        kw = max(len(k) for k, _ in params)
        for k, v in params:
            print(f"  {k:<{kw}} : {v}", file=file)
        print("─" * w, file=file)


def print_step(n: int, total: int, label: str, file=sys.stderr) -> None:
    w = _term_width()
    print(f"\n>>> Step {n}/{total} — {label}", file=file)
    print("─" * w, file=file)


def print_done(label: str, counts: List[Tuple[str, str]], file=sys.stderr) -> None:
    w = _term_width()
    print("─" * w, file=file)
    print(f"  ✓ {label}", file=file)
    for k, v in counts:
        print(f"    {k:<32} {v}", file=file)


# ── Live progress dashboard ────────────────────────────────────────────────────

class PalaeoProgress:
    """
    In-place dashboard matching the Perl NCBI-NT_Downloader TUI.
    Falls back to plain one-line-per-update when stderr is not a TTY
    or when verbose=False.
    """

    def __init__(self, total: int, label: str = "taxa", verbose: bool = True):
        self.total = total
        self.label = label
        self.verbose = verbose
        self.completed = 0
        self._recent: List[Tuple[str, str]] = []  # (name, right_label)
        self._current_name = ""
        self._current_status = ""
        self._start = time.time()
        self._updates = 0
        self._is_tty = verbose and sys.stderr.isatty()
        self._lines_drawn = 0

    def set_current(self, name: str, status: str = "processing", completed: int = -1) -> None:
        self._current_name = name
        self._current_status = status
        if completed >= 0:
            self.completed = completed
        self._draw()

    def update(self, name: str, right_label: str = "") -> None:
        self._recent.append((name, right_label))
        if len(self._recent) > RECENT_MAX:
            self._recent.pop(0)
        self.completed += 1
        self._updates += 1
        self._current_name = ""
        self._current_status = ""
        self._draw()

    def finish(self, message: str = "") -> None:
        self.completed = self.total
        self._current_name = ""
        self._draw(final=True)
        if message and self.verbose:
            print(f"  {message}", file=sys.stderr)

    def _draw(self, final: bool = False) -> None:
        if not self.verbose:
            return
        w = _term_width()
        sep = "─" * w
        elapsed = time.time() - self._start
        frac = self.completed / self.total if self.total > 0 else 0
        elapsed_str = _fmt_time(elapsed)
        if frac > 0 and self._updates >= _BURNIN:
            eta_str = _fmt_time((elapsed / frac) - elapsed)
        else:
            eta_str = "calculating..."

        pct_label = f" {frac * 100:5.1f}%"
        bar_width = max(w - len(pct_label) - 4, 20)
        filled = int(frac * bar_width)
        bar = "  [" + "█" * filled + "░" * (bar_width - filled) + "]" + pct_label

        lines: List[str] = [sep, f"  Recent {self.label}:"]

        display = [None] * max(0, RECENT_MAX - len(self._recent)) + self._recent
        for entry in display:
            if entry is None:
                lines.append("")
                continue
            name, right = entry
            mark = "  ✓ "
            name_w = max(w - len(mark) - len(right) - 4, 15)
            lines.append(f"{mark}{name[:name_w]:<{name_w}}  {right}")

        if self._current_name:
            mark = "  ▶ "
            right = f"[{self._current_status:<14}]"
            name_w = max(w - len(mark) - len(right) - 4, 15)
            lines.append(f"{mark}{self._current_name[:name_w]:<{name_w}}  {right}")
        else:
            lines.append("")

        remaining = self.total - self.completed
        lines += [
            sep,
            f"  Progress: {self.completed} / {self.total}    Remaining: {remaining}",
            bar,
            f"  Elapsed: {elapsed_str:<16} ETA: {eta_str}",
            sep,
        ]

        if self._is_tty:
            if self._lines_drawn:
                print(f"\033[{self._lines_drawn}A", end="", file=sys.stderr)
            for line in lines:
                print(f"\033[2K{line}", file=sys.stderr)
            self._lines_drawn = len(lines)
        else:
            print(
                f"  Progress: {self.completed}/{self.total} ({frac * 100:.1f}%)  "
                f"elapsed: {elapsed_str}  ETA: {eta_str}",
                file=sys.stderr,
            )
