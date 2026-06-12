"""AuthenticationPanel — per-taxon authentication scatter for Fillet inspector.

Shows damage vs reads scatter (one point per sample) for the selected taxon,
styled after metaDMG output.  No DB access required — all data comes from the
node dict returned by the tree/payload.
"""
from __future__ import annotations

import math
from typing import Dict, List, Optional

try:
    from qtpy.QtWidgets import (
        QWidget, QVBoxLayout, QHBoxLayout, QLabel, QFrame, QSizePolicy,
        QScrollArea, QGridLayout, QPushButton,
    )
    from qtpy.QtGui import QColor, QFont
    from qtpy.QtCore import Qt
except ImportError:
    raise

from .main_window import C_BG, C_BORDER, C_TEXT, C_DIM, C_GREEN, C_YELLOW, C_ORANGE, C_RED, C_BLUE

_TIER_COLORS: Dict[int, str] = {
    1: "#50b840", 2: "#3498e0", 3: "#c8a010",
    4: "#e08028", 5: "#7890a0", 0: "#e03838",
}
_TIER_LABELS: Dict[int, str] = {
    1: "T1 Excellent", 2: "T2 Strong", 3: "T3 Moderate",
    4: "T4 Weak", 5: "T5 Uncertain", 0: "T0 Rejected",
}


def _posterior_color(p: float) -> str:
    if p >= 0.90:
        return "#50b840"
    if p >= 0.70:
        return "#c8a010"
    if p >= 0.50:
        return "#e08028"
    return "#e03838"


def _kv(label: str, value: str, tip: str = "") -> QWidget:
    row = QWidget()
    lay = QHBoxLayout(row)
    lay.setContentsMargins(0, 0, 0, 0)
    lay.setSpacing(6)
    lbl = QLabel(label + ":")
    lbl.setStyleSheet(f"color: {C_DIM.name()}; font-size: 8pt;")
    lbl.setFixedWidth(120)
    if tip:
        lbl.setToolTip(tip)
    val = QLabel(value)
    val.setStyleSheet("color: #e6edf3; font-size: 8pt;")
    val.setTextInteractionFlags(Qt.TextSelectableByMouse)
    if tip:
        val.setToolTip(tip)
    lay.addWidget(lbl)
    lay.addWidget(val, 1)
    return row


class AuthenticationPanel(QWidget):
    """Per-taxon authentication scatter panel — 4th inspector tab."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self._node: Optional[dict] = None
        self._fig = None
        self._canvas = None
        self._font_scale: float = 1.0
        self._build_ui()

    def _build_ui(self) -> None:
        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.setSpacing(4)

        # Header
        hdr = QHBoxLayout()
        self._title_lbl = QLabel("Authentication summary")
        self._title_lbl.setStyleSheet(f"color: {C_DIM.name()}; font-style: italic; font-size: 9pt;")
        hdr.addWidget(self._title_lbl)
        self._tier_badge = QLabel("")
        self._tier_badge.setStyleSheet("font-weight: bold; font-size: 9pt; padding: 1px 8px;")
        self._tier_badge.hide()
        hdr.addWidget(self._tier_badge)
        hdr.addStretch()

        btn_font_up = QPushButton("A+")
        btn_font_up.setFixedSize(28, 22)
        btn_font_up.setToolTip("Increase font size")
        btn_font_up.clicked.connect(self._font_bigger)
        hdr.addWidget(btn_font_up)

        btn_font_dn = QPushButton("A-")
        btn_font_dn.setFixedSize(28, 22)
        btn_font_dn.setToolTip("Decrease font size")
        btn_font_dn.clicked.connect(self._font_smaller)
        hdr.addWidget(btn_font_dn)

        layout.addLayout(hdr)

        # Splitter: metrics on left, scatter on right
        from qtpy.QtWidgets import QSplitter
        splitter = QSplitter(Qt.Horizontal)
        splitter.setChildrenCollapsible(False)

        # Left: key metrics
        left = QScrollArea()
        left.setWidgetResizable(True)
        left.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        left.setStyleSheet("QScrollArea { border: none; background: transparent; }")
        left.setFixedWidth(260)
        self._metrics_widget = QWidget()
        self._metrics_layout = QVBoxLayout(self._metrics_widget)
        self._metrics_layout.setContentsMargins(2, 2, 2, 2)
        self._metrics_layout.setSpacing(2)
        self._metrics_layout.addStretch()
        left.setWidget(self._metrics_widget)

        # Right: matplotlib scatter
        self._canvas_container = QWidget()
        self._canvas_layout = QVBoxLayout(self._canvas_container)
        self._canvas_layout.setContentsMargins(0, 0, 0, 0)
        self._canvas_layout.setSpacing(0)
        placeholder = QLabel(
            "Select a taxon in the tree to view authentication scatter.\n"
            "Click the Load button to populate inspector panels."
        )
        placeholder.setAlignment(Qt.AlignCenter)
        placeholder.setStyleSheet(f"color: {C_DIM.name()}; font-size: 9pt;")
        placeholder.setWordWrap(True)
        self._placeholder = placeholder
        self._canvas_layout.addWidget(placeholder)

        splitter.addWidget(left)
        splitter.addWidget(self._canvas_container)
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        layout.addWidget(splitter, 1)

    def _init_canvas(self) -> None:
        if self._canvas is not None:
            return
        import matplotlib
        matplotlib.use("QtAgg")
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
        from qtpy.QtGui import QPalette as _QPalette
        self._placeholder.hide()
        self._fig = Figure(figsize=(5, 4), facecolor="#0d1220")
        self._canvas = FigureCanvas(self._fig)
        self._canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        pal = self._canvas.palette()
        pal.setColor(_QPalette.Window, QColor("#0d1220"))
        self._canvas.setPalette(pal)
        self._canvas.setAutoFillBackground(True)
        self._canvas.setStyleSheet("background: #0d1220;")
        self._canvas_layout.addWidget(self._canvas, 1)

    def _font_bigger(self) -> None:
        self._font_scale = min(2.0, self._font_scale * 1.2)
        self._redraw_scatter()

    def _font_smaller(self) -> None:
        self._font_scale = max(0.5, self._font_scale / 1.2)
        self._redraw_scatter()

    def _redraw_scatter(self) -> None:
        if self._node is not None and self._fig is not None:
            m = self._node.get("metrics") or {}
            self._draw_scatter(self._node, m)

    def clear(self) -> None:
        self._node = None
        self._title_lbl.setText("Authentication summary")
        self._title_lbl.setStyleSheet(f"color: {C_DIM.name()}; font-style: italic; font-size: 9pt;")
        self._tier_badge.hide()
        self._clear_metrics()
        if self._fig is not None:
            self._fig.clear()
            if self._canvas:
                self._canvas.draw()

    def _clear_metrics(self) -> None:
        while self._metrics_layout.count() > 0:
            item = self._metrics_layout.takeAt(0)
            w = item.widget()
            if w:
                w.deleteLater()

    def load_node(self, node: dict) -> None:
        """Populate authentication panel from a tree node dict (no DB required)."""
        self._node = node
        self._init_canvas()
        m = node.get("metrics") or {}
        name = node.get("name") or str(node.get("taxid", "?"))

        # Update header
        self._title_lbl.setText(name[:50])
        self._title_lbl.setStyleSheet("color: #e6edf3; font-size: 9pt; font-weight: bold;")
        tier = m.get("authenticity_tier")
        if tier:
            t = int(tier)
            self._tier_badge.setText(_TIER_LABELS.get(t, f"T{t}"))
            self._tier_badge.setStyleSheet(
                f"font-weight: bold; font-size: 9pt; padding: 1px 8px; color: {_TIER_COLORS.get(t, '#7890a0')};"
            )
            self._tier_badge.show()
        else:
            self._tier_badge.hide()

        # Build metrics panel
        self._clear_metrics()
        lo = self._metrics_layout

        def _hdr(text: str) -> None:
            lbl = QLabel(text)
            lbl.setStyleSheet("color: #7890a0; font-size: 7pt; font-weight: bold; "
                               "border-bottom: 1px solid #2c3a4a; margin-top: 4px;")
            lo.addWidget(lbl)

        def _row(label, val, tip="") -> None:
            lo.addWidget(_kv(label, str(val), tip))

        _hdr("Classification")
        _row("Taxon", name)
        _row("Taxid", str(node.get("taxid", "?")))
        rank = (node.get("lineage") or [{}])[-1].get("rank", "") if node.get("lineage") else ""
        if rank:
            _row("Rank", rank)

        _hdr("Authenticity")
        post = m.get("mean_posterior")
        if post is None:
            post = m.get("posterior_probability")
        _row("Tier", _TIER_LABELS.get(tier, "—") if tier else "—")
        _row("Posterior", f"{float(post):.3f}" if post is not None else "—",
             "Posterior probability that this assignment is correct (RELIC-LCA)")
        dmg = m.get("mean_damage_score")
        _row("CT damage", f"{float(dmg):.3f}" if dmg is not None else "—",
             "Mean 5' C→T misincorporation rate (proxy for aDNA damage)")
        stk = m.get("stack_concentration")
        if stk is None:
            stk = m.get("top_locus_fraction")
        _row("Stack conc.", f"{float(stk):.3f}" if stk is not None else "—",
             "Fraction of evidence reads at the single most covered position "
             "(high = reads stacked at one locus, possible false positive)")

        _hdr("Coverage")
        _row("Clade reads", f"{int(m.get('hard_reads') or 0):,}",
             "Reads assigned to this taxon plus all descendants (clade total)")
        breadth = m.get("best_reference_breadth")
        if breadth is None:
            breadth = m.get("breadth_of_coverage")
        _row("Breadth", f"{float(breadth):.4f}" if breadth is not None else "n/a",
             "Fraction of reference sequence covered by ≥1 read")

        _hdr("Blank / Decontam")
        bf = m.get("blank_fraction")
        _row("Blank fraction", f"{float(bf or 0):.3f}" if bf is not None else "—",
             "Fraction of reads in blank/negative controls (higher = potential contamination)")

        _hdr("Multi-Proxy Support")
        _row("Ecological",
             "✓" if m.get("eco_support") else "—",
             "Regional species list supports presence of this taxon (FlyGuide/Spinner)")
        _row("Palynological",
             "✓" if m.get("pal_support") else "—",
             "Pollen/spore record from same site supports this taxon")
        fos_flag = bool(m.get("fos_support"))
        fos_ev = m.get("fos_evidence_text") or ""
        fos_label = "🦴 ✓" if fos_flag else "—"
        fos_tip = (
            f"Macrofossil/archaeofaunal evidence: {fos_ev}"
            if fos_flag and fos_ev
            else "Physical macrofossil or archaeofaunal record supports this taxon at this site and time period"
            if fos_flag
            else "No fossil/archaeofaunal record for this taxon at this site in the relevant age range"
        )
        _row("Fossil record", fos_label, fos_tip)
        if fos_flag and fos_ev:
            ev_lbl = QLabel(fos_ev)
            ev_lbl.setStyleSheet(
                "color: #8bacc8; font-size: 7pt; padding-left: 8px; padding-bottom: 2px;"
            )
            ev_lbl.setWordWrap(True)
            lo.addWidget(ev_lbl)

        lo.addStretch()

        # Build scatter
        self._draw_scatter(node, m)

    def _draw_scatter(self, node: dict, m: dict) -> None:
        fs = self._font_scale  # local alias
        self._fig.clear()
        self._fig.set_facecolor("#0d1220")

        by_sample = m.get("by_sample") or {}
        if not by_sample:
            ax = self._fig.add_subplot(111)
            ax.set_facecolor("#131a28")
            ax.text(0.5, 0.5, "No per-sample data available.",
                    ha="center", va="center", transform=ax.transAxes,
                    color="#7890a0", fontsize=int(9 * fs))
            ax.set_axis_off()
            self._canvas.draw()
            return

        pts = []
        for sid, sm in by_sample.items():
            if not isinstance(sm, dict):
                continue
            r = int(sm.get("hard_reads") or 0)
            d = float(sm.get("mean_damage_score") or 0)
            p = float(sm.get("mean_posterior") or sm.get("posterior_probability") or 0)
            b = float(sm.get("best_reference_breadth") or sm.get("breadth_of_coverage") or 0)
            pts.append((sid, r, d, p, b))

        if not pts:
            ax = self._fig.add_subplot(111)
            ax.set_facecolor("#131a28")
            ax.text(0.5, 0.5, "No per-sample data available.",
                    ha="center", va="center", transform=ax.transAxes,
                    color="#7890a0", fontsize=int(9 * fs))
            ax.set_axis_off()
            self._canvas.draw()
            return

        ax1 = self._fig.add_subplot(121)
        ax2 = self._fig.add_subplot(122)
        _TEXT = "#e8f0fa"   # bright near-white for all text
        for ax in (ax1, ax2):
            ax.set_facecolor("#131a28")
            for sp in ax.spines.values():
                sp.set_color("#3a4a5e")
            ax.tick_params(colors=_TEXT, labelsize=int(7 * fs))
            ax.xaxis.label.set_color(_TEXT)
            ax.yaxis.label.set_color(_TEXT)
            ax.title.set_color(_TEXT)

        sids, reads, damages, posteriors, breadths = zip(*pts)

        # ── Left: Classification Evidence — posterior vs reads ──────────────────
        # Color = CT damage (0=blue → high=red); Size = breadth
        dmg_arr = list(damages)
        max_dmg = max(dmg_arr) if max(dmg_arr) > 0 else 0.10
        dmg_norm = [min(1.0, d / max(max_dmg, 0.01)) for d in dmg_arr]
        # Map damage to color: low=blue (#2980B9), mid=yellow (#c8a010), high=red (#C0392B)
        def _dmg_color(frac: float) -> str:
            if frac < 0.33:
                r = int(41 + (200 - 41) * frac / 0.33)
                g = int(128 + (160 - 128) * frac / 0.33)
                b = int(185 + (16 - 185) * frac / 0.33)
            elif frac < 0.66:
                f2 = (frac - 0.33) / 0.33
                r = int(200 + (192 - 200) * f2)
                g = int(160 + (57 - 160) * f2)
                b = int(16 + (43 - 16) * f2)
            else:
                f2 = (frac - 0.66) / 0.34
                r = int(192)
                g = int(57 * (1 - f2))
                b = int(43 * (1 - f2))
            return f"#{max(0,min(255,r)):02x}{max(0,min(255,g)):02x}{max(0,min(255,b)):02x}"

        colors1 = [_dmg_color(f) for f in dmg_norm]
        sizes1 = [max(25, min(220, 30 + float(b) * 250)) for b in breadths]
        xs1 = list(reads)
        ys1 = list(posteriors)
        ax1.scatter(xs1, ys1, c=colors1, s=sizes1, alpha=0.85,
                    edgecolors="#0d1220", linewidths=0.8, zorder=3)
        max_reads = max(reads, default=0)
        if max_reads > 0:
            ax1.set_xscale("log")
            import matplotlib.ticker as mticker
            ax1.xaxis.set_major_formatter(mticker.FuncFormatter(
                lambda x, _: f"{int(x):,}" if x >= 1 else f"{x:.1g}"
            ))
        else:
            ax1.set_xlim(-0.5, 1.5)
        ax1.set_ylim(-0.05, 1.10)
        ax1.set_xlabel("Reads", fontsize=int(7 * fs))
        ax1.set_ylabel("Posterior probability", fontsize=int(7 * fs))
        ax1.set_title("Classification Evidence", fontsize=int(8 * fs))
        # Threshold lines
        ax1.axhline(0.80, color="#50b840", linewidth=0.8, linestyle="--", alpha=0.6)
        ax1.axhline(0.50, color="#c8a010", linewidth=0.8, linestyle="--", alpha=0.6)
        ax1.grid(True, alpha=0.12, color="#3a4a5e")
        # Zone annotations
        ax1.text(0.02, 0.97, "Strong", transform=ax1.transAxes,
                 fontsize=int(6 * fs), color="#60d850", va="top", fontweight="bold")
        ax1.text(0.02, 0.55, "Moderate", transform=ax1.transAxes,
                 fontsize=int(6 * fs), color="#e0c040", va="top", fontweight="bold")
        ax1.text(0.02, 0.35, "Weak", transform=ax1.transAxes,
                 fontsize=int(6 * fs), color="#f05050", va="top", fontweight="bold")
        # Sample legend (instead of per-point annotation)
        for i, (sid, x, y) in enumerate(zip(sids, xs1, ys1)):
            ax1.plot([], [], 'o', color=colors1[i], markersize=4, label=sid[:15])
        _leg1 = ax1.legend(fontsize=int(5 * fs), loc='lower right', framealpha=0.4,
                           labelcolor=_TEXT, facecolor="#1c2a3a", edgecolor="#3a4a5e")
        if _leg1:
            _leg1.get_frame().set_linewidth(0.5)
        ax1.text(0.98, 0.02, "size ∝ breadth\ncolor = damage", transform=ax1.transAxes,
                 fontsize=int(5 * fs), color="#a0b8d0", ha="right", va="bottom")

        # ── Right: Authentication — damage vs breadth ───────────────────────────
        colors2 = [_posterior_color(p) for p in posteriors]
        sizes2 = [max(20, min(200, math.sqrt(r + 1) * 5)) for r in reads]
        xs2 = list(breadths)
        ys2 = list(damages)
        ax2.scatter(xs2, ys2, c=colors2, s=sizes2, alpha=0.80,
                    edgecolors="#0d1220", linewidths=0.5, zorder=3)
        ax2.set_xlabel("Breadth of coverage", fontsize=int(7 * fs))
        ax2.set_ylabel("CT damage", fontsize=int(7 * fs))
        ax2.set_title("Authentication", fontsize=int(8 * fs))
        ax2.set_xlim(-0.02, 1.05)
        # Threshold lines
        ax2.axhline(0.02, color="#c8a010", linewidth=0.8, linestyle=":", alpha=0.7)
        ax2.axhline(0.05, color="#e08028", linewidth=0.8, linestyle="--", alpha=0.6)
        ax2.axhline(0.15, color="#50b840", linewidth=0.8, linestyle="--", alpha=0.6)
        ax2.axvline(0.05, color="#7890a0", linewidth=0.8, linestyle=":", alpha=0.5)
        ax2.axvline(0.20, color="#7890a0", linewidth=0.8, linestyle="--", alpha=0.4)
        ax2.grid(True, alpha=0.12, color="#3a4a5e")
        # Zone annotations — breadth (x) close to 1 = good coverage; damage (y) > 0.05 = ancient
        # Top-right: high breadth + high damage = AUTHENTIC ✓
        ax2.text(0.98, 0.98, "Authentic", transform=ax2.transAxes,
                 fontsize=int(6 * fs), color="#60d850", va="top", ha="right", fontweight="bold")
        # Top-left: low breadth + high damage = damaged but poor coverage (uncertain)
        ax2.text(0.02, 0.98, "Damaged\npoor coverage", transform=ax2.transAxes,
                 fontsize=int(6 * fs), color="#e0c040", va="top", fontweight="bold")
        # Bottom-left: low breadth + low damage = not enough evidence
        ax2.text(0.02, 0.02, "Insufficient\nevidence", transform=ax2.transAxes,
                 fontsize=int(6 * fs), color="#f05050", va="bottom", fontweight="bold")
        # Bottom-right: high breadth + low damage = modern or contaminant
        ax2.text(0.98, 0.02, "Modern /\ncontaminant?", transform=ax2.transAxes,
                 fontsize=int(6 * fs), color="#e0c040", va="bottom", ha="right", fontweight="bold")
        ax2.text(0.50, 0.02, "size ∝ reads  ·  color = posterior", transform=ax2.transAxes,
                 fontsize=int(5 * fs), color="#a0b8d0", ha="center", va="bottom")

        # Bottom legend for authentication scatter — explicit white text for visibility
        import matplotlib.patches as mpatches
        patches = [
            mpatches.Patch(color="#60d850", label="Post ≥ 0.90"),
            mpatches.Patch(color="#e0c040", label="Post ≥ 0.70"),
            mpatches.Patch(color="#e09040", label="Post ≥ 0.50"),
            mpatches.Patch(color="#f05050", label="Post < 0.50"),
        ]
        _fig_leg = self._fig.legend(handles=patches, title="Posterior probability",
                                    fontsize=int(6 * fs), title_fontsize=int(6 * fs),
                                    loc="lower center", ncol=4,
                                    framealpha=0.5, facecolor="#1c2a3a", edgecolor="#3a4a5e",
                                    labelcolor=_TEXT, bbox_to_anchor=(0.5, 0.0))
        if _fig_leg:
            _fig_leg.get_title().set_color(_TEXT)
            _fig_leg.get_frame().set_linewidth(0.5)

        try:
            self._fig.tight_layout(pad=1.0, rect=[0, 0.08, 1, 1])
        except Exception:
            pass
        self._canvas.draw()
        self._canvas.flush_events()
        self._canvas.update()
