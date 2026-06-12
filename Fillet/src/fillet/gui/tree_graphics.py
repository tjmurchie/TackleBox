"""TaxonGraphicsTree — MEGAN7-style taxonomy tree using QGraphicsView.

Implements a Reingold-Tilford-style layout where each leaf = 1 unit of
horizontal space, internal nodes are centred over their children, and Y
depth increases downward.  Supports zoom (Ctrl+scroll), pan (middle-drag),
node selection, expand/collapse on double-click, and a context menu for
read extraction.
"""
from __future__ import annotations

import math
from typing import Dict, FrozenSet, List, Optional, Set

from qtpy.QtWidgets import (
    QGraphicsView, QGraphicsScene, QGraphicsEllipseItem,
    QGraphicsLineItem, QGraphicsTextItem, QGraphicsItem,
    QGraphicsRectItem, QGraphicsPolygonItem, QGraphicsPathItem,
    QMenu, QAction, QApplication, QWidget, QToolTip,
)
from qtpy.QtGui import (
    QPainter, QColor, QBrush, QPen, QFont, QPolygonF, QTransform,
    QWheelEvent, QMouseEvent, QContextMenuEvent, QDesktopServices,
    QPainterPath, QFontMetrics,
)
from qtpy.QtCore import Qt, QPoint, QPointF, QRect, QRectF, QUrl, Signal, QTimeLine, QTimer

from .spider_widget import metrics_to_values as _spider_metrics_to_values
from .main_window import (
    C_BG, C_BORDER, C_TEXT, C_DIM, C_GREEN, C_YELLOW, C_ORANGE, C_RED,
    C_TREE_BG, C_TREE_TEXT, C_TREE_DIM,
    SAMPLE_PALETTE, _support_color, _collect_nodes, _max_reads_in_tree,
    RANK_ORDER,
)

# ─── Layout constants ──────────────────────────────────────────────────────────
_V_SPACING = 70          # px between depth levels
_H_UNIT = 90             # px per leaf unit
_LABEL_OFFSET = 6        # px right of node circle edge
_MIN_RADIUS = 4       # zero-reads / pass-through nodes
_MIN_POS_RADIUS = 7  # minimum for any node with reads > 0
_MAX_RADIUS = 30     # maximum; large enough to show clear size differences
_BRANCH_WIDTH = 1.5

# Taxonomic ranks considered "structurally important" for label visibility.
# Nodes at these ranks always get labels even when they have no direct reads,
# so the user can orient themselves in the tree.
_MAJOR_RANKS: frozenset = frozenset({
    "genus", "family", "order", "class", "phylum",
    "kingdom", "superkingdom", "domain",
})


def _node_radius(reads: int, max_reads: int, scale_mode: str = "sqrt") -> float:
    """Map read count to a node radius in scene pixels.

    Nodes with zero reads get _MIN_RADIUS (small dot).
    Nodes with reads > 0 map from _MIN_POS_RADIUS (fewest reads) to _MAX_RADIUS
    (most reads) using the chosen scale mode, giving clear visual differentiation.
    """
    if max_reads <= 0 or reads <= 0:
        return _MIN_RADIUS
    if scale_mode == "linear":
        ratio = reads / max_reads
    elif scale_mode == "log":
        ratio = math.log1p(reads) / math.log1p(max_reads)
    else:  # sqrt (default)
        ratio = math.sqrt(reads) / math.sqrt(max_reads)
    # Map ratio to [_MIN_POS_RADIUS, _MAX_RADIUS] so any node with reads
    # is visibly larger than zero-read pass-through nodes.
    r = _MIN_POS_RADIUS + (_MAX_RADIUS - _MIN_POS_RADIUS) * ratio
    return max(_MIN_POS_RADIUS, min(_MAX_RADIUS, r))


# ─── Internal tree layout node ────────────────────────────────────────────────

class _LayoutNode:
    """Wraps a raw taxon dict with layout information."""

    def __init__(self, data: dict, depth: int, parent: Optional["_LayoutNode"] = None):
        self.data = data
        self.depth = depth
        self.parent = parent
        self.children: List[_LayoutNode] = []
        self.x: float = 0.0   # set by layout pass (leaf units)
        self.y: float = 0.0   # set by layout pass (px)
        self.collapsed: bool = False
        self.reads: int = int((data.get("metrics") or {}).get("hard_reads") or 0)
        self.taxid: str = str(data.get("taxid") or "")
        self.name: str = data.get("name") or data.get("taxid") or "?"
        self.flags: list = data.get("flags") or []
        self.by_sample: Dict[str, int] = {}   # populated during load
        m = data.get("metrics") or {}
        # direct_hard_reads: reads assigned only to this node (not descendants).
        # Absent in older SQLite files (fallback 0 is safe — "direct" label shows nothing).
        self.direct_reads: int = int(m.get("direct_hard_reads") or 0)
        _top_tier = m.get("authenticity_tier")
        _top_badge = str(m.get("authenticity_badge") or "")
        if _top_tier is None:
            # Fall back to minimum tier across samples (same logic as spider widget)
            _by_s = m.get("by_sample") or {}
            _tier_vals = [
                int((sm or {}).get("authenticity_tier"))
                for sm in _by_s.values()
                if (sm or {}).get("authenticity_tier") is not None
            ]
            _top_tier = min(_tier_vals) if _tier_vals else None
            # Badge fallback: pick badge from the sample with the min tier
            if _top_tier is not None and not _top_badge:
                for _sm in _by_s.values():
                    _b = (_sm or {}).get("authenticity_badge") or ""
                    if _b:
                        _top_badge = str(_b)
                        break
        self.authenticity_tier: Optional[int] = (
            int(_top_tier) if _top_tier is not None else None
        )
        self.authenticity_badge: str = _top_badge
        # Pre-computed metric scalars for support filtering
        _dmg = float(m.get("mean_damage_score") or 0)
        if _dmg == 0:
            by_s = m.get("by_sample") or {}
            _vals = [float((sm or {}).get("mean_damage_score") or 0) for sm in by_s.values()]
            _dmg = sum(_vals) / len(_vals) if _vals else 0
        self.mean_damage: float = _dmg
        self.mean_posterior: float = float(m.get("mean_posterior") or 0)
        self.breadth: float = float(m.get("best_reference_breadth") or 0)
        # Fossil/archaeofaunal support flag (age/site-resolved)
        _fos = bool(m.get("fos_support"))
        if not _fos:
            _by_s = m.get("by_sample") or {}
            _fos = any(bool((sm or {}).get("fos_support")) for sm in _by_s.values())
        self.fos_support: bool = _fos
        # Virtual (synthetic) node flags — set for the unclassified/excluded branch
        self.is_virtual_node: bool = bool(data.get("is_virtual_node"))
        self.virtual_category: str = str(data.get("virtual_category") or "")
        self.virtual_tooltip: str = str(data.get("virtual_tooltip") or "")

    def visible_children(self) -> List["_LayoutNode"]:
        if self.collapsed:
            return []
        return self.children

    def all_reads(self) -> int:
        """Sum reads of self and descendants (for filter pruning)."""
        if not hasattr(self, "_cache_all_reads"):
            total = self.reads
            for c in self.children:
                total += c.all_reads()
            self._cache_all_reads: int = total
        return self._cache_all_reads

    def max_node_reads(self) -> int:
        """Maximum reads at any single node in this subtree (for radius scaling)."""
        if not hasattr(self, "_cache_max_reads"):
            m = self.reads
            for c in self.children:
                m = max(m, c.max_node_reads())
            self._cache_max_reads: int = m
        return self._cache_max_reads

    def max_direct_reads(self) -> int:
        """Maximum direct_reads at any single node in this subtree."""
        if not hasattr(self, "_cache_max_direct"):
            m = self.direct_reads
            for c in self.children:
                m = max(m, c.max_direct_reads())
            self._cache_max_direct: int = m
        return self._cache_max_direct


def _subtree_taxids(data: dict) -> Set[str]:
    """Collect all taxids in the data subtree (for focus filtering)."""
    result: Set[str] = set()
    tid = str(data.get("taxid") or "")
    if tid:
        result.add(tid)
    for child in data.get("children") or []:
        result |= _subtree_taxids(child)
    return result


def _build_layout_tree(
    data: dict,
    depth: int = 0,
    parent: Optional[_LayoutNode] = None,
    selected_samples: Optional[List[str]] = None,
    focus_taxids: Optional[Set[str]] = None,
    support_filter: str = "all",
) -> Optional[_LayoutNode]:
    """Recursively build _LayoutNode tree, pruning nodes with no reads.

    focus_taxids: if set, only keep nodes whose taxid is in the set OR that
    have at least one descendant in the set (ancestor path nodes are kept to
    maintain tree structure).
    support_filter: multi-metric threshold key from _SUPPORT_LEVELS; leaf
    nodes that fail the threshold are pruned (internal nodes kept when a
    descendant passes).
    """
    node = _LayoutNode(data, depth, parent)

    # Compute per-sample reads filtered to selected_samples
    by_sample_raw = (data.get("metrics") or {}).get("by_sample") or {}
    if selected_samples is not None:
        node.by_sample = {
            sid: int((m or {}).get("hard_reads") or 0)
            for sid, m in by_sample_raw.items()
            if sid in selected_samples
        }
        node.reads = sum(node.by_sample.values())
    else:
        node.by_sample = {
            sid: int((m or {}).get("hard_reads") or 0)
            for sid, m in by_sample_raw.items()
        }
        node.reads = int((data.get("metrics") or {}).get("hard_reads") or 0)

    # When focus_taxids is set: once we're inside a focused subtree, pass None
    # so all descendants are included; otherwise propagate the filter.
    this_tid = str(data.get("taxid") or "")
    child_focus = None if (focus_taxids is None or this_tid in (focus_taxids or set())) else focus_taxids

    for child_data in (data.get("children") or []):
        child_node = _build_layout_tree(
            child_data, depth + 1, node, selected_samples, child_focus, support_filter
        )
        if child_node is not None:
            node.children.append(child_node)

    # Focus filter: prune if this node is not in focus AND has no focused descendants
    if focus_taxids is not None and this_tid not in focus_taxids:
        if not node.children:
            return None   # leaf not in focus — prune

    # Read-based prune: keep node only if it or its subtree has reads
    if node.reads == 0 and not node.children:
        return None

    # Support filter: prune leaf nodes that fail the multi-metric threshold.
    # Internal nodes (with children) are kept as long as a descendant passed.
    if support_filter != "all" and not node.children:
        if not _node_passes_support(node, support_filter):
            return None

    # Infer direct_reads for leaf nodes in old SQLite files (pre-direct_hard_reads).
    # For a leaf with no children, every read is direct by definition.
    if not node.children and node.direct_reads == 0 and node.reads > 0:
        node.direct_reads = node.reads

    return node


# ─── Reingold-Tilford leaf-counting layout ────────────────────────────────────

def _count_all_leaves(node: _LayoutNode) -> int:
    """Count leaves in the FULL tree (ignoring collapse state)."""
    if not node.children:
        return 1
    return sum(_count_all_leaves(c) for c in node.children)


def _apply_auto_collapse(root: _LayoutNode, max_depth: int = 5) -> None:
    """For large trees (> 100 leaves), collapse all internal nodes at depth >= max_depth."""
    if _count_all_leaves(root) <= 100:
        return  # Small tree — show everything expanded

    def _collapse(node: _LayoutNode) -> None:
        if node.children:
            node.collapsed = node.depth >= max_depth
            if not node.collapsed:
                for c in node.children:
                    _collapse(c)
    _collapse(root)


def _smart_collapse_node(
    node: _LayoutNode,
    min_reads: int,
    large_tree: bool,
    depth_cap: int,
) -> bool:
    """Recursively decide collapse state.  Returns True if this subtree contains
    at least one well-supported (non-critical-flag, reads >= min_reads) node."""
    # Virtual nodes (unclassified branch) are always kept expanded
    if node.is_virtual_node:
        node.collapsed = False
        for child in node.children:
            _smart_collapse_node(child, min_reads, large_tree, depth_cap)
        return True
    # Hard depth cap — always collapse below this level
    if node.depth >= depth_cap:
        if node.children:
            node.collapsed = True
        return False

    _CRITICAL = {"low_read_support", "phylo_isolated"}
    node_is_good = (
        node.reads >= min_reads
        and not any(f in _CRITICAL for f in node.flags)
    )

    if not node.children:
        return node_is_good

    # Recurse into children first
    any_child_good = False
    for child in node.children:
        child_good = _smart_collapse_node(child, min_reads, large_tree, depth_cap)
        any_child_good = any_child_good or child_good

    # Collapse if neither this node nor any descendant is well-supported
    if not node_is_good and not any_child_good:
        node.collapsed = True
    else:
        node.collapsed = False

    return node_is_good or any_child_good


def _apply_smart_collapse(root: _LayoutNode, min_reads: int = 3) -> None:
    """Collapse subtrees that contain only low-quality / critically-flagged
    assignments.  Well-supported nodes (reads >= min_reads, no critical flags)
    are kept expanded so the user can explore their children.

    Critical flags that trigger collapse consideration:
        ``low_read_support``, ``phylo_isolated``

    Non-critical (kept expanded regardless):
        ``conflicted_reads_present``, ``incongruent_with_dominant``,
        ``high_entropy_family``, etc.
    """
    too_large = _count_all_leaves(root) > 200
    depth_cap = 6 if too_large else 8
    _smart_collapse_node(root, min_reads, too_large, depth_cap)


def _count_leaves(node: _LayoutNode) -> int:
    vis = node.visible_children()
    if not vis:
        return 1
    return sum(_count_leaves(c) for c in vis)


def _assign_x(node: _LayoutNode, offset: float = 0.0) -> float:
    """Assign x positions (in leaf units). Returns total width consumed."""
    vis = node.visible_children()
    if not vis:
        node.x = offset + 0.5
        return 1.0
    x = offset
    for c in vis:
        w = _assign_x(c, x)
        x += w
    node.x = (vis[0].x + vis[-1].x) / 2.0
    return x - offset


def _assign_y(node: _LayoutNode, v_spacing: float = _V_SPACING) -> None:
    node.y = float(node.depth * v_spacing)
    for c in node.children:   # recurse all children (even collapsed, we need their coords)
        _assign_y(c, v_spacing)


def _max_depth(node: _LayoutNode) -> int:
    """Return max visible depth in this subtree."""
    vis = node.visible_children()
    if not vis:
        return node.depth
    return max(_max_depth(c) for c in vis)


def _assign_y_dendrogram(node: _LayoutNode, max_depth: int,
                          v_spacing: float = _V_SPACING) -> None:
    """Right-align leaves: leaf Y = max_depth * v_spacing, internal Y pushed up."""
    vis = node.visible_children()
    if not vis:
        node.y = float(max_depth * v_spacing)
    else:
        for c in node.children:
            _assign_y_dendrogram(c, max_depth, v_spacing)
        child_ys = [c.y for c in vis]
        node.y = min(child_ys) - v_spacing


# ─── Multi-metric support levels ──────────────────────────────────────────────
# Each preset defines: (min_reads, min_damage, min_posterior, max_tier, flagged_only)
# A leaf node is pruned from the tree when it fails ALL of the criteria.
# Internal (pass-through) nodes are kept as long as a descendant passes.

_SUPPORT_LEVELS: Dict[str, dict] = {
    "all":      {"reads": 0,   "damage": 0.00, "posterior": 0.00, "tier": 99, "flagged": False},
    "minimal":  {"reads": 3,   "damage": 0.00, "posterior": 0.00, "tier": 99, "flagged": False},
    "low":      {"reads": 10,  "damage": 0.00, "posterior": 0.00, "tier": 99, "flagged": False},
    "moderate": {"reads": 10,  "damage": 0.02, "posterior": 0.70, "tier":  5, "flagged": False},
    "confident":{"reads": 20,  "damage": 0.05, "posterior": 0.80, "tier":  4, "flagged": False},
    "high":     {"reads": 50,  "damage": 0.10, "posterior": 0.85, "tier":  3, "flagged": False},
    "flagged":  {"reads": 0,   "damage": 0.00, "posterior": 0.00, "tier": 99, "flagged": True},
}


def _node_passes_support(node: "_LayoutNode", sf: str) -> bool:
    """Return True if the node meets the minimum support threshold for level sf."""
    if sf not in _SUPPORT_LEVELS:
        return True
    lvl = _SUPPORT_LEVELS[sf]
    if lvl["flagged"]:
        return len(node.flags) > 0
    if node.reads < lvl["reads"]:
        return False
    tier = node.authenticity_tier
    if tier is not None and tier > lvl["tier"]:
        return False
    # For damage/posterior: pass if EITHER metric meets threshold (OR logic),
    # or if the threshold is 0 (no constraint).
    dmg_ok = (lvl["damage"] == 0) or (node.mean_damage >= lvl["damage"])
    post_ok = (lvl["posterior"] == 0) or (node.mean_posterior >= lvl["posterior"])
    if lvl["damage"] > 0 or lvl["posterior"] > 0:
        return dmg_ok or post_ok
    return True


# ─── Graphics items ────────────────────────────────────────────────────────────

_AUTH_TIER_COLORS = {
    1: QColor("#50b840"),   # logo green  — highest authenticity
    2: QColor("#3498e0"),   # logo blue
    3: QColor("#c8a010"),   # logo gold
    4: QColor("#e08028"),   # logo orange
    5: QColor("#7890a0"),   # blue-grey   — minimal evidence
    0: QColor("#e03838"),   # logo red    — flagged
}


def _metric_color(
    node: "_LayoutNode",
    color_by: str,
    sample_colors: Optional[Dict[str, "QColor"]] = None,
) -> "QColor":
    """Return the fill color for a node given the current color-by mode."""
    m = node.data.get("metrics") or {}
    if color_by == "damage":
        dmg = float(m.get("mean_damage_score") or 0)
        if dmg == 0:
            by_s = m.get("by_sample") or {}
            vals = [float((sm or {}).get("mean_damage_score") or 0) for sm in by_s.values()]
            dmg = sum(vals) / len(vals) if vals else 0
        if node.reads == 0:
            return C_BORDER
        if dmg >= 0.25:
            return C_GREEN
        if dmg >= 0.10:
            return C_YELLOW
        if dmg >= 0.03:
            return C_ORANGE
        return C_RED
    if color_by == "tier":
        tier = node.authenticity_tier
        return _AUTH_TIER_COLORS.get(tier, C_DIM) if tier is not None else C_DIM
    if color_by == "blank":
        blank = float(m.get("blank_fraction") or 0)
        if node.reads == 0:
            return C_BORDER
        if blank >= 0.3:
            return C_RED
        if blank >= 0.1:
            return C_ORANGE
        if blank > 0:
            return C_YELLOW
        return C_GREEN
    if color_by == "posterior":
        post = float(m.get("mean_posterior") or 0)
        if node.reads == 0:
            return C_BORDER
        if post >= 0.9:
            return C_GREEN
        if post >= 0.7:
            return C_YELLOW
        if post >= 0.4:
            return C_ORANGE
        return C_RED
    if color_by == "sample" and sample_colors:
        # Color by dominant sample (the one contributing the most reads)
        by_s = node.by_sample
        if by_s:
            dominant = max(by_s, key=lambda sid: by_s.get(sid, 0))
            return sample_colors.get(dominant, _support_color(node.flags, node.reads))
    # default: "support"
    return _support_color(node.flags, node.reads)


_CMP_A_COLOR  = QColor("#3498e0")   # logo blue   — more reads in sample A
_CMP_B_COLOR  = QColor("#e08028")   # logo orange — more reads in sample B
_CMP_EQ_COLOR = QColor("#2c3a4a")   # navy-grey   — roughly equal


def _comparison_color(node: "_LayoutNode", sample_a: str, sample_b: str) -> "QColor":
    """Return A>B / B>A / equal color for comparison mode."""
    ra = node.by_sample.get(sample_a, 0)
    rb = node.by_sample.get(sample_b, 0)
    if ra == 0 and rb == 0:
        return C_BORDER
    if ra == rb:
        return _CMP_EQ_COLOR
    ratio = ra / (ra + rb) if (ra + rb) > 0 else 0.5
    # Blue for A-dominant, orange for B-dominant; intensity scales with imbalance
    if ratio > 0.5:
        strength = (ratio - 0.5) * 2.0   # 0..1
        return _CMP_A_COLOR.lighter(int(100 + 60 * (1 - strength)))
    else:
        strength = (0.5 - ratio) * 2.0
        return _CMP_B_COLOR.lighter(int(100 + 60 * (1 - strength)))


class _NodeItem(QGraphicsEllipseItem):
    """A clickable node circle (pie chart or coxcomb for multi-sample).

    Labels are attached as child QGraphicsTextItem items so they move with
    the node during animations and receive click-through to the parent node.
    """

    def __init__(self, layout_node: _LayoutNode, radius: float,
                 color: QColor, sample_colors: Dict[str, QColor],
                 review_mark: str = "", display_mode: str = "pie",
                 ltr_layout: bool = False, dark_taxon: bool = False,
                 color_by: str = "support"):
        super().__init__(-radius, -radius, radius * 2, radius * 2)
        self.layout_node = layout_node
        self._dark_taxon = dark_taxon
        self._color = _DARK_TAXON_COLOR if dark_taxon else color
        self._radius = radius
        self._sample_colors = sample_colors
        self._review_mark = review_mark
        self._display_mode = display_mode   # "pie" | "coxcomb" | "spider"
        self._auth_badge = layout_node.authenticity_badge
        self._auth_tier = layout_node.authenticity_tier
        self._fos_support = layout_node.fos_support
        self._ltr_layout = ltr_layout
        self._is_virtual_node: bool = layout_node.is_virtual_node
        if self._is_virtual_node:
            # Override color for virtual nodes — ignore the passed-in color
            self._color = _VIRTUAL_NODE_COLOR
            self._dark_taxon = False
        self.setFlag(QGraphicsItem.ItemIsSelectable, True)
        self.setAcceptHoverEvents(True)
        self.setCursor(Qt.PointingHandCursor)
        self._selected_override = False
        self._multi_selected = False
        self._search_match = False
        self._label_item: Optional[QGraphicsTextItem] = None

        # Pre-compute spider spoke values (7 spokes matching spider_widget)
        m = layout_node.data.get("metrics") or {}
        self._spider_vals: List[float] = _spider_metrics_to_values(m)

        sample_counts = [(sid, cnt) for sid, cnt in layout_node.by_sample.items() if cnt > 0]
        self._is_pie = (color_by == "sample") and len(sample_counts) > 1 and bool(sample_colors)
        self._pie_slices = sample_counts if self._is_pie else []

        if not self._is_pie:
            self.setBrush(QBrush(color))
            self.setPen(QPen(color.darker(140), 1.2))

    def set_label(self, text: str, font: QFont, color: QColor,
                  offset_x: float, offset_y: float) -> None:
        """Attach a text label as a child item at the given node-relative offset."""
        lbl = QGraphicsTextItem(text, self)
        lbl.setFont(font)
        lbl.setDefaultTextColor(color)
        lbl.setPos(offset_x, offset_y)
        lbl.setFlag(QGraphicsItem.ItemIsSelectable, False)
        lbl.setAcceptedMouseButtons(Qt.NoButton)  # pass clicks through to scene
        self._label_item = lbl

    def hoverEnterEvent(self, event) -> None:
        """Show tooltip immediately on hover — bypasses the OS delay."""
        tip = self.toolTip()
        if tip:
            sp = event.screenPos()
            # screenPos() returns QPoint in some PyQt5 builds, QPointF in others
            pos = sp.toPoint() if hasattr(sp, "toPoint") else sp
            QToolTip.showText(pos, tip)
        super().hoverEnterEvent(event)

    def hoverLeaveEvent(self, event) -> None:
        QToolTip.hideText()
        super().hoverLeaveEvent(event)

    def set_search_match(self, matched: bool) -> None:
        if self._search_match != matched:
            self._search_match = matched
            self.update()

    def set_multi_selected(self, state: bool) -> None:
        if self._multi_selected != state:
            self._multi_selected = state
            self.update()

    def _draw_spider_node(self, painter: "QPainter", r: float) -> None:
        """Draw a mini spider/radar web node (replaces pie/coxcomb when display_mode='spider')."""
        n = len(self._spider_vals)
        if n == 0 or r < 4:
            return
        tier   = self._auth_tier or 0
        colors = {1: QColor("#1a9850"), 2: QColor("#3498e0"), 3: QColor("#c8a010"),
                  4: QColor("#e08028"), 5: QColor("#e03838"), 0: QColor("#555555")}
        base_col  = colors.get(tier, self._color)
        grid_col  = QColor("#2c3a4a")
        fill_col  = QColor(base_col); fill_col.setAlpha(70)
        edge_col  = QColor(base_col); edge_col.setAlpha(220)

        # Web rings (0.33, 0.67, 1.0)
        painter.setBrush(Qt.NoBrush)
        painter.setPen(QPen(grid_col, 0.7))
        for frac in (0.33, 0.67, 1.0):
            ring = QPolygonF()
            for i in range(n):
                a = -math.pi / 2 + 2 * math.pi * i / n
                ring.append(QPointF(frac * r * math.cos(a), frac * r * math.sin(a)))
            ring.append(ring[0])
            painter.drawPolygon(ring)

        # Spokes
        painter.setPen(QPen(grid_col, 0.6))
        for i in range(n):
            a = -math.pi / 2 + 2 * math.pi * i / n
            painter.drawLine(QPointF(0, 0), QPointF(r * math.cos(a), r * math.sin(a)))

        # Data polygon
        poly = QPolygonF()
        for i, v in enumerate(self._spider_vals):
            a   = -math.pi / 2 + 2 * math.pi * i / n
            r_v = max(0.05, v) * r
            poly.append(QPointF(r_v * math.cos(a), r_v * math.sin(a)))
        poly.append(poly[0])
        painter.setBrush(QBrush(fill_col))
        painter.setPen(QPen(edge_col, 1.2))
        painter.drawPolygon(poly)

        # Outer border circle
        painter.setBrush(Qt.NoBrush)
        painter.setPen(QPen(base_col.darker(130), 1.0))
        painter.drawEllipse(QRectF(-r, -r, r * 2, r * 2))

    def _draw_coxcomb(self, painter: "QPainter", r: float) -> None:
        """Draw equal-angle sectors whose radius is proportional to sqrt(count)."""
        n = len(self._pie_slices)
        if n == 0:
            return
        max_cnt = max(cnt for _, cnt in self._pie_slices) or 1
        sector_deg = 360.0 / n
        start_deg = 90.0  # start at top, go clockwise

        for i, (sid, cnt) in enumerate(self._pie_slices):
            sector_r = r * math.sqrt(cnt / max_cnt)
            sector_rect = QRectF(-sector_r, -sector_r, sector_r * 2, sector_r * 2)
            angle_qt = int(round((start_deg - i * sector_deg) * 16))
            span_qt = int(round(-sector_deg * 16))  # negative = clockwise in Qt
            sc = self._sample_colors.get(sid, self._color)
            painter.setBrush(QBrush(sc))
            painter.setPen(QPen(sc.darker(130), 0.6))
            painter.drawPie(sector_rect, angle_qt, span_qt)

        # Outer reference circle (bounding radius)
        painter.setBrush(Qt.NoBrush)
        painter.setPen(QPen(self._color.darker(150), 1.0))
        painter.drawEllipse(QRectF(-r, -r, r * 2, r * 2))

    def _draw_virtual_node(self, painter: "QPainter", r: float) -> None:
        """Hollow dashed circle for virtual (unclassified/excluded) branch nodes."""
        rect = QRectF(-r, -r, r * 2, r * 2)
        painter.setBrush(QBrush(_VIRTUAL_NODE_COLOR))
        pen = QPen(_VIRTUAL_NODE_BORDER, 1.4)
        pen.setStyle(Qt.DashLine)
        painter.setPen(pen)
        painter.drawEllipse(rect)
        # Small crosshair inside to distinguish from real nodes
        if r >= 5:
            arm = r * 0.35
            painter.setPen(QPen(_VIRTUAL_NODE_BORDER, 0.8))
            painter.drawLine(QPointF(-arm, 0.0), QPointF(arm, 0.0))
            painter.drawLine(QPointF(0.0, -arm), QPointF(0.0, arm))

    def paint(self, painter, option, widget=None):
        try:
            self._paint_impl(painter)
        except Exception:
            pass

    def _paint_impl(self, painter):
        painter.setRenderHint(QPainter.Antialiasing)
        r = self._radius
        rect = QRectF(-r, -r, r * 2, r * 2)

        if self._is_virtual_node:
            self._draw_virtual_node(painter, r)
            # Still draw selection/search rings on top
            if self._search_match:
                painter.setBrush(Qt.NoBrush)
                painter.setPen(QPen(QColor("#c8a010"), 2.5))
                painter.drawEllipse(rect.adjusted(-4, -4, 4, 4))
            if self.isSelected() or self._selected_override:
                painter.setBrush(Qt.NoBrush)
                painter.setPen(QPen(QColor("#ffffff"), 2.0))
                painter.drawEllipse(rect.adjusted(-2, -2, 2, 2))
            if self._multi_selected:
                pen = QPen(QColor("#00e5ff"), 2.2)
                pen.setStyle(Qt.DashLine)
                painter.setBrush(Qt.NoBrush)
                painter.setPen(pen)
                painter.drawEllipse(rect.adjusted(-3, -3, 3, 3))
            return

        if self._display_mode == "spider":
            self._draw_spider_node(painter, r)
        elif self._is_pie and self._pie_slices:
            if self._display_mode == "coxcomb":
                self._draw_coxcomb(painter, r)
            else:
                total = sum(cnt for _, cnt in self._pie_slices)
                angle = int(90 * 16)
                for sid, cnt in self._pie_slices:
                    span = int(round(360 * 16 * cnt / total))
                    sc = self._sample_colors.get(sid, self._color)
                    painter.setBrush(QBrush(sc))
                    painter.setPen(QPen(Qt.NoPen))
                    painter.drawPie(rect, angle, span)
                    angle += span
                # Border
                painter.setBrush(Qt.NoBrush)
                painter.setPen(QPen(self._color.darker(140), 1.2))
                painter.drawEllipse(rect)
        else:
            painter.setBrush(QBrush(self._color))
            if self._dark_taxon:
                painter.setPen(QPen(_DARK_TAXON_BORDER, 1.5))
            else:
                painter.setPen(QPen(self._color.darker(140), 1.2))
            painter.drawEllipse(rect)
            # Draw a hollow-circle "ghost" icon for dark taxa (no ref sequences)
            if self._dark_taxon:
                painter.setBrush(Qt.NoBrush)
                painter.setPen(QPen(_DARK_TAXON_BORDER, 1.0))
                inner_r = max(2.0, r * 0.45)
                painter.drawEllipse(QRectF(-inner_r, -inner_r, inner_r * 2, inner_r * 2))

        # Search match ring (yellow, outermost)
        if self._search_match:
            painter.setBrush(Qt.NoBrush)
            painter.setPen(QPen(QColor("#c8a010"), 2.5))
            painter.drawEllipse(rect.adjusted(-4, -4, 4, 4))

        # Primary selection ring (white solid — single-select primary)
        if self.isSelected() or self._selected_override:
            painter.setBrush(Qt.NoBrush)
            painter.setPen(QPen(QColor("#ffffff"), 2.0))
            painter.drawEllipse(rect.adjusted(-2, -2, 2, 2))

        # Multi-selection ring (cyan dashed — member of multi-select set)
        if self._multi_selected:
            pen = QPen(QColor("#00e5ff"), 2.2)
            pen.setStyle(Qt.DashLine)
            painter.setBrush(Qt.NoBrush)
            painter.setPen(pen)
            painter.drawEllipse(rect.adjusted(-3, -3, 3, 3))

        # Review mark badge (top-right of node)
        if self._review_mark in _REVIEW_COLORS:
            badge_r = max(4.0, r * 0.45)
            bx = r - badge_r * 0.5
            by = -r - badge_r * 0.5
            badge_color = _REVIEW_COLORS[self._review_mark]
            painter.setBrush(QBrush(badge_color))
            painter.setPen(QPen(QColor("#0d1117"), 1.0))
            painter.drawEllipse(QRectF(bx - badge_r, by - badge_r, badge_r * 2, badge_r * 2))
            font = QFont("Segoe UI", max(5, int(badge_r * 1.1)), QFont.Bold)
            painter.setFont(font)
            painter.setPen(QPen(QColor("#0d1117")))
            painter.drawText(
                QRectF(bx - badge_r, by - badge_r, badge_r * 2, badge_r * 2),
                Qt.AlignCenter, _REVIEW_LABELS[self._review_mark]
            )

        # Authenticity tier badge (bottom-left of node)
        if self._auth_badge and self._auth_tier is not None:
            tier_color = _AUTH_TIER_COLORS.get(self._auth_tier, QColor("#7890a0"))
            badge_r = max(4.0, r * 0.42)
            bx = -r + badge_r * 0.5
            by = r - badge_r * 0.5
            painter.setBrush(QBrush(tier_color))
            painter.setPen(QPen(QColor("#0d1117"), 0.8))
            painter.drawEllipse(QRectF(bx - badge_r, by - badge_r, badge_r * 2, badge_r * 2))
            font = QFont("Segoe UI", max(5, int(badge_r * 1.0)))
            painter.setFont(font)
            painter.setPen(QPen(QColor("#0d1117")))
            painter.drawText(
                QRectF(bx - badge_r, by - badge_r, badge_r * 2, badge_r * 2),
                Qt.AlignCenter, self._auth_badge
            )

        # Fossil support badge (bottom-right of node) — shown when fos_support=True
        if self._fos_support:
            badge_r = max(4.0, r * 0.42)
            bx = r - badge_r * 0.5
            by = r - badge_r * 0.5
            painter.setBrush(QBrush(QColor("#7b5ea7")))  # muted purple
            painter.setPen(QPen(QColor("#0d1117"), 0.8))
            painter.drawEllipse(QRectF(bx - badge_r, by - badge_r, badge_r * 2, badge_r * 2))
            font = QFont("Segoe UI", max(4, int(badge_r * 0.85)))
            painter.setFont(font)
            painter.setPen(QPen(QColor("#e8d5ff")))
            painter.drawText(
                QRectF(bx - badge_r, by - badge_r, badge_r * 2, badge_r * 2),
                Qt.AlignCenter, "\U0001F9B4"  # 🦴
            )

    def boundingRect(self) -> QRectF:
        r = self._radius + 14  # extra room for search/review/auth badge rings
        return QRectF(-r, -r, r * 2, r * 2)


class _CollapseTriangle(QGraphicsPolygonItem):
    """Small triangle indicator for collapsed nodes."""

    def __init__(self, radius: float, color: QColor, pointing_right: bool = True):
        if pointing_right:
            pts = QPolygonF([
                QPointF(radius + 3, 0),
                QPointF(radius + 10, -5),
                QPointF(radius + 10, 5),
            ])
        else:
            # Pointing down (top-to-bottom layout)
            pts = QPolygonF([
                QPointF(0, radius + 3),
                QPointF(-5, radius + 10),
                QPointF(5, radius + 10),
            ])
        super().__init__(pts)
        self.setBrush(QBrush(color))
        self.setPen(QPen(Qt.NoPen))


# ─── Dark taxon visual constants ───────────────────────────────────────────────
# These mark nodes whose taxid was in the DB build BFS expansion but ended up
# with 0 qualifying reference sequences.  Reads assigned here go to the nearest
# represented ancestor.  Visually shown as a near-black fill with a muted border.

_DARK_TAXON_COLOR  = QColor("#1a1a2e")   # near-black navy
_DARK_TAXON_BORDER = QColor("#6c6c8a")   # muted purple-grey border

_VIRTUAL_NODE_COLOR  = QColor("#0f1520")  # near-black blue-grey
_VIRTUAL_NODE_BORDER = QColor("#3a5060")  # muted teal border (dashed in paint)


def _build_virtual_tooltip(node: "_LayoutNode") -> str:
    """HTML tooltip for virtual (unclassified/excluded) branch nodes."""
    DIM = "#7890a0"
    HDR = "#3a8cc8"
    FG  = "#c8d4e0"
    reads = node.reads
    desc  = node.virtual_tooltip or ""
    # Replace newlines with HTML line breaks
    desc_html = desc.replace("\n\n", "<br/><br/>").replace("\n", "<br/>")
    tip  = f'<b style="color:{HDR};">{node.name}</b>'
    tip += f'  <span style="color:{DIM};">│  virtual category</span><br/>'
    if reads > 0:
        tip += f'<span style="color:{FG};">{reads:,} reads</span><br/>'
    else:
        tip += f'<span style="color:{DIM};">0 reads in this category</span><br/>'
    if desc_html:
        tip += f'<hr style="border:0; border-top:1px solid #2c3a4a; margin:4px 0;"/>'
        tip += f'<span style="color:{FG};">{desc_html}</span>'
    return tip


class _ShortcutsOverlay(QWidget):
    """Semi-transparent shortcut-key hint pinned to the bottom-left of the tree viewport."""

    _TEXT = (
        "Mouse shortcuts\n"
        "  Left-drag          Pan\n"
        "  Middle-drag        Pan\n"
        "  Ctrl+scroll        Zoom\n"
        "  Shift+scroll       H-stretch\n"
        "  Ctrl+Shift+scroll  V-stretch\n"
        "  Double-click       Expand / collapse\n"
        "  Right-click        Context menu\n"
        "  Ctrl+click         Multi-select\n"
        "\n"
        "Keys\n"
        "  Ctrl++   Expand all\n"
        "  Ctrl+-   Collapse all\n"
        "  Home     Fit to window\n"
        "  /        Focus search\n"
        "  Esc      Clear search\n"
    )

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setAttribute(Qt.WA_TransparentForMouseEvents, False)
        self.setAttribute(Qt.WA_TranslucentBackground)
        self.setWindowFlags(Qt.SubWindow)
        self._font = QFont("Consolas", 8)
        self._update_size()

    def _update_size(self) -> None:
        from qtpy.QtGui import QFontMetrics
        fm = QFontMetrics(self._font)
        lines = self._TEXT.split("\n")
        w = max(fm.horizontalAdvance(l) for l in lines) + 24
        h = fm.height() * len(lines) + 16
        self.setFixedSize(w, h)

    def paintEvent(self, event) -> None:
        from qtpy.QtGui import QPainter, QColor, QBrush, QPen
        p = QPainter(self)
        p.setRenderHint(QPainter.Antialiasing)
        p.setBrush(QBrush(QColor(10, 15, 30, 200)))
        p.setPen(QPen(QColor(50, 80, 120, 180), 1))
        p.drawRoundedRect(self.rect().adjusted(1, 1, -1, -1), 6, 6)
        p.setFont(self._font)
        p.setPen(QPen(QColor(160, 200, 230, 220)))
        p.drawText(self.rect().adjusted(8, 6, -8, -6), Qt.AlignLeft | Qt.AlignTop, self._TEXT)
        p.end()

    def reposition(self, viewport_size) -> None:
        margin = 8
        self.move(margin, viewport_size.height() - self.height() - margin)


class _OverlayPanel(QWidget):
    """Draggable floating overlay panel — child of the viewport."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self._drag_start: Optional[QPointF] = None

    def mousePressEvent(self, event) -> None:
        if event.button() == Qt.LeftButton:
            self._drag_start = QPointF(event.pos())
            event.accept()
        else:
            super().mousePressEvent(event)

    def mouseMoveEvent(self, event) -> None:
        if self._drag_start is not None and event.buttons() & Qt.LeftButton:
            delta = event.pos() - self._drag_start
            new_pos = self.pos() + QPoint(int(delta.x()), int(delta.y()))
            # Clamp to parent bounds
            if self.parent():
                pw, ph = self.parent().width(), self.parent().height()
                new_pos.setX(max(0, min(pw - self.width(),  new_pos.x())))
                new_pos.setY(max(0, min(ph - self.height(), new_pos.y())))
            self.move(new_pos)
            event.accept()
        else:
            super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event) -> None:
        if event.button() == Qt.LeftButton:
            self._drag_start = None
            event.accept()
        else:
            super().mouseReleaseEvent(event)


# ─── Main view ─────────────────────────────────────────────────────────────────

_REVIEW_COLORS = {
    "reviewed":     QColor("#50b840"),   # green
    "uncertain":    QColor("#c8a010"),   # yellow
    "rejected":     QColor("#e03838"),   # red
    "contaminant":  QColor("#b040c0"),   # purple
}
_REVIEW_LABELS = {
    "reviewed":    "✓",
    "uncertain":   "?",
    "rejected":    "✗",
    "contaminant": "C",
}


class TaxonGraphicsTree(QGraphicsView):
    """MEGAN7-style taxonomy tree rendered in a QGraphicsView."""

    node_selected = Signal(dict)
    selection_cleared = Signal()           # emitted when clicking on empty space to deselect
    extract_requested = Signal(list)       # list of taxid strings
    mark_changed = Signal(str, str)        # taxid, mark ("reviewed"/"uncertain"/"rejected"/"")
    multi_selection_changed = Signal(list) # list of selected taxid strings (empty = cleared)
    subtree_focus_changed = Signal(bool)   # True = focused view active, False = full tree
    show_in_list_requested = Signal(str)      # taxid — switch to list view and scroll to taxon
    rerun_requested = Signal()                # open RunAnalysisDialog
    bookmark_requested = Signal(str, str)     # taxid, name — add to bookmarks
    reclassify_requested = Signal(list, str)  # [taxids], taxon_name — re-classify subtree

    def __init__(self, parent=None):
        super().__init__(parent)
        self._scene = QGraphicsScene(self)
        self.setScene(self._scene)
        self.setRenderHint(QPainter.Antialiasing)
        self.setDragMode(QGraphicsView.NoDrag)
        self.setTransformationAnchor(QGraphicsView.AnchorUnderMouse)
        self.setResizeAnchor(QGraphicsView.AnchorViewCenter)
        self.setBackgroundBrush(QBrush(C_TREE_BG))
        self.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        self.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)

        self._samples: List[str] = []
        self._selected_samples: Optional[List[str]] = None
        self._sample_colors: Dict[str, QColor] = {}
        self._max_reads: int = 1
        self._root: Optional[_LayoutNode] = None
        self._node_items: Dict[str, _NodeItem] = {}   # taxid -> _NodeItem
        self._selected_taxid: Optional[str] = None
        self._multi_selected_taxids: Set[str] = set()   # Ctrl+Click set
        self._focus_taxids: Optional[Set[str]] = None   # subtree focus filter
        self._scale_mode: str = "sqrt"
        self._layout_mode: str = "cladogram"   # "cladogram" | "dendrogram"
        self._ltr_layout: bool = True          # left-to-right (True) or top-to-bottom (False)
        self._node_display_mode: str = "pie"   # "pie" | "coxcomb" | "spider"
        self._color_by: str = "support"        # "support" | "damage" | "tier" | "blank" | "posterior"
        self._label_mode: str = "name_reads"   # "name" | "name_reads" | "name_damage" | "name_tier"
        self._count_display: str = "cumulative"  # "cumulative" | "direct" | "both"
        self._label_visibility: str = "direct_major"  # "all" | "direct" | "direct_major"
        self._max_direct_reads: int = 1          # max direct_reads across all nodes
        self._node_size_factor: float = 1.0     # user-adjustable node radius multiplier
        self._review_marks: Dict[str, str] = {}    # taxid -> "reviewed"/"uncertain"/"rejected"
        self._dark_taxids: Set[str] = set()         # taxids with no reference sequences in DB

        # Comparison mode
        self._comparison_mode: bool = False
        self._compare_a: Optional[str] = None  # sample id A
        self._compare_b: Optional[str] = None  # sample id B

        self._branch_style: str = "square"   # "square" | "rounded"

        # Label overlap tracking (reset before each draw pass)
        self._placed_label_rects: List[QRectF] = []

        # Stretch scale factors — modified by Shift+scroll / Ctrl+Shift+scroll
        # and applied to layout spacing constants before each relayout.
        self._h_scale: float = 1.0   # horizontal (leaf spread)
        self._v_scale: float = 1.0   # vertical (depth spacing)
        self._eff_h_unit: float = _H_UNIT
        self._eff_v_spacing: float = _V_SPACING

        # Pan state — left-button: track drag vs click; middle-button: always pan
        self._left_pan_start: Optional[QPointF] = None
        self._left_pan_dragging: bool = False
        self._mid_pan_active: bool = False
        self._mid_pan_start: QPointF = QPointF()
        # Deferred deselect: fires after short delay so dragging from empty space
        # doesn't clear selection before pan is confirmed (avoids dock-resize jitter).
        self._deselect_timer = QTimer(self)
        self._deselect_timer.setSingleShot(True)
        self._deselect_timer.setInterval(80)
        self._deselect_timer.timeout.connect(self._deselect_all)

        # Search
        self._search_text: str = ""

        # Shortcut hint overlay
        self._shortcuts_overlay = _ShortcutsOverlay(self.viewport())
        self._shortcuts_overlay.hide()

        # Font (size adjustable via d-pad A+/A− buttons)
        self._font_pt: int = 9
        self._label_font = QFont("Segoe UI", self._font_pt)

        # Animation
        self._anim_timeline: Optional["QTimeLine"] = None  # node slide animation

        # Debounce timer: after zoom stops, relayout labels at the new scale
        self._zoom_relabel_timer = QTimer(self)
        self._zoom_relabel_timer.setSingleShot(True)
        self._zoom_relabel_timer.setInterval(280)
        self._zoom_relabel_timer.timeout.connect(self._relayout)

        # Cursor reset timer: restores default cursor 600ms after wheel-stretch gesture
        self._cursor_reset_timer = QTimer(self)
        self._cursor_reset_timer.setSingleShot(True)
        self._cursor_reset_timer.setInterval(600)
        self._cursor_reset_timer.timeout.connect(self.viewport().unsetCursor)

        # Viewport overlay buttons (Fit, Expand, Collapse, Redraw)
        self._setup_overlay_buttons()

    # ─── Overlay buttons ───────────────────────────────────────────────────────

    def _setup_overlay_buttons(self) -> None:
        """Floating overlay: Spacing d-pad + Pan/Zoom d-pad + action row + legend."""
        from qtpy.QtWidgets import (QPushButton, QGridLayout, QHBoxLayout,
                                    QVBoxLayout, QLabel, QFrame)

        # ── hold-to-pan / hold-to-zoom timers ─────────────────────────────
        self._pan_hold_timer = QTimer(self)
        self._pan_hold_timer.setInterval(80)
        self._pan_hold_dx = 0
        self._pan_hold_dy = 0
        self._pan_hold_timer.timeout.connect(
            lambda: self._pan(self._pan_hold_dx, self._pan_hold_dy)
        )
        self._zoom_hold_timer = QTimer(self)
        self._zoom_hold_timer.setInterval(150)  # slightly slower than pan (100→150ms)
        self._zoom_hold_factor = 1.0
        self._zoom_hold_timer.timeout.connect(self._zoom_hold_tick)

        SS = (
            "_OverlayPanel { background: rgba(10,15,28,218); border-radius: 6px; }"
            "QPushButton { background: rgba(30,45,65,200); color: #dde8f0; "
            "              border: 1px solid #2a4060; font-size: 9px; border-radius: 3px; }"
            "QPushButton:hover  { color: #ffffff; background: rgba(50,80,110,220); }"
            "QPushButton:pressed{ background: rgba(0,150,200,140); }"
            "QLabel { color: #5a7898; font-size: 7pt; background: transparent; border: none; }"
        )
        panel = _OverlayPanel(self.viewport())
        panel.setStyleSheet(SS)
        outer = QVBoxLayout(panel)
        outer.setContentsMargins(6, 6, 6, 6)
        outer.setSpacing(4)

        # ── helper: normal click button ────────────────────────────────────
        def _btn(label: str, tip: str, fn, w: int = 36, h: int = 20) -> QPushButton:
            b = QPushButton(label)
            b.setFixedSize(w, h)
            b.setToolTip(tip)
            b.clicked.connect(fn)
            return b

        # ── helper: press-and-hold button (fires fn on press, repeats via timer) ─
        def _hold_btn(label: str, tip: str, start_fn, stop_fn=None,
                      w: int = 26, h: int = 20) -> QPushButton:
            b = QPushButton(label)
            b.setFixedSize(w, h)
            b.setToolTip(tip)
            b.pressed.connect(start_fn)
            b.released.connect(stop_fn or (lambda: None))
            return b

        # ── two d-pads side by side ────────────────────────────────────────
        pads_row = QHBoxLayout()
        pads_row.setSpacing(8)

        # ── Spacing d-pad ─────────────────────────────────────────────────
        # In TTB (top-to-bottom) mode: V = depth axis (vertical), H = leaf axis (horizontal).
        # In LTR (left-to-right) mode: depth is horizontal and leaf spread is vertical, so
        # buttons labelled "V+/-" should control whatever currently appears vertical on screen.
        # _update_spacing_buttons() swaps connections when orientation changes.
        sp_frame = QFrame()
        sp_frame.setStyleSheet("QFrame{background:transparent;border:none;}")
        rg = QGridLayout(sp_frame)
        rg.setContentsMargins(0, 0, 0, 0)
        rg.setSpacing(1)
        self._sp_top_btn    = _btn("↑ +", "More vertical screen spacing", self._v_expand, 36, 20)
        self._sp_left_btn   = _btn("← −", "Less horizontal screen spacing", self._h_compress, 36, 20)
        self._sp_reset_btn  = _btn("↺", "Reset spacing to default", self._reset_scales, 26, 20)
        self._sp_right_btn  = _btn("→ +", "More horizontal screen spacing", self._h_expand, 36, 20)
        self._sp_bottom_btn = _btn("↓ −", "Less vertical screen spacing", self._v_compress, 36, 20)
        rg.addWidget(self._sp_top_btn,    0, 1)
        rg.addWidget(self._sp_left_btn,   1, 0)
        rg.addWidget(self._sp_reset_btn,  1, 1)
        rg.addWidget(self._sp_right_btn,  1, 2)
        rg.addWidget(self._sp_bottom_btn, 2, 1)
        sp_lbl = QLabel("Spacing")
        sp_lbl.setAlignment(Qt.AlignCenter)
        sp_col = QVBoxLayout()
        sp_col.setSpacing(1)
        sp_col.addWidget(sp_lbl)
        sp_col.addWidget(sp_frame)
        self._update_spacing_buttons()  # wire correct functions for current layout

        # ── Pan / Zoom d-pad ───────────────────────────────────────────────
        pz_frame = QFrame()
        pz_frame.setStyleSheet("QFrame{background:transparent;border:none;}")
        pg = QGridLayout(pz_frame)
        pg.setContentsMargins(0, 0, 0, 0)
        pg.setSpacing(1)

        def _start_pan(dx, dy):
            self._pan_hold_dx, self._pan_hold_dy = dx, dy
            self._pan(dx, dy)          # fire immediately
            self._pan_hold_timer.start()

        def _stop_pan():
            self._pan_hold_timer.stop()

        def _start_zoom(factor):
            self._zoom_hold_factor = factor
            self._zoom_hold_tick()     # fire immediately
            self._zoom_hold_timer.start()

        def _stop_zoom():
            self._zoom_hold_timer.stop()

        # Zoom-in (+) at top-left, Zoom-out (−) at bottom-right
        pg.addWidget(
            _hold_btn("Zoom+", "Zoom in (Ctrl+scroll up)",
                      lambda: _start_zoom(1.20), _stop_zoom, 36, 20), 0, 0)
        pg.addWidget(
            _hold_btn("▲", "Pan up",
                      lambda: _start_pan(0, -80), _stop_pan, 26, 20), 0, 1)
        # top-right corner: fit
        pg.addWidget(
            _btn("⊡ Fit", "Fit tree to window (Home)", self._fit_all, 36, 20), 0, 2)
        pg.addWidget(
            _hold_btn("◄", "Pan left",
                      lambda: _start_pan(-80, 0), _stop_pan, 26, 20), 1, 0)
        # center: empty spacer (use empty label for alignment)
        _spacer = QLabel("")
        _spacer.setFixedSize(26, 20)
        pg.addWidget(_spacer, 1, 1)
        pg.addWidget(
            _hold_btn("►", "Pan right",
                      lambda: _start_pan(80, 0), _stop_pan, 26, 20), 1, 2)
        pg.addWidget(
            _hold_btn("Zoom−", "Zoom out (Ctrl+scroll down)",
                      lambda: _start_zoom(1/1.20), _stop_zoom, 36, 20), 2, 0)
        pg.addWidget(
            _hold_btn("▼", "Pan down",
                      lambda: _start_pan(0, 80), _stop_pan, 26, 20), 2, 1)

        pz_lbl = QLabel("Pan / Zoom")
        pz_lbl.setAlignment(Qt.AlignCenter)
        pz_col = QVBoxLayout()
        pz_col.setSpacing(1)
        pz_col.addWidget(pz_lbl)
        pz_col.addWidget(pz_frame)

        pads_row.addLayout(sp_col)
        pads_row.addLayout(pz_col)
        outer.addLayout(pads_row)

        # ── action row ─────────────────────────────────────────────────────
        act_row = QHBoxLayout()
        act_row.setSpacing(2)

        def _expand_all():
            if self._root:
                self._set_collapsed(self._root, False)
                self._relayout()
                self._auto_compact_fit()

        act_row.addWidget(_btn("Expand All",   "Expand all nodes (Ctrl++)",
                               _expand_all, 72, 22))
        act_row.addWidget(_btn("Collapse All", "Collapse all nodes (Ctrl+-)",
                               self.collapse_all, 72, 22))
        outer.addLayout(act_row)

        # ── font size row ──────────────────────────────────────────────────
        font_row = QHBoxLayout()
        font_row.setSpacing(2)
        font_lbl = QLabel("Font:")
        font_lbl.setStyleSheet("color: #5a7898; font-size: 7pt; background: transparent;")
        font_row.addWidget(font_lbl)
        font_row.addWidget(_btn("A+", "Increase label font size", self._font_up,  32, 20))
        font_row.addWidget(_btn("A−", "Decrease label font size", self._font_down, 32, 20))
        font_row.addStretch()
        outer.addLayout(font_row)

        # ── node size row ──────────────────────────────────────────────────
        sz_row = QHBoxLayout()
        sz_row.setSpacing(2)
        sz_lbl = QLabel("Node:")
        sz_lbl.setStyleSheet("color: #5a7898; font-size: 7pt; background: transparent;")
        sz_row.addWidget(sz_lbl)
        sz_row.addWidget(_btn("N+", "Increase base node size (scales all nodes up)",
                              self._node_size_up,    32, 20))
        sz_row.addWidget(_btn("N−", "Decrease base node size (scales all nodes down)",
                              self._node_size_down,  32, 20))
        sz_row.addWidget(_btn("N↺", "Reset node size to default",
                              self._node_size_reset, 26, 20))
        sz_row.addStretch()
        outer.addLayout(sz_row)

        # ── shortcut legend ────────────────────────────────────────────────
        tips = QLabel(
            "Ctrl+scroll Zoom  ·  Shift+scroll H-space\n"
            "Ctrl+Shift+scroll V-space  ·  Home Fit\n"
            "Ctrl++ Expand  ·  Ctrl+- Collapse\n"
            "Dbl-click Expand/collapse  ·  Ctrl+click Multi-sel"
        )
        tips.setWordWrap(True)
        outer.addWidget(tips)

        panel.adjustSize()
        panel.show()
        self._overlay_panel = panel
        self._reposition_overlay()

    # ── overlay action helpers ─────────────────────────────────────────────────

    def _pan(self, dx: int, dy: int) -> None:
        self.horizontalScrollBar().setValue(self.horizontalScrollBar().value() + dx)
        self.verticalScrollBar().setValue(self.verticalScrollBar().value() + dy)

    def _zoom_hold_tick(self) -> None:
        """Called by the zoom hold timer to apply one zoom step."""
        f = self._zoom_hold_factor
        cur = self.transform().m11()
        if f > 1 and cur * f > 200.0:
            return
        if f < 1 and cur * f < 0.01:
            return
        self.scale(f, f)
        self._zoom_relabel_timer.start()

    def _h_expand(self) -> None:
        self._h_scale = min(8.0, self._h_scale * 1.25)
        self._scale_relayout()

    def _h_compress(self) -> None:
        self._h_scale = max(0.08, self._h_scale / 1.25)
        self._scale_relayout()

    def _v_expand(self) -> None:
        self._v_scale = min(8.0, self._v_scale * 1.25)
        self._scale_relayout()

    def _v_compress(self) -> None:
        self._v_scale = max(0.08, self._v_scale / 1.25)
        self._scale_relayout()

    def _reset_scales(self) -> None:
        self._h_scale = 1.0
        self._v_scale = 1.0
        self._scale_relayout()
        self._fit_all()

    def _font_up(self) -> None:
        self._font_pt = min(36, self._font_pt + 1)
        self._label_font = QFont("Segoe UI", self._font_pt)
        self._scale_relayout()

    def _font_down(self) -> None:
        self._font_pt = max(5, self._font_pt - 1)
        self._label_font = QFont("Segoe UI", self._font_pt)
        self._scale_relayout()

    def _node_size_up(self) -> None:
        self._node_size_factor = min(8.0, self._node_size_factor * 1.4)
        self._relayout()

    def _node_size_down(self) -> None:
        self._node_size_factor = max(0.1, self._node_size_factor / 1.4)
        self._relayout()

    def _node_size_reset(self) -> None:
        self._node_size_factor = 1.0
        self._relayout()

    def _update_spacing_buttons(self) -> None:
        """Wire spacing d-pad buttons so top/bottom control screen-vertical spacing
        and left/right control screen-horizontal spacing, regardless of tree orientation.

        In TTB layout: depth is vertical (v_scale) → top/bottom; leaf spread is
        horizontal (h_scale) → left/right. In LTR layout the axes are transposed.
        """
        if not hasattr(self, "_sp_top_btn"):
            return
        # Disconnect all first
        for btn, fn in [
            (self._sp_top_btn,    None),
            (self._sp_bottom_btn, None),
            (self._sp_left_btn,   None),
            (self._sp_right_btn,  None),
        ]:
            try:
                btn.clicked.disconnect()
            except RuntimeError:
                pass

        if self._ltr_layout:
            # LTR: depth is horizontal (v_scale), leaf spread is vertical (h_scale)
            # ↑ = more vertical (leaf spread); ↓ = less vertical
            self._sp_top_btn.clicked.connect(self._h_expand)
            self._sp_bottom_btn.clicked.connect(self._h_compress)
            self._sp_left_btn.clicked.connect(self._v_compress)
            self._sp_right_btn.clicked.connect(self._v_expand)
            self._sp_top_btn.setToolTip("More vertical spacing (leaf spread)\n(Shift+scroll up)")
            self._sp_bottom_btn.setToolTip("Less vertical spacing (leaf spread)\n(Shift+scroll down)")
            self._sp_left_btn.setToolTip("Less horizontal spacing (depth)\n(Ctrl+Shift+scroll down)")
            self._sp_right_btn.setToolTip("More horizontal spacing (depth)\n(Ctrl+Shift+scroll up)")
        else:
            # TTB: depth is vertical (v_scale), leaf spread is horizontal (h_scale)
            # ↑ = more vertical (depth); ↓ = less vertical
            self._sp_top_btn.clicked.connect(self._v_expand)
            self._sp_bottom_btn.clicked.connect(self._v_compress)
            self._sp_left_btn.clicked.connect(self._h_compress)
            self._sp_right_btn.clicked.connect(self._h_expand)
            self._sp_top_btn.setToolTip("More vertical spacing (depth)\n(Ctrl+Shift+scroll up)")
            self._sp_bottom_btn.setToolTip("Less vertical spacing (depth)\n(Ctrl+Shift+scroll down)")
            self._sp_left_btn.setToolTip("Less horizontal spacing (leaf spread)\n(Shift+scroll down)")
            self._sp_right_btn.setToolTip("More horizontal spacing (leaf spread)\n(Shift+scroll up)")

    def _scale_relayout(self) -> None:
        """Instant (no-animation) relayout for scale changes — prevents nodes lagging behind branches."""
        if self._root is None:
            return
        if self._anim_timeline is not None:
            self._anim_timeline.stop()
            self._anim_timeline = None
        self._eff_h_unit = _H_UNIT * self._h_scale
        self._eff_v_spacing = _V_SPACING * self._v_scale
        _assign_x(self._root, 0.0)
        self._apply_y_layout(self._root)
        self._scene.clear()
        self._node_items.clear()
        self._placed_label_rects = []
        self._draw_counter = 0
        self._scene.setItemIndexMethod(QGraphicsScene.NoIndex)
        self._draw_subtree(self._root, None)
        self._scene.setItemIndexMethod(QGraphicsScene.BspTreeIndex)
        self._scene.setSceneRect(self._scene.itemsBoundingRect().adjusted(-20, -20, 20, 20))
        self._apply_search_highlights()
        if hasattr(self, "_overlay_panel"):
            self._overlay_panel.raise_()

    def _reposition_overlay(self) -> None:
        if hasattr(self, "_overlay_panel"):
            self._overlay_panel.adjustSize()
            vp = self.viewport()
            # Pin to top-left corner; never allow the panel to slide off-screen
            pw = self._overlay_panel.width()
            ph = self._overlay_panel.height()
            x = min(8, max(0, vp.width()  - pw - 4))
            y = min(8, max(0, vp.height() - ph - 4))
            self._overlay_panel.move(x, y)
            self._overlay_panel.raise_()

    # ─── Public API ────────────────────────────────────────────────────────────

    def set_dark_taxids(self, dark_taxids: Set[str]) -> None:
        """Set the taxids that had no reference sequences in the DB build.

        Call this after loading a DB (e.g. from DarkTaxidPanel data) so that
        matching nodes are rendered with the ghost styling in _draw_subtree.
        """
        self._dark_taxids = set(dark_taxids)

    def load_tree(
        self,
        tree: dict,
        samples: List[str],
        selected_samples: Optional[List[str]] = None,
        support_filter: str = "all",
        sample_colors: Optional[Dict] = None,
    ) -> None:
        self._samples = samples
        self._selected_samples = selected_samples
        if sample_colors:
            self._sample_colors = dict(sample_colors)
        else:
            self._sample_colors = {
                sid: SAMPLE_PALETTE[i % len(SAMPLE_PALETTE)]
                for i, sid in enumerate(samples)
            }
        self._max_reads = max(1, _max_reads_in_tree(tree))
        self._support_filter = support_filter
        self._build_scene(tree)

    def set_selected_samples(self, selected: List[str]) -> None:
        self._selected_samples = selected
        if self._root is not None:
            self._rebuild_scene()

    def set_scale_mode(self, mode: str) -> None:
        """Set node radius scaling: 'sqrt', 'linear', or 'log'."""
        if mode != self._scale_mode:
            self._scale_mode = mode
            self._relayout()

    def set_layout_mode(self, mode: str) -> None:
        """Set tree layout: 'cladogram' (depth-based) or 'dendrogram' (right-aligned leaves)."""
        if mode != self._layout_mode:
            self._layout_mode = mode
            self._relayout()

    def set_ltr_layout(self, ltr: bool) -> None:
        """Set orientation: True = left-to-right (root at left), False = top-to-bottom."""
        if ltr != self._ltr_layout:
            self._ltr_layout = ltr
            self._update_spacing_buttons()
            self._relayout()

    def refresh_theme(self) -> None:
        """Re-apply current global theme colors to background and redraw scene."""
        self.setBackgroundBrush(QBrush(C_TREE_BG))
        if self._root is not None:
            self._rebuild_scene()

    def set_node_display_mode(self, mode: str) -> None:
        """Set node drawing style: 'pie' (proportional angle) or 'coxcomb' (proportional radius)."""
        if mode != self._node_display_mode:
            self._node_display_mode = mode
            self._rebuild_scene()

    def set_color_by(self, mode: str) -> None:
        """Set node color-by mode: 'support', 'damage', 'tier', 'blank', 'posterior'."""
        if mode != self._color_by:
            self._color_by = mode
            self._relayout()

    def set_label_mode(self, mode: str) -> None:
        """Set node label mode: 'name', 'name_reads', 'name_damage', 'name_tier'."""
        if mode != self._label_mode:
            self._label_mode = mode
            self._relayout()

    def set_count_display(self, mode: str) -> None:
        """Set count display mode: 'cumulative', 'direct', or 'both'.

        cumulative: label and node size use sum of reads across node + all descendants.
        direct:     label and node size use only reads directly assigned to this taxon.
        both:       label shows '<direct> / <cumulative>'; size uses cumulative.
        """
        if mode not in ("cumulative", "direct", "both"):
            return
        if mode != self._count_display:
            self._count_display = mode
            self._relayout()

    def set_cumulative_counts(self, enabled: bool) -> None:
        """Backward-compatible shim; prefer set_count_display()."""
        self.set_count_display("cumulative" if enabled else "direct")

    def set_label_visibility(self, mode: str) -> None:
        """Control which nodes display a name label.

        all:          all nodes (except zero-read pass-through internal nodes)
        direct:       only nodes with direct_reads > 0
        direct_major: nodes with direct_reads > 0, OR at a standard taxonomic
                      rank (genus/family/order/class/phylum/kingdom/domain)
        """
        if mode not in ("all", "direct", "direct_major"):
            return
        if mode != self._label_visibility:
            self._label_visibility = mode
            self._relayout()

    def _should_show_label(self, node: "_LayoutNode") -> bool:
        """Return True if this node should display a name label."""
        if node.is_virtual_node:
            return True  # always label virtual nodes
        vis = self._label_visibility
        if vis == "all":
            return not (node.reads == 0 and bool(node.children))
        if vis == "direct":
            return node.direct_reads > 0
        # "direct_major": show if direct reads OR a major structural rank
        if node.direct_reads > 0:
            return True
        rank = (node.data.get("rank") or "").lower()
        return rank in _MAJOR_RANKS

    def set_comparison_mode(
        self, enabled: bool, sample_a: Optional[str] = None, sample_b: Optional[str] = None
    ) -> None:
        """Enable/disable A-vs-B comparison overlay.

        When enabled, nodes are colored blue (A > B) or orange (B > A) based on
        per-sample read counts.  The existing color-by mode is suspended.
        """
        changed = (
            self._comparison_mode != enabled
            or self._compare_a != sample_a
            or self._compare_b != sample_b
        )
        self._comparison_mode = enabled
        self._compare_a = sample_a
        self._compare_b = sample_b
        if changed:
            self._relayout()

    @property
    def comparison_active(self) -> bool:
        return self._comparison_mode and bool(self._compare_a) and bool(self._compare_b)

    def set_review_marks(self, marks: Dict[str, str]) -> None:
        """Update the review status dict and redraw nodes."""
        self._review_marks = dict(marks)
        self._relayout()

    def set_branch_style(self, style: str) -> None:
        """Set branch drawing style: 'square' (Manhattan) or 'rounded' (S-curve bezier)."""
        if style != self._branch_style:
            self._branch_style = style
            self._relayout()

    def toggle_shortcuts_overlay(self) -> None:
        """Show/hide the shortcut hint overlay."""
        if self._shortcuts_overlay.isVisible():
            self._shortcuts_overlay.hide()
        else:
            self._shortcuts_overlay.reposition(self.viewport().size())
            self._shortcuts_overlay.show()
            self._shortcuts_overlay.raise_()

    def set_support_filter(self, sf: str) -> None:
        """Apply a multi-metric support filter and rebuild the tree.

        sf must be a key in _SUPPORT_LEVELS.  Nodes (leaves) that fail
        the threshold are removed; the tree is fully rebuilt so counts
        and layouts reflect only the passing taxa.
        """
        if sf != getattr(self, "_support_filter", "all"):
            self._support_filter = sf
            self._rebuild_scene()

    # ─── Internal scene building ────────────────────────────────────────────────

    def _build_scene(self, tree: dict) -> None:
        self._tree_data = tree
        self._rebuild_scene()

    def _rebuild_scene(self) -> None:
        if self._anim_timeline is not None:
            self._anim_timeline.stop()
            self._anim_timeline = None
        self._scene.clear()
        self._node_items.clear()
        self._placed_label_rects = []
        if not hasattr(self, "_tree_data") or self._tree_data is None:
            return

        self._root = _build_layout_tree(
            self._tree_data, depth=0, parent=None,
            selected_samples=self._selected_samples,
            focus_taxids=self._focus_taxids,
            support_filter=getattr(self, "_support_filter", "all"),
        )
        if self._root is None:
            return

        self._max_reads = max(1, self._root.max_node_reads())
        self._max_direct_reads = max(1, self._root.max_direct_reads())
        _apply_smart_collapse(self._root)

        # Compute effective layout constants from user-applied stretch scales
        self._eff_h_unit = _H_UNIT * self._h_scale
        self._eff_v_spacing = _V_SPACING * self._v_scale

        # Layout
        _assign_x(self._root, 0.0)
        self._apply_y_layout(self._root)

        # Draw — disable BSP indexing during bulk insert for faster addItem()
        self._scene.setItemIndexMethod(QGraphicsScene.NoIndex)
        self._draw_counter = 0
        self._draw_subtree(self._root, None)
        self._scene.setItemIndexMethod(QGraphicsScene.BspTreeIndex)

        # Fit
        self.fitInView(self._scene.itemsBoundingRect().adjusted(-20, -20, 20, 20),
                       Qt.KeepAspectRatio)
        self._apply_search_highlights()
        if hasattr(self, "_overlay_panel"):
            self._overlay_panel.raise_()

    def _draw_subtree(self, node: _LayoutNode,
                      parent_bar: Optional[float] = None,
                      parent_px: float = 0.0, parent_py: float = 0.0,
                      parent_r: float = 0.0) -> None:
        """Draw node and descendants.

        Square mode (default): parent_bar is the crossbar coordinate.
        Rounded mode: parent_px/py/r are the parent node's scene position and radius,
        used to draw S-curve bezier branches directly parent→child.
        """
        self._draw_counter = getattr(self, '_draw_counter', 0) + 1
        hu = self._eff_h_unit
        vs = self._eff_v_spacing
        if self._ltr_layout:
            px = node.y              # depth → X axis (grows rightward)
            py = node.x * hu         # leaf spread → Y axis (grows downward)
        else:
            px = node.x * hu         # leaf spread → X axis
            py = node.y              # depth → Y axis (grows downward)

        if self.comparison_active:
            color = _comparison_color(node, self._compare_a, self._compare_b)
        else:
            color = _metric_color(node, self._color_by, self._sample_colors)
        if self._count_display == "direct":
            radius = _node_radius(node.direct_reads, self._max_direct_reads, self._scale_mode)
        else:
            radius = _node_radius(node.reads, self._max_reads, self._scale_mode)
        radius = max(_MIN_RADIUS, radius * getattr(self, "_node_size_factor", 1.0))
        pen = QPen(color, _BRANCH_WIDTH)

        vis = node.visible_children()
        rounded = self._branch_style == "rounded"

        if rounded:
            # S-curve bezier from parent node to this node (skips crossbar/stub)
            if parent_bar is not None or (parent_px != 0.0 or parent_py != 0.0):
                path = QPainterPath()
                if self._ltr_layout:
                    mid_x = (parent_px + px) / 2
                    path.moveTo(parent_px + parent_r, parent_py)
                    path.cubicTo(mid_x, parent_py, mid_x, py, px - radius, py)
                else:
                    mid_y = (parent_py + py) / 2
                    path.moveTo(parent_px, parent_py + parent_r)
                    path.cubicTo(parent_px, mid_y, px, mid_y, px, py - radius)
                pi = QGraphicsPathItem(path)
                pi.setPen(pen)
                self._scene.addItem(pi)
            for c in vis:
                self._draw_subtree(c, None, px, py, radius)
        else:
            # Square (Manhattan) branches with crossbar
            if self._ltr_layout:
                if parent_bar is not None:
                    stub = QGraphicsLineItem(parent_bar, py, px - radius, py)
                    stub.setPen(pen)
                    self._scene.addItem(stub)
                if vis:
                    bar_x = px + vs * 0.5
                    stem = QGraphicsLineItem(px + radius, py, bar_x, py)
                    stem.setPen(pen)
                    self._scene.addItem(stem)
                    top_y    = vis[0].x  * hu
                    bottom_y = vis[-1].x * hu
                    cross = QGraphicsLineItem(bar_x, top_y, bar_x, bottom_y)
                    cross.setPen(pen)
                    self._scene.addItem(cross)
                    for c in vis:
                        self._draw_subtree(c, bar_x)
            else:
                if parent_bar is not None:
                    stub = QGraphicsLineItem(px, parent_bar, px, py - radius)
                    stub.setPen(pen)
                    self._scene.addItem(stub)
                if vis:
                    bar_y = py + vs * 0.5
                    stem = QGraphicsLineItem(px, py + radius, px, bar_y)
                    stem.setPen(pen)
                    self._scene.addItem(stem)
                    left_x  = vis[0].x  * hu
                    right_x = vis[-1].x * hu
                    bar = QGraphicsLineItem(left_x, bar_y, right_x, bar_y)
                    bar.setPen(pen)
                    self._scene.addItem(bar)
                    for c in vis:
                        self._draw_subtree(c, bar_y)

        # Node circle drawn on top of branch lines
        mark = self._review_marks.get(node.taxid, "")
        is_dark = node.taxid in self._dark_taxids
        node_item = _NodeItem(node, radius, color, self._sample_colors, review_mark=mark,
                              display_mode=self._node_display_mode, ltr_layout=self._ltr_layout,
                              dark_taxon=is_dark, color_by=self._color_by)
        node_item.setPos(px, py)
        node_item.setToolTip(self._build_tooltip(node))
        self._scene.addItem(node_item)
        self._node_items[node.taxid] = node_item

        # Collapse triangle
        if node.collapsed and node.children:
            tri = _CollapseTriangle(radius, color, pointing_right=self._ltr_layout)
            tri.setPos(px, py)
            self._scene.addItem(tri)

        if self._should_show_label(node):
            label_text = self._make_label(node)
            fm = QFontMetrics(self._label_font)
            lbl_w = fm.horizontalAdvance(label_text) + 6
            lbl_h = fm.height()
            gap = 3  # px gap between node edge and label
            # candidate positions: (scene_x, scene_y, offset_x_from_node, offset_y_from_node)
            candidates = [
                (px - lbl_w / 2, py - radius - lbl_h - gap,
                 -lbl_w / 2, -radius - lbl_h - gap),                         # above (primary)
                (px - lbl_w / 2, py + radius + gap,
                 -lbl_w / 2, radius + gap),                                   # below
                (px + radius + _LABEL_OFFSET, py - lbl_h / 2,
                 radius + _LABEL_OFFSET, -lbl_h / 2),                         # right
                (px - radius - lbl_w - _LABEL_OFFSET, py - lbl_h / 2,
                 -radius - lbl_w - _LABEL_OFFSET, -lbl_h / 2),                # left
            ]
            best_off_x = candidates[0][2]
            best_off_y = candidates[0][3]
            best_rect = QRectF(candidates[0][0], candidates[0][1], lbl_w, lbl_h)
            best_overlap = float('inf')
            for sx, sy, off_x, off_y in candidates:
                rect = QRectF(sx, sy, lbl_w, lbl_h)
                overlap = sum(
                    rect.intersected(r).width() * rect.intersected(r).height()
                    for r in self._placed_label_rects if rect.intersects(r)
                )
                if overlap < best_overlap:
                    best_overlap = overlap
                    best_off_x, best_off_y = off_x, off_y
                    best_rect = rect
                if overlap == 0:
                    break  # perfect position found
            node_item.set_label(
                label_text, self._label_font,
                color if node.reads > 0 else C_TREE_DIM,
                best_off_x, best_off_y,
            )
            self._placed_label_rects.append(best_rect)

    def _passes_support_filter(self, node: _LayoutNode) -> bool:
        if node.is_virtual_node:
            return True  # always show the unclassified branch
        sf = getattr(self, "_support_filter", "all")
        reads = node.reads
        flags = node.flags
        if sf == "all":
            return True
        if sf == "confident":
            return reads >= 10 and "low_read_support" not in flags
        if sf == "high":
            return reads >= 100 and "low_read_support" not in flags
        if sf == "flagged":
            return len(flags) > 0
        return True

    _TIER_LABELS = {
        1: "★★★ Tier 1 — High-confidence ancient DNA with strong damage signal",
        2: "★★☆ Tier 2 — Authentic ancient with moderate damage",
        3: "★☆☆ Tier 3 — Probable ancient, low or no damage detected",
        4: "◎  Tier 4 — Modern-like or uncertain authenticity",
        5: "○  Tier 5 — Low support / insufficient evidence",
    }
    _FLAG_DESCRIPTIONS = {
        "low_read_support": "fewer than 10 reads — treat as tentative",
        "low_complexity": "sequence may be repetitive or low-complexity",
        "high_blank": "elevated blank/negative control signal",
        "stack_concentration": "reads pile on one locus — possible contamination or mis-assembly",
        "conflicted": "significant fraction of reads assigned ambiguously",
        "missing_damage": "no deamination signal despite being classified as ancient",
    }

    @staticmethod
    def _build_tooltip(node: _LayoutNode) -> str:
        if node.is_virtual_node:
            return _build_virtual_tooltip(node)
        m = node.data.get("metrics") or {}
        by_sample = m.get("by_sample") or {}
        dmg = float(m.get("mean_damage_score") or 0)
        if dmg == 0 and by_sample:
            dmg_vals = [float((sm or {}).get("mean_damage_score") or 0) for sm in by_sample.values()]
            dmg = sum(dmg_vals) / len(dmg_vals) if dmg_vals else 0
        wt = float(m.get("weighted_reads") or 0)
        post = float(m.get("mean_posterior") or 0)
        blank = float(m.get("blank_fraction") or 0)
        breadth = float(m.get("best_reference_breadth") or 0)
        mean_breadth = float(m.get("mean_breadth") or 0)
        max_stack = int(m.get("max_stack_depth") or 0)
        conf = int(m.get("conflicted_reads") or 0)
        stk = float(m.get("stack_concentration") or m.get("top_locus_fraction") or 0)
        unique_refs = int(m.get("unique_refs") or 0)
        rpm = float(m.get("rpm") or 0)
        cumul = node.all_reads()

        DIM = "#7890a0"
        FG  = "#e6edf3"
        WARN = "#e08028"

        def _row(label: str, val: str, color: str = FG) -> str:
            return (f'<tr><td style="color:{DIM}; padding-right:10px;">{label}</td>'
                    f'<td style="color:{color};">{val}</td></tr>')

        auth_badge = node.authenticity_badge
        auth_tier = node.authenticity_tier

        # Core classification metrics
        reads_str = f"{node.reads:,} direct"
        if cumul > node.reads:
            reads_str += f" / {cumul:,} cumul."
        rows = [
            _row("Rank",      node.data.get("rank") or "no rank"),
            _row("Taxid",     node.taxid),
            _row("Reads",     reads_str),
            _row("Weighted",  f"{wt:.1f}"),
            _row("Posterior", f"{post:.3f}"),
        ]

        # Authenticity tier with full description
        if auth_tier is not None:
            tier_label = TaxonGraphicsTree._TIER_LABELS.get(auth_tier, f"Tier {auth_tier}")
            badge_str = f"{auth_badge}  {tier_label}" if auth_badge else tier_label
            rows.append(_row("Authenticity", badge_str))

        # Coverage / breadth
        if breadth > 0:
            rows.append(_row("Best breadth", f"{breadth:.3f}"))
        if mean_breadth > 0:
            rows.append(_row("Mean breadth", f"{mean_breadth:.3f}"))
        if max_stack > 0:
            rows.append(_row("Max stack",    str(max_stack)))
        if unique_refs > 0:
            rows.append(_row("Unique refs",  str(unique_refs)))

        # Damage & contamination signals
        if dmg > 0:
            dmg_color = "#90caf9" if dmg >= 0.10 else (WARN if dmg >= 0.03 else FG)
            rows.append(_row("Damage",       f"{dmg:.3f}", dmg_color))
        if blank > 0:
            bk_color = WARN if blank >= 0.10 else FG
            rows.append(_row("Blank frac.",  f"{blank:.4f}", bk_color))
        if stk > 0:
            rows.append(_row("Stack conc.",  f"{stk:.3f}"))
        if conf > 0:
            rows.append(_row("Conflicted",   str(conf)))
        if rpm > 0:
            rows.append(_row("RPM",          f"{rpm:.2f}"))

        # Per-sample breakdown with damage
        sample_rows = ""
        if by_sample:
            visible = [
                (sid,
                 int((sm or {}).get("hard_reads") or 0),
                 float((sm or {}).get("mean_damage_score") or 0))
                for sid, sm in sorted(by_sample.items())
                if int((sm or {}).get("hard_reads") or 0) > 0
            ]
            if visible:
                sample_rows = (
                    f'<tr><td colspan="2" style="color:{DIM}; padding-top:6px;">'
                    f'<i>Per sample:</i></td></tr>'
                )
                for sid, r, d in visible:
                    dmg_note = f"  dmg={d:.3f}" if d > 0 else ""
                    sample_rows += _row(sid, f"{r:,}{dmg_note}")

        # Flags with descriptions
        flags_row = ""
        if node.flags:
            flag_lines = []
            for f in node.flags:
                desc = TaxonGraphicsTree._FLAG_DESCRIPTIONS.get(f, f)
                flag_lines.append(f"⚑ {desc}")
            flags_row = (
                f'<tr><td colspan="2" style="color:{WARN}; padding-top:6px;">'
                + "<br>".join(flag_lines) + "</td></tr>"
            )

        return (
            f'<table style="font-family:Consolas,monospace; font-size:9pt; '
            f'white-space:nowrap;">'
            f'<tr><td colspan="2" style="color:#e6edf3; font-weight:bold; '
            f'font-size:10pt; padding-bottom:4px;">{node.name}</td></tr>'
            f'{"".join(rows)}{sample_rows}{flags_row}'
            f'</table>'
        )

    # ─── Mouse / keyboard ──────────────────────────────────────────────────────

    def wheelEvent(self, event: QWheelEvent) -> None:
        try:
            mods = event.modifiers()
            delta_y = event.angleDelta().y()
            factor = 1.15 if delta_y > 0 else 1.0 / 1.15
            if mods & Qt.ControlModifier and mods & Qt.ShiftModifier:
                # Ctrl+Shift+scroll → depth spacing (_v_scale)
                # In TTB layout v_scale is vertical; in LTR it is horizontal
                self._v_scale = max(0.15, min(8.0, self._v_scale * factor))
                self._scale_relayout()
                cur = Qt.SizeHorCursor if self._ltr_layout else Qt.SizeVerCursor
                self.viewport().setCursor(cur)
                self._cursor_reset_timer.start()
                event.accept()
            elif mods & Qt.ControlModifier:
                # Ctrl+scroll → uniform zoom (view transform, does not distort)
                cur = self.transform().m11()
                new_sx = cur * factor
                if new_sx < 0.01:
                    new_sx = 0.01
                elif new_sx > 200.0:
                    new_sx = 200.0
                if abs(new_sx - cur) > 1e-6:
                    self.scale(new_sx / cur, new_sx / cur)
                    # Debounce: relayout labels ~280ms after zoom settles
                    self._zoom_relabel_timer.start()
                event.accept()
            elif mods & Qt.ShiftModifier:
                # Shift+scroll → leaf spacing (_h_scale)
                # In TTB layout h_scale is horizontal; in LTR it is vertical
                self._h_scale = max(0.15, min(8.0, self._h_scale * factor))
                self._scale_relayout()
                cur = Qt.SizeVerCursor if self._ltr_layout else Qt.SizeHorCursor
                self.viewport().setCursor(cur)
                self._cursor_reset_timer.start()
                event.accept()
            else:
                super().wheelEvent(event)
        except Exception:
            super().wheelEvent(event)

    def mousePressEvent(self, event: QMouseEvent) -> None:
        if event.button() == Qt.MiddleButton:
            self._mid_pan_active = True
            self._mid_pan_start = event.pos()
            self.viewport().setCursor(Qt.ClosedHandCursor)
            event.accept()
            return

        if event.button() == Qt.LeftButton:
            self._left_pan_start = QPointF(event.pos())
            self._left_pan_dragging = False
            self._deselect_timer.stop()
            # Fire selection immediately on press.
            # Use items(rect) for robust hit-testing: handles label sub-items,
            # multi-select rings, and items with small bounding rects.
            ep = event.pos()
            hit_items = self.items(QRect(ep.x() - 5, ep.y() - 5, 11, 11))
            item = next((it for it in hit_items if isinstance(it, _NodeItem)), None)
            if item is None:
                # Walk parent chain for any hit subitem
                scene_pos = self.mapToScene(ep)
                top = self._scene.itemAt(scene_pos, QTransform())
                while top is not None and not isinstance(top, _NodeItem):
                    top = top.parentItem()
                item = top if isinstance(top, _NodeItem) else None
            if isinstance(item, _NodeItem):
                node = item.layout_node
                if event.modifiers() & Qt.ControlModifier:
                    self._toggle_multi_select(node)
                else:
                    self._select_node(node)
            else:
                # Defer deselect: cancelled if user starts a drag before 80ms
                self._deselect_timer.start()
                self.viewport().setCursor(Qt.OpenHandCursor)
            event.accept()
            return

        super().mousePressEvent(event)

    def mouseMoveEvent(self, event: QMouseEvent) -> None:
        if self._mid_pan_active:
            delta = event.pos() - self._mid_pan_start
            self._mid_pan_start = QPointF(event.pos())
            self.horizontalScrollBar().setValue(
                self.horizontalScrollBar().value() - int(delta.x())
            )
            self.verticalScrollBar().setValue(
                self.verticalScrollBar().value() - int(delta.y())
            )
            event.accept()
            return

        if event.buttons() & Qt.LeftButton and self._left_pan_start is not None:
            delta = event.pos() - self._left_pan_start
            if not self._left_pan_dragging and delta.manhattanLength() > 5:
                self._deselect_timer.stop()   # drag started — no deselect
                self._left_pan_dragging = True
                self.viewport().setCursor(Qt.ClosedHandCursor)
            if self._left_pan_dragging:
                self.horizontalScrollBar().setValue(
                    self.horizontalScrollBar().value() - int(delta.x())
                )
                self.verticalScrollBar().setValue(
                    self.verticalScrollBar().value() - int(delta.y())
                )
                self._left_pan_start = QPointF(event.pos())
                event.accept()
                return

        super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event: QMouseEvent) -> None:
        if event.button() == Qt.MiddleButton and self._mid_pan_active:
            self._mid_pan_active = False
            self.viewport().unsetCursor()
            event.accept()
            return
        if event.button() == Qt.LeftButton:
            self._left_pan_start = None
            self._left_pan_dragging = False
            self.viewport().unsetCursor()
            event.accept()
            return
        super().mouseReleaseEvent(event)

    def mouseDoubleClickEvent(self, event: QMouseEvent) -> None:
        if event.button() == Qt.LeftButton:
            scene_pos = self.mapToScene(event.pos())
            item = self._scene.itemAt(scene_pos, QTransform())
            while item is not None and not isinstance(item, _NodeItem):
                item = item.parentItem()
            if isinstance(item, _NodeItem):
                node = item.layout_node
                if node.children:
                    node.collapsed = not node.collapsed
                    self._relayout()
                event.accept()
                return
        super().mouseDoubleClickEvent(event)

    def resizeEvent(self, event) -> None:
        super().resizeEvent(event)
        self._reposition_overlay()
        if hasattr(self, '_shortcuts_overlay') and self._shortcuts_overlay.isVisible():
            self._shortcuts_overlay.reposition(self.viewport().size())

    def contextMenuEvent(self, event: QContextMenuEvent) -> None:
        scene_pos = self.mapToScene(event.pos())
        item = self._scene.itemAt(scene_pos, QTransform())
        while item is not None and not isinstance(item, _NodeItem):
            item = item.parentItem()

        menu = QMenu(self)
        if isinstance(item, _NodeItem):
            node = item.layout_node

            # ── Virtual (unclassified/excluded) nodes: simplified menu ────────
            if node.is_virtual_node:
                lbl = QAction(f"Category: {node.name}", self)
                lbl.setEnabled(False)
                menu.addAction(lbl)
                menu.addSeparator()
                if node.children:
                    if node.collapsed:
                        a = QAction("Expand", self)
                        a.triggered.connect(lambda: self._toggle_and_relayout(node, False))
                        menu.addAction(a)
                    else:
                        a = QAction("Collapse", self)
                        a.triggered.connect(lambda: self._toggle_and_relayout(node, True))
                        menu.addAction(a)
                    menu.addSeparator()
                if node.reads > 0:
                    a_ext = QAction(f"Extract reads ({node.reads:,})…", self)
                    a_ext.triggered.connect(lambda _c=False, tid=node.taxid: self._emit_extract([tid]))
                    menu.addAction(a_ext)
                a_copy = QAction("Copy category name", self)
                a_copy.triggered.connect(lambda _c=False, n=node.name: QApplication.clipboard().setText(n))
                menu.addAction(a_copy)
                menu.exec(event.globalPos())
                return
            # ── End virtual node fast-path ────────────────────────────────────

            if node.children:
                if node.collapsed:
                    act_expand = QAction(f"Expand '{node.name}'", self)
                    act_expand.triggered.connect(lambda: self._toggle_and_relayout(node, False))
                    menu.addAction(act_expand)
                    act_expand_sub = QAction("Expand subtree", self)
                    act_expand_sub.triggered.connect(
                        lambda: self._expand_subtree(node)
                    )
                    menu.addAction(act_expand_sub)
                else:
                    act_collapse = QAction(f"Collapse '{node.name}'", self)
                    act_collapse.triggered.connect(lambda: self._toggle_and_relayout(node, True))
                    menu.addAction(act_collapse)
                    act_collapse_sub = QAction("Collapse subtree", self)
                    act_collapse_sub.triggered.connect(
                        lambda _c=False, n=node: self._collapse_subtree(n)
                    )
                    menu.addAction(act_collapse_sub)
                menu.addSeparator()

            act_extract = QAction(f"Extract reads for '{node.name}'…", self)
            act_extract.triggered.connect(lambda: self._emit_extract([node.taxid]))
            menu.addAction(act_extract)
            act_extract_sub = QAction("Extract reads for subtree…", self)
            act_extract_sub.triggered.connect(
                lambda: self._emit_extract(self._collect_subtree_taxids(node))
            )
            menu.addAction(act_extract_sub)
            menu.addSeparator()

            act_copy_name = QAction(f"Copy name", self)
            act_copy_name.triggered.connect(
                lambda: QApplication.clipboard().setText(node.name)
            )
            menu.addAction(act_copy_name)
            act_copy_id = QAction("Copy taxid", self)
            act_copy_id.triggered.connect(
                lambda: QApplication.clipboard().setText(node.taxid)
            )
            menu.addAction(act_copy_id)
            act_ncbi = QAction("View on NCBI Taxonomy…", self)
            act_ncbi.triggered.connect(
                lambda: QDesktopServices.openUrl(
                    QUrl(f"https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={node.taxid}")
                )
            )
            menu.addAction(act_ncbi)
            act_copy_lineage = QAction("Copy full lineage", self)
            act_copy_lineage.triggered.connect(
                lambda _c=False, n=node: self._copy_lineage(n)
            )
            menu.addAction(act_copy_lineage)
            act_show_list = QAction("Show in List view", self)
            act_show_list.triggered.connect(
                lambda _c=False, tid=node.taxid: self.show_in_list_requested.emit(tid)
            )
            menu.addAction(act_show_list)
            act_bookmark = QAction("Bookmark this taxon", self)
            act_bookmark.triggered.connect(
                lambda _c=False, n=node: self.bookmark_requested.emit(n.taxid, n.name)
            )
            menu.addAction(act_bookmark)
            act_rerun = QAction("Run new analysis…", self)
            act_rerun.triggered.connect(lambda: self.rerun_requested.emit())
            menu.addAction(act_rerun)
            act_reclassify = QAction(f"Re-classify subtree of '{node.name}'…", self)
            act_reclassify.triggered.connect(
                lambda _c=False, n=node: self.reclassify_requested.emit(
                    list(_subtree_taxids(n.data)), n.name
                )
            )
            menu.addAction(act_reclassify)
            menu.addSeparator()

            # Review marks
            current_mark = self._review_marks.get(node.taxid, "")
            mark_menu = QMenu("Mark as…", self)
            for mark_key, mark_label in [
                ("reviewed",    "✓ Reviewed"),
                ("uncertain",   "? Uncertain"),
                ("rejected",    "✗ Rejected"),
                ("contaminant", "C Contaminant"),
            ]:
                act_m = QAction(mark_label, self)
                act_m.setCheckable(True)
                act_m.setChecked(current_mark == mark_key)
                act_m.triggered.connect(
                    lambda _chk, mk=mark_key, tid=node.taxid: self._set_review_mark(tid, mk)
                )
                mark_menu.addAction(act_m)
            if current_mark:
                mark_menu.addSeparator()
                act_clear = QAction("Clear mark", self)
                act_clear.triggered.connect(
                    lambda _chk=False, tid=node.taxid: self._set_review_mark(tid, "")
                )
                mark_menu.addAction(act_clear)
            menu.addMenu(mark_menu)
            menu.addSeparator()

            # ── Select submenu ────────────────────────────────────────────────
            sel_menu = QMenu("Select…", self)

            is_selected = node.taxid in self._multi_selected_taxids
            act_sel_toggle = QAction(
                "Remove from selection" if is_selected else "Add to selection (Ctrl+click)", self
            )
            act_sel_toggle.triggered.connect(
                lambda _c=False, n=node: self._toggle_multi_select(n)
            )
            sel_menu.addAction(act_sel_toggle)

            act_sel_only = QAction("Select only this node", self)
            act_sel_only.triggered.connect(
                lambda _c=False, tid=node.taxid: self.set_multi_selection([tid])
            )
            sel_menu.addAction(act_sel_only)

            sel_menu.addSeparator()

            act_sel_sub = QAction("Select subtree", self)
            act_sel_sub.triggered.connect(
                lambda _c=False, n=node: self.set_multi_selection(
                    self._collect_subtree_taxids(n)
                )
            )
            sel_menu.addAction(act_sel_sub)

            act_add_sub = QAction("Add subtree to selection", self)
            act_add_sub.triggered.connect(
                lambda _c=False, n=node: self.set_multi_selection(
                    list(set(list(self._multi_selected_taxids) + self._collect_subtree_taxids(n)))
                )
            )
            sel_menu.addAction(act_add_sub)

            act_sel_leaves = QAction("Select leaves in subtree", self)
            act_sel_leaves.triggered.connect(
                lambda _c=False, n=node: self.set_multi_selection(
                    [t for t in self._collect_subtree_taxids(n)
                     if t in self._node_items and not self._node_items[t].layout_node.children]
                )
            )
            sel_menu.addAction(act_sel_leaves)

            sel_menu.addSeparator()

            act_sel_sibs = QAction("Select siblings", self)
            def _do_sel_siblings(n=node):
                parent = n.parent
                if parent:
                    sibs = [c.taxid for c in parent.children if c.taxid != n.taxid]
                    self.set_multi_selection(sibs)
            act_sel_sibs.triggered.connect(lambda _c=False, fn=_do_sel_siblings: fn())
            sel_menu.addAction(act_sel_sibs)

            act_sel_sib_incl = QAction("Select siblings + this node", self)
            def _do_sel_siblings_incl(n=node):
                parent = n.parent
                if parent:
                    sibs = [c.taxid for c in parent.children]
                    self.set_multi_selection(sibs)
            act_sel_sib_incl.triggered.connect(lambda _c=False, fn=_do_sel_siblings_incl: fn())
            sel_menu.addAction(act_sel_sib_incl)

            sel_menu.addSeparator()

            # Global select operations available from node context
            act_sel_all = QAction("Select all nodes", self)
            act_sel_all.triggered.connect(
                lambda: self.set_multi_selection(list(self._node_items.keys()))
            )
            sel_menu.addAction(act_sel_all)

            act_sel_all_leaves = QAction("Select all leaves", self)
            act_sel_all_leaves.triggered.connect(
                lambda: self.set_multi_selection(
                    [tid for tid, itm in self._node_items.items()
                     if not itm.layout_node.children]
                )
            )
            sel_menu.addAction(act_sel_all_leaves)

            act_sel_reads = QAction("Select nodes with direct reads", self)
            act_sel_reads.triggered.connect(
                lambda: self.set_multi_selection(
                    [tid for tid, itm in self._node_items.items()
                     if itm.layout_node.direct_reads > 0]
                )
            )
            sel_menu.addAction(act_sel_reads)

            act_sel_inv = QAction("Invert selection", self)
            act_sel_inv.triggered.connect(
                lambda: self.set_multi_selection(
                    [tid for tid in self._node_items if tid not in self._multi_selected_taxids]
                )
            )
            sel_menu.addAction(act_sel_inv)

            menu.addMenu(sel_menu)
            menu.addSeparator()

        # Multi-selection actions
        if self._multi_selected_taxids:
            n_sel = len(self._multi_selected_taxids)
            act_extract_sel = QAction(f"Extract reads for {n_sel} selected taxa…", self)
            act_extract_sel.triggered.connect(
                lambda: self._emit_extract(list(self._multi_selected_taxids))
            )
            menu.addAction(act_extract_sel)

            act_focus = QAction(f"Focus: show subtree of {n_sel} selected taxa", self)
            act_focus.triggered.connect(self._focus_selected_subtree)
            menu.addAction(act_focus)

            act_clear_sel = QAction("Clear selection (Ctrl+click to add/remove)", self)
            act_clear_sel.triggered.connect(self._clear_multi_selection)
            menu.addAction(act_clear_sel)
            menu.addSeparator()

        if self._focus_taxids is not None:
            act_unfocus = QAction("Return to full tree", self)
            act_unfocus.triggered.connect(self.clear_focus)
            menu.addAction(act_unfocus)
            menu.addSeparator()

        act_fit = QAction("Fit to window", self)
        act_fit.triggered.connect(self._fit_all)
        menu.addAction(act_fit)
        act_expand_all = QAction("Expand all", self)
        act_expand_all.triggered.connect(self.expand_all)
        menu.addAction(act_expand_all)
        act_collapse_all = QAction("Collapse all", self)
        act_collapse_all.triggered.connect(self.collapse_all)
        menu.addAction(act_collapse_all)
        menu.exec(event.globalPos())

    # ─── Helpers ───────────────────────────────────────────────────────────────

    def _copy_lineage(self, node: _LayoutNode) -> None:
        """Walk the parent chain to build a semicolon-separated lineage string and copy it."""
        parts: List[str] = []
        cur: Optional[_LayoutNode] = node
        while cur is not None:
            parts.append(cur.name)
            cur = cur.parent
        parts.reverse()
        QApplication.clipboard().setText("; ".join(parts))

    def _make_label(self, node: _LayoutNode) -> str:
        """Build the node label string according to the current label mode."""
        if node.is_virtual_node:
            reads = node.reads
            if reads > 0:
                return f"[{node.name}  {reads:,}]"
            return f"[{node.name}]"
        if self.comparison_active:
            ra = node.by_sample.get(self._compare_a, 0)
            rb = node.by_sample.get(self._compare_b, 0)
            if ra > 0 or rb > 0:
                return f"{node.name}  (A:{ra:,}  B:{rb:,})"
            return node.name
        m = node.data.get("metrics") or {}
        if self._label_mode == "name":
            return node.name
        if self._label_mode == "name_damage":
            dmg = float(m.get("mean_damage_score") or 0)
            if dmg == 0:
                by_s = m.get("by_sample") or {}
                vals = [float((sm or {}).get("mean_damage_score") or 0) for sm in by_s.values()]
                dmg = sum(vals) / len(vals) if vals else 0
            if node.reads > 0 and dmg > 0:
                return f"{node.name}  ({dmg:.2f})"
            return node.name
        if self._label_mode == "name_tier":
            badge = node.authenticity_badge
            if node.reads > 0 and badge:
                return f"{node.name}  {badge}"
            return node.name
        # default: "name_reads"
        # node.reads = cumulative_hard_reads (reads at node + all descendants)
        # node.direct_reads = direct_hard_reads (reads assigned only to this node)
        cumul = node.reads
        direct = node.direct_reads
        mode = self._count_display
        if mode == "both":
            if cumul > 0:
                if direct > 0 and direct != cumul:
                    return f"{node.name}  ({direct:,} / {cumul:,})"
                return f"{node.name}  ({cumul:,})"
        elif mode == "direct":
            if direct > 0:
                return f"{node.name}  ({direct:,})"
        else:  # "cumulative" (default)
            if cumul > 0:
                return f"{node.name}  ({cumul:,})"
        return node.name

    def _set_review_mark(self, taxid: str, mark: str) -> None:
        if mark:
            self._review_marks[taxid] = mark
        else:
            self._review_marks.pop(taxid, None)
        self.mark_changed.emit(taxid, mark)
        self._relayout()

    def _select_node(self, node: _LayoutNode) -> None:
        # Regular click: clear multi-selection, select single node
        self._clear_multi_selection(emit=False)
        if self._selected_taxid and self._selected_taxid in self._node_items:
            self._node_items[self._selected_taxid]._selected_override = False
            self._node_items[self._selected_taxid].update()
        self._selected_taxid = node.taxid
        item = self._node_items.get(node.taxid)
        if item:
            item._selected_override = True
            item.update()
        self.node_selected.emit(node.data)

    def _toggle_multi_select(self, node: _LayoutNode) -> None:
        """Ctrl+Click: add or remove a node from the multi-selection set."""
        tid = node.taxid
        if tid in self._multi_selected_taxids:
            self._multi_selected_taxids.discard(tid)
            item = self._node_items.get(tid)
            if item:
                item.set_multi_selected(False)
        else:
            self._multi_selected_taxids.add(tid)
            item = self._node_items.get(tid)
            if item:
                item.set_multi_selected(True)
        # Include primary single-click selection in emitted count
        all_sel = set(self._multi_selected_taxids)
        if self._selected_taxid:
            all_sel.add(self._selected_taxid)
        self.multi_selection_changed.emit(list(all_sel))

    def _clear_multi_selection(self, emit: bool = True) -> None:
        for tid in list(self._multi_selected_taxids):
            item = self._node_items.get(tid)
            if item:
                item.set_multi_selected(False)
        self._multi_selected_taxids.clear()
        if emit:
            self.multi_selection_changed.emit([])

    def set_multi_selection(self, taxids: list) -> None:
        """Programmatically replace multi-selection with the given taxid list."""
        self._clear_multi_selection(emit=False)
        for tid in taxids:
            s = str(tid)
            self._multi_selected_taxids.add(s)
            item = self._node_items.get(s)
            if item:
                item.set_multi_selected(True)
        self.multi_selection_changed.emit(list(self._multi_selected_taxids))

    def _deselect_all(self) -> None:
        """Clear both primary and multi-selection; emit selection_cleared."""
        self._clear_multi_selection(emit=False)
        if self._selected_taxid and self._selected_taxid in self._node_items:
            self._node_items[self._selected_taxid]._selected_override = False
            self._node_items[self._selected_taxid].update()
        self._selected_taxid = None
        self.multi_selection_changed.emit([])
        self.selection_cleared.emit()

    def _focus_selected_subtree(self) -> None:
        """Filter the tree to show only selected nodes + all their descendants."""
        if not self._multi_selected_taxids or self._tree_data is None:
            return
        # Collect all descendants of each selected node
        all_subtree_ids: Set[str] = set()
        def _collect(data: dict, inside: bool) -> bool:
            tid = str(data.get("taxid") or "")
            is_focused = inside or tid in self._multi_selected_taxids
            has_focus_below = is_focused
            for child in data.get("children") or []:
                has_focus_below = _collect(child, is_focused) or has_focus_below
            if has_focus_below:
                all_subtree_ids.add(tid)
            return has_focus_below
        _collect(self._tree_data, False)

        self._focus_taxids = all_subtree_ids
        self._clear_multi_selection(emit=False)
        self._rebuild_scene()
        self.subtree_focus_changed.emit(True)

    def clear_focus(self) -> None:
        """Return to full tree view, clearing any subtree focus."""
        if self._focus_taxids is not None:
            self._focus_taxids = None
            self._rebuild_scene()
            self.subtree_focus_changed.emit(False)

    @property
    def is_focused(self) -> bool:
        """True when a subtree focus filter is active."""
        return self._focus_taxids is not None

    def _emit_extract(self, taxids: List[str]) -> None:
        self.extract_requested.emit(taxids)

    @staticmethod
    def _collect_subtree_taxids(node: _LayoutNode) -> List[str]:
        result: List[str] = [node.taxid]
        for c in node.children:
            result.extend(TaxonGraphicsTree._collect_subtree_taxids(c))
        return result

    def _toggle_and_relayout(self, node: "_LayoutNode", collapsed: bool) -> None:
        node.collapsed = collapsed
        self._relayout()

    def _expand_subtree(self, node: "_LayoutNode") -> None:
        self._set_collapsed(node, False)
        self._relayout()

    def _collapse_subtree(self, node: "_LayoutNode") -> None:
        """Collapse all children/descendants of node, leaving node itself expanded."""
        for child in node.children:
            self._set_collapsed(child, True)
        self._relayout()

    def _fit_all(self) -> None:
        r = self._scene.itemsBoundingRect()
        if r.isEmpty():
            return
        self.fitInView(r.adjusted(-20, -20, 20, 20), Qt.KeepAspectRatio)

    def _auto_compact_fit(self) -> None:
        """Adjust h/v spacing so the tree fills the viewport compactly, then fitInView.

        Always resets h_scale/v_scale to 1.0 first so repeated calls (e.g. Expand All)
        compute from the baseline and never compound-squish the tree.
        """
        if self._root is None:
            return
        if self._anim_timeline is not None:
            self._anim_timeline.stop()
            self._anim_timeline = None

        # Step 1 — draw at default scale to get an accurate bounding rect
        self._h_scale = 1.0
        self._v_scale = 1.0
        self._eff_h_unit = _H_UNIT
        self._eff_v_spacing = _V_SPACING
        _assign_x(self._root, 0.0)
        self._apply_y_layout(self._root)
        self._scene.clear()
        self._node_items.clear()
        self._placed_label_rects = []
        self._draw_counter = 0
        self._draw_subtree(self._root, None)

        vp = self.viewport().size()
        r = self._scene.itemsBoundingRect()
        if r.isEmpty() or vp.width() <= 0 or vp.height() <= 0:
            self._fit_all()
            return

        # Guard: if the viewport is suspiciously small the widget hasn't been
        # laid out yet (e.g. loaded while hidden in a QStackedWidget). Retry
        # after the event loop has processed the show event.
        if vp.width() < 50 or vp.height() < 50:
            from qtpy.QtCore import QTimer
            QTimer.singleShot(300, self._auto_compact_fit)
            return

        pad = 40
        sx = vp.width()  / (r.width()  + pad)
        sy = vp.height() / (r.height() + pad)

        # Step 2 — scale axes to fill the viewport.
        # h_scale fills the viewport width; v_scale compresses/stretches the depth.
        # v_scale is floored at a minimum that keeps labels readable (≥22 px per node).
        _MIN_V_PX = 22.0
        min_v_scale = _MIN_V_PX / _V_SPACING
        new_h = max(0.20, min(sx, 2.0))
        new_v = max(min_v_scale, min(sy, 1.5))

        if abs(new_h - 1.0) > 0.03 or abs(new_v - 1.0) > 0.03:
            self._h_scale = new_h
            self._v_scale = new_v
            self._eff_h_unit = _H_UNIT * self._h_scale
            self._eff_v_spacing = _V_SPACING * self._v_scale
            _assign_x(self._root, 0.0)
            self._apply_y_layout(self._root)
            self._scene.clear()
            self._node_items.clear()
            self._placed_label_rects = []
            self._draw_counter = 0
            self._draw_subtree(self._root, None)
            self._apply_search_highlights()

        r2 = self._scene.itemsBoundingRect()
        if r2.isEmpty():
            return
        self._scene.setSceneRect(r2.adjusted(-10, -10, 10, 10))
        scene_h = r2.height() + 20

        # Minimum scene-transform scale that keeps labels readable.
        # At 9pt base font, scale 0.70 → ~6.3pt effective — clearly legible.
        _MIN_ZOOM = 0.70

        if scene_h <= vp.height():
            # Tree fits in viewport height: use KeepAspectRatio so node circles
            # stay circular (IgnoreAspectRatio stretches x/y independently when the
            # h/v scale caps don't exactly match the viewport aspect ratio).
            self.fitInView(r2.adjusted(-10, -10, 10, 10), Qt.KeepAspectRatio)
            # Enforce minimum readable zoom — on wide/shallow trees fitInView can
            # return a scale so small that labels are invisible.
            if self.transform().m11() < _MIN_ZOOM:
                self.setTransform(QTransform.fromScale(_MIN_ZOOM, _MIN_ZOOM))
                self.verticalScrollBar().setValue(0)
                self.horizontalScrollBar().setValue(0)
        else:
            # Tree is taller than viewport at minimum readable spacing.
            # Reset to 1:1 (h_scale already fills width) and let scrollbars handle height.
            self.resetTransform()
            self.verticalScrollBar().setValue(0)

    # ─── Public helpers ────────────────────────────────────────────────────────

    def highlight_search(self, text: str) -> None:
        """Highlight tree nodes whose name contains `text` (case-insensitive)."""
        self._search_text = text.strip().lower()
        self._apply_search_highlights()

    def _apply_search_highlights(self) -> None:
        text = self._search_text
        for taxid, item in self._node_items.items():
            matched = bool(text) and text in item.layout_node.name.lower()
            item.set_search_match(matched)
            # Restore multi-selection visual state after scene rebuild
            item.set_multi_selected(taxid in self._multi_selected_taxids)

    def select_taxid(self, taxid: str) -> None:
        item = self._node_items.get(str(taxid))
        if item:
            self._select_node(item.layout_node)
            self.centerOn(item)

    def select_taxid_silent(self, taxid: str) -> None:
        """Update visual selection state without emitting node_selected.
        Used for cross-widget sync to prevent feedback loops.
        """
        item = self._node_items.get(str(taxid))
        if item is None:
            return
        node = item.layout_node
        self._clear_multi_selection(emit=False)
        if self._selected_taxid and self._selected_taxid in self._node_items:
            self._node_items[self._selected_taxid]._selected_override = False
            self._node_items[self._selected_taxid].update()
        self._selected_taxid = node.taxid
        item._selected_override = True
        item.update()
        self.centerOn(item)

    def get_selected_taxids(self) -> List[str]:
        """Return all selected taxids: multi-selection set + primary if active, else single."""
        if self._multi_selected_taxids:
            result = set(self._multi_selected_taxids)
            if self._selected_taxid:
                result.add(self._selected_taxid)
            return list(result)
        if self._selected_taxid:
            return [self._selected_taxid]
        return []

    def expand_to_taxid(self, taxid: str) -> None:
        """Expand the path from root to taxid and center the view on the node."""
        if self._root is None:
            return
        target_tid = str(taxid)
        changed = self._expand_path(self._root, target_tid)
        if changed:
            self._relayout()
        item = self._node_items.get(target_tid)
        if item:
            self.centerOn(item)

    def _expand_path(self, node: "_LayoutNode", target: str) -> bool:
        """Return True if target is in this subtree; expand the path to it."""
        if node.taxid == target:
            return True
        for c in node.children:
            if self._expand_path(c, target):
                node.collapsed = False
                return True
        return False

    def expand_all(self) -> None:
        if self._root is None:
            return
        self._set_collapsed(self._root, False)
        self._relayout()

    def collapse_all(self) -> None:
        if self._root is None:
            return
        self._set_collapsed(self._root, True)
        self._relayout()

    def collapse_to_rank(self, rank: str) -> None:
        """Collapse all nodes whose rank is below 'rank' in the taxonomy hierarchy."""
        if self._root is None:
            return
        try:
            cutoff = RANK_ORDER.index(rank)
        except ValueError:
            cutoff = len(RANK_ORDER)

        def _collapse(node: "_LayoutNode") -> None:
            if not node.children:
                return
            node_rank = node.data.get("rank") or ""
            try:
                node_rank_idx = RANK_ORDER.index(node_rank)
            except ValueError:
                node_rank_idx = -1
            # Collapse if this node's rank is at or below the cutoff rank
            node.collapsed = node_rank_idx >= cutoff
            for c in node.children:
                _collapse(c)
        _collapse(self._root)
        self._relayout()

    def collapse_to_clade(self, clade_taxid: Optional[int]) -> None:
        """Collapse everything outside the given clade's subtree.

        Nodes on the path from root → clade root are force-expanded.
        Nodes inside the clade subtree are left unchanged (rank filter applies).
        Pass None to remove the clade filter (expand all).
        """
        if self._root is None:
            return
        if clade_taxid is None:
            self._set_collapsed(self._root, False)
            self._relayout()
            return

        clade_str = str(clade_taxid)
        in_path: set = set()

        def _mark(node: "_LayoutNode") -> bool:
            this = str(node.data.get("taxid") or "")
            matched = this == clade_str
            child_hit = any(_mark(c) for c in node.children)
            if matched or child_hit:
                in_path.add(id(node))
            return matched or child_hit

        _mark(self._root)

        def _apply(node: "_LayoutNode", inside: bool = False) -> None:
            this = str(node.data.get("taxid") or "")
            if inside or this == clade_str:
                # Inside the target clade — don't change collapsed state;
                # rank filter already applied it.
                for c in node.children:
                    _apply(c, inside=True)
            elif id(node) in in_path:
                # Ancestor of the clade root — keep expanded.
                node.collapsed = False
                for c in node.children:
                    _apply(c)
            else:
                # Outside the clade — collapse this subtree.
                if node.children:
                    node.collapsed = True

        _apply(self._root)
        self._relayout()

    @staticmethod
    def _set_collapsed(node: "_LayoutNode", state: bool) -> None:
        if node.children:
            node.collapsed = state
            for c in node.children:
                TaxonGraphicsTree._set_collapsed(c, state)

    def _apply_y_layout(self, root: _LayoutNode) -> None:
        vs = self._eff_v_spacing
        if self._layout_mode == "dendrogram":
            md = _max_depth(root)
            _assign_y_dendrogram(root, md, vs)
            for c in root.children:
                _assign_y_dendrogram(c, md, vs)
        else:
            _assign_y(root, vs)

    def _relayout(self) -> None:
        """Recompute layout and redraw; preserves current zoom/pan.

        When the tree has ≤300 visible nodes, node circles slide from their
        old positions to their new ones over 220 ms (ease-in-out). Branch lines
        snap immediately to the new layout — the same approach used by MEGAN.
        """
        if self._root is None:
            return

        # Stop any running animation BEFORE clearing the scene.
        # The _tick closure holds live _NodeItem references; scene.clear() deletes the
        # underlying C++ objects, so any subsequent _tick call raises RuntimeError.
        if self._anim_timeline is not None:
            self._anim_timeline.stop()
            self._anim_timeline = None

        # Snapshot old node positions before clearing
        old_pos: Dict[str, QPointF] = {
            tid: item.pos() for tid, item in self._node_items.items()
        }

        self._eff_h_unit = _H_UNIT * self._h_scale
        self._eff_v_spacing = _V_SPACING * self._v_scale
        _assign_x(self._root, 0.0)
        self._apply_y_layout(self._root)
        self._scene.clear()
        self._node_items.clear()
        self._placed_label_rects = []
        self._draw_counter = 0
        self._draw_subtree(self._root, None)
        self._scene.setSceneRect(self._scene.itemsBoundingRect().adjusted(-20, -20, 20, 20))
        self._apply_search_highlights()

        if old_pos and len(self._node_items) <= 300:
            self._animate_nodes(old_pos)
        if hasattr(self, "_overlay_panel"):
            self._overlay_panel.raise_()
        self.viewport().update()

    def _animate_nodes(self, old_pos: Dict[str, QPointF]) -> None:
        """Slide node circles from old_pos to their current (new) positions."""
        transitions = []
        for tid, item in self._node_items.items():
            if tid in old_pos:
                frm = old_pos[tid]
                to  = item.pos()
                if (to - frm).manhattanLength() > 0.5:
                    transitions.append((item, frm, to))

        if not transitions:
            return

        # Stop any running animation
        if self._anim_timeline is not None:
            self._anim_timeline.stop()
            self._anim_timeline = None

        # Place all animated items at their start positions
        for item, frm, _to in transitions:
            item.setPos(frm)

        tl = QTimeLine(220, self)
        tl.setFrameRange(0, 100)
        tl.setCurveShape(QTimeLine.EaseInOutCurve)

        def _tick(frame: int) -> None:
            t = frame / 100.0
            for item, frm, to in transitions:
                item.setPos(
                    frm.x() + (to.x() - frm.x()) * t,
                    frm.y() + (to.y() - frm.y()) * t,
                )

        def _done() -> None:
            for item, _frm, to in transitions:
                item.setPos(to)
            self._anim_timeline = None

        tl.frameChanged.connect(_tick)
        tl.finished.connect(_done)
        tl.start()
        self._anim_timeline = tl
