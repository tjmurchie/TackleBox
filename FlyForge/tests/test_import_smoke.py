"""Smoke test: confirms FlyForge.py imports cleanly as a module (no argv
parsing or other side effects at import time)."""

import FlyForge as ff


def test_import_succeeds():
    assert hasattr(ff, "Bait")
    assert hasattr(ff, "compute_tm")
    assert hasattr(ff, "tile_sequence")
    assert hasattr(ff, "design_opool")
