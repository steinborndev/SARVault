"""Tests for scripts/profile_sar.py (pure Tanimoto, cliff tallying, guard, report)."""

import pytest

from scripts import profile_sar as p


# --- integer Tanimoto ---
def test_tanimoto_identity_disjoint_partial():
    assert p.tanimoto(0xFF, 8, 0xFF, 8) == 1.0
    assert p.tanimoto(0xF0, 4, 0x0F, 4) == 0.0
    assert p.tanimoto(0xFF, 8, 0x0F, 4) == 0.5


# --- cliff grid + top-SALI ---
def _per_target():
    # One target, three compounds:
    #   A "ff" pot 8.0, B "ff" pot 5.0 (identical fp -> T=1.0), C "0f" pot 8.0 (T=0.5 vs A/B)
    return {
        (1, "T"): [
            (1, "A", "ff", 8.0),
            (2, "B", "ff", 5.0),
            (3, "C", "0f", 8.0),
        ]
    }


def test_profile_cliffs_grid_counts():
    grid, top, budget = p.profile_cliffs(
        _per_target(),
        tan_thresholds=(0.5, 0.9),
        delta_thresholds=(1.0, 2.0),
        top_k=10,
        max_per_target=0,
    )
    # Pairs: (A,B) T=1.0 Δ=3.0 ; (A,C) T=0.5 Δ=0 ; (B,C) T=0.5 Δ=3.0
    assert grid[(0.5, 1.0)] == 2  # (A,B) and (B,C)
    assert grid[(0.5, 2.0)] == 2
    assert grid[(0.9, 1.0)] == 1  # only (A,B) clears T>=0.9
    assert grid[(0.9, 2.0)] == 1


def test_profile_cliffs_pair_budget_and_top_sali():
    _, top, budget = p.profile_cliffs(
        _per_target(), (0.5,), (1.0,), top_k=10, max_per_target=0
    )
    # 3 compounds -> 3 pairs for the single target.
    assert budget == [("T", 3, 3)]
    # Sharpest cliff first: (A,B) identical fp -> infinite SALI; then (B,C) -> 6.0.
    assert top[0][1:] == ("T", "A", "B", 1.0, 3.0)
    assert top[0][0] == float("inf")
    assert top[1][1:] == ("T", "B", "C", 0.5, 3.0)


def test_profile_cliffs_respects_max_per_target():
    # Capping to the 2 most potent keeps A (8.0) and C (8.0), drops B (5.0).
    _, top, budget = p.profile_cliffs(
        _per_target(), (0.4,), (0.0,), top_k=10, max_per_target=2
    )
    assert budget == [("T", 2, 1)]  # capped to 2 compounds -> 1 pair


# --- cheminfo guard ---
class _FakeCon:
    def __init__(self, missing):
        self._missing = missing

    def execute(self, sql):
        if any(m in sql for m in self._missing):
            raise RuntimeError("no such table")
        return self

    def fetchone(self):
        return (0,)


def test_require_cheminfo_exits_when_marts_missing():
    with pytest.raises(SystemExit):
        p.require_cheminfo(_FakeCon(missing=["mart_compound_fingerprint"]))


def test_require_cheminfo_passes_when_present():
    # Nothing missing -> no exit.
    assert p.require_cheminfo(_FakeCon(missing=[])) is None


# --- report rendering ---
def test_render_report_has_all_sections():
    grid, top, budget = p.profile_cliffs(_per_target(), (0.5,), (1.0,), 10, 0)
    volume = {
        "compounds (dim_compound)": 3,
        "fingerprinted compounds": 3,
    }
    md = p.render_report(volume, [2, 1], 0, grid, top, budget, (0.5,), (1.0,))
    for heading in (
        "## 1. Volume & survival",
        "## 2. Scaffold distribution",
        "## 3. Per-target pair budget",
        "## 4. Activity-cliff matrix",
        "## 5. Top cliffs by SALI",
    ):
        assert heading in md
    assert "100.0%" in md  # fingerprint coverage 3/3
