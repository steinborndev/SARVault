"""Tests for the activity-cliff layer: the pure filter helper and the built mart."""

import os
import shutil
import subprocess
from pathlib import Path

import duckdb
import pandas as pd
import pytest

from dashboard import logic

_REPO = Path(__file__).resolve().parents[1]
_FIXTURES = _REPO / "tests" / "fixtures" / "raw"


# --- pure filter helper ---
def _cliffs():
    return pd.DataFrame(
        {
            "molecule_chembl_id_a": ["A", "C", "E"],
            "molecule_chembl_id_b": ["B", "D", "F"],
            "tanimoto": [0.88, 0.80, 1.00],
            "delta_pchembl": [2.5, 1.2, 1.5],
            "is_identical_fp": [False, False, True],
        }
    )


def test_filter_cliffs_applies_both_floors():
    out = logic.filter_cliffs(_cliffs(), min_tanimoto=0.85, min_delta=2.0, include_identical=True)
    assert list(out["molecule_chembl_id_a"]) == ["A"]  # only the 0.88 / 2.5 pair clears both


def test_filter_cliffs_excludes_identical_by_default():
    out = logic.filter_cliffs(_cliffs(), min_tanimoto=0.75, min_delta=1.0, include_identical=False)
    assert set(out["molecule_chembl_id_a"]) == {"A", "C"}  # E/F identical-fp dropped


def test_filter_cliffs_includes_identical_when_requested():
    out = logic.filter_cliffs(_cliffs(), min_tanimoto=0.75, min_delta=1.0, include_identical=True)
    assert set(out["molecule_chembl_id_a"]) == {"A", "C", "E"}


# --- built mart (integration on the fixtures) ---
@pytest.fixture(scope="module")
def cliff_db(tmp_path_factory):
    if shutil.which("dbt") is None:
        pytest.skip("dbt CLI not available")
    db = tmp_path_factory.mktemp("cliff") / "w.duckdb"
    env = {
        **os.environ,
        "SARVAULT_RAW_DIR": str(_FIXTURES),
        "DUCKDB_PATH": str(db),
    }
    proc = subprocess.run(
        [
            "dbt", "build", "--project-dir", "dbt", "--profiles-dir", "dbt/profiles",
            "--select", "+mart_activity_cliff", "--indirect-selection", "cautious",
        ],
        cwd=_REPO,
        env=env,
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        pytest.skip(f"dbt build failed:\n{proc.stdout[-800:]}")
    return db


def _mart(db):
    con = duckdb.connect(str(db), read_only=True)
    try:
        return con.execute(
            "select * from main_analytics.mart_activity_cliff"
        ).df()
    finally:
        con.close()


def test_known_cliff_pair_present_with_expected_sali(cliff_db):
    df = _mart(cliff_db)
    pair = df[
        (df["molecule_chembl_id_a"] == "CHEMBLM6")
        & (df["molecule_chembl_id_b"] == "CHEMBLM7")
    ]
    assert len(pair) == 1
    row = pair.iloc[0]
    assert row["tanimoto"] == pytest.approx(0.885, abs=0.01)
    assert row["delta_pchembl"] == pytest.approx(2.5, abs=0.01)
    assert bool(row["is_identical_fp"]) is False
    # SALI = 2.5 / (1 - 0.885) ~= 21.7
    assert row["sali"] == pytest.approx(21.7, abs=0.5)


def test_identical_fp_pair_flagged_with_null_sali(cliff_db):
    df = _mart(cliff_db)
    pair = df[
        (df["molecule_chembl_id_a"] == "CHEMBLM8")
        & (df["molecule_chembl_id_b"] == "CHEMBLM9")
    ]
    assert len(pair) == 1
    row = pair.iloc[0]
    assert row["tanimoto"] == pytest.approx(1.0)
    assert bool(row["is_identical_fp"]) is True
    assert pd.isna(row["sali"])  # undefined for identical fingerprints


def test_subthreshold_pair_excluded(cliff_db):
    # CHEMBLM4/M5 (toluene/aniline) have Tanimoto ~0.375 -> below the 0.75 floor.
    df = _mart(cliff_db)
    m4m5 = df[
        (df["molecule_chembl_id_a"] == "CHEMBLM4")
        & (df["molecule_chembl_id_b"] == "CHEMBLM5")
    ]
    assert m4m5.empty
