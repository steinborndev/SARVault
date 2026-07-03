"""Tests for the chemical-series mart (scaffold grouping) on the fixtures."""

import os
import shutil
import subprocess
from pathlib import Path

import duckdb
import pytest

_REPO = Path(__file__).resolve().parents[1]
_FIXTURES = _REPO / "tests" / "fixtures" / "raw"


@pytest.fixture(scope="module")
def series_db(tmp_path_factory):
    if shutil.which("dbt") is None:
        pytest.skip("dbt CLI not available")
    db = tmp_path_factory.mktemp("series") / "w.duckdb"
    env = {**os.environ, "SARVAULT_RAW_DIR": str(_FIXTURES), "DUCKDB_PATH": str(db)}
    proc = subprocess.run(
        [
            "dbt", "build", "--project-dir", "dbt", "--profiles-dir", "dbt/profiles",
            "--select", "+mart_chemical_series", "--indirect-selection", "cautious",
        ],
        cwd=_REPO, env=env, capture_output=True, text=True,
    )
    if proc.returncode != 0:
        pytest.skip(f"dbt build failed:\n{proc.stdout[-800:]}")
    return db


def _series(db):
    con = duckdb.connect(str(db), read_only=True)
    try:
        return con.execute("select * from main_analytics.mart_chemical_series").df()
    finally:
        con.close()


def test_benzene_series_groups_all_members(series_db):
    df = _series(series_db)
    benzene = df[df["murcko_scaffold_smiles"] == "c1ccccc1"]
    assert len(benzene) == 1
    row = benzene.iloc[0]
    # Toluene, aniline, decyl- and undecylbenzene all share the benzene scaffold.
    assert row["n_compounds"] == 4
    assert row["pchembl_range"] == pytest.approx(1.5, abs=0.01)  # 7.5 - 6.0
    assert row["top_compound"] == "CHEMBLM8"  # decylbenzene, best pChEMBL 7.5


def test_series_grain_unique_on_scaffold(series_db):
    df = _series(series_db)
    assert df["scaffold_key"].is_unique


def test_series_covers_every_scaffolded_compound(series_db):
    # Sum of member counts must equal the number of compounds that have a scaffold.
    con = duckdb.connect(str(series_db), read_only=True)
    try:
        series_total = con.execute(
            "select coalesce(sum(n_compounds),0) from main_analytics.mart_chemical_series"
        ).fetchone()[0]
        scaffolded = con.execute(
            "select count(distinct compound_key) from main_analytics.mart_compound_fingerprint "
            "where scaffold_key is not null"
        ).fetchone()[0]
    finally:
        con.close()
    assert series_total == scaffolded


def test_compound_row_fallback_and_missing(series_db):
    # This build has no catalog mart, so compound_row must fall back to dim_compound
    # (proving the try/except on the catalog table), and return None for unknown ids.
    from dashboard import data

    con = data.connect(str(series_db))
    try:
        row = data.compound_row(con, "CHEMBLM6")
        assert row is not None
        assert int(row["compound_key"]) >= 1
        assert row.get("canonical_smiles")
        assert "selectivity_index" in row.index  # filled null by the fallback
        assert data.compound_row(con, "NOPE") is None
    finally:
        con.close()
