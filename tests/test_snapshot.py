"""Incremental fact + SCD2 snapshot across two ChEMBL "releases" (F2.2).

Drives dbt against a throwaway DuckDB file twice: once on the v1 fixture and once on
the v2 delta (``tests/fixtures/build_raw_v2.py``). Proves three things the single-run
``dbt build`` in CI cannot show on its own:

  * ``fact_activity`` refreshes incrementally — a re-run with no new activities is a
    no-op, and the v2 delta inserts only the new activity_ids.
  * ``compound_status`` records a status change as a second SCD2 window (the old row is
    closed, a new open row is written) when CHEMBLM2 advances max_phase 2 -> 4.
  * No molecule ever has more than one open validity window.
"""

from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

import duckdb
import pytest

_REPO = Path(__file__).resolve().parents[1]
_RAW_V1 = _REPO / "tests" / "fixtures" / "raw"
_RAW_V2 = _REPO / "tests" / "fixtures" / "raw_v2"

pytestmark = pytest.mark.skipif(
    shutil.which("dbt") is None, reason="dbt CLI not installed (needs the `dbt` extra)"
)


def _dbt(*args: str, raw_dir: Path, db_path: Path) -> None:
    """Run a dbt command against a specific raw fixture dir and DuckDB file."""
    env = {"SARVAULT_RAW_DIR": str(raw_dir), "DUCKDB_PATH": str(db_path)}
    result = subprocess.run(
        ["dbt", *args, "--project-dir", "dbt", "--profiles-dir", "dbt/profiles"],
        cwd=_REPO,
        env={**_os_environ(), **env},
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise AssertionError(
            f"dbt {' '.join(args)} failed:\n{result.stdout[-3000:]}\n{result.stderr[-1500:]}"
        )


def _os_environ() -> dict:
    import os

    return dict(os.environ)


def _refresh(raw_dir: Path, db_path: Path) -> None:
    """Rebuild the fact-activity lineage + re-run the snapshot for a release."""
    _dbt("run", "--select", "+fact_activity", raw_dir=raw_dir, db_path=db_path)
    _dbt("snapshot", raw_dir=raw_dir, db_path=db_path)


def _fact_ids(db_path: Path) -> set[int]:
    con = duckdb.connect(str(db_path), read_only=True)
    try:
        return {r[0] for r in con.sql("select activity_id from main_marts.fact_activity").fetchall()}
    finally:
        con.close()


def _windows(db_path: Path, molecule: str) -> list[tuple]:
    con = duckdb.connect(str(db_path), read_only=True)
    try:
        return con.sql(
            "select max_phase, is_approved_drug, (dbt_valid_to is null) as is_open "
            "from main.compound_status where molecule_chembl_id = ? "
            "order by dbt_valid_from",
            params=[molecule],
        ).fetchall()
    finally:
        con.close()


def _max_open_windows(db_path: Path) -> int:
    con = duckdb.connect(str(db_path), read_only=True)
    try:
        row = con.sql(
            "select coalesce(max(n), 0) from ("
            "  select count(*) as n from main.compound_status "
            "  where dbt_valid_to is null group by molecule_chembl_id)"
        ).fetchone()
        return int(row[0])
    finally:
        con.close()


@dataclass
class ReleaseRun:
    """Captured state from the two-release build, plus the DuckDB path to query."""

    db: Path
    ids_v1: set[int]
    ids_v1_rerun: set[int]
    ids_v2: set[int]


@pytest.fixture(scope="module")
def two_release_db(tmp_path_factory: pytest.TempPathFactory) -> ReleaseRun:
    """Build v1, re-run v1 (no-op), then apply the v2 delta — all on one DuckDB file."""
    assert _RAW_V2.exists(), "run `python tests/fixtures/build_raw_v2.py` to generate raw_v2"
    db = tmp_path_factory.mktemp("scd2") / "warehouse.duckdb"

    _refresh(_RAW_V1, db)
    ids_v1 = _fact_ids(db)

    _refresh(_RAW_V1, db)  # idempotent re-run
    ids_v1_rerun = _fact_ids(db)

    _refresh(_RAW_V2, db)
    ids_v2 = _fact_ids(db)

    return ReleaseRun(db=db, ids_v1=ids_v1, ids_v1_rerun=ids_v1_rerun, ids_v2=ids_v2)


def test_incremental_rerun_is_a_noop(two_release_db: ReleaseRun) -> None:
    # Re-running the same release must not add, drop, or duplicate any fact rows.
    assert two_release_db.ids_v1_rerun == two_release_db.ids_v1


def test_incremental_delta_inserts_only_new_activities(two_release_db: ReleaseRun) -> None:
    added = two_release_db.ids_v2 - two_release_db.ids_v1
    removed = two_release_db.ids_v1 - two_release_db.ids_v2
    # Exactly the two new activity_ids appear; nothing is lost.
    assert added == {13, 14}
    assert removed == set()


def test_scd2_captures_status_change_as_two_windows(two_release_db: ReleaseRun) -> None:
    windows = _windows(two_release_db.db, "CHEMBLM2")
    # One closed research-phase window, then one open approved window.
    assert len(windows) == 2
    (mp_old, approved_old, open_old), (mp_new, approved_new, open_new) = windows
    assert (mp_old, approved_old, open_old) == (2.0, False, False)
    assert (mp_new, approved_new, open_new) == (4.0, True, True)


def test_unchanged_compound_keeps_single_window(two_release_db: ReleaseRun) -> None:
    # CHEMBLM1 never changes status, so it must still have exactly one open window.
    windows = _windows(two_release_db.db, "CHEMBLM1")
    assert len(windows) == 1
    assert windows[0][2] is True  # still open


def test_no_compound_has_multiple_open_windows(two_release_db: ReleaseRun) -> None:
    assert _max_open_windows(two_release_db.db) == 1
