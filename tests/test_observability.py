"""Data-observability anomaly guards fire on bad data and stay quiet on good data (F2.3).

Builds the fact-activity lineage against a deliberately anomalous fixture
(``tests/fixtures/build_raw_bad.py``) and asserts every anomaly test in
``dbt/tests/assert_*.sql`` fails; then confirms the same tests pass on the good v1
fixture. Mirrors the two-release harness in test_snapshot.py.
"""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

_REPO = Path(__file__).resolve().parents[1]
_RAW_V1 = _REPO / "tests" / "fixtures" / "raw"

_ANOMALY_TESTS = (
    "assert_rowcount_floors",
    "assert_null_rate_ceiling",
    "assert_pchembl_mean_in_band",
)

pytestmark = pytest.mark.skipif(
    shutil.which("dbt") is None, reason="dbt CLI not installed (needs the `dbt` extra)"
)


def _os_environ() -> dict:
    import os

    return dict(os.environ)


def _dbt(*args: str, raw_dir: Path, db_path: Path, target_path: Path) -> subprocess.CompletedProcess:
    env = {"SARVAULT_RAW_DIR": str(raw_dir), "DUCKDB_PATH": str(db_path)}
    return subprocess.run(
        [
            "dbt", *args,
            "--project-dir", "dbt",
            "--profiles-dir", "dbt/profiles",
            "--target-path", str(target_path),
        ],
        cwd=_REPO,
        env={**_os_environ(), **env},
        capture_output=True,
        text=True,
    )


def _anomaly_statuses(raw_dir: Path, tmp: Path) -> dict[str, str]:
    """Build the fact lineage against raw_dir, run the anomaly tests, return {test: status}."""
    db_path = tmp / "warehouse.duckdb"
    target_path = tmp / "target"

    built = _dbt("run", "--select", "+fact_activity", raw_dir=raw_dir, db_path=db_path,
                 target_path=target_path)
    assert built.returncode == 0, f"dbt run failed:\n{built.stdout[-2000:]}\n{built.stderr[-800:]}"

    _dbt("test", "--select", *_ANOMALY_TESTS, raw_dir=raw_dir, db_path=db_path,
         target_path=target_path)  # exit code ignored; statuses come from run_results

    results = json.loads((target_path / "run_results.json").read_text())
    statuses = {}
    for r in results["results"]:
        name = r["unique_id"].split(".")[-1]
        if name in _ANOMALY_TESTS:
            statuses[name] = r["status"]
    return statuses


def test_anomaly_checks_fire_on_bad_data(tmp_path: Path) -> None:
    sys.path.insert(0, str(_REPO / "tests" / "fixtures"))
    from build_raw_bad import write_bad_fixture

    bad_raw = write_bad_fixture(tmp_path / "raw_bad")
    statuses = _anomaly_statuses(bad_raw, tmp_path / "bad")

    assert set(statuses) == set(_ANOMALY_TESTS), f"missing results: {statuses}"
    failed = {t for t, s in statuses.items() if s == "fail"}
    assert failed == set(_ANOMALY_TESTS), f"expected all anomaly tests to fail, got {statuses}"


def test_anomaly_checks_pass_on_good_data(tmp_path: Path) -> None:
    statuses = _anomaly_statuses(_RAW_V1, tmp_path / "good")
    assert all(s == "pass" for s in statuses.values()), f"unexpected failures: {statuses}"
