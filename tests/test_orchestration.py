"""Tests for the Dagster orchestration layer (orchestration/definitions.py).

These are collection-safe: if the dbt manifest has not been parsed yet they skip
rather than erroring at import (``@dbt_assets`` needs the manifest at import time).
CI parses the manifest before running pytest, so they execute there.
"""

from pathlib import Path

import pytest

_REPO = Path(__file__).resolve().parents[1]
_MANIFEST = _REPO / "dbt" / "target" / "manifest.json"
_FIXTURES = _REPO / "tests" / "fixtures" / "raw"


def _require_manifest():
    if not _MANIFEST.exists():
        pytest.skip(
            "dbt manifest missing; run "
            "`dbt parse --project-dir dbt --profiles-dir dbt/profiles` first"
        )


def test_definitions_load_with_full_medallion_graph():
    _require_manifest()
    from orchestration.definitions import defs

    graph = defs.resolve_asset_graph()
    keys = {k.to_user_string() for k in graph.get_all_asset_keys()}

    # Governed raw roots (published by the ChEMBL extract multi-asset).
    for table in ("activities", "molecules", "targets", "assays"):
        assert f"raw/{table}" in keys
    # One asset per dbt model, namespaced by layer.
    assert "staging/stg_activities" in keys
    assert "marts/fact_activity" in keys
    assert "analytics/mart_target_sar" in keys


def test_extract_roots_are_executable():
    _require_manifest()
    from dagster import AssetKey

    from orchestration.definitions import defs

    graph = defs.resolve_asset_graph()
    # The four ChEMBL raw tables are governed (materialisable), not external inputs.
    assert graph.get(AssetKey(["raw", "activities"])).is_executable


def test_lineage_connects_raw_to_staging_to_fact():
    _require_manifest()
    from dagster import AssetKey

    from orchestration.definitions import defs

    graph = defs.resolve_asset_graph()
    stg_parents = {
        d.to_user_string()
        for d in graph.get(AssetKey(["staging", "stg_activities"])).parent_keys
    }
    assert "raw/activities" in stg_parents

    fact_parents = {
        d.to_user_string()
        for d in graph.get(AssetKey(["marts", "fact_activity"])).parent_keys
    }
    assert "staging/stg_activities" in fact_parents


def test_jobs_and_schedule_defined():
    _require_manifest()
    from orchestration.definitions import defs

    assert defs.resolve_job_def("sarvault_pipeline") is not None
    assert defs.resolve_job_def("sarvault_transform") is not None
    schedules = list(defs.schedules) if defs.schedules else []
    assert any(s.name == "daily_refresh" for s in schedules)


def test_dbt_tests_surface_as_asset_checks():
    _require_manifest()
    from orchestration.definitions import defs

    graph = defs.resolve_asset_graph()
    # 69 dbt data tests + 2 explicit warehouse checks = 71.
    assert len(list(graph.asset_check_keys)) >= 70


def test_transform_job_builds_warehouse_and_checks_pass(tmp_path, monkeypatch):
    _require_manifest()
    # DbtCliResource runs dbt with cwd = the dbt project dir, so the raw dir must
    # be ABSOLUTE (a relative path would resolve against dbt/, not the repo root).
    monkeypatch.setenv("SARVAULT_RAW_DIR", str(_FIXTURES))
    monkeypatch.setenv("DUCKDB_PATH", str(tmp_path / "orch.duckdb"))

    from orchestration.definitions import defs

    result = defs.resolve_job_def("sarvault_transform").execute_in_process()
    assert result.success

    evals = {e.check_name: e.passed for e in result.get_asset_check_evaluations()}
    # Explicit warehouse smoke checks.
    assert evals.get("fact_activity_not_empty") is True
    assert evals.get("pchembl_within_range") is True
    # dbt data tests also ran as checks.
    assert len(evals) >= 70
