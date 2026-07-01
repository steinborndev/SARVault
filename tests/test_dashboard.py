"""M5: tests for the dashboard data-access, logic and chart builders."""

import importlib
import os
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go
import pytest

from dashboard import charts, data, logic


# --- chart builders (pure) ---
def _sar_df():
    return pd.DataFrame(
        {
            "compound_key": [1, 2],
            "target_key": [1, 1],
            "molecule_chembl_id": ["CHEMBL1", "CHEMBL2"],
            "target_chembl_id": ["CHEMBLT", "CHEMBLT"],
            "target_pref_name": ["Tubulin", "Tubulin"],
            "median_pchembl": [7.5, 6.2],
            "max_pchembl": [7.5, 6.2],
            "n_measurements": [1, 1],
            "n_assays": [1, 1],
        }
    )


def test_sar_ranking_bar_returns_figure():
    assert isinstance(charts.sar_ranking_bar(_sar_df()), go.Figure)


def test_selectivity_scatter_returns_figure():
    df = pd.DataFrame(
        {
            "molecule_chembl_id": ["CHEMBL1"],
            "best_pchembl": [8.0],
            "selectivity_index": [1.5],
            "best_target": ["CHEMBLT"],
            "n_targets": [2],
        }
    )
    assert isinstance(charts.selectivity_scatter(df), go.Figure)


def test_chemical_space_charts_handle_nulls():
    df = pd.DataFrame(
        {
            "molecule_chembl_id": ["A", "B"],
            "pref_name": ["a", "b"],
            "mw_freebase": [350.0, None],
            "alogp": [2.0, 3.0],
            "psa": [60.0, 70.0],
            "qed_weighted": [0.7, 0.6],
            "rotatable_bonds": [5, 4],
            "is_approved_drug": [True, False],
            "best_pchembl": [8.0, None],
        }
    )
    assert isinstance(charts.chemical_space_scatter(df, "mw_freebase", "alogp"), go.Figure)
    assert isinstance(charts.property_histogram(df, "mw_freebase"), go.Figure)


def test_confidence_bar_returns_figure():
    df = pd.DataFrame({"confidence_score": [5, 8, 9], "n_activities": [1221, 291, 710], "n_assays": [1, 1, 1]})
    assert isinstance(charts.confidence_bar(df), go.Figure)


# --- pure scope / metrics logic ---
def _scope_fixture():
    target_sar = pd.DataFrame(
        {
            "compound_key": [1, 1, 2, 3],
            "target_pref_name": ["Tubulin", "DNA topoisomerase 1", "Tubulin", "DNA topoisomerase 1"],
            "n_measurements": [2, 1, 3, 1],
        }
    )
    selectivity = pd.DataFrame({"compound_key": [1, 2, 3], "n_targets": [2, 1, 1]})
    chem = pd.DataFrame({"compound_key": [1, 2, 3], "is_approved_drug": [True, False, False]})
    return target_sar, selectivity, chem


def test_scope_compound_keys():
    target_sar, _, _ = _scope_fixture()
    assert logic.scope_compound_keys(target_sar, ["Tubulin"]) == {1, 2}
    assert logic.scope_compound_keys(target_sar, None) == {1, 2, 3}


def test_overview_metrics_respects_scope():
    target_sar, selectivity, chem = _scope_fixture()
    all_m = logic.overview_metrics(target_sar, selectivity, chem, None)
    assert all_m["compounds"] == 3
    assert all_m["pairs"] == 4
    assert all_m["multi_target"] == 1
    assert all_m["approved"] == 1
    tub = logic.overview_metrics(target_sar, selectivity, chem, ["Tubulin"])
    assert tub["compounds"] == 2  # compounds 1 and 2
    assert tub["activities"] == 5  # 2 + 3 measurements
    assert tub["targets"] == 1


# --- view modules import cleanly ---
def test_view_modules_import():
    for name in ("overview", "sar", "selectivity", "chemical_space", "data_quality"):
        importlib.import_module(f"dashboard.views.{name}")


# --- data access against a real (fixture-built) warehouse ---
_WAREHOUSE = Path(os.environ.get("DUCKDB_PATH", "warehouse.duckdb"))


@pytest.mark.skipif(not _WAREHOUSE.exists(), reason="warehouse not built")
def test_data_access_against_warehouse():
    con = data.connect(_WAREHOUSE)
    assert {"median_pchembl", "max_pchembl"}.issubset(data.load_target_sar(con).columns)
    assert data.list_target_names(con)
    assert {"target", "n_pairs"}.issubset(data.target_summary(con).columns)
    assert {"confidence_score", "n_activities"}.issubset(data.confidence_distribution(con).columns)
    assert {"table", "rows"}.issubset(data.layer_counts(con).columns)
    cfg = data.pipeline_config()
    assert cfg["chembl_version"] and cfg["min_confidence_score"] >= 0
