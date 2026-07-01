"""Tests for the dashboard data-access, logic, chem and chart builders."""

import importlib
import os
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go
import pytest

from dashboard import charts, chem, data, logic


# --- chart builders (pure) ---
def _sar_df():
    return pd.DataFrame(
        {
            "compound_key": [1, 2],
            "molecule_chembl_id": ["CHEMBL1", "CHEMBL2"],
            "target_pref_name": ["Tubulin", "Tubulin"],
            "median_pchembl": [7.5, 6.2],
            "max_pchembl": [7.5, 6.2],
            "n_measurements": [1, 1],
            "n_assays": [1, 1],
        }
    )


def test_sar_ranking_bar_returns_figure():
    assert isinstance(charts.sar_ranking_bar(_sar_df()), go.Figure)


def test_compound_potency_bar_returns_figure():
    df = pd.DataFrame(
        {
            "target": ["Tubulin"],
            "median_pchembl": [7.5],
            "max_pchembl": [7.5],
            "n_measurements": [1],
        }
    )
    assert isinstance(charts.compound_potency_bar(df), go.Figure)


def test_confidence_bar_returns_figure():
    df = pd.DataFrame({"confidence_score": [5, 9], "n_activities": [10, 5], "n_assays": [1, 1]})
    assert isinstance(charts.confidence_bar(df), go.Figure)


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


# --- RDKit structure rendering ---
def test_smiles_to_svg_valid():
    svg = chem.smiles_to_svg("CCO")
    assert svg and "<svg" in svg


def test_smiles_to_svg_invalid_returns_none():
    assert chem.smiles_to_svg("not-a-smiles") is None
    assert chem.smiles_to_svg("") is None


# --- pure scope logic ---
def _scope_fixtures():
    target_sar = pd.DataFrame(
        {
            "compound_key": [1, 1, 2, 3],
            "target_pref_name": [
                "Tubulin",
                "DNA topoisomerase 1",
                "Tubulin",
                "DNA topoisomerase 1",
            ],
            "n_measurements": [2, 1, 3, 1],
        }
    )
    catalog = pd.DataFrame(
        {
            "compound_key": [1, 2, 3],
            "is_approved_drug": [True, False, False],
            "best_pchembl": [9.0, 6.0, 5.0],
            "n_targets": [2, 1, 1],
        }
    )
    return target_sar, catalog


def test_resolve_scope_keys_target_and_approval():
    target_sar, catalog = _scope_fixtures()
    assert logic.resolve_scope_keys(target_sar, catalog, {}) == {1, 2, 3}
    assert logic.resolve_scope_keys(target_sar, catalog, {"targets": ["Tubulin"]}) == {1, 2}
    assert logic.resolve_scope_keys(target_sar, catalog, {"approval": "approved"}) == {1}
    assert logic.resolve_scope_keys(target_sar, catalog, {"approval": "research"}) == {2, 3}


def test_resolve_scope_keys_min_potency():
    target_sar, catalog = _scope_fixtures()
    assert logic.resolve_scope_keys(target_sar, catalog, {"min_pchembl": 6.0}) == {1, 2}


def test_ro5_breakdown():
    row = pd.Series(
        {"mw_freebase": 451.0, "alogp": 6.9, "hbd": 1, "hba": 4, "num_ro5_violations": 1}
    )
    result = logic.ro5_breakdown(row)
    assert result["violations"] == 1
    by_label = {item["label"]: item["pass"] for item in result["items"]}
    assert by_label["logP ≤ 5"] is False
    assert by_label["MW ≤ 500"] is True

    missing = logic.ro5_breakdown(
        pd.Series({"mw_freebase": None, "alogp": 2.0, "hbd": 1, "hba": 4})
    )
    assert missing["violations"] == 0
    assert missing["items"][0]["pass"] is None


def test_overview_metrics_respects_scope():
    target_sar, catalog = _scope_fixtures()
    all_m = logic.overview_metrics(target_sar, catalog, {})
    assert all_m["compounds"] == 3
    assert all_m["pairs"] == 4
    assert all_m["multi_target"] == 1
    assert all_m["approved"] == 1
    tub = logic.overview_metrics(target_sar, catalog, {"targets": ["Tubulin"]})
    assert tub["compounds"] == 2
    assert tub["activities"] == 5
    assert tub["targets"] == 1


# --- view modules import cleanly ---
def test_view_modules_import():
    for name in (
        "overview",
        "compound_library",
        "sar",
        "selectivity",
        "chemical_space",
        "data_quality",
    ):
        importlib.import_module(f"dashboard.views.{name}")


def test_compound_property_formatting_is_string():
    from dashboard.views.compound_library import _fmt

    assert _fmt(350.44) == "350.44"
    assert _fmt("N") == "N"
    assert _fmt(3) == "3"
    assert _fmt(None) == "—"
    assert _fmt(float("nan")) == "—"


# --- data access against a real (fixture-built) warehouse ---
_WAREHOUSE = Path(os.environ.get("DUCKDB_PATH", "warehouse.duckdb"))


@pytest.mark.skipif(not _WAREHOUSE.exists(), reason="warehouse not built")
def test_data_access_against_warehouse():
    con = data.connect(_WAREHOUSE)
    assert {"median_pchembl", "max_pchembl"}.issubset(data.load_target_sar(con).columns)
    cat = data.load_compound_catalog(con)
    assert {"canonical_smiles", "best_pchembl", "n_targets"}.issubset(cat.columns)
    key = int(cat["compound_key"].iloc[0])
    assert {"target", "median_pchembl"}.issubset(data.compound_target_profile(con, key).columns)
    assert {"display_name", "xref_id", "url"}.issubset(data.compound_xrefs(con, key).columns)
    assert {"ligand_code", "pdb_id", "title", "method", "year", "resolution"}.issubset(
        data.compound_pdb_entries(con, key).columns
    )
    assert data.list_target_names(con)
    cfg = data.pipeline_config()
    assert cfg["chembl_version"] and cfg["min_confidence_score"] >= 0
