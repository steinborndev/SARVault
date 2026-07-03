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


def test_sar_heatmap_returns_figure():
    assert isinstance(charts.sar_heatmap(_sar_df()), go.Figure)
    # empty scope must not raise
    assert isinstance(charts.sar_heatmap(_sar_df().iloc[0:0]), go.Figure)


def test_target_potency_violin_returns_figure():
    assert isinstance(charts.target_potency_violin(_sar_df()), go.Figure)


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


def test_standard_type_bar_returns_figure():
    df = pd.DataFrame(
        {"standard_type": ["IC50", "Ki"], "n_activities": [10, 3], "n_compounds": [8, 2]}
    )
    assert isinstance(charts.standard_type_bar(df), go.Figure)


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


def test_heavy_atom_count():
    assert chem.heavy_atom_count("CCO") == 3  # ethanol: C, C, O
    assert chem.heavy_atom_count("not-a-smiles") is None
    assert chem.heavy_atom_count("") is None


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
            "has_pdb": [True, False, True],
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


def test_resolve_scope_keys_structure_only():
    target_sar, catalog = _scope_fixtures()
    assert logic.resolve_scope_keys(target_sar, catalog, {"structure_only": True}) == {1, 3}
    assert logic.resolve_scope_keys(target_sar, catalog, {"structure_only": False}) == {1, 2, 3}
    # structure_only composes with the other facets (approved -> only compound 1)
    assert (
        logic.resolve_scope_keys(
            target_sar, catalog, {"structure_only": True, "approval": "approved"}
        )
        == {1}
    )


def test_ligand_and_lipophilic_efficiency():
    assert logic.ligand_efficiency(7.3, 20) == pytest.approx(1.37 * 7.3 / 20)
    assert logic.ligand_efficiency(None, 20) is None
    assert logic.ligand_efficiency(7.3, 0) is None
    assert logic.ligand_efficiency(7.3, None) is None
    assert logic.lipophilic_efficiency(7.0, 2.0) == pytest.approx(5.0)
    assert logic.lipophilic_efficiency(7.0, None) is None


def test_add_efficiency():
    df = pd.DataFrame(
        {
            "compound_key": [1, 2],
            "canonical_smiles": ["CCO", "bad"],
            "best_pchembl": [6.0, 8.0],
            "alogp": [1.5, None],
        }
    )
    out = logic.add_efficiency(df, {"CCO": 20, "bad": None})
    assert out.loc[0, "heavy_atoms"] == 20
    assert out.loc[0, "ligand_efficiency"] == pytest.approx(1.37 * 6.0 / 20)
    assert out.loc[0, "lipophilic_efficiency"] == pytest.approx(4.5)
    # unparseable SMILES / missing logP -> NaN in the float column, no exception
    assert pd.isna(out.loc[1, "ligand_efficiency"])
    assert pd.isna(out.loc[1, "lipophilic_efficiency"])


def test_efficiency_scatter_returns_figure():
    df = pd.DataFrame(
        {
            "molecule_chembl_id": ["A", "B"],
            "best_pchembl": [7.0, None],
            "ligand_efficiency": [0.45, None],
            "heavy_atoms": [21, 30],
            "mw_freebase": [300.0, 420.0],
            "is_approved_drug": [True, False],
        }
    )
    assert isinstance(
        charts.efficiency_scatter(df, "ligand_efficiency", "ligand efficiency"), go.Figure
    )


def test_resolve_scope_keys_structure_only_missing_column():
    # A warehouse built before has_pdb existed must not crash the toggle; the
    # structure filter is simply skipped (graceful degradation).
    target_sar, catalog = _scope_fixtures()
    stale = catalog.drop(columns=["has_pdb"])
    assert logic.resolve_scope_keys(target_sar, stale, {"structure_only": True}) == {1, 2, 3}


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


def test_preselect_first_row():
    # First open: an absent key is seeded so row 0 is marked (matches the detail below).
    state = {}
    logic.preselect_first_row(state, "lib_rows")
    assert state["lib_rows"] == {"selection": {"rows": [0], "columns": [], "cells": []}}

    # A user's existing selection is never overridden on later reruns.
    state["lib_rows"] = {"selection": {"rows": [3], "columns": [], "cells": []}}
    logic.preselect_first_row(state, "lib_rows")
    assert state["lib_rows"]["selection"]["rows"] == [3]

    # Independent keys (e.g. per-scaffold member tables) are seeded independently.
    logic.preselect_first_row(state, "members_42")
    assert state["members_42"]["selection"]["rows"] == [0]
    assert "lib_rows" in state and state["lib_rows"]["selection"]["rows"] == [3]


def test_step_selection_walks_and_clamps():
    state = {"m": {"selection": {"rows": [10]}}}
    logic.step_selection(state, "m", 1, 68)
    assert state["m"]["selection"]["rows"] == [11]
    logic.step_selection(state, "m", -1, 68)
    assert state["m"]["selection"]["rows"] == [10]

    # Clamp at both ends rather than wrapping.
    state["m"] = {"selection": {"rows": [0]}}
    logic.step_selection(state, "m", -1, 68)
    assert state["m"]["selection"]["rows"] == [0]
    state["m"] = {"selection": {"rows": [67]}}
    logic.step_selection(state, "m", 1, 68)
    assert state["m"]["selection"]["rows"] == [67]

    # Absent key defaults to row 0; an empty table is a safe no-op.
    empty = {}
    logic.step_selection(empty, "m", 1, 68)
    assert empty["m"]["selection"]["rows"] == [1]
    logic.step_selection({}, "m", 1, 0)  # no rows -> no crash


def test_smiles_to_svg_scaffold_aligned_and_highlighted():
    scaffold = "c1ccc(-c2ccccc2)cc1"  # biphenyl core
    member = "Cc1ccc(-c2ccc(Cl)cc2)cc1"  # substituted members share the core
    aligned = chem.smiles_to_svg(
        member, scaffold_smiles=scaffold, align_to_scaffold=True, highlight_scaffold=True
    )
    assert aligned and "<svg" in aligned

    # A scaffold the member does not contain must not break rendering (plain fallback).
    fallback = chem.smiles_to_svg(
        member, scaffold_smiles="C1CCCCC1CCCC", align_to_scaffold=True, highlight_scaffold=True
    )
    assert fallback and "<svg" in fallback


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
        "activity_cliffs",
        "chemical_series",
        "chemical_space",
        "data_quality",
    ):
        importlib.import_module(f"dashboard.views.{name}")


def test_compound_property_formatting_is_string():
    from dashboard.compound_detail import _fmt

    assert _fmt(350.44) == "350.44"
    assert _fmt("N") == "N"
    assert _fmt(3) == "3"
    assert _fmt(None) == "-"
    assert _fmt(float("nan")) == "-"


# --- warehouse bootstrap (cloud fetch-at-boot) ---
def test_ensure_warehouse_existing_is_untouched(tmp_path):
    db = tmp_path / "warehouse.duckdb"
    db.write_bytes(b"already here")
    # a present file is returned as-is even if a url is given (no overwrite)
    out = data.ensure_warehouse(db_path=str(db), url="https://example.invalid/w.duckdb")
    assert out == db
    assert db.read_bytes() == b"already here"


def test_ensure_warehouse_missing_without_url_no_error(tmp_path):
    db = tmp_path / "warehouse.duckdb"
    # no url -> no download attempt, caller handles the missing file downstream
    out = data.ensure_warehouse(db_path=str(db), url=None)
    assert out == db
    assert not db.exists()


def test_parse_asset_url():
    owner, repo, tag, name = data._parse_asset_url(
        "https://github.com/Knusperftw/SARVault/releases/download/warehouse-v1/warehouse.duckdb"
    )
    assert (owner, repo, tag, name) == ("Knusperftw", "SARVault", "warehouse-v1", "warehouse.duckdb")
    with pytest.raises(ValueError):
        data._parse_asset_url("https://example.com/not/a/release")


def test_layer_counts_skips_unreadable_relations():
    # A shipped snapshot may lack some relations (e.g. staging views over absent
    # raw files); layer_counts must return what it can without raising.
    import duckdb

    con = duckdb.connect(":memory:")
    con.execute("create schema main_marts")
    con.execute("create table main_marts.dim_compound as select 1 as compound_key")
    out = data.layer_counts(con)  # other relations don't exist -> skipped, no error
    assert "dim_compound" in set(out["table"])
    assert int(out.loc[out["table"] == "dim_compound", "rows"].iloc[0]) == 1
    assert "stg_activities" not in set(out["table"])


# --- data access against a real (fixture-built) warehouse ---
_WAREHOUSE = Path(os.environ.get("DUCKDB_PATH", "warehouse.duckdb"))


@pytest.mark.skipif(not _WAREHOUSE.exists(), reason="warehouse not built")
def test_data_access_against_warehouse():
    con = data.connect(_WAREHOUSE)
    assert {"median_pchembl", "max_pchembl"}.issubset(data.load_target_sar(con).columns)
    cat = data.load_compound_catalog(con)
    assert {"canonical_smiles", "best_pchembl", "n_targets", "has_pdb", "n_pdb_entries"}.issubset(
        cat.columns
    )
    assert {"standard_type", "n_activities"}.issubset(
        data.standard_type_distribution(con).columns
    )
    key = int(cat["compound_key"].iloc[0])
    assert {"target", "median_pchembl"}.issubset(data.compound_target_profile(con, key).columns)
    assert {"display_name", "xref_id", "url"}.issubset(data.compound_xrefs(con, key).columns)
    assert {"ligand_code", "pdb_id", "title", "method", "year", "resolution"}.issubset(
        data.compound_pdb_entries(con, key).columns
    )
    assert data.list_target_names(con)
    cfg = data.pipeline_config()
    assert cfg["chembl_version"] and cfg["min_confidence_score"] >= 0
