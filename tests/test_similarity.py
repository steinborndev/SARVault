"""Tests for structural similarity (Tanimoto / nearest neighbours) and substructure.

The Tanimoto and nearest-neighbour helpers are pure (no RDKit) and use tiny
hand-checkable hex fingerprints. The substructure helpers need RDKit and are
skipped when it is absent.
"""

import pandas as pd
import pytest

from dashboard import logic


# --- Tanimoto over hex fingerprints (pure) ---
def test_tanimoto_identity_and_disjoint():
    assert logic.tanimoto_hex("ff", "ff") == 1.0  # identical
    assert logic.tanimoto_hex("f0", "0f") == 0.0  # disjoint bit halves


def test_tanimoto_partial_overlap():
    # a = 1111_1111 (8 on), b = 0000_1111 (4 on); intersection 4, union 8 -> 0.5
    assert logic.tanimoto_hex("ff", "0f") == 0.5


def test_tanimoto_empty_inputs_are_zero():
    assert logic.tanimoto_hex("", "ff") == 0.0
    assert logic.tanimoto_hex("ff", "") == 0.0
    assert logic.tanimoto_hex("00", "00") == 0.0  # both empty fingerprints


# --- nearest_neighbors (pure) ---
def _fp():
    return pd.DataFrame(
        {
            "compound_key": [1, 2, 3],
            "molecule_chembl_id": ["A", "B", "C"],
            "ecfp4_hex": ["ff", "f0", "ff"],  # C is identical to the query A
        }
    )


def _catalog():
    return pd.DataFrame(
        {
            "compound_key": [1, 2, 3],
            "molecule_chembl_id": ["A", "B", "C"],
            "pref_name": ["A_name", "B_name", "C_name"],
            "best_pchembl": [7.0, 6.0, 5.0],
            "best_target": ["T1", "T1", "T2"],
        }
    )


def test_nearest_neighbors_ranks_by_similarity_and_excludes_self():
    nn = logic.nearest_neighbors(1, _fp(), _catalog())
    # Self (A) excluded by default; C (identical) ranks above B (partial).
    assert list(nn["molecule_chembl_id"]) == ["C", "B"]
    assert nn.iloc[0]["tanimoto"] == 1.0
    assert nn.iloc[1]["tanimoto"] == 0.5


def test_nearest_neighbors_delta_pchembl_relative_to_query():
    nn = logic.nearest_neighbors(1, _fp(), _catalog())
    by_id = nn.set_index("molecule_chembl_id")
    assert by_id.loc["C", "delta_pchembl"] == -2.0  # 5.0 - 7.0
    assert by_id.loc["B", "delta_pchembl"] == -1.0  # 6.0 - 7.0


def test_nearest_neighbors_self_query_returns_identity_at_one():
    nn = logic.nearest_neighbors(1, _fp(), _catalog(), include_self=True)
    top = nn.iloc[0]
    assert top["molecule_chembl_id"] == "A"
    assert top["tanimoto"] == 1.0
    assert top["delta_pchembl"] == 0.0


def test_nearest_neighbors_min_similarity_filter():
    nn = logic.nearest_neighbors(1, _fp(), _catalog(), min_similarity=0.6)
    # Only the identical C (1.0) clears 0.6; B (0.5) is dropped.
    assert list(nn["molecule_chembl_id"]) == ["C"]


def test_nearest_neighbors_missing_query_is_empty_not_error():
    nn = logic.nearest_neighbors(999, _fp(), _catalog())
    assert nn.empty
    assert list(nn.columns) == [
        "molecule_chembl_id",
        "pref_name",
        "tanimoto",
        "best_pchembl",
        "delta_pchembl",
        "best_target",
    ]


# --- substructure (RDKit) ---
rdkit = pytest.importorskip("rdkit")

from dashboard import chem  # noqa: E402


def test_is_valid_smarts():
    assert chem.is_valid_smarts("c1ccccc1") is True
    assert chem.is_valid_smarts("C(=O)N") is True
    assert chem.is_valid_smarts("[[not-smarts") is False
    assert chem.is_valid_smarts("") is False


def test_has_substructure_matches_expected_compounds():
    # Benzene ring: present in toluene and aniline, absent from ethanol.
    assert chem.has_substructure("Cc1ccccc1", "c1ccccc1") is True
    assert chem.has_substructure("Nc1ccccc1", "c1ccccc1") is True
    assert chem.has_substructure("CCO", "c1ccccc1") is False


def test_has_substructure_none_on_unusable_input():
    assert chem.has_substructure("CCO", "[[bad") is None  # bad query
    assert chem.has_substructure("not-a-smiles", "c1ccccc1") is None  # bad molecule
    assert chem.has_substructure("", "c1ccccc1") is None


def test_smiles_to_svg_highlight_renders_and_is_backward_compatible():
    plain = chem.smiles_to_svg("Cc1ccccc1")
    assert plain is not None and "<svg" in plain
    highlighted = chem.smiles_to_svg("Cc1ccccc1", highlight_smarts="c1ccccc1")
    assert highlighted is not None and "<svg" in highlighted
    # An invalid highlight SMARTS must not break rendering.
    assert chem.smiles_to_svg("Cc1ccccc1", highlight_smarts="[[bad") is not None
