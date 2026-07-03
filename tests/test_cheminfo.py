"""Tests for the RDKit cheminfo compute stage (extract/cheminfo.py).

The determinism anchor pins the ECFP4 bit pattern for a fixed SMILES: given the
pinned RDKit version, the fingerprint is a pure function of the input, so a drift
in the hash flags a fingerprint-definition change (radius, bit size, algorithm).
"""

import hashlib

import pytest

rdkit = pytest.importorskip("rdkit")

from extract import cheminfo  # noqa: E402

# Aspirin. ECFP4 (Morgan r2, 2048-bit) -> hex -> sha256, captured on RDKit 2026.03.3.
_ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"
_ASPIRIN_ECFP4_SHA256 = "3896710a7bb107c9cb95bea4b1214a0fba201c12f91a566b80e98ec83a3620ec"


def test_fingerprint_is_deterministic():
    info = cheminfo.compute_cheminfo(_ASPIRIN)
    assert info is not None
    # 2048-bit vector -> 512 hex chars.
    assert len(info["ecfp4_hex"]) == cheminfo.FP_NBITS // 4
    digest = hashlib.sha256(info["ecfp4_hex"].encode()).hexdigest()
    assert digest == _ASPIRIN_ECFP4_SHA256


def test_onbits_are_sorted_and_match_popcount():
    info = cheminfo.compute_cheminfo(_ASPIRIN)
    onbits = info["ecfp4_onbits"]
    assert onbits == sorted(onbits)
    assert len(onbits) == info["n_onbits"]
    assert all(0 <= b < cheminfo.FP_NBITS for b in onbits)


def test_fingerprint_stable_across_calls():
    a = cheminfo.compute_cheminfo(_ASPIRIN)
    b = cheminfo.compute_cheminfo(_ASPIRIN)
    assert a["ecfp4_hex"] == b["ecfp4_hex"]
    assert a["n_onbits"] == b["n_onbits"]


def test_murcko_scaffold_extracted_for_cyclic_molecule():
    info = cheminfo.compute_cheminfo(_ASPIRIN)
    # Aspirin's Bemis-Murcko scaffold is a bare benzene ring.
    assert info["murcko_scaffold_smiles"] == "c1ccccc1"
    assert info["heavy_atom_count"] == 13


def test_acyclic_molecule_has_no_scaffold():
    info = cheminfo.compute_cheminfo("CCO")  # ethanol: no ring system
    assert info is not None
    assert info["murcko_scaffold_smiles"] is None


def test_unparseable_or_empty_smiles_returns_none():
    assert cheminfo.compute_cheminfo("not-a-smiles") is None
    assert cheminfo.compute_cheminfo("") is None
    assert cheminfo.compute_cheminfo(None) is None


def test_build_cheminfo_drops_unparseable_and_stamps_provenance(tmp_path):
    import pandas as pd

    # Two parseable, one unparseable SMILES landed as a minimal raw_molecules file.
    pd.DataFrame(
        {
            "molecule_chembl_id": ["C1", "C2", "C3"],
            "molecule_structures.canonical_smiles": ["c1ccccc1", "CCO", "not-a-smiles"],
        }
    ).to_parquet(tmp_path / "raw_molecules.parquet", index=False)

    df = cheminfo.build_cheminfo(raw_dir=tmp_path)

    assert len(df) == 2  # the unparseable row is dropped
    assert set(df["molecule_chembl_id"]) == {"C1", "C2"}
    for col in cheminfo.PROVENANCE_COLUMNS:
        assert col in df.columns
    assert (df["_source"] == "rdkit").all()
    assert (df["_rdkit_version"] == cheminfo.rdkit_version()).all()
    assert df["_row_hash"].is_unique
