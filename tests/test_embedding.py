"""Tests for the UMAP embedding stage (extract/embedding.py).

The reproducibility test is the important one: with a fixed seed and single thread
the projection is a pure function of the fingerprints, so two runs must agree. The
tests use the committed cheminfo fixture (real molecule fingerprints), which forms a
well-defined similarity graph -- unlike random sparse bit vectors, which are nearly
disjoint and make any manifold method degenerate.
"""

from pathlib import Path

import duckdb
import numpy as np
import pytest

pytest.importorskip("umap")

from extract import embedding  # noqa: E402

_FIXTURE = (
    Path(__file__).resolve().parents[1]
    / "tests" / "fixtures" / "raw" / "raw_compound_cheminfo.parquet"
)


def _fixture_onbits():
    rows = duckdb.sql(
        f"select molecule_chembl_id, ecfp4_onbits from read_parquet('{_FIXTURE}') "
        "order by molecule_chembl_id"
    ).df()
    return rows["molecule_chembl_id"].tolist(), rows["ecfp4_onbits"].tolist()


def test_embedding_is_reproducible_for_a_fixed_seed():
    ids, bits = _fixture_onbits()
    a = embedding.compute_embedding(ids, bits, seed=42)
    b = embedding.compute_embedding(ids, bits, seed=42)
    assert list(a["molecule_chembl_id"]) == list(b["molecule_chembl_id"])
    coords_a = a[["umap_x", "umap_y"]].to_numpy()
    coords_b = b[["umap_x", "umap_y"]].to_numpy()
    assert np.isfinite(coords_a).all()
    assert np.allclose(coords_a, coords_b, atol=1e-6)


def test_embedding_shape_and_columns():
    ids, bits = _fixture_onbits()
    out = embedding.compute_embedding(ids, bits, seed=42)
    assert list(out.columns) == ["molecule_chembl_id", "umap_x", "umap_y"]
    assert len(out) == len(ids)
    assert out[["umap_x", "umap_y"]].notna().all().all()


def test_embedding_skips_when_too_few_compounds():
    ids, bits = _fixture_onbits()
    k = embedding.MIN_COMPOUNDS - 1
    out = embedding.compute_embedding(ids[:k], bits[:k], seed=42)
    assert out.empty


def test_build_embedding_stamps_provenance(tmp_path):
    import shutil

    shutil.copy(_FIXTURE, tmp_path / "raw_compound_cheminfo.parquet")
    df = embedding.build_embedding(raw_dir=tmp_path, seed=42)
    assert not df.empty
    for col in embedding.PROVENANCE_COLUMNS:
        assert col in df.columns
    assert (df["_source"] == "umap").all()
    assert (df["_embedding_seed"] == 42).all()
    assert df["_row_hash"].is_unique
