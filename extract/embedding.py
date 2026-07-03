"""Compute a 2-D UMAP embedding of the ECFP4 fingerprints and land it to raw.

A property scatter (MW vs logP …) shows physicochemistry; a fingerprint embedding
shows *structural* neighbourhoods — which compounds are close in chemical space and
where potency concentrates. UMAP over the ECFP4 on-bit vectors with the Jaccard
(Tanimoto) metric is the standard cheminformatics projection.

The projection is fitted once over the whole compound set in the pipeline (not per
request in the dashboard) and is deterministic given a fixed seed and single thread,
so a rebuild from the same fingerprints reproduces the same coordinates. Reads
``raw_compound_cheminfo.parquet`` (the on-bit lists) and writes
``raw_compound_embedding.parquet``.

Usage:
    python -m extract.embedding
"""

import hashlib
from datetime import datetime, timezone
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd

RAW_DIR = Path(__file__).resolve().parents[1] / "raw"

FP_NBITS = 2048
EMBED_SEED = 42
# UMAP needs at least a few points to fit a manifold; below this we skip (the
# dashboard then simply has no embedding to plot for such a tiny set).
MIN_COMPOUNDS = 4

PROVENANCE_COLUMNS = ["_fetch_ts", "_source", "_umap_version", "_embedding_seed", "_row_hash"]


def umap_version() -> str:
    import umap

    return umap.__version__


def _onbits_matrix(onbits_lists) -> np.ndarray:
    """Dense uint8 fingerprint matrix (n_compounds x FP_NBITS) from on-bit lists."""
    m = np.zeros((len(onbits_lists), FP_NBITS), dtype=np.uint8)
    for i, bits in enumerate(onbits_lists):
        if bits is not None:
            m[i, list(bits)] = 1
    return m


def compute_embedding(chembl_ids, onbits_lists, seed: int = EMBED_SEED) -> pd.DataFrame:
    """2-D UMAP (Jaccard metric) over the fingerprints; deterministic for a given seed.

    Returns a DataFrame of (molecule_chembl_id, umap_x, umap_y). If there are too few
    compounds to fit a manifold, returns an empty frame.
    """
    import umap

    n = len(chembl_ids)
    if n < MIN_COMPOUNDS:
        return pd.DataFrame(columns=["molecule_chembl_id", "umap_x", "umap_y"])

    matrix = _onbits_matrix(onbits_lists)
    reducer = umap.UMAP(
        n_components=2,
        n_neighbors=min(15, n - 1),
        metric="jaccard",
        random_state=seed,
        n_jobs=1,  # required for reproducibility
    )
    coords = reducer.fit_transform(matrix)
    return pd.DataFrame(
        {
            "molecule_chembl_id": list(chembl_ids),
            "umap_x": coords[:, 0].astype(float),
            "umap_y": coords[:, 1].astype(float),
        }
    )


def _read_cheminfo(raw_dir: Path) -> pd.DataFrame:
    path = raw_dir / "raw_compound_cheminfo.parquet"
    return duckdb.sql(
        f"select molecule_chembl_id, ecfp4_onbits from read_parquet('{path}')"
    ).df()


def _row_hash(chembl_id: str, x: float, y: float) -> str:
    return hashlib.sha256(f"{chembl_id}\t{x}\t{y}".encode()).hexdigest()


def build_embedding(raw_dir: Path | str = RAW_DIR, seed: int = EMBED_SEED) -> pd.DataFrame:
    raw_dir = Path(raw_dir)
    cheminfo = _read_cheminfo(raw_dir)
    df = compute_embedding(
        cheminfo["molecule_chembl_id"].tolist(),
        cheminfo["ecfp4_onbits"].tolist(),
        seed=seed,
    )
    print(f"embedding: {len(df)} compounds projected (seed {seed})")
    if df.empty:
        return df
    df["_fetch_ts"] = datetime.now(timezone.utc).isoformat()
    df["_source"] = "umap"
    df["_umap_version"] = umap_version()
    df["_embedding_seed"] = seed
    df["_row_hash"] = [_row_hash(c, x, y) for c, x, y in zip(df.molecule_chembl_id, df.umap_x, df.umap_y)]
    return df


def land_embedding(
    df: pd.DataFrame, raw_dir: Path | str = RAW_DIR, validate_schema: bool = True
) -> Path:
    if validate_schema and not df.empty:
        from validation.schemas import validate

        df = validate("compound_embedding", df)
    raw_dir = Path(raw_dir)
    raw_dir.mkdir(parents=True, exist_ok=True)
    out_path = raw_dir / "raw_compound_embedding.parquet"
    df.to_parquet(out_path, index=False)
    return out_path


def main(raw_dir: Path | str = RAW_DIR, seed: int = EMBED_SEED) -> Path:
    df = build_embedding(raw_dir, seed=seed)
    out = land_embedding(df, raw_dir)
    print(f"embedding: {len(df)} rows -> {out}")
    return out


if __name__ == "__main__":
    main()
