"""Generate a deliberately anomalous raw fixture to prove the observability checks fire.

`write_bad_fixture(dst)` copies the good v1 raw fixture and injects three faults, one per
anomaly guard in dbt/tests/assert_*.sql, so a build against it trips every check:

  * ROW-COUNT FLOOR  — all but two activities are dropped, so fact_activity falls below
    the per-layer floor (a stand-in for a silently-empty extract).
  * NULL-RATE CEILING — canonical_smiles is nulled for all but one molecule, blowing past
    the tolerated missing-structure rate.
  * DISTRIBUTION SHIFT — the surviving pChEMBL values are set to 99 (a unit/merge-error
    stand-in), pushing the mean far outside its plausible band.

The output is not committed; tests/test_observability.py writes it to a temp dir.
"""
from __future__ import annotations

import shutil
from pathlib import Path

import pandas as pd

SRC = Path(__file__).parent / "raw"
_SMILES_COL = "molecule_structures.canonical_smiles"


def write_bad_fixture(dst: Path) -> Path:
    """Materialise the anomalous raw fixture into ``dst`` and return the path."""
    dst.mkdir(parents=True, exist_ok=True)
    for p in sorted(SRC.glob("raw_*.parquet")):
        shutil.copy2(p, dst / p.name)

    # ROW-COUNT FLOOR + DISTRIBUTION SHIFT: keep only two surviving activities, and set
    # their pChEMBL to an implausible value.
    act = pd.read_parquet(SRC / "raw_activities.parquet")
    keep = act[act["activity_id"].isin([9, 12])].copy()  # both pass staging filters
    keep["pchembl_value"] = 99.0
    keep.to_parquet(dst / "raw_activities.parquet", index=False)

    # NULL-RATE CEILING: null out canonical_smiles for all but the first molecule.
    mol = pd.read_parquet(SRC / "raw_molecules.parquet")
    mol.loc[mol.index[1:], _SMILES_COL] = None
    mol.to_parquet(dst / "raw_molecules.parquet", index=False)

    return dst


if __name__ == "__main__":
    out = write_bad_fixture(Path(__file__).parent / "raw_bad")
    print(f"bad fixture written to {out}")
