"""Generate the raw_v2 fixture: a synthetic "next ChEMBL release" delta on raw/.

Two documented changes vs. raw/ exercise the incremental fact and the SCD2 snapshot
(F2.2). Re-run to regenerate deterministically:

    python tests/fixtures/build_raw_v2.py

Delta v1 -> v2:
  1. STATUS CHANGE (SCD2): CHEMBLM2 advances max_phase 2 -> 4 (research -> approved),
     which also flips is_approved_drug. The snapshot must open a second validity window.
  2. NEW ACTIVITIES (incremental): activity_id 13 (CHEMBLM6) and 14 (CHEMBLM7) are
     appended on the surviving target/assay (CHEMBL2095182 / CHEMBLA2), so both pass the
     staging filters and land as exactly two new fact_activity rows.
All other files are copied verbatim so the two releases differ only in these records.
"""
from __future__ import annotations

import shutil
from pathlib import Path

import pandas as pd

SRC = Path(__file__).parent / "raw"
DST = Path(__file__).parent / "raw_v2"


def main() -> None:
    DST.mkdir(exist_ok=True)
    # Copy every raw parquet verbatim first; mutate the two that change below.
    for p in sorted(SRC.glob("raw_*.parquet")):
        shutil.copy2(p, DST / p.name)

    # 1. Status change: CHEMBLM2 becomes an approved drug (max_phase 2 -> 4).
    mol = pd.read_parquet(SRC / "raw_molecules.parquet")
    mol.loc[mol["molecule_chembl_id"] == "CHEMBLM2", "max_phase"] = 4.0
    mol.to_parquet(DST / "raw_molecules.parquet", index=False)

    # 2. Two new activities appended (stable, higher activity_ids).
    act = pd.read_parquet(SRC / "raw_activities.parquet")
    new = pd.DataFrame(
        [
            {
                "activity_id": 13,
                "molecule_chembl_id": "CHEMBLM6",
                "target_chembl_id": "CHEMBL2095182",
                "assay_chembl_id": "CHEMBLA2",
                "standard_type": "IC50",
                "standard_relation": "=",
                "standard_value": 90.0,
                "standard_units": "nM",
                "pchembl_value": 7.0,
                "data_validity_comment": None,
                "document_chembl_id": "CHEMBLD13",
            },
            {
                "activity_id": 14,
                "molecule_chembl_id": "CHEMBLM7",
                "target_chembl_id": "CHEMBL2095182",
                "assay_chembl_id": "CHEMBLA2",
                "standard_type": "EC50",
                "standard_relation": "=",
                "standard_value": 150.0,
                "standard_units": "nM",
                "pchembl_value": 6.8,
                "data_validity_comment": None,
                "document_chembl_id": "CHEMBLD14",
            },
        ]
    )
    act = pd.concat([act, new[act.columns]], ignore_index=True)
    act.to_parquet(DST / "raw_activities.parquet", index=False)

    print(f"raw_v2 written to {DST} ({len(list(DST.glob('*.parquet')))} files)")


if __name__ == "__main__":
    main()
