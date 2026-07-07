"""Build the synthetic raw_cell_activities fixture for the F3.2 cellular lineage.

    python tests/fixtures/build_cell_activities.py

Writes tests/fixtures/raw/raw_cell_activities.parquet: a few reference-payload
cellular readouts across two classes and two cell lines, plus rows that must be
filtered out in staging (non-'=' relation, non-concentration unit).
"""

from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

OUT = Path(__file__).resolve().parent / "raw" / "raw_cell_activities.parquet"

_ROWS = [
    # molecule, name, class, cell target, cell name, assay, type, rel, value, units
    ("CHEMBL65", "Camptothecin", "topo1_inhibitor", "CHEMBLC1", "HeLa", "CHEMBLA1", "GI50", "=", 44.0, "nM"),
    ("CHEMBL65", "Camptothecin", "topo1_inhibitor", "CHEMBLC1", "HeLa", "CHEMBLA2", "GI50", "=", 60.0, "nM"),
    ("CHEMBL65", "Camptothecin", "topo1_inhibitor", "CHEMBLC2", "A549", "CHEMBLA3", "GI50", "=", 30.0, "nM"),
    ("CHEMBL84", "Topotecan", "topo1_inhibitor", "CHEMBLC1", "HeLa", "CHEMBLA4", "IC50", "=", 100.0, "nM"),
    ("CHEMBL4082989", "MMAE", "tubulin_inhibitor", "CHEMBLC3", "DU-145", "CHEMBLA5", "GI50", "=", 0.424, "nM"),
    # filtered out in staging: non-'=' relation and non-concentration unit
    ("CHEMBL65", "Camptothecin", "topo1_inhibitor", "CHEMBLC2", "A549", "CHEMBLA6", "GI50", ">", 1000.0, "nM"),
    ("CHEMBL84", "Topotecan", "topo1_inhibitor", "CHEMBLC2", "A549", "CHEMBLA7", "GI50", "=", 5.0, "ug.mL-1"),
]

_COLS = [
    "molecule_chembl_id", "reference_name", "payload_class", "target_chembl_id",
    "target_pref_name", "assay_chembl_id", "standard_type", "standard_relation",
    "standard_value", "standard_units",
]


def build() -> Path:
    df = pd.DataFrame(_ROWS, columns=_COLS)
    df.insert(0, "activity_id", range(1, len(df) + 1))
    df["molecule_pref_name"] = df["reference_name"].str.upper()
    df["canonical_smiles"] = "CCO"  # placeholder; the cytotox mart does not parse it
    df["_fetch_ts"] = datetime.now(timezone.utc).isoformat()
    df["_source_endpoint"] = "activity"
    df["_chembl_version"] = "36"
    df["_row_hash"] = ["cellhash%03d" % i for i in df["activity_id"]]
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(OUT, index=False)
    return OUT


if __name__ == "__main__":
    out = build()
    print(f"wrote {out} ({len(pd.read_parquet(out))} rows)")
