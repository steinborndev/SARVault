"""Resolve PDBe ligand (het) codes to co-crystal PDB entries (extract stage 3).

Compounds carry a PDBe cross-reference whose id is a chemical-component (het)
code. The PDBe REST API's ``in_pdb`` endpoint lists every PDB entry that contains
a given component, so this stage turns each ligand code into the concrete
co-crystal structures the dashboard renders in its 3D viewer.

Run after the ChEMBL/UniChem extract (it reads their raw Parquet) and before
``dbt build``:

    python -m extract.pdbe
"""

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import duckdb
import pandas as pd

from extract.chembl_client import build_session

PDBE_IN_PDB = "https://www.ebi.ac.uk/pdbe/api/pdb/compound/in_pdb"
RAW_DIR = Path(__file__).resolve().parents[1] / "raw"
PROVENANCE_COLUMNS = ["_fetch_ts", "_source_endpoint", "_chembl_version", "_row_hash"]
# Common ligands appear in hundreds of entries; cap per code to keep the raw
# layer bounded and the resolution deterministic.
MAX_PER_LIGAND = 50


def ligand_codes_from_raw(raw_dir: Path | str = RAW_DIR) -> pd.DataFrame:
    """Collect distinct (molecule_chembl_id, ligand_code) PDBe references from raw/.

    Reads PDBe cross-references from both ChEMBL's molecule ``cross_references``
    and the UniChem bulk dump, mirroring how ``mart_compound_xref`` combines the
    two. Returns columns: molecule_chembl_id, ligand_code.
    """
    raw_dir = Path(raw_dir)
    con = duckdb.connect()
    frames = []
    mol_path = raw_dir / "raw_molecules.parquet"
    if mol_path.exists():
        frames.append(
            con.execute(
                """
                with exploded as (
                    select
                        molecule_chembl_id,
                        unnest(json_extract(cross_references, '$[*]')) as ref
                    from read_parquet(?)
                    where cross_references is not null and cross_references != 'null'
                )
                select distinct
                    molecule_chembl_id,
                    json_extract_string(ref, '$.xref_id') as ligand_code
                from exploded
                where lower(json_extract_string(ref, '$.xref_src')) = 'pdbe'
                  and json_extract_string(ref, '$.xref_id') is not null
                """,
                [str(mol_path)],
            ).df()
        )
    uni_path = raw_dir / "raw_xref_unichem.parquet"
    if uni_path.exists():
        frames.append(
            con.execute(
                """
                select distinct molecule_chembl_id, xref_id as ligand_code
                from read_parquet(?)
                where source = 'pdbe' and xref_id is not null
                """,
                [str(uni_path)],
            ).df()
        )
    if not frames:
        return pd.DataFrame(columns=["molecule_chembl_id", "ligand_code"])
    return pd.concat(frames, ignore_index=True).drop_duplicates().reset_index(drop=True)


def resolve_pdb_entries(
    ligand_codes,
    session=None,
    max_per_ligand: int = MAX_PER_LIGAND,
    timeout: int = 60,
) -> pd.DataFrame:
    """Resolve each ligand (het) code to its PDB entries via the PDBe in_pdb API.

    The endpoint returns ``{"<CODE>": ["1abc", ...]}``. Entries are lowercased,
    de-duplicated, sorted and capped at ``max_per_ligand`` per code for bounded,
    deterministic output. Missing components (404) are skipped. Returns columns:
    ligand_code, pdb_id.
    """
    session = session or build_session()
    rows = []
    for code in sorted({c for c in ligand_codes if c}):
        response = session.get(f"{PDBE_IN_PDB}/{code}", timeout=timeout)
        if response.status_code == 404:
            continue
        response.raise_for_status()
        payload = response.json() or {}
        pdb_ids = sorted({str(p).lower() for p in payload.get(code, [])})
        for pdb_id in pdb_ids[:max_per_ligand]:
            rows.append({"ligand_code": code, "pdb_id": pdb_id})
    return pd.DataFrame(rows, columns=["ligand_code", "pdb_id"])


def _row_hash(record: dict) -> str:
    encoded = json.dumps(record, sort_keys=True, default=str).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def land_pdbe(df: pd.DataFrame, chembl_version: str, raw_dir: Path | str = RAW_DIR) -> Path:
    """Stamp provenance and write raw_pdbe_structures.parquet (empty df is fine).

    The raw grain mirrors UniChem: keyed by the ligand code it resolves, so the
    code column is landed as ``xref_id`` and renamed back in staging.
    """
    raw_dir = Path(raw_dir)
    raw_dir.mkdir(parents=True, exist_ok=True)
    out = df.rename(columns={"ligand_code": "xref_id"}).copy()
    out["_fetch_ts"] = datetime.now(timezone.utc).isoformat()
    out["_source_endpoint"] = "pdbe/compound/in_pdb"
    out["_chembl_version"] = str(chembl_version)
    out["_row_hash"] = [_row_hash(r) for r in out.to_dict("records")]
    out_path = raw_dir / "raw_pdbe_structures.parquet"
    out.to_parquet(out_path, index=False)
    return out_path


def main() -> None:
    from extract.config import load_config

    config = load_config()
    refs = ligand_codes_from_raw()
    structures = resolve_pdb_entries(refs["ligand_code"].tolist())
    out = land_pdbe(structures, config.chembl_version)
    print(
        f"pdbe structures: {len(structures)} rows "
        f"from {refs['ligand_code'].nunique()} ligand code(s) -> {out}"
    )


if __name__ == "__main__":
    main()
