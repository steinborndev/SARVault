"""Compute cheminformatics descriptors and land them to the raw layer.

For every compound with a parseable canonical SMILES this stage derives, with
RDKit, an ECFP4 (Morgan radius-2, 2048-bit) fingerprint and its Bemis-Murcko
scaffold (both the concrete scaffold and the generic graph framework), then lands
one row per compound to ``raw_compound_cheminfo.parquet`` with provenance columns.

It is a deterministic *compute* stage, not an API extract: given a pinned RDKit
version the output is a pure function of the input SMILES. It runs after the
molecules are landed (it reads ``raw_molecules.parquet``) and feeds the dbt
``stg_compound_cheminfo`` model, ``dim_scaffold`` and ``mart_compound_fingerprint``.

Usage:
    python -m extract.cheminfo
"""

import hashlib
from datetime import datetime, timezone
from pathlib import Path

import duckdb
import pandas as pd

RAW_DIR = Path(__file__).resolve().parents[1] / "raw"

# ECFP4 == Morgan fingerprint, radius 2. 2048 bits balances collision rate and size.
FP_RADIUS = 2
FP_NBITS = 2048

PROVENANCE_COLUMNS = ["_fetch_ts", "_source", "_rdkit_version", "_row_hash"]


def rdkit_version() -> str:
    """Pinned RDKit version, recorded per row for reproducibility."""
    import rdkit

    return rdkit.__version__


def _row_hash(molecule_chembl_id: str, smiles: str) -> str:
    """Stable SHA-256 over the identity + structure the descriptors are derived from."""
    encoded = f"{molecule_chembl_id}\t{smiles}".encode()
    return hashlib.sha256(encoded).hexdigest()


def compute_cheminfo(smiles: str) -> dict | None:
    """Derive ECFP4 + Murcko descriptors for one SMILES, or None if it won't parse.

    Returns a dict with:
      - ``ecfp4_hex``            2048-bit Morgan fingerprint as a 512-char hex string
      - ``n_onbits``            number of set bits (popcount)
      - ``heavy_atom_count``    non-hydrogen atom count
      - ``murcko_scaffold_smiles``  Bemis-Murcko scaffold (None for acyclic molecules)
      - ``murcko_generic_smiles``   generic (element/bond-agnostic) framework, or None
    """
    if not smiles:
        return None
    from rdkit import Chem
    from rdkit.Chem import rdFingerprintGenerator
    from rdkit.Chem.Scaffolds import MurckoScaffold

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    gen = rdFingerprintGenerator.GetMorganGenerator(radius=FP_RADIUS, fpSize=FP_NBITS)
    fp = gen.GetFingerprint(mol)
    bitstring = fp.ToBitString()  # deterministic length-FP_NBITS '0'/'1' string
    ecfp4_hex = f"{int(bitstring, 2):0{FP_NBITS // 4}x}"

    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    scaffold_smiles = Chem.MolToSmiles(scaffold) if scaffold.GetNumAtoms() else None
    generic_smiles = None
    if scaffold_smiles:
        try:
            generic = MurckoScaffold.MakeScaffoldGeneric(scaffold)
            generic_smiles = Chem.MolToSmiles(generic)
        except Exception:
            # A few pathological scaffolds can't be genericised; keep the concrete one.
            generic_smiles = None

    return {
        "ecfp4_hex": ecfp4_hex,
        "n_onbits": int(fp.GetNumOnBits()),
        "heavy_atom_count": int(mol.GetNumHeavyAtoms()),
        "murcko_scaffold_smiles": scaffold_smiles,
        "murcko_generic_smiles": generic_smiles,
    }


def _read_molecule_smiles(raw_dir: Path) -> pd.DataFrame:
    """Read (molecule_chembl_id, canonical_smiles) from the landed raw molecules."""
    path = raw_dir / "raw_molecules.parquet"
    return duckdb.sql(
        f"""
        select
            molecule_chembl_id,
            "molecule_structures.canonical_smiles" as canonical_smiles
        from read_parquet('{path}')
        """
    ).df()


def build_cheminfo(raw_dir: Path | str = RAW_DIR) -> pd.DataFrame:
    """Compute descriptors for every parseable molecule; return a stamped DataFrame.

    Molecules with a null or unparseable SMILES are dropped and counted; the count
    is printed so a real run surfaces how much of the set had usable structures.
    """
    raw_dir = Path(raw_dir)
    molecules = _read_molecule_smiles(raw_dir)

    rows: list[dict] = []
    skipped = 0
    for rec in molecules.itertuples(index=False):
        info = compute_cheminfo(rec.canonical_smiles)
        if info is None:
            skipped += 1
            continue
        rows.append(
            {
                "molecule_chembl_id": rec.molecule_chembl_id,
                "canonical_smiles": rec.canonical_smiles,
                **info,
                "_row_hash": _row_hash(rec.molecule_chembl_id, rec.canonical_smiles),
            }
        )

    print(f"cheminfo: {len(rows)} computed, {skipped} skipped (null/unparseable SMILES)")

    df = pd.DataFrame(rows)
    df["_fetch_ts"] = datetime.now(timezone.utc).isoformat()
    df["_source"] = "rdkit"
    df["_rdkit_version"] = rdkit_version()
    return df


def land_cheminfo(
    df: pd.DataFrame, raw_dir: Path | str = RAW_DIR, validate_schema: bool = True
) -> Path:
    """Validate (Pandera) and write the descriptor frame to raw_compound_cheminfo.parquet."""
    if validate_schema:
        from validation.schemas import validate

        df = validate("compound_cheminfo", df)
    raw_dir = Path(raw_dir)
    raw_dir.mkdir(parents=True, exist_ok=True)
    out_path = raw_dir / "raw_compound_cheminfo.parquet"
    df.to_parquet(out_path, index=False)
    return out_path


def main(raw_dir: Path | str = RAW_DIR) -> Path:
    df = build_cheminfo(raw_dir)
    out = land_cheminfo(df, raw_dir)
    print(f"cheminfo: {len(df)} rows -> {out}")
    return out


if __name__ == "__main__":
    main()
