"""Pandera schemas validating raw ChEMBL entities before they are landed.

The raw layer is intentionally permissive (most filtering happens in staging),
so these schemas assert structural essentials: the key identifier columns are
present and the provenance stamp is complete and non-null.
"""

from pandera.pandas import Column, DataFrameSchema

_PROVENANCE = {
    "_fetch_ts": Column(str, nullable=False),
    "_source_endpoint": Column(str, nullable=False),
    "_chembl_version": Column(str, nullable=False),
    "_row_hash": Column(str, nullable=False),
}


def _schema(columns: dict) -> DataFrameSchema:
    return DataFrameSchema({**columns, **_PROVENANCE}, coerce=True, strict=False)


ACTIVITIES_SCHEMA = _schema(
    {
        "molecule_chembl_id": Column(str, nullable=True),
        "target_chembl_id": Column(str, nullable=True),
        "assay_chembl_id": Column(str, nullable=True),
        "standard_type": Column(str, nullable=True),
        "pchembl_value": Column(float, nullable=True),
    }
)
MOLECULES_SCHEMA = _schema({"molecule_chembl_id": Column(str, nullable=True)})
TARGETS_SCHEMA = _schema({"target_chembl_id": Column(str, nullable=True)})
ASSAYS_SCHEMA = _schema({"assay_chembl_id": Column(str, nullable=True)})

# The cheminfo stage is an RDKit compute, not a ChEMBL fetch, so it carries its own
# provenance (source = rdkit, plus the pinned RDKit version) and stricter columns:
# every landed row must have an identity, a fingerprint and an atom count.
_CHEMINFO_PROVENANCE = {
    "_fetch_ts": Column(str, nullable=False),
    "_source": Column(str, nullable=False),
    "_rdkit_version": Column(str, nullable=False),
    "_row_hash": Column(str, nullable=False),
}
COMPOUND_CHEMINFO_SCHEMA = DataFrameSchema(
    {
        "molecule_chembl_id": Column(str, nullable=False),
        "ecfp4_hex": Column(str, nullable=False),
        "ecfp4_onbits": Column(object, nullable=False),
        "n_onbits": Column(int, nullable=False),
        "heavy_atom_count": Column(int, nullable=False),
        "murcko_scaffold_smiles": Column(str, nullable=True),
        "murcko_generic_smiles": Column(str, nullable=True),
        **_CHEMINFO_PROVENANCE,
    },
    coerce=True,
    strict=False,
)

SCHEMAS = {
    "activities": ACTIVITIES_SCHEMA,
    "molecules": MOLECULES_SCHEMA,
    "targets": TARGETS_SCHEMA,
    "assays": ASSAYS_SCHEMA,
    "compound_cheminfo": COMPOUND_CHEMINFO_SCHEMA,
}


def validate(entity: str, df):
    """Validate a raw entity DataFrame against its schema; returns the coerced df."""
    return SCHEMAS[entity].validate(df)
