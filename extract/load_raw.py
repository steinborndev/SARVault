"""Land fetched ChEMBL entities to the raw layer as Parquet.

Stamps provenance metadata on every row (_fetch_ts, _source_endpoint,
_chembl_version, _row_hash). Implemented in milestone M1 (feat/extract-raw).
"""


def load_raw(entity: str) -> None:
    """Write a fetched entity to the raw layer with provenance metadata."""
    raise NotImplementedError("Implemented in M1 (feat/extract-raw).")
