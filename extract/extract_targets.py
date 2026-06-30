"""Extract ChEMBL targets for the configured target set.

Reads config/target_set.yml, queries the ChEMBL API with pagination, and
returns raw records for the raw-load step. Implemented in milestone M1
(feat/extract-raw).
"""


def extract_targets():
    """Fetch targets records for the configured scope."""
    raise NotImplementedError("Implemented in M1 (feat/extract-raw).")
