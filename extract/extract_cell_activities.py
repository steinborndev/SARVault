"""Extract cellular cytotoxicity activities for the reference ADC payloads.

Cellular GI50/IC50 (assay_type F) is reported with a raw standard_value but no
pchembl_value, so it never enters the target-scoped fact_activity path. This
compound-based pull captures it into a separate lineage (raw_cell_activities ->
mart_compound_cytotoxicity), keyed on the curated reference payloads and stamped
with each payload's mechanism class so cytotoxicity can be grouped by class.
"""

import requests

from extract.chembl_client import chunked, fetch_all
from extract.config import ExtractConfig, load_config

_CHUNK = 20  # small: ChEMBL rejects over-long __in URLs


def build_cell_activity_params(ids: list[str], config: ExtractConfig) -> dict:
    """Build the ChEMBL /activity query for cellular readouts of the given molecules."""
    return {
        "molecule_chembl_id__in": ",".join(ids),
        "standard_type__in": ",".join(config.cell_activity.standard_types),
        "assay_type": config.cell_activity.assay_type,
    }


def extract_cell_activities(
    config: ExtractConfig | None = None,
    session: requests.Session | None = None,
) -> list[dict]:
    """Fetch cellular cytotoxicity records for the reference payloads.

    Each record is stamped with the payload's class and display name (from config)
    so the downstream mart can group cytotoxicity by mechanism class without a join.
    """
    config = config or load_config()
    class_map = config.reference_payload_class_map
    name_map = config.reference_payload_name_map
    ids = sorted(class_map)
    records: list[dict] = []
    for chunk in chunked(ids, _CHUNK):
        params = build_cell_activity_params(chunk, config)
        for rec in fetch_all("activity", params, session=session):
            mid = rec.get("molecule_chembl_id")
            rec["payload_class"] = class_map.get(mid)
            rec["reference_name"] = name_map.get(mid)
            records.append(rec)
    return records
