"""Run the M1 extract: fetch the scoped ChEMBL slice and land it to raw/.

Usage:
    python -m extract.run
"""

from extract.config import ExtractConfig, load_config
from extract.extract_activities import extract_activities
from extract.extract_assays import extract_assays
from extract.extract_cell_activities import extract_cell_activities
from extract.extract_molecules import extract_molecules
from extract.extract_targets import extract_targets
from extract.load_raw import load_raw


def _ids(records: list[dict], key: str) -> set[str]:
    return {r[key] for r in records if r.get(key)}


def run(config: ExtractConfig | None = None) -> None:
    """Extract activities, derive dimension IDs, and land every entity to raw/."""
    config = config or load_config()
    version = config.chembl_version

    activities = extract_activities(config)
    load_raw("activities", activities, version, endpoint="activity")
    print(f"activities: {len(activities)} rows")

    molecule_ids = _ids(activities, "molecule_chembl_id")
    target_ids = _ids(activities, "target_chembl_id")
    assay_ids = _ids(activities, "assay_chembl_id")

    for entity, endpoint, records in (
        ("molecules", "molecule", extract_molecules(molecule_ids)),
        ("targets", "target", extract_targets(target_ids)),
        ("assays", "assay", extract_assays(assay_ids)),
    ):
        load_raw(entity, records, version, endpoint=endpoint)
        print(f"{entity}: {len(records)} rows")

    # Cellular cytotoxicity for the reference payloads: a separate, compound-based
    # pull (no pchembl on cellular GI50/IC50), landed into its own raw entity so it
    # never enters the target-scoped fact_activity path.
    if config.reference_payloads and config.cell_activity:
        cell_activities = extract_cell_activities(config)
        load_raw("cell_activities", cell_activities, version, endpoint="activity")
        print(f"cell_activities: {len(cell_activities)} rows")


if __name__ == "__main__":
    run()
