"""Load and validate the target-set configuration (config/target_set.yml)."""

from dataclasses import dataclass
from pathlib import Path

import yaml

DEFAULT_CONFIG_PATH = Path(__file__).resolve().parents[1] / "config" / "target_set.yml"


@dataclass(frozen=True)
class Target:
    chembl_id: str
    label: str
    payload_class: str | None = None


@dataclass(frozen=True)
class ActivityFilter:
    standard_types: list[str]
    require_pchembl: bool
    min_confidence_score: int


@dataclass(frozen=True)
class ReferencePayload:
    molecule_chembl_id: str
    name: str
    payload_class: str


@dataclass(frozen=True)
class CellActivityFilter:
    standard_types: list[str]
    assay_type: str


@dataclass(frozen=True)
class ExtractConfig:
    chembl_version: str
    organism: str
    targets: list[Target]
    activity: ActivityFilter
    reference_payloads: list[ReferencePayload] = ()
    cell_activity: CellActivityFilter | None = None

    @property
    def target_ids(self) -> list[str]:
        return [t.chembl_id for t in self.targets]

    @property
    def payload_class_map(self) -> dict[str, str | None]:
        """Map each target ChEMBL id to its ADC-payload mechanism class."""
        return {t.chembl_id: t.payload_class for t in self.targets}

    @property
    def reference_payload_ids(self) -> list[str]:
        return [p.molecule_chembl_id for p in self.reference_payloads]

    @property
    def reference_payload_class_map(self) -> dict[str, str]:
        """Map each reference-payload molecule id to its payload class."""
        return {p.molecule_chembl_id: p.payload_class for p in self.reference_payloads}

    @property
    def reference_payload_name_map(self) -> dict[str, str]:
        """Map each reference-payload molecule id to its display name."""
        return {p.molecule_chembl_id: p.name for p in self.reference_payloads}


def load_config(path: Path | str = DEFAULT_CONFIG_PATH) -> ExtractConfig:
    """Parse the YAML target set into a typed, frozen config object."""
    data = yaml.safe_load(Path(path).read_text())
    targets = [
        Target(chembl_id=t["chembl_id"], label=t["label"], payload_class=t.get("payload_class"))
        for t in data["targets"]
    ]
    activity = ActivityFilter(
        standard_types=list(data["activity"]["standard_types"]),
        require_pchembl=bool(data["activity"]["require_pchembl"]),
        min_confidence_score=int(data["activity"]["min_confidence_score"]),
    )
    reference_payloads = tuple(
        ReferencePayload(
            molecule_chembl_id=p["molecule_chembl_id"],
            name=p["name"],
            payload_class=p["payload_class"],
        )
        for p in data.get("reference_payloads", [])
    )
    cell = data.get("cell_activity")
    cell_activity = (
        CellActivityFilter(
            standard_types=list(cell["standard_types"]),
            assay_type=str(cell["assay_type"]),
        )
        if cell
        else None
    )
    return ExtractConfig(
        chembl_version=str(data["chembl_version"]),
        organism=str(data["organism"]),
        targets=targets,
        activity=activity,
        reference_payloads=reference_payloads,
        cell_activity=cell_activity,
    )
