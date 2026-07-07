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
class ExtractConfig:
    chembl_version: str
    organism: str
    targets: list[Target]
    activity: ActivityFilter

    @property
    def target_ids(self) -> list[str]:
        return [t.chembl_id for t in self.targets]

    @property
    def payload_class_map(self) -> dict[str, str | None]:
        """Map each target ChEMBL id to its ADC-payload mechanism class."""
        return {t.chembl_id: t.payload_class for t in self.targets}


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
    return ExtractConfig(
        chembl_version=str(data["chembl_version"]),
        organism=str(data["organism"]),
        targets=targets,
        activity=activity,
    )
