"""M0: validate the target-set configuration is present and well-formed."""

from pathlib import Path

import yaml

CONFIG = Path(__file__).resolve().parents[1] / "config" / "target_set.yml"


def test_target_set_exists():
    assert CONFIG.exists(), "config/target_set.yml is missing"


def test_target_set_is_valid_yaml():
    data = yaml.safe_load(CONFIG.read_text())
    assert isinstance(data, dict)
    assert "chembl_version" in data
    assert isinstance(data.get("targets"), list) and data["targets"]
    assert isinstance(data.get("activity"), dict)


def test_activity_block_shape():
    data = yaml.safe_load(CONFIG.read_text())
    activity = data["activity"]
    assert isinstance(activity.get("standard_types"), list)
    assert "require_pchembl" in activity
    assert "min_confidence_score" in activity
