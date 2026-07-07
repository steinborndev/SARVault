"""F3.1: payload_class dimension.

The single source of truth is config/target_set.yml. dim_target sources the
mapping from the payload_classes dbt var (dbt/dbt_project.yml), mirrored from the
config the same way standard_types is. These tests guard that mirror against
drift and that the extract config exposes the mapping.
"""

from pathlib import Path

import yaml

from extract.config import load_config

REPO = Path(__file__).resolve().parents[1]
DBT_PROJECT = REPO / "dbt" / "dbt_project.yml"
KNOWN_CLASSES = {"tubulin_inhibitor", "topo1_inhibitor", "topo2_inhibitor"}


def test_every_target_has_a_known_payload_class():
    config = load_config()
    for target in config.targets:
        assert target.payload_class in KNOWN_CLASSES, (
            f"{target.chembl_id} has payload_class {target.payload_class!r}"
        )


def test_payload_class_map_matches_targets():
    config = load_config()
    assert config.payload_class_map == {
        t.chembl_id: t.payload_class for t in config.targets
    }


def test_dbt_var_mirrors_config():
    config = load_config()
    expected = {t.chembl_id: t.payload_class for t in config.targets}

    project = yaml.safe_load(DBT_PROJECT.read_text())
    dbt_map = {
        m["target_chembl_id"]: m["payload_class"]
        for m in project["vars"]["payload_classes"]
    }
    assert dbt_map == expected, (
        "dbt payload_classes var is out of sync with config/target_set.yml"
    )
