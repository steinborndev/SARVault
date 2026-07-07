"""F3.2: cellular cytotoxicity extract for the reference ADC payloads."""

from extract.config import load_config
from extract.extract_cell_activities import build_cell_activity_params

KNOWN_CLASSES = {"tubulin_inhibitor", "topo1_inhibitor"}


def test_reference_payloads_parsed_with_classes():
    config = load_config()
    assert config.reference_payloads, "config has no reference_payloads"
    assert set(config.reference_payload_class_map.values()) <= KNOWN_CLASSES
    for p in config.reference_payloads:
        assert p.molecule_chembl_id.startswith("CHEMBL")
        assert p.name and p.payload_class in KNOWN_CLASSES


def test_cell_activity_params_shape():
    config = load_config()
    params = build_cell_activity_params(["CHEMBL65", "CHEMBL84"], config)
    assert params["molecule_chembl_id__in"] == "CHEMBL65,CHEMBL84"
    assert params["standard_type__in"] == "GI50,IC50"
    assert params["assay_type"] == "F"
