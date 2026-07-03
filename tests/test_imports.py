"""M0: the scaffold packages import cleanly without heavy extras installed."""

import importlib

import pytest

MODULES = [
    "extract",
    "validation",
    "orchestration",
    "extract.chembl_client",
    "extract.load_raw",
    "extract.extract_molecules",
    "extract.extract_targets",
    "extract.extract_assays",
    "extract.extract_activities",
    "validation.schemas",
    # orchestration.definitions is intentionally excluded: it now requires the
    # Dagster extras and a parsed dbt manifest, so it is exercised by its own
    # suite (tests/test_orchestration.py) rather than this light-import check.
]


@pytest.mark.parametrize("module", MODULES)
def test_module_imports(module):
    importlib.import_module(module)
