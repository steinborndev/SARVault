"""M1: unit tests for the extract layer (no live API calls)."""

import pandas as pd
import pytest

from extract import load_raw as load_raw_mod
from extract.chembl_client import chunked, fetch_all
from extract.config import load_config
from extract.extract_activities import build_activity_params
from extract.load_raw import load_raw, row_hash
from validation.schemas import validate


class _FakeResponse:
    def __init__(self, payload):
        self._payload = payload

    def raise_for_status(self):
        return None

    def json(self):
        return self._payload


class _FakeSession:
    """Returns queued pages in order, recording the offsets requested."""

    def __init__(self, pages):
        self._pages = list(pages)
        self.offsets = []
        self.calls = 0

    def get(self, url, params=None, timeout=None):
        self.offsets.append(params.get("offset"))
        page = self._pages[self.calls]
        self.calls += 1
        return _FakeResponse(page)


def test_fetch_all_paginates_until_next_is_null():
    pages = [
        {"activities": [{"id": 1}, {"id": 2}], "page_meta": {"next": "/x?offset=2"}},
        {"activities": [{"id": 3}], "page_meta": {"next": None}},
    ]
    session = _FakeSession(pages)
    records = list(fetch_all("activity", {}, session=session, page_limit=2))
    assert [r["id"] for r in records] == [1, 2, 3]
    assert session.offsets == [0, 2]


def test_chunked_splits_evenly():
    assert list(chunked(["a", "b", "c", "d", "e"], 2)) == [["a", "b"], ["c", "d"], ["e"]]


def test_config_has_four_verified_targets():
    config = load_config()
    assert config.target_ids == [
        "CHEMBL2095182",
        "CHEMBL3832942",
        "CHEMBL1781",
        "CHEMBL1806",
    ]


def test_build_activity_params_encodes_filters():
    params = build_activity_params(load_config())
    assert "CHEMBL1781" in params["target_chembl_id__in"]
    assert "IC50" in params["standard_type__in"]
    assert params["pchembl_value__isnull"] == "false"


def test_row_hash_is_stable_and_order_independent():
    assert row_hash({"x": 1, "y": 2}) == row_hash({"y": 2, "x": 1})


def test_load_raw_writes_parquet_with_provenance(tmp_path):
    records = [
        {"molecule_chembl_id": "CHEMBL1", "synonyms": ["foo", "bar"]},
        {"molecule_chembl_id": "CHEMBL2", "synonyms": []},
    ]
    out = load_raw("molecules", records, "36", endpoint="molecule", raw_dir=tmp_path)
    df = pd.read_parquet(out)
    assert len(df) == 2
    for col in load_raw_mod.PROVENANCE_COLUMNS:
        assert col in df.columns
        assert df[col].notna().all()
    assert (df["_source_endpoint"] == "molecule").all()
    assert (df["_chembl_version"] == "36").all()
    assert isinstance(df["synonyms"].iloc[0], str)  # nested list serialized to JSON


def test_pandera_schema_rejects_missing_provenance():
    valid = pd.DataFrame(
        {
            "molecule_chembl_id": ["CHEMBL1"],
            "_fetch_ts": ["2026-01-01T00:00:00+00:00"],
            "_source_endpoint": ["molecule"],
            "_chembl_version": ["36"],
            "_row_hash": ["abc"],
        }
    )
    validate("molecules", valid)  # should not raise
    with pytest.raises(Exception):
        validate("molecules", valid.drop(columns=["_row_hash"]))


def test_load_unichem_mappings_filters_and_labels(tmp_path):
    """UniChem bulk reader: header dropped via id-join, filtered to our compounds."""
    import gzip

    from extract import unichem

    path = tmp_path / "src1src22.txt.gz"
    with gzip.open(path, "wt") as fh:
        fh.write("From src:'1'\tTo src:'22'\n")  # header -> dropped (not in ids)
        fh.write("CHEMBL25\t2244\n")
        fh.write("CHEMBL999\t9999\n")  # out of scope -> filtered out
    df = unichem.load_unichem_mappings(tmp_path, ["CHEMBL25"], {"pubchem": 22})
    assert list(df["molecule_chembl_id"]) == ["CHEMBL25"]
    assert list(df["source"]) == ["pubchem"]
    assert list(df["xref_id"]) == ["2244"]


def test_land_unichem_writes_provenance(tmp_path):
    import pandas as pd

    from extract import unichem

    df = pd.DataFrame(
        {"molecule_chembl_id": ["CHEMBL25"], "source": ["pubchem"], "xref_id": ["2244"]}
    )
    out = unichem.land_unichem(df, "36", raw_dir=tmp_path)
    landed = pd.read_parquet(out)
    for col in unichem.PROVENANCE_COLUMNS:
        assert col in landed.columns


class _PdbeResponse:
    def __init__(self, payload, status_code=200):
        self._payload = payload
        self.status_code = status_code

    def raise_for_status(self):
        if self.status_code >= 400:
            raise AssertionError("raise_for_status should not fire on the 404 path")

    def json(self):
        return self._payload


class _PdbeSession:
    """Serves a fixed het-code -> PDB-id-list mapping; 404s unknown codes."""

    def __init__(self, mapping):
        self.mapping = mapping

    def get(self, url, timeout=None):
        code = url.rstrip("/").split("/")[-1]
        if code not in self.mapping:
            return _PdbeResponse({}, status_code=404)
        return _PdbeResponse({code: self.mapping[code]})


def _in_pdb(pdb_ids):
    """Build the dict-per-entry shape the current PDBe in_pdb API returns."""
    return [
        {"pdb_id": p, "ligand_type": "B", "ligand_total_number": 1, "ligand_distinct_number": 1}
        for p in pdb_ids
    ]


def test_resolve_pdb_entries_extracts_pdb_id_from_dict_items():
    from extract import pdbe

    session = _PdbeSession({"DM5": _in_pdb(["198D", "1D12", "198d"])})
    df = pdbe.resolve_pdb_entries(["DM5"], session=session)
    assert set(df["pdb_id"]) == {"198d", "1d12"}  # deduped, lowercased
    assert all("{" not in pid for pid in df["pdb_id"])  # no stringified dict leaks


def test_resolve_pdb_entries_accepts_legacy_string_items():
    from extract import pdbe

    session = _PdbeSession({"AIN": ["2QQT", "1OXR", "1oxr"]})
    df = pdbe.resolve_pdb_entries(["AIN"], session=session)
    assert list(df["pdb_id"]) == ["1oxr", "2qqt"]


def test_resolve_pdb_entries_caps_per_ligand():
    from extract import pdbe

    session = _PdbeSession({"AIN": _in_pdb([f"{i:04d}" for i in range(100)])})
    df = pdbe.resolve_pdb_entries(["AIN"], session=session, max_per_ligand=10)
    assert len(df) == 10


def test_resolve_pdb_entries_skips_missing_component():
    from extract import pdbe

    df = pdbe.resolve_pdb_entries(["ZZZ"], session=_PdbeSession({}))
    assert df.empty


def test_land_pdbe_lands_xref_id_grain_with_provenance(tmp_path):
    from extract import pdbe

    df = pd.DataFrame({"ligand_code": ["AIN"], "pdb_id": ["1oxr"]})
    out = pdbe.land_pdbe(df, "36", raw_dir=tmp_path)
    landed = pd.read_parquet(out)
    assert "xref_id" in landed.columns and "pdb_id" in landed.columns
    assert "ligand_code" not in landed.columns  # renamed to xref_id at landing
    for col in pdbe.PROVENANCE_COLUMNS:
        assert col in landed.columns and landed[col].notna().all()


def test_ligand_codes_from_raw_reads_unichem_pdbe(tmp_path):
    from extract import pdbe

    pd.DataFrame(
        {"molecule_chembl_id": ["CHEMBLM1"], "source": ["pdbe"], "xref_id": ["AIN"]}
    ).to_parquet(tmp_path / "raw_xref_unichem.parquet", index=False)
    refs = pdbe.ligand_codes_from_raw(raw_dir=tmp_path)
    assert list(refs["ligand_code"]) == ["AIN"]


class _PdbeMetaSession:
    """Serves summary/experiment records per id via GET, keyed by endpoint URL."""

    def __init__(self, summaries, experiments):
        self.summaries = summaries
        self.experiments = experiments
        self.gets = []

    def get(self, url, timeout=None):
        pid = url.rstrip("/").split("/")[-1]
        self.gets.append(url)
        src = self.summaries if "summary" in url else self.experiments
        if pid not in src:
            return _PdbeResponse({}, status_code=404)
        # PDBe wraps each entry's record in a one-element list keyed by pdb id.
        return _PdbeResponse({pid: [src[pid]]})


def test_fetch_pdb_metadata_merges_summary_and_experiment():
    from extract import pdbe

    session = _PdbeMetaSession(
        summaries={
            "1oxr": {
                "title": "COX-2 with aspirin",
                "experimental_method": ["X-ray diffraction"],
                "release_date": "20040115",
            }
        },
        experiments={"1oxr": {"resolution": 2.0}},
    )
    df = pdbe.fetch_pdb_metadata(["1OXR"], session=session)
    row = df.iloc[0]
    assert row["pdb_id"] == "1oxr"  # lowercased
    assert row["title"] == "COX-2 with aspirin"
    assert row["method"] == "X-ray diffraction"
    assert row["year"] == "2004"
    assert row["resolution"] == 2.0


def test_fetch_pdb_metadata_tolerates_missing_entries():
    from extract import pdbe

    df = pdbe.fetch_pdb_metadata(["ZZZZ"], session=_PdbeMetaSession({}, {}))
    row = df.iloc[0]
    assert row["pdb_id"] == "zzzz"
    assert row["title"] is None
    assert row["method"] is None
    assert row["year"] is None
    assert row["resolution"] is None


def test_land_pdbe_summary_writes_provenance(tmp_path):
    from extract import pdbe

    df = pd.DataFrame(
        [
            {
                "pdb_id": "1oxr",
                "title": "t",
                "method": "X-ray diffraction",
                "year": "2004",
                "resolution": 2.0,
            }
        ]
    )
    out = pdbe.land_pdbe_summary(df, "36", raw_dir=tmp_path)
    landed = pd.read_parquet(out)
    assert {"pdb_id", "title", "method", "year", "resolution"}.issubset(landed.columns)
    for col in pdbe.PROVENANCE_COLUMNS:
        assert col in landed.columns and landed[col].notna().all()
