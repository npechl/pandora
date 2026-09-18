from pathlib import Path

import httpx
import pytest

from pandora.ingestion.mmcif import (
    _pdbe_revision_date,
    _revision_date_from_cif_file,
    ingest_local_mmcif,
)
from pandora.ingestion.policy import load_policy
from pandora.ingestion.search import search_pdbe, search_rcsb

POLICY_PATH = (
    Path(__file__).parent.parent / "datasets" / "canonicalisation.yaml"
)


def test_load_policy():
    policy = load_policy(str(POLICY_PATH))

    assert policy.policy_id == "overview-remap"
    assert policy.identifier_rules.chain_id.strategy == "remap"
    assert policy.altloc_rules.strategy == "select_best_occupancy"
    assert policy.ligand_rules.keep_waters is False


def test_revision_date_from_cif_file_takes_latest_row(tmp_path):
    cif_path = tmp_path / "with_history.cif"
    cif_path.write_text(
        "data_TEST\n"
        "loop_\n"
        "_pdbx_audit_revision_history.ordinal\n"
        "_pdbx_audit_revision_history.revision_date\n"
        "1 2020-01-01\n"
        "2 2021-06-15\n"
    )

    assert _revision_date_from_cif_file(cif_path) == "2021-06-15"


def test_revision_date_from_cif_file_missing_category(tmp_path):
    cif_path = tmp_path / "no_history.cif"
    cif_path.write_text("data_TEST\n_entry.id TEST\n")

    assert _revision_date_from_cif_file(cif_path) is None


def test_pdbe_revision_date_normalizes_and_handles_failure(monkeypatch):
    class FakeResponse:
        def raise_for_status(self):
            pass

        def json(self):
            return {"104m": [{"revision_date": "20260812"}]}

    monkeypatch.setattr(httpx, "get", lambda *a, **k: FakeResponse())
    assert _pdbe_revision_date("104m") == "2026-08-12"

    def raise_network_error(*a, **k):
        raise httpx.ConnectError("no network")

    monkeypatch.setattr(httpx, "get", raise_network_error)
    assert _pdbe_revision_date("104m") is None


def test_ingest_local_mmcif_defaults_source_uri_to_path(tmp_path):
    cif_path = tmp_path / "with_history.cif"
    cif_path.write_text(
        "data_TEST\n"
        "loop_\n"
        "_pdbx_audit_revision_history.ordinal\n"
        "_pdbx_audit_revision_history.revision_date\n"
        "1 2020-01-01\n"
    )

    provenance = ingest_local_mmcif(cif_path)
    assert provenance.provider == "local"
    assert provenance.source_uri == str(cif_path)
    assert provenance.from_cache is False
    assert provenance.revision_date == "2020-01-01"

    labeled = ingest_local_mmcif(cif_path, source_uri="pdb_snapshot_2024-01-29")
    assert labeled.source_uri == "pdb_snapshot_2024-01-29"


def test_ingest_local_mmcif_missing_path_raises(tmp_path):
    with pytest.raises(ValueError):
        ingest_local_mmcif(tmp_path / "does_not_exist.cif")


class _FakeResponse:
    def __init__(self, status_code=200, json_data=None, text=""):
        self.status_code = status_code
        self._json_data = json_data
        self.text = text

    def raise_for_status(self):
        if self.status_code >= 400:
            raise httpx.HTTPStatusError("error", request=None, response=self)

    def json(self):
        return self._json_data


def test_search_rcsb_returns_identifiers(monkeypatch):
    captured = {}

    def fake_post(url, json=None, timeout=None):
        captured["json"] = json
        return _FakeResponse(
            json_data={"result_set": [{"identifier": "104M", "score": 1.0}]}
        )

    monkeypatch.setattr(httpx, "post", fake_post)
    ids = search_rcsb({"type": "terminal"}, rows=10, start=0)

    assert ids == ["104M"]
    assert captured["json"]["request_options"]["paginate"] == {
        "start": 0,
        "rows": 10,
    }


def test_search_rcsb_no_matches_returns_empty(monkeypatch):
    monkeypatch.setattr(
        httpx, "post", lambda *a, **k: _FakeResponse(status_code=204)
    )
    assert search_rcsb({"type": "terminal"}) == []


def test_search_rcsb_bad_request_raises(monkeypatch):
    monkeypatch.setattr(
        httpx,
        "post",
        lambda *a, **k: _FakeResponse(status_code=400, text="bad attribute"),
    )
    with pytest.raises(RuntimeError, match="bad attribute"):
        search_rcsb({"type": "terminal"})


def test_search_pdbe_returns_identifiers(monkeypatch):
    captured = {}

    def fake_get(url, params=None, timeout=None):
        captured["params"] = params
        return _FakeResponse(
            json_data={"response": {"docs": [{"pdb_id": "8wa8"}]}}
        )

    monkeypatch.setattr(httpx, "get", fake_get)
    ids = search_pdbe("resolution:[0 TO 1.5]", rows=5)

    assert ids == ["8wa8"]
    assert captured["params"]["q"] == "resolution:[0 TO 1.5]"
    assert captured["params"]["fl"] == "pdb_id"


def test_search_pdbe_dedupes_preserving_order(monkeypatch):
    monkeypatch.setattr(
        httpx,
        "get",
        lambda *a, **k: _FakeResponse(
            json_data={
                "response": {
                    "docs": [
                        {"pdb_id": "7oju"},
                        {"pdb_id": "7o34"},
                        {"pdb_id": "7oju"},
                    ]
                }
            }
        ),
    )
    assert search_pdbe("resolution:[0 TO 1.5]") == ["7oju", "7o34"]


def test_search_pdbe_bad_request_raises(monkeypatch):
    monkeypatch.setattr(
        httpx,
        "get",
        lambda *a, **k: _FakeResponse(status_code=500, text="solr error"),
    )
    with pytest.raises(RuntimeError, match="solr error"):
        search_pdbe("bad query")
