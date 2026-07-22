#! /usr/bin/env python3
"""Offline tests for NCBI `datasets` request hygiene and chunking.

Covers the two halves of the initialization failure where every
`dataset_report` POST was refused before it left the machine:

    net/http: invalid header field value for "Api-Key"
    net/http: invalid header field value for "X-Datasets-Client-Cmd"

An api key travels as an HTTP header and `datasets` echoes its own argv into a
second one, so a key holding a newline poisons both. The chunking half keeps a
16,000-accession initialization from riding on a single request.
"""
import inspect
import json
import zipfile

import pytest

from mycotools.lib.dbtools import clean_api_key
from mycotools.mtdb import update as mod


# --------------------------------------------------------------------------- #
# api key sanitation
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize(
    "raw,expect",
    [
        ("0123456789abcdef", "0123456789abcdef"),  # untouched
        ("0123456789abcdef\n", "0123456789abcdef"),  # pasted with its newline
        ("0123456789abcdef\r\n", "0123456789abcdef"),  # CRLF store
        ("  0123456789abcdef  ", "0123456789abcdef"),  # padded
        ("0123\n4567", "01234567"),  # embedded control character
    ],
)
def test_clean_api_key_yields_legal_header(raw, expect):
    assert clean_api_key(raw) == expect


def test_clean_api_key_preserves_falsy_types():
    """`None` and `""` mean different things to callers that persist the key."""
    assert clean_api_key(None) is None
    assert clean_api_key("") == ""


@pytest.mark.parametrize("raw", ["k\n", "k\r", "k\x00", "k\x7f", "k\x1f"])
def test_cleaned_keys_hold_no_control_characters(raw):
    """Go rejects any byte below 0x20 (bar space/tab) and 0x7f in a header."""
    cleaned = clean_api_key(raw)
    assert not any(ord(c) < 0x20 or ord(c) == 0x7F for c in cleaned)


def test_run_datasets_mutes_stdout_but_reports_failures(monkeypatch, tmp_path, caplog):
    """Silenced so it cannot fight the chunk bar -- but a failure still speaks."""
    import subprocess

    from mycotools.download import ncbi

    seen = {}

    def fake_run(cmd, **kwargs):
        seen.update(kwargs)
        return subprocess.CompletedProcess(
            cmd, 1, stdout=None, stderr="Collecting 5 records [--]\rdone\nError: nope\n"
        )

    monkeypatch.setattr(ncbi.subprocess, "run", fake_run)
    with caplog.at_level("ERROR"):
        assert ncbi.run_datasets(None, "accs.txt", str(tmp_path) + "/", True) == 1

    assert seen["stdout"] == subprocess.DEVNULL
    assert "Error: nope" in caplog.text
    # the progress bar's redraw frames are collapsed, not replayed
    assert "[--]" not in caplog.text


def test_mute_stderr_keeps_datasets_off_the_terminal(monkeypatch, tmp_path, caplog):
    """A caller drawing its own bar silences datasets without losing the text."""
    import subprocess

    from mycotools.download import ncbi

    def fake_run(cmd, **kwargs):
        return subprocess.CompletedProcess(cmd, 1, stdout=None, stderr="Error: nope\n")

    monkeypatch.setattr(ncbi.subprocess, "run", fake_run)
    with caplog.at_level("ERROR"):
        ncbi.run_datasets(
            None, "accs.txt", str(tmp_path) + "/", True, mute_stderr=True
        )
    assert caplog.text == ""

    caplog.clear()
    with caplog.at_level("DEBUG"):
        ncbi.run_datasets(
            None, "accs.txt", str(tmp_path) + "/", True, mute_stderr=True
        )
    assert "Error: nope" in caplog.text  # retained for -v, not shown by default


def test_chunks_mute_datasets_stderr(tmp_path, monkeypatch):
    """dwnld_data_reports always mutes; there is no CLI knob for it."""
    from mycotools.download import ncbi

    seen = []

    def fake_run(include, accs_file, output_path, annotated, **kwargs):
        from pathlib import Path

        seen.append(kwargs)
        with open(accs_file, "r") as raw:
            _fake_report(Path(output_path), [x for x in raw.read().split("\n") if x])
        return 0

    monkeypatch.setattr(mod, "run_datasets", fake_run)
    mod.dwnld_data_reports(
        [f"GCA_{i:09d}.1" for i in range(20)], str(tmp_path) + "/", chunk=10
    )

    assert seen and all(kw["mute_stderr"] is True for kw in seen)
    assert all(kw["verbose"] is False for kw in seen)
    # hardcoded rather than a parameter, so there is nothing for the CLI to
    # thread through and no way for a user to turn the bar back into noise
    assert "mute_stderr" not in inspect.signature(mod.dwnld_data_reports).parameters


def test_run_datasets_never_passes_an_illegal_key(monkeypatch, tmp_path):
    """The argv `datasets` reflects into X-Datasets-Client-Cmd must be clean."""
    from mycotools.download import ncbi

    import subprocess

    seen = {}

    def fake_run(cmd, **kwargs):
        seen["cmd"] = cmd
        return subprocess.CompletedProcess(cmd, 0, stdout=None, stderr="")

    monkeypatch.setattr(ncbi.subprocess, "run", fake_run)
    ncbi.run_datasets(
        None, str(tmp_path / "accs.txt"), str(tmp_path) + "/", True, api="APIKEY\n"
    )
    assert "APIKEY" in seen["cmd"]
    assert not any("\n" in arg or "\r" in arg for arg in seen["cmd"])


# --------------------------------------------------------------------------- #
# chunked data report acquisition
# --------------------------------------------------------------------------- #
def _fake_report(out_dir, accs):
    """Write the ncbi_dataset.zip `datasets` would leave for a chunk, so the
    caller's own unzip/parse path is the one under test."""
    report = "\n".join(
        json.dumps(
            {
                "accession": acc,
                "organism": {"organismName": f"Aspergillus {acc.lower()} X1"},
                "assemblyInfo": {"assemblyName": acc},
                "annotationInfo": {"provider": "test"},
            }
        )
        for acc in accs
    )
    with zipfile.ZipFile(out_dir / "ncbi_dataset.zip", "w") as zf:
        zf.writestr("ncbi_dataset/data/assembly_data_report.jsonl", report + "\n")


@pytest.fixture
def stub_datasets(monkeypatch):
    """Replace `datasets` with a local writer, recording each invocation."""
    calls = []

    def fake_run(include, accs_file, output_path, annotated, **kwargs):
        from pathlib import Path

        with open(accs_file, "r") as raw:
            accs = [x for x in raw.read().split("\n") if x]
        calls.append(accs)
        _fake_report(Path(output_path), accs)
        return 0

    monkeypatch.setattr(mod, "run_datasets", fake_run)
    return calls


def test_chunk_progress_is_reported(tmp_path, stub_datasets, monkeypatch):
    """One bar over the whole acquisition, advanced per chunk."""
    bars = []

    class FakeTqdm:
        def __init__(self, iterable, **kwargs):
            self.iterable = iterable
            bars.append(kwargs)

        def __iter__(self):
            return iter(self.iterable)

    monkeypatch.setattr(mod, "tqdm", FakeTqdm)
    mod.dwnld_data_reports(
        [f"GCA_{i:09d}.1" for i in range(30)], str(tmp_path) + "/", chunk=10
    )

    assert len(bars) == 1
    assert bars[0]["total"] == 3
    assert bars[0]["disable"] is False


def test_single_chunk_draws_no_bar(tmp_path, stub_datasets, monkeypatch):
    """A bar that fills in one step is noise."""
    bars = []
    monkeypatch.setattr(
        mod, "tqdm", lambda iterable, **kw: bars.append(kw) or iterable
    )
    mod.dwnld_data_reports(["GCA_000000000.1"], str(tmp_path) + "/", chunk=100)
    assert bars[0]["disable"] is True


def test_reports_are_chunked(tmp_path, stub_datasets):
    accs = [f"GCA_{i:09d}.1" for i in range(250)]
    acc2org, acc2meta, failed = mod.dwnld_data_reports(
        accs, str(tmp_path) + "/", chunk=100
    )

    assert [len(c) for c in stub_datasets] == [100, 100, 50]
    assert len(acc2org) == 250
    assert set(acc2org) == set(accs)


def test_chunking_does_not_change_the_result(tmp_path, stub_datasets):
    """Whatever the chunk size, the merged report is the same."""
    accs = [f"GCA_{i:09d}.1" for i in range(60)]
    whole, _, _ = mod.dwnld_data_reports(accs, str(tmp_path / "a") + "/", chunk=1000)
    split, _, _ = mod.dwnld_data_reports(accs, str(tmp_path / "b") + "/", chunk=7)
    assert whole == split


def test_completed_chunks_are_not_redownloaded(tmp_path, stub_datasets):
    """A resumed run picks up where it stopped rather than starting over."""
    accs = [f"GCA_{i:09d}.1" for i in range(30)]
    out = str(tmp_path) + "/"
    mod.dwnld_data_reports(accs, out, chunk=10)
    assert len(stub_datasets) == 3

    stub_datasets.clear()
    acc2org, _, _ = mod.dwnld_data_reports(accs, out, chunk=10)
    assert stub_datasets == []
    assert len(acc2org) == 30


def test_shifted_chunk_boundaries_invalidate_the_cache(tmp_path, stub_datasets):
    """A cached chunk is reused only if it covers exactly the same accessions."""
    out = str(tmp_path) + "/"
    mod.dwnld_data_reports([f"GCA_{i:09d}.1" for i in range(20)], out, chunk=10)
    stub_datasets.clear()

    grown = [f"GCA_{i:09d}.1" for i in range(20, 45)]
    acc2org, _, _ = mod.dwnld_data_reports(grown, out, chunk=10)
    assert stub_datasets  # re-acquired rather than trusting stale chunk dirs
    assert set(acc2org) == set(grown)


def test_retries_stay_off_the_terminal(tmp_path, monkeypatch, caplog):
    """Retry chatter would overdraw the bar; one summary replaces it."""

    def fake_run(include, accs_file, output_path, annotated, **kwargs):
        return 1  # never writes an archive, so every attempt is spent

    monkeypatch.setattr(mod, "run_datasets", fake_run)
    accs = [f"GCA_{i:09d}.1" for i in range(20)]
    with caplog.at_level("INFO"):
        mod.dwnld_data_reports(accs, str(tmp_path) + "/", chunk=10)

    assert "Reattempting" not in caplog.text
    assert "datasets failed - " not in caplog.text
    # the failure itself is still reported, once, after the bar
    assert caplog.text.count("returned no genomes") == 1
    assert "2/2" in caplog.text


def test_a_failed_chunk_does_not_discard_the_others(tmp_path, monkeypatch):
    """One empty chunk costs its own accessions, not the whole download."""

    def fake_run(include, accs_file, output_path, annotated, **kwargs):
        from pathlib import Path

        with open(accs_file, "r") as raw:
            accs = [x for x in raw.read().split("\n") if x]
        if accs and accs[0].endswith("0000010.1"):  # the second chunk
            return 1
        _fake_report(Path(output_path), accs)
        return 0

    monkeypatch.setattr(mod, "run_datasets", fake_run)
    accs = [f"GCA_{i:09d}.1" for i in range(30)]
    acc2org, _, _ = mod.dwnld_data_reports(accs, str(tmp_path) + "/", chunk=10)

    assert len(acc2org) == 20
    assert "GCA_000000015.1" not in acc2org


def test_no_accessions_makes_no_requests(tmp_path, stub_datasets):
    acc2org, acc2meta, failed = mod.dwnld_data_reports([], str(tmp_path) + "/")
    assert stub_datasets == []
    assert (acc2org, acc2meta, failed) == ({}, {}, [])
