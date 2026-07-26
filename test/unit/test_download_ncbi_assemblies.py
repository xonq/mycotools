#! /usr/bin/env python3
"""Offline tests for assembly acquisition surviving NCBI dropping a transfer.

Covers the failure where an initialization died partway through the assembly
downloads that follow MycoCosm:

    Downloading: ncbi_dataset.zip    724MB error
    Error: Download error: stream error: stream ID 5; INTERNAL_ERROR
    ERROR:                  ncbiDwnld failed 3 attempts

NCBI resets these streams, and the odds of a reset scale with the size of the
archive -- so retrying a multi-GB chunk whole retries it at exactly the size
that was already failing. These tests pin the three behaviors that keep one
reset from ending a run: the partial archive is cleared before a retry, a chunk
that never completes is halved rather than retried whole, and a chunk that is
genuinely unavailable is reported instead of discarding the chunks that worked.
"""
import json
import subprocess
import zipfile

import pytest

from mycotools.download import ncbi


REQ_FILES = {"fna"}


@pytest.fixture(autouse=True)
def no_backoff(monkeypatch):
    """The retry pause is real time; nothing here is testing the clock."""
    monkeypatch.setattr(ncbi.time, "sleep", lambda *a, **kw: None)


def _fake_download(out_dir, accs):
    """Write the ncbi_dataset.zip `datasets` leaves on success, so the caller's
    own unzip/parse path is the one under test."""
    catalog = {
        "assemblies": [
            {
                "accession": acc,
                "files": [
                    {
                        "fileType": "GENOMIC_NUCLEOTIDE_FASTA",
                        "filePath": f"{acc}/{acc}.fna",
                    }
                ],
            }
            for acc in accs
        ]
    }
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
    with zipfile.ZipFile(out_dir + "ncbi_dataset.zip", "w") as zf:
        zf.writestr("ncbi_dataset/data/dataset_catalog.json", json.dumps(catalog))
        zf.writestr("ncbi_dataset/data/assembly_data_report.jsonl", report + "\n")
        for acc in accs:
            zf.writestr(f"ncbi_dataset/data/{acc}/{acc}.fna", ">ctg\nACGT\n")


def _truncated_download(out_dir):
    """What a dropped stream leaves behind: a prefix of an archive."""
    with open(out_dir + "ncbi_dataset.zip", "wb") as out:
        out.write(b"PK\x03\x04" + b"\x00" * 512)


def _stub(monkeypatch, behavior):
    """Route run_datasets to `behavior(accs, out_dir)`, recording each call."""
    calls = []

    def fake_run(include, accs_file, output_path, annotated, **kwargs):
        with open(accs_file, "r") as raw:
            accs = [x for x in raw.read().split("\n") if x]
        calls.append(accs)
        return behavior(accs, output_path)

    monkeypatch.setattr(ncbi, "run_datasets", fake_run)
    return calls


def _download(tmp_path, accs, **kwargs):
    out_dir = str(tmp_path) + "/"
    return ncbi.download_datasets(
        accs,
        out_dir + "assembly_accs.txt",
        "genome",
        REQ_FILES,
        False,
        out_dir,
        **kwargs,
    )


# --------------------------------------------------------------------------- #
# a dropped transfer is retried, not fatal
# --------------------------------------------------------------------------- #
def test_dropped_transfer_is_retried(tmp_path, monkeypatch):
    """The reset that killed the run is transient; the second attempt gets it."""
    accs = [f"GCA_{i:09d}.1" for i in range(5)]

    def behavior(chunk, out_dir):
        if len(calls) == 1:  # first call
            _truncated_download(out_dir)
            return 1
        _fake_download(out_dir, chunk)
        return 0

    calls = _stub(monkeypatch, behavior)
    acc2files, acc2org, failed = _download(tmp_path, accs)

    assert sorted(acc2files) == accs
    assert len(calls) == 2


def test_partial_archive_never_reaches_the_next_attempt(tmp_path, monkeypatch):
    """datasets leaves the truncated file behind; parsing it fails the retry
    before it starts, and the run reported three attempts against one archive."""
    accs = [f"GCA_{i:09d}.1" for i in range(5)]
    saw_stale = []

    def behavior(chunk, out_dir):
        saw_stale.append(ncbi.Path(out_dir + "ncbi_dataset.zip").is_file())
        if len(calls) < 3:
            _truncated_download(out_dir)
            return 1
        _fake_download(out_dir, chunk)
        return 0

    calls = _stub(monkeypatch, behavior)
    acc2files, _, _ = _download(tmp_path, accs)

    assert sorted(acc2files) == accs
    assert saw_stale == [False, False, False]


def test_a_call_that_writes_nothing_is_not_a_crash(tmp_path, monkeypatch):
    """datasets can die before it opens the archive -- an absent file is the
    same failure as a truncated one, not a FileNotFoundError out of the parser."""
    accs = [f"GCA_{i:09d}.1" for i in range(5)]

    def behavior(chunk, out_dir):
        if len(calls) == 1:
            return 1  # no archive at all
        _fake_download(out_dir, chunk)
        return 0

    calls = _stub(monkeypatch, behavior)
    acc2files, _, _ = _download(tmp_path, accs)

    assert sorted(acc2files) == accs


# --------------------------------------------------------------------------- #
# a chunk that never completes is halved
# --------------------------------------------------------------------------- #
def test_oversized_chunk_is_split_rather_than_retried_whole(tmp_path, monkeypatch):
    """The reset is a function of transfer size, so the recovery has to shrink
    the transfer; retrying 100 accessions three times only repeats the size that
    was failing."""
    accs = [f"GCA_{i:09d}.1" for i in range(100)]

    def behavior(chunk, out_dir):
        if len(chunk) > 25:  # what NCBI would not hold open long enough to send
            _truncated_download(out_dir)
            return 1
        _fake_download(out_dir, chunk)
        return 0

    calls = _stub(monkeypatch, behavior)
    acc2files, acc2org, failed = _download(tmp_path, accs)

    assert sorted(acc2files) == accs  # nothing lost to the resets
    assert sorted(acc2org) == accs
    assert failed == []
    assert max(len(c) for c in calls if len(c) <= 25) <= 25


def test_a_half_that_works_survives_the_half_that_does_not(tmp_path, monkeypatch):
    """Losing a chunk wholesale is what made one bad accession expensive."""
    accs = [f"GCA_{i:09d}.1" for i in range(40)]
    poison = accs[0]

    def behavior(chunk, out_dir):
        if poison in chunk:
            _truncated_download(out_dir)
            return 1
        _fake_download(out_dir, chunk)
        return 0

    _stub(monkeypatch, behavior)
    acc2files, _, _ = _download(tmp_path, accs, min_accs=1)

    assert acc2files is not False
    assert poison not in acc2files
    assert len(acc2files) >= len(accs) // 2


def test_splitting_stops_rather_than_hammering_a_dead_repository(
    tmp_path, monkeypatch
):
    """When nothing is being served, splitting only multiplies the calls."""
    accs = [f"GCA_{i:09d}.1" for i in range(100)]

    def behavior(chunk, out_dir):
        _truncated_download(out_dir)
        return 1

    calls = _stub(monkeypatch, behavior)
    acc2files, acc2org, failed = _download(tmp_path, accs, min_accs=10)

    assert (acc2files, acc2org, failed) == (False, False, False)
    # a repository that serves nothing at any size is down rather than
    # overloaded: it gives up one halving past the depth of the split tree,
    # rather than walking the whole tree and spending ~90 calls on the same answer
    # 6 is the floor here: at min_accs it is abandoned rather than split, so the
    # last two batches are the halves of 12
    assert [len(c) for c in calls] == sum(
        ([n] * 3 for n in (100, 50, 25, 12, 6, 6)), []
    )


# --------------------------------------------------------------------------- #
# one dead chunk does not discard the chunks that succeeded
# --------------------------------------------------------------------------- #
def _ncbi_df(accs):
    import pandas as pd

    return pd.DataFrame({"assembly_acc": accs, "version": ["20240101"] * len(accs)})


def test_one_dead_chunk_leaves_the_rest_of_the_run_intact(tmp_path, monkeypatch):
    """The run died on chunk 3 and threw away the ~4GB chunks 1 and 2 had
    already retrieved; a chunk NCBI will not serve is a reported failure."""
    accs = [f"GCA_{i:09d}.1" for i in range(30)]
    dead = set(accs[10:20])
    # the reattempt against the other repository has to fail too, or the chunk
    # is not actually lost and the test proves nothing
    dead_either_repo = dead | {a.replace("GCA_", "GCF_") for a in dead}

    def fake_download(chunk, *args, **kwargs):
        if set(chunk) <= dead_either_repo:
            return False, False, False
        return (
            {a: {"fna": f"/x/{a}.fna"} for a in chunk},
            {a: {"genus": "Aspergillus", "species": "sp.", "strain": ""} for a in chunk},
            [],
        )

    monkeypatch.setattr(ncbi, "download_datasets", fake_download)
    new_df, rep_failed = ncbi.main(
        ncbi_df=_ncbi_df(accs),
        gff3=False,
        output_path=str(tmp_path) + "/",
        chunk=10,
    )

    assert sorted(new_df["assembly_acc"]) == sorted(set(accs) - dead)
    assert {x[0] for x in rep_failed} == dead


def test_consecutive_dead_chunks_still_stop_the_run(tmp_path, monkeypatch):
    """Continuing past one chunk is recovery; continuing past every chunk would
    silently mark a whole database as permanently failed."""
    accs = [f"GCA_{i:09d}.1" for i in range(50)]

    monkeypatch.setattr(
        ncbi, "download_datasets", lambda *a, **kw: (False, False, False)
    )
    with pytest.raises(SystemExit) as exc:
        ncbi.main(
            ncbi_df=_ncbi_df(accs),
            gff3=False,
            output_path=str(tmp_path) + "/",
            chunk=10,
        )
    assert exc.value.code == 10


# --------------------------------------------------------------------------- #
# nothing leaves the run unaccounted for
# --------------------------------------------------------------------------- #
def _served(accs):
    """What download_datasets returns for accessions NCBI supplied in full."""
    return (
        {a: {"fna": f"/x/{a}.fna"} for a in accs},
        {a: {"genus": "Aspergillus", "species": "sp.", "strain": ""} for a in accs},
        [],
    )


def test_an_unserved_accession_is_reported_rather_than_dropped(
    tmp_path, monkeypatch, caplog
):
    """An accession absent from the archive is absent from the catalog too, so
    it is never named in the parser's failed list -- only subtracting what came
    back finds it. Rebuilding the failure list from the reattempt alone dropped
    214 of 525 genomes out of an initialization with no entry in failed.tsv,
    which left them invisible to --rerun as well."""
    accs = [f"GCA_{i:09d}.1" for i in range(6)]
    served = set(accs[:4])

    monkeypatch.setattr(
        ncbi, "download_datasets",
        lambda chunk, *a, **kw: _served([x for x in chunk if x in served]),
    )
    monkeypatch.setattr(ncbi, "resolve_paired_accs", lambda *a, **kw: {})

    with caplog.at_level("DEBUG", logger=ncbi.logger.name):
        new_df, rep_failed = ncbi.main(
            ncbi_df=_ncbi_df(accs),
            gff3=False,
            output_path=str(tmp_path) + "/",
            chunk=10,
        )

    assert sorted(new_df["assembly_acc"]) == sorted(served)
    # every requested accession is either retrieved or reported, never neither
    assert {x[0] for x in rep_failed} == set(accs) - served
    # and the failure list is what reported them: the catch-all exists to make
    # a miss impossible, so it firing here would mean the accounting is back to
    # relying on it. Pinning that keeps the two fixes independently covered
    assert "neither retrieved nor reported" not in caplog.text


def test_a_lowercase_accession_is_retrieved_rather_than_lost(tmp_path, monkeypatch):
    """The accession column is uppercased because NCBI's is, and the index is
    taken from that column. Normalizing after the index was derived left the two
    disagreeing, and a row whose index matched neither the downloads nor the
    failure list dropped out of the run entirely."""
    accs = ["gca_000000001.1", "GCA_000000002.1"]

    monkeypatch.setattr(
        ncbi, "download_datasets", lambda chunk, *a, **kw: _served(chunk)
    )
    monkeypatch.setattr(ncbi, "resolve_paired_accs", lambda *a, **kw: {})

    new_df, rep_failed = ncbi.main(
        ncbi_df=_ncbi_df(accs),
        gff3=False,
        output_path=str(tmp_path) + "/",
        chunk=10,
    )

    assert sorted(new_df["assembly_acc"]) == ["GCA_000000001.1", "GCA_000000002.1"]
    assert rep_failed == []


def test_requests_carry_the_normalized_accession(tmp_path, monkeypatch):
    """Uppercasing the column but requesting the raw value would ask NCBI for an
    accession it does not recognize, so the normalization has to reach the
    request rather than only the bookkeeping."""
    requested = []

    def fake_download(chunk, *a, **kw):
        requested.append(list(chunk))
        return _served(chunk)

    monkeypatch.setattr(ncbi, "download_datasets", fake_download)
    monkeypatch.setattr(ncbi, "resolve_paired_accs", lambda *a, **kw: {})

    ncbi.main(
        ncbi_df=_ncbi_df(["gca_000000001.1"]),
        gff3=False,
        output_path=str(tmp_path) + "/",
        chunk=10,
    )

    assert requested == [["GCA_000000001.1"]]


def test_a_frame_indexed_on_another_column_reports_rather_than_drops(
    tmp_path, monkeypatch
):
    """The catch-all only means something if something can reach it. Indexing on
    a column other than the accession is a misconfiguration rather than a
    supported mode -- but it must surface as every row being reported, not as
    every row disappearing."""
    accs = ["GCA_000000001.1", "GCA_000000002.1"]
    ncbi_df = _ncbi_df(accs)
    ncbi_df["biosample"] = ["SAMN00000001", "SAMN00000002"]

    monkeypatch.setattr(
        ncbi, "download_datasets", lambda chunk, *a, **kw: _served(chunk)
    )
    monkeypatch.setattr(ncbi, "resolve_paired_accs", lambda *a, **kw: {})

    new_df, rep_failed = ncbi.main(
        ncbi_df=ncbi_df,
        gff3=False,
        output_path=str(tmp_path) + "/",
        chunk=10,
        column="biosample",
    )

    assert {x[0] for x in rep_failed} == set(ncbi_df["biosample"])
    assert len(new_df) == 0


def test_the_counterpart_is_the_one_ncbi_holds(tmp_path, monkeypatch):
    """GenBank and RefSeq version independently: GCA_017499595.2 pairs with
    GCF_017499595.1, so swapping the prefix and keeping the version asks for a
    GCF_017499595.2 that does not exist and the reattempt cannot ever land."""
    accs = ["GCA_017499595.2"]
    requested = []

    def fake_download(chunk, *a, **kw):
        requested.append(list(chunk))
        if list(chunk) == ["GCF_017499595.1"]:
            return _served(chunk)
        return {}, {}, []

    monkeypatch.setattr(ncbi, "download_datasets", fake_download)
    monkeypatch.setattr(
        ncbi, "resolve_paired_accs",
        lambda *a, **kw: {"GCA_017499595.2": "GCF_017499595.1"},
    )

    new_df, rep_failed = ncbi.main(
        ncbi_df=_ncbi_df(accs),
        gff3=False,
        output_path=str(tmp_path) + "/",
        chunk=10,
    )

    assert requested[1] == ["GCF_017499595.1"]
    assert list(new_df["assembly_acc"]) == accs
    assert rep_failed == []


def test_an_unmirrored_accession_is_not_reattempted(tmp_path, monkeypatch):
    """NCBI reporting no counterpart is the answer; requesting one regardless
    spends a call to be told the same thing."""
    accs = [f"GCA_{i:09d}.1" for i in range(3)]
    requested = []

    def fake_download(chunk, *a, **kw):
        requested.append(list(chunk))
        return {}, {}, []

    monkeypatch.setattr(ncbi, "download_datasets", fake_download)
    monkeypatch.setattr(ncbi, "resolve_paired_accs", lambda *a, **kw: {})

    new_df, rep_failed = ncbi.main(
        ncbi_df=_ncbi_df(accs),
        gff3=False,
        output_path=str(tmp_path) + "/",
        chunk=10,
    )

    assert len(requested) == 1  # the primary pass only
    assert {x[0] for x in rep_failed} == set(accs)


# --------------------------------------------------------------------------- #
# resolving the counterpart
# --------------------------------------------------------------------------- #
def _stub_summary(monkeypatch, stdout="", returncode=0, raises=None):
    calls = []

    def fake_run(cmd, **kwargs):
        calls.append(cmd)
        if raises is not None:
            raise raises
        return subprocess.CompletedProcess(cmd, returncode, stdout=stdout, stderr="")

    monkeypatch.setattr(ncbi.subprocess, "run", fake_run)
    return calls


def test_only_accessions_with_a_counterpart_are_returned(tmp_path, monkeypatch):
    """NCBI omits paired_accession for an assembly it never mirrored."""
    lines = "\n".join(
        json.dumps(x)
        for x in (
            {"accession": "GCA_017499595.2", "paired_accession": "GCF_017499595.1"},
            {"accession": "GCA_013368295.1"},
            {"accession": "GCA_002938375.1", "paired_accession": None},
        )
    )
    _stub_summary(monkeypatch, stdout=lines)

    pairs = ncbi.resolve_paired_accs(
        ["GCA_017499595.2", "GCA_013368295.1", "GCA_002938375.1"],
        str(tmp_path / "accs.txt"),
        str(tmp_path) + "/",
    )

    assert pairs == {"GCA_017499595.2": "GCF_017499595.1"}


def test_pair_lookup_is_chunked(tmp_path, monkeypatch):
    """These are metadata records rather than genomes, so the chunk sits well
    above the genome chunk for the same reason the data reports do."""
    calls = _stub_summary(monkeypatch, stdout="")

    ncbi.resolve_paired_accs(
        [f"GCA_{i:09d}.1" for i in range(1200)],
        str(tmp_path / "accs.txt"),
        str(tmp_path) + "/",
        summary_chunk=500,
    )

    assert len(calls) == 3


def test_an_unresolvable_chunk_leaves_its_accessions_failed(tmp_path, monkeypatch):
    """Losing the lookup must not lose the accessions: without a counterpart
    they keep the failure they already had, which is the outcome regardless."""
    _stub_summary(monkeypatch, returncode=1)

    pairs = ncbi.resolve_paired_accs(
        ["GCA_017499595.2"], str(tmp_path / "accs.txt"), str(tmp_path) + "/"
    )

    assert pairs == {}


def test_a_missing_datasets_does_not_end_the_run(tmp_path, monkeypatch):
    """The reattempt is an extra chance, not a step the run depends on; an
    uninstalled datasets would otherwise raise past everything downloaded."""
    _stub_summary(monkeypatch, raises=FileNotFoundError(2, "No such file", "datasets"))

    pairs = ncbi.resolve_paired_accs(
        ["GCA_017499595.2"], str(tmp_path / "accs.txt"), str(tmp_path) + "/"
    )

    assert pairs == {}


def test_the_lookup_returns_to_the_original_directory(tmp_path, monkeypatch):
    """resolve_paired_accs chdirs to run datasets; a failure that skipped the
    chdir back would silently relocate the rest of the run."""
    _stub_summary(monkeypatch, raises=FileNotFoundError(2, "No such file", "datasets"))
    before = str(ncbi.Path.cwd())

    ncbi.resolve_paired_accs(
        ["GCA_017499595.2"], str(tmp_path / "accs.txt"), str(tmp_path) + "/"
    )

    assert str(ncbi.Path.cwd()) == before
