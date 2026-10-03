#! /usr/bin/env python3
"""Offline tests for acquiring the MycoCosm table without caching an outage.

Covers the failure where an update died parsing the table JGI had supposedly
supplied:

    pandas.errors.ParserError: Error tokenizing data. C error:
    Expected 1 fields in line 5, saw 2

The file was not the table: JGI was down, and answers a failed request with a
404 HTML page whose fifth line is

    <meta name="viewport" content="width=device-width, initial-scale=1">

-- one comma, hence two fields where the parser had inferred one. curl reports
that page as a successful download, and the presence of the output file is what
marks the download done, so the error page was saved and every later run reread
it. These tests pin that an error page is never accepted as the table, that one
already on disk is discarded rather than parsed, and that a real table is still
read whatever its encoding.
"""
import subprocess

import pytest

import mycotools.mtdb.update as mod


NOT_FOUND_PAGE = (
    "<!DOCTYPE html>\n"
    '<!-- saved from url=(0031)https://jgi.doe.gov/nopage.html -->\n'
    '<html lang="en-US"><head><meta http-equiv="Content-Type" '
    'content="text/html; charset=UTF-8">\n'
    "\n"
    '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
    "<title>Page not found - DOE Joint Genome Institute</title>\n"
)

HEADER = (
    "##,name,portal,NCBI Taxon,assembly length,#of genes,is restricted,"
    "is public,is published,is superseded,superseded by,publication(s),"
    "pubmed id(s),doi id(s)\n"
)
ROW = "1,Amoeboaphelidium occidentale KS120 v1.0,Amoocc1,1498855,21443760,7883,N,Y,N,N,,,,\n"


@pytest.fixture(autouse=True)
def no_backoff(monkeypatch):
    """The retry pause is real time; nothing here is testing the clock."""
    monkeypatch.setattr(mod.time, "sleep", lambda *a, **kw: None)


def _table(rows=1, encoding="utf-8", publication="Grigoriev et al."):
    # the trailing ,,,, is superseded by/publication(s)/pubmed id(s)/doi id(s)
    body = "".join(
        ROW.replace("1,", f"{i + 1},", 1).replace(",,,,", f",,{publication},,", 1)
        for i in range(rows)
    )
    return (HEADER + body).encode(encoding)


def _stub_curl(monkeypatch, behavior):
    """Route the curl subprocess to `behavior(out_path)`, recording each call."""
    calls = []

    def fake_run(cmd, **kwargs):
        out_path = cmd[cmd.index("-o") + 1]
        calls.append(out_path)
        payload = behavior(len(calls))
        if payload is None:  # what --fail does to an HTTP error: no file, rc 1
            raise subprocess.CalledProcessError(22, cmd)
        with open(out_path, "wb") as out:
            out.write(payload)

    monkeypatch.setattr(mod.subprocess, "run", fake_run)
    monkeypatch.setattr(mod, "find_execs", lambda *a, **kw: ["curl"])
    return calls


# --------------------------------------------------------------------------- #
# an error page is not a table
# --------------------------------------------------------------------------- #
def test_the_404_page_is_rejected_rather_than_parsed(tmp_path):
    """The exact payload that produced the ParserError, named for what it is."""
    table = tmp_path / "mycocosm.csv"
    table.write_text(NOT_FOUND_PAGE)

    with pytest.raises(ValueError):
        mod.read_mycocosm(str(table))


def test_a_served_error_page_is_never_saved_as_the_table(tmp_path, monkeypatch):
    """A 200 is not proof the table was served; only parsing it is."""
    out_file = str(tmp_path / "mycocosm.csv")
    _stub_curl(monkeypatch, lambda n: NOT_FOUND_PAGE.encode())

    with pytest.raises(SystemExit) as exc:
        mod.dwnld_mycocosm(out_file, max_attempts=2)

    assert exc.value.code == 24
    # the run must not leave behind a file that later runs would take as done
    assert not mod.Path(out_file).is_file()
    assert not mod.Path(out_file + ".tmp").is_file()


def test_a_cached_error_page_is_discarded_and_refetched(tmp_path, monkeypatch):
    """The bug outlived the outage: the presence of the file marks the download
    done, so a saved error page fails every subsequent run identically."""
    out_file = str(tmp_path / "mycocosm.csv")
    with open(out_file, "w") as out:
        out.write(NOT_FOUND_PAGE)
    calls = _stub_curl(monkeypatch, lambda n: _table(rows=3))

    jgi_df = mod.dwnld_mycocosm(out_file)

    assert len(calls) == 1  # it refetched rather than reread
    assert list(jgi_df["portal"]) == ["Amoocc1"] * 3


def test_a_cached_table_is_not_refetched(tmp_path, monkeypatch):
    """Discarding a bad cache must not discard a good one -- this table is one
    download per run and the reason the file is kept at all."""
    out_file = str(tmp_path / "mycocosm.csv")
    with open(out_file, "wb") as out:
        out.write(_table(rows=2))
    calls = _stub_curl(monkeypatch, lambda n: pytest.fail("refetched a good table"))

    jgi_df = mod.dwnld_mycocosm(out_file)

    assert calls == []
    assert len(jgi_df) == 2


def test_an_outage_that_lifts_is_retried(tmp_path, monkeypatch):
    """JGI being down is transient; the next attempt gets the table."""
    out_file = str(tmp_path / "mycocosm.csv")
    calls = _stub_curl(monkeypatch, lambda n: None if n < 2 else _table(rows=4))

    jgi_df = mod.dwnld_mycocosm(out_file)

    assert len(calls) == 2
    assert len(jgi_df) == 4
    assert mod.Path(out_file).is_file()


def test_a_truncated_table_is_not_accepted(tmp_path, monkeypatch):
    """A dropped transfer can leave a prefix that still parses as CSV; without
    the header there is nothing to key the JGI merge on."""
    out_file = str(tmp_path / "mycocosm.csv")
    _stub_curl(monkeypatch, lambda n: ROW.encode())  # rows, no header

    with pytest.raises(SystemExit):
        mod.dwnld_mycocosm(out_file, max_attempts=1)


# --------------------------------------------------------------------------- #
# a real table is read whatever its encoding
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("encoding", ["utf-8", "cp1252", "latin1"])
def test_the_table_is_read_in_any_encoding_jgi_serves(tmp_path, encoding):
    """MycoCosm declares no encoding and its publication and strain fields
    carry non-ASCII names."""
    table = tmp_path / "mycocosm.csv"
    with open(table, "wb") as out:
        out.write(_table(publication="Grünwald et al.", encoding=encoding))

    jgi_df = mod.read_mycocosm(str(table))

    assert list(jgi_df["portal"]) == ["Amoocc1"]


def test_utf8_is_tried_before_the_single_byte_codecs(tmp_path):
    """cp1252 and latin1 map nearly every byte, so they never raise on input
    that is not theirs -- they decode it to mojibake. Trying utf-8 first is
    what makes the fallback able to choose rather than always take the first."""
    table = tmp_path / "mycocosm.csv"
    with open(table, "wb") as out:
        out.write(_table(publication="Grünwald et al.", encoding="utf-8"))

    jgi_df = mod.read_mycocosm(str(table))

    assert "Grünwald et al." in list(jgi_df["publication(s)"])
    assert "GrÃ¼nwald et al." not in list(jgi_df["publication(s)"])


def test_undecodable_bytes_fall_through_to_latin1(tmp_path):
    """One malformed field must not cost the whole table; latin1 is the last
    resort precisely because it maps every byte."""
    table = tmp_path / "mycocosm.csv"
    with open(table, "wb") as out:
        out.write(HEADER.encode() + b"1,\xff\xfe\x00\x00bad,Amoocc1\n")

    # latin1 maps every byte, so something is always returned; the point is
    # that it is the table, not an exception
    jgi_df = mod.read_mycocosm(str(table))
    assert "portal" in jgi_df.columns


def test_quoted_headers_are_stripped(tmp_path):
    """JGI has served the header row quoted; the column names are looked up
    literally downstream."""
    table = tmp_path / "mycocosm.csv"
    with open(table, "w") as out:
        out.write('"##","name","portal"\n1,Amoeboaphelidium occidentale,Amoocc1\n')

    jgi_df = mod.read_mycocosm(str(table))

    assert "portal" in jgi_df.columns
    assert "name" in jgi_df.columns
