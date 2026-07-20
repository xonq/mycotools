#! /usr/bin/env python3
"""Full-script tests for ``mycotools.mtdb.update`` user-argument scenarios.

Scope and rationale
--------------------
``mtdb.update`` orchestrates MycotoolsDB (MTDB) initialization/updating. The
heavy lifting (``ref_update``/``rogue_update``/``taxonomy_update``) depends on
authenticated JGI/NCBI access (Entrez, MycoCosm, ``datasets``) and downloads
whole genomes, so it is deliberately *not* exercised here. Instead these tests
cover the layers a user actually drives through command-line arguments, which
run before any login/network is touched:

  * ``control_flow`` argument validation (the mutually-exclusive / required-flag
    guards that ``sys.exit`` with distinct codes). ``loginCheck`` is guarded by
    ``if not ncbi_email`` and every validation guard fires *before* that line,
    so passing a dummy ``ncbi_email`` keeps these tests fully offline.
  * The argparse CLI surface (``python -m mycotools.mtdb.update``): help,
    unknown flags, and type-checked options all exit at the parser, before the
    ``datasets`` dependency check or ``control_flow``.
  * The login-free helper "modules" the script is built from, driven with the
    committed ``test/ust.mtdb`` fungal database and an in-code lineage filter
    (as permitted by the task): config generation, internal dereplication,
    genus-only lineage constraint extraction, and the ledger/sidecar parsers.

Run with::

    micromamba run -n mycotools pytest test/update_mtdb/ -v

All tests are quick (no downloads, no genome assimilation) and require no JGI or
NCBI credentials.
"""
import os
import sys
import subprocess
import tempfile
from pathlib import Path

import pandas as pd
import pytest

import mycotools.mtdb.update as u
from mycotools.lib.dbtools import db2df


# --------------------------------------------------------------------------- #
# Paths / fixtures
# --------------------------------------------------------------------------- #
REPO_ROOT = Path(__file__).resolve().parents[2]
UST_MTDB = REPO_ROOT / "test" / "ust.mtdb"


@pytest.fixture(autouse=True)
def offline_env(monkeypatch):
    """Force the login-free branches for every test.

    * Removing ``MYCODB`` makes ``primaryDB()`` return ``None`` and keeps
      ``control_flow`` out of the "already initialized" config path.
    * ``MYCOFNA``/``MYCOFAA``/``MYCOGFF3`` are only ever string-concatenated by
      the pandas MTDB machinery (``pd2mtdb``/``db2df``); dummy prefixes satisfy
      them without any filesystem or network access.
    """
    monkeypatch.delenv("MYCODB", raising=False)
    monkeypatch.setenv("MYCOFNA", "/tmp/mtdb_test_fna/")
    monkeypatch.setenv("MYCOFAA", "/tmp/mtdb_test_faa/")
    monkeypatch.setenv("MYCOGFF3", "/tmp/mtdb_test_gff3/")


@pytest.fixture
def tmp_ledger_dir():
    with tempfile.TemporaryDirectory() as d:
        yield Path(d)


# --------------------------------------------------------------------------- #
# control_flow argument validation
# --------------------------------------------------------------------------- #
# control_flow positional signature (see mtdb.update.control_flow):
#   init, update, reference, add, taxonomy, predb, save, nonpublished,
#   ncbi_only, lineage, rank, kingdom, failed, forbidden, resume, no_md5, cpu
_CF_DEFAULTS = dict(
    init=None,
    update=False,
    reference=None,
    add=None,
    taxonomy=False,
    predb=None,
    save=False,
    nonpublished=False,
    ncbi_only=False,
    lineage=None,
    rank=None,
    kingdom="fungi",
    failed=False,
    forbidden=False,
    resume=None,
    no_md5=False,
    cpu=1,
    # dummy credential -> skips loginCheck() entirely
    ncbi_email="tester@example.com",
)


def _run_control_flow(**overrides):
    """Invoke control_flow with defaults + overrides; return its SystemExit code
    (or ``None`` if it returned without exiting)."""
    kwargs = dict(_CF_DEFAULTS)
    kwargs.update(overrides)
    try:
        u.control_flow(**kwargs)
    except SystemExit as exc:
        return exc.code
    return None


# Each case is a documented "user argument scenario" -> expected exit code.
# The exit codes are the script's own contract (see control_flow).
VALIDATION_CASES = [
    # id                              overrides                                              exit
    ("no-mode",                       dict(),                                                15),
    # --predb is NOT counted as a standalone mode by the first guard, so
    # `--predb X` on its own exits 15 (no mode), not 18.
    ("predb-alone",                   dict(predb="/tmp/p.tsv"),                              15),
    ("reference-without-init",        dict(reference=str(UST_MTDB)),                         14),
    ("lineage-without-rank",          dict(init="/tmp/x", lineage="Ceraceosorus"),           16),
    ("lineage-without-init",          dict(update=True, lineage="Ceraceosorus",
                                           rank="genus"),                                     17),
    # --predb without --init, but with another mode set so we pass the first guard
    ("predb-without-init",            dict(update=True, predb="/tmp/p.tsv"),                 18),
    ("predb-with-lineage",            dict(init="/tmp/x", predb="/tmp/p.tsv",
                                           lineage="Foo", rank="genus"),                      20),
    ("reference-plus-add",            dict(init="/tmp/x", reference=str(UST_MTDB),
                                           add="/tmp/a.mtdb"),                                13),
    ("reference-plus-predb",          dict(init="/tmp/x", reference=str(UST_MTDB),
                                           predb="/tmp/p.tsv"),                               19),
    ("invalid-kingdom",               dict(update=True, kingdom="xyz"),                      431),
    ("lineage-rank-length-mismatch",  dict(init="/tmp/x", lineage="A,B", rank="genus"),      18),
    ("invalid-rank",                  dict(init="/tmp/x", lineage="Foo", rank="badrank"),    22),
]


@pytest.mark.parametrize(
    "overrides,expected",
    [pytest.param(o, e, id=i) for (i, o, e) in VALIDATION_CASES],
)
def test_control_flow_validation_exit_codes(overrides, expected):
    """Bad argument combinations exit early with their documented codes,
    before any login/network is attempted."""
    assert _run_control_flow(**overrides) == expected


@pytest.mark.parametrize(
    "kingdom", ["f", "a", "b", "p", "r", "fungi", "animals",
                "bacteria", "plants", "archaea", "FUNGI"],
)
def test_control_flow_accepts_valid_kingdoms(kingdom):
    """Valid --kingdom values (abbreviations, full names, any case) pass the
    kingdom gate; with no mode set they fall through to the mode guard (15),
    proving they were NOT rejected as invalid (431)."""
    code = _run_control_flow(kingdom=kingdom)
    assert code != 431
    assert code == 15


@pytest.mark.parametrize("kingdom", ["xyz", "fungee", "protist", ""])
def test_control_flow_rejects_invalid_kingdoms(kingdom):
    assert _run_control_flow(update=True, kingdom=kingdom) == 431


# --------------------------------------------------------------------------- #
# argparse CLI surface (subprocess -> real parser in main())
# --------------------------------------------------------------------------- #
def _cli(*args, timeout=60):
    """Run `python -m mycotools.mtdb.update <args>` offline (stdin closed)."""
    return subprocess.run(
        [sys.executable, "-m", "mycotools.mtdb.update", *args],
        stdin=subprocess.DEVNULL,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        timeout=timeout,
    )


def test_cli_help_exits_zero():
    res = _cli("--help")
    assert res.returncode == 0
    assert "usage" in res.stdout.lower()
    # a couple of real options should be advertised
    assert "--init" in res.stdout
    assert "--lineage" in res.stdout


def test_cli_unknown_flag_is_rejected():
    res = _cli("--definitely-not-a-flag")
    assert res.returncode == 2  # argparse usage error


@pytest.mark.parametrize("bad", [("--cpu", "notanint"), ("--resume", "notadate")])
def test_cli_type_checked_options_reject_bad_values(bad):
    """--cpu and --resume are ``type=int``; non-integers fail at the parser."""
    res = _cli(*bad)
    assert res.returncode == 2


# --------------------------------------------------------------------------- #
# Login-free helper "modules"
# --------------------------------------------------------------------------- #
def test_gen_config_shape_and_values():
    cfg = u.gen_config(
        branch="fungi",
        forbidden="$MYCODB/log/forbidden.tsv",
        repo="git@example:repo",
        rogue=True,
        nonpublished="yes",
        jgi=True,
        rank2lineages={"genus": ["ceraceosorus"]},
    )
    assert cfg == {
        "forbidden": "$MYCODB/log/forbidden.tsv",
        "repository": "git@example:repo",
        "branch": "fungi",
        "nonpublished": "yes",
        "rogue": True,
        "jgi": True,
        "lineage_constraints": {"genus": ["ceraceosorus"]},
    }


def test_ust_mtdb_loads():
    """Sanity: the reference test database loads into the expected shape."""
    df = db2df(str(UST_MTDB))
    assert len(df) == 42
    assert set(df["source"]) == {"jgi", "ncbi"}
    assert "assembly_acc" in df.columns


def test_internal_redundancy_check_noop_on_clean_db():
    """ust.mtdb has no redundant assembly accessions -> nothing is dropped."""
    df = db2df(str(UST_MTDB))
    out = u.internal_redundancy_check(df)
    assert isinstance(out, pd.DataFrame)
    assert set(out["ome"]) == set(df["ome"])


def test_internal_redundancy_check_keeps_highest_ncbi_version():
    """Two versions of the same NCBI assembly base collapse to the higher one."""
    df = db2df(str(UST_MTDB))
    template = df[df["source"] == "ncbi"].iloc[0].copy()

    lo = template.copy()
    lo["ome"], lo["assembly_acc"] = "dupomelo", "GCA_999999999.1"
    hi = template.copy()
    hi["ome"], hi["assembly_acc"] = "dupomehi", "GCA_999999999.2"

    df_dup = pd.concat([df, pd.DataFrame([lo, hi])], ignore_index=True)
    out = u.internal_redundancy_check(df_dup)
    kept = set(out["ome"])

    assert "dupomehi" in kept          # higher version retained
    assert "dupomelo" not in kept      # lower version dereplicated
    assert len(out) == len(df) + 1     # exactly one of the two added rows kept


def test_extract_constraint_lineages_genus_only_is_offline():
    """A genus-only constraint filters purely on the 'genus' column and returns
    before any Entrez taxonomy query (ncbi_api=None would otherwise be used for
    a network call)."""
    df = db2df(str(UST_MTDB))
    tax_dicts, filtered = u.extract_constraint_lineages(
        df.copy(), None, "fungi", {"genus": ["ceraceosorus"]}, {}, "/tmp/unused.tsv"
    )
    assert set(filtered["genus"]) == {"Ceraceosorus"}
    assert 0 < len(filtered) < len(df)
    # nothing was queried, so no taxonomy dict was accumulated
    assert tax_dicts == {}


# --- ledger / sidecar parsers (round-trips) -------------------------------- #
def test_read_ledger_skips_comments_and_blanks(tmp_ledger_dir):
    p = tmp_ledger_dir / "ledger.txt"
    p.write_text("# a comment\n\nline1\n   \nline2\n")
    assert u._read_ledger(str(p)) == ["line1", "line2"]


def test_read_ledger_missing_file_is_empty(tmp_ledger_dir):
    assert u._read_ledger(str(tmp_ledger_dir / "nope.txt")) == []


def test_forbid_omes_round_trip(tmp_ledger_dir):
    p = str(tmp_ledger_dir / "relics.txt")
    u.write_forbid_omes({"acaing1", "antflo1"}, p)
    assert u.acq_forbid_omes(p) == {"acaing1", "antflo1"}
    # a second write unions rather than overwrites
    u.write_forbid_omes({"cerbom1"}, p)
    assert u.acq_forbid_omes(p) == {"acaing1", "antflo1", "cerbom1"}


def test_parse_failed_creates_header_then_round_trips(tmp_ledger_dir):
    p = str(tmp_ledger_dir / "failed.tsv")
    assert u.parse_failed(file_path=p) == {}          # creates the header file
    u.add_failed("acaing1", "jgi", "1.0", "20260101", p)
    assert u.parse_failed(file_path=p) == {
        "acaing1": {"source": "jgi", "version": "1.0", "attempt_date": "20260101"}
    }


def test_jgi2ncbi_round_trip(tmp_ledger_dir):
    p = str(tmp_ledger_dir / "jgi2ncbi.tsv")
    assert u.parse_jgi2ncbi(p) == {}
    u.add_jgi2ncbi({"acain1": "GCA_000417875.1"}, file_path=p)
    assert u.parse_jgi2ncbi(p) == {"acain1": "GCA_000417875.1"}


def test_true_ncbi_round_trip(tmp_ledger_dir):
    p = str(tmp_ledger_dir / "supported_ncbi.tsv")
    assert u.parse_true_ncbi(p) == set()
    u.add_true_ncbi({"GCA_1", "GCA_2"}, file_path=p)
    assert u.parse_true_ncbi(p) == {"GCA_1", "GCA_2"}


def test_parse_dups(tmp_ledger_dir):
    p = str(tmp_ledger_dir / "duplicates.tsv")
    p_path = Path(p)
    p_path.write_text("# header comment\ndupcode\tA\tB\tC\n")
    assert u.parse_dups(p) == {"dupcode": ["A", "B", "C"]}
