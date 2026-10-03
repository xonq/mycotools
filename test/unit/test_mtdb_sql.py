#! /usr/bin/env python3
"""Tests for the SQLite MTDB backend.

The governing property is *equivalence*: a database read from SQLite must be
indistinguishable from the same database read from its `.mtdb` interchange file.
Most tests below assert that directly, against the committed 42-entry
Ustilaginomycotina database, so the backend cannot drift from the flat format it
has to keep exporting.
"""
import json
import os
import sqlite3
import subprocess
import sys
from pathlib import Path

import pytest

from mycotools.lib import mtdb_sql
from mycotools.lib.dbtools import db_stem, load_omes, mtdb, omes_from_accessions

REPO_ROOT = Path(__file__).resolve().parents[2]


@pytest.fixture
def mtdb_root(tmp_path, monkeypatch, ust_mtdb):
    """An MTDB directory hierarchy holding the committed test database."""
    root = tmp_path / "mycotoolsdb"
    for sub in ("mtdb", "config", "log", "data/fna", "data/faa", "data/gff3"):
        (root / sub).mkdir(parents=True)
    flat = root / "mtdb" / "20240101.mtdb"
    flat.write_text(Path(ust_mtdb).read_text())
    (root / "config" / "mtdb.json").write_text(json.dumps({"branch": "fungi"}))

    monkeypatch.setenv("MYCODB", str(root / "mtdb") + "/")
    monkeypatch.setenv("MYCOFNA", str(root / "data/fna") + "/")
    monkeypatch.setenv("MYCOFAA", str(root / "data/faa") + "/")
    monkeypatch.setenv("MYCOGFF3", str(root / "data/gff3") + "/")
    return root


@pytest.fixture
def flat_db(mtdb_root):
    return str(mtdb_root / "mtdb" / "20240101.mtdb")


@pytest.fixture
def sql_db(mtdb_root, flat_db):
    """The same database, converted to the SQLite backend."""
    target = str(mtdb_root / "mtdb" / mtdb_sql.PRIMARY_DB_NAME)
    mtdb(flat_db).to_sql(target)
    return target


# --------------------------------------------------------------------------- #
# format detection
# --------------------------------------------------------------------------- #
def test_is_sqlite_discriminates(flat_db, sql_db):
    assert mtdb_sql.is_sqlite(sql_db)
    assert not mtdb_sql.is_sqlite(flat_db)


def test_is_sqlite_tolerates_missing_and_empty(tmp_path):
    assert not mtdb_sql.is_sqlite(None)
    assert not mtdb_sql.is_sqlite(tmp_path / "nope")
    empty = tmp_path / "empty"
    empty.write_text("")
    assert not mtdb_sql.is_sqlite(empty)


def test_schema_version_recorded(sql_db):
    conn = sqlite3.connect(sql_db)
    assert conn.execute("PRAGMA user_version").fetchone()[0] == mtdb_sql.SCHEMA_VERSION
    conn.close()


def test_future_schema_is_refused(sql_db):
    conn = sqlite3.connect(sql_db)
    conn.execute(f"PRAGMA user_version = {mtdb_sql.SCHEMA_VERSION + 1}")
    conn.commit()
    conn.close()
    with pytest.raises(mtdb_sql.MTDBSchemaError):
        mtdb_sql.connect(sql_db)


# --------------------------------------------------------------------------- #
# equivalence with the flat interchange format
# --------------------------------------------------------------------------- #
def test_sqlite_read_equals_flat_read(flat_db, sql_db):
    flat, sql = mtdb(flat_db), mtdb(sql_db)
    assert len(flat["ome"]) == len(sql["ome"])
    flat_rows = flat.set_index("ome")
    sql_rows = sql.set_index("ome")
    assert set(flat_rows) == set(sql_rows)
    for ome in flat_rows:
        assert sql_rows[ome] == flat_rows[ome], ome


def test_round_trip_is_byte_identical(tmp_path, flat_db, sql_db):
    """flat -> SQLite -> flat must reproduce the flat export exactly."""
    a, b = tmp_path / "a.mtdb", tmp_path / "b.mtdb"
    mtdb(flat_db).df2db(str(a))
    mtdb(sql_db).df2db(str(b))
    assert a.read_text() == b.read_text()


def test_taxonomy_is_normalized_not_duplicated(sql_db):
    conn = sqlite3.connect(sql_db)
    genomes = conn.execute("SELECT COUNT(*) FROM genome").fetchone()[0]
    lineages = conn.execute("SELECT COUNT(*) FROM taxonomy").fetchone()[0]
    genera = conn.execute("SELECT COUNT(DISTINCT genus) FROM genome").fetchone()[0]
    conn.close()
    assert lineages == genera
    assert lineages < genomes  # the point of normalizing


def test_genome_table_has_no_taxonomy_column(sql_db):
    conn = sqlite3.connect(sql_db)
    cols = {r[1] for r in conn.execute("PRAGMA table_info(genome)")}
    conn.close()
    assert "taxonomy" not in cols
    assert "ome" in cols and "assembly_acc" in cols


def test_paths_stored_abbreviated_survive_a_moved_root(mtdb_root, sql_db, monkeypatch):
    """An abbreviated path follows $MYCO* rather than pinning the old root."""
    moved = mtdb_root.parent / "elsewhere"
    monkeypatch.setenv("MYCOFNA", str(moved / "fna") + "/")
    monkeypatch.setenv("MYCOFAA", str(moved / "faa") + "/")
    monkeypatch.setenv("MYCOGFF3", str(moved / "gff3") + "/")
    row = mtdb(sql_db).set_index("ome")["ustmay1"]
    assert row["fna"] == str(moved / "fna" / "ustmay1.fna")


def test_explicit_absolute_paths_are_preserved(tmp_path, mtdb_root, flat_db):
    """A standalone genome's explicit path must not be rewritten to $MYCOFNA."""
    lines = Path(flat_db).read_text().splitlines()
    fields = lines[0].split("\t")
    fields[11:14] = ["/elsewhere/x.fna", "/elsewhere/x.faa", "/elsewhere/x.gff3"]
    standalone = tmp_path / "standalone.mtdb"
    standalone.write_text("\t".join(fields) + "\n")

    sql_target = str(tmp_path / "standalone.db")
    mtdb(str(standalone)).to_sql(sql_target)
    row = mtdb(sql_target).set_index("ome")[fields[0]]
    assert row["fna"] == "/elsewhere/x.fna"
    assert row["gff3"] == "/elsewhere/x.gff3"


# --------------------------------------------------------------------------- #
# indexed access
# --------------------------------------------------------------------------- #
def test_select_omes_matches_full_read(sql_db):
    full = mtdb(sql_db).set_index("ome")
    wanted = sorted(full)[:3]
    subset = mtdb(mtdb_sql.select_omes(sql_db, wanted)).set_index("ome")
    assert sorted(subset) == wanted
    for ome in wanted:
        assert subset[ome] == full[ome]


def test_select_omes_ignores_unknown_and_empty(sql_db):
    got = mtdb_sql.select_omes(sql_db, ["ustmay1", "", "not_an_ome"])
    assert got["ome"] == ["ustmay1"]


def test_select_omes_beyond_sqlite_parameter_cap(sql_db):
    """Chunking: more omes than SQLite's per-statement host-parameter limit."""
    real = mtdb_sql.omes(sql_db)
    padded = real + [f"absent{i}" for i in range(2000)]
    assert mtdb_sql.select_omes(sql_db, padded)["ome"] == sorted(real)


def test_select_ome_prefix_finds_versioned_ome(sql_db):
    """`antflo2` must resolve to the versioned row `antflo2.1`."""
    assert "antflo2.1" in mtdb_sql.select_ome_prefix(sql_db, "antflo2")["ome"]


def test_select_column_by_assembly_acc(sql_db):
    full = mtdb(sql_db).set_index("ome")
    acc = full["ustmay1"]["assembly_acc"]
    got = mtdb_sql.select_column(sql_db, "assembly_acc", [acc])
    assert got["ome"] == ["ustmay1"]


def test_select_column_rejects_unknown_column(sql_db):
    with pytest.raises(KeyError):
        mtdb_sql.select_column(sql_db, "taxonomy", ["x"])


def test_load_omes_agrees_across_backends(flat_db, sql_db):
    wanted = {"ustmay1", "acaing1"}
    from_flat = load_omes(flat_db, wanted).set_index("ome")
    from_sql = load_omes(sql_db, wanted).set_index("ome")
    assert set(from_flat) == set(from_sql) == wanted
    for ome in wanted:
        assert from_flat[ome] == from_sql[ome]


def test_infer_rank_agrees_with_in_memory(sql_db):
    db = mtdb(sql_db).set_index("ome")
    for lineage in ("Ustilaginales", "Basidiomycota", "Ustilago"):
        assert mtdb_sql.infer_rank(sql_db, lineage) == db.infer_rank(lineage.lower())


def test_genera_for_lineages_narrows_correctly(sql_db):
    genera = mtdb_sql.genera_for_lineages(sql_db, ["Ustilaginales"])
    narrowed = mtdb(mtdb_sql.select_column(sql_db, "genus", genera)).set_index("ome")
    expected = mtdb(sql_db).set_index("ome").extract_tax(["ustilaginales"])
    assert set(narrowed) == set(expected)


def test_omes_and_count(sql_db, flat_db):
    assert mtdb_sql.count(sql_db) == len(mtdb(flat_db)["ome"])
    assert mtdb_sql.omes(sql_db) == sorted(mtdb(flat_db)["ome"])


# --------------------------------------------------------------------------- #
# writing
# --------------------------------------------------------------------------- #
def test_write_is_atomic_on_failure(tmp_path, flat_db, monkeypatch):
    """A write that dies midway must leave the previous database intact."""
    target = tmp_path / "primary.db"
    mtdb(flat_db).to_sql(str(target))
    before = target.read_bytes()

    def boom(*a, **k):
        raise RuntimeError("interrupted")

    monkeypatch.setattr(mtdb_sql, "split_taxonomy", boom)
    with pytest.raises(RuntimeError):
        mtdb(flat_db).to_sql(str(target))
    assert target.read_bytes() == before
    assert not (tmp_path / "primary.db.tmp").exists()


def test_upsert_and_delete(sql_db):
    row = mtdb(sql_db).set_index("ome")["ustmay1"]
    row = {**row, "ome": "newome1", "strain": "NEW"}
    mtdb_sql.upsert_rows(sql_db, [row])
    assert mtdb(mtdb_sql.select_omes(sql_db, ["newome1"]))["strain"] == ["NEW"]
    mtdb_sql.delete_omes(sql_db, ["newome1"])
    assert mtdb_sql.select_omes(sql_db, ["newome1"])["ome"] == []


def test_conflicting_lineages_are_reported(caplog):
    """A genus with two different lineages means a corrupt source, not a merge."""
    columns = {
        "ome": ["a1", "a2"],
        "genus": ["Ustilago", "Ustilago"],
        "taxonomy": [{"order": "Ustilaginales"}, {"order": "Somethingelse"}],
    }
    with caplog.at_level("WARNING"):
        mtdb_sql.split_taxonomy(
            columns["ome"], columns["genus"], columns["taxonomy"]
        )
    assert "conflicting lineages" in caplog.text


# --------------------------------------------------------------------------- #
# installing a new primary (the `mtdb update` write path)
# --------------------------------------------------------------------------- #
def test_write_primary_installs_sqlite_and_archives_flat(mtdb_root, flat_db):
    from mycotools.mtdb.update import write_primary

    log_dir = mtdb_root / "log" / "20240202"
    log_dir.mkdir(parents=True)
    db = mtdb(flat_db)

    new_path = write_primary(db, "20240202", str(log_dir) + "/")

    # the primary is now SQLite and holds every genome
    assert new_path == str(mtdb_root / "mtdb" / mtdb_sql.PRIMARY_DB_NAME)
    assert mtdb_sql.count(new_path) == len(db["ome"])
    # the flat primary it replaced was archived, then removed from $MYCODB
    assert (log_dir / "20240101.mtdb").is_file()
    assert not Path(flat_db).exists()
    # and a readable snapshot of the new primary sits beside the archive
    snapshot = log_dir / "20240202.mtdb"
    assert snapshot.is_file()
    assert snapshot.read_text().startswith("#ome\t")


def test_write_primary_snapshot_reloads_identically(mtdb_root, flat_db):
    from mycotools.mtdb.update import write_primary

    log_dir = mtdb_root / "log" / "20240202"
    log_dir.mkdir(parents=True)
    write_primary(mtdb(flat_db), "20240202", str(log_dir) + "/")

    primary = mtdb(str(mtdb_root / "mtdb" / mtdb_sql.PRIMARY_DB_NAME)).set_index("ome")
    snapshot = mtdb(str(log_dir / "20240202.mtdb")).set_index("ome")
    assert set(primary) == set(snapshot)
    for ome in primary:
        assert primary[ome] == snapshot[ome], ome


def test_write_primary_replaces_an_existing_sqlite_primary(mtdb_root, sql_db):
    from mycotools.mtdb.update import write_primary

    log_dir = mtdb_root / "log" / "20240303"
    log_dir.mkdir(parents=True)
    db = mtdb(sql_db).set_index("ome")
    del db["ustmay1"]
    write_primary(db.reset_index(), "20240303", str(log_dir) + "/")

    assert "ustmay1" not in set(mtdb_sql.omes(sql_db))
    assert (log_dir / mtdb_sql.PRIMARY_DB_NAME).is_file()  # prior primary archived


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize(
    "path,expected",
    [
        ("/db/mtdb/20240101.mtdb", "20240101"),
        ("/db/mtdb/mtdb.db", "mtdb"),
        ("/db/mtdb/custom", "custom"),
        (None, "mtdb"),
    ],
)
def test_db_stem(path, expected):
    assert db_stem(path) == expected


def test_omes_from_accessions():
    assert omes_from_accessions(["ustmay1_ABC", "acaing1_1", "nounderscore"]) == {
        "ustmay1",
        "acaing1",
    }


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def _cli(root, *args):
    env = dict(os.environ)
    env.update(
        MYCODB=str(root / "mtdb") + "/",
        MYCOFNA=str(root / "data/fna") + "/",
        MYCOFAA=str(root / "data/faa") + "/",
        MYCOGFF3=str(root / "data/gff3") + "/",
        PYTHONPATH=str(REPO_ROOT),
    )
    return subprocess.run(
        [sys.executable, "-m", "mycotools.mtdb", *args],
        stdin=subprocess.DEVNULL,
        capture_output=True,
        text=True,
        timeout=120,
        env=env,
    )


def test_cli_primary_prefers_sqlite(mtdb_root, flat_db, sql_db):
    assert _cli(mtdb_root).stdout.strip() == sql_db


def test_cli_primary_falls_back_to_flat(mtdb_root, flat_db):
    assert _cli(mtdb_root).stdout.strip() == flat_db


def test_cli_ome_lookup_identical_across_backends(mtdb_root, flat_db, sql_db):
    from_sql = _cli(mtdb_root, "ustmay1").stdout
    Path(sql_db).unlink()
    from_flat = _cli(mtdb_root, "ustmay1").stdout
    assert from_sql == from_flat


def test_cli_migrate(mtdb_root, flat_db):
    result = _cli(mtdb_root, "manage", "--migrate", "-y")
    assert result.returncode == 0, result.stderr
    assert mtdb_sql.is_sqlite(mtdb_root / "mtdb" / mtdb_sql.PRIMARY_DB_NAME)
    # migrating twice is a no-op, not an error
    assert _cli(mtdb_root, "manage", "--migrate", "-y").returncode == 0


def test_cli_extract_identical_across_backends(mtdb_root, flat_db, sql_db):
    from_sql = _cli(mtdb_root, "extract", "-n").stdout
    Path(sql_db).unlink()
    from_flat = _cli(mtdb_root, "extract", "-n").stdout
    assert from_sql == from_flat and from_sql


def test_cli_extract_lineage_identical_across_backends(mtdb_root, flat_db, sql_db):
    from_sql = _cli(mtdb_root, "extract", "-n", "-l", "Ustilaginales").stdout
    Path(sql_db).unlink()
    from_flat = _cli(mtdb_root, "extract", "-n", "-l", "Ustilaginales").stdout
    assert from_sql == from_flat and from_sql


def test_cli_extract_stdin_round_trip(mtdb_root, sql_db):
    dumped = _cli(mtdb_root, "extract", "-n").stdout
    env = dict(os.environ)
    env.update(
        MYCODB=str(mtdb_root / "mtdb") + "/",
        MYCOFNA=str(mtdb_root / "data/fna") + "/",
        MYCOFAA=str(mtdb_root / "data/faa") + "/",
        MYCOGFF3=str(mtdb_root / "data/gff3") + "/",
        PYTHONPATH=str(REPO_ROOT),
    )
    piped = subprocess.run(
        [sys.executable, "-m", "mycotools.mtdb", "extract", "-n", "-d", "-"],
        input=dumped,
        capture_output=True,
        text=True,
        timeout=120,
        env=env,
    )
    assert piped.returncode == 0, piped.stderr
    assert piped.stdout == dumped


def test_cli_extract_seed_is_reproducible(mtdb_root, sql_db):
    a = _cli(mtdb_root, "extract", "-n", "-r", "genus", "-a", "1", "--seed", "7").stdout
    b = _cli(mtdb_root, "extract", "-n", "-r", "genus", "-a", "1", "--seed", "7").stdout
    assert a == b and a
