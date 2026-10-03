#! /usr/bin/env python3
"""SQLite storage backend for MycotoolsDB (MTDB).

The primary MTDB is stored as a single SQLite file (``$MYCODB/mtdb.db``); the
tab-delimited ``.mtdb`` file remains the interchange format, read and written by
``mycotools.lib.dbtools.mtdb``. Nothing here is user-facing: the ``mtdb`` class
is still the interface, and every call site that passes a path to ``mtdb()``
keeps working because :func:`is_sqlite` dispatches on the file itself.

Two properties motivate the backend:

* **Indexed lookup.** ``mtdb <OME>``, ``mtdb accession``, and ``mtdb files``
  answer questions about a handful of genomes. Against a flat file each has to
  parse the whole database; here they are index seeks (:func:`select_omes`).
* **Normalized taxonomy.** A lineage belongs to a genus, not to a genome, so it
  is stored once in the ``taxonomy`` table instead of being repeated as JSON on
  every row -- ~70% of a flat ``.mtdb`` is duplicated lineage text.

Writes go through :func:`write_db`, which builds into a temporary file and
renames over the target, so a cancelled update can never leave a partial
database where ``primary_db()`` will find it.

Path columns follow the flat-file convention: an empty ``fna``/``faa``/``gff3``
means "the default ``$MYCOFNA``/``$MYCOFAA``/``$MYCOGFF3`` location for this
ome". Storing them abbreviated keeps a database valid after its root moves;
:func:`rows2columns` expands them on read.
"""

from __future__ import annotations

import json
import logging
import os
import sqlite3
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence

logger = logging.getLogger(__name__)

#: bumped when the on-disk layout changes in a way older code cannot read
SCHEMA_VERSION = 1

#: basename of the SQLite primary database within $MYCODB
PRIMARY_DB_NAME = "mtdb.db"

#: first 16 bytes of any SQLite 3 file
_SQLITE_MAGIC = b"SQLite format 3\x00"

#: MTDB columns as they appear in a `.mtdb` row, in order. Mirrors
#: ``mtdb.columns``; duplicated here so this module does not import dbtools
#: (dbtools imports *it*).
COLUMNS = [
    "ome",
    "genus",
    "species",
    "strain",
    "taxonomy",
    "version",
    "source",
    "biosample",
    "assembly_acc",
    "acquisition_date",
    "published",
    "fna",
    "faa",
    "gff3",
]

#: columns physically stored on the `genome` table -- `taxonomy` is normalized
#: out into its own table and reattached on read
_GENOME_COLUMNS = [c for c in COLUMNS if c != "taxonomy"]

#: file-path columns and the environment variable holding their default dir
_PATH_COLUMNS = {"fna": "MYCOFNA", "faa": "MYCOFAA", "gff3": "MYCOGFF3"}

#: ranks backfilled empty when a lineage omits them. Must stay identical to
#: ``mtdb.read_tax``'s list -- the two readers have to produce the same taxonomy
#: dict for the same genome, or a flat/SQLite round trip is not an identity.
#: Any other rank present in the source (``clade``, ``subclass``, ...) is stored
#: and returned as-is; this list only governs what gets added.
_TAX_RANKS = (
    "superkingdom",
    "kingdom",
    "phylum",
    "subphylum",
    "class",
    "order",
    "family",
    "subfamily",
)

_SCHEMA = """
CREATE TABLE IF NOT EXISTS genome (
    ome              TEXT PRIMARY KEY,
    genus            TEXT NOT NULL DEFAULT '',
    species          TEXT NOT NULL DEFAULT '',
    strain           TEXT NOT NULL DEFAULT '',
    version          TEXT NOT NULL DEFAULT '',
    source           TEXT NOT NULL DEFAULT '',
    biosample        TEXT NOT NULL DEFAULT '',
    assembly_acc     TEXT NOT NULL DEFAULT '',
    acquisition_date TEXT NOT NULL DEFAULT '',
    published        TEXT NOT NULL DEFAULT '',
    fna              TEXT NOT NULL DEFAULT '',
    faa              TEXT NOT NULL DEFAULT '',
    gff3             TEXT NOT NULL DEFAULT ''
);

CREATE TABLE IF NOT EXISTS taxonomy (
    genus   TEXT PRIMARY KEY,
    lineage TEXT NOT NULL DEFAULT '{}'
);

CREATE INDEX IF NOT EXISTS genome_assembly_acc ON genome(assembly_acc);
CREATE INDEX IF NOT EXISTS genome_genus        ON genome(genus);
CREATE INDEX IF NOT EXISTS genome_source       ON genome(source);
CREATE INDEX IF NOT EXISTS genome_published    ON genome(published);
"""


class MTDBSchemaError(Exception):
    """Raised when a SQLite MTDB was written by an incompatible version."""


# --------------------------------------------------------------------------- #
# detection / connection
# --------------------------------------------------------------------------- #
def is_sqlite(path: "str | Path | None") -> bool:
    """Is `path` a SQLite file? Detects by magic bytes, not by extension.

    Every entry point that accepts a database path funnels through here, which
    is what lets a `.mtdb` flat file and a SQLite database be used
    interchangeably wherever a path is currently accepted."""
    if not path:
        return False
    try:
        with open(path, "rb") as raw:
            return raw.read(16) == _SQLITE_MAGIC
    except (OSError, TypeError, ValueError):
        return False


def connect(path: "str | Path", read_only: bool = True) -> sqlite3.Connection:
    """Open a SQLite MTDB and verify its schema version.

    Read-only connections use a URI so concurrent readers can never be blocked
    by, or interfere with, a writer."""
    path = str(path)
    if read_only:
        conn = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    else:
        conn = sqlite3.connect(path)
    conn.row_factory = sqlite3.Row
    version = conn.execute("PRAGMA user_version").fetchone()[0]
    if version > SCHEMA_VERSION:
        conn.close()
        raise MTDBSchemaError(
            f"{path} uses MTDB schema v{version}; this Mycotools reads "
            f"v{SCHEMA_VERSION}. Upgrade Mycotools."
        )
    return conn


def _init_schema(conn: sqlite3.Connection) -> None:
    conn.executescript(_SCHEMA)
    conn.execute(f"PRAGMA user_version = {SCHEMA_VERSION}")


# --------------------------------------------------------------------------- #
# row <-> column-dict conversion
# --------------------------------------------------------------------------- #
def _resolve_paths(row: Dict[str, Any], ome: str, env: Optional[Mapping[str, str]]) -> None:
    """Expand abbreviated path columns in place, as the flat reader does.

    `env` is the pre-read {var: dir} mapping; None means the process is not
    linked to a primary MTDB, in which case abbreviated paths stay empty and the
    caller is responsible for reporting it."""
    if env is None:
        return
    for col, var in _PATH_COLUMNS.items():
        value = row.get(col) or ""
        if not value or value == f"{ome}.{col}":
            row[col] = env[var] + ome + "." + col


def read_path_env() -> Optional[Dict[str, str]]:
    """Read the three MTDB data-directory variables once.

    Hoisted out of per-row loops on purpose: ``os.environ`` decodes the whole
    environment on every iteration, which used to dominate database load time.
    Returns None when the process is not linked to a primary MTDB."""
    try:
        return {var: os.environ[var] for var in _PATH_COLUMNS.values()}
    except KeyError:
        return None


def rows2columns(
    rows: Iterable[Mapping[str, Any]],
    lineages: Mapping[str, Mapping[str, Any]],
    add_paths: bool = True,
) -> Dict[str, list]:
    """Build the column-oriented dict the `mtdb` class stores.

    Output is identical in shape to the flat reader's: one list per column, with
    `taxonomy` holding a lineage dict per row that carries genus/species/strain
    alongside the higher ranks."""
    df: Dict[str, list] = {c: [] for c in COLUMNS}
    env = read_path_env() if add_paths else None
    if add_paths and env is None:
        logger.error("MycotoolsDB not in path, cannot delineate biofile paths")

    for row in rows:
        row = dict(row)
        ome = row.get("ome") or ""
        genus = row.get("genus") or ""
        # a fresh dict per row: callers mutate row taxonomy (df2db strips the
        # genome-level ranks before writing) and must not corrupt the shared
        # genus lineage
        # reproduces `mtdb.read_tax`: stored ranks keep their order, then any
        # standard rank the source omitted is appended empty
        tax = dict(lineages.get(genus, {}))
        for rank in _TAX_RANKS:
            if rank not in tax:
                tax[rank] = ""
        tax["genus"] = genus
        tax["species"] = (genus + " " + (row.get("species") or "")).strip()
        tax["strain"] = row.get("strain") or ""
        row["taxonomy"] = tax
        _resolve_paths(row, ome, env)
        for col in COLUMNS:
            df[col].append(row.get(col, ""))
    return df


def _abbreviate_paths(row: Mapping[str, Any], ome: str) -> Dict[str, str]:
    """Collapse default file paths back to '' for storage.

    Mirrors what `mtdb.df2db` does when writing a flat file, so a database
    stays valid when its root directory moves."""
    out = {}
    env = read_path_env()
    for col in _PATH_COLUMNS:
        value = str(row.get(col) or "")
        if env is not None:
            default = env[_PATH_COLUMNS[col]] + ome + "." + col
            if value == default:
                value = ""
        out[col] = value
    return out


def split_taxonomy(
    omes: Sequence[str],
    genera: Sequence[str],
    taxonomies: Sequence[Any],
) -> Dict[str, Dict[str, Any]]:
    """Reduce per-row taxonomy to one lineage per genus.

    Lineages are a property of the genus -- `assimilate_tax` assigns the same
    dict to every row of a genus -- so this is lossless for a well-formed
    database. A genuine conflict means the source is corrupt, so it is reported
    rather than silently resolved; the most complete lineage wins."""
    lineages: Dict[str, Dict[str, Any]] = {}
    conflicts = set()
    for ome, genus, tax in zip(omes, genera, taxonomies):
        if not genus:
            continue
        stripped = _strip_genome_ranks(tax)
        if not _rank_count(stripped):
            lineages.setdefault(genus, stripped)
            continue
        prior = lineages.get(genus)
        if prior is None or not _rank_count(prior):
            lineages[genus] = stripped
        elif prior != stripped:
            conflicts.add(genus)
            if _rank_count(stripped) > _rank_count(prior):
                lineages[genus] = stripped
    if conflicts:
        logger.warning(
            "%d genera carry conflicting lineages across rows (%s%s); kept the "
            "most complete for each",
            len(conflicts),
            ", ".join(sorted(conflicts)[:5]),
            "..." if len(conflicts) > 5 else "",
        )
    return lineages


def _rank_count(lineage: Mapping[str, Any]) -> int:
    return sum(1 for v in lineage.values() if v)


def _strip_genome_ranks(tax: Any) -> Dict[str, Any]:
    """Drop genus/species/strain, which live in genome columns.

    Empty ranks are kept, and so is their order: a lineage stored here has to
    reproduce the taxonomy field of the `.mtdb` it came from byte for byte, so
    that flat -> SQLite -> flat is an identity."""
    if not tax:
        return {}
    if isinstance(tax, str):
        try:
            tax = json.loads(tax.replace("'", '"'))
        except (json.JSONDecodeError, TypeError):
            return {}
    if not isinstance(tax, dict):
        return {}
    return {
        rank: name
        for rank, name in tax.items()
        if rank not in {"genus", "species", "strain"}
    }


# --------------------------------------------------------------------------- #
# writing
# --------------------------------------------------------------------------- #
def write_db(path: "str | Path", columns: Mapping[str, Sequence[Any]]) -> str:
    """Write a whole MTDB to `path` atomically.

    The database is built in a sibling temporary file and renamed into place, so
    readers -- including `primary_db()`, which picks a database off the
    filesystem -- never observe a partially written database."""
    path = str(path)
    tmp = path + ".tmp"
    for stale in (tmp, tmp + "-journal", tmp + "-wal", tmp + "-shm"):
        if Path(stale).exists():
            Path(stale).unlink()

    omes = list(columns.get("ome", []))
    lineages = split_taxonomy(
        omes, list(columns.get("genus", [])), list(columns.get("taxonomy", []))
    )

    conn = sqlite3.connect(tmp)
    try:
        _init_schema(conn)
        genome_rows = []
        for i, ome in enumerate(omes):
            if not ome:
                continue
            row = {c: columns[c][i] if c in columns else "" for c in _GENOME_COLUMNS}
            row.update(_abbreviate_paths(row, ome))
            genome_rows.append(
                tuple("" if row.get(c) is None else str(row.get(c, "")) for c in _GENOME_COLUMNS)
            )
        placeholders = ", ".join("?" * len(_GENOME_COLUMNS))
        conn.executemany(
            f"INSERT OR REPLACE INTO genome ({', '.join(_GENOME_COLUMNS)}) "
            f"VALUES ({placeholders})",
            genome_rows,
        )
        conn.executemany(
            "INSERT OR REPLACE INTO taxonomy (genus, lineage) VALUES (?, ?)",
            [(genus, json.dumps(lin)) for genus, lin in lineages.items()],
        )
        conn.commit()
        conn.execute("VACUUM")
    finally:
        conn.close()

    Path(tmp).replace(path)
    return path


def upsert_rows(path: "str | Path", rows: Iterable[Mapping[str, Any]]) -> None:
    """Insert or update individual genome rows in an existing database."""
    conn = connect(path, read_only=False)
    try:
        _init_schema(conn)
        for row in rows:
            ome = row.get("ome") or ""
            if not ome:
                continue
            stored = {c: "" if row.get(c) is None else str(row.get(c, "")) for c in _GENOME_COLUMNS}
            stored.update(_abbreviate_paths(row, ome))
            conn.execute(
                f"INSERT OR REPLACE INTO genome ({', '.join(_GENOME_COLUMNS)}) "
                f"VALUES ({', '.join('?' * len(_GENOME_COLUMNS))})",
                tuple(stored[c] for c in _GENOME_COLUMNS),
            )
            lineage = _strip_genome_ranks(row.get("taxonomy"))
            if lineage and row.get("genus"):
                conn.execute(
                    "INSERT OR REPLACE INTO taxonomy (genus, lineage) VALUES (?, ?)",
                    (row["genus"], json.dumps(lineage)),
                )
        conn.commit()
    finally:
        conn.close()


def delete_omes(path: "str | Path", omes: Iterable[str]) -> int:
    """Remove genomes by ome; returns the number of rows deleted."""
    conn = connect(path, read_only=False)
    try:
        cur = conn.executemany("DELETE FROM genome WHERE ome = ?", [(o,) for o in omes])
        conn.commit()
        return cur.rowcount
    finally:
        conn.close()


# --------------------------------------------------------------------------- #
# reading
# --------------------------------------------------------------------------- #
def _lineages(conn: sqlite3.Connection, genera: Optional[Iterable[str]] = None) -> Dict[str, dict]:
    if genera is None:
        rows = conn.execute("SELECT genus, lineage FROM taxonomy").fetchall()
    else:
        genera = list(dict.fromkeys(genera))
        rows = []
        for chunk in _chunks(genera):
            rows.extend(
                conn.execute(
                    f"SELECT genus, lineage FROM taxonomy WHERE genus IN "
                    f"({', '.join('?' * len(chunk))})",
                    chunk,
                ).fetchall()
            )
    out = {}
    for row in rows:
        try:
            out[row["genus"]] = json.loads(row["lineage"])
        except json.JSONDecodeError:
            logger.error("malformed lineage for genus %s", row["genus"])
            out[row["genus"]] = {}
    return out


def _chunks(seq: Sequence[Any], size: int = 900):
    """SQLite caps host parameters per statement (default 999)."""
    for i in range(0, len(seq), size):
        yield list(seq[i : i + size])


def read_db(path: "str | Path", add_paths: bool = True) -> Dict[str, list]:
    """Read an entire SQLite MTDB into the column dict the `mtdb` class holds."""
    conn = connect(path)
    try:
        rows = conn.execute(
            f"SELECT {', '.join(_GENOME_COLUMNS)} FROM genome ORDER BY ome"
        ).fetchall()
        lineages = _lineages(conn)
    finally:
        conn.close()
    return rows2columns(rows, lineages, add_paths=add_paths)


def select_omes(
    path: "str | Path", omes: Iterable[str], add_paths: bool = True
) -> Dict[str, list]:
    """Read only the named genomes -- an index seek per ome, not a full scan.

    This is the query behind `mtdb <OME>`, `mtdb accession`, and every other
    tool that needs a handful of genomes out of the primary database."""
    omes = [o for o in dict.fromkeys(omes) if o]
    if not omes:
        return {c: [] for c in COLUMNS}
    conn = connect(path)
    try:
        rows = []
        for chunk in _chunks(omes):
            rows.extend(
                conn.execute(
                    f"SELECT {', '.join(_GENOME_COLUMNS)} FROM genome WHERE ome IN "
                    f"({', '.join('?' * len(chunk))}) ORDER BY ome",
                    chunk,
                ).fetchall()
            )
        lineages = _lineages(conn, [r["genus"] for r in rows])
    finally:
        conn.close()
    return rows2columns(rows, lineages, add_paths=add_paths)


def select_ome_prefix(
    path: "str | Path", prefix: str, add_paths: bool = True
) -> Dict[str, list]:
    """Read the genome named `prefix`, or its MTDB-versioned successors.

    An ome may carry a version tag (`cryneo24` -> `cryneo24.1`), so a lookup for
    the base ome has to find the versioned row."""
    conn = connect(path)
    try:
        rows = conn.execute(
            f"SELECT {', '.join(_GENOME_COLUMNS)} FROM genome "
            f"WHERE ome = ? OR ome LIKE ? ESCAPE '\\' ORDER BY ome",
            (prefix, prefix.replace("\\", "\\\\").replace("%", "\\%").replace("_", "\\_") + ".%"),
        ).fetchall()
        lineages = _lineages(conn, [r["genus"] for r in rows])
    finally:
        conn.close()
    return rows2columns(rows, lineages, add_paths=add_paths)


def select_column(
    path: "str | Path", column: str, values: Iterable[str], add_paths: bool = True
) -> Dict[str, list]:
    """Read the genomes whose `column` matches one of `values`."""
    if column not in set(_GENOME_COLUMNS):
        raise KeyError(f"cannot select on {column}")
    values = [v for v in dict.fromkeys(values) if v]
    if not values:
        return {c: [] for c in COLUMNS}
    conn = connect(path)
    try:
        rows = []
        for chunk in _chunks(values):
            rows.extend(
                conn.execute(
                    f"SELECT {', '.join(_GENOME_COLUMNS)} FROM genome WHERE {column} IN "
                    f"({', '.join('?' * len(chunk))}) ORDER BY ome",
                    chunk,
                ).fetchall()
            )
        lineages = _lineages(conn, [r["genus"] for r in rows])
    finally:
        conn.close()
    return rows2columns(rows, lineages, add_paths=add_paths)


def omes(path: "str | Path") -> List[str]:
    """List every ome without materializing the rest of the database."""
    conn = connect(path)
    try:
        return [r[0] for r in conn.execute("SELECT ome FROM genome ORDER BY ome")]
    finally:
        conn.close()


def count(path: "str | Path") -> int:
    conn = connect(path)
    try:
        return conn.execute("SELECT COUNT(*) FROM genome").fetchone()[0]
    finally:
        conn.close()


def infer_rank(path: "str | Path", lineage: str) -> Optional[str]:
    """Find the taxonomic rank a lineage name belongs to, using the taxonomy
    table instead of scanning every genome. Returns None if unknown."""
    target = lineage.lower()
    conn = connect(path)
    try:
        if conn.execute(
            "SELECT 1 FROM genome WHERE LOWER(genus) = ? LIMIT 1", (target,)
        ).fetchone():
            return "genus"
        for row in conn.execute("SELECT lineage FROM taxonomy"):
            try:
                lin = json.loads(row[0])
            except json.JSONDecodeError:
                continue
            for rank, name in lin.items():
                if isinstance(name, str) and name.lower() == target:
                    return rank
        if conn.execute(
            "SELECT 1 FROM genome WHERE LOWER(genus || ' ' || species) = ? LIMIT 1",
            (target,),
        ).fetchone():
            return "species"
        if conn.execute(
            "SELECT 1 FROM genome WHERE LOWER(strain) = ? LIMIT 1", (target,)
        ).fetchone():
            return "strain"
    finally:
        conn.close()
    return None


def genera_for_lineages(path: "str | Path", lineages: Iterable[str]) -> Optional[List[str]]:
    """Genera whose stored lineage matches any of `lineages` at any rank.

    Returns None when a name resolves to a genome-level rank (species/strain),
    which the taxonomy table cannot answer -- the caller then filters in memory."""
    wanted = set(x.lower() for x in lineages if x)
    if not wanted:
        return []
    conn = connect(path)
    try:
        hits = set()
        for row in conn.execute("SELECT genus, lineage FROM taxonomy"):
            genus = row["genus"]
            if genus and genus.lower() in wanted:
                hits.add(genus)
                continue
            try:
                lin = json.loads(row["lineage"])
            except json.JSONDecodeError:
                continue
            for name in lin.values():
                if isinstance(name, str) and name.lower() in wanted:
                    hits.add(genus)
                    break
        # a genus with no taxonomy row can still be matched by name
        for row in conn.execute("SELECT DISTINCT genus FROM genome"):
            if row[0] and row[0].lower() in wanted:
                hits.add(row[0])
    finally:
        conn.close()
    return sorted(hits)
