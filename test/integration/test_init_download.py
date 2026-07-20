#! /usr/bin/env python3
"""LIVE integration test: initialize a primary MycotoolsDB from test/ust.mtdb.

Unlike the offline suites (test/unit, test/update_mtdb), this test actually
downloads JGI (MycoCosm) and NCBI genome data and gathers NCBI taxonomy, using
the *no-password* credential store (``~/.mycotools/mtdb_credentials.json``, see
dbtools.store_login). It therefore requires:

  * network access to JGI and NCBI,
  * the `datasets` executable on PATH (ships in the `mycotools` conda env), and
  * stored, accessible NCBI/JGI credentials.

Per request, a missing credential store is a hard ERROR (not a skip): the test
raises so the absence is loud.

JGI downloads use the current JGI Data Portal API (see mycotools.jgiDwnld). JGI
archives most files to tape (file_status PURGED); a first download restores them
to disk, which "typically takes less than an hour but can take up to a night".
To keep this test deterministic and fast, its fixture pre-warms the restore for
the JGI genome (bounded); if the restore has not completed within that bound the
test skips (rerun once JGI has staged the data) rather than block on tape I/O.

Isolation: the run happens in a throwaway ``HOME`` that contains only a copy of
the credential file - never the real ``~/.mycotools/config.json`` - so the
user's linked/active database is never modified. update_mtdb writes its new
config into the temp HOME, which is deleted on teardown.

Run:

    micromamba run -n mycotools pytest test/integration/ -v -s

To bound runtime, only a small subset of ust.mtdb is initialized (see
N_PER_SOURCE); one JGI + one NCBI entry still exercises the full download and
curation paths for both sources.
"""
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from mycotools.lib.dbtools import read_plain_login

REPO_ROOT = Path(__file__).resolve().parents[2]
UST_MTDB = REPO_ROOT / "test" / "ust.mtdb"
CREDS_PATH = Path.home() / ".mycotools" / "mtdb_credentials.json"

# entries per source (jgi/ncbi) to pull from ust.mtdb into the reference; small
# by default so the live download stays quick.
N_PER_SOURCE = 1

# column indices in an .mtdb row (0-based): ome, genus, species, strain,
# taxonomy, version, source, biosample, assembly_acc, ...
_SOURCE_COL = 6
_ASSEMBLY_ACC_COL = 8

# how long to wait for JGI to restore archived (on-tape) files before skipping
JGI_RESTORE_TIMEOUT = 1200


def require_credentials():
    """Return (ncbi_email, ncbi_api, jgi_email, jgi_pwd) from the no-password
    store, raising if the store or any required field is absent."""
    if not CREDS_PATH.is_file():
        raise RuntimeError(
            f"no-password credential store not found at {CREDS_PATH}. Store "
            "credentials first with `mtdb manage --store` (or `mtdb m -s`)."
        )
    ncbi_email, ncbi_api, jgi_email, jgi_pwd = read_plain_login(str(CREDS_PATH))
    required = {"ncbi_email": ncbi_email, "jgi_email": jgi_email, "jgi_pwd": jgi_pwd}
    missing = sorted(name for name, val in required.items() if not val)
    if missing:
        raise RuntimeError(
            f"credential store {CREDS_PATH} is missing required field(s): {missing}"
        )
    return ncbi_email, ncbi_api, jgi_email, jgi_pwd


def _iter_rows(source=None):
    """Yield the tab-split fields of each data row in ust.mtdb, optionally
    filtered to a single source (jgi/ncbi)."""
    for line in UST_MTDB.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if source is None or fields[_SOURCE_COL].strip().lower() == source:
            yield fields


def build_reference_subset(dest: Path, n_per_source: int = N_PER_SOURCE) -> Path:
    """Write a small reference .mtdb containing n JGI + n NCBI entries taken from
    ust.mtdb (preserving its exact columns)."""
    rows = []
    for src in ("jgi", "ncbi"):
        rows.extend(["\t".join(f) for f in list(_iter_rows(src))[:n_per_source]])
    assert any("\tjgi\t" in r or r.split("\t")[_SOURCE_COL] == "jgi" for r in rows)
    assert any(r.split("\t")[_SOURCE_COL] == "ncbi" for r in rows)
    dest.write_text("\n".join(rows) + "\n")
    return dest


def jgi_portal_ids(n_per_source: int = N_PER_SOURCE):
    """The JGI portal ids (assembly_acc) that build_reference_subset will use."""
    return [f[_ASSEMBLY_ACC_COL].strip() for f in list(_iter_rows("jgi"))[:n_per_source]]


def warm_jgi_restore(portal_ids, jgi_email, jgi_pwd, timeout=JGI_RESTORE_TIMEOUT):
    """Ensure the assembly + gff3 for each JGI portal id are RESTORED (on disk),
    requesting a tape restore and polling if needed. Dogfoods the jgiDwnld API
    helpers so the warmed selection matches what update_mtdb will download.
    Returns the set of portal ids whose essential files are ready."""
    from mycotools.jgiDwnld import (
        jgi_api_login, search_organism, select_file, request_restore,
        poll_restore, _mycocosm_ids, _is_restored,
    )

    session, token = jgi_api_login(jgi_email, jgi_pwd)
    ready = set()
    for portal_id in portal_ids:
        org, files = search_organism(session, portal_id)
        if not org:
            continue
        selected = [
            f for f in (select_file(files, "fna", masked=True),
                        select_file(files, "gff3", masked=True))
            if f is not None
        ]
        if len(selected) < 2:  # need both an assembly and an annotation
            continue
        purged = [f["_id"] for f in selected if not _is_restored(f)]
        if purged:
            status_url = request_restore(
                session, token,
                _mycocosm_ids(org.get("id"), (org.get("top_hit") or {}).get("_id"),
                              org.get("mycocosm_portal_id") or portal_id, purged),
            )
            if not poll_restore(session, status_url, timeout=timeout, interval=20):
                continue
        ready.add(portal_id)
    return ready


@pytest.fixture
def isolated_home(tmp_path):
    """A throwaway HOME containing only a copy of the real credential store."""
    require_credentials()  # hard error if creds are missing
    home = tmp_path / "home"
    (home / ".mycotools").mkdir(parents=True)
    dst = home / ".mycotools" / "mtdb_credentials.json"
    shutil.copy(CREDS_PATH, dst)
    os.chmod(dst, 0o600)
    return home


def run_update_mtdb(home: Path, init_dir: Path, reference: Path, timeout=1800):
    """Invoke `update_mtdb --init <init_dir> --reference <reference>` in an
    isolated HOME, with MYCODB stripped so no linked DB is inherited."""
    env = dict(os.environ)
    env["HOME"] = str(home)
    for var in ("MYCODB", "MYCOFNA", "MYCOFAA", "MYCOGFF3"):
        env.pop(var, None)
    return subprocess.run(
        [sys.executable, "-m", "mycotools.mtdb.update",
         "--init", str(init_dir), "--reference", str(reference)],
        stdin=subprocess.DEVNULL,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        timeout=timeout,
        env=env,
    )


def read_primary_sources(init_dir: Path):
    """Return the list of `source` values in the produced primary .mtdb."""
    mtdb_files = list(init_dir.glob("**/mtdb/*.mtdb"))
    assert mtdb_files, "no <date>.mtdb produced"
    primary = max(mtdb_files, key=lambda p: p.stat().st_size)
    sources = []
    for line in primary.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) > _SOURCE_COL:
            sources.append(fields[_SOURCE_COL].strip().lower())
    return primary, sources


@pytest.mark.integration
def test_credentials_present():
    """The no-password credential store exists and has NCBI + JGI logins."""
    ncbi_email, ncbi_api, jgi_email, jgi_pwd = require_credentials()
    assert "@" in ncbi_email
    assert "@" in jgi_email
    assert jgi_pwd


@pytest.mark.integration
def test_init_primary_db_downloads_jgi_and_ncbi(isolated_home, tmp_path):
    """Initialize a primary MTDB from a ust.mtdb subset (1 JGI + 1 NCBI) and
    confirm genome data from *both* sources is downloaded and curated.

    JGI downloads use the current JGI Data Portal API; the fixture warms the
    tape restore first so the download is deterministic. Both the JGI and NCBI
    genomes must appear in the finished primary database, each with sequence and
    annotation files, and the init must never die on an XML-parse error."""
    _, _, jgi_email, jgi_pwd = require_credentials()

    portal_ids = jgi_portal_ids()
    ready = warm_jgi_restore(portal_ids, jgi_email, jgi_pwd)
    if set(portal_ids) - ready:
        pytest.skip(
            "JGI has not finished restoring "
            f"{sorted(set(portal_ids) - ready)} from tape; rerun once staged."
        )

    reference = build_reference_subset(tmp_path / "reference_subset.mtdb")
    init_dir = isolated_home / "mtdb_init"   # non-existent -> becomes the DB root

    result = run_update_mtdb(isolated_home, init_dir, reference)

    # surface the full log on failure so download / curation errors are visible
    assert result.returncode == 0, (
        f"update_mtdb --init exited {result.returncode}\n"
        f"----- output -----\n{result.stdout}"
    )

    # regression guard: the JGI directory data must never crash the run
    assert "Traceback (most recent call last)" not in result.stdout, result.stdout
    assert "ParseError" not in result.stdout, result.stdout

    # the finished primary database contains BOTH sources
    primary, sources = read_primary_sources(init_dir)
    assert sources, f"primary MTDB is empty\n{result.stdout}"
    assert "jgi" in sources, f"JGI genome missing from primary DB\n{result.stdout}"
    assert "ncbi" in sources, f"NCBI genome missing from primary DB\n{result.stdout}"

    # downloaded + curated sequence/annotation data (>=2 genomes: 1 jgi + 1 ncbi)
    fna = list(init_dir.glob("**/data/fna/*.fna"))
    faa = list(init_dir.glob("**/data/faa/*.faa"))
    gff3 = list(init_dir.glob("**/data/gff3/*.gff3"))
    assert len(fna) >= 2, f"expected >=2 genome assemblies (.fna)\n{result.stdout}"
    assert len(faa) >= 2, f"expected >=2 proteomes (.faa)\n{result.stdout}"
    assert len(gff3) >= 2, f"expected >=2 annotations (.gff3)\n{result.stdout}"
