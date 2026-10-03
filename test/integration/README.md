# Live integration tests

These tests **hit the network** and **download real genome data**. They are the
opposite of the offline suites (`test/unit`, `test/update_mtdb`) and are marked
`@pytest.mark.integration`.

## Requirements

- The `mycotools` conda env (pandas, biopython, and the `datasets` executable on
  PATH).
- Network access to NCBI (and JGI, when its API works).
- **Stored no-password credentials** at `~/.mycotools/mtdb_credentials.json`
  (set via `mtdb-manage -s`). A missing store is a hard **error**, not a skip.

## Run

```bash
micromamba run -n mycotools pytest test/integration/ -v -s
```

`test_init_primary_db_downloads_jgi_and_ncbi` initializes a throwaway primary
MTDB from a 1-JGI + 1-NCBI subset of `test/ust.mtdb` and asserts genome data
from **both** sources is downloaded and curated (the finished primary `.mtdb`
must contain a `jgi`- and an `ncbi`-sourced row, each with `.fna`/`.faa`/`.gff3`
files). It runs in an **isolated `HOME`** (a temp dir holding only a copy of the
credential file) so the real `~/.mycotools/config.json` and your linked/active
database are never touched.

## JGI download API + tape restores

JGI genome downloads use the current **JGI Data Portal API** (see
`mycotools/jgiDwnld.py`): `mycocosm_file_list` (search) → `request_archived_files`
(restore) → `download_files` (zip stream), authorized with a signon session
token. The old `get-directory` XML endpoint was retired by JGI.

JGI keeps most files in tape archive (`file_status: PURGED`); the first download
of a genome requests a restore to disk, which "typically takes less than an hour
but can take up to a night". To stay deterministic, the test's fixture pre-warms
the restore for its JGI genome (bounded by `JGI_RESTORE_TIMEOUT`). If the restore
has not completed within that bound the test **skips** with a message to rerun
once JGI has staged the data, rather than blocking on tape I/O. Once a genome's
files are on disk they remain downloadable (no wait) until they re-purge.
