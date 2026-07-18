# mycotools unit-test scaffold

Offline, quick unit tests for the mycotools scripts. No JGI/NCBI login, no
downloads, no genome assimilation — every test exercises argument handling,
argparse surfaces, or login-free helper functions.

## Run

```bash
micromamba run -n mycotools pytest test/unit/ -q          # this scaffold
micromamba run -n mycotools pytest test/ -q               # scaffold + update_mtdb suite
```

(~9 s.)

## Layout

| File | What it is |
| --- | --- |
| `conftest.py` | Shared fixtures: `offline_env` (autouse; strips `MYCODB`, sets dummy `MYCO*` prefixes), path fixtures (`ust_mtdb`, `reference_mtdb`, `repo_root`, `test_data_dir`), and a `run_cli` subprocess helper. |
| `test_cli_smoke.py` | Data-driven **import** + **`--help` doesn't crash** checks for every entry-point script. The backbone: one parametrized case per script. |
| `test_lib_biotools.py` | A concrete, worked example of deeper unit tests (pure `fa2dict`/`dict2fa` parsers) — the pattern to copy for real logic. |
| `_template.py` | Copy to `test_<script>.py` and fill in the TODOs to grow a per-script suite. Not collected (leading underscore). |

The fuller, fully-implemented example lives one level up in
[`test/update_mtdb/`](../update_mtdb/) — argument-validation exit codes, CLI
parsing, and login-free helpers driven by `test/ust.mtdb`.

## How to extend a script

1. `cp test/unit/_template.py test/unit/test_<script>.py`
2. Point the import at the module, fill in the `--help` module name.
3. Add argument-validation tests (assert on `SystemExit.code` for bad flag
   combinations) and unit tests for the module's login-free helpers.
4. Use `ust_mtdb` / `reference_mtdb` fixtures for data; constrain heavier paths
   with a lineage or a tiny `.mtdb` to stay fast.

## Currently broken scripts (xfailed, not hidden)

`test_cli_smoke.py` marks these `xfail` so the suite stays green while
documenting the breakage. Fixing one turns its test **XPASS**, signalling that
the entry should be dropped from `KNOWN_BROKEN`:

| Script | Failure |
| --- | --- |
| `fa2hmmer2fa` | `ImportError: cannot import name 'compAcc2fa' from mycotools.db2search` |
| `s2subs` | registered as a console entry point in `pyproject.toml`, but `mycotools/s2subs.py` does not exist (`ModuleNotFoundError`) |
| `treetools` | `TypeError` on a `type | None` union hint — that syntax needs Python ≥ 3.10; the pinned env is 3.9 |

`test_entry_point_modules_match_pyproject` fails if a new `[project.scripts]`
entry is added without a scaffold entry, keeping the list from drifting.
