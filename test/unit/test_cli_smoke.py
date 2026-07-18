#! /usr/bin/env python3
"""CLI smoke tests across every mycotools entry-point script.

Two cheap, offline checks per script:

  * ``test_module_imports``   - the module imports without error (catches
    syntax errors, missing symbols, and version-incompatible type hints).
  * ``test_cli_help``         - ``python -m mycotools.<script> --help`` runs to
    completion without a Python traceback (catches import/parse-time crashes on
    the CLI path).

These are deliberately shallow: they prove each script is *invocable*, not that
its logic is correct. Per-script logic tests belong in dedicated
``test_<script>.py`` modules (see ``_template.py`` and the fuller
``test/update_mtdb/`` suite for the pattern).

Scripts that are currently broken are marked ``xfail`` so the suite stays green
while still documenting the breakage; when one is fixed its test will XPASS,
signalling that the entry should be removed from ``KNOWN_BROKEN``.
"""
import importlib
import re
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]

# --------------------------------------------------------------------------- #
# The user-facing CLI surface.
# --------------------------------------------------------------------------- #
# Modules registered as console entry points in pyproject.toml [project.scripts].
# Kept in sync with pyproject by test_entry_point_modules_match_pyproject below.
ENTRY_POINT_MODULES = [
    "mtdb",
    "acc2fa",
    "acc2gbk",
    "acc2gff",
    "acc2locus",
    "add2gff",
    "annotationStats",
    "assemblyStats",
    "bioreform",
    "coords2fa",
    "crap",
    "db2files",
    "db2hgs",
    "db2microsyntree",
    "db2search",
    "extract_mtdb",
    "fa2clus",
    "fa2hmmer2fa",
    "fa2mass",
    "fa2tree",
    "fna2faa",
    "gff2seq",
    "gff2svg",
    "jgiDwnld",
    "manage_mtdb",
    "ncbiAcc2fa",
    "ncbiDwnld",
    "ome2name",
    "predb2mtdb",
    "s2subs",
    "update_mtdb",
]

# Modules that ship a cli()/main() but are NOT registered as entry points.
EXTRA_CLI_MODULES = [
    "acc2fq",
    "ncbi_dwnld_fallback",
    "treetools",
]

ALL_MODULES = ENTRY_POINT_MODULES + EXTRA_CLI_MODULES

# module -> reason it currently fails to import / invoke. See test/README notes.
KNOWN_BROKEN = {
    "fa2hmmer2fa": "ImportError: cannot import name 'compAcc2fa' from mycotools.db2search",
    "s2subs": "entry point registered in pyproject but mycotools/s2subs.py does not exist",
    "treetools": "TypeError on `type | None` union hint; requires Python >= 3.10 (env is 3.9)",
}


def _params(modules):
    """Wrap module ids in pytest.param, attaching xfail marks for known breakage."""
    out = []
    for m in modules:
        marks = []
        if m in KNOWN_BROKEN:
            marks.append(pytest.mark.xfail(reason=KNOWN_BROKEN[m], strict=False))
        out.append(pytest.param(m, id=m, marks=marks))
    return out


@pytest.mark.parametrize("module", _params(ALL_MODULES))
def test_module_imports(module):
    """Every CLI module imports cleanly."""
    importlib.import_module(f"mycotools.{module}")


@pytest.mark.parametrize("module", _params(ALL_MODULES))
def test_cli_help(module, run_cli):
    """``--help`` runs without a Python traceback or a module-load failure.

    Some older scripts parse ``sys.argv`` by hand and do not implement argparse's
    ``--help`` (they may exit non-zero), so exit code is intentionally not
    asserted here - only that the process did not crash with a traceback.
    """
    res = run_cli(module, "--help")
    combined = res.stdout + res.stderr
    assert "Traceback (most recent call last)" not in res.stderr, combined
    assert "No module named" not in res.stderr, combined


def _registered_entry_point_modules():
    """Parse the ``[project.scripts]`` block of pyproject.toml without a TOML
    library (the test env is Python 3.9, which predates ``tomllib``)."""
    text = (REPO_ROOT / "pyproject.toml").read_text()
    block = re.search(r"^\[project\.scripts\](.*?)^\[", text, re.S | re.M)
    body = block.group(1) if block else text.split("[project.scripts]", 1)[-1]
    modules = set()
    for line in body.splitlines():
        m = re.match(r'\s*[\w.-]+\s*=\s*"mycotools\.([\w]+):', line)
        if m:
            modules.add(m.group(1))
    return modules


def test_entry_point_modules_match_pyproject():
    """Guard against drift: every pyproject console script is scaffolded here."""
    registered = _registered_entry_point_modules()
    assert registered, "failed to parse [project.scripts] from pyproject.toml"
    missing = registered - set(ENTRY_POINT_MODULES)
    assert not missing, f"new pyproject entry points not scaffolded: {sorted(missing)}"
