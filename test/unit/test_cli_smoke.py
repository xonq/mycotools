#! /usr/bin/env python3
"""CLI smoke tests across every mycotools entry-point and tool module.

Two cheap, offline checks per module:

  * ``test_module_imports``   - the module imports without error (catches
    syntax errors, missing symbols, and version-incompatible type hints).
  * ``test_cli_help``         - ``python -m mycotools.<module> --help`` runs to
    completion without a Python traceback (catches import/parse-time crashes on
    the CLI path).

These are deliberately shallow: they prove each tool is *invocable*, not that
its logic is correct. Per-script logic tests belong in dedicated
``test_<script>.py`` modules (see ``_template.py`` and the fuller
``test/update_mtdb/`` suite for the pattern).

The CLI surface is nested (``mycotools <GROUP> <TOOL>`` / ``mtdb <SUBCOMMAND>``),
so the module lists below mirror that shape: the two console scripts, the group
dispatchers they route through, and the leaf tool modules those dispatchers
import and call ``cli()`` on. Three drift guards keep the lists honest -
``test_entry_point_modules_match_pyproject``,
``test_deprecated_shims_target_scaffolded_modules``, and
``test_dispatch_targets_are_scaffolded`` - so a tool added to a dispatch table
or pyproject without being scaffolded here fails rather than going untested.

Modules that are currently broken are marked ``xfail`` via ``KNOWN_BROKEN`` so
the suite stays green while still documenting the breakage; when one is fixed
its test XPASSes, signalling that the entry should be removed. Modules gated on
an optional third-party dependency are listed in ``OPTIONAL_DEPS`` and skip
(rather than fail) when it is not installed.
"""
import functools
import importlib
import importlib.util
import re
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]

# --------------------------------------------------------------------------- #
# The user-facing CLI surface.
# --------------------------------------------------------------------------- #
# The two console scripts in pyproject.toml [project.scripts]: `mtdb` is the
# database lifecycle, `mycotools` (mycotools.cli) the downstream analysis
# toolkit. Kept in sync with pyproject by the drift guards below.
PRIMARY_MODULES = [
    "mtdb",
    "cli",
]

# Group dispatchers: each parses a single positional plus REMAINDER and hands
# the rest of argv to a leaf module (see mycotools/lib/subcmd.py). They are not
# console scripts of their own, but each is importable and `python -m`-invocable.
GROUP_MODULES = [
    "download",
    "homology",
    "cluster",
    "phylo",
    "stats",
    "seq",
    "gff",
    "mtdb.acc2",
]

# Leaf tools: the modules a dispatcher actually imports and calls `cli()` on.
LEAF_MODULES = [
    # mycotools <GROUP> <TOOL>
    "download.jgi",
    "download.ncbi",
    "homology.db",
    "homology.fasta",
    "cluster.db",
    "cluster.fasta",
    "phylo.crap",
    "phylo.tree",
    "phylo.synteny",
    "phylo.tools",
    "stats.annotation",
    "stats.assembly",
    "seq.translate",
    "seq.coords",
    "seq.gff",
    "seq.convert",
    "seq.mass",
    "gff.add",
    "gff.svg",
    "rename",
    # mtdb <SUBCOMMAND>
    "mtdb.extract",
    "mtdb.update",
    "mtdb.predb",
    "mtdb.configure",
    "mtdb.manage",
    "mtdb.files",
    # mtdb accession <FORMAT>
    "mtdb.acc2.fa",
    "mtdb.acc2.gff",
    "mtdb.acc2.gbk",
    "mtdb.acc2.locus",
]

# Everything with a CLI behind it.
CLI_MODULES = PRIMARY_MODULES + GROUP_MODULES + LEAF_MODULES

# `mycotools.deprecated` holds the pre-v2 flat console-script shims. It exposes
# no CLI of its own (each shim warns and forwards to a leaf module's `cli()`),
# so it is import-checked but not `--help`-checked.
ALL_MODULES = CLI_MODULES + ["deprecated"]

# module -> reason it currently fails to import / invoke. See test/README notes.
KNOWN_BROKEN = {}

# module -> third-party package it imports at module scope. An optional
# dependency that is absent (or, like ete3 on Python >= 3.13, installed but not
# importable) skips rather than fails: the module is fine, the environment is
# not the full runtime.
OPTIONAL_DEPS = {
    "phylo.crap": "ete3",
}


@functools.lru_cache(maxsize=None)
def _dep_unusable(dep):
    """Return a reason string if `dep` cannot be imported here, else None.

    Presence is not enough: ete3 ships a `cgi` import that Python 3.13 removed,
    so it resolves as installed and then raises on import."""
    if importlib.util.find_spec(dep) is None:
        return f"optional dependency {dep!r} is not installed"
    try:
        importlib.import_module(dep)
    except Exception as error:  # noqa: BLE001 - any import failure disqualifies it
        return f"optional dependency {dep!r} is not importable: {error}"
    return None


def _params(modules):
    """Wrap module ids in pytest.param, attaching xfail/skipif marks."""
    out = []
    for m in modules:
        marks = []
        if m in KNOWN_BROKEN:
            marks.append(pytest.mark.xfail(reason=KNOWN_BROKEN[m], strict=False))
        dep = OPTIONAL_DEPS.get(m)
        if dep is not None:
            unusable = _dep_unusable(dep)
            marks.append(
                pytest.mark.skipif(unusable is not None, reason=unusable or "")
            )
        out.append(pytest.param(m, id=m, marks=marks))
    return out


@pytest.mark.parametrize("module", _params(ALL_MODULES))
def test_module_imports(module):
    """Every CLI module imports cleanly."""
    importlib.import_module(f"mycotools.{module}")


@pytest.mark.parametrize("module", _params(CLI_MODULES))
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


@pytest.mark.parametrize("module", _params(LEAF_MODULES))
def test_leaf_modules_expose_cli(module):
    """Dispatchers and deprecation shims invoke leaves as ``module.cli()``.

    A leaf that loses its ``cli`` still imports and still passes ``--help``
    (``python -m`` runs ``__main__``/module scope), so only this catches it -
    and it breaks every nested invocation of the tool."""
    mod = importlib.import_module(f"mycotools.{module}")
    assert callable(getattr(mod, "cli", None)), f"mycotools.{module} has no cli()"


# --------------------------------------------------------------------------- #
# Drift guards: the lists above must cover the real, declared CLI surface.
# --------------------------------------------------------------------------- #
_SCRIPT_RE = re.compile(r'^\s*([\w.-]+)\s*=\s*"mycotools\.([\w.]+):(\w+)"')

# console script -> (module within `mycotools`, attribute) for the two real
# entrypoints. Everything else in [project.scripts] is a deprecated flat alias.
PRIMARY_ENTRY_POINTS = {
    "mtdb": ("mtdb", "main"),
    "mycotools": ("cli", "main"),
}


def _registered_entry_points():
    """Parse ``[project.scripts]`` into ``{script: (module, attr)}``.

    Done with a regex rather than a TOML library so the suite still runs on the
    Python 3.9 reference environment, which predates ``tomllib``."""
    text = (REPO_ROOT / "pyproject.toml").read_text()
    block = re.search(r"^\[project\.scripts\](.*?)^\[", text, re.S | re.M)
    body = block.group(1) if block else text.split("[project.scripts]", 1)[-1]
    scripts = {}
    for line in body.splitlines():
        m = _SCRIPT_RE.match(line)
        if m:
            scripts[m.group(1)] = (m.group(2), m.group(3))
    return scripts


def test_entry_point_modules_match_pyproject():
    """Guard against drift: every pyproject console script resolves to a module
    scaffolded here - the two real entrypoints, or the deprecation shim module."""
    registered = _registered_entry_points()
    assert registered, "failed to parse [project.scripts] from pyproject.toml"
    unexpected = {m for m, _ in registered.values()} - {"mtdb", "cli", "deprecated"}
    assert not unexpected, f"new pyproject entry point modules not scaffolded: {sorted(unexpected)}"


def test_primary_entry_points_resolve():
    """`mtdb` and `mycotools` point at a callable that actually exists."""
    registered = _registered_entry_points()
    for script, target in PRIMARY_ENTRY_POINTS.items():
        assert registered.get(script) == target, (
            f"{script} should be mycotools.{target[0]}:{target[1]}, "
            f"got {registered.get(script)}"
        )
        module, attr = target
        assert callable(getattr(importlib.import_module(f"mycotools.{module}"), attr))


def test_deprecated_entry_points_resolve():
    """Every flat alias resolves to a shim that exists in mycotools.deprecated.

    A console script naming a missing attribute installs fine and only fails when
    the user runs it, so nothing but this catches a stale alias."""
    deprecated = importlib.import_module("mycotools.deprecated")
    aliases = {
        script: attr
        for script, (module, attr) in _registered_entry_points().items()
        if module == "deprecated"
    }
    assert aliases, "no deprecated flat aliases parsed from [project.scripts]"
    for script, attr in sorted(aliases.items()):
        shim = getattr(deprecated, attr, None)
        assert callable(shim), f"`{script}` resolves to a missing shim {attr!r}"


def _shim_target(shim):
    """Return the module a deprecation shim delegates to.

    `_make_shim` closes over the target module path; reading it back off the
    closure keeps the shim table itself as the single source of truth."""
    free = dict(zip(shim.__code__.co_freevars, shim.__closure__ or ()))
    cell = free.get("target")
    return cell.cell_contents if cell is not None else None


def test_deprecated_shims_target_scaffolded_modules():
    """Each flat alias forwards to a relocated module that is smoke-tested here."""
    deprecated = importlib.import_module("mycotools.deprecated")
    scaffolded = {f"mycotools.{m}" for m in ALL_MODULES}
    for script, (module, attr) in sorted(_registered_entry_points().items()):
        if module != "deprecated":
            continue
        target = _shim_target(getattr(deprecated, attr))
        assert target, f"could not resolve the delegate target of `{script}`"
        assert target in scaffolded, f"`{script}` forwards to unscaffolded {target}"


def _dispatch_targets():
    """Resolve every nested subcommand to the module it dispatches to."""
    targets = set()
    # `mycotools <GROUP>` and `mtdb <SUBCOMMAND>` are the two roots; each group
    # dispatcher then carries its own table of leaves.
    roots = [("mycotools", "cli"), ("mycotools.mtdb", "mtdb")]
    for package, module in roots:
        table = importlib.import_module(f"mycotools.{module}").SUBCOMMANDS
        targets.update(f"{package}.{sub}" for sub in table.values())
    for group in GROUP_MODULES:
        table = importlib.import_module(f"mycotools.{group}").SUBCOMMANDS
        targets.update(f"mycotools.{group}.{sub}" for sub in table.values())
    return targets


def test_dispatch_targets_are_scaffolded():
    """Guard against drift: a tool added to a dispatch table is tested here.

    The nested dispatchers import their leaves lazily, so a subcommand pointing
    at a module that does not exist stays invisible until someone runs it."""
    scaffolded = {f"mycotools.{m}" for m in ALL_MODULES}
    missing = _dispatch_targets() - scaffolded
    assert not missing, f"dispatch targets not scaffolded here: {sorted(missing)}"
