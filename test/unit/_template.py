#! /usr/bin/env python3
"""COPY-ME template for a per-script mycotools unit suite.

This file is NOT collected by pytest (leading underscore). To scaffold tests for
a script ``mycotools/<script>.py``:

    cp test/unit/_template.py test/unit/test_<script>.py

then fill in the TODOs. Shared fixtures (``offline_env`` autouse, ``run_cli``,
``ust_mtdb``, ``reference_mtdb``, ``repo_root``, ``test_data_dir``) come from
``test/unit/conftest.py`` - just declare them as arguments.

Guiding principles
------------------
* Stay OFFLINE. Do not exercise anything that reaches JGI/NCBI login or the
  network (Entrez queries, downloads, the ``datasets`` utility). Test argument
  handling and login-free helper functions only.
* Prefer calling functions directly over the subprocess CLI when you can, and
  assert on ``SystemExit.code`` for argument-validation paths. See
  ``test/update_mtdb/test_update_mtdb.py`` for a fully worked example.
* Use the committed ``test/ust.mtdb`` / ``test/reference.mtdb`` databases as
  fixtures; constrain any heavier path with a lineage or a tiny .mtdb so tests
  stay fast.
"""
import pytest

import mycotools.SCRIPT as mod  # TODO: rename import to the module under test


# --------------------------------------------------------------------------- #
# CLI surface
# --------------------------------------------------------------------------- #
def test_imports():
    """Module imports cleanly (already covered by test_cli_smoke, but keep a
    local anchor so this file is meaningful on its own)."""
    assert mod is not None


def test_cli_help(run_cli):
    res = run_cli("SCRIPT", "--help")  # TODO: module name
    assert "Traceback (most recent call last)" not in res.stderr


# --------------------------------------------------------------------------- #
# Argument validation (login-free)
# --------------------------------------------------------------------------- #
@pytest.mark.skip(reason="TODO: scaffold - assert exit codes for bad argument combos")
def test_argument_validation():
    """Example shape:

    with pytest.raises(SystemExit) as exc:
        mod.some_entry(bad_args...)
    assert exc.value.code == EXPECTED
    """
    ...


# --------------------------------------------------------------------------- #
# Pure helper functions (login-free)
# --------------------------------------------------------------------------- #
@pytest.mark.skip(reason="TODO: scaffold - unit test the module's pure helpers")
def test_pure_helpers(ust_mtdb): ...
