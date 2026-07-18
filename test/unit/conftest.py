#! /usr/bin/env python3
"""Shared fixtures for the mycotools unit-test scaffold.

Everything here keeps tests *offline*: no JGI/NCBI login, no downloads. The
per-script suites built on top of these fixtures should only exercise argument
handling, argparse surfaces, and login-free helper functions.
"""
import os
import subprocess
import sys
from pathlib import Path

import pytest

# test/  (this file lives in test/unit/)
TEST_DIR = Path(__file__).resolve().parents[1]
REPO_ROOT = TEST_DIR.parent


@pytest.fixture(autouse=True)
def offline_env(monkeypatch):
    """Force the login-free code paths for every test.

    * Dropping ``MYCODB`` makes ``dbtools.primaryDB()`` return ``None`` instead
      of resolving a real, linked database.
    * ``MYCOFNA``/``MYCOFAA``/``MYCOGFF3`` are only string-concatenated by the
      pandas MTDB helpers, so dummy prefixes satisfy them without touching disk.
    """
    monkeypatch.delenv("MYCODB", raising=False)
    monkeypatch.setenv("MYCOFNA", "/tmp/mtdb_unittest_fna/")
    monkeypatch.setenv("MYCOFAA", "/tmp/mtdb_unittest_faa/")
    monkeypatch.setenv("MYCOGFF3", "/tmp/mtdb_unittest_gff3/")


@pytest.fixture(scope="session")
def repo_root():
    return REPO_ROOT


@pytest.fixture(scope="session")
def test_data_dir():
    return TEST_DIR


@pytest.fixture(scope="session")
def ust_mtdb():
    """Path to the committed 42-entry Ustilaginomycotina test database."""
    return TEST_DIR / "ust.mtdb"


@pytest.fixture(scope="session")
def reference_mtdb():
    """Path to the committed 9-entry reference test database."""
    return TEST_DIR / "reference.mtdb"


@pytest.fixture
def run_cli():
    """Return a helper that runs ``python -m mycotools.<module> <args>`` offline.

    The child process gets ``stdin`` closed (so an unexpected ``input()`` prompt
    fails fast instead of hanging) and an environment with ``MYCODB`` stripped.
    """
    def _run(module, *args, timeout=60):
        env = dict(os.environ)
        env.pop("MYCODB", None)
        env.setdefault("MYCOFNA", "/tmp/mtdb_unittest_fna/")
        env.setdefault("MYCOFAA", "/tmp/mtdb_unittest_faa/")
        env.setdefault("MYCOGFF3", "/tmp/mtdb_unittest_gff3/")
        return subprocess.run(
            [sys.executable, "-m", f"mycotools.{module}", *args],
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=timeout,
            env=env,
        )

    return _run
