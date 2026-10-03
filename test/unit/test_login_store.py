#! /usr/bin/env python3
"""Tests for the no-password credential store (dbtools.store_login / login_check).

All tests use temp paths and monkeypatch ``dbtools.PLAIN_LOGIN_PATH`` so the real
``~/.mycotools`` credential files are never read or written.
"""
import io
import json
import os
import stat

import pytest

import mycotools.lib.dbtools as dbtools


# a structurally representative JGI Data Portal token; not a live credential
TOKEN = "eyJhbGciOiJSUzI1NiJ9.eyJleHAiOjF9.c2lnbmF0dXJl"


@pytest.fixture
def paths(tmp_path):
    return {
        "plain": str(tmp_path / "mtdb_credentials.json"),
        "enc": str(tmp_path / "mtdb_key"),
    }


def test_store_login_round_trip(paths):
    dbtools.store_login(
        "APIKEY",
        TOKEN,
        info_path=paths["plain"],
        encrypted_path=paths["enc"],
    )
    assert dbtools.read_plain_login(paths["plain"]) == ("APIKEY", TOKEN)


def test_store_login_strips_bearer_prefix(paths):
    """A token is copied out of a JGI session with its scheme attached; the
    store holds the bare credential so callers concatenate exactly one."""
    dbtools.store_login(
        "APIKEY",
        "Bearer " + TOKEN,
        info_path=paths["plain"],
        encrypted_path=paths["enc"],
    )
    assert dbtools.read_plain_login(paths["plain"]) == ("APIKEY", TOKEN)


@pytest.mark.parametrize(
    "raw",
    [
        "Bearer " + TOKEN,
        "bearer " + TOKEN,
        "  Bearer   " + TOKEN + "  ",
        TOKEN + "\n",
    ],
)
def test_clean_jgi_token_normalizes(raw):
    assert dbtools.clean_jgi_token(raw) == TOKEN


def test_clean_jgi_token_passes_falsy_through():
    # callers distinguish None (never supplied) from "" (supplied empty)
    assert dbtools.clean_jgi_token(None) is None
    assert dbtools.clean_jgi_token("") == ""


def test_store_login_sets_owner_only_permissions(paths):
    dbtools.store_login(
        "K",
        TOKEN,
        info_path=paths["plain"],
        encrypted_path=paths["enc"],
    )
    assert stat.S_IMODE(os.stat(paths["plain"]).st_mode) == 0o600


def test_store_login_writes_valid_json_with_all_fields(paths):
    dbtools.store_login(
        "",
        TOKEN,
        info_path=paths["plain"],
        encrypted_path=paths["enc"],
    )
    with open(paths["plain"]) as fh:
        data = json.load(fh)
    assert set(data) == {"ncbi_api", "jgi_token"}
    assert data["ncbi_api"] == ""  # empty field preserved, not dropped


def test_store_login_removes_existing_encrypted_key(paths):
    with open(paths["enc"], "wb") as fh:
        fh.write(b"encrypted junk")
    dbtools.store_login(
        "K",
        TOKEN,
        info_path=paths["plain"],
        encrypted_path=paths["enc"],
    )
    assert not os.path.exists(paths["enc"])  # one source of truth


def test_logincheck_prefers_encrypted_key_over_plain(monkeypatch, paths):
    """When both stores exist, the encrypted key branch runs first. Its junk
    bytes fail to decrypt, and that failure proves the plaintext store (which
    would otherwise return cleanly) was never reached."""
    dbtools.store_login(
        "P",
        TOKEN,
        info_path=paths["plain"],
        encrypted_path=paths["enc"],
    )
    monkeypatch.setattr(dbtools, "PLAIN_LOGIN_PATH", paths["plain"])
    with open(paths["enc"], "wb") as fh:
        fh.write(b"not a valid fernet token")
    # feed a password through both the tty and non-tty code paths
    monkeypatch.setattr(dbtools.getpass, "getpass", lambda prompt="": "pw")
    monkeypatch.setattr("sys.stdin", io.StringIO("pw\n"))
    with pytest.raises(Exception):
        dbtools.login_check(info_path=paths["enc"])


def test_logincheck_uses_plain_store_when_no_encrypted_key(
    monkeypatch, tmp_path, paths
):
    dbtools.store_login(
        "APIKEY",
        TOKEN,
        info_path=paths["plain"],
        encrypted_path=paths["enc"],
    )
    monkeypatch.setattr(dbtools, "PLAIN_LOGIN_PATH", paths["plain"])
    missing_key = str(tmp_path / "no_such_key")
    assert dbtools.login_check(info_path=missing_key) == ("APIKEY", TOKEN)


def test_logincheck_prompts_when_no_store(monkeypatch, tmp_path):
    monkeypatch.setattr(dbtools, "PLAIN_LOGIN_PATH", str(tmp_path / "absent.json"))
    calls = {}

    def fake_getLogin(ncbi, jgi):
        calls["args"] = (ncbi, jgi)
        return ("K", TOKEN)

    monkeypatch.setattr(dbtools, "get_login", fake_getLogin)
    result = dbtools.login_check(
        info_path=str(tmp_path / "no_key"), ncbi=True, jgi=False
    )
    assert result == ("K", TOKEN)
    assert calls["args"] == (True, False)  # ncbi/jgi flags forwarded


def test_encrypt_pw_round_trips_token(monkeypatch, paths):
    """The encrypted key holds the same two fields the plain store does."""
    pytest.importorskip("cryptography")
    monkeypatch.setattr(dbtools, "PLAIN_LOGIN_PATH", paths["plain"])
    monkeypatch.setattr(dbtools.getpass, "getpass", lambda prompt="": "hunter2")
    dbtools.encrypt_pw("K", TOKEN, info_path=paths["enc"])
    monkeypatch.setattr("sys.stdin", io.StringIO("hunter2\n"))
    assert dbtools.login_check(info_path=paths["enc"]) == ("K", TOKEN)


def test_encrypt_pw_removes_plain_store(monkeypatch, paths):
    """Setting a password drops the unencrypted store, keeping one source."""
    pytest.importorskip("cryptography")
    dbtools.store_login(
        "K",
        TOKEN,
        info_path=paths["plain"],
        encrypted_path=paths["enc"],
    )
    monkeypatch.setattr(dbtools, "PLAIN_LOGIN_PATH", paths["plain"])
    monkeypatch.setattr(dbtools.getpass, "getpass", lambda prompt="": "hunter2")
    dbtools.encrypt_pw("K", TOKEN, info_path=paths["enc"])
    assert os.path.exists(paths["enc"])  # encrypted key written
    assert not os.path.exists(paths["plain"])  # plaintext store removed
