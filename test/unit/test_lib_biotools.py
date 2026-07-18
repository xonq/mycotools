#! /usr/bin/env python3
"""Unit tests for the pure parsers in ``mycotools.lib.biotools``.

These are login-free, filesystem-light functions - ideal unit-test targets and a
worked example of how to grow the per-script scaffolds beyond CLI smoke tests.
"""
import pytest

from mycotools.lib.biotools import fa2dict, dict2fa


@pytest.fixture
def sample_fasta(tmp_path):
    fa = tmp_path / "sample.fa"
    fa.write_text(
        ">seq1 a description here\n"
        "ACGTACGT\n"
        "AAAA\n"          # multi-line sequence should be concatenated
        ">seq2\n"
        "TTTT\n"
    )
    return fa


def test_fa2dict_parses_headers_and_sequences(sample_fasta):
    parsed = fa2dict(str(sample_fasta))
    assert set(parsed) == {"seq1", "seq2"}
    assert parsed["seq1"]["sequence"] == "ACGTACGTAAAA"   # wrapped lines joined
    assert parsed["seq1"]["description"] == "a description here"
    assert parsed["seq2"]["sequence"] == "TTTT"
    assert parsed["seq2"]["description"] == ""             # no description


def test_dict2fa_emits_expected_format(sample_fasta):
    parsed = fa2dict(str(sample_fasta))
    out = dict2fa(parsed)
    assert ">seq1 a description here" in out
    assert ">seq2" in out
    assert "ACGTACGTAAAA" in out


def test_fa2dict_dict2fa_round_trip(sample_fasta, tmp_path):
    parsed = fa2dict(str(sample_fasta))
    roundtrip_path = tmp_path / "roundtrip.fa"
    roundtrip_path.write_text(dict2fa(parsed))
    assert fa2dict(str(roundtrip_path)) == parsed


# --------------------------------------------------------------------------- #
# TODO(scaffold): extend coverage of the other pure biotools parsers.
# --------------------------------------------------------------------------- #
@pytest.mark.skip(reason="TODO: scaffold - add gff2list / list2gff round-trip tests")
def test_gff2list_list2gff_round_trip():
    ...
