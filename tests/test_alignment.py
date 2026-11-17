
import pytest
from app.alignment import find_matches


def test_forward_strand():
    result = find_matches("ACGTACGTAGG", "ACGTACGT", 0)
    assert len(result) == 1
    assert result[0]["strand"] == "+"
    assert result[0]["pam"] == "AGG"


def test_reverse_strand():
    result = find_matches("CCTACGTACGT", "ACGTACGT", 0)
    assert len(result) == 1
    assert result[0]["strand"] == "-"


def test_zero_mismatches():
    result = find_matches("ACGTACGTAGG", "ACGTACGT", 0)
    assert result[0]["mismatches"] == 0


def test_one_mismatch():
    result = find_matches("ACGTTCGTAGG", "ACGTACGT", 1)
    assert result[0]["mismatches"] == 1


def test_invalid_guide():
    with pytest.raises(ValueError):
        find_matches("ACGTACGTAGG", "ACGTXYZ", 0)
