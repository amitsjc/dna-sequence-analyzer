"""Tests for the ORF finder."""

import pytest
from dna_analyzer.orf_finder import translate, find_orfs, longest_orf, ORF
from dna_analyzer.parser import Sequence


def make_seq(dna: str, seq_id: str = "test") -> Sequence:
    return Sequence(seq_id, "", dna)


# --- translate() ---

def test_translate_simple():
    # ATG (M) + AAA (K) + TAA (stop)
    assert translate("ATGAAATAA") == "MK"


def test_translate_stops_at_stop_codon():
    # Should not include anything after the stop
    assert translate("ATGAAATAAGGG") == "MK"


def test_translate_full_protein():
    # ATG-TTT-GGG-TAA → M-F-G-stop
    assert translate("ATGTTTGGGTAA") == "MFG"


def test_translate_unknown_codon():
    # NNN is not in the table → "?"
    assert "?" in translate("ATGNNNTAA")


# --- find_orfs() ---

def test_finds_orf_on_forward_strand():
    # Simple ORF: ATG...TAA in frame +1
    dna = "ATGAAAGGGCCCTAA" + "N" * 90  # pad to pass min_length=0
    seq = make_seq(dna)
    orfs = find_orfs(seq, min_length=0)
    assert any(o.strand == "+" for o in orfs)


def test_orf_has_correct_protein():
    # ATG-TTT-GGG-TAA → MFG, padded to exceed default min_length
    core = "ATGTTTGGGTAA"
    dna = core + "A" * 200
    seq = make_seq(dna)
    orfs = find_orfs(seq, min_length=0)
    forward_orfs = [o for o in orfs if o.strand == "+"]
    assert any(o.protein == "MFG" for o in forward_orfs)


def test_no_orf_when_no_start():
    seq = make_seq("TTTCCCGGG" * 20)
    orfs = find_orfs(seq, min_length=0)
    assert orfs == []


def test_min_length_filters_short_orfs():
    # Short ORF: 9 bp — should be excluded with min_length=100
    seq = make_seq("ATGAAATAA" + "N" * 200)
    orfs = find_orfs(seq, min_length=100)
    assert all(len(o) >= 100 for o in orfs)


def test_orfs_sorted_longest_first():
    seq = make_seq("ATGAAATAA" * 5 + "A" * 200)
    orfs = find_orfs(seq, min_length=0)
    lengths = [len(o) for o in orfs]
    assert lengths == sorted(lengths, reverse=True)


def test_longest_orf_returns_none_when_empty():
    seq = make_seq("TTTCCC" * 50)
    assert longest_orf(seq, min_length=0) is None


def test_orf_dataclass_repr():
    orf = ORF(start=0, end=12, frame=1, strand="+", nucleotides="ATGAAATAA", protein="MK")
    assert "frame=+1" in repr(orf)
