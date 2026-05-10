"""Tests for GC content analysis."""

from dna_analyzer.gc_content import gc_content, gc_content_windows, gc_skew, nucleotide_counts, summarize
from dna_analyzer.parser import Sequence


def test_gc_content_all_gc():
    assert gc_content("GGCC") == 100.0


def test_gc_content_all_at():
    assert gc_content("AATT") == 0.0


def test_gc_content_mixed():
    # ATGC → 2 GC out of 4 = 50%
    assert gc_content("ATGC") == 50.0


def test_gc_content_empty():
    assert gc_content("") == 0.0


def test_gc_content_case_insensitive():
    assert gc_content("atgc") == gc_content("ATGC")


def test_gc_skew_positive():
    # More G than C → positive skew
    assert gc_skew("GGGC") > 0


def test_gc_skew_negative():
    assert gc_skew("GCCC") < 0


def test_gc_skew_balanced():
    assert gc_skew("GGCC") == 0.0


def test_gc_skew_no_gc():
    assert gc_skew("AAAA") == 0.0


def test_nucleotide_counts():
    counts = nucleotide_counts("AATGCCN")
    assert counts == {"A": 2, "T": 1, "G": 1, "C": 2, "N": 1}


def test_sliding_window_count():
    # 20bp sequence, window=10, step=5 → windows at 0,5,10
    seq = "A" * 10 + "G" * 10
    windows = gc_content_windows(seq, window_size=10, step=5)
    assert len(windows) == 3
    assert windows[0]["gc_percent"] == 0.0   # all A
    assert windows[2]["gc_percent"] == 100.0  # all G


def test_summarize_fields():
    seq = Sequence("test", "Test sequence", "ATGCATGC")
    result = summarize(seq)
    assert result["id"] == "test"
    assert result["length"] == 8
    assert result["gc_percent"] == 50.0
    assert result["at_percent"] == 50.0
    assert "nucleotide_counts" in result
