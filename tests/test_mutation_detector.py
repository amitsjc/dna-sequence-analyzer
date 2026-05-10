"""Tests for the mutation detector."""

from dna_analyzer.mutation_detector import compare, _classify_mutation, SNP, Indel
from dna_analyzer.parser import Sequence


def make_seq(dna: str, seq_id: str = "test") -> Sequence:
    return Sequence(seq_id, "", dna)


# --- _classify_mutation ---

def test_transition_purine_to_purine():
    assert _classify_mutation("A", "G") == "transition"
    assert _classify_mutation("G", "A") == "transition"


def test_transition_pyrimidine_to_pyrimidine():
    assert _classify_mutation("C", "T") == "transition"
    assert _classify_mutation("T", "C") == "transition"


def test_transversion_purine_to_pyrimidine():
    assert _classify_mutation("A", "C") == "transversion"
    assert _classify_mutation("G", "T") == "transversion"


# --- compare(): identical sequences ---

def test_identical_sequences_no_mutations():
    seq = "ATGCATGCATGC"
    ref = make_seq(seq, "ref")
    query = make_seq(seq, "query")
    report = compare(ref, query)
    assert report.snps == []
    assert report.indels == []
    assert report.identity_percent == 100.0


# --- compare(): SNP detection ---

def test_single_snp_detected():
    ref   = make_seq("ATGCATGC", "ref")
    query = make_seq("ATGGATGC", "query")  # C→G at position 3
    report = compare(ref, query)
    assert len(report.snps) == 1
    assert report.snps[0].ref_base == "C"
    assert report.snps[0].alt_base == "G"


def test_snp_classified_correctly():
    ref   = make_seq("AAAAAAAA", "ref")
    query = make_seq("AAAGAAAA", "query")  # A→G = transition
    report = compare(ref, query)
    assert report.transitions == 1
    assert report.transversions == 0


def test_multiple_snps():
    ref   = make_seq("ATGATGATG", "ref")
    query = make_seq("ACGACGACG", "query")  # T→C at positions 1,4,7
    report = compare(ref, query)
    assert len(report.snps) == 3


# --- compare(): indels ---

def test_insertion_detected():
    ref   = make_seq("ATGCATGC", "ref")
    query = make_seq("ATGCAATGC", "query")  # extra A inserted
    report = compare(ref, query)
    assert any(i.kind == "insertion" for i in report.indels)


def test_deletion_detected():
    ref   = make_seq("ATGCATGC", "ref")
    query = make_seq("ATGATGC", "query")  # one C deleted
    report = compare(ref, query)
    assert any(i.kind == "deletion" for i in report.indels)


# --- report metadata ---

def test_report_ids_and_lengths():
    ref   = make_seq("ATGCATGC", "refseq")
    query = make_seq("ATGCATGC", "queryseq")
    report = compare(ref, query)
    assert report.ref_id == "refseq"
    assert report.query_id == "queryseq"
    assert report.ref_length == 8


def test_ts_tv_ratio_infinite_when_no_transversions():
    ref   = make_seq("AAAAAA", "ref")
    query = make_seq("AGAAAA", "query")  # only transition
    report = compare(ref, query)
    assert report.ts_tv_ratio == float("inf")


def test_summary_contains_key_fields():
    ref   = make_seq("ATGC", "ref")
    query = make_seq("AGGC", "query")
    report = compare(ref, query)
    summary = report.summary()
    assert "Identity" in summary
    assert "SNPs" in summary
    assert "Ts/Tv" in summary
