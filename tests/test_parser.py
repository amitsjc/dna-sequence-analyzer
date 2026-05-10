"""Tests for the FASTA parser."""

import tempfile
import os
import pytest
from dna_analyzer.parser import parse_fasta, Sequence


SAMPLE_FASTA = """>seq1 Test sequence one
ATGCATGCATGC
TTTTAAAA
>seq2 Test sequence two
GGGGCCCC
"""


def write_temp_fasta(content: str) -> str:
    f = tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False)
    f.write(content)
    f.close()
    return f.name


def test_parse_returns_correct_count():
    path = write_temp_fasta(SAMPLE_FASTA)
    seqs = parse_fasta(path)
    assert len(seqs) == 2
    os.unlink(path)


def test_sequence_ids_and_descriptions():
    path = write_temp_fasta(SAMPLE_FASTA)
    seqs = parse_fasta(path)
    assert seqs[0].id == "seq1"
    assert seqs[0].description == "Test sequence one"
    assert seqs[1].id == "seq2"
    os.unlink(path)


def test_multiline_sequence_joined():
    path = write_temp_fasta(SAMPLE_FASTA)
    seqs = parse_fasta(path)
    assert seqs[0].sequence == "ATGCATGCATGCTTTTAAAA"
    os.unlink(path)


def test_sequence_uppercased():
    path = write_temp_fasta(">s1 desc\natgcATGC\n")
    seqs = parse_fasta(path)
    assert seqs[0].sequence == "ATGCATGC"
    os.unlink(path)


def test_reverse_complement():
    seq = Sequence("s1", "", "ATGC")
    assert seq.reverse_complement() == "GCAT"


def test_is_valid_dna():
    assert Sequence("s1", "", "ATGC").is_valid()
    assert not Sequence("s1", "", "ATGX").is_valid()


def test_file_not_found():
    with pytest.raises(FileNotFoundError):
        parse_fasta("/nonexistent/path/file.fasta")
