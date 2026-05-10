"""Smoke tests for the CLI entry points."""

import subprocess
import sys
import os

PYTHON = sys.executable
FASTA = os.path.join(os.path.dirname(__file__), "..", "data", "example.fasta")


def run(*args):
    result = subprocess.run(
        [PYTHON, "-m", "dna_analyzer", *args],
        capture_output=True, text=True,
        cwd=os.path.join(os.path.dirname(__file__), "..")
    )
    return result


def test_analyze_exits_zero():
    r = run("analyze", FASTA, "--min-orf", "0")
    assert r.returncode == 0


def test_analyze_output_contains_gc():
    r = run("analyze", FASTA, "--min-orf", "0")
    assert "GC Content" in r.stdout


def test_orfs_command():
    r = run("orfs", FASTA, "--min-orf", "0")
    assert r.returncode == 0
    assert "Open Reading Frames" in r.stdout


def test_compare_command():
    r = run("compare", FASTA, FASTA, "--query-index", "1")
    assert r.returncode == 0
    assert "Identity" in r.stdout
    assert "99.51" in r.stdout


def test_missing_file_exits_nonzero():
    r = run("analyze", "nonexistent.fasta")
    assert r.returncode != 0


def test_help_flag():
    r = run("--help")
    assert r.returncode == 0
    assert "analyze" in r.stdout
    assert "compare" in r.stdout
