"""Command-line interface — three subcommands: analyze, orfs, compare."""

import argparse
import sys
from dna_analyzer.parser import parse_fasta
from dna_analyzer.gc_content import summarize, gc_content_windows
from dna_analyzer.orf_finder import find_orfs
from dna_analyzer.mutation_detector import compare
from dna_analyzer import report


def cmd_analyze(args):
    """Full analysis: GC content + ORFs for every sequence in a FASTA file."""
    report.print_welcome()
    sequences = parse_fasta(args.fasta)
    print(f"  Loaded {len(sequences)} sequence(s) from {args.fasta}\n")

    for seq in sequences:
        gc_summary = summarize(seq)
        report.print_gc_report(gc_summary)

        if args.windows:
            windows = gc_content_windows(seq.sequence, window_size=args.window_size, step=args.step)
            print(f"\n  Sliding window GC% (size={args.window_size}, step={args.step}):")
            for w in windows:
                bar_len = int(w['gc_percent'] / 5)
                bar = "█" * bar_len
                print(f"    {w['start']:>5}–{w['end']:<5} {w['gc_percent']:>5.1f}%  {bar}")

        orfs = find_orfs(seq, min_length=args.min_orf)
        report.print_orf_report(seq.id, orfs)

    print()


def cmd_orfs(args):
    """Find and display ORFs in a FASTA file."""
    sequences = parse_fasta(args.fasta)
    for seq in sequences:
        orfs = find_orfs(seq, min_length=args.min_orf)
        report.print_orf_report(seq.id, orfs)
    print()


def cmd_compare(args):
    """Compare two sequences from separate FASTA files."""
    ref_seqs = parse_fasta(args.ref)
    query_seqs = parse_fasta(args.query)

    ref = ref_seqs[args.ref_index]
    query = query_seqs[args.query_index]

    mutation_report = compare(ref, query)
    report.print_mutation_report(mutation_report)
    print()


def main():
    parser = argparse.ArgumentParser(
        prog="dna_analyzer",
        description="DNA Sequence Analyzer — bioinformatics toolkit",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # --- analyze ---
    p_analyze = subparsers.add_parser("analyze", help="Full analysis: GC content + ORFs")
    p_analyze.add_argument("fasta", help="Path to input FASTA file")
    p_analyze.add_argument("--min-orf", type=int, default=100, metavar="BP",
                           help="Minimum ORF length in bp (default: 100)")
    p_analyze.add_argument("--windows", action="store_true",
                           help="Show sliding window GC content")
    p_analyze.add_argument("--window-size", type=int, default=60, metavar="N",
                           help="Sliding window size (default: 60)")
    p_analyze.add_argument("--step", type=int, default=10, metavar="N",
                           help="Sliding window step size (default: 10)")
    p_analyze.set_defaults(func=cmd_analyze)

    # --- orfs ---
    p_orfs = subparsers.add_parser("orfs", help="Find Open Reading Frames")
    p_orfs.add_argument("fasta", help="Path to input FASTA file")
    p_orfs.add_argument("--min-orf", type=int, default=100, metavar="BP",
                        help="Minimum ORF length in bp (default: 100)")
    p_orfs.set_defaults(func=cmd_orfs)

    # --- compare ---
    p_compare = subparsers.add_parser("compare", help="Compare two sequences and detect mutations")
    p_compare.add_argument("ref", help="Reference FASTA file")
    p_compare.add_argument("query", help="Query FASTA file")
    p_compare.add_argument("--ref-index", type=int, default=0, metavar="N",
                           help="Which sequence to use from ref file (default: 0)")
    p_compare.add_argument("--query-index", type=int, default=0, metavar="N",
                           help="Which sequence to use from query file (default: 0)")
    p_compare.set_defaults(func=cmd_compare)

    args = parser.parse_args()
    try:
        args.func(args)
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
