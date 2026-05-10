# DNA Sequence Analyzer

A Python command-line tool for analyzing DNA sequences. Built for bioinformatics exploration — useful for students and researchers working with genomic data.

## Features

- **FASTA Parser** — load single or multi-sequence `.fasta` files
- **GC Content Calculator** — compute GC% per sequence and genome-wide
- **ORF Finder** — detect Open Reading Frames across all 6 reading frames
- **Mutation Detector** — compare two sequences and highlight SNPs and indels
- **Summary Report** — generate a printable analysis report from the CLI

## Requirements

- Python 3.8+
- See `requirements.txt`

## Installation

```bash
git clone https://github.com/amitsjc/dna-sequence-analyzer.git
cd dna-sequence-analyzer
pip install -r requirements.txt
```

## Usage

```bash
# Analyze a FASTA file
python -m dna_analyzer analyze data/example.fasta

# Find ORFs in a sequence
python -m dna_analyzer orfs data/example.fasta

# Compare two sequences
python -m dna_analyzer compare data/seq1.fasta data/seq2.fasta
```

## Project Structure

```
dna-sequence-analyzer/
├── dna_analyzer/       # Core analysis modules
├── data/               # Sample FASTA files
├── tests/              # Unit tests
├── requirements.txt
└── README.md
```

## Background

DNA sequences are stored as strings of four nucleotides: **A** (adenine), **T** (thymine), **G** (guanine), and **C** (cytosine). This tool implements common bioinformatics algorithms used in genomics research, including GC content analysis (used to identify genome regions and species), open reading frame detection (locating protein-coding genes), and mutation comparison (the basis of SNP analysis in genetics studies).
