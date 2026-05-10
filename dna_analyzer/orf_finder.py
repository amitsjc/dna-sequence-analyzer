"""ORF finder — detects Open Reading Frames across all 6 reading frames."""

from dataclasses import dataclass
from dna_analyzer.parser import Sequence

# Standard genetic code (codon → amino acid, single-letter)
CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

START_CODON = "ATG"
STOP_CODONS = {"TAA", "TAG", "TGA"}


@dataclass
class ORF:
    """Represents a single Open Reading Frame."""
    start: int         # 0-based start position on the original sequence
    end: int           # 0-based end position (exclusive)
    frame: int         # reading frame: +1, +2, +3, -1, -2, -3
    strand: str        # "+" or "-"
    nucleotides: str   # the raw DNA of this ORF
    protein: str       # translated amino acid sequence (without stop *)

    def __len__(self):
        return self.end - self.start

    def __repr__(self):
        return (
            f"ORF(frame={self.frame:+d}, start={self.start}, end={self.end}, "
            f"length={len(self)}bp, protein_len={len(self.protein)}aa)"
        )


def translate(sequence: str) -> str:
    """
    Translate a DNA sequence into a protein sequence.
    Reads in triplets from position 0. Stops at first stop codon.
    Returns single-letter amino acid string without the stop codon symbol.
    """
    seq = sequence.upper()
    protein = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        aa = CODON_TABLE.get(codon, "?")
        if aa == "*":
            break
        protein.append(aa)
    return "".join(protein)


def _find_orfs_in_strand(sequence: str, strand: str, offset_in_original: int = 0) -> list[ORF]:
    """Find all ORFs in all 3 reading frames of a single strand."""
    orfs = []
    seq = sequence.upper()
    length = len(seq)

    for frame_offset in range(3):
        frame_num = frame_offset + 1
        i = frame_offset
        orf_start = None

        while i <= length - 3:
            codon = seq[i:i+3]
            if codon == START_CODON and orf_start is None:
                orf_start = i
            elif codon in STOP_CODONS and orf_start is not None:
                orf_end = i + 3
                orf_seq = seq[orf_start:orf_end]
                protein = translate(orf_seq)
                if strand == "+":
                    abs_start = offset_in_original + orf_start
                    abs_end = offset_in_original + orf_end
                else:
                    # Reverse complement coords: map back to forward strand
                    abs_end = offset_in_original - orf_start
                    abs_start = offset_in_original - orf_end
                orfs.append(ORF(
                    start=abs_start,
                    end=abs_end,
                    frame=frame_num if strand == "+" else -frame_num,
                    strand=strand,
                    nucleotides=orf_seq,
                    protein=protein,
                ))
                orf_start = None
            i += 3

    return orfs


def find_orfs(seq: Sequence, min_length: int = 100) -> list[ORF]:
    """
    Find all ORFs in all 6 reading frames of a sequence.

    Searches 3 frames on the forward (+) strand and 3 frames on the
    reverse complement (-) strand. Returns ORFs sorted by length
    (longest first), filtered to those >= min_length nucleotides.
    """
    forward = seq.sequence
    reverse = seq.reverse_complement()

    orfs = []
    orfs += _find_orfs_in_strand(forward, strand="+", offset_in_original=0)
    orfs += _find_orfs_in_strand(reverse, strand="-", offset_in_original=len(forward))

    orfs = [o for o in orfs if len(o) >= min_length]
    orfs.sort(key=len, reverse=True)
    return orfs


def longest_orf(seq: Sequence, min_length: int = 100) -> ORF | None:
    """Return the single longest ORF, or None if none found."""
    orfs = find_orfs(seq, min_length=min_length)
    return orfs[0] if orfs else None
