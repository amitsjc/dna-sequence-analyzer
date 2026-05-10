"""FASTA file parser — reads single and multi-sequence .fasta files."""

from pathlib import Path


class Sequence:
    """Represents a single DNA sequence with its metadata."""

    VALID_BASES = set("ATGCatgcNn")

    def __init__(self, seq_id: str, description: str, sequence: str):
        self.id = seq_id
        self.description = description
        self.sequence = sequence.upper()

    def __len__(self):
        return len(self.sequence)

    def __repr__(self):
        return f"Sequence(id={self.id!r}, length={len(self)})"

    def is_valid(self) -> bool:
        """Return True if the sequence contains only valid DNA bases."""
        return all(base in self.VALID_BASES for base in self.sequence)

    def complement(self) -> str:
        """Return the complementary strand (5'→3')."""
        table = str.maketrans("ATGC", "TACG")
        return self.sequence.translate(table)

    def reverse_complement(self) -> str:
        """Return the reverse complement — the other strand read 5'→3'."""
        return self.complement()[::-1]


def parse_fasta(filepath: str) -> list[Sequence]:
    """
    Parse a FASTA file and return a list of Sequence objects.

    FASTA format: lines starting with '>' are headers; all other lines
    are sequence data that may be split across multiple lines.
    """
    path = Path(filepath)
    if not path.exists():
        raise FileNotFoundError(f"FASTA file not found: {filepath}")
    if not path.suffix.lower() in (".fasta", ".fa", ".fna"):
        raise ValueError(f"Expected a .fasta/.fa/.fna file, got: {path.suffix}")

    sequences = []
    current_id = None
    current_desc = ""
    current_seq_parts = []

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(";"):
                continue
            if line.startswith(">"):
                # Save previous sequence before starting a new one
                if current_id is not None:
                    sequences.append(
                        Sequence(current_id, current_desc, "".join(current_seq_parts))
                    )
                    current_seq_parts = []
                header = line[1:].split(None, 1)
                current_id = header[0]
                current_desc = header[1] if len(header) > 1 else ""
            else:
                current_seq_parts.append(line)

    # Save the final sequence
    if current_id is not None:
        sequences.append(
            Sequence(current_id, current_desc, "".join(current_seq_parts))
        )

    if not sequences:
        raise ValueError(f"No sequences found in {filepath}")

    return sequences
