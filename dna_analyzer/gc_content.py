"""GC content analysis — computes GC%, sliding window, and GC skew."""

from dna_analyzer.parser import Sequence


def gc_content(sequence: str) -> float:
    """
    Return GC content as a percentage (0–100).

    GC pairs have 3 hydrogen bonds vs AT's 2, so GC-rich regions are
    more thermally stable. Used to characterize genomes and identify
    promoter regions and CpG islands.
    """
    seq = sequence.upper()
    if not seq:
        return 0.0
    gc = seq.count("G") + seq.count("C")
    return round((gc / len(seq)) * 100, 2)


def gc_content_windows(sequence: str, window_size: int = 100, step: int = 10) -> list[dict]:
    """
    Compute GC content across a sliding window along the sequence.

    Returns a list of dicts with keys: start, end, gc_percent.
    Useful for spotting GC-rich promoter regions and CpG islands.
    """
    seq = sequence.upper()
    results = []
    for start in range(0, len(seq) - window_size + 1, step):
        window = seq[start : start + window_size]
        results.append({
            "start": start,
            "end": start + window_size,
            "gc_percent": gc_content(window),
        })
    return results


def gc_skew(sequence: str) -> float:
    """
    Compute GC skew: (G - C) / (G + C).

    Values range from -1 to +1. Positive skew = more G than C.
    Sign shifts around the replication origin and terminus, making
    this useful for finding replication origins in bacterial genomes.
    """
    seq = sequence.upper()
    g = seq.count("G")
    c = seq.count("C")
    if g + c == 0:
        return 0.0
    return round((g - c) / (g + c), 4)


def nucleotide_counts(sequence: str) -> dict:
    """Return raw counts of each nucleotide A, T, G, C, and N."""
    seq = sequence.upper()
    return {
        "A": seq.count("A"),
        "T": seq.count("T"),
        "G": seq.count("G"),
        "C": seq.count("C"),
        "N": seq.count("N"),
    }


def summarize(seq: Sequence) -> dict:
    """Return a full GC summary for a Sequence object."""
    counts = nucleotide_counts(seq.sequence)
    return {
        "id": seq.id,
        "length": len(seq),
        "gc_percent": gc_content(seq.sequence),
        "gc_skew": gc_skew(seq.sequence),
        "nucleotide_counts": counts,
        "at_percent": round(100 - gc_content(seq.sequence), 2),
    }
