"""Mutation detector — compares two DNA sequences and classifies differences."""

from dataclasses import dataclass
from dna_analyzer.parser import Sequence

# Transitions: purine↔purine (A↔G) or pyrimidine↔pyrimidine (C↔T)
# Transversions: purine↔pyrimidine (A/G ↔ C/T)
PURINES = {"A", "G"}
PYRIMIDINES = {"C", "T"}


@dataclass
class SNP:
    """A single nucleotide polymorphism at a specific position."""
    position: int     # 0-based
    ref_base: str     # base in the reference sequence
    alt_base: str     # base in the query sequence
    mutation_type: str  # "transition" or "transversion"

    def __repr__(self):
        return (
            f"SNP(pos={self.position}, {self.ref_base}→{self.alt_base}, "
            f"{self.mutation_type})"
        )


@dataclass
class Indel:
    """An insertion or deletion relative to the reference."""
    position: int
    kind: str       # "insertion" or "deletion"
    base: str       # the inserted or deleted base

    def __repr__(self):
        return f"Indel(pos={self.position}, {self.kind}, base={self.base!r})"


@dataclass
class MutationReport:
    """Full comparison result between a reference and query sequence."""
    ref_id: str
    query_id: str
    ref_length: int
    query_length: int
    snps: list[SNP]
    indels: list[Indel]
    identity_percent: float     # % of aligned positions that are identical
    transitions: int
    transversions: int

    @property
    def ts_tv_ratio(self) -> float:
        """Transition/transversion ratio. >2 suggests natural divergence."""
        return round(self.transitions / self.transversions, 2) if self.transversions else float("inf")

    def summary(self) -> str:
        lines = [
            f"Reference : {self.ref_id} ({self.ref_length} bp)",
            f"Query     : {self.query_id} ({self.query_length} bp)",
            f"Identity  : {self.identity_percent:.2f}%",
            f"SNPs      : {len(self.snps)}",
            f"Indels    : {len(self.indels)}",
            f"Transitions    : {self.transitions}",
            f"Transversions  : {self.transversions}",
            f"Ts/Tv ratio    : {self.ts_tv_ratio}",
        ]
        return "\n".join(lines)


def _classify_mutation(ref: str, alt: str) -> str:
    """Return 'transition' or 'transversion' for a SNP."""
    if (ref in PURINES and alt in PURINES) or (ref in PYRIMIDINES and alt in PYRIMIDINES):
        return "transition"
    return "transversion"


def _align_pairwise(ref: str, query: str) -> list[tuple[str | None, str | None]]:
    """
    Needleman-Wunsch global alignment — returns list of (ref_base, query_base) pairs.
    Gap is represented as None. Simple linear gap penalty.
    """
    match_score = 1
    mismatch = -1
    gap = -2

    m, n = len(ref), len(query)
    dp = [[0] * (n + 1) for _ in range(m + 1)]

    for i in range(m + 1):
        dp[i][0] = i * gap
    for j in range(n + 1):
        dp[0][j] = j * gap

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            diag = dp[i-1][j-1] + (match_score if ref[i-1] == query[j-1] else mismatch)
            delete = dp[i-1][j] + gap
            insert = dp[i][j-1] + gap
            dp[i][j] = max(diag, delete, insert)

    # Traceback
    alignment = []
    i, j = m, n
    while i > 0 or j > 0:
        if i > 0 and j > 0:
            diag = dp[i-1][j-1] + (match_score if ref[i-1] == query[j-1] else mismatch)
            if dp[i][j] == diag:
                alignment.append((ref[i-1], query[j-1]))
                i -= 1
                j -= 1
                continue
        if i > 0 and dp[i][j] == dp[i-1][j] + gap:
            alignment.append((ref[i-1], None))
            i -= 1
        else:
            alignment.append((None, query[j-1]))
            j -= 1

    return list(reversed(alignment))


def compare(ref: Sequence, query: Sequence) -> MutationReport:
    """
    Align two sequences and return a full MutationReport.

    Uses Needleman-Wunsch global alignment so insertions and deletions
    are properly handled rather than masked as mismatches.
    """
    alignment = _align_pairwise(ref.sequence, query.sequence)

    snps = []
    indels = []
    identical = 0
    aligned_positions = 0

    for pos, (r, q) in enumerate(alignment):
        if r is None:
            indels.append(Indel(position=pos, kind="insertion", base=q))
        elif q is None:
            indels.append(Indel(position=pos, kind="deletion", base=r))
        else:
            aligned_positions += 1
            if r == q:
                identical += 1
            else:
                mutation_type = _classify_mutation(r, q)
                snps.append(SNP(position=pos, ref_base=r, alt_base=q, mutation_type=mutation_type))

    identity = round((identical / aligned_positions) * 100, 2) if aligned_positions else 0.0
    transitions = sum(1 for s in snps if s.mutation_type == "transition")
    transversions = sum(1 for s in snps if s.mutation_type == "transversion")

    return MutationReport(
        ref_id=ref.id,
        query_id=query.id,
        ref_length=len(ref),
        query_length=len(query),
        snps=snps,
        indels=indels,
        identity_percent=identity,
        transitions=transitions,
        transversions=transversions,
    )
