"""Formatted terminal output for analysis results."""

from colorama import Fore, Style, init

init(autoreset=True)

DIVIDER = "─" * 60


def _header(title: str) -> str:
    return f"\n{Fore.CYAN}{Style.BRIGHT}{DIVIDER}\n  {title}\n{DIVIDER}{Style.RESET_ALL}"


def _label(key: str, value) -> str:
    return f"  {Fore.YELLOW}{key:<22}{Style.RESET_ALL}{value}"


def print_gc_report(summary: dict) -> None:
    print(_header(f"GC Content — {summary['id']}"))
    print(_label("Length", f"{summary['length']} bp"))
    print(_label("GC Content", f"{summary['gc_percent']}%"))
    print(_label("AT Content", f"{summary['at_percent']}%"))
    print(_label("GC Skew", summary['gc_skew']))
    counts = summary['nucleotide_counts']
    print(_label("Nucleotides", f"A={counts['A']}  T={counts['T']}  G={counts['G']}  C={counts['C']}  N={counts['N']}"))


def print_orf_report(seq_id: str, orfs: list) -> None:
    print(_header(f"Open Reading Frames — {seq_id}"))
    if not orfs:
        print(f"  {Fore.RED}No ORFs found above minimum length threshold.{Style.RESET_ALL}")
        return
    print(_label("ORFs found", len(orfs)))
    for i, orf in enumerate(orfs[:10], 1):
        strand_color = Fore.GREEN if orf.strand == "+" else Fore.MAGENTA
        print(
            f"  {Fore.WHITE}{i:>2}.{Style.RESET_ALL} "
            f"frame {strand_color}{orf.frame:+d}{Style.RESET_ALL}  "
            f"pos {orf.start}–{orf.end}  "
            f"{len(orf)} bp  "
            f"protein: {Fore.WHITE}{orf.protein[:30]}{'...' if len(orf.protein) > 30 else ''}{Style.RESET_ALL}"
        )
    if len(orfs) > 10:
        print(f"  {Fore.YELLOW}... and {len(orfs) - 10} more ORFs{Style.RESET_ALL}")


def print_mutation_report(report) -> None:
    print(_header(f"Mutation Report"))
    print(_label("Reference", f"{report.ref_id} ({report.ref_length} bp)"))
    print(_label("Query", f"{report.query_id} ({report.query_length} bp)"))

    identity_color = Fore.GREEN if report.identity_percent >= 99 else Fore.YELLOW if report.identity_percent >= 95 else Fore.RED
    print(_label("Identity", f"{identity_color}{report.identity_percent:.2f}%{Style.RESET_ALL}"))
    print(_label("SNPs", len(report.snps)))
    print(_label("Indels", len(report.indels)))
    print(_label("Transitions", report.transitions))
    print(_label("Transversions", report.transversions))
    ratio = f"{report.ts_tv_ratio}" if report.ts_tv_ratio != float('inf') else "∞ (no transversions)"
    print(_label("Ts/Tv ratio", ratio))

    if report.snps:
        print(f"\n  {Fore.CYAN}SNP details:{Style.RESET_ALL}")
        for snp in report.snps[:20]:
            ts_color = Fore.BLUE if snp.mutation_type == "transition" else Fore.RED
            print(
                f"    pos {snp.position:>5}  "
                f"{Fore.WHITE}{snp.ref_base}→{snp.alt_base}{Style.RESET_ALL}  "
                f"{ts_color}{snp.mutation_type}{Style.RESET_ALL}"
            )
        if len(report.snps) > 20:
            print(f"    {Fore.YELLOW}... and {len(report.snps) - 20} more SNPs{Style.RESET_ALL}")


def print_welcome() -> None:
    print(f"""
{Fore.CYAN}{Style.BRIGHT}
  ██████╗ ███╗   ██╗ █████╗      █████╗ ███╗   ██╗ █████╗ ██╗  ██╗
  ██╔══██╗████╗  ██║██╔══██╗    ██╔══██╗████╗  ██║██╔══██╗██║  ██║
  ██║  ██║██╔██╗ ██║███████║    ███████║██╔██╗ ██║███████║███████║
  ██║  ██║██║╚██╗██║██╔══██║    ██╔══██║██║╚██╗██║██╔══██║██╔══██║
  ██████╔╝██║ ╚████║██║  ██║    ██║  ██║██║ ╚████║██║  ██║██║  ██║
  ╚═════╝ ╚═╝  ╚═══╝╚═╝  ╚═╝    ╚═╝  ╚═╝╚═╝  ╚═══╝╚═╝  ╚═╝╚═╝  ╚═╝
{Style.RESET_ALL}
  {Fore.WHITE}DNA Sequence Analyzer — bioinformatics toolkit{Style.RESET_ALL}
""")
