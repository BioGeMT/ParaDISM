#!/usr/bin/env python3
"""
Homology Mapper Pipeline - Unified Entry Point

This script handles:
1. Interactive mode with rich UI (no arguments)
2. Argument-driven mode (with --read1/--read2/--reference)
3. Complete pipeline execution (MSA → alignment → mapping → output)
"""

import argparse
import atexit
import sys
from pathlib import Path

SRC_DIR = Path(__file__).parent / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

def _lazy_imports():
    """Import heavy pipeline modules only when needed (not for liftover)."""
    from pipeline.executor import SimpleParaDISMExecutor
    from ui.interactive import interactive_mode
    from ui.ui_components import console
    return SimpleParaDISMExecutor, interactive_mode, console


def _restore_cursor() -> None:
    """Ensure terminal cursor is visible when the process exits."""
    sys.stderr.write("\033[?25h")
    sys.stderr.flush()


atexit.register(_restore_cursor)


def run_with_arguments(args: argparse.Namespace) -> None:
    """Run mapper using explicit command-line arguments (supports both paired-end and single-end)."""
    SimpleParaDISMExecutor, _, console = _lazy_imports()

    if not Path(args.read1).exists():
        console.print(f"[red]✗ R1 file not found: {args.read1}[/red]")
        sys.exit(1)
    if args.read2 and not Path(args.read2).exists():
        console.print(f"[red]✗ R2 file not found: {args.read2}[/red]")
        sys.exit(1)
    if not Path(args.reference).exists():
        console.print(f"[red]✗ Reference file not found: {args.reference}[/red]")
        sys.exit(1)
    if args.sam and not Path(args.sam).exists():
        console.print(f"[red]✗ SAM file not found: {args.sam}[/red]")
        sys.exit(1)

    if args.aligner == "minimap2" and not args.minimap2_profile:
        console.print("[red]✗ --minimap2-profile must be provided when --aligner minimap2[/red]")
        sys.exit(1)

    profile = args.minimap2_profile or "short"

    executor = SimpleParaDISMExecutor(
        output_dir=args.output_dir,
        prefix=args.prefix,
        min_alternate_count=args.min_alternate_count,
        add_quality_filters=args.add_quality_filters,
        qual_threshold=args.qual_threshold,
        dp_threshold=args.dp_threshold,
        af_threshold=args.af_threshold,
    )
    executor.run_pipeline(
        r1=args.read1,
        r2=args.read2,
        ref=args.reference,
        aligner=args.aligner,
        threads=args.threads,
        sam=args.sam,
        minimap2_profile=profile,
        show_header=True,
        iterations=args.iterations,
        threshold=args.threshold,
    )


class CustomHelpFormatter(argparse.RawDescriptionHelpFormatter):
    """Help formatter with wider columns for option descriptions."""

    def __init__(self, prog: str) -> None:
        super().__init__(prog, max_help_position=40)
        self._width = 120


def build_parser() -> argparse.ArgumentParser:
    """Construct the CLI argument parser."""

    parser = argparse.ArgumentParser(
        description="Read-mapping and refinement workflow for highly homologous regions (supports both paired-end and single-end)",
        formatter_class=CustomHelpFormatter,
        epilog="""
Examples:
  # Interactive mode (scans current directory)
  python paradism.py

  # Paired-end mode
  python paradism.py --read1 r1.fq --read2 r2.fq --reference ref.fa

  # Single-end mode
  python paradism.py --read1 reads.fq --reference ref.fa

  # Liftover VCF to chromosomal coordinates
  python paradism.py liftover --vcf variants.vcf --positions positions.txt -o lifted.vcf

  # Liftover BED to chromosomal coordinates
  python paradism.py liftover --bed exons.bed --positions positions.txt -o lifted.bed
        """,
    )

    # Liftover subcommand
    subparsers = parser.add_subparsers(dest="command")
    liftover_parser = subparsers.add_parser(
        "liftover",
        help="Convert VCF/BED from gene-local to chromosomal coordinates",
        formatter_class=CustomHelpFormatter,
    )
    liftover_parser.add_argument("--positions", required=True, help="Gene positions file (e.g., PKD1_b38_pseudogene_positions.txt)")
    liftover_parser.add_argument("--output", "-o", required=True, help="Output file path")
    liftover_group = liftover_parser.add_mutually_exclusive_group(required=True)
    liftover_group.add_argument("--vcf", help="Input VCF file to liftover")
    liftover_group.add_argument("--bed", help="Input BED file to liftover")

    required = parser.add_argument_group(
        "Required arguments (argument-driven mode)",
        "Provide these together to skip the interactive prompts.",
    )
    required.add_argument("--read1", metavar="READ1", help="R1 FASTQ file (or single-end reads)")
    required.add_argument("--read2", metavar="READ2", help="R2 FASTQ file (optional, for paired-end mode)")
    required.add_argument("--reference", metavar="REFERENCE", help="Reference FASTA file")

    optional = parser.add_argument_group("Optional arguments")
    optional.add_argument("--sam", metavar="SAM", help="Existing SAM file (skips new alignment)")
    optional.add_argument(
        "--aligner",
        metavar="ALIGNER",
        default="bowtie2",
        choices=["bwa-mem2", "bowtie2", "minimap2"],
        help=(
            "Read aligner [bwa-mem2|bowtie2|minimap2] (default: Bowtie2 for short reads)"
        ),
    )
    optional.add_argument("--threads", metavar="THREADS", type=int, default=4, help="Threads to use (default: 4)")
    optional.add_argument(
        "--minimap2-profile",
        metavar="PROFILE",
        choices=["short", "pacbio-hifi", "pacbio-clr", "ont-q20", "ont-standard"],
        help="Minimap2 profile (required when --aligner minimap2): "
             "short (short-read), pacbio-hifi (PacBio HiFi/CCS), pacbio-clr (PacBio CLR), "
             "ont-q20 (Nanopore Q20+), ont-standard (Nanopore standard).",
    )
    optional.add_argument(
        "--output-dir",
        metavar="OUTPUT_DIR",
        default="./output",
        help="Output directory (default: ./output)",
    )
    optional.add_argument(
        "--prefix",
        metavar="PREFIX",
        default=None,
        help="Prefix for output files (default: derived from output directory name)",
    )
    optional.add_argument(
        "--iterations",
        metavar="N",
        type=int,
        default=1,
        help="Number of ParaDISM runs (default: 1 = no refinement). "
             "iterations=1 runs ParaDISM once, iterations=2 runs twice (1 refinement iteration), etc. "
             "Each refinement iteration calls variants, updates the reference, and re-runs ParaDISM.",
    )
    optional.add_argument(
        "--input-dir",
        metavar="INPUT_DIR",
        default=".",
        help="Input directory for interactive mode file scanning (default: current directory)",
    )
    optional.add_argument(
        "--threshold",
        metavar="THRESHOLD",
        default=None,
        help="Minimum alignment score threshold. "
             "For bwa-mem2/minimap2: integer score (e.g., 240 for 150bp, 160 for 100bp). "
             "For bowtie2: score function (e.g., 'G,40,40' for 150bp, 'G,30,30' for 100bp). "
             "Default: 240 for bwa-mem2/minimap2, 'G,40,40' for bowtie2.",
    )
    optional.add_argument(
        "--min-alternate-count",
        metavar="N",
        type=int,
        default=5,
        help="Minimum alternate allele count for FreeBayes variant calling (default: 5). "
             "Lower values (e.g., 2-3) increase recall but may decrease precision.",
    )
    optional.add_argument(
        "--add-quality-filters",
        action="store_true",
        help="Apply quality filters during variant calling iterations",
    )
    optional.add_argument(
        "--qual-threshold",
        metavar="N",
        type=int,
        default=20,
        help="Minimum QUAL score for quality filtering (default: 20)",
    )
    optional.add_argument(
        "--dp-threshold",
        metavar="N",
        type=int,
        default=10,
        help="Minimum depth (DP) for quality filtering (default: 10)",
    )
    optional.add_argument(
        "--af-threshold",
        metavar="F",
        type=float,
        default=0.05,
        help="Minimum allele frequency (AF) for quality filtering (default: 0.05)",
    )

    return parser


def main() -> None:
    """Entry point for the mapper CLI."""

    parser = build_parser()
    args = parser.parse_args()

    # Handle liftover subcommand
    if args.command == "liftover":
        from liftover import run_liftover
        run_liftover(args)
        return

    # Check if we have minimum required arguments for CLI mode
    if args.read1 and args.reference:
        run_with_arguments(args)
        return

    _, interactive_mode, console = _lazy_imports()

    # Partial arguments provided
    if any([args.read1, args.reference]):
        console.print(
            "[yellow]⚠ Partial arguments provided. Minimum required: --read1 and --reference[/yellow]"
        )
        console.print("[yellow]  Launching interactive mode instead...[/yellow]\n")

    # Validate input directory if provided
    input_dir = Path(args.input_dir)
    if not input_dir.exists():
        console.print(f"[red]✗ Input directory not found: {args.input_dir}[/red]")
        sys.exit(1)
    if not input_dir.is_dir():
        console.print(f"[red]✗ Not a directory: {args.input_dir}[/red]")
        sys.exit(1)

    interactive_mode(input_dir=str(input_dir.resolve()), output_dir=args.output_dir)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n✗ Pipeline cancelled by user")
        sys.exit(1)
