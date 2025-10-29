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

from pipeline.executor import PipelineExecutor
from ui.interactive import interactive_mode
from ui.ui_components import console


def _restore_cursor() -> None:
    """Ensure terminal cursor is visible when the process exits."""
    sys.stderr.write("\033[?25h")
    sys.stderr.flush()


atexit.register(_restore_cursor)


def run_with_arguments(args: argparse.Namespace) -> None:
    """Run mapper using explicit command-line arguments (supports both paired-end and single-end)."""

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

    # Show mode info
    mode = "paired-end" if args.read2 else "single-end"
    console.print(f"[cyan]Running in {mode} mode[/cyan]")
    if not args.read2:
        console.print("[yellow]⚠ Single-end mode has lower mapping confidence (no mate validation)[/yellow]")

    profile = args.minimap2_profile or "short"

    executor = PipelineExecutor(output_dir=args.output_dir)
    executor.run_pipeline(
        r1=args.read1,
        r2=args.read2,
        ref=args.reference,
        aligner=args.aligner,
        threads=args.threads,
        sam=args.sam,
        minimap2_profile=profile,
        show_header=True,
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
  python mapper.py

  # Interactive mode with custom input directory
  python mapper.py --input-dir /path/to/data

  # Paired-end mode
  python mapper.py --read1 r1.fq --read2 r2.fq --reference ref.fa

  # Single-end mode
  python mapper.py --read1 reads.fq --reference ref.fa

  # Use existing SAM file
  python mapper.py --read1 r1.fq --read2 r2.fq --reference ref.fa --sam mapped.sam

  # Specify aligner and threads
  python mapper.py --read1 r1.fq --read2 r2.fq --reference ref.fa --aligner bwa-mem2 --threads 16

  # PacBio HiFi single-end
  python mapper.py --read1 pacbio_hifi.fq --reference ref.fa --aligner minimap2 --minimap2-profile pacbio-hifi

  # Nanopore Q20+ single-end
  python mapper.py --read1 ont_q20.fq --reference ref.fa --aligner minimap2 --minimap2-profile ont-q20
        """,
    )

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
        default="bwa-mem2",
        choices=["bowtie2", "bwa-mem2", "minimap2"],
        help="Read aligner [bowtie2|bwa-mem2|minimap2] (default: bwa-mem2)",
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
        "--input-dir",
        metavar="INPUT_DIR",
        default=".",
        help="Input directory for interactive mode file scanning (default: current directory)",
    )

    return parser


def main() -> None:
    """Entry point for the mapper CLI."""

    parser = build_parser()
    args = parser.parse_args()

    # Check if we have minimum required arguments for CLI mode
    if args.read1 and args.reference:
        run_with_arguments(args)
        return

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
        console.print("\n[red]✗ Pipeline cancelled by user[/red]")
        sys.exit(1)
