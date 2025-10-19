#!/usr/bin/env python3
"""Interactive workflow for the homologous-region mapper."""

import sys
from pathlib import Path

from pipeline.executor import PipelineExecutor
from ui.ui_components import (
    display_file_pairs,
    display_validation_results,
    display_pipeline_config,
    print_header,
    print_section,
    console,
)
from utils.file_scanner import (
    scan_fastq_metadata,
    scan_fasta_metadata,
    scan_sam_metadata,
    find_paired_fastqs,
    find_references,
    find_sam_files,
)
from utils.validators import (
    validate_fastq_pair,
    validate_fasta,
    validate_sam,
)


def interactive_mode():
    """Run mapper in interactive mode with rich UI."""

    print_header()

    print_section("Input Files")

    pairs = find_paired_fastqs(".")
    fastq_pairs = []
    for r1_path, r2_path in pairs:
        r1_meta = scan_fastq_metadata(r1_path)
        r2_meta = scan_fastq_metadata(r2_path)
        fastq_pairs.append((r1_path, r2_path, r1_meta, r2_meta))

    ref_files = find_references(".")
    references = []
    for ref_path in ref_files:
        ref_meta = scan_fasta_metadata(ref_path)
        references.append((ref_path, ref_meta))

    sam_file_paths = find_sam_files(".")
    sam_files = []
    for sam_path in sam_file_paths:
        sam_meta = scan_sam_metadata(sam_path)
        sam_files.append((sam_path, sam_meta))

    display_file_pairs(fastq_pairs, references, sam_files)

    all_fastq_files = sorted(list(Path(".").glob("*.fq")) + list(Path(".").glob("*.fastq")))

    if not fastq_pairs and len(all_fastq_files) < 2:
        console.print("\n[red]✗ No FASTQ files found![/red]")
        console.print("Please ensure you have FASTQ files in the current directory.\n")
        sys.exit(1)

    if not references:
        console.print("\n[red]✗ No reference FASTA files found![/red]")
        console.print("Please ensure you have a reference FASTA file in the current directory.\n")
        sys.exit(1)

    console.print()

    if not fastq_pairs:
        console.print("[yellow]⚠ No auto-paired FASTQ files detected.[/yellow]")
        console.print("[yellow]  Please manually select R1 and R2 files.[/yellow]\n")

        console.print("[cyan]Available FASTQ files:[/cyan]")
        for idx, fq_file in enumerate(all_fastq_files, 1):
            fq_meta = scan_fastq_metadata(str(fq_file))
            console.print(f"  [green]{idx}[/green]. {fq_file.name} [dim]({fq_meta['size_human']}, {fq_meta['read_count']:,} reads)[/dim]")

        console.print()

        while True:
            choice = console.input(f"[green]Select R1 file [1-{len(all_fastq_files)}]:[/green] ").strip()
            try:
                idx = int(choice) - 1
                if 0 <= idx < len(all_fastq_files):
                    r1_path = str(all_fastq_files[idx])
                    r1_meta = scan_fastq_metadata(r1_path)
                    console.print(f"[green]✓[/green] R1: [cyan]{Path(r1_path).name}[/cyan]")
                    break
                else:
                    console.print(f"[red]Invalid selection. Please enter 1-{len(all_fastq_files)}[/red]")
            except ValueError:
                console.print(f"[red]Invalid selection. Please enter 1-{len(all_fastq_files)}[/red]")

        while True:
            choice = console.input(f"[green]Select R2 file [1-{len(all_fastq_files)}]:[/green] ").strip()
            try:
                idx = int(choice) - 1
                if 0 <= idx < len(all_fastq_files):
                    r2_path = str(all_fastq_files[idx])
                    if r2_path == r1_path:
                        console.print("[red]R2 cannot be the same as R1. Please select a different file.[/red]")
                        continue
                    r2_meta = scan_fastq_metadata(r2_path)
                    console.print(f"[green]✓[/green] R2: [cyan]{Path(r2_path).name}[/cyan]")
                    break
                else:
                    console.print(f"[red]Invalid selection. Please enter 1-{len(all_fastq_files)}[/red]")
            except ValueError:
                console.print(f"[red]Invalid selection. Please enter 1-{len(all_fastq_files)}[/red]")
    else:
        while True:
            choice = console.input(f"[green]Select FASTQ pair [1-{len(fastq_pairs)}] or 'm' for manual selection:[/green] ").strip().lower()
            if choice == 'm':
                console.print()
                console.print("[cyan]Available FASTQ files:[/cyan]")
                for idx, fq_file in enumerate(all_fastq_files, 1):
                    fq_meta = scan_fastq_metadata(str(fq_file))
                    console.print(f"  [green]{idx}[/green]. {fq_file.name} [dim]({fq_meta['size_human']}, {fq_meta['read_count']:,} reads)[/dim]")
                console.print()

                while True:
                    r1_choice = console.input(f"[green]Select R1 file [1-{len(all_fastq_files)}]:[/green] ").strip()
                    try:
                        idx = int(r1_choice) - 1
                        if 0 <= idx < len(all_fastq_files):
                            r1_path = str(all_fastq_files[idx])
                            r1_meta = scan_fastq_metadata(r1_path)
                            console.print(f"[green]✓[/green] R1: [cyan]{Path(r1_path).name}[/cyan]")
                            break
                        else:
                            console.print(f"[red]Invalid selection. Please enter 1-{len(all_fastq_files)}[/red]")
                    except ValueError:
                        console.print(f"[red]Invalid selection. Please enter 1-{len(all_fastq_files)}[/red]")

                while True:
                    r2_choice = console.input(f"[green]Select R2 file [1-{len(all_fastq_files)}]:[/green] ").strip()
                    try:
                        idx = int(r2_choice) - 1
                        if 0 <= idx < len(all_fastq_files):
                            r2_path = str(all_fastq_files[idx])
                            if r2_path == r1_path:
                                console.print("[red]R2 cannot be the same as R1. Please select a different file.[/red]")
                                continue
                            r2_meta = scan_fastq_metadata(r2_path)
                            console.print(f"[green]✓[/green] R2: [cyan]{Path(r2_path).name}[/cyan]")
                            break
                        else:
                            console.print(f"[red]Invalid selection. Please enter 1-{len(all_fastq_files)}[/red]")
                    except ValueError:
                        console.print(f"[red]Invalid selection. Please enter 1-{len(all_fastq_files)}[/red]")
                break
            else:
                choice = choice if choice else "1"
                try:
                    idx = int(choice) - 1
                    if 0 <= idx < len(fastq_pairs):
                        r1_path, r2_path, r1_meta, r2_meta = fastq_pairs[idx]
                        console.print(f"[green]✓[/green] Selected: [cyan]{Path(r1_path).name}[/cyan] / [cyan]{Path(r2_path).name}[/cyan]")
                        break
                    else:
                        console.print(f"[red]Invalid selection. Please enter 1-{len(fastq_pairs)} or 'm' for manual[/red]")
                except ValueError:
                    console.print(f"[red]Invalid selection. Please enter 1-{len(fastq_pairs)} or 'm' for manual[/red]")

    console.print()

    all_ref_files = sorted(
        list(Path(".").glob("*.fa")) +
        list(Path(".").glob("*.fasta")) +
        list(Path(".").glob("*.fas")) +
        list(Path(".").glob("*.fna"))
    )

    if len(references) == 1:
        ref_path, ref_meta = references[0]
        console.print(f"[green]✓[/green] Auto-selected: [cyan]{Path(ref_path).name}[/cyan]")
        console.print()
    else:
        while True:
            choice = console.input(
                f"[green]Select reference [1-{len(references)}], 'm' to see all FASTA files, or ENTER for first:[/green] "
            ).strip().lower()

            if choice == 'm':
                console.print()
                console.print("[cyan]Available FASTA files:[/cyan]")
                for idx, ref_file in enumerate(all_ref_files, 1):
                    ref_meta_temp = scan_fasta_metadata(str(ref_file))
                    console.print(f"  [green]{idx}[/green]. {ref_file.name} [dim]({ref_meta_temp['num_sequences']} sequences, {ref_meta_temp['size_human']})[/dim]")
                console.print()

                while True:
                    ref_choice = console.input(f"[green]Select reference [1-{len(all_ref_files)}]:[/green] ").strip()
                    try:
                        idx = int(ref_choice) - 1
                        if 0 <= idx < len(all_ref_files):
                            ref_path = str(all_ref_files[idx])
                            ref_meta = scan_fasta_metadata(ref_path)
                            console.print(f"[green]✓[/green] Selected: [cyan]{Path(ref_path).name}[/cyan]")
                            break
                        else:
                            console.print(f"[red]Invalid selection. Please enter 1-{len(all_ref_files)}[/red]")
                    except ValueError:
                        console.print(f"[red]Invalid selection. Please enter 1-{len(all_ref_files)}[/red]")
                break
            else:
                choice = choice if choice else "1"
                try:
                    idx = int(choice) - 1
                    if 0 <= idx < len(references):
                        ref_path, ref_meta = references[idx]
                        console.print(f"[green]✓[/green] Selected: [cyan]{Path(ref_path).name}[/cyan]")
                        break
                    else:
                        console.print(f"[red]Invalid selection. Please enter 1-{len(references)} or 'm' to see all files[/red]")
                except ValueError:
                    console.print(f"[red]Invalid selection. Please enter 1-{len(references)} or 'm' to see all files[/red]")

    sam_path = None
    if sam_files:
        console.print()
        while True:
            choice = console.input("[green]Use existing SAM alignment? (y/n):[/green] ").strip().lower()
            if choice in ['y', 'yes']:
                if len(sam_files) == 1:
                    sam_path, sam_meta = sam_files[0]
                    console.print(f"[green]✓[/green] Auto-selected: [cyan]{Path(sam_path).name}[/cyan]")
                else:
                    while True:
                        sam_choice = console.input(f"[green]Select SAM file [1-{len(sam_files)}]:[/green] ").strip()
                        try:
                            idx = int(sam_choice) - 1
                            if 0 <= idx < len(sam_files):
                                sam_path, sam_meta = sam_files[idx]
                                console.print(f"[green]✓[/green] Selected: [cyan]{Path(sam_path).name}[/cyan]")
                                break
                            else:
                                console.print(f"[red]Invalid selection. Please enter 1-{len(sam_files)}[/red]")
                        except ValueError:
                            console.print(f"[red]Invalid selection. Please enter 1-{len(sam_files)}[/red]")
                break
            elif choice in ['', 'n', 'no']:
                break
            else:
                console.print("[yellow]Please answer y or n[/yellow]")

    aligner = "bwa-mem2"
    minimap2_profile = "short"
    threads = 4

    if not sam_path:
        from rich import box as rbox
        from rich.panel import Panel
        from rich.table import Table

        console.print()
        print_section("Aligner Selection")

        aligner_table = Table(show_header=False, box=None, padding=(0, 1), expand=False)
        aligner_table.add_column("Number", style="dim", width=4)
        aligner_table.add_column("Aligner", style="bright_white")
        aligner_table.add_column("Description", style="cyan")
        aligner_table.add_row("1", "BWA-MEM2", "(default, fast for short reads)")
        aligner_table.add_row("2", "Bowtie2", "(alternative for short reads)")
        aligner_table.add_row("3", "Minimap2", "(versatile, supports long reads)")

        aligner_panel = Panel(
            aligner_table,
            title="[bold green]Select Aligner[/bold green]",
            title_align="left",
            border_style="green",
            box=rbox.ROUNDED,
            expand=False,
            padding=(1, 2)
        )
        console.print(aligner_panel)
        console.print()

        while True:
            choice = console.input("[green]Select [1-3] or ENTER for BWA-MEM2:[/green] ").strip()
            choice = choice if choice else "1"

            if choice == "1":
                aligner = "bwa-mem2"
                console.print("[green]✓[/green] Selected: [cyan]BWA-MEM2[/cyan]")
                break
            elif choice == "2":
                aligner = "bowtie2"
                console.print("[green]✓[/green] Selected: [cyan]Bowtie2[/cyan]")
                break
            elif choice == "3":
                aligner = "minimap2"

                console.print()
                console.print("[cyan]Select Minimap2 Profile:[/cyan]")
                console.print("  [green]1[/green] - short (Illumina paired-end)")
                console.print("  [green]2[/green] - pacbio (PacBio HiFi/CLR)")
                console.print("  [green]3[/green] - nanopore (Oxford Nanopore)")
                console.print()

                while True:
                    profile_choice = console.input("[green]Select [1-3] or ENTER for short:[/green] ").strip()
                    profile_choice = profile_choice if profile_choice else "1"

                    if profile_choice == "1":
                        minimap2_profile = "short"
                        break
                    elif profile_choice == "2":
                        minimap2_profile = "pacbio"
                        break
                    elif profile_choice == "3":
                        minimap2_profile = "nanopore"
                        break
                    else:
                        console.print("[red]Invalid selection. Please enter 1-3[/red]")

                console.print(f"[green]✓[/green] Selected: [cyan]Minimap2 ({minimap2_profile})[/cyan]")
                break
            else:
                console.print("[red]Invalid selection. Please enter 1-3[/red]")

        console.print()

        while True:
            choice = console.input("[green]Number of threads (ENTER for 4):[/green] ").strip()
            if not choice:
                break
            try:
                threads = int(choice)
                if threads < 1:
                    raise ValueError
                console.print(f"[green]✓[/green] Threads: [cyan]{threads}[/cyan]")
                break
            except ValueError:
                console.print("[red]Invalid input. Please enter a positive integer.[/red]")

    console.print()

    validations = []

    if fastq_pairs:
        validations.extend(validate_fastq_pair(r1_path, r2_path, r1_meta, r2_meta))
    if ref_meta:
        validations.extend(validate_fasta(ref_path, ref_meta))
    if sam_path:
        validations.extend(validate_sam(sam_path, sam_meta))

    is_safe = display_validation_results(validations)

    if not is_safe:
        console.print("[red]✗ Blocking errors detected. Please fix the issues above and try again.[/red]")
        sys.exit(1)

    console.print()
    print_section("Ready to Run")

    display_pipeline_config(
        r1_file=Path(r1_path).name,
        r2_file=Path(r2_path).name,
        ref_file=Path(ref_path).name,
        read_count=r1_meta['read_count'],
        ref_sequences=ref_meta['num_sequences'],
        ref_size_kb=ref_meta['total_length'] / 1000,
        aligner=aligner,
        threads=threads,
        sam_file=Path(sam_path).name if sam_path else None,
        minimap2_profile=minimap2_profile if aligner == "minimap2" else None
    )

    console.print()
    console.input("[green]Press ENTER to start pipeline or Ctrl+C to cancel...[/green]")

    executor = PipelineExecutor()
    executor.run_pipeline(
        r1=r1_path,
        r2=r2_path,
        ref=ref_path,
        aligner=aligner,
        threads=threads,
        sam=sam_path,
        minimap2_profile=minimap2_profile
    )

    console.print("[green]═══════════════════════════════════════════════════════[/green]")
    console.print("[green]Pipeline Complete![/green]")
    console.print("[green]═══════════════════════════════════════════════════════[/green]")
    console.print()
    console.print("[cyan]Outputs:[/cyan]")
    console.print("  • Unique mappings: [green]./output/unique_mappings.tsv[/green]")
    console.print("  • Gene-specific FASTQs: [green]./output/fastq/[/green]")
    console.print("  • Gene-specific BAMs: [green]./output/bam/[/green]")
    console.print()
