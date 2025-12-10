#!/usr/bin/env python3
"""Interactive workflow for the ParaDISM mapper."""

import sys
import os
from pathlib import Path

from pipeline.executor import SimpleParaDISMExecutor
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
    format_bytes,
)
from utils.validators import (
    validate_fastq_pair,
    validate_fasta,
    validate_sam,
)


def interactive_mode(input_dir: str = ".", output_dir: str = "./output"):
    input_path = Path(input_dir)
    if not input_path.exists():
        console.print(f"[red]✗ Input directory not found: {input_dir}[/red]")
        sys.exit(1)
    if not input_path.is_dir():
        console.print(f"[red]✗ Not a directory: {input_dir}[/red]")
        sys.exit(1)
    
    input_dir_resolved = str(input_path.resolve())

    print_header()

    # Sequencing Mode Selection - FIRST, before file detection
    from rich import box as rbox
    from rich.panel import Panel
    from rich.table import Table

    print_section("Sequencing Mode Selection")

    mode_table = Table(show_header=False, box=None, padding=(0, 1), expand=False)
    mode_table.add_column("Number", style="dim", width=4)
    mode_table.add_column("Mode", style="bright_white")
    mode_table.add_row("1", "Paired-End")
    mode_table.add_row("2", "Single-End")

    mode_panel = Panel(
        mode_table,
        title="[bold green]Select Sequencing Mode[/bold green]",
        title_align="left",
        border_style="green",
        box=rbox.ROUNDED,
        expand=False,
        padding=(1, 2)
    )
    console.print(mode_panel)
    console.print()

    is_paired = True
    while True:
        choice = console.input("[green]Select [1-2]:[/green] ").strip()
        if choice == "1":
            is_paired = True
            console.print("[green]✓[/green] Selected: [cyan]Paired-End Mode[/cyan]")
            break
        elif choice == "2":
            is_paired = False
            console.print("[green]✓[/green] Selected: [cyan]Single-End Mode[/cyan]")
            break
        else:
            console.print("[red]Invalid selection. Please enter 1 or 2[/red]")

    console.print()

    # Now detect files based on mode
    print_section("Input Files")
    
    if input_dir_resolved != ".":
        console.print(f"[cyan]Scanning directory: {input_dir_resolved}[/cyan]\n")

    all_fastq_files = sorted(list(input_path.glob("*.fq")) + list(input_path.glob("*.fastq")))

    # Only auto-detect pairs if in paired-end mode
    fastq_pairs = []
    if is_paired:
        pairs = find_paired_fastqs(input_dir_resolved)
        for r1_path, r2_path in pairs:
            r1_size = os.path.getsize(r1_path)
            r2_size = os.path.getsize(r2_path)
            fastq_pairs.append((r1_path, r2_path, r1_size, r2_size))

    ref_files = find_references(input_dir_resolved)
    references = []
    for ref_path in ref_files:
        ref_size = os.path.getsize(ref_path)
        references.append((ref_path, ref_size))

    sam_file_paths = find_sam_files(input_dir_resolved)
    sam_files = []
    for sam_path in sam_file_paths:
        sam_size = os.path.getsize(sam_path)
        sam_files.append((sam_path, sam_size))

    # Display files appropriately based on mode
    if is_paired:
        display_file_pairs(fastq_pairs, references, sam_files)
    else:
        # For single-end, show individual FASTQ files (no pairing)
        from ui.ui_components import display_single_end_files
        display_single_end_files(all_fastq_files, references, sam_files)

    if not all_fastq_files:
        console.print(f"\n[red]✗ No FASTQ files found in {input_dir_resolved}![/red]")
        console.print(f"Please ensure you have FASTQ files in the specified directory.\n")
        sys.exit(1)

    if not references:
        console.print(f"\n[red]✗ No reference FASTA files found in {input_dir_resolved}![/red]")
        console.print(f"Please ensure you have a reference FASTA file in the specified directory.\n")
        sys.exit(1)

    console.print()

    # FASTQ file selection based on mode
    if not is_paired:
        # === SINGLE-END MODE ===
        # Show all files and let user select one
        while True:
            choice = console.input(f"[green]Select FASTQ file [1-{len(all_fastq_files)}]:[/green] ").strip()
            try:
                idx = int(choice) - 1
                if 0 <= idx < len(all_fastq_files):
                    r1_path = str(all_fastq_files[idx])
                    r1_meta = scan_fastq_metadata(r1_path)
                    console.print(f"[green]✓[/green] Selected: [cyan]{Path(r1_path).name}[/cyan]")
                    r2_path = None
                    r2_meta = None
                    break
                else:
                    console.print(f"[red]Invalid selection. Please enter 1-{len(all_fastq_files)}[/red]")
            except ValueError:
                console.print(f"[red]Invalid selection. Please enter 1-{len(all_fastq_files)}[/red]")

    elif fastq_pairs:
        # === PAIRED-END MODE with auto-detected pairs ===
        while True:
            choice = console.input(f"[green]Select FASTQ pair [1-{len(fastq_pairs)}] or 'm' for manual selection:[/green] ").strip().lower()
            if choice == 'm':
                console.print()
                console.print("[cyan]Available FASTQ files:[/cyan]")
                for idx, fq_file in enumerate(all_fastq_files, 1):
                    fq_size = os.path.getsize(str(fq_file))
                    fq_size_str = format_bytes(fq_size)
                    console.print(f"  [green]{idx}[/green]. {fq_file.name} [dim]({fq_size_str})[/dim]")
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
                try:
                    idx = int(choice) - 1
                    if 0 <= idx < len(fastq_pairs):
                        r1_path, r2_path, r1_size, r2_size = fastq_pairs[idx]
                        console.print(f"[green]✓[/green] Selected: [cyan]{Path(r1_path).name}[/cyan] / [cyan]{Path(r2_path).name}[/cyan]")
                        r1_meta = scan_fastq_metadata(r1_path)
                        r2_meta = scan_fastq_metadata(r2_path)
                        break
                    else:
                        console.print(f"[red]Invalid selection. Please enter 1-{len(fastq_pairs)} or 'm' for manual[/red]")
                except ValueError:
                    console.print(f"[red]Invalid selection. Please enter 1-{len(fastq_pairs)} or 'm' for manual[/red]")

    else:
        # === PAIRED-END MODE without auto-detected pairs ===
        console.print("[yellow]⚠ No auto-paired FASTQ files detected.[/yellow]")
        console.print("[yellow]  Please manually select R1 and R2 files.[/yellow]\n")

        console.print("[cyan]Available FASTQ files:[/cyan]")
        for idx, fq_file in enumerate(all_fastq_files, 1):
            fq_size = os.path.getsize(str(fq_file))
            fq_size_str = format_bytes(fq_size)
            console.print(f"  [green]{idx}[/green]. {fq_file.name} [dim]({fq_size_str})[/dim]")
        console.print()

        # Select R1
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

        # Select R2
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

    console.print()

    all_ref_files = sorted(
        list(input_path.glob("*.fa")) +
        list(input_path.glob("*.fasta")) +
        list(input_path.glob("*.fas")) +
        list(input_path.glob("*.fna"))
    )

    if len(references) == 1:
        ref_path, ref_size = references[0]
        console.print(f"[green]✓[/green] Auto-selected: [cyan]{Path(ref_path).name}[/cyan]")
        ref_meta = scan_fasta_metadata(ref_path)
        console.print()
    else:
        while True:
            choice = console.input(
                f"[green]Select reference [1-{len(references)}] or 'm' to see all FASTA files:[/green] "
            ).strip().lower()

            if choice == 'm':
                console.print()
                console.print("[cyan]Available FASTA files:[/cyan]")
                for idx, ref_file in enumerate(all_ref_files, 1):
                    ref_size = os.path.getsize(str(ref_file))
                    ref_size_str = format_bytes(ref_size)
                    console.print(f"  [green]{idx}[/green]. {ref_file.name} [dim]({ref_size_str})[/dim]")
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
                try:
                    idx = int(choice) - 1
                    if 0 <= idx < len(references):
                        ref_path, ref_size = references[idx]
                        ref_meta = scan_fasta_metadata(ref_path)
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
                    sam_path, sam_size = sam_files[0]
                    console.print(f"[green]✓[/green] Auto-selected: [cyan]{Path(sam_path).name}[/cyan]")
                    sam_meta = scan_sam_metadata(sam_path)
                else:
                    while True:
                        sam_choice = console.input(f"[green]Select SAM file [1-{len(sam_files)}]:[/green] ").strip()
                        try:
                            idx = int(sam_choice) - 1
                            if 0 <= idx < len(sam_files):
                                sam_path, sam_size = sam_files[idx]
                                sam_meta = scan_sam_metadata(sam_path)
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

    aligner = "bowtie2"
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
        aligner_table.add_row("1", "Bowtie2", "(default for short reads)")
        aligner_table.add_row("2", "BWA-MEM2", "(alternative for short reads)")
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
            choice = console.input("[green]Select [1-3]:[/green] ").strip()

            if choice == "1":
                aligner = "bowtie2"
                console.print("[green]✓[/green] Selected: [cyan]Bowtie2[/cyan]")
                break
            elif choice == "2":
                aligner = "bwa-mem2"
                console.print("[green]✓[/green] Selected: [cyan]BWA-MEM2[/cyan]")
                break
            elif choice == "3":
                aligner = "minimap2"

                console.print()
                print_section("Minimap2 Profile Selection")

                profile_table = Table(show_header=False, box=None, padding=(0, 1), expand=False)
                profile_table.add_column("Number", style="dim", width=4)
                profile_table.add_column("Profile", style="bright_white")
                profile_table.add_column("Description", style="cyan")
                profile_table.add_row("1", "Short", "(Short-Read Sequencing)")
                profile_table.add_row("2", "PacBio-HiFi", "(PacBio HiFi/CCS)")
                profile_table.add_row("3", "PacBio-CLR", "(PacBio Continuous)")
                profile_table.add_row("4", "ONT-Q20", "(Nanopore Q20+)")
                profile_table.add_row("5", "ONT-Standard", "(Nanopore Standard)")

                profile_panel = Panel(
                    profile_table,
                    title="[bold green]Select Minimap2 Profile[/bold green]",
                    title_align="left",
                    border_style="green",
                    box=rbox.ROUNDED,
                    expand=False,
                    padding=(1, 2)
                )
                console.print(profile_panel)
                console.print()

                while True:
                    profile_choice = console.input("[green]Select [1-5]:[/green] ").strip()

                    if profile_choice == "1":
                        minimap2_profile = "short"
                        display_name = "Short"
                        break
                    elif profile_choice == "2":
                        minimap2_profile = "pacbio-hifi"
                        display_name = "PacBio-HiFi"
                        break
                    elif profile_choice == "3":
                        minimap2_profile = "pacbio-clr"
                        display_name = "PacBio-CLR"
                        break
                    elif profile_choice == "4":
                        minimap2_profile = "ont-q20"
                        display_name = "ONT-Q20"
                        break
                    elif profile_choice == "5":
                        minimap2_profile = "ont-standard"
                        display_name = "ONT-Standard"
                        break
                    else:
                        console.print("[red]Invalid selection. Please enter 1-5[/red]")

                console.print(f"[green]✓[/green] Selected: [cyan]Minimap2 ({display_name})[/cyan]")
                break
            else:
                console.print("[red]Invalid selection. Please enter 1-3[/red]")

        console.print()

        max_threads = os.cpu_count() or 1
        while True:
            choice = console.input(f"[green]Number of threads ({max_threads} available):[/green] ").strip()
            if not choice:
                console.print("[red]Please enter a number of threads.[/red]")
                continue
            try:
                threads = int(choice)
                if threads < 1:
                    raise ValueError
                if threads > max_threads:
                    console.print(f"[yellow]⚠ Warning: {threads} threads requested but only {max_threads} available[/yellow]")
                console.print(f"[green]✓[/green] Threads: [cyan]{threads}[/cyan]")
                break
            except ValueError:
                console.print("[red]Invalid input. Please enter a positive integer.[/red]")

        console.print()

        # Ask for iterations (optional) - after thread selection
        iterations = 1
        while True:
            choice = console.input("[green]Number of ParaDISM runs (1 = no refinement, 2 = 1 refinement iteration, default: 1):[/green] ").strip()
            if not choice:
                iterations = 1
                break
            try:
                iterations = int(choice)
                if iterations < 1:
                    raise ValueError
                if iterations > 1:
                    console.print(f"[green]✓[/green] Iterations: [cyan]{iterations}[/cyan] (will run ParaDISM {iterations} times)")
                    console.print("[yellow]⚠ Note: Each refinement iteration calls variants, updates reference, and re-runs ParaDISM[/yellow]")
                else:
                    console.print(f"[green]✓[/green] Iterations: [cyan]1[/cyan] (no refinement, single run)")
                break
            except ValueError:
                console.print("[red]Invalid input. Please enter a positive integer (>= 1).[/red]")

    console.print()

    validations = []

    # Validate FASTQ files
    if is_paired and r2_path:
        validations.extend(validate_fastq_pair(r1_path, r2_path, r1_meta, r2_meta))
    else:
        # Single-end validation (just R1)
        from utils.validators import validate_fastq_single
        validations.extend(validate_fastq_single(r1_path, r1_meta))

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
        r2_file=Path(r2_path).name if r2_path else None,
        ref_file=Path(ref_path).name,
        read_count=r1_meta['read_count'],
        ref_sequences=ref_meta['num_sequences'],
        ref_size_kb=ref_meta['total_length'] / 1000,
        aligner=aligner,
        threads=threads,
        sam_file=Path(sam_path).name if sam_path else None,
        minimap2_profile=minimap2_profile if aligner == "minimap2" else None,
        output_dir=output_dir,
        iterations=iterations,
    )

    console.print()
    console.input("[green]Press ENTER to start pipeline or Ctrl+C to cancel...[/green]")

    # Create output directory if it doesn't exist
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    executor = SimpleParaDISMExecutor(output_dir=output_dir)
    executor.run_pipeline(
        r1=r1_path,
        r2=r2_path,
        ref=ref_path,
        aligner=aligner,
        threads=threads,
        sam=sam_path,
        minimap2_profile=minimap2_profile,
        show_header=True,
        iterations=iterations,
    )

    console.print("[green]═══════════════════════════════════════════════════════[/green]")
    console.print("[green]Pipeline Complete![/green]")
    console.print("[green]═══════════════════════════════════════════════════════[/green]")
    console.print()
    console.print("[cyan]Outputs:[/cyan]")
    output_path = Path(output_dir)
    prefix_name = output_path.name if output_path.name != "." else "output"
    console.print(f"  • Gene-specific FASTQs: [green]{output_path}/final_outputs/{prefix_name}_fastq/[/green]")
    console.print(f"  • Gene-specific BAMs: [green]{output_path}/final_outputs/{prefix_name}_bam/[/green]")
    console.print()
