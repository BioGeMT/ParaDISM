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


def _find_fastq_files(input_path: Path) -> list[Path]:
    return sorted(list(input_path.glob("*.fq")) + list(input_path.glob("*.fastq")))


def _find_reference_files(input_path: Path) -> list[Path]:
    reference_files: list[Path] = []
    for pattern in ("*.fa", "*.fasta", "*.fas", "*.fna"):
        reference_files.extend(input_path.glob(pattern))
    return sorted(reference_files)


def _is_generated_output_path(path: Path, root_path: Path) -> bool:
    try:
        relative_parts = path.relative_to(root_path).parts
    except ValueError:
        return False
    return (
        "final_outputs" in relative_parts
        or any(part.startswith("iteration_") for part in relative_parts)
    )


def _candidate_reference_files(search_root: Path, reads_path: Path) -> list[Path]:
    candidates: list[Path] = []
    seen: set[Path] = set()

    for path in _find_reference_files(reads_path):
        resolved_path = path.resolve()
        candidates.append(resolved_path)
        seen.add(resolved_path)

    recursive_paths: list[Path] = []
    for pattern in ("*.fa", "*.fasta", "*.fas", "*.fna"):
        recursive_paths.extend(search_root.rglob(pattern))
    for path in sorted(recursive_paths):
        resolved_path = path.resolve()
        if resolved_path in seen:
            continue
        if _is_generated_output_path(path.parent, search_root):
            continue
        candidates.append(resolved_path)
        seen.add(resolved_path)

    return candidates


def _candidate_input_directories(input_path: Path) -> list[Path]:
    return sorted(
        {
            path.parent
            for pattern in ("*.fq", "*.fastq")
            for path in input_path.rglob(pattern)
            if path.parent != input_path and not _is_generated_output_path(path.parent, input_path)
        }
    )


def _resolve_user_path(raw_path: str, base_path: Path) -> Path:
    selected_path = Path(raw_path).expanduser()
    if not selected_path.is_absolute():
        selected_path = base_path / selected_path
    return selected_path.resolve()


def _select_input_directory(input_path: Path) -> Path:
    while True:
        if _find_fastq_files(input_path):
            return input_path

        console.print(f"\n[yellow]No FASTQ files found in {input_path.resolve()}[/yellow]")
        candidates = _candidate_input_directories(input_path)
        if len(candidates) == 1:
            candidate = candidates[0]
            selected_path = candidate.resolve()
            fastq_count = len(_find_fastq_files(selected_path))
            ref_count = len(find_references(str(selected_path)))
            console.print(
                f"[green]✓[/green] Auto-selected reads directory: "
                f"[cyan]{candidate.relative_to(input_path)}[/cyan] "
                f"[dim]({fastq_count} FASTQ, {ref_count} FASTA)[/dim]"
            )
            return selected_path
        if candidates:
            console.print("[cyan]Reads directories with FASTQ files:[/cyan]")
            for index, candidate in enumerate(candidates, 1):
                fastq_count = len(_find_fastq_files(candidate))
                ref_count = len(find_references(str(candidate)))
                console.print(
                    f"  [green]{index}[/green]. {candidate.relative_to(input_path)} "
                    f"[dim]({fastq_count} FASTQ, {ref_count} FASTA)[/dim]"
                )
            console.print()
            prompt = f"[green]Select reads directory [1-{len(candidates)}], enter a reads directory path, or 'q' to quit:[/green] "
        else:
            prompt = "[green]Enter a reads directory path, or 'q' to quit:[/green] "

        choice = console.input(prompt).strip()
        if choice.lower() in {"q", "quit", "exit"}:
            sys.exit(1)

        try:
            index = int(choice) - 1
        except ValueError:
            selected_path = _resolve_user_path(choice, input_path)
        else:
            if not (0 <= index < len(candidates)):
                console.print(f"[red]Invalid selection. Please enter 1-{len(candidates)} or a reads directory path[/red]")
                continue
            selected_path = candidates[index].resolve()

        if not selected_path.exists():
            console.print(f"[red]Reads directory not found: {selected_path}[/red]")
            continue
        if not selected_path.is_dir():
            console.print(f"[red]Not a directory: {selected_path}[/red]")
            continue
        input_path = selected_path


def _select_reference_path(input_path: Path) -> Path:
    console.print(f"\n[yellow]No reference FASTA files found in {input_path.resolve()}[/yellow]")
    while True:
        choice = console.input("[green]Enter reference FASTA path, or 'q' to quit:[/green] ").strip()
        if choice.lower() in {"q", "quit", "exit"}:
            sys.exit(1)

        selected_path = _resolve_user_path(choice, input_path)
        if not selected_path.exists():
            console.print(f"[red]Reference file not found: {selected_path}[/red]")
            continue
        if not selected_path.is_file():
            console.print(f"[red]Reference is not a file: {selected_path}[/red]")
            continue
        return selected_path


def interactive_mode(input_dir: str = ".", output_dir: str = "./output", reference: str | None = None):
    input_path = Path(input_dir)
    if not input_path.exists():
        console.print(f"[red]✗ Input directory not found: {input_dir}[/red]")
        sys.exit(1)
    if not input_path.is_dir():
        console.print(f"[red]✗ Not a directory: {input_dir}[/red]")
        sys.exit(1)

    provided_reference = None
    if reference:
        provided_reference = Path(reference)
        if not provided_reference.exists():
            console.print(f"[red]✗ Reference file not found: {reference}[/red]")
            sys.exit(1)
        if not provided_reference.is_file():
            console.print(f"[red]✗ Reference is not a file: {reference}[/red]")
            sys.exit(1)
        provided_reference = provided_reference.resolve()
    
    search_root = input_path.resolve()
    input_dir_resolved = str(search_root)

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

    input_path = _select_input_directory(input_path)
    if input_dir_resolved != ".":
        console.print(f"[cyan]Scanning directory: {input_path.resolve()}[/cyan]\n")
    input_dir_resolved = str(input_path.resolve())

    all_fastq_files = _find_fastq_files(input_path)

    # Only auto-detect pairs if in paired-end mode
    fastq_pairs = []
    if is_paired:
        pairs = find_paired_fastqs(input_dir_resolved)
        for r1_path, r2_path in pairs:
            r1_size = os.path.getsize(r1_path)
            r2_size = os.path.getsize(r2_path)
            fastq_pairs.append((r1_path, r2_path, r1_size, r2_size))

    references = []
    ref_files = _candidate_reference_files(search_root, input_path)
    for ref_path in ref_files:
        if provided_reference and ref_path == provided_reference:
            continue
        references.append((str(ref_path), os.path.getsize(ref_path)))
    if not references and not provided_reference:
        provided_reference = _select_reference_path(input_path)
    if provided_reference:
        references.insert(0, (str(provided_reference), os.path.getsize(provided_reference)))

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

    all_ref_files = [Path(ref_path) for ref_path, _ref_size in references]

    if provided_reference:
        ref_path = str(provided_reference)
        console.print(f"[green]✓[/green] Using provided reference: [cyan]{provided_reference.name}[/cyan]")
        ref_meta = scan_fasta_metadata(ref_path)
        console.print()
    elif len(references) == 1:
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
    threads = 4
    workers = 1
    threshold = None
    iterations = 1
    n_anchors = 1
    min_alternate_count = 5
    add_quality_filters = False
    qual_threshold = 20
    dp_threshold = 10
    af_threshold = 0.05

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
        aligner_table.add_row("3", "Minimap2", "(short-read mode)")

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
                console.print("[green]✓[/green] Selected: [cyan]Minimap2[/cyan]")
                break
            else:
                console.print("[red]Invalid selection. Please enter 1-3[/red]")

        console.print()

        # Ask for threshold (optional)
        threshold = None
        if aligner == "bowtie2":
            default_threshold = "G,40,40"
            threshold_prompt = f"[green]Alignment score threshold (default: {default_threshold} for 150bp, use G,30,30 for 100bp):[/green] "
        else:
            default_threshold = "240"
            threshold_prompt = f"[green]Alignment score threshold (default: {default_threshold} for 150bp, use 160 for 100bp):[/green] "

        while True:
            choice = console.input(threshold_prompt).strip()
            if not choice:
                threshold = None  # Use default
                console.print(f"[green]✓[/green] Threshold: [cyan]{default_threshold}[/cyan] (default)")
                break
            else:
                threshold = choice
                console.print(f"[green]✓[/green] Threshold: [cyan]{threshold}[/cyan]")
                break

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

    max_workers = os.cpu_count() or 1
    while True:
        choice = console.input("[green]Read-assignment worker processes (default: 1):[/green] ").strip()
        if not choice:
            workers = 1
            console.print("[green]✓[/green] Assignment workers: [cyan]1[/cyan] (default)")
            break
        try:
            workers = int(choice)
            if workers < 1:
                raise ValueError
            if workers > max_workers:
                console.print(f"[yellow]⚠ Warning: {workers} workers requested but only {max_workers} CPUs available[/yellow]")
            console.print(f"[green]✓[/green] Assignment workers: [cyan]{workers}[/cyan]")
            break
        except ValueError:
            console.print("[red]Invalid input. Please enter a positive integer.[/red]")

    console.print()

    # Ask for iterations (optional) - after thread/worker selection
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

    # Variant calling options (only relevant if iterations > 1)
    if iterations > 1:
        console.print()
        print_section("Variant Calling Options")

        # Min alternate count
        while True:
            choice = console.input("[green]Minimum alternate allele count for variant calling (default: 5):[/green] ").strip()
            if not choice:
                min_alternate_count = 5
                console.print(f"[green]✓[/green] Min alternate count: [cyan]5[/cyan] (default)")
                break
            try:
                min_alternate_count = int(choice)
                if min_alternate_count < 1:
                    raise ValueError
                console.print(f"[green]✓[/green] Min alternate count: [cyan]{min_alternate_count}[/cyan]")
                break
            except ValueError:
                console.print("[red]Invalid input. Please enter a positive integer.[/red]")

        # Quality filters
        console.print()
        while True:
            choice = console.input("[green]Apply quality filters? (y/n, default: n):[/green] ").strip().lower()
            if choice in ['', 'n', 'no']:
                add_quality_filters = False
                console.print(f"[green]✓[/green] Quality filters: [cyan]disabled[/cyan]")
                break
            elif choice in ['y', 'yes']:
                add_quality_filters = True
                console.print(f"[green]✓[/green] Quality filters: [cyan]enabled[/cyan]")

                # QUAL threshold
                console.print()
                while True:
                    choice = console.input("[green]  Minimum QUAL score (default: 20):[/green] ").strip()
                    if not choice:
                        qual_threshold = 20
                        console.print(f"  [green]✓[/green] QUAL threshold: [cyan]20[/cyan] (default)")
                        break
                    try:
                        qual_threshold = int(choice)
                        if qual_threshold < 0:
                            raise ValueError
                        console.print(f"  [green]✓[/green] QUAL threshold: [cyan]{qual_threshold}[/cyan]")
                        break
                    except ValueError:
                        console.print("  [red]Invalid input. Please enter a non-negative integer.[/red]")

                # DP threshold
                while True:
                    choice = console.input("[green]  Minimum depth DP (default: 10):[/green] ").strip()
                    if not choice:
                        dp_threshold = 10
                        console.print(f"  [green]✓[/green] DP threshold: [cyan]10[/cyan] (default)")
                        break
                    try:
                        dp_threshold = int(choice)
                        if dp_threshold < 0:
                            raise ValueError
                        console.print(f"  [green]✓[/green] DP threshold: [cyan]{dp_threshold}[/cyan]")
                        break
                    except ValueError:
                        console.print("  [red]Invalid input. Please enter a non-negative integer.[/red]")

                # AF threshold
                while True:
                    choice = console.input("[green]  Minimum allele frequency AF (default: 0.05):[/green] ").strip()
                    if not choice:
                        af_threshold = 0.05
                        console.print(f"  [green]✓[/green] AF threshold: [cyan]0.05[/cyan] (default)")
                        break
                    try:
                        af_threshold = float(choice)
                        if af_threshold < 0 or af_threshold > 1:
                            raise ValueError
                        console.print(f"  [green]✓[/green] AF threshold: [cyan]{af_threshold}[/cyan]")
                        break
                    except ValueError:
                        console.print("  [red]Invalid input. Please enter a number between 0 and 1.[/red]")
                break
            else:
                console.print("[red]Please answer y or n[/red]")

    console.print()
    while True:
        choice = console.input("[green]Minimum C1 anchor positions for assignment (default: 1):[/green] ").strip()
        if not choice:
            n_anchors = 1
            console.print("[green]✓[/green] C1 anchors: [cyan]1[/cyan] (default)")
            break
        try:
            n_anchors = int(choice)
            if n_anchors < 1:
                raise ValueError
            console.print(f"[green]✓[/green] C1 anchors: [cyan]{n_anchors}[/cyan]")
            break
        except ValueError:
            console.print("[red]Invalid input. Please enter a positive integer.[/red]")

    console.print()

    choice = console.input(f"[green]Output directory (default: {output_dir}):[/green] ").strip()
    if choice:
        output_dir = str(Path(choice).expanduser())
    console.print(f"[green]✓[/green] Output directory: [cyan]{output_dir}[/cyan]")
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
        workers=workers,
        sam_file=Path(sam_path).name if sam_path else None,
        output_dir=output_dir,
        iterations=iterations,
        n_anchors=n_anchors,
        min_alternate_count=min_alternate_count,
        add_quality_filters=add_quality_filters,
        qual_threshold=qual_threshold,
        dp_threshold=dp_threshold,
        af_threshold=af_threshold,
    )

    console.print()
    console.input("[green]Press ENTER to start pipeline or Ctrl+C to cancel...[/green]")

    # Create output directory if it doesn't exist
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    executor = SimpleParaDISMExecutor(
        output_dir=output_dir,
        min_alternate_count=min_alternate_count,
        add_quality_filters=add_quality_filters,
        qual_threshold=qual_threshold,
        dp_threshold=dp_threshold,
        af_threshold=af_threshold,
    )
    executor.run_pipeline(
        r1=r1_path,
        r2=r2_path,
        ref=ref_path,
        aligner=aligner,
        threads=threads,
        sam=sam_path,
        show_header=True,
        iterations=iterations,
        threshold=threshold,
        n_anchors=n_anchors,
        workers=workers,
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
