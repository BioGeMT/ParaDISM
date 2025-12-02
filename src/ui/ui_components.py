#!/usr/bin/env python3
"""Rich UI components for file selection, validation, and config display."""

from typing import Dict, List, Optional, Tuple
import os
from pathlib import Path
from rich.console import Console
from rich.table import Table
from rich.panel import Panel
from rich.text import Text
from rich import box


console = Console()


def format_bytes(bytes_size: int) -> str:
    """Return human-readable byte size."""
    for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
        if bytes_size < 1024.0:
            return f"{bytes_size:.1f} {unit}"
        bytes_size /= 1024.0
    return f"{bytes_size:.1f} PB"


def display_single_end_files(
    fastq_files: List,
    references: List[Tuple[str, int]],
    sam_files: List[Tuple[str, int]]
) -> None:
    """Show single-end FASTQs with reference/SAM summaries."""
    from pathlib import Path

    panels = []

    # FASTQ Files Panel (individual files, not pairs)
    if fastq_files:
        table = Table(show_header=False, box=None, padding=(0, 1), expand=False)
        table.add_column("Number", style="dim", width=4)
        table.add_column("File", no_wrap=False)

        for idx, fq_file in enumerate(fastq_files, 1):
            fq_name = fq_file.name if isinstance(fq_file, Path) else Path(fq_file).name
            fq_size = os.path.getsize(str(fq_file))
            fq_info = format_bytes(fq_size)

            table.add_row(
                f"{idx}.",
                f"[bright_white]{fq_name}[/bright_white]\n[cyan]{fq_info}[/cyan]"
            )

        panels.append(Panel(
            table,
            title="[bold green]FASTQ Files (Single-End)[/bold green]",
            title_align="left",
            border_style="green",
            box=box.ROUNDED,
            expand=False
        ))

    # Reference Files Panel
    if references:
        table = Table(show_header=False, box=None, padding=(0, 1), expand=False)
        table.add_column("Number", style="dim", width=4)
        table.add_column("File", no_wrap=False)

        for idx, (ref_path, ref_size) in enumerate(references, 1):
            ref_name = ref_path.split('/')[-1]
            ref_info = format_bytes(ref_size)

            table.add_row(
                f"{idx}.",
                f"[bright_white]{ref_name}[/bright_white]\n[cyan]{ref_info}[/cyan]"
            )

        panels.append(Panel(
            table,
            title="[bold green]Reference Files[/bold green]",
            title_align="left",
            border_style="green",
            box=box.ROUNDED,
            expand=False
        ))

    # SAM Files Panel
    if sam_files:
        table = Table(show_header=False, box=None, padding=(0, 1), expand=False)
        table.add_column("Number", style="dim", width=4)
        table.add_column("File", no_wrap=False)

        for idx, (sam_path, sam_size) in enumerate(sam_files, 1):
            sam_name = sam_path.split('/')[-1]
            sam_info = format_bytes(sam_size)

            table.add_row(
                f"{idx}.",
                f"[bright_white]{sam_name}[/bright_white]\n[cyan]{sam_info}[/cyan]"
            )

        panels.append(Panel(
            table,
            title="[bold green]SAM Files (Optional)[/bold green]",
            title_align="left",
            border_style="green",
            box=box.ROUNDED,
            expand=False
        ))

    # Display all panels
    if panels:
        if len(panels) > 1:
            # Check if terminal width is sufficient for side-by-side layout
            if console.width >= 140:
                # Use Table to place panels side-by-side
                layout_table = Table(show_header=False, box=None, padding=(0, 1), expand=False)
                for _ in panels:
                    layout_table.add_column(no_wrap=True)
                layout_table.add_row(*panels)
                content = layout_table
            else:
                # Stack panels vertically for narrow terminals
                from rich.console import Group
                content = Group(*panels)
        else:
            content = panels[0]

        main_panel = Panel(
            content,
            title="[bold blue]DETECTED FILES[/bold blue]",
            title_align="left",
            border_style="blue",
            box=box.DOUBLE,
            padding=(1, 2),
            expand=False
        )

        console.print(main_panel)
    else:
        console.print("[yellow]No files detected in current directory[/yellow]")


def display_file_pairs(
    fastq_pairs: List[Tuple[str, str, int, int]],
    references: List[Tuple[str, int]],
    sam_files: List[Tuple[str, int]]
) -> None:
    """Show paired-end FASTQs with reference/SAM summaries."""
    panels = []

    # Paired Reads Panel
    if fastq_pairs:
        table = Table(show_header=False, box=None, padding=(0, 1), expand=False)
        table.add_column("Number", style="dim", width=4)
        table.add_column("R1", no_wrap=False)
        table.add_column("Arrow", width=3, justify="center")
        table.add_column("R2", no_wrap=False)

        for idx, (r1_path, r2_path, r1_size, r2_size) in enumerate(fastq_pairs, 1):
            r1_name = r1_path.split('/')[-1]
            r2_name = r2_path.split('/')[-1]
            r1_info = format_bytes(r1_size)
            r2_info = format_bytes(r2_size)

            table.add_row(
                f"{idx}.",
                f"[bright_white]{r1_name}[/bright_white]\n[cyan]{r1_info}[/cyan]",
                "←→",
                f"[bright_white]{r2_name}[/bright_white]\n[cyan]{r2_info}[/cyan]"
            )

        panels.append(Panel(
            table,
            title="[bold green]Paired Reads[/bold green]",
            title_align="left",
            border_style="green",
            box=box.ROUNDED,
            expand=False
        ))

    # Reference Files Panel
    if references:
        table = Table(show_header=False, box=None, padding=(0, 1), expand=False)
        table.add_column("Number", style="dim", width=4)
        table.add_column("File", no_wrap=False)

        for idx, (ref_path, ref_size) in enumerate(references, 1):
            ref_name = ref_path.split('/')[-1]
            ref_info = format_bytes(ref_size)

            table.add_row(
                f"{idx}.",
                f"[bright_white]{ref_name}[/bright_white]\n[cyan]{ref_info}[/cyan]"
            )

        panels.append(Panel(
            table,
            title="[bold green]Reference Files[/bold green]",
            title_align="left",
            border_style="green",
            box=box.ROUNDED,
            expand=False
        ))

    # SAM Files Panel
    if sam_files:
        table = Table(show_header=False, box=None, padding=(0, 1), expand=False)
        table.add_column("Number", style="dim", width=4)
        table.add_column("File", no_wrap=False)

        for idx, (sam_path, sam_size) in enumerate(sam_files, 1):
            sam_name = sam_path.split('/')[-1]
            sam_info = format_bytes(sam_size)

            table.add_row(
                f"{idx}.",
                f"[bright_white]{sam_name}[/bright_white]\n[cyan]{sam_info}[/cyan]"
            )

        panels.append(Panel(
            table,
            title="[bold green]SAM Files (Optional)[/bold green]",
            title_align="left",
            border_style="green",
            box=box.ROUNDED,
            expand=False
        ))

    # Display all panels in a single container
    if panels:
        if len(panels) > 1:
            # Check if terminal width is sufficient for side-by-side layout
            # Use 140 characters as threshold (comfortable for two panels)
            if console.width >= 140:
                # Use Table to place panels side-by-side without expansion
                layout_table = Table(show_header=False, box=None, padding=(0, 1), expand=False)
                for _ in panels:
                    layout_table.add_column(no_wrap=True)
                layout_table.add_row(*panels)
                content = layout_table
            else:
                # Stack panels vertically for narrow terminals
                from rich.console import Group
                content = Group(*panels)
        else:
            content = panels[0]

        # Create outer panel - should not expand if content doesn't expand
        main_panel = Panel(
            content,
            title="[bold blue]DETECTED FILES[/bold blue]",
            title_align="left",
            border_style="blue",
            box=box.DOUBLE,
            padding=(1, 2),
            expand=False
        )

        console.print(main_panel)
    else:
        console.print("[yellow]No files detected in current directory[/yellow]")


def display_validation_results(validations: List[Dict]) -> bool:
    """Render validation results; return True if no blocking errors."""
    if not validations:
        return True

    lines = []
    has_errors = False

    for result in validations:
        status = result["status"]
        message = result["message"]
        blocking = result.get("blocking", False)

        # Messages already include symbols (✓, ⚠, ✗), just apply white color
        lines.append(f"[bright_white]{message}[/bright_white]")

        if status == "error" and blocking:
            has_errors = True

    panel = Panel(
        "\n".join(lines),
        title="[bold green]Validation Results[/bold green]",
        title_align="left",
        border_style="green",
        box=box.ROUNDED,
        padding=(1, 2),
        expand=False
    )

    console.print(panel)

    return not has_errors


def display_pipeline_config(
    r1_file: str,
    r2_file: Optional[str],
    ref_file: str,
    read_count: int,
    ref_sequences: int,
    ref_size_kb: float,
    aligner: str,
    threads: int,
    sam_file: Optional[str] = None,
    minimap2_profile: Optional[str] = None,
    output_dir: str = "./output",
    iterations: int = 0,
) -> None:
    """Display pipeline configuration summary (paired or single-end)."""
    lines = []

    # INPUT FILES section
    lines.append("[bold]INPUT FILES[/bold]")

    # Display reads based on mode
    if r2_file:
        # Paired-end
        lines.append(f"  Reads:      [cyan]{r1_file}[/cyan] / [cyan]{r2_file}[/cyan] [dim]({read_count:,} pairs)[/dim]")
    else:
        # Single-end
        lines.append(f"  Reads:      [cyan]{r1_file}[/cyan] [dim]({read_count:,} reads, single-end)[/dim]")

    lines.append(f"  Reference:  [cyan]{ref_file}[/cyan] [dim]({ref_sequences} sequences, {ref_size_kb:.1f} Kbp)[/dim]")

    if sam_file:
        lines.append(f"  SAM:        [cyan]{sam_file}[/cyan] [dim](using existing alignment)[/dim]")
    else:
        # Match capitalization from aligner selection UI
        aligner_map = {
            "bwa-mem2": "BWA-MEM2",
            "bowtie2": "Bowtie2",
            "minimap2": "Minimap2"
        }
        profile_map = {
            "short": "Short",
            "pacbio-hifi": "PacBio-HiFi",
            "pacbio-clr": "PacBio-CLR",
            "ont-q20": "ONT-Q20",
            "ont-standard": "ONT-Standard",
        }
        aligner_display = aligner_map.get(aligner, aligner.upper())
        if minimap2_profile:
            profile_display = profile_map.get(minimap2_profile, minimap2_profile.upper())
            aligner_display += f" ({profile_display})"
        lines.append(f"  Aligner:    [cyan]{aligner_display}[/cyan]")
        lines.append(f"  Threads:    [cyan]{threads}[/cyan] cores")
    
    if iterations > 0:
        lines.append(f"  Iterations: [cyan]{iterations}[/cyan] [dim](iterative refinement enabled)[/dim]")
    else:
        lines.append(f"  Iterations: [dim]0 (disabled)[/dim]")

    lines.append("")

    # PIPELINE STEPS section
    lines.append("[bold]PIPELINE STEPS[/bold]")

    # Use the same display name for index build step
    aligner_display_name = aligner_map.get(aligner, aligner.upper())
    steps = [
        ("1", "Multiple Sequence Alignment"),
        ("2", f"{aligner_display_name} Index Build"),
        ("3", "Read Alignment"),
        ("4", "Mapping of Reference Sequences to MSA"),
        ("5", "Mapping of Reads to Reference Sequences"),
        ("6", "Mapping Refinement"),
        ("7", "Output File Creation")
    ]

    # Gray out steps 2-3 if using existing SAM
    for num, step_name in steps:
        if sam_file and num in ["2", "3"]:
            lines.append(f"  [dim]{num}. ⚙  {step_name} (skipped)[/dim]")
        else:
            lines.append(f"  {num}. [cyan]⚙[/cyan]  {step_name}")

    lines.append("")

    # OUTPUT DIRECTORY section
    lines.append("[bold]OUTPUT DIRECTORY[/bold]")
    output_path = Path(output_dir)
    output_dir_display = str(output_path.resolve())
    lines.append(f"  [cyan]{output_dir_display}/[/cyan]")
    lines.append("  ├─ [green]unique_mappings.tsv[/green]")
    lines.append(f"  ├─ [green]fastq/PKD1*.fq[/green] [dim]({ref_sequences} files)[/dim]")
    lines.append(f"  ├─ [green]bam/PKD1*.bam + .bai[/green] [dim]({ref_sequences * 2} files)[/dim]")
    lines.append("  └─ [dim]pipeline_YYYYMMDD_HHMMSS.log[/dim]")

    # Create panel
    panel = Panel(
        "\n".join(lines),
        title="[bold blue]PIPELINE CONFIGURATION[/bold blue]",
        title_align="left",
        border_style="blue",
        box=box.DOUBLE,
        padding=(1, 2),
        expand=False
    )

    console.print(panel)


def print_header():
    """Print the mapper ASCII art header."""
    console.clear()
    header = Text("""
██████╗  █████╗ ██████╗  █████╗ ██████╗ ██╗███████╗███╗   ███╗
██╔══██╗██╔══██╗██╔══██╗██╔══██╗██╔══██╗██║██╔════╝████╗ ████║
██████╔╝███████║██████╔╝███████║██║  ██║██║███████╗██╔████╔██║
██╔═══╝ ██╔══██║██╔══██╗██╔══██║██║  ██║██║╚════██║██║╚██╔╝██║
██║     ██║  ██║██║  ██║██║  ██║██████╔╝██║███████║██║ ╚═╝ ██║
╚═╝     ╚═╝  ╚═╝╚═╝  ╚═╝╚═╝  ╚═╝╚═════╝ ╚═╝╚══════╝╚═╝     ╚═╝
    """, style="bold blue")
    console.print(header)


def print_section(title: str):
    """Print a section header."""
    console.print()
    console.print(f"[bold cyan]▶ {title}[/bold cyan]")
    console.print()


def print_success(message: str):
    """Print a success message."""
    console.print(f"[green]✓[/green] {message}")


def print_error(message: str):
    """Print an error message."""
    console.print(f"[red]✗[/red] {message}")


def print_warning(message: str):
    """Print a warning message."""
    console.print(f"[yellow]⚠[/yellow] {message}")
