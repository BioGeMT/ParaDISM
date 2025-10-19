#!/usr/bin/env python3
"""
PKD1 Mapper Pipeline - Unified Entry Point

This script handles:
1. Interactive mode with rich UI (no arguments)
2. Direct CLI mode (with arguments)
3. Complete pipeline execution (MSA → alignment → mapping → output)

Usage:
    python mapper.py                           # Interactive mode
    python mapper.py --read1 R1 --read2 R2 --reference REF --aligner ALIGNER --threads N
"""

import sys
import os
import argparse
import subprocess
import tempfile
from pathlib import Path
from datetime import datetime

# Add scripts directory to path
SCRIPT_DIR = Path(__file__).parent.absolute()
SCRIPTS_DIR = SCRIPT_DIR / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))

from file_scanner import (
    scan_fastq_metadata,
    scan_fasta_metadata,
    scan_sam_metadata,
    find_paired_fastqs,
    find_references,
    find_sam_files
)
from validators import (
    validate_fastq_pair,
    validate_fasta,
    validate_sam
)
from ui_components import (
    display_file_pairs,
    display_validation_results,
    display_pipeline_config,
    print_header,
    print_section,
    console
)


class PipelineExecutor:
    """Executes the PKD1 mapper pipeline steps."""

    def __init__(self, output_dir="./output", scripts_dir=None):
        self.output_dir = Path(output_dir)
        self.scripts_dir = Path(scripts_dir) if scripts_dir else SCRIPTS_DIR
        self.output_dir.mkdir(exist_ok=True)

        # Create log file
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        self.log_file = self.output_dir / f"pipeline_{timestamp}.log"

        # Spinner characters
        self.spinner_chars = '⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏'
        self.spinner_thread = None
        self.spinner_stop = None

    def log(self, message):
        """Write message to log file."""
        with open(self.log_file, 'a') as f:
            timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            f.write(f"\n{'='*48}\n")
            f.write(f"[{timestamp}] {message}\n")
            f.write(f"{'='*48}\n")

    def run_spinner(self, message):
        """Run spinner animation in background thread."""
        import threading
        import time

        self.spinner_stop = threading.Event()

        def spin():
            i = 0
            while not self.spinner_stop.is_set():
                char = self.spinner_chars[i % len(self.spinner_chars)]
                # Use carriage return to overwrite line
                print(f"\r{char} {message}", end="", file=sys.stderr)
                sys.stderr.flush()
                time.sleep(0.1)
                i += 1

        self.spinner_thread = threading.Thread(target=spin, daemon=True)
        self.spinner_thread.start()

    def stop_spinner(self, message, success=True):
        """Stop spinner and show final status."""
        if self.spinner_stop:
            self.spinner_stop.set()
        if self.spinner_thread:
            self.spinner_thread.join()

        if success:
            print(f"\r\033[0;36m✓\033[0m \033[0;36m{message}\033[0m", file=sys.stderr)
        else:
            print(f"\r\033[0;31m✗\033[0m \033[0;31m{message}\033[0m", file=sys.stderr)

    def run_with_spinner(self, cmd, message, shell=False):
        """Run command with spinner animation."""
        self.log(message)
        self.run_spinner(message)

        try:
            result = subprocess.run(
                cmd,
                shell=shell,
                capture_output=True,
                text=True,
                check=True
            )

            # Log output
            with open(self.log_file, 'a') as f:
                if result.stdout:
                    f.write(result.stdout)
                if result.stderr:
                    f.write(result.stderr)

            self.stop_spinner(message, success=True)
            return result

        except subprocess.CalledProcessError as e:
            self.stop_spinner(message, success=False)
            print("Error occurred. Check log: " + str(self.log_file), file=sys.stderr)
            print("Last 20 lines:", file=sys.stderr)
            with open(self.log_file, 'r') as f:
                lines = f.readlines()
                for line in lines[-20:]:
                    print("  " + line.rstrip(), file=sys.stderr)
            sys.exit(1)

    def run_with_progress(self, cmd, message):
        """Run command with progress output (for Python scripts)."""
        self.log(message)

        import codecs
        import select
        import time

        spinner_interval = 0.1
        spinner_index = 0
        last_refresh = 0.0
        display_progress = ""
        progress_fragment = ""

        def render_frame(spinner_char, progress_line):
            """Redraw spinner and progress lines without flooding the terminal."""
            formatted_progress = f"  {progress_line}" if progress_line else ""
            sys.stderr.write("\033[2A\r")
            sys.stderr.write(f"{spinner_char} {message}\033[K\n")
            if formatted_progress:
                sys.stderr.write(f"{formatted_progress}\033[K\n")
            else:
                sys.stderr.write("\033[K\n")
            sys.stderr.flush()

        # Prepare screen space for spinner + progress
        sys.stderr.write("\033[?25l")  # hide cursor
        sys.stderr.write(f"{self.spinner_chars[0]} {message}\n\n")
        sys.stderr.flush()

        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )

        stdout_fd = process.stdout.fileno() if process.stdout else None
        stderr_fd = process.stderr.fileno() if process.stderr else None
        stdout_done = stdout_fd is None
        stderr_done = stderr_fd is None

        stdout_decoder = codecs.getincrementaldecoder("utf-8")() if not stdout_done else None
        stderr_decoder = codecs.getincrementaldecoder("utf-8")() if not stderr_done else None

        exit_code = None

        try:
            with open(self.log_file, 'a') as log_handle:
                while True:
                    read_fds = []
                    if not stdout_done:
                        read_fds.append(stdout_fd)
                    if not stderr_done:
                        read_fds.append(stderr_fd)

                    if read_fds:
                        ready, _, _ = select.select(read_fds, [], [], spinner_interval)
                    else:
                        ready = []
                        time.sleep(spinner_interval)

                    now = time.time()
                    if now - last_refresh >= spinner_interval:
                        spinner_char = self.spinner_chars[spinner_index % len(self.spinner_chars)]
                        spinner_index += 1
                        render_frame(spinner_char, display_progress)
                        last_refresh = now

                    for fd in ready:
                        if fd == stdout_fd and not stdout_done:
                            chunk = os.read(stdout_fd, 4096)
                            if not chunk:
                                stdout_done = True
                                remaining = stdout_decoder.decode(b"", final=True) if stdout_decoder else ""
                                if remaining:
                                    log_handle.write(remaining)
                                    log_handle.flush()
                                    for ch in remaining:
                                        if ch in ("\r", "\n"):
                                            if progress_fragment:
                                                display_progress = progress_fragment
                                            progress_fragment = ""
                                        else:
                                            progress_fragment += ch
                                    if progress_fragment:
                                        display_progress = progress_fragment
                                        progress_fragment = ""
                            else:
                                text = stdout_decoder.decode(chunk)
                                if text:
                                    log_handle.write(text)
                                    log_handle.flush()
                                    for ch in text:
                                        if ch in ("\r", "\n"):
                                            if progress_fragment:
                                                display_progress = progress_fragment
                                            progress_fragment = ""
                                        else:
                                            progress_fragment += ch
                                    if progress_fragment:
                                        display_progress = progress_fragment
                        elif fd == stderr_fd and not stderr_done:
                            chunk = os.read(stderr_fd, 4096)
                            if not chunk:
                                stderr_done = True
                                remaining = stderr_decoder.decode(b"", final=True) if stderr_decoder else ""
                                if remaining:
                                    log_handle.write(remaining)
                                    log_handle.flush()
                            else:
                                text = stderr_decoder.decode(chunk)
                                if text:
                                    log_handle.write(text)
                                    log_handle.flush()

                    if process.poll() is not None and stdout_done and stderr_done:
                        break

            # Flush any trailing fragments into display
            if progress_fragment:
                display_progress = progress_fragment

            exit_code = process.wait()

            status_line = (
                f"\033[0;36m✓\033[0m \033[0;36m{message}\033[0m"
                if exit_code == 0
                else f"\033[0;31m✗\033[0m \033[0;31m{message}\033[0m"
            )

            formatted_progress_line = f"  {display_progress}" if exit_code != 0 and display_progress else ""

            sys.stderr.write("\033[2A\r")
            sys.stderr.write(f"{status_line}\033[K\n")
            if formatted_progress_line:
                sys.stderr.write(f"{formatted_progress_line}\033[K\n")
            else:
                sys.stderr.write("\033[K")
                sys.stderr.write("\r")
            sys.stderr.flush()

        finally:
            sys.stderr.write("\033[?25h")  # show cursor
            sys.stderr.flush()

        if exit_code != 0:
            print("Error occurred. Check log: " + str(self.log_file), file=sys.stderr)
            print("Last 20 lines:", file=sys.stderr)
            with open(self.log_file, 'r') as f:
                lines = f.readlines()
                for line in lines[-20:]:
                    print("  " + line.rstrip(), file=sys.stderr)
            sys.exit(1)

    def run_pipeline(self, r1, r2, ref, aligner="bwa-mem2", threads=4,
                     sam=None, minimap2_profile="short", show_header=True):
        """Execute the complete pipeline."""

        if show_header:
            print("\n\033[0;36mRunning mapper pipeline...\033[0m\n", file=sys.stderr)

        # Step 1: MSA
        msa_output = self.output_dir / "ref_seq_msa.aln"
        self.run_with_spinner(
            f"mafft --auto '{ref}' > '{msa_output}'",
            "Creating Multiple Sequence Alignment",
            shell=True
        )

        # Step 2-3: Alignment (or copy SAM)
        sam_output = self.output_dir / "mapped_reads.sam"

        if sam:
            self.run_with_spinner(
                ["cp", str(sam), str(sam_output)],
                "Copying provided SAM file"
            )
        else:
            if aligner == "bowtie2":
                index_base = self.output_dir / "ref_index"
                self.run_with_spinner(
                    ["bowtie2-build", str(ref), str(index_base)],
                    "Building Bowtie2 index"
                )
                self.run_with_spinner(
                    f"bowtie2 -p {threads} -x '{index_base}' -1 '{r1}' -2 '{r2}' -S '{sam_output}'",
                    "Aligning reads with Bowtie2",
                    shell=True
                )
            elif aligner == "bwa-mem2":
                index_base = self.output_dir / "ref_index"
                self.run_with_spinner(
                    ["bwa-mem2", "index", "-p", str(index_base), str(ref)],
                    "Building BWA-MEM2 index"
                )
                self.run_with_spinner(
                    f"bwa-mem2 mem -t {threads} '{index_base}' '{r1}' '{r2}' > '{sam_output}'",
                    "Aligning reads with BWA-MEM2",
                    shell=True
                )
            elif aligner == "minimap2":
                index_file = self.output_dir / "ref_index.mmi"

                preset_map = {
                    "short": "sr",
                    "pacbio": "map-pb",
                    "nanopore": "map-ont"
                }
                preset = preset_map.get(minimap2_profile, "sr")

                self.run_with_spinner(
                    ["minimap2", "-x", preset, "-d", str(index_file), str(ref)],
                    f"Building minimap2 ({minimap2_profile}) index"
                )
                self.run_with_spinner(
                    f"minimap2 -ax {preset} --MD -t {threads} '{index_file}' '{r1}' '{r2}' > '{sam_output}'",
                    f"Aligning reads with minimap2 ({minimap2_profile})",
                    shell=True
                )

        # Step 4: Map reference to MSA
        ref_msa_tsv = self.output_dir / "ref_seq_msa.tsv"
        self.run_with_spinner(
            [
                "python", str(self.scripts_dir / "ref_2_msa.py"),
                "--reference_fasta", str(ref),
                "--msa_file", str(msa_output),
                "--output", str(ref_msa_tsv)
            ],
            "Mapping reference sequences to MSA"
        )

        # Step 5: Map reads to genes
        mapped_reads_tsv = self.output_dir / "mapped_reads.tsv"
        self.run_with_progress(
            [
                "python", str(self.scripts_dir / "read_2_gene.py"),
                "--sam", str(sam_output),
                "--output", str(mapped_reads_tsv),
                "--fastq", str(r1)
            ],
            "Mapping reads to reference sequences"
        )

        # Step 6: Unique mapping algorithm
        unique_mappings_tsv = self.output_dir / "unique_mappings.tsv"
        self.run_with_progress(
            [
                "python", str(self.scripts_dir / "mapper_algo_snp_only.py"),
                "--read_map", str(mapped_reads_tsv),
                "--msa", str(ref_msa_tsv),
                "--output_file", str(unique_mappings_tsv)
            ],
            "Refining mapping"
        )

        # Step 7: Generate outputs
        fastq_dir = self.output_dir / "fastq"
        bam_dir = self.output_dir / "bam"
        self.run_with_progress(
            [
                "python", str(self.scripts_dir / "output.py"),
                "--tsv", str(unique_mappings_tsv),
                "--r1", str(r1),
                "--r2", str(r2),
                "--ref", str(ref),
                "--fastq-dir", str(fastq_dir),
                "--bam-dir", str(bam_dir),
                "--aligner", aligner,
                "--threads", str(threads),
                "--minimap2-profile", minimap2_profile
            ],
            "Writing output files"
        )

        # Brief delay to let output buffers flush
        import time
        time.sleep(0.2)

        # Cleanup intermediate files
        self.log("Cleaning up intermediate files")
        for pattern in ["ref_index.*", "ref_seq_msa.aln", "ref_seq_msa.tsv", "mapped_reads.tsv"]:
            for f in self.output_dir.glob(pattern):
                f.unlink()

        print("\033[0;36m✓\033[0m \033[0;36mCleaning up intermediate files\033[0m", file=sys.stderr)

        print(f"\nPipeline complete. Outputs saved to: {self.output_dir}\n", file=sys.stderr)


def interactive_mode():
    """Run mapper in interactive mode with rich UI."""

    print_header()

    # ========================================================================
    # File Detection
    # ========================================================================

    print_section("Input Files")

    # Find files
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

    # Display detected files
    display_file_pairs(fastq_pairs, references, sam_files)

    # Get all individual FASTQ files for manual selection fallback
    all_fastq_files = sorted(list(Path(".").glob("*.fq")) + list(Path(".").glob("*.fastq")))

    # Check minimum requirements
    if not fastq_pairs and len(all_fastq_files) < 2:
        console.print("\n[red]✗ No FASTQ files found![/red]")
        console.print("Please ensure you have FASTQ files in the current directory.\n")
        sys.exit(1)

    if not references:
        console.print("\n[red]✗ No reference FASTA files found![/red]")
        console.print("Please ensure you have a reference FASTA file in the current directory.\n")
        sys.exit(1)

    # ========================================================================
    # File Selection
    # ========================================================================

    console.print()

    # Select FASTQ pair
    if not fastq_pairs:
        # Manual selection mode
        console.print("[yellow]⚠ No auto-paired FASTQ files detected.[/yellow]")
        console.print("[yellow]  Please manually select R1 and R2 files.[/yellow]\n")

        console.print("[cyan]Available FASTQ files:[/cyan]")
        for idx, fq_file in enumerate(all_fastq_files, 1):
            fq_meta = scan_fastq_metadata(str(fq_file))
            console.print(f"  [green]{idx}[/green]. {fq_file.name} [dim]({fq_meta['size_human']}, {fq_meta['read_count']:,} reads)[/dim]")

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

    elif len(fastq_pairs) == 1:
        r1_path, r2_path, r1_meta, r2_meta = fastq_pairs[0]
        console.print(f"[green]✓[/green] Auto-selected: [cyan]{Path(r1_path).name}[/cyan] / [cyan]{Path(r2_path).name}[/cyan]")
        console.print()  # Empty line after selection
    else:
        # Multiple auto-paired options
        while True:
            choice = console.input(f"[green]Select read pair [1-{len(fastq_pairs)}], 'm' for manual, or ENTER for first:[/green] ").strip().lower()

            if choice == 'm':
                # Switch to manual mode
                console.print()
                console.print("[cyan]Available FASTQ files:[/cyan]")
                for idx, fq_file in enumerate(all_fastq_files, 1):
                    fq_meta = scan_fastq_metadata(str(fq_file))
                    console.print(f"  [green]{idx}[/green]. {fq_file.name} [dim]({fq_meta['size_human']}, {fq_meta['read_count']:,} reads)[/dim]")
                console.print()

                # Select R1
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

                # Select R2
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
                # Select from auto-paired list
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

    console.print()  # Empty line after FASTQ selection

    # Select reference
    # Get all reference files for manual selection
    all_ref_files = sorted(list(Path(".").glob("*.fa")) + list(Path(".").glob("*.fasta")) +
                          list(Path(".").glob("*.fas")) + list(Path(".").glob("*.fna")))

    if len(references) == 1:
        ref_path, ref_meta = references[0]
        console.print(f"[green]✓[/green] Auto-selected: [cyan]{Path(ref_path).name}[/cyan]")
        console.print()  # Empty line after selection
    else:
        while True:
            choice = console.input(f"[green]Select reference [1-{len(references)}], 'm' to see all FASTA files, or ENTER for first:[/green] ").strip().lower()

            if choice == 'm':
                # Show all FASTA files (not just those in references list)
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

    # SAM file selection (optional)
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

    # ========================================================================
    # Aligner Selection (if not using SAM)
    # ========================================================================

    aligner = "bwa-mem2"
    minimap2_profile = "short"
    threads = 4

    if not sam_path:
        from rich.table import Table
        from rich.panel import Panel
        from rich import box as rbox
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

        # Thread selection
        max_threads = os.cpu_count() or 4
        console.print()
        while True:
            choice = console.input(f"[green]Number of threads [1-{max_threads}] or ENTER for 4:[/green] ").strip()
            choice = choice if choice else "4"

            try:
                threads = int(choice)
                if 1 <= threads <= max_threads:
                    console.print(f"[green]✓[/green] Using [cyan]{threads}[/cyan] threads")
                    break
                else:
                    console.print(f"[red]Invalid thread count. Please enter 1-{max_threads}[/red]")
            except ValueError:
                console.print(f"[red]Invalid thread count. Please enter 1-{max_threads}[/red]")

    # ========================================================================
    # Validation
    # ========================================================================

    console.print()
    print_section("Input Validation")

    all_validations = []
    all_validations.extend(validate_fastq_pair(r1_path, r2_path, r1_meta, r2_meta))
    all_validations.extend(validate_fasta(ref_path, ref_meta))
    if sam_path:
        all_validations.extend(validate_sam(sam_path, scan_sam_metadata(sam_path)))

    can_proceed = display_validation_results(all_validations)

    if not can_proceed:
        console.print("\n[red]✗ Validation failed! Cannot proceed.[/red]\n")
        sys.exit(1)

    # ========================================================================
    # Configuration Display & Confirmation
    # ========================================================================

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

    # ========================================================================
    # Run Pipeline
    # ========================================================================

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

    # Success message
    console.print("[green]═══════════════════════════════════════════════════════[/green]")
    console.print("[green]Pipeline Complete![/green]")
    console.print("[green]═══════════════════════════════════════════════════════[/green]")
    console.print()
    console.print("[cyan]Outputs:[/cyan]")
    console.print("  • Unique mappings: [green]./output/unique_mappings.tsv[/green]")
    console.print("  • Gene-specific FASTQs: [green]./output/fastq/[/green]")
    console.print("  • Gene-specific BAMs: [green]./output/bam/[/green]")
    console.print()


def cli_mode(args):
    """Run mapper in CLI mode with provided arguments."""

    # Validate files exist
    if not Path(args.read1).exists():
        console.print(f"[red]✗ R1 file not found: {args.read1}[/red]")
        sys.exit(1)
    if not Path(args.read2).exists():
        console.print(f"[red]✗ R2 file not found: {args.read2}[/red]")
        sys.exit(1)
    if not Path(args.reference).exists():
        console.print(f"[red]✗ Reference file not found: {args.reference}[/red]")
        sys.exit(1)
    if args.sam and not Path(args.sam).exists():
        console.print(f"[red]✗ SAM file not found: {args.sam}[/red]")
        sys.exit(1)

    # Run pipeline directly
    executor = PipelineExecutor(output_dir=args.output_dir)
    executor.run_pipeline(
        r1=args.read1,
        r2=args.read2,
        ref=args.reference,
        aligner=args.aligner,
        threads=args.threads,
        sam=args.sam,
        minimap2_profile=args.minimap2_profile,
        show_header=True
    )


def main():
    parser = argparse.ArgumentParser(
        description="PKD1 Mapper Pipeline - Read mapping and gene-specific output generation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Interactive mode (recommended)
  python mapper.py

  # Direct CLI mode
  python mapper.py --read1 r1.fq --read2 r2.fq --reference ref.fa

  # Use existing SAM file
  python mapper.py --read1 r1.fq --read2 r2.fq --reference ref.fa --sam mapped.sam

  # Specify aligner and threads
  python mapper.py --read1 r1.fq --read2 r2.fq --reference ref.fa --aligner bowtie2 --threads 8
        """
    )

    parser.add_argument("--read1", help="R1 FASTQ file")
    parser.add_argument("--read2", help="R2 FASTQ file")
    parser.add_argument("--reference", help="Reference FASTA file")
    parser.add_argument("--sam", help="Existing SAM file (optional, skips alignment)")
    parser.add_argument("--aligner", default="bwa-mem2", choices=["bowtie2", "bwa-mem2", "minimap2"],
                        help="Read aligner (default: bwa-mem2)")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads (default: 4)")
    parser.add_argument("--minimap2-profile", default="short", choices=["short", "pacbio", "nanopore"],
                        help="Minimap2 profile (default: short)")
    parser.add_argument("--output-dir", default="./output", help="Output directory (default: ./output)")

    args = parser.parse_args()

    # If all required arguments provided, run in CLI mode
    if args.read1 and args.read2 and args.reference:
        cli_mode(args)
    else:
        # Otherwise run interactive mode
        if any([args.read1, args.read2, args.reference]):
            console.print("[yellow]⚠ Partial arguments provided. Use --read1, --read2, and --reference together for CLI mode.[/yellow]")
            console.print("[yellow]  Launching interactive mode instead...[/yellow]\n")

        interactive_mode()


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        console.print("\n\n[yellow]Pipeline cancelled by user[/yellow]\n")
        sys.exit(130)
