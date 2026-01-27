#!/usr/bin/env python3
"""Simple ParaDISM executor: alignment → MSA → ParaDISM algorithm."""

from __future__ import annotations

import sys
import time
import subprocess
import threading
import shutil
import gzip
from pathlib import Path
from typing import Sequence
from Bio import SeqIO

from utils.logger import PipelineLogger
from utils.progress import ProgressRunner
from .paradism_algo import load_msa, process_sam_to_dict, write_fastq_outputs, create_bam_files

PIPELINE_DIR = Path(__file__).resolve().parent


class SimpleParaDISMExecutor:
    """Simple executor that runs alignment, MSA, and ParaDISM algorithm."""

    def __init__(self, output_dir: str | Path = "./output", pipeline_dir: str | Path | None = None, prefix: str | None = None, min_alternate_count: int = 5, add_quality_filters: bool = False, qual_threshold: int = 20, dp_threshold: int = 10, af_threshold: float = 0.05) -> None:
        self.output_dir = Path(output_dir)
        self.pipeline_dir = Path(pipeline_dir) if pipeline_dir else PIPELINE_DIR
        self.output_dir.mkdir(exist_ok=True)
        self.min_alternate_count = min_alternate_count
        self.add_quality_filters = add_quality_filters
        self.qual_threshold = qual_threshold
        self.dp_threshold = dp_threshold
        self.af_threshold = af_threshold

        # Use provided prefix, or extract from output directory name
        if prefix is not None:
            self.prefix = prefix
        else:
            self.prefix = self.output_dir.name if self.output_dir.name != "." else "output"

        timestamp = time.strftime('%Y%m%d_%H%M%S')
        self.log_file = self.output_dir / f"{self.prefix}_pipeline_{timestamp}.log"
        self.logger = PipelineLogger(self.log_file)
        self.progress = ProgressRunner(self.logger)

    def _run_spinner(self, command: Sequence[str] | str | callable, message: str, *, shell: bool = False):
        """Run command or callable with spinner. Returns result if callable."""
        if callable(command):
            # It's a Python function - run with spinner wrapper
            stop_event = threading.Event()
            
            def spin():
                index = 0
                while not stop_event.is_set():
                    char = self.progress.spinner_chars[index % len(self.progress.spinner_chars)]
                    print(f"\r  {char} {message}", end="", file=sys.stderr)
                    sys.stderr.flush()
                    index += 1
                    time.sleep(0.1)
            
            spinner_thread = threading.Thread(target=spin, daemon=True)
            spinner_thread.start()
            
            try:
                result = command()
                stop_event.set()
                spinner_thread.join()
                print(f"\r  \033[0;36m✓ {message}\033[0m", file=sys.stderr)
                return result
            except Exception:
                stop_event.set()
                spinner_thread.join()
                print(f"\r  \033[0;31m✗ {message}\033[0m", file=sys.stderr)
                raise
        else:
            # It's a command - use progress runner
            self.progress.run_with_spinner(command, message, shell=shell)
            return None

    def _extract_none_reads_from_assignments(self, assignments: dict[str, str]) -> set[str]:
        """Extract read IDs that mapped to NONE from assignments dict."""
        return {read_name for read_name, gene in assignments.items() if gene == "NONE"}

    def _extract_reads_from_fastq(self, fastq_path: Path, read_ids: set[str], output_path: Path) -> int:
        """Extract reads from FASTQ file based on read IDs. Returns count of extracted reads."""
        output_path.parent.mkdir(parents=True, exist_ok=True)
        count = 0
        
        # Handle compressed input files
        fastq_handle = gzip.open(fastq_path, "rt") if str(fastq_path).endswith(".gz") else open(fastq_path, "rt")
        
        with open(output_path, 'w') as outfile:
            try:
                for rec in SeqIO.parse(fastq_handle, "fastq"):
                    # Normalize read ID (remove /1 or /2 suffix for matching)
                    base_id = rec.id
                    if base_id.endswith("/1") or base_id.endswith("/2"):
                        base_id = base_id[:-2]
                    
                    # Check if this read should be extracted
                    if base_id in read_ids or rec.id in read_ids:
                        SeqIO.write(rec, outfile, "fastq")
                        count += 1
            finally:
                fastq_handle.close()
        return count

    def _merge_assignments(self, previous_assignments: dict[str, str], new_assignments: dict[str, str]) -> dict[str, str]:
        """Merge assignments: keep previous non-NONE mappings + new mappings from NONE reads."""
        merged = {}
        
        # Keep all previous non-NONE mappings
        for read_name, gene in previous_assignments.items():
            if gene != "NONE":
                merged[read_name] = gene
        
        # Add new mappings for reads that were NONE in previous iteration
        previous_none = {read_name for read_name, gene in previous_assignments.items() if gene == "NONE"}
        for read_name, gene in new_assignments.items():
            if read_name in previous_none:
                merged[read_name] = gene
        
        return merged

    def _run_single_iteration_refinement(
        self,
        r1: str,
        r2: str | None,
        previous_ref: Path,
        previous_bam_dir: Path,
        previous_assignments: dict[str, str],
        aligner: str,
        threads: int,
        minimap2_profile: str,
        is_paired: bool,
        iteration: int,
        bowtie2_score_min: str = "G,40,40",
        bwa_min_score: int = 240,
        minimap2_min_score: int = 240,
    ) -> tuple[Path, dict[str, str], bool]:
        """
        Run one iteration of refinement on NONE reads only.
        
        Returns:
            Tuple of (updated_reference_path, merged_assignments, converged)
            where converged=True indicates no progress was made
        """
        print(f"\n  \033[0;33m=== Iteration {iteration} ===\033[0m\n", file=sys.stderr)
        
        # Create iteration directory
        iter_output_dir = self.output_dir / f"iteration_{iteration}"
        iter_output_dir.mkdir(exist_ok=True)
        
        # Log iteration start to main log
        self.logger.section(f"Iteration {iteration}")
        
        # 1. Identify NONE reads
        none_reads = self._extract_none_reads_from_assignments(previous_assignments)
        
        if len(none_reads) == 0:
            print(f"  No NONE reads to refine. Pipeline has converged.", file=sys.stderr)
            return previous_ref, previous_assignments, True
        
        # 2. Extract NONE reads from original FASTQ files
        none_r1_path = iter_output_dir / "none_reads_r1.fq"
        none_r2_path = iter_output_dir / "none_reads_r2.fq" if is_paired else None
        
        r1_count = self._extract_reads_from_fastq(Path(r1), none_reads, none_r1_path)
        r2_count = 0
        if is_paired and r2:
            r2_count = self._extract_reads_from_fastq(Path(r2), none_reads, none_r2_path)
        
        if r1_count == 0:
            print(f"  No reads extracted. Pipeline has converged.", file=sys.stderr)
            return previous_ref, previous_assignments, True
        
        # 3. Call variants from successfully mapped reads and apply to reference
        variant_calling_dir = iter_output_dir / "variant_calling"
        variant_calling_dir.mkdir(exist_ok=True)
        per_gene_vcf_dir = variant_calling_dir / "per_gene_vcfs"
        per_gene_vcf_dir.mkdir(exist_ok=True)
        
        vcf_output = variant_calling_dir / "variants.vcf"
        updated_ref = variant_calling_dir / "updated_ref.fa"
        
        variant_cmd = [
            "python",
            str(self.pipeline_dir / "variant_refinement.py"),
            "--bam-dir", str(previous_bam_dir),
            "--reference", str(previous_ref),
            "--output-vcf", str(vcf_output),
            "--output-ref", str(updated_ref),
            "--per-gene-vcf-dir", str(per_gene_vcf_dir),
            "--min-alternate-count", str(self.min_alternate_count),
            "--qual-threshold", str(self.qual_threshold),
            "--dp-threshold", str(self.dp_threshold),
            "--af-threshold", str(self.af_threshold),
        ]
        if self.add_quality_filters:
            variant_cmd.append("--add-qfilters")
        
        def _call_variants():
            result = subprocess.run(variant_cmd, capture_output=True, text=True)
            self.logger.write(result.stdout or "")
            self.logger.write(result.stderr or "")
            
            if result.returncode == 2:
                # No variants found - converged
                return None  # Signal convergence
            elif result.returncode != 0:
                raise RuntimeError("Variant calling failed")
            return True
        
        variant_result = self._run_spinner(_call_variants, "Calling variants from mapped reads")
        if variant_result is None:
            print(f"  No variants found. Stopping refinement at iteration {iteration}.", file=sys.stderr)
            return previous_ref, previous_assignments, True
        
        # 4. Re-align NONE reads with updated reference
        iter_fastq_dir = iter_output_dir / f"{self.prefix}_fastq"
        iter_bam_dir = iter_output_dir / f"{self.prefix}_bam"
        
        # Create MSA for updated reference
        iter_msa = iter_output_dir / "ref_seq_msa.aln"
        self._run_spinner(
            f"mafft --auto '{updated_ref}' > '{iter_msa}'",
            "Creating Multiple Sequence Alignment",
            shell=True,
        )
        
        # Align NONE reads
        iter_sam = iter_output_dir / "mapped_reads.sam"
        if aligner == "bowtie2":
            iter_index = iter_output_dir / "ref_index"
            self._run_spinner(
                ["bowtie2-build", str(updated_ref), str(iter_index)],
                "Building Bowtie2 index",
            )
            if is_paired:
                self._run_spinner(
                    f"bowtie2 --local --score-min {bowtie2_score_min} -p {threads} -x '{iter_index}' -1 '{none_r1_path}' -2 '{none_r2_path}' -S '{iter_sam}'",
                    "Aligning reads with Bowtie2",
                    shell=True,
                )
            else:
                self._run_spinner(
                    f"bowtie2 --local --score-min {bowtie2_score_min} -p {threads} -x '{iter_index}' -U '{none_r1_path}' -S '{iter_sam}'",
                    "Aligning reads with Bowtie2",
                    shell=True,
                )
        elif aligner == "bwa-mem2":
            iter_index = iter_output_dir / "ref_index"
            self._run_spinner(
                ["bwa-mem2", "index", "-p", str(iter_index), str(updated_ref)],
                "Building BWA-MEM2 index",
            )
            awk_filter = f"awk '/^@/{{print;next}} $3==\"*\"{{print;next}} {{for(i=12;i<=NF;i++)if($i~/^AS:i:/){{split($i,a,\":\");if(a[3]>={bwa_min_score})print;next}}}}'"
            if is_paired:
                self._run_spinner(
                    f"bwa-mem2 mem -A 2 -B 8 -T {bwa_min_score} -t {threads} '{iter_index}' '{none_r1_path}' '{none_r2_path}' | {awk_filter} > '{iter_sam}'",
                    "Aligning reads with BWA-MEM2",
                    shell=True,
                )
            else:
                self._run_spinner(
                    f"bwa-mem2 mem -A 2 -B 8 -T {bwa_min_score} -t {threads} '{iter_index}' '{none_r1_path}' | {awk_filter} > '{iter_sam}'",
                    "Aligning reads with BWA-MEM2",
                    shell=True,
                )
        elif aligner == "minimap2":
            iter_index = iter_output_dir / "ref_index.mmi"
            preset_map = {
                "short": "sr",
                "pacbio-hifi": "map-hifi",
                "pacbio-clr": "map-pb",
                "ont-q20": "lr:hq",
                "ont-standard": "map-ont",
            }
            preset = preset_map.get(minimap2_profile, "sr")
            score_threshold = f"-s {minimap2_min_score}" if preset == "sr" else ""
            self._run_spinner(
                ["minimap2", "-x", preset, "-d", str(iter_index), str(updated_ref)],
                f"Building minimap2 ({minimap2_profile}) index",
            )
            if is_paired:
                self._run_spinner(
                    f"minimap2 -ax {preset} --MD {score_threshold} -t {threads} '{iter_index}' '{none_r1_path}' '{none_r2_path}' > '{iter_sam}'",
                    f"Aligning reads with minimap2 ({minimap2_profile})",
                    shell=True,
                )
            else:
                self._run_spinner(
                    f"minimap2 -ax {preset} --MD {score_threshold} -t {threads} '{iter_index}' '{none_r1_path}' > '{iter_sam}'",
                    f"Aligning reads with minimap2 ({minimap2_profile})",
                    shell=True,
                )
        
        # Run ParaDISM on NONE reads
        def _run_paradism_iteration():
            iter_msa_obj, iter_seq_to_aln, iter_gene_names = load_msa(str(iter_msa))
            new_assignments = process_sam_to_dict(str(iter_sam), iter_msa_obj, iter_seq_to_aln, iter_gene_names)
            
            # Write outputs for this iteration
            iter_genes = self._write_fastq_outputs(new_assignments, str(none_r1_path), str(none_r2_path) if none_r2_path else None, iter_fastq_dir)
            if iter_genes:
                create_bam_files(iter_genes, str(updated_ref), str(iter_fastq_dir), str(iter_bam_dir), aligner, threads, minimap2_profile, self.prefix, bowtie2_score_min, bwa_min_score, minimap2_min_score)
            return new_assignments
        
        new_assignments = self._run_spinner(_run_paradism_iteration, "Running ParaDISM algorithm")
        
        # 5. Merge assignments
        merged_assignments = self._merge_assignments(previous_assignments, new_assignments)
        
        # 6. Check convergence
        new_none_reads = self._extract_none_reads_from_assignments(merged_assignments)
        reads_rescued = len(none_reads) - len(new_none_reads)
        
        if reads_rescued == 0:
            print(f"  No reads were rescued from NONE. Pipeline has converged.", file=sys.stderr)
            return updated_ref, merged_assignments, True
        
        return updated_ref, merged_assignments, False

    def _write_fastq_outputs(self, assignments: dict[str, str], r1_path: str, r2_path: str | None, fastq_dir: Path) -> list[str]:
        """Write FASTQ outputs from assignments dict. Returns list of genes."""
        return write_fastq_outputs(assignments, r1_path, r2_path, str(fastq_dir), self.prefix)

    def run_pipeline(
        self,
        r1: str | Path,
        r2: str | Path | None,
        ref: str | Path,
        aligner: str = "bwa-mem2",
        threads: int = 4,
        sam: str | Path | None = None,
        minimap2_profile: str | None = None,
        show_header: bool = True,
        iterations: int = 1,
        threshold: str | None = None,
    ) -> None:
        """Execute the ParaDISM pipeline with optional iterative refinement.

        Args:
            iterations: Number of ParaDISM runs (1 = no refinement, 2+ = refinement iterations)
            threshold: Alignment score threshold. For bwa-mem2/minimap2: integer (e.g., "240").
                      For bowtie2: score function (e.g., "G,40,40"). Default based on aligner.
        """

        is_paired = r2 is not None
        r1 = str(r1)
        r2 = str(r2) if r2 else None
        ref = str(ref)
        original_ref = Path(ref)
        sam = str(sam) if sam else None

        # Validate minimap2 profile is provided when using minimap2
        if aligner == "minimap2" and not minimap2_profile:
            raise ValueError("--minimap2-profile must be provided when --aligner minimap2")

        # Set default profile for other aligners
        if not minimap2_profile:
            minimap2_profile = "short"

        # Set default thresholds based on aligner if not provided
        if threshold is None:
            bowtie2_score_min = "G,40,40"
            bwa_min_score = 240
            minimap2_min_score = 240
        else:
            if aligner == "bowtie2":
                bowtie2_score_min = threshold
                bwa_min_score = 240
                minimap2_min_score = 240
            else:
                bowtie2_score_min = "G,40,40"
                bwa_min_score = int(threshold)
                minimap2_min_score = int(threshold)

        if show_header:
            mode_text = "paired-end" if is_paired else "single-end"
            print(f"\n  Running in {mode_text} mode\n", file=sys.stderr)
            print("  Running ParaDISM pipeline...\n", file=sys.stderr)

        # 1. Determine output directories and print iteration header
        original_output_dir = self.output_dir
        if iterations > 1:
            # Create iteration_1 subdirectory for multiple iterations
            iter1_output_dir = self.output_dir / "iteration_1"
            iter1_output_dir.mkdir(exist_ok=True)
            msa_output_dir = iter1_output_dir
            sam_output_dir = iter1_output_dir
            if show_header:
                print(f"\n  \033[0;33m=== Iteration 1 ===\033[0m\n", file=sys.stderr)
            self.logger.section("Iteration 1")
        else:
            msa_output_dir = self.output_dir
            sam_output_dir = self.output_dir
        
        # 2. Create MSA
        msa_output = msa_output_dir / "ref_seq_msa.aln"
        self._run_spinner(
            f"mafft --auto '{ref}' > '{msa_output}'",
            "Creating Multiple Sequence Alignment",
            shell=True,
        )

        # 3. Align reads
        
        sam_output = sam_output_dir / "mapped_reads.sam"

        if sam:
            self._run_spinner(
                ["cp", sam, str(sam_output)],
                "Copying provided SAM file",
            )
        else:
            if aligner == "bowtie2":
                index_base = self.output_dir / "ref_index"
                self._run_spinner(
                    ["bowtie2-build", ref, str(index_base)],
                    "Building Bowtie2 index",
                )
                if is_paired:
                    self._run_spinner(
                        f"bowtie2 --local --score-min {bowtie2_score_min} -p {threads} -x '{index_base}' -1 '{r1}' -2 '{r2}' -S '{sam_output}'",
                        "Aligning reads with Bowtie2",
                        shell=True,
                    )
                else:
                    self._run_spinner(
                        f"bowtie2 --local --score-min {bowtie2_score_min} -p {threads} -x '{index_base}' -U '{r1}' -S '{sam_output}'",
                        "Aligning reads with Bowtie2",
                        shell=True,
                    )
            elif aligner == "bwa-mem2":
                index_base = self.output_dir / "ref_index"
                self._run_spinner(
                    ["bwa-mem2", "index", "-p", str(index_base), ref],
                    "Building BWA-MEM2 index",
                )
                awk_filter = f"awk '/^@/{{print;next}} $3==\"*\"{{print;next}} {{for(i=12;i<=NF;i++)if($i~/^AS:i:/){{split($i,a,\":\");if(a[3]>={bwa_min_score})print;next}}}}'"
                if is_paired:
                    self._run_spinner(
                        f"bwa-mem2 mem -A 2 -B 8 -T {bwa_min_score} -t {threads} '{index_base}' '{r1}' '{r2}' | {awk_filter} > '{sam_output}'",
                        "Aligning reads with BWA-MEM2",
                        shell=True,
                    )
                else:
                    self._run_spinner(
                        f"bwa-mem2 mem -A 2 -B 8 -T {bwa_min_score} -t {threads} '{index_base}' '{r1}' | {awk_filter} > '{sam_output}'",
                        "Aligning reads with BWA-MEM2",
                        shell=True,
                    )
            elif aligner == "minimap2":
                index_file = self.output_dir / "ref_index.mmi"
                preset_map = {
                    "short": "sr",
                    "pacbio-hifi": "map-hifi",
                    "pacbio-clr": "map-pb",
                    "ont-q20": "lr:hq",
                    "ont-standard": "map-ont",
                }
                preset = preset_map.get(minimap2_profile, "sr")
                score_threshold = f"-s {minimap2_min_score}" if preset == "sr" else ""

                self._run_spinner(
                    ["minimap2", "-x", preset, "-d", str(index_file), ref],
                    f"Building minimap2 ({minimap2_profile}) index",
                )
                if is_paired:
                    self._run_spinner(
                        f"minimap2 -ax {preset} --MD {score_threshold} -t {threads} '{index_file}' '{r1}' '{r2}' > '{sam_output}'",
                        f"Aligning reads with minimap2 ({minimap2_profile})",
                        shell=True,
                    )
                else:
                    self._run_spinner(
                        f"minimap2 -ax {preset} --MD {score_threshold} -t {threads} '{index_file}' '{r1}' > '{sam_output}'",
                        f"Aligning reads with minimap2 ({minimap2_profile})",
                        shell=True,
                    )

        # 3. Run initial ParaDISM algorithm
        if iterations > 1:
            # Switch to iteration_1 directory for processing
            self.output_dir = iter1_output_dir
        
        fastq_dir = self.output_dir / f"{self.prefix}_fastq"
        bam_dir = self.output_dir / f"{self.prefix}_bam"

        # Run ParaDISM algorithm directly (not via subprocess) to get assignments dict
        def _run_paradism():
            msa_obj, seq_to_aln, gene_names = load_msa(str(msa_output))
            assignments = process_sam_to_dict(str(sam_output), msa_obj, seq_to_aln, gene_names)
            genes = write_fastq_outputs(assignments, r1, r2, str(fastq_dir), self.prefix)
            if genes:
                create_bam_files(genes, ref, str(fastq_dir), str(bam_dir), aligner, threads, minimap2_profile, self.prefix, bowtie2_score_min, bwa_min_score, minimap2_min_score)
            return assignments

        current_assignments = self._run_spinner(
            _run_paradism,
            "Running ParaDISM algorithm",
        )
        current_ref = Path(ref)
        current_bam_dir = bam_dir

        # Store iteration outputs (only if we'll have multiple iterations)
        iteration_outputs = []
        if iterations > 1:
            iteration_outputs.append({
                'iteration': 1,
                'reference': current_ref,
                'output_dir': self.output_dir,
                'assignments': current_assignments,
                'bam_dir': current_bam_dir,
            })
            self.output_dir = original_output_dir

        # Iterative refinement
        refinement_iterations = iterations - 1
        if refinement_iterations > 0:
            for iteration in range(2, iterations + 1):
                current_ref, current_assignments, converged = self._run_single_iteration_refinement(
                    r1=r1,
                    r2=r2,
                    previous_ref=iteration_outputs[-1]['reference'],
                    previous_bam_dir=iteration_outputs[-1]['bam_dir'],
                    previous_assignments=iteration_outputs[-1]['assignments'],
                    aligner=aligner,
                    threads=threads,
                    minimap2_profile=minimap2_profile,
                    is_paired=is_paired,
                    iteration=iteration,
                    bowtie2_score_min=bowtie2_score_min,
                    bwa_min_score=bwa_min_score,
                    minimap2_min_score=minimap2_min_score,
                )

                if converged:
                    print(f"\n  \033[0;34m✓ Pipeline converged.\033[0m\n", file=sys.stderr)
                    break
                
                iter_output_dir = self.output_dir / f"iteration_{iteration}"
                iter_bam_dir = iter_output_dir / f"{self.prefix}_bam"
                
                iteration_outputs.append({
                    'iteration': iteration,
                    'reference': current_ref,
                    'output_dir': iter_output_dir,
                    'assignments': current_assignments,
                    'bam_dir': iter_bam_dir,
                })

            # Write final outputs from converged/merged assignments
            final_output = iteration_outputs[-1]
            final_outputs_dir = original_output_dir / "final_outputs"
            final_fastq_dir = final_outputs_dir / f"{self.prefix}_fastq"
            final_bam_dir = final_outputs_dir / f"{self.prefix}_bam"
            
            def _write_final_outputs():
                final_genes = self._write_fastq_outputs(
                    final_output['assignments'],
                    r1,
                    r2,
                    str(final_fastq_dir)
                )
                if final_genes:
                    create_bam_files(
                        final_genes,
                        str(original_ref),
                        str(final_fastq_dir),
                        str(final_bam_dir),
                        aligner,
                        threads,
                        minimap2_profile,
                        self.prefix,
                        bowtie2_score_min,
                        bwa_min_score,
                        minimap2_min_score
                    )

            print(f"\n  Iterative refinement complete.\n", file=sys.stderr)
            
            self._run_spinner(_write_final_outputs, "Writing final outputs")
        
        # Write final outputs for single iteration case
        if iterations == 1:
            final_outputs_dir = original_output_dir / "final_outputs"
            final_fastq_dir = final_outputs_dir / f"{self.prefix}_fastq"
            final_bam_dir = final_outputs_dir / f"{self.prefix}_bam"
            
            def _write_final_outputs():
                final_genes = self._write_fastq_outputs(
                    current_assignments,
                    r1,
                    r2,
                    str(final_fastq_dir)
                )
                if final_genes:
                    create_bam_files(
                        final_genes,
                        str(current_ref),
                        str(final_fastq_dir),
                        str(final_bam_dir),
                        aligner,
                        threads,
                        minimap2_profile,
                        self.prefix,
                        bowtie2_score_min,
                        bwa_min_score,
                        minimap2_min_score
                    )

            self._run_spinner(_write_final_outputs, "Writing final outputs")

        # 4. Cleanup intermediate files
        time.sleep(0.2)
        self.logger.section("Cleaning up intermediate files")
        
        # Cleanup top-level intermediate files
        for pattern in [
            "ref_index.*",
            "ref_seq_msa.aln",
        ]:
            for file_path in original_output_dir.glob(pattern):
                file_path.unlink()
        
        # Cleanup intermediate files from iteration directories
        for iteration_dir in original_output_dir.glob("iteration_*"):
            if iteration_dir.is_dir():
                for pattern in [
                    "ref_index.*",
                    "ref_seq_msa.aln",
                    "none_reads_r1.fq",
                    "none_reads_r2.fq",
                ]:
                    for file_path in iteration_dir.glob(pattern):
                        file_path.unlink()

        print("  \033[0;36m✓ Cleaning up intermediate files\033[0m", file=sys.stderr)
        final_outputs_path = original_output_dir / "final_outputs"
        print(f"\n  \033[0;36m✓ Pipeline complete. Final outputs in: {final_outputs_path}\033[0m", file=sys.stderr)
