#!/usr/bin/env python3
"""Pipeline execution utilities for the homologous-region mapper."""

from __future__ import annotations

import sys
import time
from pathlib import Path
from typing import Sequence
from Bio import SeqIO

from utils.logger import PipelineLogger
from utils.progress import ProgressRunner

PIPELINE_DIR = Path(__file__).resolve().parent


class PipelineExecutor:
    """Executes the homologous-region mapper pipeline steps."""

    def __init__(self, output_dir: str | Path = "./output", pipeline_dir: str | Path | None = None, prefix: str | None = None) -> None:
        self.output_dir = Path(output_dir)
        self.pipeline_dir = Path(pipeline_dir) if pipeline_dir else PIPELINE_DIR
        self.output_dir.mkdir(exist_ok=True)

        # Use provided prefix, or extract from output directory name
        if prefix is not None:
            self.prefix = prefix
        else:
            self.prefix = self.output_dir.name if self.output_dir.name != "." else "output"

        timestamp = time.strftime('%Y%m%d_%H%M%S')
        self.log_file = self.output_dir / f"{self.prefix}_pipeline_{timestamp}.log"
        self.logger = PipelineLogger(self.log_file)
        self.progress = ProgressRunner(self.logger)
        
        # Track iteration outputs for analysis
        self.iteration_outputs = []

    # ------------------------------------------------------------------ #
    # Helpers
    # ------------------------------------------------------------------ #
    def _run_spinner(self, command: Sequence[str] | str, message: str, *, shell: bool = False) -> None:
        self.progress.run_with_spinner(command, message, shell=shell)

    def _run_progress(self, command: Sequence[str], message: str) -> None:
        self.progress.run_with_progress(command, message)

    def _extract_none_reads_from_tsv(self, tsv_path: Path) -> set[str]:
        """Extract read IDs that mapped to NONE from a unique_mappings.tsv file."""
        none_reads = set()
        with open(tsv_path, 'r') as f:
            next(f)  # Skip header
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = line.split('\t')
                read_name = parts[0]
                gene = parts[1] if len(parts) > 1 else "NONE"
                # Remove strand suffixes (+/-) if present
                read_name = read_name.rstrip('+-')
                # Remove (plus)/(minus) tags from gene if present
                if " (plus)" in gene or " (minus)" in gene:
                    gene = gene.replace(" (plus)", "").replace(" (minus)", "")
                if gene == "NONE":
                    none_reads.add(read_name)
        return none_reads
    
    def _extract_reads_from_fastq(self, fastq_path: Path, read_ids: set[str], output_path: Path, is_paired: bool = False) -> int:
        """Extract reads from FASTQ file based on read IDs. Returns count of extracted reads."""
        output_path.parent.mkdir(parents=True, exist_ok=True)
        count = 0
        with open(output_path, 'w') as outfile:
            for rec in SeqIO.parse(fastq_path, "fastq"):
                # Normalize read ID (remove /1 or /2 suffix for matching)
                base_id = rec.id
                if base_id.endswith("/1") or base_id.endswith("/2"):
                    base_id = base_id[:-2]
                
                # Check if this read should be extracted
                if base_id in read_ids or rec.id in read_ids:
                    SeqIO.write(rec, outfile, "fastq")
                    count += 1
        return count
    
    def _merge_tsv_files(self, previous_tsv: Path, new_tsv: Path, output_tsv: Path) -> None:
        """Merge two TSV files, keeping non-NONE mappings from previous and new mappings from NONE reads."""
        previous_mappings = {}
        previous_was_none = set()
        
        # Read previous TSV - keep all mappings, track which were NONE
        with open(previous_tsv, 'r') as f:
            header = next(f).strip()
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = line.split('\t')
                read_name = parts[0]
                gene = parts[1] if len(parts) > 1 else "NONE"
                # Remove strand suffixes
                read_name = read_name.rstrip('+-')
                # Remove strand tags from gene for comparison
                gene_clean = gene.replace(" (plus)", "").replace(" (minus)", "")
                
                # Store all mappings
                previous_mappings[read_name] = line
                if gene_clean == "NONE":
                    previous_was_none.add(read_name)
        
        # Read new TSV - these are mappings from re-aligned NONE reads
        new_mappings = {}
        with open(new_tsv, 'r') as f:
            next(f)  # Skip header
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = line.split('\t')
                read_name = parts[0]
                gene = parts[1] if len(parts) > 1 else "NONE"
                # Remove strand suffixes
                read_name = read_name.rstrip('+-')
                new_mappings[read_name] = line
        
        # Merge strategy:
        # 1. Keep all previous non-NONE mappings (they were successful, don't change them)
        # 2. For reads that were NONE in previous iteration, use new mapping (even if still NONE)
        merged_mappings = {}
        
        # First, add all previous non-NONE mappings
        for read_name, line in previous_mappings.items():
            if read_name not in previous_was_none:
                merged_mappings[read_name] = line
        
        # Then, add new mappings for reads that were NONE
        for read_name, new_line in new_mappings.items():
            if read_name in previous_was_none:
                # This read was NONE in previous iteration, use new mapping
                merged_mappings[read_name] = new_line
            elif read_name not in merged_mappings:
                # New read not in previous (shouldn't happen, but handle it)
                merged_mappings[read_name] = new_line
        
        # Write merged TSV
        output_tsv.parent.mkdir(parents=True, exist_ok=True)
        with open(output_tsv, 'w') as f:
            f.write(header + '\n')
            for read_name in sorted(merged_mappings.keys()):
                f.write(merged_mappings[read_name] + '\n')
    
    def _run_single_iteration_refinement(
        self,
        r1: str,
        r2: str | None,
        previous_ref: Path,
        previous_bam_dir: Path,
        previous_mappings_tsv: Path,
        aligner: str,
        threads: int,
        minimap2_profile: str,
        is_paired: bool,
        iteration: int,
    ) -> tuple[Path, Path]:
        """
        Run one iteration of refinement on NONE reads only:
        1. Identify reads that mapped to NONE in previous iteration
        2. Call variants from successfully mapped reads (non-NONE)
        3. Apply variants to reference
        4. Re-align only the NONE reads with updated reference
        5. Merge results: keep previous successful mappings + new mappings from NONE reads
        
        Returns:
            Tuple of (updated_reference_path, iteration_output_directory)
        """
        print(f"\n\033[0;33m=== Iteration {iteration} (refining NONE reads only) ===\033[0m", file=sys.stderr)
        
        # Create iteration directory with all outputs organized together
        iter_output_dir = self.output_dir / f"iteration_{iteration}"
        iter_output_dir.mkdir(exist_ok=True)
        
        # 1. Identify NONE reads from previous iteration
        print(f"  Identifying NONE reads from previous iteration...", file=sys.stderr)
        none_reads = self._extract_none_reads_from_tsv(previous_mappings_tsv)
        print(f"  Found {len(none_reads)} reads mapped to NONE", file=sys.stderr)
        
        if len(none_reads) == 0:
            print(f"  No NONE reads to refine. Skipping iteration {iteration}.", file=sys.stderr)
            # Copy previous mappings as final output
            final_mappings_tsv = iter_output_dir / f"{self.prefix}_unique_mappings.tsv"
            import shutil
            shutil.copy(previous_mappings_tsv, final_mappings_tsv)
            return previous_ref, iter_output_dir
        
        # 2. Extract NONE reads from original FASTQ files
        none_r1_path = iter_output_dir / "none_reads_r1.fq"
        none_r2_path = iter_output_dir / "none_reads_r2.fq" if is_paired else None
        
        print(f"  Extracting {len(none_reads)} NONE reads from FASTQ files...", file=sys.stderr)
        r1_count = self._extract_reads_from_fastq(Path(r1), none_reads, none_r1_path, is_paired)
        r2_count = 0
        if is_paired and r2:
            r2_count = self._extract_reads_from_fastq(Path(r2), none_reads, none_r2_path, is_paired)
        
        print(f"  Extracted {r1_count} R1 reads and {r2_count} R2 reads", file=sys.stderr)
        
        # Check if we actually extracted any reads
        if r1_count == 0:
            print(f"  No reads extracted. Skipping iteration {iteration}.", file=sys.stderr)
            final_mappings_tsv = iter_output_dir / f"{self.prefix}_unique_mappings.tsv"
            import shutil
            shutil.copy(previous_mappings_tsv, final_mappings_tsv)
            return previous_ref, iter_output_dir
        
        # 3. Call variants from successfully mapped reads (non-NONE) and apply to reference
        variant_calling_dir = iter_output_dir / "variant_calling"
        variant_calling_dir.mkdir(exist_ok=True)
        
        per_gene_vcf_dir = variant_calling_dir / "per_gene_vcfs"
        per_gene_vcf_dir.mkdir(exist_ok=True)
        
        vcf_output = variant_calling_dir / "variants.vcf"
        updated_ref = variant_calling_dir / "updated_ref.fa"
        
        self._run_progress(
            [
                "python",
                str(self.pipeline_dir / "variant_refinement.py"),
                "--bam-dir", str(previous_bam_dir),
                "--reference", str(previous_ref),
                "--output-vcf", str(vcf_output),
                "--output-ref", str(updated_ref),
                "--per-gene-vcf-dir", str(per_gene_vcf_dir),
            ],
            f"Calling variants from mapped reads (iteration {iteration})",
        )
        
        # 4. Re-align only NONE reads with updated reference
        print(f"  Re-aligning {len(none_reads)} NONE reads with refined reference...", file=sys.stderr)
        
        # Create temporary executor for this iteration
        iter_executor = PipelineExecutor(
            output_dir=iter_output_dir,
            prefix=self.prefix,
        )
        
        # Run pipeline on NONE reads only
        iter_executor.run_pipeline(
            r1=str(none_r1_path),
            r2=str(none_r2_path) if none_r2_path else None,
            ref=str(updated_ref),
            aligner=aligner,
            threads=threads,
            sam=None,  # Always realign with new reference
            minimap2_profile=minimap2_profile,
            show_header=False,
            iterations=0,  # No nested iterations
        )
        
        # 5. Merge results: keep previous successful mappings + new mappings from NONE reads
        new_mappings_tsv = iter_output_dir / f"{self.prefix}_unique_mappings.tsv"
        final_mappings_tsv = iter_output_dir / f"{self.prefix}_unique_mappings.tsv"
        
        print(f"  Merging previous successful mappings with new NONE read mappings...", file=sys.stderr)
        self._merge_tsv_files(previous_mappings_tsv, new_mappings_tsv, final_mappings_tsv)
        
        # Update BAM directory to include merged results (we'll need to regenerate BAMs)
        # For now, we'll keep the new BAM directory, but ideally we'd merge BAMs too
        # This is a limitation - BAM merging would require more complex logic
        
        return updated_ref, iter_output_dir

    # ------------------------------------------------------------------ #
    # Public API
    # ------------------------------------------------------------------ #
    def run_pipeline(
        self,
        r1: str | Path,
        r2: str | Path | None,
        ref: str | Path,
        aligner: str = "bwa-mem2",
        threads: int = 4,
        sam: str | Path | None = None,
        minimap2_profile: str = "short",
        show_header: bool = True,
        iterations: int = 1,
    ) -> None:
        """Execute the complete pipeline (supports both paired-end and single-end).
        
        Args:
            iterations: Number of ParaDISM runs (1 = no refinement, 2 = 1 refinement iteration, etc.).
                        Each refinement iteration calls variants, updates reference, and re-runs ParaDISM.
        """

        if show_header:
            print("\n\033[0;36mRunning mapper pipeline...\033[0m\n", file=sys.stderr)

        r1 = str(r1)
        r2 = str(r2) if r2 else None
        ref = str(ref)
        sam = str(sam) if sam else None
        is_paired = r2 is not None

        # For iteration 1 (initial run), create iteration_1 directory to keep structure consistent
        # Only create iteration_1 if iterations > 0 (not when called from refinement with iterations=0)
        original_output_dir = self.output_dir
        if iterations > 0:
            iter1_output_dir = self.output_dir / "iteration_1"
            iter1_output_dir.mkdir(exist_ok=True)
            # Temporarily redirect self.output_dir for initial run
            self.output_dir = iter1_output_dir
        else:
            # When iterations=0, use output_dir directly (no subdirectory)
            iter1_output_dir = self.output_dir
        
        msa_output = self.output_dir / "ref_seq_msa.aln"
        self._run_spinner(
            f"mafft --auto '{ref}' > '{msa_output}'",
            "Creating Multiple Sequence Alignment",
            shell=True,
        )

        sam_output = self.output_dir / "mapped_reads.sam"

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
                        f"bowtie2 --local --score-min G,40,40 -p {threads} -x '{index_base}' -1 '{r1}' -2 '{r2}' -S '{sam_output}'",
                        "Aligning reads with Bowtie2",
                        shell=True,
                    )
                else:
                    self._run_spinner(
                        f"bowtie2 --local --score-min G,40,40 -p {threads} -x '{index_base}' -U '{r1}' -S '{sam_output}'",
                        "Aligning reads with Bowtie2",
                        shell=True,
                    )
            elif aligner == "bwa-mem2":
                index_base = self.output_dir / "ref_index"
                self._run_spinner(
                    ["bwa-mem2", "index", "-p", str(index_base), ref],
                    "Building BWA-MEM2 index",
                )
                if is_paired:
                    self._run_spinner(
                        f"bwa-mem2 mem -A 2 -B 8 -T 240 -t {threads} '{index_base}' '{r1}' '{r2}' > '{sam_output}'",
                        "Aligning reads with BWA-MEM2",
                        shell=True,
                    )
                else:
                    self._run_spinner(
                        f"bwa-mem2 mem -A 2 -B 8 -T 240 -t {threads} '{index_base}' '{r1}' > '{sam_output}'",
                        "Aligning reads with BWA-MEM2",
                        shell=True,
                    )
            elif aligner == "minimap2":
                index_file = self.output_dir / "ref_index.mmi"
                # Modern minimap2 presets (2024-2025)
                preset_map = {
                    "short": "sr",
                    "pacbio-hifi": "map-hifi",
                    "pacbio-clr": "map-pb",
                    "ont-q20": "lr:hq",
                    "ont-standard": "map-ont",
                }
                preset = preset_map.get(minimap2_profile, "sr")

                # Add stringent score threshold only for short reads
                score_threshold = "-s 240" if preset == "sr" else ""

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

        ref_msa_tsv = self.output_dir / "ref_seq_msa.tsv"
        self._run_spinner(
            [
                "python",
                str(self.pipeline_dir / "ref_2_msa.py"),
                "--reference_fasta",
                ref,
                "--msa_file",
                str(msa_output),
                "--output_file",
                str(ref_msa_tsv),
            ],
            "Mapping reference sequences to MSA",
        )

        mapped_reads_tsv = self.output_dir / "mapped_reads.tsv"
        read_2_gene_cmd = [
            "python",
            str(self.pipeline_dir / "read_2_gene.py"),
            "--sam",
            str(sam_output),
            "--output",
            str(mapped_reads_tsv),
            "--fastq",
            r1,
        ]
        if not is_paired:
            read_2_gene_cmd.append("--single-end")

        self._run_progress(
            read_2_gene_cmd,
            "Mapping reads to reference sequences",
        )

        unique_mappings_tsv = self.output_dir / f"{self.prefix}_unique_mappings.tsv"
        self._run_progress(
            [
                "python",
                str(self.pipeline_dir / "mapper_algo_snp_only.py"),
                "--read_map",
                str(mapped_reads_tsv),
                "--msa",
                str(ref_msa_tsv),
                "--output_file",
                str(unique_mappings_tsv),
                "--batch_size",
                "10000",
            ]
            + ([] if is_paired else ["--single-end"]),
            "Refining mapping",
        )

        fastq_dir = self.output_dir / f"{self.prefix}_fastq"
        bam_dir = self.output_dir / f"{self.prefix}_bam"

        # Build command with optional R2
        output_cmd = [
            "python",
            str(self.pipeline_dir / "output.py"),
            "--tsv",
            str(unique_mappings_tsv),
            "--r1",
            r1,
        ]
        if r2:
            output_cmd.extend(["--r2", r2])
        output_cmd.extend([
            "--ref",
            ref,
            "--fastq-dir",
            str(fastq_dir),
            "--bam-dir",
            str(bam_dir),
            "--aligner",
            aligner,
            "--threads",
            str(threads),
            "--minimap2-profile",
            minimap2_profile,
            "--prefix",
            self.prefix,
        ])

        self._run_progress(
            output_cmd,
            "Writing output files",
        )

        time.sleep(0.2)

        self.logger.section("Cleaning up intermediate files")
        for pattern in [
            "ref_index.*",
            "ref_seq_msa.aln",
            "ref_seq_msa.tsv",
            "mapped_reads.tsv",
            # "mapped_reads.sam",  # Keep SAM file
        ]:
            for file_path in self.output_dir.glob(pattern):
                file_path.unlink()

        print("\033[0;36m✓\033[0m \033[0;36mCleaning up intermediate files\033[0m", file=sys.stderr)
        
        # Store iteration 1 outputs (initial run)
        # Only store if this was iteration 1 (not when called from refinement)
        if iterations > 0:
            self.iteration_outputs.append({
                'iteration': 1,
                'reference': Path(ref),
                'output_dir': self.output_dir,
                'mappings_tsv': unique_mappings_tsv,
                'bam_dir': bam_dir,
            })
            
            # Restore original output_dir for subsequent operations
            self.output_dir = original_output_dir
        else:
            # When iterations=0 (called from refinement), don't store in iteration_outputs
            # The outputs are already in the correct directory
            pass
        
        # Iterative refinement
        # iterations=1 means run once (no refinement)
        # iterations=2 means run twice (1 refinement iteration)
        # So refinement_iterations = iterations - 1
        refinement_iterations = iterations - 1
        if refinement_iterations > 0:
            print(f"\n\033[0;36mStarting iterative refinement ({refinement_iterations} iteration(s))...\033[0m\n", file=sys.stderr)
            for iteration in range(2, iterations + 1):
                iter_ref, iter_output_dir = self._run_single_iteration_refinement(
                    r1=r1,
                    r2=r2,
                    previous_ref=self.iteration_outputs[-1]['reference'],
                    previous_bam_dir=self.iteration_outputs[-1]['bam_dir'],
                    previous_mappings_tsv=self.iteration_outputs[-1]['mappings_tsv'],
                    aligner=aligner,
                    threads=threads,
                    minimap2_profile=minimap2_profile,
                    is_paired=is_paired,
                    iteration=iteration,
                )
                
                iter_mappings_tsv = iter_output_dir / f"{self.prefix}_unique_mappings.tsv"
                iter_bam_dir = iter_output_dir / f"{self.prefix}_bam"
                
                self.iteration_outputs.append({
                    'iteration': iteration,
                    'reference': iter_ref,
                    'output_dir': iter_output_dir,
                    'mappings_tsv': iter_mappings_tsv,
                    'bam_dir': iter_bam_dir,
                })
            
            # Update output_dir to final iteration
            final_output = self.iteration_outputs[-1]
            self.output_dir = final_output['output_dir']
            print(f"\n\033[0;36mIterative refinement complete. Final outputs from iteration {iterations}.\033[0m\n", file=sys.stderr)
        
        print(f"\nPipeline complete. Outputs saved to: {self.output_dir}\n", file=sys.stderr)
