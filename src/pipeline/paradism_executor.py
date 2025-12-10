#!/usr/bin/env python3
"""Streamlined ParaDISM executor that writes outputs directly without intermediate TSV files.

This executor uses paradism_algo.py's process_sam_to_dict() function to get assignments
as a dictionary, then writes FASTQ and BAM outputs directly without creating intermediate
TSV files. A TSV is optionally written at the end for reference.

Key differences from executor.py:
- No intermediate TSV step (assignments stored in dict)
- Direct FASTQ writing from assignments dict
- Simpler, more efficient workflow
"""

from __future__ import annotations

import sys
import time
import shutil
import threading
import subprocess
from pathlib import Path
from typing import Sequence
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from utils.logger import PipelineLogger
from utils.progress import ProgressRunner
from paradism_algo import load_msa, process_sam_to_dict
from output import create_bam_files

PIPELINE_DIR = Path(__file__).resolve().parent


class ParaDISMExecutor:
    """Streamlined executor for ParaDISM algorithm that writes outputs directly."""

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

    def _run_step_with_spinner(self, message: str, func, *args, **kwargs):
        """Run a Python callable with a spinner and plain success mark."""
        self.logger.section(message)

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
            result = func(*args, **kwargs)
            stop_event.set()
            spinner_thread.join()
            print(f"\r  \033[0;36m✓ {message}\033[0m", file=sys.stderr)
            return result
        except Exception:
            stop_event.set()
            spinner_thread.join()
            print(f"\r  \033[0;31m✗ {message}\033[0m", file=sys.stderr)
            raise

    def _run_spinner(self, command: Sequence[str] | str, message: str, *, shell: bool = False) -> None:
        self.progress.run_with_spinner(command, message, shell=shell)

    def _run_progress(self, command: Sequence[str], message: str) -> None:
        self.progress.run_with_progress(command, message)

    def _base_read_id(self, read_id: str) -> str:
        """Normalize read ID to match SAM/TSV names."""
        if read_id.endswith("/1") or read_id.endswith("/2"):
            return read_id[:-2]
        return read_id

    def _write_fastq_outputs(self, assignments: dict[str, str], r1_path: str, r2_path: str | None, fastq_dir: Path) -> list[str]:
        """Write per-gene FASTQ files directly from assignments dict."""
        fastq_dir.mkdir(exist_ok=True)
        
        # Load FASTQ records
        r1_records = {}
        for rec in SeqIO.parse(r1_path, "fastq"):
            r1_records[self._base_read_id(rec.id)] = rec

        r2_records = {}
        if r2_path:
            for rec in SeqIO.parse(r2_path, "fastq"):
                r2_records[self._base_read_id(rec.id)] = rec

        # Organize reads by gene
        gene_collection = defaultdict(list)
        is_paired = r2_path is not None

        for read_name, gene in assignments.items():
            if gene == "NONE":
                continue

            # Add R1
            r1 = r1_records.get(read_name)
            if r1 is not None:
                if is_paired:
                    gene_collection[gene].append(SeqRecord(
                        r1.seq,
                        id=f"{read_name}/1",
                        description="",
                        letter_annotations=r1.letter_annotations
                    ))
                else:
                    gene_collection[gene].append(SeqRecord(
                        r1.seq,
                        id=read_name,
                        description="",
                        letter_annotations=r1.letter_annotations
                    ))

            # Add R2 if paired-end
            if is_paired:
                r2 = r2_records.get(read_name)
                if r2 is not None:
                    gene_collection[gene].append(SeqRecord(
                        r2.seq,
                        id=f"{read_name}/2",
                        description="",
                        letter_annotations=r2.letter_annotations
                    ))

        # Write FASTQ files
        processed_genes = []
        for gene, records in gene_collection.items():
            if not records:
                continue

            filename = f"{self.prefix}_{gene}.fq"
            output_file = fastq_dir / filename

            with open(output_file, "w") as out_f:
                SeqIO.write(records, out_f, "fastq")

            print(f"Created {filename} with {len(records)} reads", file=sys.stderr)
            processed_genes.append(gene)

        return processed_genes

    def _create_bam_files(self, genes: list[str], ref_fasta: str, fastq_dir: Path, bam_dir: Path, 
                          aligner: str, threads: int, minimap2_profile: str, is_paired: bool) -> None:
        """Create BAM files for each gene."""
        create_bam_files(
            genes=genes,
            ref_fasta=ref_fasta,
            fastq_dir=str(fastq_dir),
            output_dir=str(bam_dir),
            aligner=aligner,
            threads=threads,
            minimap2_profile=minimap2_profile,
            is_paired=is_paired,
            prefix=self.prefix,
        )

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
    ) -> None:
        """Execute the complete ParaDISM pipeline with direct output writing."""

        is_paired = r2 is not None
        r1 = str(r1)
        r2 = str(r2) if r2 else None
        ref = str(ref)
        sam = str(sam) if sam else None

        if show_header:
            mode_text = "paired-end" if is_paired else "single-end"
            print(f"\n  Running in {mode_text} mode\n", file=sys.stderr)
            print("  Running ParaDISM pipeline...\n", file=sys.stderr)

        # 1. Create MSA
        msa_output = self.output_dir / "ref_seq_msa.aln"
        self._run_spinner(
            f"mafft --auto '{ref}' > '{msa_output}'",
            "Creating Multiple Sequence Alignment",
            shell=True,
        )

        # 2. Align reads
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
                bwa_min_score = 240
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

        # 3. Load MSA and process SAM to get assignments dict
        assignments = self._run_step_with_spinner(
            "Running ParaDISM algorithm",
            self._process_assignments,
            str(sam_output),
            str(msa_output),
        )

        # 4. Write outputs directly from dict
        fastq_dir = self.output_dir / f"{self.prefix}_fastq"
        bam_dir = self.output_dir / f"{self.prefix}_bam"

        genes = self._run_step_with_spinner(
            "Writing FASTQ files",
            self._write_fastq_outputs,
            assignments,
            r1,
            r2,
            fastq_dir,
        )

        if genes:
            self._run_step_with_spinner(
                "Creating BAM files",
                self._create_bam_files,
                genes,
                ref,
                fastq_dir,
                bam_dir,
                aligner,
                threads,
                minimap2_profile,
                is_paired,
            )

        # 5. Cleanup
        time.sleep(0.2)
        self.logger.section("Cleaning up intermediate files")
        for pattern in [
            "ref_index.*",
            "ref_seq_msa.aln",
        ]:
            for file_path in self.output_dir.glob(pattern):
                file_path.unlink()

        print("  \033[0;36m✓ Cleaning up intermediate files\033[0m", file=sys.stderr)
        print(f"\n  \033[0;36m✓ Pipeline complete. Outputs in: {self.output_dir}\033[0m", file=sys.stderr)

    def _process_assignments(self, sam_path: str, msa_path: str) -> dict[str, str]:
        """Process SAM file and return assignments dict."""
        msa, seq_to_aln, gene_names = load_msa(msa_path)
        assignments = process_sam_to_dict(sam_path, msa, seq_to_aln, gene_names)
        return assignments
