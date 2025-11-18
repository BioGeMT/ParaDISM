#!/usr/bin/env python3
"""Pipeline execution utilities for the homologous-region mapper."""

from __future__ import annotations

import sys
import time
from pathlib import Path
from typing import Sequence

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

    # ------------------------------------------------------------------ #
    # Helpers
    # ------------------------------------------------------------------ #
    def _run_spinner(self, command: Sequence[str] | str, message: str, *, shell: bool = False) -> None:
        self.progress.run_with_spinner(command, message, shell=shell)

    def _run_progress(self, command: Sequence[str], message: str) -> None:
        self.progress.run_with_progress(command, message)

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
    ) -> None:
        """Execute the complete pipeline (supports both paired-end and single-end)."""

        if show_header:
            print("\n\033[0;36mRunning mapper pipeline...\033[0m\n", file=sys.stderr)

        r1 = str(r1)
        r2 = str(r2) if r2 else None
        ref = str(ref)
        sam = str(sam) if sam else None
        is_paired = r2 is not None

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
                        f"bowtie2 -p {threads} -x '{index_base}' -1 '{r1}' -2 '{r2}' -S '{sam_output}'",
                        "Aligning reads with Bowtie2",
                        shell=True,
                    )
                else:
                    self._run_spinner(
                        f"bowtie2 -p {threads} -x '{index_base}' -U '{r1}' -S '{sam_output}'",
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
                        f"bwa-mem2 mem -t {threads} '{index_base}' '{r1}' '{r2}' > '{sam_output}'",
                        "Aligning reads with BWA-MEM2",
                        shell=True,
                    )
                else:
                    self._run_spinner(
                        f"bwa-mem2 mem -t {threads} '{index_base}' '{r1}' > '{sam_output}'",
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

                self._run_spinner(
                    ["minimap2", "-x", preset, "-d", str(index_file), ref],
                    f"Building minimap2 ({minimap2_profile}) index",
                )
                if is_paired:
                    self._run_spinner(
                        f"minimap2 -ax {preset} --MD -t {threads} '{index_file}' '{r1}' '{r2}' > '{sam_output}'",
                        f"Aligning reads with minimap2 ({minimap2_profile})",
                        shell=True,
                    )
                else:
                    self._run_spinner(
                        f"minimap2 -ax {preset} --MD -t {threads} '{index_file}' '{r1}' > '{sam_output}'",
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
        print(f"\nPipeline complete. Outputs saved to: {self.output_dir}\n", file=sys.stderr)
