#!/usr/bin/env python3
"""Pipeline execution utilities for mapper2 (per-gene alignment mode)."""

from __future__ import annotations

import concurrent.futures
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Sequence

import subprocess
from pipeline.mapper_algo_snp_only import load_msa_mapping
from utils.logger import PipelineLogger
from utils.progress import ProgressRunner

PIPELINE_DIR = Path(__file__).resolve().parent


class PipelineExecutor:
    """Executes the mapper2 pipeline steps."""

    def __init__(
        self,
        output_dir: str | Path = "./output",
        pipeline_dir: str | Path | None = None,
        prefix: str | None = None,
    ) -> None:
        self.output_dir = Path(output_dir)
        self.pipeline_dir = Path(pipeline_dir) if pipeline_dir else PIPELINE_DIR
        self.output_dir.mkdir(parents=True, exist_ok=True)

        if prefix is not None:
            self.prefix = prefix
        else:
            self.prefix = self.output_dir.name if self.output_dir.name != "." else "output"

        timestamp = time.strftime('%Y%m%d_%H%M%S')
        self.log_file = self.output_dir / f"{self.prefix}_pipeline_{timestamp}.log"
        self.logger = PipelineLogger(self.log_file)
        self.progress = ProgressRunner(self.logger)
        self.per_gene_dir = self.output_dir / "per_gene"
        self.per_gene_dir.mkdir(exist_ok=True)

    # ------------------------------------------------------------------ #
    # Helpers
    # ------------------------------------------------------------------ #
    def _run_spinner(self, command: Sequence[str] | str, message: str, *, shell: bool = False) -> None:
        self.progress.run_with_spinner(command, message, shell=shell)

    def _run_progress(self, command: Sequence[str], message: str) -> None:
        self.progress.run_with_progress(command, message)

    def _run_command(self, cmd: Sequence[str], message: str, gene: str | None = None, stdout=None) -> None:
        prefix = f"[{gene}] " if gene else ""
        print(f"{prefix}{message}", file=sys.stderr)
        run_kwargs: Dict = {"stderr": subprocess.PIPE, "text": True, "check": True}
        if stdout is None:
            run_kwargs["stdout"] = subprocess.PIPE
        else:
            run_kwargs["stdout"] = stdout
        try:
            result = subprocess.run(cmd, **run_kwargs)
        except subprocess.CalledProcessError as exc:
            stdout_data = getattr(exc, "stdout", None)
            stderr_data = getattr(exc, "stderr", None)
            if stdout_data:
                self.logger.write(stdout_data)
            if stderr_data:
                self.logger.write(stderr_data)
            raise
        else:
            if stdout is None and result.stdout:
                self.logger.write(result.stdout)
            if result.stderr:
                self.logger.write(result.stderr)

    def _write_slice_file(
        self,
        read_map_path: Path,
        gene_name: str,
        sequence_positions: Dict[str, Dict[int, int]],
        output_path: Path,
    ) -> None:
        mapping = sequence_positions.get(gene_name, {})
        per_read: Dict[str, int] = {}
        with open(read_map_path, "r") as handle:
            next(handle)
            for line in handle:
                parts = line.rstrip().split("\t")
                if len(parts) < 7:
                    continue
                read_name = parts[0].rstrip("+-")
                ref_pos = parts[2]
                ref_gene = parts[6]
                if ref_gene != gene_name or ref_pos == "-":
                    continue
                if read_name in per_read:
                    continue
                try:
                    msa_pos = mapping.get(int(ref_pos))
                except ValueError:
                    continue
                if msa_pos is None:
                    continue
                per_read[read_name] = msa_pos
        with open(output_path, "w") as out:
            out.write("Read_Name\tMSA_Start\n")
            for read, start in per_read.items():
                out.write(f"{read}\t{start}\n")

    def _process_gene(
        self,
        gene_name: str,
        gene_sequence: str,
        r1: str,
        r2: str | None,
        aligner: str,
        threads: int,
        minimap2_profile: str,
        is_paired: bool,
        sequence_positions: Dict[str, Dict[int, int]],
        msa_tsv: Path,
    ) -> Dict[str, Path]:
        gene_dir = self.per_gene_dir / gene_name
        gene_dir.mkdir(parents=True, exist_ok=True)

        gene_fasta = gene_dir / f"{gene_name}.fa"
        with open(gene_fasta, "w") as handle:
            handle.write(f">{gene_name}\n{gene_sequence}\n")

        sam_output = gene_dir / f"{gene_name}.sam"
        mapped_reads_tsv = gene_dir / f"{gene_name}_mapped_reads.tsv"
        unique_mappings = gene_dir / f"{gene_name}_unique_mappings.tsv"
        slice_file = gene_dir / f"{gene_name}_slice_starts.tsv"

        if aligner == "bowtie2":
            index_base = gene_dir / f"{gene_name}_index"
            self._run_command(["bowtie2-build", str(gene_fasta), str(index_base)], "Building Bowtie2 index", gene=gene_name)
            if is_paired:
                if not r2:
                    raise ValueError("Paired-end mode requires both R1 and R2")
                cmd = [
                    "bowtie2",
                    "-p",
                    str(threads),
                    "-x",
                    str(index_base),
                    "-1",
                    r1,
                    "-2",
                    r2,
                    "-S",
                    str(sam_output),
                ]
            else:
                cmd = ["bowtie2", "-p", str(threads), "-x", str(index_base), "-U", r1, "-S", str(sam_output)]
            self._run_command(cmd, "Aligning reads with Bowtie2", gene=gene_name)
        elif aligner == "bwa-mem2":
            index_base = gene_dir / f"{gene_name}_index"
            self._run_command(["bwa-mem2", "index", "-p", str(index_base), str(gene_fasta)], "Building BWA-MEM2 index", gene=gene_name)
            cmd = ["bwa-mem2", "mem", "-t", str(threads), str(index_base), r1]
            if is_paired:
                if not r2:
                    raise ValueError("Paired-end mode requires both R1 and R2")
                cmd.append(r2)
            with open(sam_output, "w") as sam_handle:
                self._run_command(cmd, "Aligning reads with BWA-MEM2", gene=gene_name, stdout=sam_handle)
        elif aligner == "minimap2":
            preset_map = {
                "short": "sr",
                "pacbio-hifi": "map-hifi",
                "pacbio-clr": "map-pb",
                "ont-q20": "lr:hq",
                "ont-standard": "map-ont",
            }
            preset = preset_map.get(minimap2_profile, "sr")
            index_file = gene_dir / f"{gene_name}.mmi"
            self._run_command(
                ["minimap2", "-x", preset, "-d", str(index_file), str(gene_fasta)],
                f"Building minimap2 ({minimap2_profile}) index",
                gene=gene_name,
            )
            cmd = ["minimap2", "-ax", preset, "--MD", "-t", str(threads), str(index_file), r1]
            if is_paired:
                if not r2:
                    raise ValueError("Paired-end mode requires both R1 and R2")
                cmd.append(r2)
            with open(sam_output, "w") as sam_handle:
                self._run_command(cmd, f"Aligning reads with minimap2 ({minimap2_profile})", gene=gene_name, stdout=sam_handle)
        else:
            raise ValueError(f"Unsupported aligner: {aligner}")

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
        self._run_command(read_2_gene_cmd, "Converting SAM to TSV", gene=gene_name)

        mapper_cmd = [
            "python",
            str(self.pipeline_dir / "mapper_algo_snp_only.py"),
            "--read_map",
            str(mapped_reads_tsv),
            "--msa",
            str(msa_tsv),
            "--output_file",
            str(unique_mappings),
        ]
        if not is_paired:
            mapper_cmd.append("--single-end")
        self._run_command(mapper_cmd, "Running refinement", gene=gene_name)

        self._write_slice_file(mapped_reads_tsv, gene_name, sequence_positions, slice_file)

        return {
            "gene": gene_name,
            "unique": unique_mappings,
            "read_map": mapped_reads_tsv,
            "slice": slice_file,
        }

    def _reconcile_results(
        self,
        gene_results: List[Dict[str, Path]],
        gene_order: List[str],
    ) -> tuple[Path, Dict[str, int], set[str]]:
        final_unique = self.output_dir / f"{self.prefix}_unique_mappings.tsv"
        per_gene_calls: Dict[str, Dict[str, str]] = defaultdict(dict)
        all_reads: set[str] = set()

        for result in gene_results:
            gene = result["gene"]
            with open(result["unique"], "r") as handle:
                next(handle)
                for line in handle:
                    parts = line.rstrip().split("\t")
                    if len(parts) < 2:
                        continue
                    read, assignment = parts[0], parts[1]
                    per_gene_calls[read][gene] = assignment
                    all_reads.add(read)

        consensus_counts = {"single_gene": 0, "all_none": 0, "conflict": 0}
        with open(final_unique, "w") as out:
            out.write("Read_Name\tUniquely_Mapped\n")
            for read in sorted(all_reads):
                votes = [per_gene_calls[read].get(gene, "NONE") for gene in gene_order]
                non_none = [v for v in votes if v != "NONE"]
                if not non_none:
                    final_call = "NONE"
                    category = "all_none"
                elif len(set(non_none)) == 1:
                    final_call = non_none[0]
                    category = "single_gene"
                else:
                    final_call = "NONE"
                    category = "conflict"
                out.write(f"{read}\t{final_call}\n")
                consensus_counts[category] += 1

        return final_unique, consensus_counts, all_reads

    def _build_slice_lookup(self, gene_results: List[Dict[str, Path]]) -> Dict[str, Dict[str, int]]:
        lookup: Dict[str, Dict[str, int]] = {}
        for result in gene_results:
            gene = result["gene"]
            lookup[gene] = {}
            with open(result["slice"], "r") as handle:
                next(handle)
                for line in handle:
                    read, start = line.rstrip().split("\t")
                    lookup[gene][read] = int(start)
        return lookup

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
        minimap2_profile: str = "short",
        show_header: bool = True,
    ) -> None:
        from Bio import SeqIO

        if show_header:
            print("\n\033[0;36mRunning mapper2 pipeline...\033[0m\n", file=sys.stderr)

        r1 = str(r1)
        r2 = str(r2) if r2 else None
        ref = str(ref)
        is_paired = r2 is not None

        msa_output = self.output_dir / "ref_seq_msa.aln"
        self._run_spinner(
            f"mafft --auto '{ref}' > '{msa_output}'",
            "Creating Multiple Sequence Alignment",
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

        sequence_positions, _, gene_names = load_msa_mapping(str(ref_msa_tsv))
        if not gene_names:
            raise RuntimeError("No gene names detected in MSA mapping")

        records = []
        gene_order_index = {name: idx for idx, name in enumerate(gene_names)}
        for record in SeqIO.parse(ref, "fasta"):
            gene_name = record.id.split()[0]
            if gene_name in gene_order_index:
                records.append((gene_name, str(record.seq)))
        records.sort(key=lambda item: gene_order_index[item[0]])
        if not records:
            raise RuntimeError("Reference FASTA does not match MSA gene names")

        print(f"Launching per-gene alignments ({len(records)} genes, {threads} threads each)", file=sys.stderr)
        futures = []
        results: List[Dict[str, Path]] = []
        with concurrent.futures.ThreadPoolExecutor(max_workers=len(records)) as executor:
            for gene_name, seq in records:
                futures.append(
                    executor.submit(
                        self._process_gene,
                        gene_name,
                        seq,
                        r1,
                        r2,
                        aligner,
                        threads,
                        minimap2_profile,
                        is_paired,
                        sequence_positions,
                        ref_msa_tsv,
                    )
                )
            for future in concurrent.futures.as_completed(futures):
                results.append(future.result())

        slice_lookup = self._build_slice_lookup(results)
        final_unique, consensus_counts, all_reads = self._reconcile_results(results, gene_names)

        slice_table = self.output_dir / f"{self.prefix}_slice_positions.tsv"
        with open(slice_table, "w") as out:
            header = ["Read_Name"] + [f"{gene}_MSA_Start" for gene in gene_names]
            out.write("\t".join(header) + "\n")
            for read in sorted(all_reads):
                row = [read]
                for gene in gene_names:
                    start = slice_lookup.get(gene, {}).get(read)
                    row.append(str(start) if start is not None else "NA")
                out.write("\t".join(row) + "\n")

        identical = 0
        disagreement = 0
        for read in sorted(all_reads):
            starts = []
            missing = False
            for gene in gene_names:
                start = slice_lookup.get(gene, {}).get(read)
                if start is None:
                    missing = True
                    break
                starts.append(start)
            if not missing and len(set(starts)) == 1:
                identical += 1
            else:
                disagreement += 1

        slice_summary = self.output_dir / f"{self.prefix}_slice_overlap_summary.txt"
        with open(slice_summary, "w") as out:
            out.write(f"Reads with identical slices across all genes: {identical}\n")
            out.write(f"Reads with any slice disagreement:            {disagreement}\n")

        consensus_summary = self.output_dir / f"{self.prefix}_consensus_summary.tsv"
        with open(consensus_summary, "w") as out:
            out.write("Outcome\tCount\n")
            out.write(f"Single gene consensus\t{consensus_counts['single_gene']}\n")
            out.write(f"Unanimous NONE\t{consensus_counts['all_none']}\n")
            out.write(f"Conflict -> NONE\t{consensus_counts['conflict']}\n")

        fastq_dir = self.output_dir / f"{self.prefix}_fastq"
        bam_dir = self.output_dir / f"{self.prefix}_bam"
        output_cmd = [
            "python",
            str(self.pipeline_dir / "output.py"),
            "--tsv",
            str(final_unique),
            "--r1",
            r1,
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
        ]
        if r2:
            output_cmd.extend(["--r2", r2])
        self._run_progress(output_cmd, "Writing output files")

        print(f"\nPipeline complete. Outputs saved to: {self.output_dir}\n", file=sys.stderr)
