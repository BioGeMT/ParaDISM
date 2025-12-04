#!/usr/bin/env python3
"""
Variant calling and reference refinement for iterative ParaDISM pipeline.

This script:
1. Calls variants from per-gene BAM files using FreeBayes
2. Filters to SNPs only (matching existing ParaDISM approach)
3. Merges per-gene VCFs
4. Applies variants to reference using bcftools consensus
"""

import argparse
import subprocess
import sys
from pathlib import Path
from typing import Optional

# Hardcoded defaults (matching existing ParaDISM variant calling)
VARIANT_MIN_ALT = 5
VARIANT_PLOIDY = 2
VARIANT_FILTER_SNPS_ONLY = True
# Quality filtering thresholds
VARIANT_MIN_QUAL = 20  # Minimum quality score (Phred-scaled)
VARIANT_MIN_DP = 10    # Minimum total depth
VARIANT_MIN_AF = 0.05  # Minimum allele frequency (5%)


def check_tool(tool_name: str) -> bool:
    """Check if a required tool is available."""
    try:
        subprocess.run(
            [tool_name, "--version"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=True,
        )
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False


def call_variants_from_bams(
    bam_dir: Path,
    reference: Path,
    output_vcf: Path,
    per_gene_vcf_dir: Path,
) -> None:
    """
    Call variants from per-gene BAM files using FreeBayes.
    
    Args:
        bam_dir: Directory containing per-gene BAM files (*.sorted.bam)
        reference: Reference FASTA file
        output_vcf: Output path for combined VCF
        per_gene_vcf_dir: Directory to save per-gene VCFs
    """
    if not check_tool("freebayes"):
        print("Error: freebayes not found. Please install FreeBayes.", file=sys.stderr)
        sys.exit(1)
    
    if not check_tool("bcftools"):
        print("Error: bcftools not found. Please install bcftools.", file=sys.stderr)
        sys.exit(1)
    
    bam_dir = Path(bam_dir)
    reference = Path(reference)
    per_gene_vcf_dir = Path(per_gene_vcf_dir)
    output_vcf = Path(output_vcf)
    
    per_gene_vcf_dir.mkdir(parents=True, exist_ok=True)
    
    # Find all per-gene BAM files
    bam_files = list(bam_dir.glob("*.sorted.bam"))
    if not bam_files:
        print(f"Warning: No BAM files found in {bam_dir}", file=sys.stderr)
        # Create empty VCF
        output_vcf.parent.mkdir(parents=True, exist_ok=True)
        with open(output_vcf, "w") as f:
            f.write("##fileformat=VCFv4.2\n")
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        return
    
    per_gene_vcfs = []
    
    # Call variants for each gene BAM
    for bam_file in sorted(bam_files):
        # Extract gene name from filename (e.g., prefix_GENE.sorted.bam -> GENE)
        gene_name = bam_file.stem.replace(".sorted", "")
        # Remove prefix if present (everything before last underscore)
        if "_" in gene_name:
            gene_name = gene_name.rsplit("_", 1)[-1]
        
        gene_vcf = per_gene_vcf_dir / f"{gene_name}.vcf"
        gene_vcf_tmp = per_gene_vcf_dir / f"{gene_name}.tmp.vcf"
        
        print(f"  Calling variants for {gene_name}...", file=sys.stderr)
        
        # Call FreeBayes - fail fast if it fails
        with open(gene_vcf_tmp, "w") as vcf_out, open(
            per_gene_vcf_dir / f"{gene_name}.log", "w"
        ) as log_out:
            try:
                subprocess.run(
                    [
                        "freebayes",
                        "--bam", str(bam_file),
                        "--fasta-reference", str(reference),
                        "--ploidy", str(VARIANT_PLOIDY),
                        "--min-alternate-count", str(VARIANT_MIN_ALT),
                        "--region", gene_name,
                        "-i",
                    ],
                    stdout=vcf_out,
                    stderr=log_out,
                    check=True,
                )
            except subprocess.CalledProcessError as e:
                # Read log file for error details
                log_path = per_gene_vcf_dir / f"{gene_name}.log"
                error_details = ""
                if log_path.exists():
                    with open(log_path, "r") as log_f:
                        error_details = log_f.read()
                print(f"Error: FreeBayes failed for {gene_name}", file=sys.stderr)
                if error_details:
                    print(f"  FreeBayes stderr: {error_details}", file=sys.stderr)
                sys.exit(1)
        
        if VARIANT_FILTER_SNPS_ONLY:
            snps_vcf = gene_vcf_tmp.with_suffix(".snps.vcf")
            with open(snps_vcf, "w") as vcf_out:
                result = subprocess.run(
                    ["bcftools", "view", "-v", "snps", str(gene_vcf_tmp)],
                    stdout=vcf_out,
                    stderr=subprocess.PIPE,
                    check=True,
                )
            gene_vcf_tmp.unlink()
            snps_vcf.rename(gene_vcf)
        else:
            gene_vcf_tmp.rename(gene_vcf)
        
        # Only keep VCFs that contain at least one variant line
        has_variants = False
        with open(gene_vcf, "r") as vf:
            for line in vf:
                if line.startswith("#"):
                    continue
                has_variants = True
                break
        if has_variants:
            per_gene_vcfs.append(gene_vcf)
        else:
            # Remove empty VCF to avoid merge failures; keep log for reference
            gene_vcf.unlink(missing_ok=True)
    
    # Merge all per-gene VCFs
    if not per_gene_vcfs:
        print("  No variants found in any gene", file=sys.stderr)
        output_vcf.parent.mkdir(parents=True, exist_ok=True)
        with open(output_vcf, "w") as f:
            f.write("##fileformat=VCFv4.2\n")
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        return

    print(f"  Merging {len(per_gene_vcfs)} per-gene VCFs...", file=sys.stderr)

    output_vcf.parent.mkdir(parents=True, exist_ok=True)

    # If only one VCF, just copy it (no need to merge)
    if len(per_gene_vcfs) == 1:
        print(f"  Only one gene with variants, copying VCF directly...", file=sys.stderr)
        import shutil
        shutil.copy(per_gene_vcfs[0], output_vcf)

        # Still need to compress and index for bcftools consensus
        vcf_gz = output_vcf.with_suffix(".vcf.gz")
        print(f"  Compressing VCF for bcftools consensus...", file=sys.stderr)
        try:
            subprocess.run(
                ["bgzip", "-c", str(output_vcf)],
                stdout=open(vcf_gz, "wb"),
                stderr=subprocess.PIPE,
                check=True,
            )
            subprocess.run(
                ["bcftools", "index", str(vcf_gz)],
                stderr=subprocess.PIPE,
                check=True,
            )
            print(f"  Compressed VCF saved to: {vcf_gz}", file=sys.stderr)
        except subprocess.CalledProcessError as e:
            print(f"Error: bgzip/index failed: {e.stderr.decode() if e.stderr else 'Unknown error'}", file=sys.stderr)
            sys.exit(1)

        print(f"  Combined VCF saved to: {output_vcf}", file=sys.stderr)
        return

    # bcftools merge requires bgzipped VCFs - compress per-gene VCFs first
    print(f"  Compressing per-gene VCFs for merging...", file=sys.stderr)
    compressed_vcfs = []
    for vcf_file in per_gene_vcfs:
        vcf_gz = vcf_file.with_suffix(".vcf.gz")
        try:
            subprocess.run(
                ["bgzip", "-c", str(vcf_file)],
                stdout=open(vcf_gz, "wb"),
                stderr=subprocess.PIPE,
                check=True,
            )
            # Index compressed VCF
            subprocess.run(
                ["bcftools", "index", str(vcf_gz)],
                stderr=subprocess.PIPE,
                check=True,
            )
            compressed_vcfs.append(vcf_gz)
        except subprocess.CalledProcessError as e:
            print(f"Error: Failed to compress {vcf_file}: {e.stderr.decode() if e.stderr else 'Unknown error'}", file=sys.stderr)
            sys.exit(1)

    # Merge VCFs using bcftools merge - fail fast if it fails
    merged_vcf = output_vcf.with_suffix(".merged.vcf")
    try:
        with open(merged_vcf, "w") as vcf_out:
            subprocess.run(
                ["bcftools", "merge", "-m", "all", "--force-samples"] + [str(vcf) for vcf in compressed_vcfs],
                stdout=vcf_out,
                stderr=subprocess.PIPE,
                check=True,
            )
    except subprocess.CalledProcessError as e:
        print(f"Error: bcftools merge failed: {e.stderr.decode() if e.stderr else 'Unknown error'}", file=sys.stderr)
        sys.exit(1)

    # Sort VCF - fail fast if it fails
    try:
        with open(output_vcf, "w") as sorted_out:
            subprocess.run(
                ["bcftools", "sort", str(merged_vcf)],
                stdout=sorted_out,
                stderr=subprocess.PIPE,
                check=True,
            )
        merged_vcf.unlink()  # Remove temporary merged file
    except subprocess.CalledProcessError as e:
        print(f"Error: bcftools sort failed: {e.stderr.decode() if e.stderr else 'Unknown error'}", file=sys.stderr)
        sys.exit(1)
    
    # Compress and index VCF for bcftools consensus - require bgzip, fail fast
    if not check_tool("bgzip"):
        print("Error: bgzip not found. bgzip is required for variant refinement.", file=sys.stderr)
        sys.exit(1)
    
    vcf_gz = output_vcf.with_suffix(".vcf.gz")
    print(f"  Compressing VCF for bcftools consensus...", file=sys.stderr)
    try:
        subprocess.run(
            ["bgzip", "-c", str(output_vcf)],
            stdout=open(vcf_gz, "wb"),
            stderr=subprocess.PIPE,
            check=True,
        )
        # Index the compressed VCF
        subprocess.run(
            ["bcftools", "index", str(vcf_gz)],
            stderr=subprocess.PIPE,
            check=True,
        )
        print(f"  Compressed VCF saved to: {vcf_gz}", file=sys.stderr)
    except subprocess.CalledProcessError as e:
        print(f"Error: bgzip/index failed: {e.stderr.decode() if e.stderr else 'Unknown error'}", file=sys.stderr)
        sys.exit(1)
    
    print(f"  Combined VCF saved to: {output_vcf}", file=sys.stderr)


def apply_variants_to_reference(
    reference: Path,
    vcf: Path,
    output_ref: Path,
) -> bool:
    """
    Apply variants to reference using bcftools consensus.

    Args:
        reference: Input reference FASTA
        vcf: Variants VCF file (uncompressed)
        output_ref: Output path for updated reference

    Returns:
        True if variants were applied, False if no variants found
    """
    if not check_tool("bcftools"):
        print("Error: bcftools not found. Please install bcftools.", file=sys.stderr)
        sys.exit(1)
    
    if not check_tool("samtools"):
        print("Error: samtools not found. Please install samtools.", file=sys.stderr)
        sys.exit(1)
    
    reference = Path(reference)
    vcf = Path(vcf)
    output_ref = Path(output_ref)
    
    # Check if VCF has variants
    has_variants = False
    if vcf.exists():
        with open(vcf, "r") as f:
            for line in f:
                if not line.startswith("#"):
                    has_variants = True
                    break
    
    if not has_variants:
        print("  No variants to apply; copying reference unchanged", file=sys.stderr)
        import shutil
        output_ref.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(reference, output_ref)
        # Index the reference
        subprocess.run(
            ["samtools", "faidx", str(output_ref)],
            stderr=subprocess.DEVNULL,
            check=False,
        )
        return False
    
    print(f"  Applying variants to reference...", file=sys.stderr)
    
    output_ref.parent.mkdir(parents=True, exist_ok=True)
    
    # Require compressed VCF - fail fast if not found
    vcf_gz = vcf.with_suffix(".vcf.gz")
    if not vcf_gz.exists():
        print(f"Error: Compressed VCF not found: {vcf_gz}", file=sys.stderr)
        print(f"  Variant refinement requires bgzipped VCF.", file=sys.stderr)
        sys.exit(1)
    
    # Apply variants using bcftools consensus - fail fast if it fails
    try:
        with open(output_ref, "w") as ref_out:
            subprocess.run(
                ["bcftools", "consensus", "-f", str(reference), str(vcf_gz)],
                stdout=ref_out,
                stderr=subprocess.PIPE,
                check=True,
            )
        print(f"  Successfully applied variants from {vcf_gz.name}", file=sys.stderr)
    except subprocess.CalledProcessError as e:
        print(f"Error: bcftools consensus failed: {e.stderr.decode() if e.stderr else 'Unknown error'}", file=sys.stderr)
        sys.exit(1)

    # Index the new reference - fail fast if it fails
    try:
        subprocess.run(
            ["samtools", "faidx", str(output_ref)],
            stderr=subprocess.PIPE,
            check=True,
        )
    except subprocess.CalledProcessError as e:
        print(f"Error: samtools faidx failed: {e.stderr.decode() if e.stderr else 'Unknown error'}", file=sys.stderr)
        sys.exit(1)

    print(f"  Updated reference saved to: {output_ref}", file=sys.stderr)
    return True


def main():
    """CLI entry point for variant refinement."""
    parser = argparse.ArgumentParser(
        description="Call variants from per-gene BAMs and apply to reference"
    )
    parser.add_argument(
        "--bam-dir",
        required=True,
        help="Directory containing per-gene BAM files (*.sorted.bam)",
    )
    parser.add_argument(
        "--reference",
        required=True,
        help="Reference FASTA file",
    )
    parser.add_argument(
        "--output-vcf",
        required=True,
        help="Output path for combined VCF file",
    )
    parser.add_argument(
        "--output-ref",
        required=True,
        help="Output path for updated reference FASTA",
    )
    parser.add_argument(
        "--per-gene-vcf-dir",
        required=True,
        help="Directory to save per-gene VCF files",
    )
    
    args = parser.parse_args()

    call_variants_from_bams(
        bam_dir=Path(args.bam_dir),
        reference=Path(args.reference),
        output_vcf=Path(args.output_vcf),
        per_gene_vcf_dir=Path(args.per_gene_vcf_dir),
    )

    has_variants = apply_variants_to_reference(
        reference=Path(args.reference),
        vcf=Path(args.output_vcf),
        output_ref=Path(args.output_ref),
    )

    # Exit with status 2 if no variants found (for caller to detect convergence)
    if not has_variants:
        sys.exit(2)


if __name__ == "__main__":
    main()
