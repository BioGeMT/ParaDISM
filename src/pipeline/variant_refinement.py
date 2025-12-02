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
        
        # Call FreeBayes
        try:
            with open(gene_vcf_tmp, "w") as vcf_out, open(
                per_gene_vcf_dir / f"{gene_name}.log", "w"
            ) as log_out:
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
        except subprocess.CalledProcessError:
            print(f"    Warning: FreeBayes failed for {gene_name}; creating empty VCF", file=sys.stderr)
            # Create empty VCF
            with open(gene_vcf, "w") as f:
                f.write("##fileformat=VCFv4.2\n")
                f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            gene_vcf_tmp.unlink(missing_ok=True)
            continue
        
        # Filter to SNPs only (matching existing ParaDISM approach)
        if VARIANT_FILTER_SNPS_ONLY:
            try:
                # Use bcftools to filter SNPs
                with open(gene_vcf, "w") as vcf_out:
                    subprocess.run(
                        ["bcftools", "view", "-v", "snps", str(gene_vcf_tmp)],
                        stdout=vcf_out,
                        stderr=subprocess.DEVNULL,
                        check=True,
                    )
                gene_vcf_tmp.unlink()
            except subprocess.CalledProcessError:
                # Fallback to awk filter (matching existing code)
                print(f"    Using awk filter for {gene_name}...", file=sys.stderr)
                with open(gene_vcf, "w") as vcf_out:
                    with open(gene_vcf_tmp, "r") as vcf_in:
                        for line in vcf_in:
                            if line.startswith("#"):
                                vcf_out.write(line)
                            else:
                                parts = line.strip().split("\t")
                                if len(parts) >= 5:
                                    ref = parts[3]
                                    alt = parts[4]
                                    # SNPs only: single nucleotide substitutions
                                    if len(ref) == 1 and len(alt) == 1:
                                        vcf_out.write(line)
                gene_vcf_tmp.unlink()
        else:
            gene_vcf_tmp.rename(gene_vcf)
        
        if gene_vcf.exists() and gene_vcf.stat().st_size > 100:  # Non-empty VCF
            per_gene_vcfs.append(gene_vcf)
    
    # Merge all per-gene VCFs
    if not per_gene_vcfs:
        print("  No variants found in any gene", file=sys.stderr)
        output_vcf.parent.mkdir(parents=True, exist_ok=True)
        with open(output_vcf, "w") as f:
            f.write("##fileformat=VCFv4.2\n")
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        return
    
    print(f"  Merging {len(per_gene_vcfs)} per-gene VCFs...", file=sys.stderr)
    
    # Use bcftools merge if available, otherwise concatenate
    try:
        output_vcf.parent.mkdir(parents=True, exist_ok=True)
        with open(output_vcf, "w") as vcf_out:
            subprocess.run(
                ["bcftools", "merge", "-m", "all"] + [str(vcf) for vcf in per_gene_vcfs],
                stdout=vcf_out,
                stderr=subprocess.DEVNULL,
                check=True,
            )
    except subprocess.CalledProcessError:
        # Fallback: simple concatenation (matching existing code)
        print("  Using simple VCF concatenation...", file=sys.stderr)
        output_vcf.parent.mkdir(parents=True, exist_ok=True)
        with open(output_vcf, "w") as out_f:
            # Write header from first VCF
            with open(per_gene_vcfs[0], "r") as first_vcf:
                for line in first_vcf:
                    if line.startswith("#"):
                        out_f.write(line)
                    else:
                        break
            
            # Write variants from all VCFs
            for vcf_file in per_gene_vcfs:
                with open(vcf_file, "r") as vcf_in:
                    for line in vcf_in:
                        if not line.startswith("#"):
                            out_f.write(line)
        
        # Sort by chromosome and position
        try:
            sorted_vcf = output_vcf.with_suffix(".sorted.vcf")
            with open(sorted_vcf, "w") as sorted_out:
                subprocess.run(
                    ["bcftools", "sort", str(output_vcf)],
                    stdout=sorted_out,
                    stderr=subprocess.DEVNULL,
                    check=True,
                )
            sorted_vcf.rename(output_vcf)
        except subprocess.CalledProcessError:
            # If sorting fails, keep unsorted VCF
            pass
    
    print(f"  Combined VCF saved to: {output_vcf}", file=sys.stderr)


def apply_variants_to_reference(
    reference: Path,
    vcf: Path,
    output_ref: Path,
) -> None:
    """
    Apply variants to reference using bcftools consensus.
    
    Args:
        reference: Input reference FASTA
        vcf: Variants VCF file
        output_ref: Output path for updated reference
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
        return
    
    print(f"  Applying variants to reference...", file=sys.stderr)
    
    output_ref.parent.mkdir(parents=True, exist_ok=True)
    
    # Apply variants using bcftools consensus
    try:
        with open(output_ref, "w") as ref_out:
            subprocess.run(
                ["bcftools", "consensus", "-f", str(reference), str(vcf)],
                stdout=ref_out,
                stderr=subprocess.PIPE,
                check=True,
            )
    except subprocess.CalledProcessError as e:
        print(f"  Error applying variants: {e.stderr.decode() if e.stderr else 'Unknown error'}", file=sys.stderr)
        # Fallback: copy reference unchanged
        import shutil
        shutil.copy(reference, output_ref)
        print("  Warning: Variant application failed; using original reference", file=sys.stderr)
    
    # Index the new reference
    subprocess.run(
        ["samtools", "faidx", str(output_ref)],
        stderr=subprocess.DEVNULL,
        check=False,
    )
    
    print(f"  Updated reference saved to: {output_ref}", file=sys.stderr)


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
    
    apply_variants_to_reference(
        reference=Path(args.reference),
        vcf=Path(args.output_vcf),
        output_ref=Path(args.output_ref),
    )


if __name__ == "__main__":
    main()
