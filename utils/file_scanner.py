#!/usr/bin/env python3
"""
File scanning utilities for the homologous-region mapper pipeline.
Extract metadata from FASTQ, FASTA, and SAM files for UI display.
"""

import os
import subprocess
from typing import Dict, List, Tuple, Optional
from pathlib import Path
import re


def format_bytes(bytes_size: int) -> str:
    """Convert bytes to human-readable format."""
    for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
        if bytes_size < 1024.0:
            return f"{bytes_size:.1f} {unit}"
        bytes_size /= 1024.0
    return f"{bytes_size:.1f} PB"


def scan_fastq_metadata(fq_path: str) -> Dict:
    """
    Quick scan of FASTQ file for display metadata.

    Args:
        fq_path: Path to FASTQ file

    Returns:
        Dict with size, read count, and quality encoding info
    """
    try:
        # Get file size
        size_bytes = os.path.getsize(fq_path)
        size_human = format_bytes(size_bytes)

        # Estimate read count from first 100K lines
        result = subprocess.run(
            f"head -n 100000 '{fq_path}' | wc -l",
            shell=True,
            capture_output=True,
            text=True,
            timeout=5
        )
        sample_lines = int(result.stdout.strip())

        # Get total line count (fast approximation)
        result = subprocess.run(
            f"wc -l < '{fq_path}'",
            shell=True,
            capture_output=True,
            text=True,
            timeout=10
        )
        total_lines = int(result.stdout.strip())

        # FASTQ has 4 lines per read
        read_count = total_lines // 4

        # Detect quality encoding from first few reads
        result = subprocess.run(
            f"head -n 400 '{fq_path}' | sed -n '4~4p' | head -c 1000",
            shell=True,
            capture_output=True,
            text=True,
            timeout=5
        )
        quality_chars = result.stdout

        # Check ASCII range
        min_qual = min(ord(c) for c in quality_chars if c.strip())
        if min_qual < 59:  # Phred+33 uses ASCII 33-73
            quality_encoding = "Phred+33"
        else:  # Phred+64 uses ASCII 64-104
            quality_encoding = "Phred+64"

        # Estimate mean read length from first 100 reads
        result = subprocess.run(
            f"head -n 400 '{fq_path}' | sed -n '2~4p' | head -n 100 | awk '{{sum+=length($0)}} END {{print int(sum/NR)}}'",
            shell=True,
            capture_output=True,
            text=True,
            timeout=5
        )
        read_length_mean = int(result.stdout.strip()) if result.stdout.strip() else 0

        return {
            "size_bytes": size_bytes,
            "size_human": size_human,
            "read_count": read_count,
            "read_length_mean": read_length_mean,
            "quality_encoding": quality_encoding,
            "exists": True,
            "error": None
        }

    except Exception as e:
        return {
            "size_bytes": 0,
            "size_human": "0 B",
            "read_count": 0,
            "read_length_mean": 0,
            "quality_encoding": "Unknown",
            "exists": False,
            "error": str(e)
        }


def scan_fasta_metadata(fa_path: str) -> Dict:
    """
    Scan FASTA file for sequence information.

    Args:
        fa_path: Path to FASTA file

    Returns:
        Dict with size, sequence count, and sequence IDs
    """
    try:
        from Bio import SeqIO

        size_bytes = os.path.getsize(fa_path)
        size_human = format_bytes(size_bytes)

        sequences = list(SeqIO.parse(fa_path, "fasta"))
        num_sequences = len(sequences)
        total_length = sum(len(seq.seq) for seq in sequences)
        sequence_ids = [seq.id for seq in sequences]

        # Find longest sequence
        longest_seq = max(sequences, key=lambda s: len(s.seq)) if sequences else None
        longest_name = longest_seq.id if longest_seq else ""
        longest_length = len(longest_seq.seq) if longest_seq else 0

        return {
            "size_bytes": size_bytes,
            "size_human": size_human,
            "num_sequences": num_sequences,
            "total_length": total_length,
            "sequence_ids": sequence_ids,
            "longest_name": longest_name,
            "longest_length": longest_length,
            "exists": True,
            "error": None
        }

    except Exception as e:
        return {
            "size_bytes": 0,
            "size_human": "0 B",
            "num_sequences": 0,
            "total_length": 0,
            "sequence_ids": [],
            "longest_name": "",
            "longest_length": 0,
            "exists": False,
            "error": str(e)
        }


def scan_sam_metadata(sam_path: str) -> Dict:
    """
    Validate and extract SAM file metadata.

    Args:
        sam_path: Path to SAM file

    Returns:
        Dict with size, read counts, and validation status
    """
    try:
        size_bytes = os.path.getsize(sam_path)
        size_human = format_bytes(size_bytes)

        # Check for valid header
        result = subprocess.run(
            f"samtools view -H '{sam_path}' 2>/dev/null | wc -l",
            shell=True,
            capture_output=True,
            text=True,
            timeout=10
        )
        has_header = int(result.stdout.strip()) > 0

        # Count aligned reads
        result = subprocess.run(
            f"samtools view -c -F 4 '{sam_path}' 2>/dev/null",
            shell=True,
            capture_output=True,
            text=True,
            timeout=30
        )
        aligned_reads = int(result.stdout.strip()) if result.stdout.strip() else 0

        # Count total reads
        result = subprocess.run(
            f"samtools view -c '{sam_path}' 2>/dev/null",
            shell=True,
            capture_output=True,
            text=True,
            timeout=30
        )
        total_reads = int(result.stdout.strip()) if result.stdout.strip() else 0

        # Check for MD tags
        result = subprocess.run(
            f"samtools view '{sam_path}' 2>/dev/null | head -n 100 | grep -c 'MD:Z:' || echo 0",
            shell=True,
            capture_output=True,
            text=True,
            timeout=10
        )
        has_md_tags = int(result.stdout.strip()) > 0

        # Get reference names from header
        result = subprocess.run(
            f"samtools view -H '{sam_path}' 2>/dev/null | grep '^@SQ' | cut -f2 | sed 's/SN://'",
            shell=True,
            capture_output=True,
            text=True,
            timeout=10
        )
        references = [ref.strip() for ref in result.stdout.strip().split('\n') if ref.strip()]

        return {
            "size_bytes": size_bytes,
            "size_human": size_human,
            "has_header": has_header,
            "aligned_reads": aligned_reads,
            "total_reads": total_reads,
            "alignment_rate": (aligned_reads / total_reads * 100) if total_reads > 0 else 0,
            "has_md_tags": has_md_tags,
            "references": references,
            "exists": True,
            "error": None
        }

    except Exception as e:
        return {
            "size_bytes": 0,
            "size_human": "0 B",
            "has_header": False,
            "aligned_reads": 0,
            "total_reads": 0,
            "alignment_rate": 0,
            "has_md_tags": False,
            "references": [],
            "exists": False,
            "error": str(e)
        }


def find_paired_fastqs(directory: str = ".") -> List[Tuple[str, str]]:
    """
    Auto-detect R1/R2 FASTQ pairs in a directory.

    Matches files that differ by exactly one character (1→2).

    Args:
        directory: Directory to search (default: current)

    Returns:
        List of (r1_path, r2_path) tuples
    """
    pairs = []
    fastq_files = sorted(list(Path(directory).glob("*.fq")) + list(Path(directory).glob("*.fastq")))
    seen = set()

    for i, file1 in enumerate(fastq_files):
        if str(file1) in seen:
            continue

        name1 = file1.name

        # Try to find a matching pair
        for file2 in fastq_files[i+1:]:
            if str(file2) in seen:
                continue

            name2 = file2.name

            # Files must have same length
            if len(name1) != len(name2):
                continue

            # Count character differences
            diff_count = sum(c1 != c2 for c1, c2 in zip(name1, name2))

            # Must differ by exactly 1 character
            if diff_count != 1:
                continue

            # Find the differing position and check if it's 1→2
            diff_pos = next(i for i, (c1, c2) in enumerate(zip(name1, name2)) if c1 != c2)
            char1, char2 = name1[diff_pos], name2[diff_pos]

            # Valid pair if char1='1' and char2='2'
            if char1 == '1' and char2 == '2':
                pairs.append((str(file1), str(file2)))
                seen.add(str(file1))
                seen.add(str(file2))
                break

    return pairs


def find_references(directory: str = ".") -> List[str]:
    """
    Find FASTA reference files in directory.

    Args:
        directory: Directory to search

    Returns:
        List of reference file paths
    """
    ref_files = []
    for ext in ["*.fa", "*.fasta", "*.fas", "*.fna"]:
        ref_files.extend(Path(directory).glob(ext))

    return sorted([str(f) for f in ref_files])


def find_sam_files(directory: str = ".") -> List[str]:
    """
    Find SAM files in directory.

    Args:
        directory: Directory to search

    Returns:
        List of SAM file paths
    """
    sam_files = list(Path(directory).glob("*.sam"))
    return sorted([str(f) for f in sam_files])
