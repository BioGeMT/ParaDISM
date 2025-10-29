#!/usr/bin/env python3
"""Fast metadata scanners for FASTQ/FASTA/SAM (UI display)."""

import os
import subprocess
from typing import Dict, List, Tuple, Optional
from pathlib import Path


def format_bytes(bytes_size: int) -> str:
    for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
        if bytes_size < 1024.0:
            return f"{bytes_size:.1f} {unit}"
        bytes_size /= 1024.0
    return f"{bytes_size:.1f} PB"


def scan_fastq_metadata(fq_path: str) -> Dict:
    size_bytes = os.path.getsize(fq_path)
    size_human = format_bytes(size_bytes)
    estimated_read_count = int(size_bytes / 500)
    
    try:
        result = subprocess.run(
            f"wc -l < '{fq_path}'",
            shell=True,
            capture_output=True,
            text=True,
            timeout=5
        )
        if result.returncode == 0:
            total_lines = int(result.stdout.strip())
            read_count = total_lines // 4
        else:
            read_count = estimated_read_count
    except (subprocess.TimeoutExpired, ValueError):
        read_count = estimated_read_count

    result = subprocess.run(
        f"head -n 400 '{fq_path}' 2>/dev/null | sed -n '4~4p' | head -c 1000",
        shell=True,
        capture_output=True,
        text=True,
        timeout=5
    )
    quality_chars = result.stdout
    if quality_chars:
        min_qual = min(ord(c) for c in quality_chars if c.strip())
        quality_encoding = "Phred+33" if min_qual < 59 else "Phred+64"
    else:
        quality_encoding = "Unknown"

    result = subprocess.run(
        f"head -n 400 '{fq_path}' 2>/dev/null | sed -n '2~4p' | head -n 100 | awk '{{sum+=length($0)}} END {{if(NR>0) print int(sum/NR); else print 100}}'",
        shell=True,
        capture_output=True,
        text=True,
        timeout=5
    )
    read_length_mean = int(result.stdout.strip()) if result.stdout.strip() else 100

    return {
        "size_bytes": size_bytes,
        "size_human": size_human,
        "read_count": read_count,
        "read_length_mean": read_length_mean,
        "quality_encoding": quality_encoding,
        "exists": True,
        "error": None
    }


def scan_fasta_metadata(fa_path: str) -> Dict:
    size_bytes = os.path.getsize(fa_path)
    size_human = format_bytes(size_bytes)
    sequence_ids = []
    total_length = 0
    longest_name = ""
    longest_length = 0
    
    with open(fa_path, 'r') as f:
        current_id = None
        current_length = 0
        
        for line in f:
            if line.startswith('>'):
                if current_id is not None:
                    sequence_ids.append(current_id)
                    total_length += current_length
                    if current_length > longest_length:
                        longest_length = current_length
                        longest_name = current_id
                
                current_id = line[1:].strip().split()[0]
                current_length = 0
            else:
                current_length += len(line.strip())
        
        if current_id is not None:
            sequence_ids.append(current_id)
            total_length += current_length
            if current_length > longest_length:
                longest_length = current_length
                longest_name = current_id
    
    return {
        "size_bytes": size_bytes,
        "size_human": size_human,
        "num_sequences": len(sequence_ids),
        "total_length": total_length,
        "sequence_ids": sequence_ids,
        "longest_name": longest_name,
        "longest_length": longest_length,
        "exists": True,
        "error": None
    }


def scan_sam_metadata(sam_path: str) -> Dict:
    size_bytes = os.path.getsize(sam_path)
    size_human = format_bytes(size_bytes)

    result = subprocess.run(
        f"samtools view -H '{sam_path}' 2>/dev/null | wc -l",
        shell=True,
        capture_output=True,
        text=True,
        timeout=10
    )
    has_header = int(result.stdout.strip()) > 0

    result = subprocess.run(
        f"samtools view -c -F 4 '{sam_path}' 2>/dev/null",
        shell=True,
        capture_output=True,
        text=True,
        timeout=30
    )
    aligned_reads = int(result.stdout.strip()) if result.stdout.strip() else 0

    result = subprocess.run(
        f"samtools view -c '{sam_path}' 2>/dev/null",
        shell=True,
        capture_output=True,
        text=True,
        timeout=30
    )
    total_reads = int(result.stdout.strip()) if result.stdout.strip() else 0

    result = subprocess.run(
        f"samtools view '{sam_path}' 2>/dev/null | head -n 100 | grep -c 'MD:Z:' || echo 0",
        shell=True,
        capture_output=True,
        text=True,
        timeout=10
    )
    has_md_tags = int(result.stdout.strip()) > 0

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


def find_paired_fastqs(directory: str = ".") -> List[Tuple[str, str]]:
    pairs = []
    fastq_files = sorted(list(Path(directory).glob("*.fq")) + list(Path(directory).glob("*.fastq")))
    seen = set()

    for i, file1 in enumerate(fastq_files):
        if str(file1) in seen:
            continue

        name1 = file1.name
        for file2 in fastq_files[i+1:]:
            if str(file2) in seen:
                continue

            name2 = file2.name
            if len(name1) != len(name2):
                continue

            diff_count = sum(c1 != c2 for c1, c2 in zip(name1, name2))
            if diff_count != 1:
                continue

            diff_pos = next(i for i, (c1, c2) in enumerate(zip(name1, name2)) if c1 != c2)
            char1, char2 = name1[diff_pos], name2[diff_pos]
            if char1 == '1' and char2 == '2':
                pairs.append((str(file1), str(file2)))
                seen.add(str(file1))
                seen.add(str(file2))
                break

    return pairs


def find_references(directory: str = ".") -> List[str]:
    ref_files = []
    for ext in ["*.fa", "*.fasta", "*.fas", "*.fna"]:
        ref_files.extend(Path(directory).glob(ext))
    return sorted([str(f) for f in ref_files])


def find_sam_files(directory: str = ".") -> List[str]:
    sam_files = list(Path(directory).glob("*.sam"))
    return sorted([str(f) for f in sam_files])
