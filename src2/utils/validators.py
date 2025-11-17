#!/usr/bin/env python3
"""Validation helpers for basic FASTQ/FASTA/SAM checks."""

from typing import Dict, List
from pathlib import Path
import subprocess


def validate_fastq_pair(r1_path: str, r2_path: str, r1_metadata: Dict, r2_metadata: Dict) -> List[Dict]:
    """Validate compatibility of an R1/R2 FASTQ pair."""
    results = []

    # Check if files exist
    if not r1_metadata["exists"]:
        results.append({
            "status": "error",
            "message": f"✗ R1 file not found or cannot be read: {r1_path}",
            "blocking": True
        })
        return results

    if not r2_metadata["exists"]:
        results.append({
            "status": "error",
            "message": f"✗ R2 file not found or cannot be read: {r2_path}",
            "blocking": True
        })
        return results

    # Check read counts match
    r1_count = r1_metadata["read_count"]
    r2_count = r2_metadata["read_count"]

    if r1_count != r2_count:
        diff_pct = abs(r1_count - r2_count) / max(r1_count, r2_count) * 100
        if diff_pct > 0.1:  # More than 0.1% difference is an error
            results.append({
                "status": "error",
                "message": f"✗ R1 and R2 have different read counts: {r1_count:,} vs {r2_count:,}",
                "blocking": True
            })
        else:
            results.append({
                "status": "warning",
                "message": f"⚠ R1 and R2 read counts slightly different: {r1_count:,} vs {r2_count:,} (within tolerance)",
                "blocking": False
            })
    else:
        results.append({
            "status": "success",
            "message": f"✓ R1 and R2 have matching read counts ({r1_count:,} pairs)",
            "blocking": False
        })

    # Check quality encoding compatibility
    if r1_metadata["quality_encoding"] != r2_metadata["quality_encoding"]:
        results.append({
            "status": "warning",
            "message": f"⚠ R1 and R2 have different quality encodings: {r1_metadata['quality_encoding']} vs {r2_metadata['quality_encoding']}",
            "blocking": False
        })
    else:
        results.append({
            "status": "success",
            "message": f"✓ Quality encoding detected: {r1_metadata['quality_encoding']}",
            "blocking": False
        })

    # Check read lengths are similar
    r1_len = r1_metadata["read_length_mean"]
    r2_len = r2_metadata["read_length_mean"]
    if r1_len > 0 and r2_len > 0:
        if abs(r1_len - r2_len) > 10:  # More than 10bp difference
            results.append({
                "status": "warning",
                "message": f"⚠ R1 and R2 have different mean read lengths: {r1_len}bp vs {r2_len}bp",
                "blocking": False
            })
        else:
            results.append({
                "status": "success",
                "message": f"✓ Mean read length: {r1_len}bp (R1), {r2_len}bp (R2)",
                "blocking": False
            })

    return results


def validate_fastq_single(r1_path: str, r1_metadata: Dict) -> List[Dict]:
    """Validate a single-end FASTQ file."""
    results = []

    # Check if file exists
    if not r1_metadata["exists"]:
        results.append({
            "status": "error",
            "message": f"✗ FASTQ file not found or cannot be read: {r1_path}",
            "blocking": True
        })
        return results

    # Show read count
    read_count = r1_metadata["read_count"]
    results.append({
        "status": "success",
        "message": f"✓ FASTQ contains {read_count:,} reads",
        "blocking": False
    })

    # Show quality encoding
    results.append({
        "status": "success",
        "message": f"✓ Quality encoding detected: {r1_metadata['quality_encoding']}",
        "blocking": False
    })

    # Show read length
    read_len = r1_metadata["read_length_mean"]
    if read_len > 0:
        results.append({
            "status": "success",
            "message": f"✓ Mean read length: {read_len}bp",
            "blocking": False
        })

    return results


def validate_fasta(fa_path: str, fa_metadata: Dict) -> List[Dict]:
    """Validate a reference FASTA file."""
    results = []

    # Check if file exists
    if not fa_metadata["exists"]:
        results.append({
            "status": "error",
            "message": f"✗ Reference file not found or cannot be read: {fa_path}",
            "blocking": True
        })
        return results

    # Check if file has sequences
    if fa_metadata["num_sequences"] == 0:
        results.append({
            "status": "error",
            "message": "✗ Reference file contains no sequences",
            "blocking": True
        })
        return results

    results.append({
        "status": "success",
        "message": f"✓ Reference contains {fa_metadata['num_sequences']} sequence(s)",
        "blocking": False
    })

    # Show sequence IDs
    seq_ids = fa_metadata["sequence_ids"]
    if len(seq_ids) <= 10:
        results.append({
            "status": "success",
            "message": f"✓ Sequences: {', '.join(seq_ids)}",
            "blocking": False
        })
    else:
        results.append({
            "status": "success",
            "message": f"✓ Sequences: {', '.join(seq_ids[:10])} ... and {len(seq_ids) - 10} more",
            "blocking": False
        })

    # Check total length
    if fa_metadata["total_length"] > 0:
        total_kb = fa_metadata["total_length"] / 1000
        results.append({
            "status": "success",
            "message": f"✓ Total reference length: {total_kb:.1f} Kbp",
            "blocking": False
        })

    return results


def validate_sam(sam_path: str, sam_metadata: Dict) -> List[Dict]:
    """Validate SAM for header, alignment rate, MD tags, and refs."""
    results = []

    # Check if file exists
    if not sam_metadata["exists"]:
        results.append({
            "status": "error",
            "message": f"✗ SAM file not found or cannot be read: {sam_path}",
            "blocking": True
        })
        return results

    # Check for valid header
    if not sam_metadata["has_header"]:
        results.append({
            "status": "error",
            "message": "✗ SAM file has no valid header",
            "blocking": True
        })
    else:
        results.append({
            "status": "success",
            "message": "✓ SAM file has valid header",
            "blocking": False
        })

    # Check for aligned reads
    aligned = sam_metadata["aligned_reads"]
    total = sam_metadata["total_reads"]
    alignment_rate = sam_metadata["alignment_rate"]

    if aligned == 0:
        results.append({
            "status": "error",
            "message": "✗ SAM file contains no aligned reads",
            "blocking": True
        })
    else:
        results.append({
            "status": "success",
            "message": f"✓ SAM contains {aligned:,} aligned reads ({alignment_rate:.1f}% alignment rate)",
            "blocking": False
        })

    # Check for MD tags (required for mapper pipeline)
    if not sam_metadata["has_md_tags"]:
        results.append({
            "status": "error",
            "message": "✗ SAM file missing required MD tags (rerun alignment with --MD flag or calmd)",
            "blocking": True
        })
    else:
        results.append({
            "status": "success",
            "message": "✓ SAM file has required MD tags",
            "blocking": False
        })

    # Show reference sequences
    refs = sam_metadata["references"]
    if refs:
        if len(refs) <= 10:
            results.append({
                "status": "success",
                "message": f"✓ Reference sequences: {', '.join(refs)}",
                "blocking": False
            })
        else:
            results.append({
                "status": "success",
                "message": f"✓ Reference sequences: {', '.join(refs[:10])} ... and {len(refs) - 10} more",
                "blocking": False
            })

    return results


def has_blocking_errors(validation_results: List[Dict]) -> bool:
    """Return True if any result is a blocking error."""
    return any(result["status"] == "error" and result["blocking"] for result in validation_results)
