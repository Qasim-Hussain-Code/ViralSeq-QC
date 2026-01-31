"""
ViralSeq-QC: A zero-dependency toolkit for viral genome quality control.

This package provides comprehensive quality control functions for viral
consensus sequences, including:

- Streamed FASTA parsing (memory efficient)
- Length and N-content filtering
- Advanced metrics (GC content, sequence complexity, homopolymer detection)
- Multiple output formats (TSV, JSON, FASTA)

Example usage:
    >>> from src import parse_fasta, is_high_quality, calculate_gc_content
    >>> for header, seq in parse_fasta("input.fasta"):
    ...     if is_high_quality(seq, min_length=1000, max_n_content=1.0):
    ...         print(f"{header}: GC={calculate_gc_content(seq):.1f}%")
"""

# Core I/O functions
from .input_output import (
    parse_fasta,
    write_tsv_report,
    write_json_report,
    write_fasta,
    validate_fasta_format,
    FastaParseError,
    InvalidFastaError,
)

# QC functions
from .qc import (
    # Core metrics
    calculate_gc_content,
    calculate_n_content,
    check_length,
    is_high_quality,
    # Advanced metrics
    detect_homopolymer_runs,
    calculate_sequence_complexity,
    check_terminal_ns,
    validate_nucleotides,
    get_sequence_metrics,
    # Batch processing
    process_batch,
)

__all__ = [
    # I/O
    "parse_fasta",
    "write_tsv_report",
    "write_json_report",
    "write_fasta",
    "validate_fasta_format",
    "FastaParseError",
    "InvalidFastaError",
    # Core QC
    "calculate_gc_content",
    "calculate_n_content",
    "check_length",
    "is_high_quality",
    # Advanced QC
    "detect_homopolymer_runs",
    "calculate_sequence_complexity",
    "check_terminal_ns",
    "validate_nucleotides",
    "get_sequence_metrics",
    # Batch
    "process_batch",
]

__version__ = "2.0.0"
__author__ = "Qasim Hussain"