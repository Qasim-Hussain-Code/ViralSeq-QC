"""
ViralSeq-QC: Input/Output Module

This module handles file parsing and writing operations, including:
- FASTA file parsing (standard and gzip-compressed)
- Output writing in multiple formats (TSV, JSON, FASTA)
- Robust error handling for malformed input files
"""

import gzip
import json
import logging
from pathlib import Path
from typing import Dict, Generator, List, Optional, TextIO, Tuple

logger = logging.getLogger(__name__)


# =============================================================================
# Custom Exceptions
# =============================================================================

class FastaParseError(Exception):
    """Exception raised for errors during FASTA parsing."""
    pass


class InvalidFastaError(FastaParseError):
    """Exception raised when FASTA format is invalid."""
    pass


# =============================================================================
# FASTA Parsing
# =============================================================================

def parse_fasta(file_path: str) -> Generator[Tuple[str, str], None, None]:
    """
    Read a FASTA file and yield (header, sequence) tuples one by one.

    Supports both standard and gzip-compressed (.gz) FASTA files.
    Handles multiline sequences using a buffer strategy for memory efficiency.

    Args:
        file_path: Path to the input FASTA file (.fasta, .fa, .fna, or .gz).

    Yields:
        Tuple of (header, sequence) for each record.

    Raises:
        FastaParseError: If file cannot be read or is malformed.

    Example:
        >>> for header, seq in parse_fasta("sequences.fasta"):
        ...     print(f"{header}: {len(seq)} bp")
    """
    path = Path(file_path)

    if not path.exists():
        raise FastaParseError(f"File not found: {file_path}")

    # Determine if gzip compressed
    is_gzipped = path.suffix.lower() == '.gz'

    try:
        if is_gzipped:
            file_handle = gzip.open(file_path, 'rt', encoding='utf-8')
        else:
            file_handle = open(file_path, encoding='utf-8')

        yield from _parse_fasta_handle(file_handle, file_path)
        file_handle.close()

    except gzip.BadGzipFile as e:
        raise FastaParseError(f"Invalid gzip file: {file_path}") from e
    except UnicodeDecodeError as e:
        raise FastaParseError(f"Encoding error in {file_path}: {e}") from e
    except Exception as e:
        raise FastaParseError(f"Error reading {file_path}: {e}") from e


def _parse_fasta_handle(
    file_handle: TextIO,
    file_path: str
) -> Generator[Tuple[str, str], None, None]:
    """
    Internal function to parse FASTA from a file handle.

    Args:
        file_handle: Open file handle (text mode).
        file_path: Original file path (for error messages).

    Yields:
        Tuple of (header, sequence) for each record.
    """
    header = None
    sequence_buffer: List[str] = []
    line_number = 0
    record_count = 0

    for line in file_handle:
        line_number += 1
        line = line.strip()

        # Skip empty lines
        if not line:
            continue

        # Skip comment lines (some FASTA files have these)
        if line.startswith(';'):
            continue

        if line.startswith(">"):
            # Yield previous record if exists
            if header is not None:
                full_sequence = "".join(sequence_buffer)
                record_count += 1
                logger.debug(f"Parsed record {record_count}: {header[:50]}...")
                yield header, full_sequence

            # Start new record
            header = line[1:].strip()  # Remove '>' and whitespace
            if not header:
                logger.warning(f"Empty header at line {line_number} in {file_path}")
                header = f"unnamed_sequence_{line_number}"
            sequence_buffer = []

        else:
            # Accumulate sequence lines
            if header is None:
                raise InvalidFastaError(
                    f"Sequence data before header at line {line_number} in {file_path}"
                )
            # Clean sequence line (remove spaces, validate characters optionally)
            clean_line = line.replace(' ', '').replace('\t', '')
            sequence_buffer.append(clean_line)

    # Yield the last record
    if header is not None:
        full_sequence = "".join(sequence_buffer)
        record_count += 1
        logger.debug(f"Parsed record {record_count}: {header[:50]}...")
        yield header, full_sequence

    logger.info(f"Finished parsing {file_path}: {record_count} records")


def validate_fasta_format(file_path: str) -> Tuple[bool, str]:
    """
    Validate that a file is properly formatted FASTA.

    Args:
        file_path: Path to the file to validate.

    Returns:
        Tuple of (is_valid, message).
    """
    try:
        record_count = 0
        for header, sequence in parse_fasta(file_path):
            record_count += 1
            if not sequence:
                return False, f"Empty sequence for header: {header}"

        if record_count == 0:
            return False, "No FASTA records found in file"

        return True, f"Valid FASTA file with {record_count} records"

    except FastaParseError as e:
        return False, str(e)


# =============================================================================
# Output Writing
# =============================================================================

def write_tsv_report(
    results: List[Dict],
    output_path: str,
    columns: Optional[List[str]] = None
) -> None:
    """
    Write QC results to a TSV (Tab-Separated Values) file.

    Args:
        results: List of dictionaries containing QC results.
        output_path: Path for the output TSV file.
        columns: Optional list of columns to include (default: all columns).

    Raises:
        IOError: If file cannot be written.
    """
    if not results:
        logger.warning("No results to write to TSV")
        return

    # Determine columns
    if columns is None:
        columns = list(results[0].keys())

    try:
        with open(output_path, 'w', encoding='utf-8') as f:
            # Write header
            f.write('\t'.join(columns) + '\n')

            # Write data rows
            for row in results:
                values = [str(row.get(col, '')) for col in columns]
                f.write('\t'.join(values) + '\n')

        logger.info(f"TSV report written to: {output_path}")

    except OSError as e:
        raise OSError(f"Failed to write TSV file: {e}") from e


def write_json_report(
    results: List[Dict],
    output_path: str,
    include_summary: bool = True
) -> None:
    """
    Write QC results to a JSON file.

    Args:
        results: List of dictionaries containing QC results.
        output_path: Path for the output JSON file.
        include_summary: If True, include summary statistics.

    Raises:
        IOError: If file cannot be written.
    """
    output = {
        'results': results
    }

    if include_summary and results:
        passed = sum(1 for r in results if r.get('status') == 'PASS')
        failed = len(results) - passed
        output['summary'] = {
            'total_sequences': len(results),
            'passed': passed,
            'failed': failed,
            'pass_rate': round(passed / len(results) * 100, 2) if results else 0
        }

    try:
        with open(output_path, 'w', encoding='utf-8') as f:
            json.dump(output, f, indent=2)

        logger.info(f"JSON report written to: {output_path}")

    except OSError as e:
        raise OSError(f"Failed to write JSON file: {e}") from e


def write_fasta(
    sequences: List[Tuple[str, str]],
    output_path: str,
    line_width: int = 80
) -> None:
    """
    Write sequences to a FASTA file.

    Args:
        sequences: List of (header, sequence) tuples.
        output_path: Path for the output FASTA file.
        line_width: Maximum line width for sequences (default: 80).

    Raises:
        IOError: If file cannot be written.
    """
    try:
        with open(output_path, 'w', encoding='utf-8') as f:
            for header, sequence in sequences:
                f.write(f">{header}\n")

                # Wrap sequence to specified line width
                for i in range(0, len(sequence), line_width):
                    f.write(sequence[i:i + line_width] + '\n')

        logger.info(f"FASTA file written to: {output_path} ({len(sequences)} sequences)")

    except OSError as e:
        raise OSError(f"Failed to write FASTA file: {e}") from e
