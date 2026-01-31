"""
ViralSeq-QC: Unit Tests for Input/Output Module

Comprehensive tests for FASTA parsing and output writing functions.
"""

import gzip
import json
import os
import tempfile

import pytest

from src.input_output import (
    FastaParseError,
    parse_fasta,
    validate_fasta_format,
    write_fasta,
    write_json_report,
    write_tsv_report,
)


class TestParseFasta:
    """Tests for parse_fasta function."""

    def test_parse_standard_fasta(self, temp_fasta_file):
        """Parse a standard FASTA file."""
        records = list(parse_fasta(temp_fasta_file))
        assert len(records) > 0

        # Check that each record is a tuple of (header, sequence)
        for header, sequence in records:
            assert isinstance(header, str)
            assert isinstance(sequence, str)
            assert len(header) > 0

    def test_parse_multiline_sequences(self, temp_multiline_fasta):
        """Parse FASTA with multiline sequences."""
        records = list(parse_fasta(temp_multiline_fasta))

        # Should correctly join multiline sequences
        assert len(records) == 2

        # First record should have all 3 lines joined
        header, seq = records[0]
        assert header == "multiline_seq"
        assert len(seq) == 76 * 3  # 3 lines of 76 characters

    def test_file_not_found(self):
        """Should raise FastaParseError for missing file."""
        with pytest.raises(FastaParseError):
            list(parse_fasta("nonexistent_file.fasta"))

    def test_empty_file(self, empty_fasta_file):
        """Empty file should yield no records."""
        records = list(parse_fasta(empty_fasta_file))
        assert len(records) == 0

    def test_malformed_fasta(self, malformed_fasta_file):
        """Malformed FASTA should raise FastaParseError."""
        with pytest.raises(FastaParseError):
            list(parse_fasta(malformed_fasta_file))


class TestParseFastaGzip:
    """Tests for gzip-compressed FASTA parsing."""

    def test_parse_gzipped_fasta(self):
        """Parse a gzip-compressed FASTA file."""
        content = b">test_seq\nATGCATGCATGC\n"

        with tempfile.NamedTemporaryFile(suffix='.fasta.gz', delete=False) as f:
            with gzip.open(f.name, 'wb') as gz:
                gz.write(content)
            temp_path = f.name

        try:
            records = list(parse_fasta(temp_path))
            assert len(records) == 1
            assert records[0][0] == "test_seq"
            assert records[0][1] == "ATGCATGCATGC"
        finally:
            os.unlink(temp_path)

    def test_invalid_gzip_file(self):
        """Invalid gzip file should raise FastaParseError."""
        with tempfile.NamedTemporaryFile(suffix='.gz', delete=False) as f:
            f.write(b"not a gzip file")
            temp_path = f.name

        try:
            with pytest.raises(FastaParseError):
                list(parse_fasta(temp_path))
        finally:
            os.unlink(temp_path)


class TestValidateFastaFormat:
    """Tests for validate_fasta_format function."""

    def test_valid_fasta(self, temp_fasta_file):
        """Valid FASTA should return True."""
        is_valid, message = validate_fasta_format(temp_fasta_file)
        assert is_valid is True
        assert "Valid FASTA" in message

    def test_invalid_fasta(self, malformed_fasta_file):
        """Invalid FASTA should return False with error message."""
        is_valid, message = validate_fasta_format(malformed_fasta_file)
        assert is_valid is False


class TestWriteTsvReport:
    """Tests for write_tsv_report function."""

    def test_write_basic_report(self, temp_output_dir):
        """Write a basic TSV report."""
        output_path = os.path.join(temp_output_dir, "report.tsv")

        results = [
            {'sequence_id': 'seq1', 'length': 1000, 'status': 'PASS'},
            {'sequence_id': 'seq2', 'length': 500, 'status': 'FAIL'},
        ]

        write_tsv_report(results, output_path)

        assert os.path.exists(output_path)

        with open(output_path) as f:
            lines = f.readlines()

        assert len(lines) == 3  # Header + 2 data rows
        assert 'sequence_id' in lines[0]
        assert 'seq1' in lines[1]
        assert 'seq2' in lines[2]

    def test_write_with_custom_columns(self, temp_output_dir):
        """Write TSV with specific columns."""
        output_path = os.path.join(temp_output_dir, "report.tsv")

        results = [
            {'sequence_id': 'seq1', 'length': 1000, 'gc': 50.0, 'status': 'PASS'},
        ]

        write_tsv_report(results, output_path, columns=['sequence_id', 'status'])

        with open(output_path) as f:
            header = f.readline()

        assert 'sequence_id' in header
        assert 'status' in header
        assert 'length' not in header  # Excluded column

    def test_write_empty_results(self, temp_output_dir):
        """Empty results should not create file content."""
        output_path = os.path.join(temp_output_dir, "report.tsv")
        write_tsv_report([], output_path)
        # Function should handle gracefully (warning logged)


class TestWriteJsonReport:
    """Tests for write_json_report function."""

    def test_write_json_with_summary(self, temp_output_dir):
        """Write JSON report with summary statistics."""
        output_path = os.path.join(temp_output_dir, "report.json")

        results = [
            {'sequence_id': 'seq1', 'status': 'PASS'},
            {'sequence_id': 'seq2', 'status': 'PASS'},
            {'sequence_id': 'seq3', 'status': 'FAIL'},
        ]

        write_json_report(results, output_path, include_summary=True)

        assert os.path.exists(output_path)

        with open(output_path) as f:
            data = json.load(f)

        assert 'results' in data
        assert 'summary' in data
        assert data['summary']['total_sequences'] == 3
        assert data['summary']['passed'] == 2
        assert data['summary']['failed'] == 1
        assert data['summary']['pass_rate'] == pytest.approx(66.67, rel=0.01)

    def test_write_json_without_summary(self, temp_output_dir):
        """Write JSON report without summary."""
        output_path = os.path.join(temp_output_dir, "report.json")

        results = [{'sequence_id': 'seq1', 'status': 'PASS'}]

        write_json_report(results, output_path, include_summary=False)

        with open(output_path) as f:
            data = json.load(f)

        assert 'results' in data
        assert 'summary' not in data


class TestWriteFasta:
    """Tests for write_fasta function."""

    def test_write_fasta_file(self, temp_output_dir):
        """Write sequences to FASTA format."""
        output_path = os.path.join(temp_output_dir, "output.fasta")

        sequences = [
            ('seq1', 'ATGCATGCATGC'),
            ('seq2', 'GCTAGCTAGCTA'),
        ]

        write_fasta(sequences, output_path)

        assert os.path.exists(output_path)

        # Read back and verify
        records = list(parse_fasta(output_path))
        assert len(records) == 2
        assert records[0] == ('seq1', 'ATGCATGCATGC')
        assert records[1] == ('seq2', 'GCTAGCTAGCTA')

    def test_write_fasta_with_line_wrap(self, temp_output_dir):
        """Long sequences should be wrapped."""
        output_path = os.path.join(temp_output_dir, "output.fasta")

        long_seq = "ATGC" * 100  # 400 bp
        sequences = [('long_seq', long_seq)]

        write_fasta(sequences, output_path, line_width=80)

        with open(output_path) as f:
            lines = f.readlines()

        # Header + 5 wrapped lines (400 / 80 = 5)
        assert len(lines) == 6

        # Verify sequence is recoverable
        records = list(parse_fasta(output_path))
        assert records[0][1] == long_seq

    def test_write_empty_fasta(self, temp_output_dir):
        """Empty sequence list should create empty file."""
        output_path = os.path.join(temp_output_dir, "output.fasta")
        write_fasta([], output_path)

        assert os.path.exists(output_path)
        with open(output_path) as f:
            content = f.read()
        assert content == ""
