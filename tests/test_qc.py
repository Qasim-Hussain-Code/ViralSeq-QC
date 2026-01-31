"""
ViralSeq-QC: Unit Tests for QC Module

Comprehensive tests for all quality control functions including:
- Core metrics (GC content, N content, length checks)
- Advanced metrics (homopolymers, complexity, terminal N's)
- Composite QC functions
"""

import pytest
from src.qc import (
    calculate_gc_content,
    calculate_n_content,
    check_length,
    is_high_quality,
    detect_homopolymer_runs,
    calculate_sequence_complexity,
    check_terminal_ns,
    validate_nucleotides,
    get_sequence_metrics,
    process_batch,
)


class TestCalculateGCContent:
    """Tests for calculate_gc_content function."""
    
    def test_balanced_sequence(self):
        """50% GC content sequence."""
        assert calculate_gc_content("ATGCATGC") == 50.0
    
    def test_high_gc_sequence(self):
        """100% GC content sequence."""
        assert calculate_gc_content("GGCCGGCC") == 100.0
    
    def test_zero_gc_sequence(self):
        """0% GC content sequence (only A and T)."""
        assert calculate_gc_content("ATATATATAT") == 0.0
    
    def test_empty_sequence(self):
        """Empty sequence should return 0.0."""
        assert calculate_gc_content("") == 0.0
    
    def test_case_insensitive(self):
        """Should handle mixed case sequences."""
        assert calculate_gc_content("atgcATGC") == 50.0
    
    def test_lowercase_sequence(self):
        """Should handle fully lowercase sequences."""
        assert calculate_gc_content("atgcatgc") == 50.0
    
    def test_single_base(self):
        """Single base sequences."""
        assert calculate_gc_content("G") == 100.0
        assert calculate_gc_content("A") == 0.0
        
    def test_sequence_with_n(self):
        """Sequence with N bases (N's not counted as GC)."""
        result = calculate_gc_content("ATGCNNNN")  # 2 GC out of 8
        assert result == 25.0


class TestCalculateNContent:
    """Tests for calculate_n_content function."""
    
    def test_no_n_content(self):
        """Sequence with no N's."""
        assert calculate_n_content("ATGCATGC") == 0.0
    
    def test_all_n_content(self):
        """Sequence that is all N's."""
        assert calculate_n_content("NNNNNNNN") == 100.0
    
    def test_partial_n_content(self):
        """Sequence with some N's."""
        result = calculate_n_content("ATGCNNNN")  # 4 N out of 8
        assert result == 50.0
    
    def test_empty_sequence(self):
        """Empty sequence should return 0.0."""
        assert calculate_n_content("") == 0.0
    
    def test_case_insensitive(self):
        """Should handle lowercase n's."""
        result = calculate_n_content("ATGCnnnn")
        assert result == 50.0


class TestCheckLength:
    """Tests for check_length function."""
    
    def test_sequence_at_threshold(self):
        """Sequence exactly at threshold should pass."""
        seq = "A" * 200
        assert check_length(seq, min_length=200) is True
    
    def test_sequence_above_threshold(self):
        """Sequence above threshold should pass."""
        seq = "A" * 500
        assert check_length(seq, min_length=200) is True
    
    def test_sequence_below_threshold(self):
        """Sequence below threshold should fail."""
        seq = "A" * 100
        assert check_length(seq, min_length=200) is False
    
    def test_empty_sequence(self):
        """Empty sequence should fail."""
        assert check_length("", min_length=200) is False
        assert check_length("", min_length=0) is False
    
    def test_custom_threshold(self):
        """Test with custom threshold values."""
        seq = "A" * 29000
        assert check_length(seq, min_length=29000) is True
        assert check_length(seq, min_length=30000) is False


class TestDetectHomopolymerRuns:
    """Tests for detect_homopolymer_runs function."""
    
    def test_long_homopolymer(self):
        """Detect a single long homopolymer."""
        seq = "ATGC" + "A" * 15 + "TGCA"
        runs = detect_homopolymer_runs(seq, min_length=10)
        assert len(runs) == 1
        assert runs[0]['base'] == 'A'
        assert runs[0]['length'] == 15
    
    def test_multiple_homopolymers(self):
        """Detect multiple homopolymer runs."""
        seq = "A" * 12 + "TGCA" + "T" * 12
        runs = detect_homopolymer_runs(seq, min_length=10)
        assert len(runs) == 2
    
    def test_no_homopolymers(self):
        """No homopolymers of sufficient length."""
        seq = "ATGCATGCATGCATGCATGC"
        runs = detect_homopolymer_runs(seq, min_length=10)
        assert len(runs) == 0
    
    def test_empty_sequence(self):
        """Empty sequence should return empty list."""
        assert detect_homopolymer_runs("", min_length=10) == []
    
    def test_n_homopolymer(self):
        """Detect N homopolymer runs."""
        seq = "ATGC" + "N" * 15 + "TGCA"
        runs = detect_homopolymer_runs(seq, min_length=10)
        assert len(runs) == 1
        assert runs[0]['base'] == 'N'


class TestCalculateSequenceComplexity:
    """Tests for calculate_sequence_complexity function."""
    
    def test_repetitive_sequence(self):
        """Highly repetitive sequence should have low complexity."""
        seq = "ATG" * 100  # Very repetitive
        complexity = calculate_sequence_complexity(seq, k=3)
        # Only 3 unique 3-mers: ATG, TGA, GAT
        assert complexity < 0.1
    
    def test_high_complexity_sequence(self):
        """Diverse sequence should have higher complexity."""
        seq = "ATGCAGTCGATCGTAGCTAG" * 10
        complexity = calculate_sequence_complexity(seq, k=3)
        assert complexity > 0.2
    
    def test_short_sequence(self):
        """Sequence shorter than k should return 0."""
        assert calculate_sequence_complexity("AT", k=3) == 0.0
    
    def test_empty_sequence(self):
        """Empty sequence should return 0."""
        assert calculate_sequence_complexity("", k=3) == 0.0


class TestCheckTerminalNs:
    """Tests for check_terminal_ns function."""
    
    def test_no_terminal_ns(self):
        """Sequence with no terminal N's."""
        passes, details = check_terminal_ns("ATGCATGCATGC", max_terminal_n=5)
        assert passes is True
        assert details['5_prime_ns'] == 0
        assert details['3_prime_ns'] == 0
    
    def test_acceptable_terminal_ns(self):
        """Sequence with acceptable terminal N count."""
        seq = "NNN" + "ATGCATGC" + "NN"
        passes, details = check_terminal_ns(seq, max_terminal_n=5)
        assert passes is True
        assert details['5_prime_ns'] == 3
        assert details['3_prime_ns'] == 2
    
    def test_excessive_terminal_ns(self):
        """Sequence with too many terminal N's should fail."""
        seq = "N" * 20 + "ATGCATGC" + "N" * 20
        passes, details = check_terminal_ns(seq, max_terminal_n=10)
        assert passes is False
        assert details['5_prime_ns'] == 20
        assert details['3_prime_ns'] == 20
    
    def test_empty_sequence(self):
        """Empty sequence should pass."""
        passes, details = check_terminal_ns("", max_terminal_n=5)
        assert passes is True


class TestValidateNucleotides:
    """Tests for validate_nucleotides function."""
    
    def test_valid_sequence(self):
        """Standard ATGCN sequence should be valid."""
        is_valid, details = validate_nucleotides("ATGCNATGC")
        assert is_valid is True
        assert details['valid'] is True
    
    def test_invalid_characters(self):
        """Sequence with invalid characters should fail."""
        is_valid, details = validate_nucleotides("ATGCXYZ")
        assert is_valid is False
        assert 'X' in details['invalid_chars']
        assert 'Y' in details['invalid_chars']
        assert 'Z' in details['invalid_chars']
    
    def test_empty_sequence(self):
        """Empty sequence should be valid."""
        is_valid, details = validate_nucleotides("")
        assert is_valid is True


class TestIsHighQuality:
    """Tests for is_high_quality composite function."""
    
    def test_high_quality_sequence(self, high_quality_sequence):
        """High quality sequence should pass."""
        assert is_high_quality(
            high_quality_sequence, 
            min_length=200, 
            max_n_content=5.0
        ) is True
    
    def test_fails_length(self):
        """Short sequence should fail."""
        assert is_high_quality("ATGC", min_length=200) is False
    
    def test_fails_n_content(self):
        """High N content sequence should fail."""
        seq = "ATGC" * 50 + "N" * 100  # > 5% N
        assert is_high_quality(seq, max_n_content=5.0) is False
    
    def test_empty_sequence(self):
        """Empty sequence should fail."""
        assert is_high_quality("") is False
    
    def test_custom_thresholds(self):
        """Test with viral genome thresholds (SARS-CoV-2 ~29.9kb)."""
        long_seq = "ATGC" * 7500  # 30,000 bp
        assert is_high_quality(long_seq, min_length=29000, max_n_content=1.0) is True


class TestGetSequenceMetrics:
    """Tests for get_sequence_metrics function."""
    
    def test_complete_metrics(self, high_quality_sequence):
        """Should return all expected metrics."""
        metrics = get_sequence_metrics(high_quality_sequence)
        
        assert 'length' in metrics
        assert 'gc_content' in metrics
        assert 'n_content' in metrics
        assert 'complexity' in metrics
        assert 'homopolymer_count' in metrics
        assert 'valid_nucleotides' in metrics
        assert '5_prime_ns' in metrics
        assert '3_prime_ns' in metrics
        
        assert metrics['length'] > 0
        assert 0 <= metrics['gc_content'] <= 100
        assert metrics['n_content'] == 0.0  # No N's in test sequence
    
    def test_empty_sequence_metrics(self):
        """Empty sequence should return zero metrics."""
        metrics = get_sequence_metrics("")
        assert metrics['length'] == 0
        assert metrics['gc_content'] == 0.0


class TestProcessBatch:
    """Tests for process_batch function."""
    
    def test_mixed_batch(self, sample_sequences):
        """Process batch with mixed quality sequences."""
        sequences = [
            sample_sequences['high_quality'],
            sample_sequences['short_fragment'],
            sample_sequences['mixed_case'],
        ]
        
        passed, failed = process_batch(sequences, min_length=100)
        
        assert len(passed) == 2  # high_quality and mixed_case
        assert len(failed) == 1  # short_fragment
    
    def test_all_pass(self):
        """Batch where all pass."""
        sequences = ["ATGC" * 100, "GCTA" * 100]
        passed, failed = process_batch(sequences, min_length=100)
        
        assert len(passed) == 2
        assert len(failed) == 0
    
    def test_all_fail(self):
        """Batch where all fail."""
        sequences = ["ATG", "GCT", ""]
        passed, failed = process_batch(sequences, min_length=100)
        
        assert len(passed) == 0
        assert len(failed) == 3
    
    def test_empty_batch(self):
        """Empty batch should return empty lists."""
        passed, failed = process_batch([])
        assert passed == []
        assert failed == []
