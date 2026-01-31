"""
ViralSeq-QC: Quality Control Module

This module provides comprehensive quality control functions for viral genome sequences,
including standard metrics (GC content, length, N-content) and advanced bioinformatics
metrics (homopolymer detection, sequence complexity, terminal N analysis).
"""

import re
from collections import Counter
from typing import Dict, List, Tuple

# =============================================================================
# Core QC Metrics
# =============================================================================

def calculate_gc_content(sequence: str) -> float:
    """
    Calculate the GC content percentage of a viral sequence.

    GC content is an important quality metric that can indicate:
    - Sequencing bias (extreme values may suggest contamination)
    - Taxonomic classification (different viral families have characteristic GC%)

    Args:
        sequence: The DNA sequence (case-insensitive).

    Returns:
        GC percentage (0-100). Returns 0.0 for empty sequences.

    Example:
        >>> calculate_gc_content("ATGCATGC")
        50.0
    """
    if not sequence:
        return 0.0

    seq_upper = sequence.upper()
    g_count = seq_upper.count('G')
    c_count = seq_upper.count('C')
    total_len = len(seq_upper)

    if total_len == 0:
        return 0.0

    return ((g_count + c_count) / total_len) * 100.0


def calculate_n_content(sequence: str) -> float:
    """
    Calculate the percentage of ambiguous bases (N) in a sequence.

    High N-content indicates:
    - Low coverage regions during sequencing
    - Assembly failures or masked regions
    - Poor sequence quality

    Args:
        sequence: The DNA sequence (case-insensitive).

    Returns:
        N percentage (0-100). Returns 0.0 for empty sequences.

    Example:
        >>> calculate_n_content("ATGNNATGC")
        22.22
    """
    if not sequence:
        return 0.0

    seq_upper = sequence.upper()
    n_count = seq_upper.count('N')
    total_len = len(sequence)

    if total_len == 0:
        return 0.0

    return (n_count / total_len) * 100.0


def check_length(sequence: str, min_length: int = 200) -> bool:
    """
    Validate if a sequence meets the minimum length threshold.

    Short sequences may indicate:
    - Fragmented assembly artifacts
    - Primer dimers or adapter contamination
    - Incomplete viral genomes

    Args:
        sequence: The DNA sequence.
        min_length: Minimum acceptable length in base pairs.

    Returns:
        True if sequence length >= min_length, False otherwise.
    """
    if not sequence:
        return False
    return len(sequence) >= min_length


# =============================================================================
# Advanced Bioinformatics Metrics
# =============================================================================

def detect_homopolymer_runs(sequence: str, min_length: int = 10) -> List[Dict]:
    """
    Detect potentially problematic homopolymer stretches in a sequence.

    Long homopolymer runs can indicate:
    - Sequencing errors (especially with nanopore/PacBio)
    - Assembly artifacts
    - Real biological features (rare in coding regions)

    Args:
        sequence: The DNA sequence (case-insensitive).
        min_length: Minimum homopolymer length to report.

    Returns:
        List of dictionaries with homopolymer details:
        - 'base': The repeated nucleotide
        - 'length': Length of the homopolymer run
        - 'start': Start position (0-indexed)
        - 'end': End position (0-indexed, exclusive)

    Example:
        >>> detect_homopolymer_runs("ATGAAAAAAAATGC", min_length=5)
        [{'base': 'A', 'length': 9, 'start': 3, 'end': 12}]
    """
    if not sequence:
        return []

    seq_upper = sequence.upper()
    homopolymers = []

    # Regex pattern to find runs of the same nucleotide
    pattern = re.compile(r'([ATGCN])\1{' + str(min_length - 1) + r',}')

    for match in pattern.finditer(seq_upper):
        homopolymers.append({
            'base': match.group(1),
            'length': len(match.group(0)),
            'start': match.start(),
            'end': match.end()
        })

    return homopolymers


def calculate_sequence_complexity(sequence: str, k: int = 3) -> float:
    """
    Calculate linguistic complexity using k-mer diversity.

    Low complexity sequences may indicate:
    - Repetitive elements
    - Sequencing/assembly artifacts
    - Low-complexity regions that should be masked

    Complexity is calculated as: (observed k-mers) / (maximum possible k-mers)

    Args:
        sequence: The DNA sequence (case-insensitive).
        k: K-mer size for complexity calculation (default: 3 for trinucleotides).

    Returns:
        Complexity score (0.0 to 1.0). Higher values indicate more complex sequences.
        Returns 0.0 for sequences shorter than k.

    Example:
        >>> calculate_sequence_complexity("ATGATGATGATG", k=3)  # Repetitive
        0.25
        >>> calculate_sequence_complexity("ATGCATGCNNNN", k=3)  # More diverse
        0.7
    """
    if not sequence or len(sequence) < k:
        return 0.0

    seq_upper = sequence.upper()

    # Count all k-mers
    kmers = [seq_upper[i:i+k] for i in range(len(seq_upper) - k + 1)]
    unique_kmers = len(set(kmers))

    # Maximum possible k-mers (limited by sequence length or 4^k for DNA)
    max_possible = min(len(kmers), 4 ** k)

    if max_possible == 0:
        return 0.0

    return unique_kmers / max_possible


def check_terminal_ns(sequence: str, max_terminal_n: int = 10) -> Tuple[bool, Dict]:
    """
    Check for excessive N's at sequence termini (common assembly artifact).

    Terminal N stretches often indicate:
    - Low coverage at contig ends
    - Assembly graph issues
    - Primer binding regions with poor sequencing

    Args:
        sequence: The DNA sequence (case-insensitive).
        max_terminal_n: Maximum allowed N's at each terminus.

    Returns:
        Tuple of (passes_check, details_dict):
        - passes_check: True if terminal N content is acceptable
        - details_dict: Contains '5_prime_ns' and '3_prime_ns' counts

    Example:
        >>> check_terminal_ns("NNNATGCATGCNNN", max_terminal_n=5)
        (True, {'5_prime_ns': 3, '3_prime_ns': 3})
    """
    if not sequence:
        return True, {'5_prime_ns': 0, '3_prime_ns': 0}

    seq_upper = sequence.upper()

    # Count leading N's (5' end)
    five_prime_ns = 0
    for base in seq_upper:
        if base == 'N':
            five_prime_ns += 1
        else:
            break

    # Count trailing N's (3' end)
    three_prime_ns = 0
    for base in reversed(seq_upper):
        if base == 'N':
            three_prime_ns += 1
        else:
            break

    details = {
        '5_prime_ns': five_prime_ns,
        '3_prime_ns': three_prime_ns
    }

    passes = five_prime_ns <= max_terminal_n and three_prime_ns <= max_terminal_n
    return passes, details


def validate_nucleotides(sequence: str) -> Tuple[bool, Dict]:
    """
    Validate that sequence contains only valid nucleotide characters.

    Valid characters: A, T, G, C, N (and their lowercase equivalents)
    Standard IUPAC ambiguity codes are flagged as warnings.

    Args:
        sequence: The DNA sequence.

    Returns:
        Tuple of (is_valid, details_dict):
        - is_valid: True if only ATGCN characters present
        - details_dict: Contains character counts and any invalid characters

    Example:
        >>> validate_nucleotides("ATGCATGC")
        (True, {'valid': True, 'base_counts': {'A': 2, 'T': 2, 'G': 2, 'C': 2}})
    """
    if not sequence:
        return True, {'valid': True, 'base_counts': {}, 'invalid_chars': []}

    seq_upper = sequence.upper()
    valid_bases = set('ATGCN')

    base_counts = Counter(seq_upper)
    invalid_chars = [char for char in base_counts.keys() if char not in valid_bases]

    details = {
        'valid': len(invalid_chars) == 0,
        'base_counts': dict(base_counts),
        'invalid_chars': invalid_chars
    }

    return len(invalid_chars) == 0, details


# =============================================================================
# Composite QC Functions
# =============================================================================

def is_high_quality(
    sequence: str,
    min_length: int = 200,
    max_n_content: float = 5.0,
    min_complexity: float = 0.0,
    max_terminal_n: int = 50
) -> bool:
    """
    Determine if a sequence passes all QC thresholds.

    This is the primary QC function that combines multiple checks:
    1. Minimum sequence length
    2. Maximum N-content percentage
    3. Minimum sequence complexity (optional)
    4. Terminal N content (optional)

    Args:
        sequence: Viral DNA sequence.
        min_length: Minimum acceptable length (bp).
        max_n_content: Maximum allowed N percentage.
        min_complexity: Minimum complexity score (0-1, default 0 = disabled).
        max_terminal_n: Maximum N's at each terminus (default 50 = lenient).

    Returns:
        True if sequence passes all checks, False otherwise.
    """
    # Check 1: Length
    if not check_length(sequence, min_length):
        return False

    # Check 2: N-content (Ambiguity)
    n_percent = calculate_n_content(sequence)
    if n_percent > max_n_content:
        return False

    # Check 3: Sequence complexity (if enabled)
    if min_complexity > 0:
        complexity = calculate_sequence_complexity(sequence)
        if complexity < min_complexity:
            return False

    # Check 4: Terminal N's (if threshold is set)
    if max_terminal_n < 50:  # Only check if not using default lenient value
        passes_terminal, _ = check_terminal_ns(sequence, max_terminal_n)
        if not passes_terminal:
            return False

    return True


def get_sequence_metrics(sequence: str) -> Dict:
    """
    Calculate all QC metrics for a sequence in a single call.

    This is useful for generating comprehensive QC reports.

    Args:
        sequence: The DNA sequence.

    Returns:
        Dictionary containing all calculated metrics.
    """
    if not sequence:
        return {
            'length': 0,
            'gc_content': 0.0,
            'n_content': 0.0,
            'complexity': 0.0,
            'homopolymer_count': 0,
            'valid_nucleotides': True,
            '5_prime_ns': 0,
            '3_prime_ns': 0
        }

    _, terminal_details = check_terminal_ns(sequence)
    _, nucleotide_details = validate_nucleotides(sequence)

    return {
        'length': len(sequence),
        'gc_content': calculate_gc_content(sequence),
        'n_content': calculate_n_content(sequence),
        'complexity': calculate_sequence_complexity(sequence),
        'homopolymer_count': len(detect_homopolymer_runs(sequence)),
        'valid_nucleotides': nucleotide_details['valid'],
        '5_prime_ns': terminal_details['5_prime_ns'],
        '3_prime_ns': terminal_details['3_prime_ns']
    }


def process_batch(
    sequences: List[str],
    min_length: int = 200,
    max_n_content: float = 5.0
) -> Tuple[List[str], List[str]]:
    """
    Process a batch of sequences and separate them into passed and failed lists.

    Args:
        sequences: A list of raw viral sequences.
        min_length: Threshold for length check.
        max_n_content: Maximum allowed N percentage.

    Returns:
        Tuple of (passed_sequences, failed_sequences)
    """
    passed = []
    failed = []

    for seq in sequences:
        if is_high_quality(seq, min_length=min_length, max_n_content=max_n_content):
            passed.append(seq)
        else:
            failed.append(seq)

    return passed, failed
