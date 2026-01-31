"""
ViralSeq-QC: Test Suite Configuration

Pytest fixtures and shared test resources for ViralSeq-QC testing.
"""

import os
import tempfile

import pytest

# =============================================================================
# Test Data Fixtures
# =============================================================================

@pytest.fixture
def sample_sequences():
    """Provide a dictionary of test sequences with known properties."""
    return {
        'high_quality': "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC" * 10,  # 540bp, 50% GC
        'short_fragment': "ATGCATGC",  # 8bp - should fail length check
        'high_n_content': "ATGCNNNNNNNNNNNNNNNNNNNATGC",  # High N content
        'low_complexity': "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" * 5,  # Low complexity
        'homopolymer': "ATGC" + "A" * 20 + "TGCATGC" * 10,  # Contains homopolymer run
        'terminal_ns': "NNNNN" + "ATGCATGC" * 30 + "NNNNN",  # Terminal N's
        'empty': "",
        'mixed_case': "AtGcAtGcAtGcAtGcAtGcAtGcAtGc" * 10,
    }


@pytest.fixture
def high_quality_sequence():
    """A high-quality viral genome fragment that should pass all QC checks."""
    # 1000bp of realistic viral sequence (high complexity, low N, no long homopolymers)
    return (
        "ATGCAGTCGATCGATCGATCGTAGCTAGCTAGCTAGCTAGCATGCATGC"
        "GCTAGCTAGCTAGCATGCATGCATGCATGCGATCGATCGATCGATCGAT"
        "CGATCGATCGATCGATCGTAGCTAGCTAGCTAGCTAGCATGCATGCATG"
        "CATGCATGCATGCGATCGATCGATCGATCGATCGATCGATCGATCGTAG"
        "CTAGCTAGCTAGCTAGCATGCATGCATGCATGCATGCATGCGATCGATC"
        "GATCGATCGATCGATCGATCGATCGTAGCTAGCTAGCTAGCTAGCATGC"
        "ATGCATGCATGCATGCATGCGATCGATCGATCGATCGATCGATCGATCG"
        "ATCGTAGCTAGCTAGCTAGCTAGCATGCATGCATGCATGCATGCATGCG"
        "ATCGATCGATCGATCGATCGATCGATCGATCGTAGCTAGCTAGCTAGCT"
        "AGCATGCATGCATGCATGCATGCATGCGATCGATCGATCGATCGATCGA"
    )


@pytest.fixture
def temp_fasta_file(sample_sequences):
    """Create a temporary FASTA file with test sequences."""
    content = ""
    for name, seq in sample_sequences.items():
        if seq:  # Skip empty sequences
            content += f">{name}\n{seq}\n"

    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(content)
        temp_path = f.name

    yield temp_path

    # Cleanup
    if os.path.exists(temp_path):
        os.unlink(temp_path)


@pytest.fixture
def temp_multiline_fasta():
    """Create a temporary FASTA file with multiline sequences."""
    content = """>multiline_seq
ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
>single_line_seq
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT
"""

    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(content)
        temp_path = f.name

    yield temp_path

    # Cleanup
    if os.path.exists(temp_path):
        os.unlink(temp_path)


@pytest.fixture
def temp_output_dir():
    """Create a temporary directory for output files."""
    temp_dir = tempfile.mkdtemp()
    yield temp_dir

    # Cleanup
    import shutil
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)


@pytest.fixture
def empty_fasta_file():
    """Create an empty FASTA file."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        temp_path = f.name

    yield temp_path

    if os.path.exists(temp_path):
        os.unlink(temp_path)


@pytest.fixture
def malformed_fasta_file():
    """Create a malformed FASTA file (sequence before header)."""
    content = """ATGCATGCATGC
>proper_header
GCTAGCTAGCTA
"""

    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(content)
        temp_path = f.name

    yield temp_path

    if os.path.exists(temp_path):
        os.unlink(temp_path)
