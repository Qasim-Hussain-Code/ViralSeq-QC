from typing import Optional

def calculate_gc_content(sequence: str) -> float:
    """
    Calculates the GC content percentage of a viral sequence.
    
    Args:
        sequence (str): The DNA sequence (case-insensitive).
        
    Returns:
        float: GC percentage (0-100). Returns 0.0 for empty sequences.
    """
    if not sequence:
        return 0.0
        
    # Standardize to uppercase to handle mixed input
    seq_upper = sequence.upper()
    
    # Calculate counts 
    g_count = seq_upper.count('G')
    c_count = seq_upper.count('C')
    total_len = len(seq_upper)
    
    if total_len == 0:
        return 0.0
        
    return ((g_count + c_count) / total_len) * 100.0

def check_length(sequence: str, min_length: int = 200) -> bool:
    """
    Validates if a sequence meets the minimum length threshold.
    Useful for filtering out fragmented assembly artifacts.
    """
    if not sequence:
        return False
    return len(sequence) >= min_length

def is_high_quality(sequence: str, min_length: int = 200, max_n_content: float = 5.0) -> bool:
    """
    Determines if a sequence passes all QC thresholds.
    
    Args:
        sequence (str): Viral DNA sequence.
        min_length (int): Minimum acceptable length (bp).
        max_n_content (float): Maximum allowed N percentage.
        
    Returns:
        bool: True if sequence passes all checks, False otherwise.
    """
    # Check 1: Length
    if not check_length(sequence, min_length):
        return False
        
    # Check 2: N-content (Ambiguity)
    seq_upper = sequence.upper()
    n_count = seq_upper.count('N')
    total_len = len(sequence)
    
    # Avoid division by zero
    if total_len == 0:
        return False
        
    n_percent = (n_count / total_len) * 100.0
    
    # Decision logic
    if n_percent > max_n_content:
        return False
        
    return True

from typing import List, Tuple

def process_batch(sequences: List[str], min_length: int = 200) -> Tuple[List[str], List[str]]:
    """
    Process a batch of sequences and separate them into passed and failed lists.
    
    Args:
        sequences (List[str]): A list of raw viral sequences.
        min_length (int): Threshold for length check.
        
    Returns:
        Tuple[List[str], List[str]]: (passed_sequences, failed_sequences)
    """
    passed = []
    failed = []
    
    # Loop through every sequence in the input list
    for seq in sequences:
        # Apply the decision logic
        if is_high_quality(seq, min_length=min_length):
            passed.append(seq)
        else:
            failed.append(seq)
            
    return passed, failed