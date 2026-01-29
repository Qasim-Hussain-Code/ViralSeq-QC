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