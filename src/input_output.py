from typing import Generator, Tuple

def parse_fasta(file_path: str) -> Generator[Tuple[str, str], None, None]:
    """
    Reads a FASTA file and yields (header, sequence) tuples one by one.
    Handles multiline sequences using a buffer strategy.
    
    Args:
        file_path (str): Path to the input FASTA file.
        
    Yields:
        Tuple[str, str]: (Header, DNA_Sequence)
    """
    header = None
    sequence_buffer = []

    try:
        with open(file_path, 'r') as file:
            for line in file:
                line = line.strip()
                if not line:
                    continue  # Skip empty lines

                if line.startswith(">"):
                    # If we have a previous record, yield it now
                    if header:
                        full_sequence = "".join(sequence_buffer)
                        yield header, full_sequence
                    
                    # Reset for the new record
                    header = line[1:]  # Remove the '>'
                    sequence_buffer = []
                else:
                    # Accumulate DNA lines
                    sequence_buffer.append(line)

            # yield the last record (The "Orphan" logic)
            if header:
                full_sequence = "".join(sequence_buffer)
                yield header, full_sequence

    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.")