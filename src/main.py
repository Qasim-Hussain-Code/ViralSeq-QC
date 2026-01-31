import argparse
import sys
from src.input_output import parse_fasta
from src.qc import is_high_quality

def main():
    # 1. Setup the Argument Parser
    parser = argparse.ArgumentParser(
        description="ViralSeq-QC: A toolkit for viral genome quality control."
    )
    
    # 2. Define the rules (What flags do we accept?)
    parser.add_argument("--input", required=True, help="Path to input FASTA file")
    parser.add_argument("--output", required=True, help="Path to output QC report (TSV)")
    parser.add_argument("--min_length", type=int, default=200, help="Minimum sequence length (default: 200)")
    parser.add_argument("--max_n", type=float, default=5.0, help="Maximum % N-content (default: 5.0)")
    
    # 3. Read the user's orders
    args = parser.parse_args()
    
    print(f"--- Starting ViralSeq-QC ---")
    print(f"Processing: {args.input}")
    print(f"Filters: Length >= {args.min_length}, N-content <= {args.max_n}%")
    
    # 4. The Workflow Loop
    passed_count = 0
    total_count = 0
    
    try:
        with open(args.output, "w") as out_file:
            # Write Header
            out_file.write("Sequence_ID\tLength\tGC_Content\tStatus\n")
            
            # Stream the data
            for header, sequence in parse_fasta(args.input):
                total_count += 1
                
                # Ask the QC module for a verdict
                passed = is_high_quality(
                    sequence, 
                    min_length=args.min_length, 
                    max_n_content=args.max_n
                )
                
                # Write result logic
                status = "PASS" if passed else "FAIL"
                if passed:
                    passed_count += 1
                
                # Simple GC calc for report
                gc = (sequence.count("G") + sequence.count("C")) / len(sequence) * 100 if len(sequence) > 0 else 0
                
                out_file.write(f"{header}\t{len(sequence)}\t{gc:.2f}\t{status}\n")
                
    except FileNotFoundError:
        print(f"Error: Input file '{args.input}' not found.")
        sys.exit(1)
        
    # 5. Final Summary
    print(f"--- Complete ---")
    print(f"Total Sequences: {total_count}")
    print(f"Passed QC: {passed_count}")
    print(f"Report saved to: {args.output}")

if __name__ == "__main__":
    main()