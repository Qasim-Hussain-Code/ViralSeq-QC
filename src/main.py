#!/usr/bin/env python3
"""
ViralSeq-QC: Command-Line Interface

A zero-dependency CLI toolkit for viral genome quality control.
Supports multiple output formats and comprehensive QC metrics.

Usage:
    python -m src.main --input sequences.fasta --output qc_report.tsv
    viralseq-qc --input sequences.fasta --output qc_report.json --json
"""

import argparse
import logging
import sys
import json
from pathlib import Path
from typing import List, Dict, Tuple, Optional

from src.input_output import (
    parse_fasta, 
    write_tsv_report, 
    write_json_report, 
    write_fasta,
    FastaParseError
)
from src.qc import (
    is_high_quality, 
    calculate_gc_content, 
    calculate_n_content,
    get_sequence_metrics,
    detect_homopolymer_runs,
    calculate_sequence_complexity
)
from src import __version__


def setup_logging(verbose: bool = False) -> None:
    """Configure logging based on verbosity level."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


def create_parser() -> argparse.ArgumentParser:
    """Create and configure the argument parser."""
    parser = argparse.ArgumentParser(
        prog='viralseq-qc',
        description=(
            "ViralSeq-QC: A zero-dependency toolkit for viral genome quality control. "
            "Filters sequences based on length, N-content, and other quality metrics."
        ),
        epilog=(
            "Examples:\n"
            "  viralseq-qc --input raw.fasta --output qc_report.tsv\n"
            "  viralseq-qc --input raw.fasta.gz --output report.json --json\n"
            "  viralseq-qc --input raw.fasta --output report.tsv --fasta-out passed.fasta\n"
            "  viralseq-qc --input raw.fasta --output report.tsv --min-length 29000 --max-n 1.0"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Required arguments
    parser.add_argument(
        '--input', '-i',
        required=True, 
        help="Path to input FASTA file (supports .gz compression)"
    )
    parser.add_argument(
        '--output', '-o',
        required=True, 
        help="Path to output QC report (TSV by default)"
    )
    
    # QC threshold arguments
    qc_group = parser.add_argument_group('QC Thresholds')
    qc_group.add_argument(
        '--min-length', 
        type=int, 
        default=200, 
        metavar='BP',
        help="Minimum sequence length in base pairs (default: 200)"
    )
    qc_group.add_argument(
        '--max-n', 
        type=float, 
        default=5.0, 
        metavar='PERCENT',
        help="Maximum allowed N-content percentage (default: 5.0)"
    )
    qc_group.add_argument(
        '--min-complexity',
        type=float,
        default=0.0,
        metavar='SCORE',
        help="Minimum sequence complexity score 0-1 (default: 0 = disabled)"
    )
    
    # Output format arguments
    output_group = parser.add_argument_group('Output Options')
    output_group.add_argument(
        '--json', 
        action='store_true',
        help="Output report in JSON format instead of TSV"
    )
    output_group.add_argument(
        '--fasta-out',
        metavar='FILE',
        help="Write passing sequences to a new FASTA file"
    )
    output_group.add_argument(
        '--failed-out',
        metavar='FILE',
        help="Write failing sequences to a separate FASTA file"
    )
    output_group.add_argument(
        '--detailed',
        action='store_true',
        help="Include all QC metrics in the output report"
    )
    
    # Utility arguments
    util_group = parser.add_argument_group('Utility')
    util_group.add_argument(
        '--verbose', '-v',
        action='store_true',
        help="Enable verbose output with debug information"
    )
    util_group.add_argument(
        '--quiet', '-q',
        action='store_true',
        help="Suppress all output except errors"
    )
    util_group.add_argument(
        '--version',
        action='version',
        version=f'%(prog)s {__version__}'
    )
    
    return parser


def process_sequences(
    input_path: str,
    min_length: int,
    max_n: float,
    min_complexity: float,
    detailed: bool,
    logger: logging.Logger
) -> Tuple[List[Dict], List[Tuple[str, str]], List[Tuple[str, str]]]:
    """
    Process all sequences from the input file and generate QC results.
    
    Returns:
        Tuple of (results_list, passed_sequences, failed_sequences)
    """
    results = []
    passed_seqs = []
    failed_seqs = []
    
    for header, sequence in parse_fasta(input_path):
        # Calculate core metrics
        length = len(sequence)
        gc_content = calculate_gc_content(sequence)
        n_content = calculate_n_content(sequence)
        
        # Apply QC logic
        passed = is_high_quality(
            sequence, 
            min_length=min_length, 
            max_n_content=max_n,
            min_complexity=min_complexity
        )
        status = "PASS" if passed else "FAIL"
        
        # Build result record
        result = {
            'sequence_id': header,
            'length': length,
            'gc_content': round(gc_content, 2),
            'n_content': round(n_content, 2),
            'status': status
        }
        
        # Add detailed metrics if requested
        if detailed:
            complexity = calculate_sequence_complexity(sequence)
            homopolymers = detect_homopolymer_runs(sequence)
            metrics = get_sequence_metrics(sequence)
            
            result['complexity'] = round(complexity, 3)
            result['homopolymer_count'] = len(homopolymers)
            result['5_prime_ns'] = metrics['5_prime_ns']
            result['3_prime_ns'] = metrics['3_prime_ns']
        
        results.append(result)
        
        # Track passed/failed sequences for optional FASTA output
        if passed:
            passed_seqs.append((header, sequence))
        else:
            failed_seqs.append((header, sequence))
        
        logger.debug(f"Processed: {header[:40]}... -> {status}")
    
    return results, passed_seqs, failed_seqs


def print_summary(
    results: List[Dict], 
    output_path: str,
    quiet: bool,
    logger: logging.Logger
) -> None:
    """Print a summary of QC results."""
    if quiet:
        return
    
    total = len(results)
    passed = sum(1 for r in results if r['status'] == 'PASS')
    failed = total - passed
    pass_rate = (passed / total * 100) if total > 0 else 0
    
    logger.info("=" * 50)
    logger.info("ViralSeq-QC Complete")
    logger.info("=" * 50)
    logger.info(f"Total Sequences:  {total}")
    logger.info(f"Passed QC:        {passed} ({pass_rate:.1f}%)")
    logger.info(f"Failed QC:        {failed} ({100-pass_rate:.1f}%)")
    logger.info(f"Report saved to:  {output_path}")
    logger.info("=" * 50)


def main() -> int:
    """
    Main entry point for ViralSeq-QC CLI.
    
    Returns:
        Exit code (0 for success, 1 for error).
    """
    parser = create_parser()
    args = parser.parse_args()
    
    # Setup logging
    if args.quiet:
        logging.basicConfig(level=logging.ERROR)
    else:
        setup_logging(verbose=args.verbose)
    
    logger = logging.getLogger(__name__)
    
    # Validate input file exists
    input_path = Path(args.input)
    if not input_path.exists():
        logger.error(f"Input file not found: {args.input}")
        return 1
    
    logger.info(f"Starting ViralSeq-QC v{__version__}")
    logger.info(f"Input: {args.input}")
    logger.info(f"Thresholds: length >= {args.min_length} bp, N-content <= {args.max_n}%")
    
    if args.min_complexity > 0:
        logger.info(f"Complexity filter: >= {args.min_complexity}")
    
    try:
        # Process sequences
        results, passed_seqs, failed_seqs = process_sequences(
            args.input,
            args.min_length,
            args.max_n,
            args.min_complexity,
            args.detailed,
            logger
        )
        
        if not results:
            logger.warning("No sequences found in input file")
            return 1
        
        # Write main output report
        if args.json:
            write_json_report(results, args.output, include_summary=True)
        else:
            columns = ['sequence_id', 'length', 'gc_content', 'n_content', 'status']
            if args.detailed:
                columns.extend(['complexity', 'homopolymer_count', '5_prime_ns', '3_prime_ns'])
            write_tsv_report(results, args.output, columns=columns)
        
        # Write optional FASTA outputs
        if args.fasta_out and passed_seqs:
            write_fasta(passed_seqs, args.fasta_out)
            logger.info(f"Passed sequences saved to: {args.fasta_out}")
        
        if args.failed_out and failed_seqs:
            write_fasta(failed_seqs, args.failed_out)
            logger.info(f"Failed sequences saved to: {args.failed_out}")
        
        # Print summary
        print_summary(results, args.output, args.quiet, logger)
        
        return 0
        
    except FastaParseError as e:
        logger.error(f"FASTA parsing error: {e}")
        return 1
    except IOError as e:
        logger.error(f"I/O error: {e}")
        return 1
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())