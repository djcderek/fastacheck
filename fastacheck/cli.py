"""
FASTA Check - Command Line Interface
A tool for validating and analyzing FASTA files
"""

import argparse
import sys
import json
from pathlib import Path
from typing import Optional

from .parser import FastaParser
from .stats import BasicStats, AdvancedStats


def create_parser() -> argparse.ArgumentParser:
    """Create and configure the argument parser."""
    parser = argparse.ArgumentParser(
        prog='fastacheck',
        description='Validate and analyze FASTA files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  fastacheck validate input.fasta
  fastacheck analyze input.fasta --stats
  fastacheck validate input.fasta.gz --streaming
  fastacheck analyze input.fasta --output results.json
        """
    )
    
    parser.add_argument(
        '--version',
        action='version',
        version='%(prog)s 1.0.0'
    )
    
    # Subcommands
    subparsers = parser.add_subparsers(
        dest='command',
        help='Available commands',
        required=True
    )
    
    # Validate command
    validate_parser = subparsers.add_parser(
        'validate',
        help='Validate FASTA file format'
    )
    validate_parser.add_argument(
        'file',
        help='FASTA file to validate'
    )
    validate_parser.add_argument(
        '--streaming',
        action='store_true',
        help='Use streaming mode for large files'
    )
    validate_parser.add_argument(
        '--quiet',
        action='store_true',
        help='Only show validation result (pass/fail)'
    )
    
    # Analyze command
    analyze_parser = subparsers.add_parser(
        'analyze',
        help='Analyze FASTA file and show statistics'
    )
    analyze_parser.add_argument(
        'file',
        help='FASTA file to analyze'
    )
    analyze_parser.add_argument(
        '--streaming',
        action='store_true',
        help='Use streaming mode for large files'
    )
    analyze_parser.add_argument(
        '--detailed',
        action='store_true',
        help='Show detailed length and GC statistics'
    )
    analyze_parser.add_argument(
        '--assembly',
        action='store_true',
        help='Show genome assembly metrics (N50, auN, etc.)'
    )
    analyze_parser.add_argument(
        '--gene-set',
        action='store_true',
        help='Show gene/protein set metrics'
    )
    analyze_parser.add_argument(
        '--outliers',
        choices=['iqr', 'zscore'],
        help='Detect outlier sequences using specified method'
    )
    analyze_parser.add_argument(
        '--threshold',
        type=float,
        default=1.5,
        help='Threshold for outlier detection (default: 1.5)'
    )
    analyze_parser.add_argument(
        '--output',
        help='Output results to file'
    )
    analyze_parser.add_argument(
        '--format',
        choices=['json', 'txt', 'csv'],
        default='txt',
        help='Output format (default: txt)'
    )
    
    return parser


def validate_file(file_path: str, streaming: bool = False, quiet: bool = False) -> int:
    """
    Validate a FASTA file.
    
    Returns:
        0 if valid, 1 if invalid, 2 if error
    """
    try:
        parser = FastaParser(file_path, streaming=streaming)
        
        if not quiet:
            print(f"Validating: {file_path}")
        
        validation_result = parser.validate_format()
        
        if validation_result['is_valid']:
            if not quiet:
                print(f"✓ Valid FASTA file ({validation_result['sequence_count']} sequences)")
                
                if validation_result['warnings']:
                    print("\nWarnings:")
                    for warning in validation_result['warnings']:
                        print(f"  ⚠ {warning}")
            else:
                print("VALID")
            return 0
        else:
            if not quiet:
                print("✗ Invalid FASTA file")
                print("\nErrors:")
                for error in validation_result['errors']:
                    print(f"  ✗ {error}")
                    
                if validation_result['warnings']:
                    print("\nWarnings:")
                    for warning in validation_result['warnings']:
                        print(f"  ⚠ {warning}")
            else:
                print("INVALID")
            return 1
            
    except FileNotFoundError:
        print(f"Error: File not found: {file_path}", file=sys.stderr)
        return 2
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return 2


def analyze_file(file_path: str, streaming: bool = False, detailed: bool = False,
                assembly: bool = False, gene_set: bool = False, outliers: Optional[str] = None,
                threshold: float = 1.5, output_file: Optional[str] = None, 
                output_format: str = 'txt') -> int:
    """
    Analyze a FASTA file and show statistics.
    
    Returns:
        0 if successful, 1 if invalid, 2 if error
    """
    try:
        parser = FastaParser(file_path, streaming=streaming)
        
        print(f"Analyzing: {file_path}")
        
        # First validate
        validation_result = parser.validate_format()
        
        if not validation_result['is_valid']:
            print("✗ Invalid FASTA file - cannot analyze")
            print("\nErrors:")
            for error in validation_result['errors']:
                print(f"  ✗ {error}")
            return 1
        
        # Parse sequences and build statistics
        basic_stats = BasicStats()
        
        print("Processing sequences...")
        for header, sequence in parser.parse_sequences():
            # Extract header string for BasicStats
            header_str = header.get('sequence_id', '')
            basic_stats.add_sequence(sequence, header_str)
        
        # Get basic summary
        summary = basic_stats.get_summary()
        
        # Basic results
        print(f"\n✓ Analysis complete")
        print(f"Sequences: {summary['basic']['sequence_count']:,}")
        print(f"Total length: {summary['basic']['total_length']:,} bp")
        print(f"Average length: {summary['basic']['mean_length']:.1f} bp")
        print(f"Shortest: {summary['basic']['min_length']:,} bp")
        print(f"Longest: {summary['basic']['max_length']:,} bp")
        
        # Assembly metrics (N50, etc.)
        if summary['assembly_stats']['n50'] > 0:
            print(f"N50: {summary['assembly_stats']['n50']:,} bp")
            print(f"L50: {summary['assembly_stats']['l50']:,} sequences")
        
        # GC content
        if summary['gc_stats']:
            print(f"Overall GC content: {summary['gc_stats']['overall_gc_percent']:.1f}%")
        
        # Show detailed statistics if requested
        if detailed:
            print(f"\n=== Detailed Length Statistics ===")
            length_stats = summary['length_stats']
            print(f"Median length: {length_stats.get('median_length', 0):,.0f} bp")
            print(f"Standard deviation: {length_stats.get('std_dev', 0):,.1f}")
            print(f"Q1 (25th percentile): {length_stats.get('q1', 0):,.0f} bp")
            print(f"Q3 (75th percentile): {length_stats.get('q3', 0):,.0f} bp")
            
            if summary['gc_stats']:
                print(f"\n=== GC Content Statistics ===")
                gc_stats = summary['gc_stats']
                print(f"Mean GC content: {gc_stats.get('mean_gc_percent', 0):.1f}%")
                print(f"Median GC content: {gc_stats.get('median_gc_percent', 0):.1f}%")
                print(f"GC range: {gc_stats.get('min_gc_percent', 0):.1f}% - {gc_stats.get('max_gc_percent', 0):.1f}%")
                print(f"GC std dev: {gc_stats.get('std_dev_gc', 0):.1f}%")
            
            if summary['n_stats'] and summary['n_stats'].get('total_n_bases', 0) > 0:
                print(f"\n=== N Base Statistics ===")
                n_stats = summary['n_stats']
                print(f"Total N bases: {n_stats.get('total_n_bases', 0):,}")
                print(f"Sequences with Ns: {n_stats.get('sequences_with_n', 0):,} ({n_stats.get('percent_sequences_with_n', 0):.1f}%)")
                print(f"Overall N content: {n_stats.get('overall_n_percent', 0):.2f}%")
        
        # Assembly-specific metrics
        if assembly:
            advanced_stats = AdvancedStats(basic_stats)
            assembly_metrics = advanced_stats.calculate_genome_assembly_metrics()
            
            print(f"\n=== Genome Assembly Metrics ===")
            print(f"N25: {assembly_metrics.get('n25', 0):,} bp (L25: {assembly_metrics.get('l25', 0):,})")
            print(f"N50: {assembly_metrics.get('n50', 0):,} bp (L50: {assembly_metrics.get('l50', 0):,})")
            print(f"N75: {assembly_metrics.get('n75', 0):,} bp (L75: {assembly_metrics.get('l75', 0):,})")
            print(f"N90: {assembly_metrics.get('n90', 0):,} bp (L90: {assembly_metrics.get('l90', 0):,})")
            print(f"N95: {assembly_metrics.get('n95', 0):,} bp (L95: {assembly_metrics.get('l95', 0):,})")
            print(f"auN (Area under Nx curve): {assembly_metrics.get('auN', 0):,.1f}")
            
            print(f"\n=== Contig Size Distribution ===")
            print(f"Large contigs (≥10kb): {assembly_metrics.get('large_contigs', 0):,}")
            print(f"Medium contigs (1-10kb): {assembly_metrics.get('medium_contigs', 0):,}")
            print(f"Small contigs (100bp-1kb): {assembly_metrics.get('small_contigs', 0):,}")
            print(f"Very small contigs (<100bp): {assembly_metrics.get('very_small_contigs', 0):,}")
        
        # Gene set metrics
        if gene_set:
            advanced_stats = AdvancedStats(basic_stats)
            gene_metrics = advanced_stats.calculate_gene_set_metrics()
            
            print(f"\n=== Gene/Protein Set Metrics ===")
            print(f"Short sequences (<300bp): {gene_metrics.get('short_genes', 0):,}")
            print(f"Medium sequences (300-1500bp): {gene_metrics.get('medium_genes', 0):,}")
            print(f"Long sequences (≥1500bp): {gene_metrics.get('long_genes', 0):,}")
            print(f"Mean gene length: {gene_metrics.get('mean_gene_length', 0):,.1f} bp")
            print(f"Median gene length: {gene_metrics.get('median_gene_length', 0):,.1f} bp")
        
        # Outlier detection
        if outliers:
            outlier_sequences = basic_stats.get_outliers(method=outliers, threshold=threshold)
            
            print(f"\n=== Outlier Detection ({outliers.upper()}, threshold={threshold}) ===")
            if outlier_sequences:
                print(f"Found {len(outlier_sequences)} outlier sequences:")
                for idx, length in outlier_sequences[:10]:  # Show first 10
                    print(f"  Sequence {idx+1}: {length:,} bp")
                if len(outlier_sequences) > 10:
                    print(f"  ... and {len(outlier_sequences) - 10} more")
            else:
                print("No outlier sequences detected")
        
        # Handle output file
        if output_file:
            output_data = {
                'file_path': file_path,
                'validation': validation_result,
                'summary': summary
            }
            
            # Add advanced metrics if requested
            if assembly or gene_set:
                advanced_stats = AdvancedStats(basic_stats)
                if assembly:
                    output_data['assembly_metrics'] = advanced_stats.calculate_genome_assembly_metrics()
                if gene_set:
                    output_data['gene_metrics'] = advanced_stats.calculate_gene_set_metrics()
            
            # Add outliers if detected
            if outliers:
                output_data['outliers'] = {
                    'method': outliers,
                    'threshold': threshold,
                    'sequences': basic_stats.get_outliers(method=outliers, threshold=threshold)
                }
            
            with open(output_file, 'w') as f:
                if output_format == 'json':
                    json.dump(output_data, f, indent=2)
                elif output_format == 'csv':
                    # Simple CSV format for basic stats
                    f.write("metric,value\n")
                    f.write(f"sequence_count,{summary['basic']['sequence_count']}\n")
                    f.write(f"total_length,{summary['basic']['total_length']}\n")
                    f.write(f"mean_length,{summary['basic']['mean_length']:.1f}\n")
                    f.write(f"min_length,{summary['basic']['min_length']}\n")
                    f.write(f"max_length,{summary['basic']['max_length']}\n")
                    if summary['assembly_stats']['n50'] > 0:
                        f.write(f"n50,{summary['assembly_stats']['n50']}\n")
                        f.write(f"l50,{summary['assembly_stats']['l50']}\n")
                    if summary['gc_stats']:
                        f.write(f"gc_percent,{summary['gc_stats']['overall_gc_percent']:.1f}\n")
                else:
                    # Text format
                    f.write(f"FASTA Analysis Results\n")
                    f.write(f"=======================\n")
                    f.write(f"File: {file_path}\n")
                    f.write(f"Sequences: {summary['basic']['sequence_count']:,}\n")
                    f.write(f"Total length: {summary['basic']['total_length']:,} bp\n")
                    f.write(f"Average length: {summary['basic']['mean_length']:.1f} bp\n")
                    f.write(f"N50: {summary['assembly_stats']['n50']:,} bp\n")
                    if summary['gc_stats']:
                        f.write(f"GC content: {summary['gc_stats']['overall_gc_percent']:.1f}%\n")
            
            print(f"\nResults saved to: {output_file}")
        
        return 0
        
    except FileNotFoundError:
        print(f"Error: File not found: {file_path}", file=sys.stderr)
        return 2
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return 2


def main():
    """Main CLI entry point."""
    parser = create_parser()
    args = parser.parse_args()
    
    if args.command == 'validate':
        return validate_file(
            args.file,
            streaming=args.streaming,
            quiet=args.quiet
        )
    
    elif args.command == 'analyze':
        return analyze_file(
            args.file,
            streaming=args.streaming,
            detailed=args.detailed,
            assembly=args.assembly,
            gene_set=getattr(args, 'gene_set', False),
            outliers=args.outliers,
            threshold=args.threshold,
            output_file=args.output,
            output_format=args.format
        )
    
    else:
        parser.print_help()
        return 1


if __name__ == '__main__':
    sys.exit(main())