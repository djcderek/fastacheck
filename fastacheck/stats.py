"""
FastaCheck Statistics Module

This module provides comprehensive statistical analysis for FASTA sequences,
including basic metrics, N50 calculations, and distribution analysis.
"""

import statistics
from typing import List, Dict, Tuple, Optional, Any
import numpy as np
from collections import Counter


class BasicStats:
    """
    Core statistics calculator for FASTA sequences.
    
    Calculates fundamental metrics like sequence counts, lengths, N50, and
    provides statistical summaries of sequence collections.
    """
    
    def __init__(self):
        """Initialize statistics tracking."""
        self.sequence_count = 0
        self.total_length = 0
        self.lengths = []
        self.min_length = float('inf')
        self.max_length = 0
        self.gc_counts = []
        self.n_counts = []
        self.headers = []
        
    def add_sequence(self, sequence: str, header: str = ""):
        """
        Add a sequence to the statistics calculation.
        
        Args:
            sequence: The DNA/RNA/protein sequence string
            header: Optional header information
        """
        length = len(sequence)
        
        self.sequence_count += 1
        self.total_length += length
        self.lengths.append(length)
        
        self.min_length = min(self.min_length, length)
        self.max_length = max(self.max_length, length)

        gc_count = sequence.upper().count('G') + sequence.upper().count('C')
        self.gc_counts.append(gc_count)
        
        n_count = sequence.upper().count('N')
        self.n_counts.append(n_count)
        
        if header:
            self.headers.append(header)
    
    def calculate_n50(self):
        """
        Calculate N50 and L50 values.
        
        N50: The length of the shortest sequence in the set of longest sequences
             that together represent at least 50% of the total length.
        L50: The number of sequences in the N50 set.
        
        Returns:
            Tuple of (N50 value, L50 count)
        """
        if not self.lengths:
            return 0, 0
            
        sorted_lengths = sorted(self.lengths, reverse=True)
        total_length = sum(sorted_lengths)
        target = total_length / 2
        
        running_sum = 0
        for i, length in enumerate(sorted_lengths):
            running_sum += length
            if running_sum >= target:
                return length, i + 1
        
        return sorted_lengths[-1], len(sorted_lengths)
    
    def calculate_nx(self, x: int):
        """
        Calculate Nx and Lx values for any percentage x.
        """
        if not self.lengths or x < 0 or x > 100:
            return 0, 0
            
        sorted_lengths = sorted(self.lengths, reverse=True)
        total_length = sum(sorted_lengths)
        target = total_length * (x / 100)
        
        running_sum = 0
        for i, length in enumerate(sorted_lengths):
            running_sum += length
            if running_sum >= target:
                return length, i + 1
        
        return sorted_lengths[-1], len(sorted_lengths)
    
    def get_length_statistics(self):
        """
        Calculate comprehensive length statistics.
        """
        if not self.lengths:
            return {}
            
        return {
            'count': len(self.lengths),
            'total_length': sum(self.lengths),
            'min_length': min(self.lengths),
            'max_length': max(self.lengths),
            'mean_length': statistics.mean(self.lengths),
            'median_length': statistics.median(self.lengths),
            'std_dev': statistics.stdev(self.lengths) if len(self.lengths) > 1 else 0,
            'variance': statistics.variance(self.lengths) if len(self.lengths) > 1 else 0,
            'q1': self._percentile(self.lengths, 25),
            'q3': self._percentile(self.lengths, 75)
        }
    
    def get_gc_statistics(self):
        """
        Calculate GC content statistics.
        """
        if not self.gc_counts or not self.lengths:
            return {}
        
        # Calculate GC percentages
        gc_percentages = [
            (gc_count / length) * 100 if length > 0 else 0
            for gc_count, length in zip(self.gc_counts, self.lengths)
        ]
        
        total_gc = sum(self.gc_counts)
        overall_gc = (total_gc / self.total_length) * 100 if self.total_length > 0 else 0
        
        return {
            'overall_gc_percent': overall_gc,
            'mean_gc_percent': statistics.mean(gc_percentages),
            'median_gc_percent': statistics.median(gc_percentages),
            'min_gc_percent': min(gc_percentages),
            'max_gc_percent': max(gc_percentages),
            'std_dev_gc': statistics.stdev(gc_percentages) if len(gc_percentages) > 1 else 0
        }
    
    def get_n_statistics(self):
        """
        Calculate statistics for N bases (ambiguous nucleotides).
        """
        if not self.n_counts:
            return {}
        
        total_n = sum(self.n_counts)
        sequences_with_n = sum(1 for count in self.n_counts if count > 0)
        
        n_percentages = [
            (n_count / length) * 100 if length > 0 else 0
            for n_count, length in zip(self.n_counts, self.lengths)
        ]
        
        return {
            'total_n_bases': total_n,
            'sequences_with_n': sequences_with_n,
            'percent_sequences_with_n': (sequences_with_n / self.sequence_count) * 100 if self.sequence_count > 0 else 0,
            'overall_n_percent': (total_n / self.total_length) * 100 if self.total_length > 0 else 0,
            'mean_n_percent': statistics.mean(n_percentages),
            'max_n_percent': max(n_percentages) if n_percentages else 0
        }
    
    def get_length_distribution(self, bins: int = 20):
        """
        Calculate length distribution for histogram generation.
        """
        if not self.lengths:
            return {'bin_edges': [], 'counts': []}
        
        # Use numpy for histogram calculation
        counts, bin_edges = np.histogram(self.lengths, bins=bins)
        
        return {
            'bin_edges': bin_edges.tolist(),
            'counts': counts.tolist()
        }
    
    def get_summary(self):
        """
        Generate comprehensive summary of all statistics.
        """
        if self.sequence_count == 0:
            return {
                'basic': {'sequence_count': 0, 'total_length': 0},
                'length_stats': {},
                'gc_stats': {},
                'n_stats': {},
                'assembly_stats': {}
            }
        
        # Fix min_length if no sequences were added
        if self.min_length == float('inf'):
            self.min_length = 0
        
        n50, l50 = self.calculate_n50()
        n90, l90 = self.calculate_nx(90)
        
        return {
            'basic': {
                'sequence_count': self.sequence_count,
                'total_length': self.total_length,
                'min_length': self.min_length,
                'max_length': self.max_length,
                'mean_length': self.total_length / self.sequence_count
            },
            'length_stats': self.get_length_statistics(),
            'gc_stats': self.get_gc_statistics(),
            'n_stats': self.get_n_statistics(),
            'assembly_stats': {
                'n50': n50,
                'l50': l50,
                'n90': n90,
                'l90': l90
            }
        }
    
    def get_outliers(self, method: str = 'iqr', threshold: float = 1.5):
        """
        Identify outlier sequences based on length.
        """
        if not self.lengths:
            return []
        
        outliers = []
        
        if method == 'iqr':
            q1 = self._percentile(self.lengths, 25)
            q3 = self._percentile(self.lengths, 75)
            iqr = q3 - q1
            lower_bound = q1 - threshold * iqr
            upper_bound = q3 + threshold * iqr
            
            for i, length in enumerate(self.lengths):
                if length < lower_bound or length > upper_bound:
                    outliers.append((i, length))
        
        elif method == 'zscore':
            mean_length = statistics.mean(self.lengths)
            std_dev = statistics.stdev(self.lengths) if len(self.lengths) > 1 else 0
            
            if std_dev > 0:
                for i, length in enumerate(self.lengths):
                    z_score = abs(length - mean_length) / std_dev
                    if z_score > threshold:
                        outliers.append((i, length))
        
        return outliers
    
    def _percentile(self, data: List[float], percentile: float):
        """
        Calculate percentile of a dataset.
        """
        if not data:
            return 0
        
        sorted_data = sorted(data)
        n = len(sorted_data)
        
        if n == 1:
            return sorted_data[0]
        
        # Calculate index
        index = (percentile / 100) * (n - 1)
        
        # Interpolate if necessary
        if index.is_integer():
            return sorted_data[int(index)]
        else:
            lower_index = int(index)
            upper_index = lower_index + 1
            if upper_index >= n:
                return sorted_data[-1]
            
            # Linear interpolation
            weight = index - lower_index
            return sorted_data[lower_index] * (1 - weight) + sorted_data[upper_index] * weight


class AdvancedStats:
    """
    Advanced statistical analysis for specialized use cases.
    
    Provides additional metrics for genome assemblies, gene sets,
    and other specialized sequence collections.
    """
    
    def __init__(self, basic_stats: BasicStats):
        """
        Initialize with basic statistics.
        
        Args:
            basic_stats: BasicStats instance with calculated statistics
        """
        self.basic_stats = basic_stats
    
    def calculate_auN(self):
        """
        Calculate Area Under the Nx curve (auN).
        
        This metric provides a single number that summarizes the entire
        Nx curve, useful for comparing assemblies.
        """
        if not self.basic_stats.lengths:
            return 0.0
        
        # Calculate N values for x from 1 to 100
        n_values = []
        for x in range(1, 101):
            nx, _ = self.basic_stats.calculate_nx(x)
            n_values.append(nx)
        
        # Calculate area under curve using trapezoidal rule
        return sum(n_values) / 100
    
    def calculate_genome_assembly_metrics(self):
        """
        Calculate metrics specifically for genome assemblies.
        """
        if not self.basic_stats.lengths:
            return {}
        
        # Sort lengths for calculations
        sorted_lengths = sorted(self.basic_stats.lengths, reverse=True)
        
        nx_values = {}
        for x in [25, 50, 75, 90, 95]:
            nx, lx = self.basic_stats.calculate_nx(x)
            nx_values[f'n{x}'] = nx
            nx_values[f'l{x}'] = lx
        
        size_categories = {
            'large_contigs': len([l for l in sorted_lengths if l >= 10000]),
            'medium_contigs': len([l for l in sorted_lengths if 1000 <= l < 10000]),
            'small_contigs': len([l for l in sorted_lengths if 100 <= l < 1000]),
            'very_small_contigs': len([l for l in sorted_lengths if l < 100])
        }
        
        return {
            **nx_values,
            'auN': self.calculate_auN(),
            'largest_contig': max(sorted_lengths),
            'total_contigs': len(sorted_lengths),
            **size_categories
        }
    
    def calculate_gene_set_metrics(self):
        """
        Calculate metrics for gene/protein sequence collections.
        """
        if not self.basic_stats.lengths:
            return {}
        
        # Typical gene/protein length categories
        gene_categories = {
            'short_genes': len([l for l in self.basic_stats.lengths if l < 300]),
            'medium_genes': len([l for l in self.basic_stats.lengths if 300 <= l < 1500]),
            'long_genes': len([l for l in self.basic_stats.lengths if l >= 1500])
        }
        
        return {
            **gene_categories,
            'mean_gene_length': statistics.mean(self.basic_stats.lengths),
            'median_gene_length': statistics.median(self.basic_stats.lengths)
        }