# FastaCheck

A comprehensive Python tool for validating and analyzing FASTA files with support for compressed formats and streaming processing.

## Features

- **Format Validation**: Validates FASTA file format and identifies common errors
- **Statistical Analysis**: Comprehensive sequence statistics including N50, GC content, and length distributions
- **Compressed File Support**: Handles `.gz` and `.bz2` compressed files automatically
- **Streaming Mode**: Memory-efficient processing for large files
- **Assembly Metrics**: Specialized metrics for genome assemblies (N50, auN, contig size distributions)
- **Gene Set Analysis**: Metrics tailored for gene/protein sequence collections
- **Outlier Detection**: Identifies unusual sequences using IQR or Z-score methods
- **Multiple Output Formats**: Results in JSON, CSV, or plain text formats

## Installation

```bash
pip install fastacheck
```

## Quick Start

### Command Line Usage

**Validate a FASTA file:**
```bash
fastacheck validate input.fasta
```

**Analyze sequences with basic statistics:**
```bash
fastacheck analyze input.fasta
```

**Detailed analysis with assembly metrics:**
```bash
fastacheck analyze genome.fasta --detailed --assembly --output results.json
```

**Process compressed files in streaming mode:**
```bash
fastacheck analyze large_file.fasta.gz --streaming
```

### Python API Usage

```python
from fastacheck.parser import FastaParser
from fastacheck.stats import BasicStats, AdvancedStats

# Parse and validate a FASTA file
parser = FastaParser("sequences.fasta")

# Validate format
validation = parser.validate_format()
if validation['is_valid']:
    print("Valid FASTA file")

# Calculate statistics
stats = BasicStats()
for header, sequence in parser.parse_sequences():
    stats.add_sequence(sequence, header['sequence_id'])

# Get comprehensive summary
summary = stats.get_summary()
print(f"Total sequences: {summary['basic']['sequence_count']}")
print(f"N50: {summary['assembly_stats']['n50']} bp")
```

## Command Line Options

### Validation
```bash
fastacheck validate [OPTIONS] FILE

Options:
  --streaming     Use streaming mode for large files
  --quiet         Only show validation result (pass/fail)
```

### Analysis
```bash
fastacheck analyze [OPTIONS] FILE

Options:
  --streaming                Use streaming mode for large files
  --detailed                 Show detailed length and GC statistics
  --assembly                 Show genome assembly metrics (N50, auN, etc.)
  --gene-set                 Show gene/protein set metrics
  --outliers {iqr,zscore}    Detect outlier sequences
  --threshold FLOAT          Threshold for outlier detection (default: 1.5)
  --output FILE              Output results to file
  --format {json,txt,csv}    Output format (default: txt)
```

## Statistics Provided

### Basic Statistics
- Sequence count and total length
- Length distribution (min, max, mean, median, standard deviation)
- GC content analysis
- N base (ambiguous nucleotide) statistics

### Assembly Metrics
- **N50, N90, N95**: Standard assembly quality metrics
- **L50, L90, L95**: Number of sequences in Nx sets
- **auN**: Area under the Nx curve
- **Contig size distribution**: Categorized by length ranges

### Gene Set Metrics
- Length distribution categories for genes/proteins
- Mean and median gene lengths

## Examples

**Basic validation:**
```bash
$ fastacheck validate sequences.fasta
Validating: sequences.fasta
✓ Valid FASTA file (1,234 sequences)
```

**Detailed analysis output:**
```bash
$ fastacheck analyze genome.fasta --detailed --assembly
Analyzing: genome.fasta
Processing sequences...

✓ Analysis complete
Sequences: 2,847
Total length: 125,834,567 bp
Average length: 44,201.2 bp
Shortest: 201 bp
Longest: 2,847,392 bp
N50: 156,847 bp
L50: 287 sequences
Overall GC content: 41.2%

=== Genome Assembly Metrics ===
N25: 284,751 bp (L25: 143)
N50: 156,847 bp (L50: 287)
N75: 89,234 bp (L75: 523)
N90: 45,123 bp (L90: 891)
auN (Area under Nx curve): 145,234.7
```

## Modules

### Parser Module (`fastacheck.parser`)
- `FastaParser`: Main parsing class with validation and sequence extraction
- Handles compressed files automatically
- Streaming mode support for memory efficiency
- Comprehensive error reporting

### Statistics Module (`fastacheck.stats`)
- `BasicStats`: Core statistical calculations
- `AdvancedStats`: Specialized metrics for assemblies and gene sets
- N50/Nx calculations, GC analysis, outlier detection
- Length distribution analysis

## Requirements

- Python ≥ 3.7
- NumPy (for statistical calculations)

## License

MIT License

## Contributing

Contributions are welcome! Please feel free to submit issues, feature requests, or pull requests.