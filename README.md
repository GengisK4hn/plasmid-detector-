# Enhanced PlasmidFinder

A comprehensive Python tool for detecting plasmid sequences in genomic data with enhanced reporting and detailed hit information.

## Features

- **Direct Abricate Integration**: Fast analysis without loading entire files into memory
- **Rich Hit Metadata**: Captures position, coverage percentage, and identity percentage for each marker
- **Enhanced Error Handling**: Robust timeout handling and graceful fallbacks
- **Flexible Environment Support**: Works with or without conda environments
- **Detailed Confidence Scoring**: HIGH (3+ hits), MODERATE (1-2 hits), NONE
- **Comprehensive Reporting**: Hit-level details with positions and quality metrics

## Installation

### Requirements

- Python 3.7+
- Biopython
- pandas

### Install Dependencies

```bash
pip install biopython pandas
```

### Install Abricate (Recommended)

For accurate plasmid detection, install Abricate with the PlasmidFinder database:

```bash
# Using conda (recommended)
conda install -c bioconda abricate

# Or using the environment of your choice
conda install -c bioconda -n YOUR_ENV_NAME abricate
```

## Usage

### Basic Usage

```bash
python plasmid_finder.py sequence.fasta
```

### With Conda Environment

```bash
python plasmid_finder.py sequence.fasta myenv
```

### Python API

```python
from plasmid_finder import PlasmidEnhancedDetector

# Initialize detector
detector = PlasmidEnhancedDetector()

# Analyze a single file
result = detector.analyze_fasta_file_abricate('sequence.fasta')

# Check results
if result['plasmid_detected']:
    print(f"Plasmid detected with {result['hit_count']} hits")
    for hit in result['hits']:
        print(f"  {hit['gene']}: positions {hit['start']}-{hit['end']}")
```

### Batch Analysis

```python
from plasmid_finder import analyze_fasta_files

# Analyze multiple files
fasta_files = ['seq1.fasta', 'seq2.fasta', 'seq3.fasta']
results_df = analyze_fasta_files(fasta_files)
```

## Output

The tool provides detailed information about detected plasmids:

- **Plasmid Detection**: Boolean indicating if plasmid markers were found
- **Confidence Level**: HIGH, MODERATE, or NONE based on hit count and quality
- **Markers Found**: List of plasmid replication genes detected
- **Hit Details**: For each marker
  - Gene name
  - Start and end positions
  - Coverage percentage (minimum 80%)
  - Identity percentage (minimum 80%)
- **Method Used**: Abricate, PlasmidFinder, or Built-in fallback

## Methods

The tool uses a priority-based approach:

1. **Abricate** (most accurate): Uses the PlasmidFinder database via Abricate
2. **PlasmidFinder CLI**: Direct use of PlasmidFinder tool
3. **Built-in detection** (fallback): Simplified pattern matching with limited accuracy

## Examples

See the `examples/` directory for detailed usage examples:

- `basic_usage.py`: Simple single-file analysis
- `batch_analysis.py`: Analyzing multiple files
- `custom_env.py`: Using with specific conda environments

## Requirements

```
biopython>=1.79
pandas>=1.3.0
```

## License

MIT License - see LICENSE file for details

## Contributing

Contributions are welcome! Please feel free to submit issues or pull requests.

## Citation

If you use this tool in your research, please cite:

```
Enhanced PlasmidFinder: Comprehensive plasmid detection with detailed reporting
https://github.com/YOURUSERNAME/enhanced-plasmid-finder
```

## Acknowledgments

- PlasmidFinder: https://bitbucket.org/genomicepidemiology/plasmidfinder/src/master/
- Abricate: https://github.com/tseemann/abricate
