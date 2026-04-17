# Enhanced PlasmidFinder

A comprehensive Python tool for detecting plasmid sequences in genomic data with enhanced reporting and detailed hit information.

## Features

- **Direct Abricate Integration**: Fast analysis without loading entire files into memory
- **Rich Hit Metadata**: Captures position, coverage percentage, and identity percentage for each marker
- **Enhanced Error Handling**: Robust timeout handling and graceful fallbacks
- **Flexible Environment Support**: Works with or without conda environments
- **Detailed Confidence Scoring**: HIGH (3+ hits), MODERATE (1-2 hits), NONE
- **Comprehensive Reporting**: Hit-level details with positions and quality metrics
- **Vaccine Plasmid Signature Detection**: Identifies vaccine plasmid signatures (Pfizer, Moderna, etc.)
- **SV40 & Promoter Element Detection**: Detects viral enhancers and promoter motifs

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

### Basic Plasmid Detection

```bash
# Detect bacterial plasmid origins
python plasmid_finder.py sequence.fasta
```

### SV40 & Promoter Element Detection

```bash
# Detect SV40 enhancer elements and promoter motifs
python verify_sv40_findings.py sequence.fasta
```

### Python API

#### Plasmid Detection

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

#### SV40/Promoter Detection

```python
from Bio import SeqIO

# Load sequence
record = SeqIO.read("sequence.fasta", "fasta")
sequence = str(record.seq).upper()

# Search for SV40 72bp enhancer
sv40_enhancer = "GGTGTGGAAAGTCCCCAGGCTCCC"
count = sequence.count(sv40_enhancer)
print(f"SV40 72bp enhancer: {count} hits")

# Search for ColE1 origin
colE1 = "AAGGATCTAGGTGAAGATCCTTTTTGATAATCTCATGACCAAAATCCCTTAACGTGAGTTTTCGTTCCACTGAGCGTCAGACCCCGT"
count = sequence.count(colE1)
print(f"ColE1 origin: {count} hits")
```

## What Gets Detected

### By PlasmidFinder (via ABRICATE):
- Bacterial plasmid replication origins
- Plasmid types (rep families)
- Host range predictions
- Example: ColE1, rep5, rep7, rep16, etc.

### By SV40/Promoter Detector:
- **SV40 72bp Enhancer** (Simian Virus 40 transcriptional enhancer)
- **SV40 GC Boxes** (Sp1 transcription factor binding sites)
- **SV40 6bp Repeats** (transcriptional regulation)
- **TATA Boxes** (core promoter elements)
- **ColE1 Origin** (bacterial origin + cryptic mammalian promoter)
- **Sp1 Sites** (transcription factor binding)

### Key Difference:

**ABRICATE/PlasmidFinder** searches for:
- Bacterial plasmid replication origins
- Purpose: Identify plasmid types and bacterial hosts
- Database: Curated bacterial replication systems

**SV40/Promoter Detector** searches for:
- Viral enhancers and promoter elements
- Transcriptional regulatory motifs
- Purpose: Identify mammalian gene expression systems
- Database: Scientific literature and motif databases

## Output

### Plasmid Detection Output:

- **Plasmid Detection**: Boolean indicating if plasmid markers were found
- **Confidence Level**: HIGH, MODERATE, or NONE based on hit count
- **Markers Found**: List of plasmid replication genes detected
- **Hit Details**: Position, coverage, identity for each marker
- **Method Used**: Abricate, PlasmidFinder, or Built-in fallback

### SV40/Promoter Detection Output:

- **Element Name**: Type of element detected (e.g., SV40_72bp_Enhancer)
- **Motif Sequence**: Exact sequence matched
- **Position**: Location in sequence (bp)
- **Percentage**: Position as percentage of sequence length
- **Hit Count**: Number of copies found

## Findings Documentation

### SV40 & Promoter Detection Results:

Detailed findings from analysis of Pfizer BNT162b2 expression vector (GenBank OR134577.1) are documented in:

**`SV40_PROMOTER_FINDINGS.md`**

Contains:
- Complete detection results with positions
- Motif sequences and counts
- Methodology documentation
- Reproducibility instructions
- Data-only presentation (no speculation)

## Methods

### Plasmid Detection Priority:

1. **Abricate** (most accurate): Uses PlasmidFinder database
2. **PlasmidFinder CLI**: Direct use of PlasmidFinder tool
3. **Built-in detection** (fallback): Simplified pattern matching

### SV40/Promoter Detection:

1. Load sequence using BioPython
2. Convert to uppercase
3. Search for exact motif matches
4. Record positions and counts
5. Map locations relative to sequence length

### Motif Database Sources:

- **SV40 enhancer**: SV40 virus genome (NCBI Reference Sequence)
- **Promoter motifs**: TRANSFAC database
- **Plasmid origins**: PlasmidFinder database

## Examples

See the `examples/` directory for detailed usage examples:

- `basic_usage.py`: Simple single-file analysis
- `batch_analysis.py`: Analyzing multiple files
- `custom_env.py`: Using with specific conda environments

## Verification

To verify the SV40 and promoter findings:

```bash
python verify_sv40_findings.py pfizer_bnt162b2.fasta
```

Expected results for Pfizer BNT162b2 (OR134577.1):
- SV40_72bp_Enhancer: 2 hits
- SV40_GC_Box: 2 hits
- SV40_GC_Box_Reverse: 6 hits
- SV40_6bp_Repeat: 3 hits
- TATA_Box: 1 hit
- TATA_Box_Variant: 1 hit
- ColE1_Origin: 1 hit

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
https://github.com/GengisK4hn/plasmid-detector-
```

## Acknowledgments

- PlasmidFinder: https://bitbucket.org/genomicepidemiology/plasmidfinder/src/master/
- Abricate: https://github.com/tseemann/abricate
- SV40 Reference: NCBI GenBank (SV40 virus genome)
- TRANSFAC: Transcription Factor database
