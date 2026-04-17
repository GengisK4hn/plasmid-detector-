# Enhanced Plasmid Finder Suite - Master Controller

## Overview

The Enhanced Plasmid Finder Suite provides comprehensive plasmid and vector analysis with full integration into the main run_everything.sh pipeline. It features error correction, automatic data discovery, and multi-tool validation.

## Key Features

### 1. Error Correction & Recovery
- Automatic error detection and logging
- Attempted corrections for common failures
- Comprehensive error tracking and reporting
- Graceful degradation when tools unavailable

### 2. Data Discovery & Scanning
- Automatic discovery of sequence files across multiple directories
- Validation of file formats (FASTA, GenBank)
- Multi-source support (local, external drives, network paths)
- File integrity checking

### 3. Pipeline Integration
- Calls main run_everything.sh (20+ analyzers, 9 phases)
- ABRICATE plasmid detection
- SV40 validation (BLAST-confirmed)
- FCS (Furin Cleavage Site) detection
- Cross-tool validation

## Usage

### Basic Usage

```bash
python ENHANCED_PLASMID_FINDER_SUITE.py <input_directory> [output_directory]
```

### Example

```bash
python ENHANCED_PLASMID_FINDER_SUITE.py data/sequences suite_results
```

## Phases

### Phase 1: Data Discovery
- Scans configured directories for sequence files
- Validates file formats and integrity
- Reports discovery statistics

### Phase 2: Main Pipeline
- Calls run_everything.sh
- Runs 20+ analyzers
- 9 comprehensive analysis phases

### Phase 3: Plasmid Detection (ABRICATE)
- ABRICATE integration
- Plasmid replication system detection
- Confidence scoring

### Phase 4: SV40 Validation
- Exact 72bp enhancer search
- BLAST validation
- Position mapping

### Phase 5: Cross-Validation
- Correlates findings across tools
- Cross-referenced results
- Statistical summary

### Phase 6: Report Generation
- JSON output (machine-readable)
- Text report (human-readable)
- Error correction log

## Data Sources

The suite automatically scans:
- `data/sequences/` - Main sequence directory
- `data/sequences/additional/` - Additional references
- `data/references/` - Vaccine references
- `/media/external_drive/mckernan_supracode_data/` - External data
- `../supracode-tool/data/sequences/` - Parent project data

## Outputs

### JSON Report
```json
{
  "timestamp": "2026-04-17T...",
  "input_directory": "data/sequences",
  "phases": {
    "data_discovery": {...},
    "plasmid_detection": {...},
    "sv40_validation": {...}
  },
  "summary": {...}
}
```

### Text Report
Human-readable summary with:
- Phase results
- Error correction statistics
- Cross-validated findings

## Error Handling

The suite includes comprehensive error handling:
- Automatic retry for transient failures
- Graceful degradation when tools unavailable
- Detailed error logging
- Summary of corrected vs uncorrected errors

## Integration with Main Pipeline

The suite integrates with the main run_everything.sh by:
1. Calling the main script as Phase 2
2. Adding specialized plasmid analysis (ABRICATE)
3. Providing SV40 validation (BLAST-confirmed)
4. Cross-referencing all findings

## Requirements

- Python 3.7+
- Biopython
- ABRICATE (optional, for plasmid detection)
- BLAST+ (optional, for SV40 validation)
- Access to run_everything.sh (for main pipeline)

## Status

✅ Data discovery and validation working
✅ Error correction mechanisms in place
✅ Plasmid detection with ABRICATE
✅ SV40 validation pipeline
✅ Cross-tool validation
✅ Multi-format reporting

## Example Output

```
╔════════════════════════════════════════════════════════════════════╗
║     ENHANCED PLASMID FINDER SUITE - MASTER CONTROLLER            ║
║     (Error Correction + Data Scanning + Full Integration)        ║
╚════════════════════════════════════════════════════════════════════╝

======================================================================
  PHASE 1: DATA DISCOVERY
======================================================================

📁 Scanning: data/sequences
   Found 7 sequence files

📊 DISCOVERY SUMMARY:
   Total files found: 13
   Valid files: 10

======================================================================
  PHASE 3: PLASMID DETECTION (ABRICATE)
======================================================================

🔬 pfizer_bnt162b2_OR134577.1.fasta
   ✅ PLASMID: ColRNAI_1

🔬 moderna_mrna_1273_OR134578.1.fasta
   ✅ PLASMID: ColRNAI_1

📊 Found 7 plasmids in 10 files

✅ All phases executed
```

## Notes

- The suite continues execution even when individual tools fail
- Error correction is attempted but not guaranteed
- Results are cross-validated for reliability
- All analysis is logged for reproducibility
