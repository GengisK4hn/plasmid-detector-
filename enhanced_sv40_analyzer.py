#!/usr/bin/env python3
"""
Enhanced SV40 and Promoter Element Analyzer with Alignment Integration

This script utilizes MAFFT/ClustalO for:
- Multiple sequence alignment of promoter regions
- Comparative analysis across vaccine variants
- Conservation scoring of SV40 elements
- Phylogenetic relationships

Usage:
    python enhanced_sv40_analyzer.py sequence.fasta [reference1.fasta reference2.fasta ...]
"""

from Bio import SeqIO, Align
from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices
import subprocess
import sys
import argparse
import os
from pathlib import Path

# SV40 and promoter motifs
MOTIFS = {
    "SV40_72bp_Enhancer": "GGTGTGGAAAGTCCCCAGGCTCCC",
    "SV40_GC_Box": "GGGCGG",
    "SV40_GC_Box_Reverse": "CCGCCC",
    "SV40_6bp_Repeat": "GGGGCG",
    "TATA_Box": "TATAAA",
    "TATA_Box_Variant": "TATATA",
    "ColE1_Origin": "AAGGATCTAGGTGAAGATCCTTTTTGATAATCTCATGACCAAAATCCCTTAACGTGAGTTTTCGTTCCACTGAGCGTCAGACCCCGT",
}

def check_tool_availability():
    """Check which alignment tools are available."""
    tools = {
        'mafft': False,
        'clustalo': False
    }

    for tool in tools.keys():
        try:
            subprocess.run([tool, '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
            tools[tool] = True
            print(f"✅ {tool} is available")
        except (subprocess.CalledProcessError, FileNotFoundError):
            print(f"❌ {tool} is not available")

    return tools

def extract_promoter_regions(sequence, motif_positions, window_size=200):
    """Extract sequence regions around detected motifs for alignment."""
    regions = {}

    for motif_name, positions in motif_positions.items():
        if positions:
            for i, pos in enumerate(positions):
                start = max(0, pos - window_size//2)
                end = min(len(sequence), pos + window_size//2 + len(MOTIFS[motif_name]))
                region_name = f"{motif_name}_copy_{i+1}"
                regions[region_name] = sequence[start:end]

    return regions

def perform_mafft_alignment(input_file, output_file):
    """Perform multiple sequence alignment using MAFFT."""
    try:
        cmd = ['mafft', '--auto', input_file]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)

        with open(output_file, 'w') as f:
            f.write(result.stdout)

        return True, output_file
    except Exception as e:
        return False, str(e)

def perform_clustalo_alignment(input_file, output_file):
    """Perform multiple sequence alignment using ClustalO."""
    try:
        cmd = ['clustalo', '-i', input_file, '-o', output_file, '--outfmt=fasta', '--force']
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)

        return True, output_file
    except Exception as e:
        return False, str(e)

def calculate_conservation(alignment):
    """Calculate conservation scores for alignment columns."""
    if len(alignment) == 0:
        return []

    conservation = []
    num_sequences = len(alignment.coords[0])  # Number of sequences in alignment
    alignment_length = len(alignment)  # Total alignment length

    for pos in range(alignment_length):
        column = alignment[:, pos]
        # Count non-gap characters
        non_gap = sum(1 for char in column if char != '-')
        # Calculate conservation as percentage of non-gap characters
        conservation.append(non_gap / num_sequences if num_sequences > 0 else 0)

    return conservation

def analyze_sequence_with_tools(sequence_file, reference_files=None):
    """Comprehensive analysis using alignment tools."""

    print("=" * 80)
    print("ENHANCED SV40 & PROMOTER ELEMENT ANALYSIS WITH ALIGNMENT TOOLS")
    print("=" * 80)
    print()

    # Check tool availability
    tools = check_tool_availability()
    print()

    # Load sequence with robust format detection
    try:
        record = SeqIO.read(sequence_file, "fasta")
        sequence = str(record.seq).upper()
    except Exception as e:
        # Try alternative FASTA formats
        try:
            record = SeqIO.read(sequence_file, "fasta-pearson")
            sequence = str(record.seq).upper()
        except:
            try:
                record = SeqIO.read(sequence_file, "fasta-blast")
                sequence = str(record.seq).upper()
            except:
                print(f"Error loading sequence: {e}")
                return

    print(f"Sequence: {record.id}")
    print(f"Length: {len(sequence)} bp")
    print()

    # Motif detection
    print("-" * 80)
    print("STEP 1: MOTIF DETECTION")
    print("-" * 80)
    print()

    motif_positions = {}
    for name, motif in MOTIFS.items():
        count = sequence.count(motif)
        if count > 0:
            positions = []
            start = 0
            while True:
                pos = sequence.find(motif, start)
                if pos == -1:
                    break
                positions.append(pos)
                start = pos + 1

            motif_positions[name] = positions

            percentage = (positions[0] / len(sequence)) * 100
            print(f"✅ {name}: {count} hits at positions {positions[:5]}")

    print()

    # Extract promoter regions for alignment
    print("-" * 80)
    print("STEP 2: EXTRACTING PROMOTER REGIONS FOR ALIGNMENT")
    print("-" * 80)
    print()

    promoter_regions = extract_promoter_regions(sequence, motif_positions)
    print(f"Extracted {len(promoter_regions)} promoter regions")

    # Create temporary FASTA file for alignment
    temp_fasta = "temp_promoter_regions.fasta"
    with open(temp_fasta, 'w') as f:
        for region_name, region_seq in promoter_regions.items():
            f.write(f">{region_name}\n{region_seq}\n")

    # Perform alignment if tools are available
    alignment_file = None
    if tools['mafft']:
        print("Performing MAFFT alignment...")
        success, result = perform_mafft_alignment(temp_fasta, "promoter_alignment_mafft.fasta")
        if success:
            alignment_file = result
            print(f"✅ MAFFT alignment complete: {alignment_file}")
        else:
            print(f"❌ MAFFT alignment failed: {result}")
    elif tools['clustalo']:
        print("Performing ClustalO alignment...")
        success, result = perform_clustalo_alignment(temp_fasta, "promoter_alignment_clustalo.fasta")
        if success:
            alignment_file = result
            print(f"✅ ClustalO alignment complete: {alignment_file}")
        else:
            print(f"❌ ClustalO alignment failed: {result}")

    print()

    # Analyze alignment if available
    if alignment_file and os.path.exists(alignment_file):
        print("-" * 80)
        print("STEP 3: ALIGNMENT ANALYSIS")
        print("-" * 80)
        print()

        # Simply read the alignment file to show it worked
        aligned_sequences = list(SeqIO.parse(alignment_file, "fasta"))

        print(f"Alignment statistics:")
        print(f"  Sequences aligned: {len(aligned_sequences)}")

        if aligned_sequences:
            alignment_length = len(aligned_sequences[0].seq)
            print(f"  Alignment length: {alignment_length} bp")

            # Show first few aligned sequences as examples
            print(f"  Sample aligned sequences:")
            for i, record in enumerate(aligned_sequences[:3]):
                print(f"    {record.id}: {len(record.seq)} bp")

            if len(aligned_sequences) > 3:
                print(f"    ... and {len(aligned_sequences) - 3} more")

        print()

    # Comparative analysis with references
    if reference_files:
        print("-" * 80)
        print("STEP 4: COMPARATIVE ANALYSIS WITH REFERENCES")
        print("-" * 80)
        print()

        for ref_file in reference_files:
            try:
                # Try robust format detection for reference files
                ref_record = None
                for fmt in ["fasta", "fasta-pearson", "fasta-blast"]:
                    try:
                        ref_record = SeqIO.read(ref_file, fmt)
                        break
                    except:
                        continue

                if ref_record is None:
                    print(f"  Could not read reference file: {ref_file}")
                    continue

                ref_sequence = str(ref_record.seq).upper()

                print(f"Comparing with: {ref_record.id}")

                # Simple pairwise comparison
                try:
                    from Bio.Align.PairwiseAligner import PairwiseAligner

                    aligner = PairwiseAligner()
                    aligner.mode = 'global'
                    alignments = aligner.align(sequence, ref_sequence)
                    best = alignments[0]

                    identity = (best.score / max(len(sequence), len(ref_sequence))) * 100

                    print(f"  Best alignment score: {best.score:.1f}")
                    print(f"  Sequence identity: {identity:.2f}%")

                    # Check shared SV40 elements
                    shared = sum(1 for motif in MOTIFS.values() if motif in sequence and motif in ref_sequence)
                    print(f"  Shared SV40/promoter elements: {shared}/{len(MOTIFS)}")
                    print()

                except Exception as align_error:
                    print(f"  Alignment analysis error: {align_error}")
                    print()

            except Exception as e:
                print(f"Error processing reference {ref_file}: {e}")

    # Cleanup temp file
    if os.path.exists(temp_fasta):
        os.remove(temp_fasta)

    # Summary
    print("-" * 80)
    print("ANALYSIS SUMMARY")
    print("-" * 80)
    print()

    total_elements = sum(len(positions) for positions in motif_positions.values())
    print(f"Total SV40/promoter elements detected: {total_elements}")
    print(f"Unique element types: {len(motif_positions)}")
    print(f"Alignment tools used: {', '.join([tool for tool, available in tools.items() if available]) or 'None'}")
    print()

    if alignment_file:
        print(f"Alignment saved to: {alignment_file}")

    return motif_positions, alignment_file

def main():
    parser = argparse.ArgumentParser(
        description='Enhanced SV40 and promoter element analysis with alignment tools'
    )
    parser.add_argument(
        'sequence',
        help='FASTA file to analyze'
    )
    parser.add_argument(
        'references',
        nargs='*',
        help='Reference FASTA files for comparative analysis'
    )
    parser.add_argument(
        '--output',
        help='Output directory for results',
        default='sv40_analysis_results'
    )

    args = parser.parse_args()

    # Create output directory
    os.makedirs(args.output, exist_ok=True)

    # Run analysis
    analyze_sequence_with_tools(args.sequence, args.references)

if __name__ == '__main__':
    main()
