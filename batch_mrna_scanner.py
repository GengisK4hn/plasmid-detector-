#!/usr/bin/env python3
"""
Batch mRNA Scanner with MAFFT Alignment

Scans all mRNA sequences in a directory and performs:
- SV40/promoter element detection
- MAFFT multiple sequence alignment
- Comparative analysis across sequences

Usage:
    python batch_mrna_scanner.py data/sequences/
"""

from Bio import SeqIO, Align
try:
    from Bio.Align import PairwiseAligner
except ImportError:
    PairwiseAligner = None
import subprocess
import sys
import os
import argparse
from pathlib import Path
import glob

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
    tools = {'mafft': False, 'clustalo': False}

    for tool in tools.keys():
        try:
            subprocess.run([tool, '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
            tools[tool] = True
        except (subprocess.CalledProcessError, FileNotFoundError):
            pass

    return tools

def load_sequence_robust(filepath):
    """Load sequence with multiple format attempts."""
    for fmt in ["fasta", "fasta-pearson", "fasta-blast", "genbank"]:
        try:
            record = SeqIO.read(filepath, fmt)
            return str(record.seq).upper(), record.id, len(record.seq)
        except:
            continue
    return None, None, None

def detect_motifs(sequence, seq_id):
    """Detect all SV40/promoter motifs in a sequence."""
    results = {
        'sequence_id': seq_id,
        'motifs_found': {},
        'total_elements': 0
    }

    for motif_name, motif_seq in MOTIFS.items():
        count = sequence.count(motif_seq)
        if count > 0:
            positions = []
            start = 0
            while True:
                pos = sequence.find(motif_seq, start)
                if pos == -1:
                    break
                positions.append(pos)
                start = pos + 1

            results['motifs_found'][motif_name] = {
                'count': count,
                'positions': positions
            }
            results['total_elements'] += count

    return results

def perform_mafft_alignment(input_file, output_file):
    """Perform MAFFT multiple sequence alignment."""
    try:
        cmd = ['mafft', '--auto', input_file]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)

        with open(output_file, 'w') as f:
            f.write(result.stdout)

        return True, output_file
    except Exception as e:
        return False, str(e)

def scan_directory(sequence_dir):
    """Scan all FASTA files in directory."""

    print("=" * 80)
    print("BATCH mRNA SCANNER WITH MAFFT ALIGNMENT")
    print("=" * 80)
    print()

    # Check tools
    tools = check_tool_availability()
    print(f"Available alignment tools:")
    for tool, available in tools.items():
        status = "✅" if available else "❌"
        print(f"  {status} {tool}")
    print()

    # Find all FASTA files
    fasta_files = []
    for ext in ['*.fasta', '*.fa', '*.fna']:
        fasta_files.extend(glob.glob(os.path.join(sequence_dir, ext)))

    if not fasta_files:
        print(f"No FASTA files found in {sequence_dir}")
        return

    print(f"Found {len(fasta_files)} FASTA files to scan")
    print()

    # Scan each file
    all_results = []
    sequences_with_sv40 = []

    print("-" * 80)
    print("STEP 1: SCANNING FOR SV40/PROMOTER ELEMENTS")
    print("-" * 80)
    print()

    for i, fasta_file in enumerate(fasta_files, 1):
        print(f"[{i}/{len(fasta_files)}] {os.path.basename(fasta_file)}")

        sequence, seq_id, length = load_sequence_robust(fasta_file)
        if sequence is None:
            print(f"  ❌ Could not read file")
            continue

        print(f"  Length: {length} bp")

        # Detect motifs
        results = detect_motifs(sequence, seq_id)

        if results['total_elements'] > 0:
            print(f"  ✅ SV40/promoter elements found: {results['total_elements']} total")
            for motif_name, motif_data in results['motifs_found'].items():
                print(f"     - {motif_name}: {motif_data['count']} copies")
            sequences_with_sv40.append(results)
        else:
            print(f"  ❌ No SV40/promoter elements detected")

        all_results.append(results)
        print()

    # Summary
    print("-" * 80)
    print("SCAN SUMMARY")
    print("-" * 80)
    print()

    print(f"Total sequences scanned: {len(all_results)}")
    print(f"Sequences with SV40/promoter elements: {len(sequences_with_sv40)}")

    if sequences_with_sv40:
        print()
        print("Top sequences by element count:")

        # Sort by total elements
        sorted_results = sorted(sequences_with_sv40, key=lambda x: x['total_elements'], reverse=True)

        for i, result in enumerate(sorted_results[:10], 1):
            print(f"  {i}. {result['sequence_id']}: {result['total_elements']} elements")

        # Create MAFFT alignment if we have sequences with SV40
        if len(sequences_with_sv40) >= 2 and tools['mafft']:
            print()
            print("-" * 80)
            print("STEP 2: MAFFT MULTIPLE SEQUENCE ALIGNMENT")
            print("-" * 80)
            print()

            # Create combined FASTA for alignment
            combined_fasta = "temp_combined_sequences.fasta"
            with open(combined_fasta, 'w') as f:
                for result in sorted_results[:10]:  # Align top 10 sequences
                    # Reload the full sequence
                    for fasta_file in fasta_files:
                        seq, seq_id, _ = load_sequence_robust(fasta_file)
                        if seq_id == result['sequence_id']:
                            f.write(f">{seq_id}\n{seq}\n")
                            break

            # Perform MAFFT alignment
            print("Performing MAFFT alignment...")
            success, result = perform_mafft_alignment(combined_fasta, "mRNA_alignment_mafft.fasta")

            if success:
                print(f"✅ MAFFT alignment complete: mRNA_alignment_mafft.fasta")

                # Read alignment stats
                aligned_sequences = list(SeqIO.parse(result, "fasta"))
                print(f"   Sequences aligned: {len(aligned_sequences)}")
                if aligned_sequences:
                    print(f"   Alignment length: {len(aligned_sequences[0].seq)} bp")
            else:
                print(f"❌ MAFFT alignment failed: {result}")

            # Cleanup
            if os.path.exists(combined_fasta):
                os.remove(combined_fasta)

        # Comparative analysis
        print()
        print("-" * 80)
        print("STEP 3: COMPARATIVE ANALYSIS")
        print("-" * 80)
        print()

        if len(sequences_with_sv40) >= 2:
            # Compare pairwise similarity
            print("Pairwise sequence similarity (top sequences):")

            for i in range(min(3, len(sorted_results))):
                for j in range(i + 1, min(3, len(sorted_results))):
                    seq1_id = sorted_results[i]['sequence_id']
                    seq2_id = sorted_results[j]['sequence_id']

                    # Find and load sequences
                    seq1, seq2 = None, None
                    for fasta_file in fasta_files:
                        seq, seq_id, _ = load_sequence_robust(fasta_file)
                        if seq_id == seq1_id:
                            seq1 = seq
                        elif seq_id == seq2_id:
                            seq2 = seq
                        if seq1 and seq2:
                            break

                    if seq1 and seq2:
                        if PairwiseAligner:
                            aligner = PairwiseAligner()
                            aligner.mode = 'global'
                            alignments = aligner.align(seq1, seq2)
                            best = alignments[0]

                            identity = (best.score / max(len(seq1), len(seq2))) * 100
                        else:
                            # Simple identity calculation
                            min_len = min(len(seq1), len(seq2))
                            matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
                            identity = (matches / max(len(seq1), len(seq2))) * 100

                        shared = sum(1 for motif in MOTIFS.values()
                                  if motif in seq1 and motif in seq2)

                        print(f"  {seq1_id} vs {seq2_id}:")
                        print(f"     Identity: {identity:.2f}%")
                        print(f"     Shared motifs: {shared}/{len(MOTIFS)}")
                        print()

    print()
    print("=" * 80)
    print("SCAN COMPLETE")
    print("=" * 80)

def main():
    parser = argparse.ArgumentParser(
        description='Batch mRNA scanner with MAFFT alignment'
    )
    parser.add_argument(
        'directory',
        help='Directory containing FASTA files to scan'
    )

    args = parser.parse_args()

    if not os.path.isdir(args.directory):
        print(f"Error: {args.directory} is not a valid directory")
        sys.exit(1)

    scan_directory(args.directory)

if __name__ == '__main__':
    main()
