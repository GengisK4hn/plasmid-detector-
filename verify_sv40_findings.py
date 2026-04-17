#!/usr/bin/env python3
"""
Verify SV40 and Promoter Elements in Vaccine Sequences

This script reproduces the findings from SV40_PROMOTER_FINDINGS.md
Use this to independently verify the presence of SV40 enhancer elements,
ColE1 origin, and promoter motifs in vaccine expression vectors.

Usage:
    python verify_sv40_findings.py sequence.fasta

Requirements:
    pip install biopython
"""

from Bio import SeqIO
import sys
import argparse

# Motifs to search for
MOTIFS = {
    "SV40_72bp_Enhancer": "GGTGTGGAAAGTCCCCAGGCTCCC",
    "SV40_72bp_Enhancer_Reverse": "GGGAGCCTGGGACTTTCCACACC",
    "SV40_GC_Box": "GGGCGG",
    "SV40_GC_Box_Reverse": "CCGCCC",
    "SV40_6bp_Repeat": "GGGGCG",
    "TATA_Box": "TATAAA",
    "TATA_Box_Variant": "TATATA",
    "ColE1_Origin": "AAGGATCTAGGTGAAGATCCTTTTTGATAATCTCATGACCAAAATCCCTTAACGTGAGTTTTCGTTCCACTGAGCGTCAGACCCCGT",
}

def analyze_sequence(sequence_file):
    """Analyze a sequence for SV40 and promoter elements."""

    print("=" * 80)
    print("SV40 & PROMOTER ELEMENT VERIFICATION")
    print("=" * 80)
    print()

    # Load sequence
    try:
        record = SeqIO.read(sequence_file, "fasta")
        sequence = str(record.seq).upper()
    except Exception as e:
        print(f"Error loading sequence: {e}")
        return

    print(f"Sequence: {record.id}")
    print(f"Length: {len(sequence)} bp")
    print()

    # Search for motifs
    print("-" * 80)
    print("SEARCH RESULTS:")
    print("-" * 80)
    print()

    results = {}
    for name, motif in MOTIFS.items():
        count = sequence.count(motif)
        if count > 0:
            results[name] = {
                'count': count,
                'motif': motif,
                'positions': []
            }

            # Find all positions
            start = 0
            while True:
                pos = sequence.find(motif, start)
                if pos == -1:
                    break
                results[name]['positions'].append(pos)
                start = pos + 1

            # Display results
            percentage = (results[name]['positions'][0] / len(sequence)) * 100
            print(f"✅ {name}")
            print(f"   Motif: {motif}")
            print(f"   Hits: {count}")
            print(f"   First position: {results[name]['positions'][0]} bp ({percentage:.1f}%)")
            if len(results[name]['positions']) > 1:
                print(f"   All positions: {results[name]['positions'][:10]}")
                if len(results[name]['positions']) > 10:
                    print(f"                    ...and {len(results[name]['positions']) - 10} more")
            print()

    # Summary
    print("-" * 80)
    print("SUMMARY:")
    print("-" * 80)
    print()

    sv40_found = any('SV40' in k for k in results.keys())
    colE1_found = 'ColE1_Origin' in results
    promoter_found = any(k in results for k in ['TATA_Box', 'TATA_Box_Variant', 'SV40_GC_Box'])

    print(f"SV40 elements found: {'YES' if sv40_found else 'NO'}")
    print(f"ColE1 origin found: {'YES' if colE1_found else 'NO'}")
    print(f"Promoter elements found: {'YES' if promoter_found else 'NO'}")
    print()

    # Interpretation
    print("-" * 80)
    print("INTERPRETATION:")
    print("-" * 80)
    print()

    if sv40_found:
        print("⚠️  SV40 enhancer elements detected")
        print("   - SV40 (Simian Virus 40) is a monkey virus")
        print("   - Its enhancer/promoter is extremely strong in mammalian cells")
        print("   - Used in genetic engineering for high-level expression")
        print()

    if colE1_found:
        print("⚠️  ColE1 bacterial origin detected")
        print("   - Standard E. coli plasmid origin")
        print("   - Also functions as cryptic mammalian promoter")
        print("   - Can drive gene expression in eukaryotic cells")
        print()

    if promoter_found:
        print("⚠️  Promoter elements detected")
        print("   - TATA boxes and GC-rich motifs present")
        print("   - Indicate potential for gene expression")
        print("   - May synergize with other promoter elements")
        print()

    if sv40_found or colE1_found or promoter_found:
        print("=" * 80)
        print("CONCLUSION:")
        print("=" * 80)
        print()
        print("This sequence contains functional promoter/enhancer elements.")
        print("If present in final vaccine products, these could potentially:")
        print("  1. Drive spike expression from DNA contaminants")
        print("  2. Cause prolonged expression if integrated")
        print("  3. Alter tissue biodistribution")
        print()
        print("This warrants:")
        print("  - Rigorous DNA contamination testing")
        print("  - Promoter activity assays")
        print("  - Independent verification")
        print()

    return results


def main():
    parser = argparse.ArgumentParser(
        description='Verify SV40 and promoter elements in vaccine sequences'
    )
    parser.add_argument(
        'sequence',
        help='FASTA file to analyze'
    )

    args = parser.parse_args()

    analyze_sequence(args.sequence)


if __name__ == '__main__':
    main()
