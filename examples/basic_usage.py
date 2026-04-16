#!/usr/bin/env python3
"""
Basic usage example for Enhanced PlasmidFinder

This example demonstrates how to analyze a single FASTA file
for plasmid content.
"""

import sys
import os

# Add parent directory to path to import the module
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from plasmid_finder import PlasmidEnhancedDetector


def main():
    """Demonstrate basic usage of the plasmid detector"""

    # Initialize the detector
    detector = PlasmidEnhancedDetector()

    # Analyze a single file
    fasta_file = "example.fasta"

    print(f"Analyzing {fasta_file}...")
    print("=" * 70)

    result = detector.analyze_fasta_file_abricate(fasta_file)

    if result and result.get('analysis_complete'):
        print("\nResults:")
        print(f"  File: {result['file']}")
        print(f"  Plasmid Detected: {result['plasmid_detected']}")
        print(f"  Confidence: {result.get('confidence', 'N/A')}")
        print(f"  Method: {result.get('method', 'N/A')}")

        if result['plasmid_detected']:
            print(f"  Hit Count: {result.get('hit_count', 0)}")
            print("\nDetailed Hits:")
            for hit in result.get('hits', []):
                print(f"    - {hit['gene']}")
                print(f"      Position: {hit['start']}-{hit['end']}")
                print(f"      Coverage: {hit['coverage']:.1f}%")
                print(f"      Identity: {hit['identity']:.1f}%")
    else:
        print("Analysis failed or returned no results")


if __name__ == "__main__":
    main()
