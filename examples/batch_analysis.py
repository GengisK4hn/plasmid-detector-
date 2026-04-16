#!/usr/bin/env python3
"""
Batch analysis example for Enhanced PlasmidFinder

This example demonstrates how to analyze multiple FASTA files
and generate a comprehensive report.
"""

import sys
import os
from pathlib import Path

# Add parent directory to path to import the module
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from plasmid_finder import analyze_fasta_files


def main():
    """Demonstrate batch analysis of multiple FASTA files"""

    # Get all FASTA files from a directory
    data_dir = Path("data/sequences")

    if not data_dir.exists():
        print(f"Directory {data_dir} not found. Creating example usage...")
        fasta_files = [
            "example1.fasta",
            "example2.fasta",
            "example3.fasta"
        ]
        print("Replace with your actual file paths:")
        for f in fasta_files:
            print(f"  {f}")
        return

    # Find all FASTA files
    fasta_files = list(data_dir.glob("*.fasta")) + list(data_dir.glob("*.fa"))

    if not fasta_files:
        print(f"No FASTA files found in {data_dir}")
        return

    print(f"Found {len(fasta_files)} FASTA files")
    print("=" * 70)

    # Analyze all files
    results_df = analyze_fasta_files([str(f) for f in fasta_files])

    # Save results to CSV
    output_file = "plasmid_analysis_results.csv"
    results_df.to_csv(output_file, index=False)
    print(f"\nResults saved to {output_file}")

    # Print summary statistics
    print("\nSummary Statistics:")
    print(f"  Total files analyzed: {len(results_df)}")
    print(f"  Plasmids detected: {sum(results_df['Plasmid_Detected'] == 'YES')}")
    print(f"  High confidence: {sum(results_df['Confidence'] == 'HIGH')}")
    print(f"  Moderate confidence: {sum(results_df['Confidence'] == 'MODERATE')}")


if __name__ == "__main__":
    main()
