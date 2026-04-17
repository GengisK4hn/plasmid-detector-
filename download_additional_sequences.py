#!/usr/bin/env python3
"""
Download additional sequences for comprehensive SV40 analysis

This script downloads:
1. SARS-CoV-2 natural sequences (negative controls)
2. SV40 positive controls
3. Additional vaccine sequences
4. Bacterial plasmid backbones

Author: Enhanced Plasmid Finder Suite
Date: April 2026
"""

import requests
import time
from pathlib import Path
import sys

class SequenceDownloader:
    """Download sequence data from NCBI and other sources."""

    def __init__(self, output_dir="data/sequences/additional"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True, parents=True)
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'Mozilla/5.0 (compatible; EnhancedPlasmidFinder/2.0)'
        })

    def download_fasta(self, accession, url=None, description=""):
        """Download FASTA sequence from NCBI or other source."""
        try:
            # Construct NCBI FASTA URL if not provided
            if url is None or "ncbi" in url.lower():
                fasta_url = f"https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&id={accession}&report=fasta&log$=seqview"
            else:
                fasta_url = url

            print(f"  Downloading {accession}...")
            if description:
                print(f"    ({description})")

            response = self.session.get(fasta_url, timeout=30)
            response.raise_for_status()

            # Save to file
            output_file = self.output_dir / f"{accession}.fasta"
            with open(output_file, 'w') as f:
                f.write(response.text)

            # Check if we got actual sequence data
            if len(response.text) < 100 or "Error" in response.text:
                print(f"    ⚠️  Warning: Possible download issue")
            else:
                print(f"    ✅ Success: {len(response.text)} bytes")

            time.sleep(1)  # Be nice to NCBI
            return True

        except Exception as e:
            print(f"    ❌ Failed: {e}")
            return False

    def download_priority_1_sars_cov2(self):
        """Download SARS-CoV-2 natural sequences (negative controls)."""
        print("\n" + "="*70)
        print("PRIORITY 1: SARS-CoV-2 Natural Sequences (Negative Controls)")
        print("="*70)

        sequences = {
            "NC_045512.2": ("https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2",
                          "Wuhan-Hu-1 reference (complete genome)"),
            "MT020880.1": ("https://www.ncbi.nlm.nih.gov/nuccore/MT020880.1",
                          "Early Wuhan spike protein"),
            "OMX095706.1": ("https://www.ncbi.nlm.nih.gov/nuccore/OMX095706.1",
                           "Delta variant spike"),
            "OMX067679.1": ("https://www.ncbi.nlm.nih.gov/nuccore/OMX067679.1",
                           "Omicron BA.1 spike"),
            "OMX067680.1": ("https://www.ncbi.nlm.nih.gov/nuccore/OMX067680.1",
                           "Omicron BA.2 spike"),
        }

        success_count = 0
        for accession, (url, desc) in sequences.items():
            if self.download_fasta(accession, url, desc):
                success_count += 1

        print(f"\n✅ Downloaded {success_count}/{len(sequences)} SARS-CoV-2 sequences")
        return success_count

    def download_priority_3_sv40_controls(self):
        """Download SV40 positive controls."""
        print("\n" + "="*70)
        print("PRIORITY 2: SV40 Positive Controls")
        print("="*70)

        sequences = {
            "NC_001669.1": ("https://www.ncbi.nlm.nih.gov/nuccore/NC_001669.1",
                           "SV40 virus complete genome"),
        }

        success_count = 0
        for accession, (url, desc) in sequences.items():
            if self.download_fasta(accession, url, desc):
                success_count += 1

        print(f"\n✅ Downloaded {success_count}/{len(sequences)} SV40 control sequences")
        return success_count

    def download_priority_4_other_vaccines(self):
        """Download other vaccine sequences for comparison."""
        print("\n" + "="*70)
        print("PRIORITY 3: Other Vaccine Platforms (Comparators)")
        print("="*70)

        # Note: Many of these may not be publicly available
        sequences = {
            # Add other vaccine sequences as they become available
            # "JX233432.1": ("", "AstraZeneca ChAdOx1"),  # Example
        }

        if not sequences:
            print("  ℹ️  No additional vaccine sequences available for download")
            print("     (Many proprietary sequences not publicly available)")
            return 0

        success_count = 0
        for accession, (url, desc) in sequences.items():
            if self.download_fasta(accession, url, desc):
                success_count += 1

        print(f"\n✅ Downloaded {success_count}/{len(sequences)} additional vaccine sequences")
        return success_count

    def download_all(self):
        """Download all priority sequences."""
        print("\n" + "="*70)
        print("ENHANCED PLASMID FINDER - SEQUENCE DOWNLOAD SUITE")
        print("="*70)
        print(f"Output directory: {self.output_dir.absolute()}")

        total = 0
        total += self.download_priority_1_sars_cov2()
        total += self.download_priority_3_sv40_controls()
        total += self.download_priority_4_other_vaccines()

        print("\n" + "="*70)
        print(f"TOTAL: Downloaded {total} sequences")
        print("="*70)

        print("\nNext steps:")
        print("1. Run batch analysis: python batch_mrna_scanner.py")
        print("2. Check results in: data/sequences/additional/")
        print("3. Compare with your Pfizer/Moderna analysis")

        return total


def main():
    """Main download function."""
    downloader = SequenceDownloader()
    downloader.download_all()


if __name__ == "__main__":
    main()
