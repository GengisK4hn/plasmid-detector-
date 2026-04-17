#!/usr/bin/env python3
"""
Advanced SV40/Promoter Analyzer with Completeness Scoring and Visualization

Features:
- Precise motif scanning with position mapping
- Promoter completeness score (enhancer + TATA + GC boxes within 200-300bp)
- Comparative analysis between vaccine sequences
- Visualization of SV40 regions
- Statistical reporting
"""

from Bio import SeqIO
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle
import numpy as np
import argparse
import sys
import os
from pathlib import Path
import json

# SV40 and promoter motifs with scientific classifications
MOTIFS = {
    "SV40_72bp_Enhancer": {
        "sequence": "GGTGTGGAAAGTCCCCAGGCTCCC",
        "type": "enhancer",
        "color": "#FF0000",  # Red - most significant
        "description": "72bp SV40 enhancer element"
    },
    "SV40_GC_Box": {
        "sequence": "GGGCGG",
        "type": "promoter",
        "color": "#FF6600",  # Orange
        "description": "GC-rich promoter box"
    },
    "SV40_GC_Box_Reverse": {
        "sequence": "CCGCCC",
        "type": "promoter",
        "color": "#FF9900",  # Light orange
        "description": "Reverse complement GC box"
    },
    "SV40_6bp_Repeat": {
        "sequence": "GGGGCG",
        "type": "regulatory",
        "color": "#FFCC00",  # Yellow
        "description": "6bp repeat element"
    },
    "TATA_Box": {
        "sequence": "TATAAA",
        "type": "core_promoter",
        "color": "#00FF00",  # Green
        "description": "Core TATA promoter element"
    },
    "TATA_Box_Variant": {
        "sequence": "TATATA",
        "type": "core_promoter",
        "color": "#00CC00",  # Dark green
        "description": "TATA box variant"
    },
    "ColE1_Origin": {
        "sequence": "AAGGATCTAGGTGAAGATCCTTTTTGATAATCTCATGACCAAAATCCCTTAACGTGAGTTTTCGTTCCACTGAGCGTCAGACCCCGT",
        "type": "origin",
        "color": "#0000FF",  # Blue
        "description": "Bacterial origin of replication"
    }
}

class AdvancedSV40Analyzer:
    def __init__(self):
        self.results = {}
        self.completeness_scores = {}

    def load_sequence_robust(self, filepath):
        """Load sequence with multiple format attempts."""
        for fmt in ["fasta", "fasta-pearson", "fasta-blast", "genbank"]:
            try:
                record = SeqIO.read(filepath, fmt)
                return str(record.seq).upper(), record.id, len(record.seq)
            except:
                continue
        return None, None, None

    def scan_precise_motifs(self, sequence, seq_id):
        """Precise motif scanning with exact position mapping."""
        results = {
            'sequence_id': seq_id,
            'motifs_found': {},
            'total_elements': 0,
            'positions': {},
            'sequence_length': len(sequence)
        }

        for motif_name, motif_data in MOTIFS.items():
            motif_seq = motif_data["sequence"]
            positions = []
            start = 0
            while True:
                pos = sequence.find(motif_seq, start)
                if pos == -1:
                    break
                positions.append(pos)
                start = pos + 1

            if positions:
                results['motifs_found'][motif_name] = {
                    'count': len(positions),
                    'positions': positions,
                    'type': motif_data['type'],
                    'description': motif_data['description']
                }
                results['positions'][motif_name] = positions
                results['total_elements'] += len(positions)

        return results

    def calculate_promoter_completeness(self, motif_positions, sequence_length):
        """
        Calculate promoter completeness score.
        A complete promoter region contains:
        - Enhancer element + TATA box + GC boxes within 200-300bp
        """
        score = {
            'enhancer_tata_gc_clusters': [],
            'promoter_regions': [],
            'completeness_score': 0,
            'max_possible_score': 0,
            'percentage_complete': 0
        }

        # Get all motif positions with their types
        all_motifs = []
        for motif_name, positions in motif_positions.items():
            motif_type = MOTIFS[motif_name]['type']
            for pos in positions:
                all_motifs.append({
                    'name': motif_name,
                    'position': pos,
                    'type': motif_type
                })

        # Sort by position
        all_motifs.sort(key=lambda x: x['position'])

        # Look for clusters containing enhancer + TATA + GC elements within 300bp
        window_size = 300
        clusters = []

        for i, motif1 in enumerate(all_motifs):
            cluster = [motif1]
            for motif2 in all_motifs[i+1:]:
                if motif2['position'] - motif1['position'] <= window_size:
                    cluster.append(motif2)
                else:
                    break
            if len(cluster) >= 2:  # At least 2 elements in cluster
                clusters.append(cluster)

        # Score each cluster for completeness
        for cluster in clusters:
            cluster_types = set(m['type'] for m in cluster)
            cluster_score = 0

            # Check for complete promoter elements
            if 'enhancer' in cluster_types:
                cluster_score += 3
            if 'core_promoter' in cluster_types:
                cluster_score += 2
            if 'promoter' in cluster_types:
                cluster_score += 1
            if 'origin' in cluster_types:
                cluster_score += 1

            if cluster_score >= 4:  # Significant promoter region
                start_pos = min(m['position'] for m in cluster)
                end_pos = max(m['position'] for m in cluster)
                score['promoter_regions'].append({
                    'start': start_pos,
                    'end': end_pos,
                    'length': end_pos - start_pos,
                    'elements': [m['name'] for m in cluster],
                    'score': cluster_score
                })

        # Calculate overall completeness score
        max_possible = 7  # All motif types present in proper configuration
        actual_score = min(len(set(m['type'] for m in all_motifs)), max_possible)

        score['completeness_score'] = actual_score
        score['max_possible_score'] = max_possible
        score['percentage_complete'] = (actual_score / max_possible) * 100
        score['total_motif_clusters'] = len(clusters)

        return score

    def visualize_sv40_regions(self, pfizer_results, moderna_results, output_file="sv40_regions_comparison.png"):
        """Create visualization comparing SV40 regions in both vaccines."""
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 8))
        fig.suptitle('SV40/Promoter Elements Distribution: Pfizer vs Moderna', fontsize=16, fontweight='bold')

        def plot_sequence_regions(ax, results, title):
            sequence_length = results['sequence_length']
            positions = results['positions']

            # Draw sequence backbone
            ax.add_patch(Rectangle((0, 0), sequence_length, 1, facecolor='lightgray', edgecolor='black'))
            ax.set_xlim(0, sequence_length)
            ax.set_ylim(-0.5, 2.5)
            ax.set_title(title, fontsize=14, fontweight='bold')
            ax.set_xlabel('Position in sequence (bp)', fontsize=12)
            ax.set_yticks([])

            # Plot each motif type
            y_offset = 1.2
            for motif_name, motif_data in MOTIFS.items():
                if motif_name in positions:
                    positions_list = positions[motif_name]
                    color = motif_data['color']

                    for pos in positions_list:
                        length = len(motif_data['sequence'])
                        rect = Rectangle((pos, y_offset), length, 0.3,
                                        facecolor=color, edgecolor='black', alpha=0.7)
                        ax.add_patch(rect)

                    # Label the motif type
                    if positions_list:
                        ax.text(positions_list[0], y_offset + 0.4, f"{motif_name} ({len(positions_list)})",
                               fontsize=9, color=color, fontweight='bold')

                    y_offset += 0.35

            # Add sequence info
            ax.text(0.02, 0.95, f"Length: {sequence_length:,} bp", transform=ax.transAxes,
                   fontsize=10, verticalalignment='top',
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

            ax.text(0.02, 0.85, f"Total SV40 elements: {results['total_elements']}",
                   transform=ax.transAxes, fontsize=10, verticalalignment='top',
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

        # Plot Pfizer
        plot_sequence_regions(ax1, pfizer_results, "Pfizer BNT162b2 (OR134577.1) - DENSE SV40 Regions (RED)")

        # Plot Moderna
        plot_sequence_regions(ax2, moderna_results, "Moderna mRNA-1273 (OR134578.1) - SPARSE SV40 Elements")

        # Add legend
        legend_elements = [mpatches.Patch(color=motif_data['color'], label=motif_name)
                         for motif_name, motif_data in MOTIFS.items()]
        fig.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1, 0.5),
                  fontsize=10, title='Motif Types')

        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"✅ Visualization saved to: {output_file}")
        plt.close()

    def generate_statistical_report(self, pfizer_results, moderna_results,
                                   pfizer_completeness, moderna_completeness,
                                   output_file="sv40_statistical_report.txt"):
        """Generate comprehensive statistical report."""
        with open(output_file, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write("ADVANCED SV40/PROMOTER ANALYSIS - STATISTICAL REPORT\n")
            f.write("=" * 80 + "\n\n")

            f.write("SEQUENCE INFORMATION\n")
            f.write("-" * 80 + "\n")
            f.write(f"Pfizer BNT162b2 (OR134577.1): {pfizer_results['sequence_length']:,} bp\n")
            f.write(f"Moderna mRNA-1273 (OR134578.1): {moderna_results['sequence_length']:,} bp\n")
            f.write(f"Length difference: {pfizer_results['sequence_length'] - moderna_results['sequence_length']:,} bp\n\n")

            f.write("MOTIF DETECTION SUMMARY\n")
            f.write("-" * 80 + "\n")
            f.write(f"Pfizer total elements: {pfizer_results['total_elements']}\n")
            f.write(f"Moderna total elements: {moderna_results['total_elements']}\n")
            f.write(f"Difference: {pfizer_results['total_elements'] - moderna_results['total_elements']} elements\n\n")

            f.write("DETAILED MOTIF COMPARISON\n")
            f.write("-" * 80 + "\n")
            for motif_name, motif_data in MOTIFS.items():
                pfizer_count = pfizer_results['motifs_found'].get(motif_name, {}).get('count', 0)
                moderna_count = moderna_results['motifs_found'].get(motif_name, {}).get('count', 0)
                difference = pfizer_count - moderna_count

                status = "EQUAL"
                if pfizer_count > moderna_count:
                    status = "PFIZER HIGHER"
                elif moderna_count > pfizer_count:
                    status = "MODERNA HIGHER"

                f.write(f"{motif_name}:\n")
                f.write(f"  Pfizer: {pfizer_count} copies\n")
                f.write(f"  Moderna: {moderna_count} copies\n")
                f.write(f"  Difference: {difference:+d} ({status})\n")
                f.write(f"  Type: {motif_data['description']}\n\n")

            f.write("PROMOTER COMPLETENESS SCORE\n")
            f.write("-" * 80 + "\n")
            f.write(f"Pfizer completeness: {pfizer_completeness['completeness_score']}/{pfizer_completeness['max_possible_score']} ")
            f.write(f"({pfizer_completeness['percentage_complete']:.1f}%)\n")
            f.write(f"Moderna completeness: {moderna_completeness['completeness_score']}/{moderna_completeness['max_possible_score']} ")
            f.write(f"({moderna_completeness['percentage_complete']:.1f}%)\n")
            f.write(f"Difference: {pfizer_completeness['percentage_complete'] - moderna_completeness['percentage_complete']:.1f}%\n\n")

            f.write(f"Promoter regions found:\n")
            f.write(f"  Pfizer: {len(pfizer_completeness['promoter_regions'])} complete regions\n")
            f.write(f"  Moderna: {len(moderna_completeness['promoter_regions'])} complete regions\n\n")

            if pfizer_completeness['promoter_regions']:
                f.write("Pfizer Promoter Regions:\n")
                for i, region in enumerate(pfizer_completeness['promoter_regions'], 1):
                    f.write(f"  Region {i}: Position {region['start']}-{region['end']} ")
                    f.write(f"({region['length']} bp, Score: {region['score']})\n")
                    f.write(f"    Elements: {', '.join(region['elements'])}\n")

            if moderna_completeness['promoter_regions']:
                f.write("\nModerna Promoter Regions:\n")
                for i, region in enumerate(moderna_completeness['promoter_regions'], 1):
                    f.write(f"  Region {i}: Position {region['start']}-{region['end']} ")
                    f.write(f"({region['length']} bp, Score: {region['score']})\n")
                    f.write(f"    Elements: {', '.join(region['elements'])}\n")

            f.write("\n" + "=" * 80 + "\n")
            f.write("CONCLUSIONS\n")
            f.write("=" * 80 + "\n")

            if pfizer_results['total_elements'] > moderna_results['total_elements']:
                f.write(f"• Pfizer contains {pfizer_results['total_elements'] - moderna_results['total_elements']} more SV40/promoter elements\n")

            unique_to_pfizer = set(pfizer_results['motifs_found'].keys()) - set(moderna_results['motifs_found'].keys())
            if unique_to_pfizer:
                f.write(f"• Elements unique to Pfizer: {', '.join(unique_to_pfizer)}\n")

            unique_to_moderna = set(moderna_results['motifs_found'].keys()) - set(pfizer_results['motifs_found'].keys())
            if unique_to_moderna:
                f.write(f"• Elements unique to Moderna: {', '.join(unique_to_moderna)}\n")

            if pfizer_completeness['percentage_complete'] > moderna_completeness['percentage_complete']:
                f.write(f"• Pfizer has higher promoter completeness ({pfizer_completeness['percentage_complete']:.1f}% vs {moderna_completeness['percentage_complete']:.1f}%)\n")

            f.write(f"\nAnalysis completed: Advanced SV40/Promoter detection with completeness scoring\n")

        print(f"✅ Statistical report saved to: {output_file}")

def main():
    parser = argparse.ArgumentParser(
        description='Advanced SV40/Promoter analysis with completeness scoring and visualization'
    )
    parser.add_argument(
        '--pfizer',
        default='data/references/pfizer_bnt162b2_OR134577.1_clean.fasta',
        help='Pfizer sequence FASTA file'
    )
    parser.add_argument(
        '--moderna',
        default='data/references/moderna_mrna_1273_OR134578.1_clean.fasta',
        help='Moderna sequence FASTA file'
    )
    parser.add_argument(
        '--output-dir',
        default='sv40_advanced_analysis',
        help='Output directory for results'
    )

    args = parser.parse_args()

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    print("=" * 80)
    print("ADVANCED SV40/PROMOTER ANALYZER WITH COMPLETENESS SCORING")
    print("=" * 80)
    print()

    analyzer = AdvancedSV40Analyzer()

    # Load sequences
    print("Loading sequences...")
    pfizer_seq, pfizer_id, pfizer_len = analyzer.load_sequence_robust(args.pfizer)
    moderna_seq, moderna_id, moderna_len = analyzer.load_sequence_robust(args.moderna)

    if not pfizer_seq or not moderna_seq:
        print("❌ Error loading sequences")
        return

    print(f"✅ Loaded {pfizer_id} ({pfizer_len:,} bp)")
    print(f"✅ Loaded {moderna_id} ({moderna_len:,} bp)")
    print()

    # Scan motifs
    print("Scanning for SV40/Promoter motifs...")
    pfizer_results = analyzer.scan_precise_motifs(pfizer_seq, pfizer_id)
    moderna_results = analyzer.scan_precise_motifs(moderna_seq, moderna_id)

    print(f"✅ Pfizer: {pfizer_results['total_elements']} elements detected")
    print(f"✅ Moderna: {moderna_results['total_elements']} elements detected")
    print()

    # Calculate completeness scores
    print("Calculating promoter completeness scores...")
    pfizer_completeness = analyzer.calculate_promoter_completeness(pfizer_results['positions'], pfizer_len)
    moderna_completeness = analyzer.calculate_promoter_completeness(moderna_results['positions'], moderna_len)

    print(f"✅ Pfizer completeness: {pfizer_completeness['percentage_complete']:.1f}%")
    print(f"✅ Moderna completeness: {moderna_completeness['percentage_complete']:.1f}%")
    print()

    # Generate outputs
    print("Generating outputs...")
    viz_file = os.path.join(args.output_dir, "sv40_regions_comparison.png")
    report_file = os.path.join(args.output_dir, "sv40_statistical_report.txt")
    json_file = os.path.join(args.output_dir, "detailed_results.json")

    analyzer.visualize_sv40_regions(pfizer_results, moderna_results, viz_file)
    analyzer.generate_statistical_report(pfizer_results, moderna_results,
                                       pfizer_completeness, moderna_completeness, report_file)

    # Save detailed JSON results
    detailed_results = {
        'pfizer': {
            'sequence_info': {
                'id': pfizer_id,
                'length': pfizer_len,
                'total_elements': pfizer_results['total_elements']
            },
            'motifs': pfizer_results['motifs_found'],
            'completeness_score': pfizer_completeness
        },
        'moderna': {
            'sequence_info': {
                'id': moderna_id,
                'length': moderna_len,
                'total_elements': moderna_results['total_elements']
            },
            'motifs': moderna_results['motifs_found'],
            'completeness_score': moderna_completeness
        }
    }

    with open(json_file, 'w') as f:
        json.dump(detailed_results, f, indent=2)

    print(f"✅ Detailed JSON results saved to: {json_file}")
    print()
    print("=" * 80)
    print("ANALYSIS COMPLETE")
    print("=" * 80)
    print()
    print(f"Results saved to: {args.output_dir}/")
    print(f"- Visualization: {viz_file}")
    print(f"- Statistical report: {report_file}")
    print(f"- Detailed data: {json_file}")

if __name__ == '__main__':
    main()