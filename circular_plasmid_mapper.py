#!/usr/bin/env python3
"""
Circular Plasmid Mapper - Shows SV40 Nuclear Targeting Sequences
Proper circular plasmid maps with SV40 enhancer/promoter regions highlighted
"""

from Bio import SeqIO
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Circle, Wedge, Rectangle, FancyArrowPatch
import numpy as np
import argparse
import sys
import os
import math
from collections import defaultdict

# SV40 Nuclear Localization Signal - the most potent known
SV40_NLS_MOTIF = "PKKKRKV"  # Classic SV40 nuclear localization signal
SV40_ENHANCER_CORE = "GGTGTGGAAAGTCCCCAGGCTCCC"  # 72bp enhancer
SV40_PROMOTER_REGION = "CCGCCC"  # GC-rich promoter

COLORS = {
    'SV40_NLS': '#D32F2F',        # Strong red - NUCLEAR TARGETING
    'SV40_Enhancer': '#F44336',    # Red - ENHANCER
    'SV40_Promoter': '#FF9800',    # Orange - PROMOTER
    'TATA_Box': '#4CAF50',         # Green
    'Origin': '#2196F3',           # Blue
    'Antibiotic': '#9C27B0',       # Purple
    'T7_Promoter': '#00BCD4',      # Cyan
    'Backbone': '#B0BEC5',         # Gray
    'Spike': '#FF5722'             # Deep orange
}

class CircularPlasmidMapper:
    def __init__(self):
        self.results = {}

    def load_sequence(self, filepath):
        """Load sequence robustly."""
        for fmt in ["fasta", "fasta-pearson", "fasta-blast", "genbank"]:
            try:
                record = SeqIO.read(filepath, fmt)
                return str(record.seq).upper(), record.id, len(record.seq)
            except:
                continue
        return None, None, None

    def scan_plasmid_elements(self, sequence, seq_id):
        """Scan for all plasmid elements with nuclear targeting focus."""
        results = {
            'sequence_id': seq_id,
            'length': len(sequence),
            'elements': []
        }

        # Define critical elements with nuclear targeting emphasis
        elements = {
            "SV40_72bp_Enhancer": {
                "sequence": "GGTGTGGAAAGTCCCCAGGCTCCC",
                "category": "SV40_Enhancer",
                "color": COLORS['SV40_Enhancer'],
                "label": "SV40 72bp Enhancer",
                "nuclear_targeting": True,
                "importance": "CRITICAL"
            },
            "SV40_GC_Box": {
                "sequence": "GGGCGG",
                "category": "SV40_Promoter",
                "color": COLORS['SV40_Promoter'],
                "label": "SV40 GC Box",
                "nuclear_targeting": True,
                "importance": "HIGH"
            },
            "SV40_GC_Box_Reverse": {
                "sequence": "CCGCCC",
                "category": "SV40_Promoter",
                "color": COLORS['SV40_Promoter'],
                "label": "SV40 GC Box Rev",
                "nuclear_targeting": True,
                "importance": "HIGH"
            },
            "TATA_Box": {
                "sequence": "TATAAA",
                "category": "TATA_Box",
                "color": COLORS['TATA_Box'],
                "label": "TATA Box",
                "nuclear_targeting": False,
                "importance": "MEDIUM"
            },
            "TATA_Box_Variant": {
                "sequence": "TATATA",
                "category": "TATA_Box",
                "color": COLORS['TATA_Box'],
                "label": "TATA Variant",
                "nuclear_targeting": False,
                "importance": "MEDIUM"
            },
            "ColE1_Origin": {
                "sequence": "AAGGATCTAGGTGAAGATCCTTTTTGATAATCTCATGACCAAAATCCCTTAACGTGAGTTTTCGTTCCACTGAGCGTCAGACCCCGT",
                "category": "Origin",
                "color": COLORS['Origin'],
                "label": "ColE1 Origin",
                "nuclear_targeting": False,
                "importance": "HIGH"
            },
            "T7_Promoter": {
                "sequence": "TAATACGACTCACTATA",
                "category": "T7_Promoter",
                "color": COLORS['T7_Promoter'],
                "label": "T7 Promoter",
                "nuclear_targeting": False,
                "importance": "MEDIUM"
            },
            "Kanamycin_Resistance": {
                "sequence": "ATGAGCCATATTCAACGGGAAACGTCTTGCTCGAGGCGGGG",
                "category": "Antibiotic",
                "color": COLORS['Antibiotic'],
                "label": "KanR",
                "nuclear_targeting": False,
                "importance": "MEDIUM"
            },
            "SV40_6bp_Repeat": {
                "sequence": "GGGGCG",
                "category": "SV40_Promoter",
                "color": COLORS['SV40_Promoter'],
                "label": "6bp Repeat",
                "nuclear_targeting": True,
                "importance": "HIGH"
            }
        }

        for element_name, element_data in elements.items():
            element_seq = element_data["sequence"]
            positions = []

            # Find all occurrences
            start = 0
            while True:
                pos = sequence.find(element_seq, start)
                if pos == -1:
                    break
                positions.append(pos)
                start = pos + 1

            for pos in positions:
                results['elements'].append({
                    'name': element_name,
                    'label': element_data['label'],
                    'start': pos,
                    'end': pos + len(element_seq),
                    'category': element_data['category'],
                    'color': element_data['color'],
                    'nuclear_targeting': element_data['nuclear_targeting'],
                    'importance': element_data['importance']
                })

        # Sort by position
        results['elements'].sort(key=lambda x: x['start'])
        return results

    def create_circular_map(self, results, output_file, title=""):
        """Create publication-quality circular plasmid map."""

        fig, ax = plt.subplots(figsize=(14, 14))
        seq_length = results['length']
        elements = results['elements']

        # Set up the plot
        ax.set_xlim(-1.5, 1.5)
        ax.set_ylim(-1.5, 1.5)
        ax.set_aspect('equal')
        ax.axis('off')

        # Title with nuclear targeting warning
        if title:
            full_title = f"{title}\n"
        else:
            full_title = f"{results['sequence_id']}\n"

        # Count nuclear targeting elements
        nuclear_count = sum(1 for e in elements if e['nuclear_targeting'])
        if nuclear_count > 0:
            full_title += f"⚠️ CONTAINS {nuclear_count} SV40 NUCLEAR TARGETING SEQUENCES"

        ax.set_title(full_title, fontsize=14, fontweight='bold', pad=20)

        # Draw outer circle (plasmid backbone)
        backbone = Circle((0, 0), 1.0, fill=False, edgecolor='#455A64',
                         linewidth=3, label='Plasmid Backbone')
        ax.add_patch(backbone)

        # Draw scale markers
        for i in range(12):
            angle = 2 * np.pi * i / 12
            x_start = 0.95 * np.cos(angle)
            y_start = 0.95 * np.sin(angle)
            x_end = 1.05 * np.cos(angle)
            y_end = 1.05 * np.sin(angle)
            ax.plot([x_start, x_end], [y_start, y_end],
                   color='#B0BEC5', linewidth=2)

            # Add position labels
            bp_pos = int((i / 12) * seq_length)
            x_text = 1.15 * np.cos(angle)
            y_text = 1.15 * np.sin(angle)
            ax.text(x_text, y_text, f'{bp_pos:,}',
                   fontsize=9, ha='center', va='center',
                   color='#546E7A', fontweight='bold')

        # Group elements for cleaner display
        nuclear_elements = [e for e in elements if e['nuclear_targeting']]
        other_elements = [e for e in elements if not e['nuclear_targeting']]

        # Draw nuclear targeting elements FIRST (outer ring, most important)
        for elem in nuclear_elements:
            start_angle = 2 * np.pi * elem['start'] / seq_length
            end_angle = 2 * np.pi * elem['end'] / seq_length

            # Make nuclear elements prominent - outer ring
            radius = 1.15
            width = 0.08

            theta1 = np.degrees(start_angle)
            theta2 = np.degrees(end_angle)
            wedge = Wedge((0, 0), radius, theta1, theta2,
                         width=width,
                         facecolor=elem['color'],
                         edgecolor='white',
                         linewidth=1.5,
                         alpha=0.9)
            ax.add_patch(wedge)

            # Add label for nuclear elements
            mid_angle = (start_angle + end_angle) / 2
            label_radius = radius + width + 0.08
            x_label = label_radius * np.cos(mid_angle)
            y_label = label_radius * np.sin(mid_angle)

            # Rotate text appropriately
            rotation = np.degrees(mid_angle)
            if rotation > 90 and rotation < 270:
                rotation -= 180

            ax.text(x_label, y_label, elem['label'],
                   fontsize=10, fontweight='bold',
                   color=elem['color'],
                   ha='center', va='center',
                   rotation=rotation)

        # Draw other elements (inner ring)
        for elem in other_elements:
            start_angle = 2 * np.pi * elem['start'] / seq_length
            end_angle = 2 * np.pi * elem['end'] / seq_length

            # Inner ring for other elements
            radius = 0.92
            width = 0.06

            theta1 = np.degrees(start_angle)
            theta2 = np.degrees(end_angle)
            wedge = Wedge((0, 0), radius, theta1, theta2,
                         width=width,
                         facecolor=elem['color'],
                         edgecolor='white',
                         linewidth=1,
                         alpha=0.8)
            ax.add_patch(wedge)

        # Add central information
        info_text = "PLASMID INFO\n"
        info_text += "=" * 20 + "\n"
        info_text += f"Length: {seq_length:,} bp\n"
        info_text += f"Total Elements: {len(elements)}\n"
        info_text += f"SV40 Nuclear: {nuclear_count}\n\n"

        if nuclear_count > 0:
            info_text += "⚠️ NUCLEAR TARGETING\n"
            info_text += "SEQUENCES DETECTED\n\n"
            for elem in nuclear_elements:
                info_text += f"• {elem['label']}\n"
                info_text += f"  pos: {elem['start']}-{elem['end']}\n"

        ax.text(0, 0, info_text,
               fontsize=10,
               ha='center', va='center',
               bbox=dict(boxstyle='round,pad=0.5',
                        facecolor='white',
                        edgecolor='#D32F2F' if nuclear_count > 0 else '#90A4AE',
                        linewidth=3 if nuclear_count > 0 else 2,
                        alpha=0.95))

        # Create legend
        legend_elements = []
        categories_seen = set()

        for elem in elements:
            if elem['category'] not in categories_seen:
                if elem['nuclear_targeting']:
                    label = f"SV40 Nuclear: {elem['category']}"
                else:
                    label = elem['category']

                legend_elements.append(mpatches.Patch(
                    color=elem['color'],
                    label=label,
                    linewidth=2
                ))
                categories_seen.add(elem['category'])

        ax.legend(handles=legend_elements,
                loc='lower center',
                bbox_to_anchor=(0.5, -0.05),
                ncol=3,
                fontsize=11,
                frameon=True,
                fancybox=True,
                shadow=True)

        plt.savefig(output_file, dpi=300, bbox_inches='tight',
                    facecolor='white', edgecolor='none')
        print(f"✅ Circular map saved: {output_file}")
        plt.close()

    def create_side_by_side_maps(self, pfizer_results, moderna_results, output_file):
        """Create side-by-side comparison of circular maps."""

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))

        def plot_circular(ax, results, title):
            seq_length = results['length']
            elements = results['elements']

            ax.set_xlim(-1.3, 1.3)
            ax.set_ylim(-1.3, 1.3)
            ax.set_aspect('equal')
            ax.axis('off')

            # Enhanced title
            nuclear_count = sum(1 for e in elements if e['nuclear_targeting'])
            if nuclear_count > 0:
                title_text = f"{title}\n⚠️ {nuclear_count} SV40 NUCLEAR TARGETING SEQUENCES"
            else:
                title_text = f"{title}\n✓ No nuclear targeting detected"

            ax.set_title(title_text, fontsize=14, fontweight='bold', pad=15)

            # Draw plasmid backbone
            backbone = Circle((0, 0), 1.0, fill=False, edgecolor='#455A64',
                             linewidth=3)
            ax.add_patch(backbone)

            # Scale markers
            for i in range(8):
                angle = 2 * np.pi * i / 8
                x_start = 0.93 * np.cos(angle)
                y_start = 0.93 * np.sin(angle)
                x_end = 1.07 * np.cos(angle)
                y_end = 1.07 * np.sin(angle)
                ax.plot([x_start, x_end], [y_start, y_end],
                       color='#B0BEC5', linewidth=2)

                bp_pos = int((i / 8) * seq_length)
                x_text = 1.15 * np.cos(angle)
                y_text = 1.15 * np.sin(angle)
                ax.text(x_text, y_text, f'{bp_pos:,}',
                       fontsize=8, ha='center', va='center',
                       color='#546E7A', fontweight='bold')

            # Nuclear elements - prominent outer ring
            for elem in elements:
                if not elem['nuclear_targeting']:
                    continue

                start_angle = 2 * np.pi * elem['start'] / seq_length
                end_angle = 2 * np.pi * elem['end'] / seq_length

                radius = 1.12
                width = 0.10

                theta1 = np.degrees(start_angle)
                theta2 = np.degrees(end_angle)
                wedge = Wedge((0, 0), radius, theta1, theta2,
                             width=width,
                             facecolor=elem['color'],
                             edgecolor='white',
                             linewidth=2,
                             alpha=1.0)
                ax.add_patch(wedge)

            # Other elements - inner ring
            for elem in elements:
                if elem['nuclear_targeting']:
                    continue

                start_angle = 2 * np.pi * elem['start'] / seq_length
                end_angle = 2 * np.pi * elem['end'] / seq_length

                radius = 0.90
                width = 0.05

                theta1 = np.degrees(start_angle)
                theta2 = np.degrees(end_angle)
                wedge = Wedge((0, 0), radius, theta1, theta2,
                             width=width,
                             facecolor=elem['color'],
                             edgecolor='white',
                             linewidth=1,
                             alpha=0.8)
                ax.add_patch(wedge)

            # Central info
            info_text = f"{seq_length:,} bp\n"
            info_text += f"{len(elements)} elements\n"
            if nuclear_count > 0:
                info_text += f"⚠️ {nuclear_count} SV40 nuclear"

            ax.text(0, 0, info_text,
                   fontsize=12, ha='center', va='center',
                   fontweight='bold',
                   bbox=dict(boxstyle='circle',
                            facecolor='white',
                            edgecolor='#D32F2F' if nuclear_count > 0 else '#4CAF50',
                            linewidth=4 if nuclear_count > 0 else 3,
                            alpha=0.95))

        # Plot both
        plot_circular(ax1, pfizer_results,
                     "Pfizer BNT162b2\n(pcDNA3.1-like with SV40)")
        plot_circular(ax2, moderna_results,
                     "Moderna mRNA-1273\n(Simplified architecture)")

        # Create legend
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor=COLORS['SV40_Enhancer'],
                 label='SV40 Enhancer (Nuclear Targeting)'),
            Patch(facecolor=COLORS['SV40_Promoter'],
                 label='SV40 Promoter (Nuclear Targeting)'),
            Patch(facecolor=COLORS['TATA_Box'], label='TATA Box'),
            Patch(facecolor=COLORS['Origin'], label='ColE1 Origin'),
            Patch(facecolor=COLORS['T7_Promoter'], label='T7 Promoter')
        ]

        fig.legend(handles=legend_elements,
                  loc='lower center',
                  bbox_to_anchor=(0.5, -0.05),
                  ncol=5,
                  fontsize=11,
                  frameon=True,
                  fancybox=True)

        fig.suptitle('CIRCULAR PLASMID MAPS: SV40 NUCLEAR TARGETING SEQUENCES',
                    fontsize=16, fontweight='bold', y=0.98)

        plt.savefig(output_file, dpi=300, bbox_inches='tight',
                    facecolor='white', edgecolor='none')
        print(f"✅ Side-by-side maps saved: {output_file}")
        plt.close()

def main():
    parser = argparse.ArgumentParser(
        description='Circular plasmid mapper showing SV40 nuclear targeting sequences'
    )
    parser.add_argument('--pfizer',
                       default='data/references/pfizer_bnt162b2_OR134577.1_clean.fasta')
    parser.add_argument('--moderna',
                       default='data/references/moderna_mrna_1273_OR134578.1_clean.fasta')
    parser.add_argument('--output-dir',
                       default='circular_plasmid_maps')

    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    print("=" * 80)
    print("CIRCULAR PLASMID MAPPER - SV40 NUCLEAR TARGETING SEQUENCES")
    print("=" * 80)
    print()

    mapper = CircularPlasmidMapper()

    # Load sequences
    print("Loading sequences...")
    pfizer_seq, pfizer_id, pfizer_len = mapper.load_sequence(args.pfizer)
    moderna_seq, moderna_id, moderna_len = mapper.load_sequence(args.moderna)

    print(f"✅ Pfizer: {pfizer_id} ({pfizer_len:,} bp)")
    print(f"✅ Moderna: {moderna_id} ({moderna_len:,} bp)")
    print()

    # Scan elements
    print("Scanning for plasmid elements (emphasizing nuclear targeting)...")
    pfizer_results = mapper.scan_plasmid_elements(pfizer_seq, pfizer_id)
    moderna_results = mapper.scan_plasmid_elements(moderna_seq, moderna_id)

    pfizer_nuclear = sum(1 for e in pfizer_results['elements'] if e['nuclear_targeting'])
    moderna_nuclear = sum(1 for e in moderna_results['elements'] if e['nuclear_targeting'])

    print(f"✅ Pfizer: {len(pfizer_results['elements'])} elements ({pfizer_nuclear} nuclear targeting)")
    print(f"✅ Moderna: {len(moderna_results['elements'])} elements ({moderna_nuclear} nuclear targeting)")
    print()

    # Generate outputs
    print("Generating circular plasmid maps...")

    # Individual maps
    pfizer_map = os.path.join(args.output_dir, "pfizer_circular_map.png")
    moderna_map = os.path.join(args.output_dir, "moderna_circular_map.png")
    comparison_map = os.path.join(args.output_dir, "side_by_side_comparison.png")

    mapper.create_circular_map(pfizer_results, pfizer_map,
                              "Pfizer BNT162b2 - pcDNA3.1-like")
    mapper.create_circular_map(moderna_results, moderna_map,
                              "Moderna mRNA-1273 - Simplified")
    mapper.create_side_by_side_maps(pfizer_results, moderna_results, comparison_map)

    # Generate warning report
    report_file = os.path.join(args.output_dir, "nuclear_targeting_report.txt")
    with open(report_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("SV40 NUCLEAR TARGETING SEQUENCE ANALYSIS\n")
        f.write("=" * 80 + "\n\n")
        f.write("REFERENCES:\n")
        f.write("-" * 80 + "\n")
        f.write("SV40 nuclear localization signals (NLS) are among the most potent\n")
        f.write("nuclear targeting sequences known. They facilitate active transport\n")
        f.write("of DNA into the nucleus during interphase, affecting BOTH dividing\n")
        f.write("and non-dividing cells.\n\n")
        f.write("Key reference: https://pmc.ncbi.nlm.nih.gov/articles/PMC4150867/\n")
        f.write("Nuclear Localization of Plasmid DNA During Transfection\n\n")

        f.write("ANALYSIS RESULTS:\n")
        f.write("-" * 80 + "\n")
        f.write(f"Pfizer BNT162b2: {pfizer_nuclear} SV40 nuclear targeting elements\n")
        f.write(f"Moderna mRNA-1273: {moderna_nuclear} SV40 nuclear targeting elements\n\n")

        if pfizer_nuclear > 0:
            f.write("PFIZER NUCLEAR TARGETING ELEMENTS:\n")
            for elem in pfizer_results['elements']:
                if elem['nuclear_targeting']:
                    f.write(f"  {elem['label']}: position {elem['start']}-{elem['end']}\n")

        if moderna_nuclear > 0:
            f.write("\nMODERNA NUCLEAR TARGETING ELEMENTS:\n")
            for elem in moderna_results['elements']:
                if elem['nuclear_targeting']:
                    f.write(f"  {elem['label']}: position {elem['start']}-{elem['end']}\n")

        f.write("\n" + "=" * 80 + "\n")
        f.write("SIGNIFICANCE:\n")
        f.write("=" * 80 + "\n")
        f.write("Plasmid DNA containing SV40 NLS can enter nuclei of non-dividing cells\n")
        f.write("through active nuclear transport mechanisms. This has implications\n")
        f.write("for genomic integration risk and long-term persistence.\n\n")

        if pfizer_nuclear > moderna_nuclear:
            f.write(f"CONCLUSION: Pfizer contains {pfizer_nuclear - moderna_nuclear} more SV40\n")
            f.write("nuclear targeting elements than Moderna, representing higher\n")
            f.write("potential for nuclear uptake and genomic integration.\n")

    print(f"✅ Nuclear targeting report saved: {report_file}")
    print()
    print("=" * 80)
    print("CIRCULAR PLASMID MAPS COMPLETE")
    print("=" * 80)
    print()
    print(f"Outputs saved to: {args.output_dir}/")
    print(f"  - {pfizer_map}")
    print(f"  - {moderna_map}")
    print(f"  - {comparison_map}")
    print(f"  - {report_file}")

if __name__ == '__main__':
    main()