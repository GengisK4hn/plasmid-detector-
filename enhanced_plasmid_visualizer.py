#!/usr/bin/env python3
"""
Enhanced Plasmid Architecture Visualizer
Clean, publication-quality visualization with proper layout and clarity
"""

from Bio import SeqIO
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle, FancyBboxPatch
import numpy as np
import argparse
import sys
import os
import json
from collections import defaultdict

# Color scheme - clean and distinguishable
COLORS = {
    'SV40_Enhancer': '#D32F2F',      # Strong red
    'SV40_Promoter': '#F44336',       # Light red
    'GC_Box': '#FF9800',             # Orange
    'TATA': '#4CAF50',               # Green
    'Origin': '#2196F3',             # Blue
    'Antibiotic': '#9C27B0',         # Purple
    'T7_Promoter': '#00BCD4',        # Cyan
    'Spike': '#FF5722',              # Deep orange
    'Backbone': '#607D8B',           # Blue gray
    'PolyA': '#9E9E9E'               # Gray
}

ELEMENT_DATA = {
    "SV40_72bp_Enhancer": {
        "sequence": "GGTGTGGAAAGTCCCCAGGCTCCC",
        "category": "SV40_Enhancer",
        "color": COLORS['SV40_Enhancer'],
        "label": "72bp Enhancer",
        "height": 1.0
    },
    "SV40_GC_Box": {
        "sequence": "GGGCGG",
        "category": "GC_Box",
        "color": COLORS['GC_Box'],
        "label": "GC Box",
        "height": 0.6
    },
    "SV40_GC_Box_Reverse": {
        "sequence": "CCGCCC",
        "category": "GC_Box",
        "color": COLORS['GC_Box'],
        "label": "GC Box Rev",
        "height": 0.6
    },
    "TATA_Box": {
        "sequence": "TATAAA",
        "category": "TATA",
        "color": COLORS['TATA'],
        "label": "TATA Box",
        "height": 0.5
    },
    "TATA_Box_Variant": {
        "sequence": "TATATA",
        "category": "TATA",
        "color": COLORS['TATA'],
        "label": "TATA Var",
        "height": 0.5
    },
    "ColE1_Origin": {
        "sequence": "AAGGATCTAGGTGAAGATCCTTTTTGATAATCTCATGACCAAAATCCCTTAACGTGAGTTTTCGTTCCACTGAGCGTCAGACCCCGT",
        "category": "Origin",
        "color": COLORS['Origin'],
        "label": "ColE1 Ori",
        "height": 1.2
    },
    "T7_Promoter": {
        "sequence": "TAATACGACTCACTATA",
        "category": "T7_Promoter",
        "color": COLORS['T7_Promoter'],
        "label": "T7 Promoter",
        "height": 0.8
    },
    "Kanamycin_Resistance": {
        "sequence": "ATGAGCCATATTCAACGGGAAACGTCTTGCTCGAGGCGGGG",
        "category": "Antibiotic",
        "color": COLORS['Antibiotic'],
        "label": "KanR",
        "height": 1.0
    }
}

def load_sequence(filepath):
    """Load sequence robustly."""
    for fmt in ["fasta", "fasta-pearson", "fasta-blast", "genbank"]:
        try:
            record = SeqIO.read(filepath, fmt)
            return str(record.seq).upper(), record.id, len(record.seq)
        except:
            continue
    return None, None, None

def scan_elements(sequence, seq_id):
    """Scan for all elements and return organized results."""
    results = {
        'sequence_id': seq_id,
        'length': len(sequence),
        'elements': []
    }

    for element_name, element_data in ELEMENT_DATA.items():
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
                'height': element_data['height']
            })

    # Sort by position
    results['elements'].sort(key=lambda x: x['start'])
    return results

def create_clean_visualization(pfizer_results, moderna_results, output_file):
    """Create clean, publication-quality visualization."""

    fig = plt.figure(figsize=(20, 12))
    gs = fig.add_gridspec(3, 1, height_ratios=[1, 1, 0.3], hspace=0.3)

    ax_pfizer = fig.add_subplot(gs[0])
    ax_moderna = fig.add_subplot(gs[1])
    ax_legend = fig.add_subplot(gs[2])

    def plot_plasmid(ax, results, title):
        """Plot a single plasmid with clean layout."""
        seq_length = results['length']
        elements = results['elements']

        # Set up the plot
        ax.set_xlim(0, seq_length)
        ax.set_ylim(-0.5, 3.5)
        ax.set_title(title, fontsize=16, fontweight='bold', pad=20)
        ax.set_xlabel('Position (base pairs)', fontsize=12, fontweight='bold')
        ax.set_yticks([])

        # Draw main sequence backbone
        backbone = Rectangle((0, 1.0), seq_length, 0.3,
                            facecolor='#ECEFF1', edgecolor='#455A64',
                            linewidth=2)
        ax.add_patch(backbone)

        # Draw scale bar
        for i in range(0, seq_length, 1000):
            ax.axvline(x=i, ymin=0.1, ymax=0.9, color='#B0BEC5',
                      linestyle='--', linewidth=0.5, alpha=0.5)
            if i % 2000 == 0:
                ax.text(i, -0.3, f'{i:,}', fontsize=8,
                       ha='center', color='#546E7A')

        # Group elements by category for cleaner display
        category_groups = defaultdict(list)
        for elem in elements:
            category_groups[elem['category']].append(elem)

        # Draw elements with better organization
        y_positions = {
            'SV40_Enhancer': 2.5,
            'GC_Box': 2.2,
            'TATA': 1.9,
            'Origin': 1.6,
            'T7_Promoter': 1.3,
            'Antibiotic': 1.0
        }

        # Draw elements
        for category, elems in category_groups.items():
            y_pos = y_positions.get(category, 1.0)

            for elem in elems:
                width = elem['end'] - elem['start']
                height = elem['height']

                # Use fancy box for better appearance
                box = FancyBboxPatch((elem['start'], y_pos - height/2),
                                   width, height,
                                   boxstyle="round,pad=0.02",
                                   facecolor=elem['color'],
                                   edgecolor='white',
                                   linewidth=1.5,
                                   alpha=0.9)
                ax.add_patch(box)

                # Add label for first occurrence only
                if elem == elems[0]:
                    ax.text(elem['start'], y_pos + height/2 + 0.1,
                           elem['label'],
                           fontsize=9, fontweight='bold',
                           color=elem['color'],
                           ha='left')

        # Add statistics box
        stats_text = f"Length: {seq_length:,} bp\n"
        stats_text += f"Elements: {len(elements)}\n"

        # Count by category
        category_counts = defaultdict(int)
        for elem in elements:
            category_counts[elem['category']] += 1

        stats_text += "\nCategory Breakdown:\n"
        for cat, count in sorted(category_counts.items(), key=lambda x: -x[1]):
            stats_text += f"  {cat}: {count}\n"

        ax.text(0.98, 0.95, stats_text,
               transform=ax.transAxes,
               fontsize=10,
               verticalalignment='top',
               horizontalalignment='right',
               bbox=dict(boxstyle='round',
                        facecolor='white',
                        edgecolor='#90A4AE',
                        alpha=0.95,
                        linewidth=2))

    # Plot Pfizer
    pfizer_title = ("Pfizer BNT162b2 (OR134577.1)\n"
                   "Dense SV40 regulatory regions + Complete backbone")
    plot_plasmid(ax_pfizer, pfizer_results, pfizer_title)

    # Plot Moderna
    moderna_title = ("Moderna mRNA-1273 (OR134578.1)\n"
                    "Simplified architecture + Reduced regulatory elements")
    plot_plasmid(ax_moderna, moderna_results, moderna_title)

    # Create clean legend
    ax_legend.axis('off')

    legend_elements = []
    categories_seen = set()

    for elem in pfizer_results['elements'] + moderna_results['elements']:
        if elem['category'] not in categories_seen:
            legend_elements.append(mpatches.Patch(
                color=elem['color'],
                label=elem['category'].replace('_', ' '),
                linewidth=2
            ))
            categories_seen.add(elem['category'])

    ax_legend.legend(handles=legend_elements,
                    loc='center',
                    ncol=6,
                    fontsize=11,
                    frameon=True,
                    fancybox=True,
                    shadow=True,
                    facecolor='white',
                    edgecolor='#90A4AE')

    # Add main title
    fig.suptitle('Comprehensive Plasmid Architecture Comparison',
                fontsize=18, fontweight='bold', y=0.98)

    plt.savefig(output_file, dpi=300, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    print(f"✅ Enhanced visualization saved: {output_file}")
    plt.close()

def create_summary_statistics(pfizer_results, moderna_results, output_file):
    """Create comprehensive statistics summary."""

    with open(output_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("COMPREHENSIVE PLASMID ARCHITECTURE ANALYSIS\n")
        f.write("=" * 80 + "\n\n")

        f.write("SEQUENCE METRICS\n")
        f.write("-" * 80 + "\n")
        f.write(f"Pfizer: {pfizer_results['length']:,} bp\n")
        f.write(f"Moderna: {moderna_results['length']:,} bp\n")
        f.write(f"Difference: {pfizer_results['length'] - moderna_results['length']:,} bp\n\n")

        f.write("ELEMENT ANALYSIS\n")
        f.write("-" * 80 + "\n")
        f.write(f"Pfizer: {len(pfizer_results['elements'])} elements\n")
        f.write(f"Moderna: {len(moderna_results['elements'])} elements\n")
        f.write(f"Difference: {len(pfizer_results['elements']) - len(moderna_results['elements'])} elements\n\n")

        f.write("CATEGORY BREAKDOWN\n")
        f.write("-" * 80 + "\n")

        # Get all categories
        all_categories = set()
        for elem in pfizer_results['elements']:
            all_categories.add(elem['category'])
        for elem in moderna_results['elements']:
            all_categories.add(elem['category'])

        for category in sorted(all_categories):
            pfizer_count = sum(1 for e in pfizer_results['elements'] if e['category'] == category)
            moderna_count = sum(1 for e in moderna_results['elements'] if e['category'] == category)
            diff = pfizer_count - moderna_count

            f.write(f"\n{category.replace('_', ' ')}:\n")
            f.write(f"  Pfizer: {pfizer_count}\n")
            f.write(f"  Moderna: {moderna_count}\n")
            if diff > 0:
                f.write(f"  Pfizer has {diff} more\n")
            elif diff < 0:
                f.write(f"  Moderna has {-diff} more\n")
            else:
                f.write(f"  Equal\n")

        f.write("\n" + "=" * 80 + "\n")
        f.write("POSITION MAP\n")
        f.write("=" * 80 + "\n\n")

        f.write("PFIZER ELEMENTS:\n")
        for elem in pfizer_results['elements']:
            f.write(f"  {elem['label']:15} {elem['category']:20} "
                   f"pos {elem['start']:5}-{elem['end']:5}\n")

        f.write("\nMODERNA ELEMENTS:\n")
        for elem in moderna_results['elements']:
            f.write(f"  {elem['label']:15} {elem['category']:20} "
                   f"pos {elem['start']:5}-{elem['end']:5}\n")

    print(f"✅ Statistics saved: {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Enhanced plasmid visualizer')
    parser.add_argument('--pfizer', default='data/references/pfizer_bnt162b2_OR134577.1_clean.fasta')
    parser.add_argument('--moderna', default='data/references/moderna_mrna_1273_OR134578.1_clean.fasta')
    parser.add_argument('--output-dir', default='enhanced_plasmid_visualization')

    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    print("=" * 80)
    print("ENHANCED PLASMID ARCHITECTURE VISUALIZER")
    print("=" * 80)
    print()

    # Load sequences
    print("Loading sequences...")
    pfizer_seq, pfizer_id, pfizer_len = load_sequence(args.pfizer)
    moderna_seq, moderna_id, moderna_len = load_sequence(args.moderna)

    print(f"✅ Pfizer: {pfizer_id} ({pfizer_len:,} bp)")
    print(f"✅ Moderna: {moderna_id} ({moderna_len:,} bp)")
    print()

    # Scan elements
    print("Scanning plasmid elements...")
    pfizer_results = scan_elements(pfizer_seq, pfizer_id)
    moderna_results = scan_elements(moderna_seq, moderna_id)

    print(f"✅ Pfizer: {len(pfizer_results['elements'])} elements")
    print(f"✅ Moderna: {len(moderna_results['elements'])} elements")
    print()

    # Generate outputs
    viz_file = os.path.join(args.output_dir, "clean_plasmid_architecture.png")
    stats_file = os.path.join(args.output_dir, "plasmid_statistics.txt")

    print("Generating enhanced visualization...")
    create_clean_visualization(pfizer_results, moderna_results, viz_file)

    print("Generating statistics...")
    create_summary_statistics(pfizer_results, moderna_results, stats_file)

    print()
    print("=" * 80)
    print("VISUALIZATION COMPLETE")
    print("=" * 80)
    print()
    print(f"Outputs saved to: {args.output_dir}/")
    print(f"  - {viz_file}")
    print(f"  - {stats_file}")

if __name__ == '__main__':
    main()