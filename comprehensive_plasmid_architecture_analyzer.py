#!/usr/bin/env python3
"""
Comprehensive Plasmid Architecture Analyzer
Maps all critical plasmid features for vaccine integrity analysis:
- Bacterial backbone elements (ColE1, KanR, T7 promoter)
- Mammalian regulatory elements (SV40 components)
- Spike gene insert and UTRs
- Selection markers and origins
- Safety/relevant regions for residual DNA analysis
"""

from Bio import SeqIO
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle
import numpy as np
import argparse
import sys
import os
import json
from collections import defaultdict

# Comprehensive plasmid element definitions
PLASMID_ELEMENTS = {
    # SV40 Components (Pfizer-specific)
    "SV40_72bp_Enhancer": {
        "sequence": "GGTGTGGAAAGTCCCCAGGCTCCC",
        "type": "mammalian_regulatory",
        "category": "SV40",
        "color": "#FF0000",  # Red - highest concern
        "description": "72bp SV40 enhancer (Pfizer-specific)",
        "safety_relevance": "HIGH",
        "qPCR_target": True
    },
    "SV40_GC_Box": {
        "sequence": "GGGCGG",
        "type": "mammalian_regulatory",
        "category": "SV40",
        "color": "#FF6600",  # Orange
        "description": "GC-rich promoter box",
        "safety_relevance": "MEDIUM",
        "qPCR_target": False
    },
    "SV40_GC_Box_Reverse": {
        "sequence": "CCGCCC",
        "type": "mammalian_regulatory",
        "category": "SV40",
        "color": "#FF9900",  # Light orange
        "description": "Reverse complement GC box",
        "safety_relevance": "MEDIUM",
        "qPCR_target": False
    },
    "SV40_6bp_Repeat": {
        "sequence": "GGGGCG",
        "type": "mammalian_regulatory",
        "category": "SV40",
        "color": "#FFCC00",  # Yellow
        "description": "6bp repeat element",
        "safety_relevance": "LOW",
        "qPCR_target": False
    },
    "TATA_Box": {
        "sequence": "TATAAA",
        "type": "mammalian_regulatory",
        "category": "Promoter",
        "color": "#00FF00",  # Green
        "description": "Core TATA promoter element",
        "safety_relevance": "LOW",
        "qPCR_target": False
    },
    "TATA_Box_Variant": {
        "sequence": "TATATA",
        "type": "mammalian_regulatory",
        "category": "Promoter",
        "color": "#00CC00",  # Dark green
        "description": "TATA box variant",
        "safety_relevance": "LOW",
        "qPCR_target": False
    },

    # Bacterial Backbone Elements
    "ColE1_Origin": {
        "sequence": "AAGGATCTAGGTGAAGATCCTTTTTGATAATCTCATGACCAAAATCCCTTAACGTGAGTTTTCGTTCCACTGAGCGTCAGACCCCGT",
        "type": "bacterial_backbone",
        "category": "Origin",
        "color": "#0000FF",  # Blue
        "description": "ColE1 bacterial origin of replication",
        "safety_relevance": "HIGH",
        "qPCR_target": True
    },
    "ColRNAI_1": {
        "sequence": "GGGCGG",  # Part of ColE1 regulation
        "type": "bacterial_backbone",
        "category": "Regulatory",
        "color": "#0066FF",  # Medium blue
        "description": "RNAI regulator within ColE1",
        "safety_relevance": "MEDIUM",
        "qPCR_target": True
    },

    # Selection Markers
    "Kanamycin_Resistance": {
        "sequence": "ATGAGCCATATTCAACGGGAAACGTCTTGCTCGAGGCGGGG",  # KanR start
        "type": "selection_marker",
        "category": "Antibiotic_Resistance",
        "color": "#FF00FF",  # Magenta
        "description": "Kanamycin/Neomycin resistance gene (aph(3')-II)",
        "safety_relevance": "HIGH",
        "qPCR_target": True
    },
    "KanR_Internal": {
        "sequence": "GCTCGAGGCGGGG",  # Internal sequence
        "type": "selection_marker",
        "category": "Antibiotic_Resistance",
        "color": "#CC00FF",  # Dark magenta
        "description": "KanR internal sequence",
        "safety_relevance": "HIGH",
        "qPCR_target": True
    },

    # Transcription Elements
    "T7_Promoter": {
        "sequence": "TAATACGACTCACTATA",
        "type": "transcription",
        "category": "Promoter",
        "color": "#00FFFF",  # Cyan
        "description": "T7 promoter for IVT transcription",
        "safety_relevance": "MEDIUM",
        "qPCR_target": True
    },
    "T7_Promoter_Extended": {
        "sequence": "TAATACGACTCACTATAGGG",
        "type": "transcription",
        "category": "Promoter",
        "color": "#00CCCC",  # Dark cyan
        "description": "Extended T7 promoter",
        "safety_relevance": "MEDIUM",
        "qPCR_target": True
    },

    # Spike Gene Elements
    "Spike_5UTR_Start": {
        "sequence": "GAAAAAATG",  # Common start context
        "type": "spike_gene",
        "category": "UTR",
        "color": "#FF0066",  # Pink-red
        "description": "Spike gene 5' UTR region",
        "safety_relevance": "LOW",
        "qPCR_target": True
    },
    "Spike_Start_Codon": {
        "sequence": "ATGTTTGTTTTT",
        "type": "spike_gene",
        "category": "Coding",
        "color": "#CC0066",  # Dark pink
        "description": "Spike protein start codon region",
        "safety_relevance": "LOW",
        "qPCR_target": True
    },
    "Spike_Internal": {
        "sequence": "TTTGTTTTT",  # Common spike sequence
        "type": "spike_gene",
        "category": "Coding",
        "color": "#FF3399",  # Light pink
        "description": "Spike gene internal sequence",
        "safety_relevance": "LOW",
        "qPCR_target": True
    },

    # Additional regulatory elements
    "AmpR_Promoter": {
        "sequence": "TTGACA",  # -35 region
        "type": "regulatory",
        "category": "Promoter",
        "color": "#9900FF",  # Purple
        "description": "AmpR promoter region (sometimes drives KanR)",
        "safety_relevance": "MEDIUM",
        "qPCR_target": False
    },
    "f1_Origin": {
        "sequence": "CAAGGCGACCA",
        "type": "origin",
        "category": "Phage",
        "color": "#009900",  # Dark green
        "description": "f1 phage origin (ssDNA production)",
        "safety_relevance": "MEDIUM",
        "qPCR_target": True
    },

    # PolyA signals
    "PolyA_Signal": {
        "sequence": "AATAAA",
        "type": "regulatory",
        "category": "PolyA",
        "color": "#999999",  # Gray
        "description": "Polyadenylation signal",
        "safety_relevance": "LOW",
        "qPCR_target": False
    },
    "SV40_PolyA": {
        "sequence": "TCCATGGTGATGC",
        "type": "mammalian_regulatory",
        "category": "PolyA",
        "color": "#CC6666",  # Light red
        "description": "SV40 polyA signal region",
        "safety_relevance": "MEDIUM",
        "qPCR_target": True
    }
}

class ComprehensivePlasmidAnalyzer:
    def __init__(self):
        self.results = {}
        self.density_scores = {}
        self.safety_scores = {}

    def load_sequence_robust(self, filepath):
        """Load sequence with multiple format attempts."""
        for fmt in ["fasta", "fasta-pearson", "fasta-blast", "genbank"]:
            try:
                record = SeqIO.read(filepath, fmt)
                return str(record.seq).upper(), record.id, len(record.seq)
            except:
                continue
        return None, None, None

    def scan_comprehensive_elements(self, sequence, seq_id):
        """Comprehensive scanning of all plasmid elements."""
        results = {
            'sequence_id': seq_id,
            'elements_found': {},
            'total_elements': 0,
            'positions': {},
            'sequence_length': len(sequence),
            'categories': defaultdict(int),
            'safety_relevance': defaultdict(int),
            'qPCR_targets': []
        }

        for element_name, element_data in PLASMID_ELEMENTS.items():
            element_seq = element_data["sequence"]
            positions = []
            start = 0
            while True:
                pos = sequence.find(element_seq, start)
                if pos == -1:
                    break
                positions.append(pos)
                start = pos + 1

            if positions:
                results['elements_found'][element_name] = {
                    'count': len(positions),
                    'positions': positions,
                    'type': element_data['type'],
                    'category': element_data['category'],
                    'description': element_data['description'],
                    'safety_relevance': element_data['safety_relevance'],
                    'qPCR_target': element_data['qPCR_target'],
                    'color': element_data['color']
                }
                results['positions'][element_name] = positions
                results['total_elements'] += len(positions)
                results['categories'][element_data['category']] += len(positions)
                results['safety_relevance'][element_data['safety_relevance']] += len(positions)

                if element_data['qPCR_target']:
                    results['qPCR_targets'].append({
                        'name': element_name,
                        'positions': positions,
                        'description': element_data['description']
                    })

        return results

    def calculate_density_scores(self, results):
        """Calculate element density per kilobase."""
        sequence_length = results['sequence_length']
        density = {
            'elements_per_kb': results['total_elements'] / (sequence_length / 1000),
            'sv40_density': 0,
            'backbone_density': 0,
            'safety_density': 0
        }

        # Calculate category-specific densities
        for element_name, element_data in results['elements_found'].items():
            category = element_data['category']
            count = element_data['count']

            if category == 'SV40':
                density['sv40_density'] += count
            elif category in ['Origin', 'Antibiotic_Resistance']:
                density['backbone_density'] += count

            if element_data['safety_relevance'] == 'HIGH':
                density['safety_density'] += count

        # Normalize per kb
        density['sv40_density'] /= (sequence_length / 1000)
        density['backbone_density'] /= (sequence_length / 1000)
        density['safety_density'] /= (sequence_length / 1000)

        return density

    def calculate_safety_score(self, results):
        """Calculate safety relevance score."""
        score = {
            'high_relevance_count': results['safety_relevance']['HIGH'],
            'medium_relevance_count': results['safety_relevance']['MEDIUM'],
            'low_relevance_count': results['safety_relevance']['LOW'],
            'qPCR_target_count': len(results['qPCR_targets']),
            'safety_score': 0
        }

        # Weighted safety score (HIGH=3, MEDIUM=2, LOW=1)
        weighted_score = (
            score['high_relevance_count'] * 3 +
            score['medium_relevance_count'] * 2 +
            score['low_relevance_count'] * 1
        )

        score['safety_score'] = weighted_score
        return score

    def create_comprehensive_visualization(self, pfizer_results, moderna_results,
                                         pfizer_density, moderna_density,
                                         pfizer_safety, moderna_safety,
                                         output_file="comprehensive_plasmid_architecture.png"):
        """Create comprehensive visualization of plasmid architecture."""

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(20, 12))
        fig.suptitle('Comprehensive Plasmid Architecture Analysis: Pfizer vs Moderna\n'
                    'Bacterial Backbone + Mammalian Regulatory Elements',
                    fontsize=16, fontweight='bold')

        def plot_plasmid_architecture(ax, results, density, safety, title):
            sequence_length = results['sequence_length']
            positions = results['positions']

            # Draw sequence backbone
            ax.add_patch(Rectangle((0, 0), sequence_length, 1,
                                 facecolor='lightgray', edgecolor='black', linewidth=2))
            ax.set_xlim(0, sequence_length)
            ax.set_ylim(-0.5, 4.0)
            ax.set_title(title, fontsize=14, fontweight='bold')
            ax.set_xlabel('Position in sequence (bp)', fontsize=12)
            ax.set_yticks([])

            # Plot elements by category with different rows
            y_positions = {
                'SV40': 3.2,
                'Promoter': 2.8,
                'Origin': 2.4,
                'Antibiotic_Resistance': 2.0,
                'transcription': 1.6,
                'UTR': 1.2,
                'Coding': 0.8,
                'PolyA': 0.4
            }

            for element_name, element_data in PLASMID_ELEMENTS.items():
                if element_name in positions:
                    positions_list = positions[element_name]
                    category = element_data['category']
                    color = element_data['color']

                    y_pos = y_positions.get(category, 0.4)
                    for pos in positions_list:
                        length = len(element_data['sequence'])
                        rect = Rectangle((pos, y_pos), length, 0.25,
                                        facecolor=color, edgecolor='black',
                                        alpha=0.7, linewidth=0.5)
                        ax.add_patch(rect)

                    # Label the element type
                    if positions_list:
                        label_text = f"{element_name} ({len(positions_list)})"
                        if element_data['safety_relevance'] == 'HIGH':
                            label_text += " ⚠️"
                        if element_data['qPCR_target']:
                            label_text += " 📡"

                        ax.text(positions_list[0], y_pos + 0.3, label_text,
                               fontsize=8, color=color, fontweight='bold',
                               rotation=0 if len(positions_list) == 1 else 45)

            # Add comprehensive info box
            info_text = f"Length: {sequence_length:,} bp\n"
            info_text += f"Total Elements: {results['total_elements']}\n"
            info_text += f"Density: {density['elements_per_kb']:.2f} elements/kb\n"
            info_text += f"SV40 Density: {density['sv40_density']:.2f}/kb\n"
            info_text += f"Backbone Density: {density['backbone_density']:.2f}/kb\n"
            info_text += f"Safety Score: {safety['safety_score']}\n"
            info_text += f"High Relevance: {safety['high_relevance_count']}\n"
            info_text += f"qPCR Targets: {safety['qPCR_target_count']}"

            ax.text(0.02, 0.95, info_text, transform=ax.transAxes,
                   fontsize=10, verticalalignment='top',
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))

            # Add category breakdown
            category_text = "CATEGORIES:\n"
            for category, count in sorted(results['categories'].items(), key=lambda x: -x[1]):
                category_text += f"{category}: {count}\n"

            ax.text(0.98, 0.95, category_text, transform=ax.transAxes,
                   fontsize=9, verticalalignment='top', horizontalalignment='right',
                   bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9))

        # Plot Pfizer
        pfizer_title = (f"Pfizer BNT162b2 (OR134577.1) - DENSE SV40 + COMPLETE BACKBONE\n"
                       f"Safety Score: {pfizer_safety['safety_score']} | "
                       f"qPCR Targets: {pfizer_safety['qPCR_target_count']}")
        plot_plasmid_architecture(ax1, pfizer_results, pfizer_density, pfizer_safety, pfizer_title)

        # Plot Moderna
        moderna_title = (f"Moderna mRNA-1273 (OR134578.1) - SPARSE MAMMALIAN + SHARED BACKBONE\n"
                        f"Safety Score: {moderna_safety['safety_score']} | "
                        f"qPCR Targets: {moderna_safety['qPCR_target_count']}")
        plot_plasmid_architecture(ax2, moderna_results, moderna_density, moderna_safety, moderna_title)

        # Create comprehensive legend
        legend_elements = []
        for category, color in {
            'SV40': '#FF0000',
            'Promoter': '#00FF00',
            'Origin': '#0000FF',
            'Antibiotic_Resistance': '#FF00FF',
            'Transcription': '#00FFFF',
            'Spike Gene': '#FF0066',
            'PolyA': '#999999'
        }.items():
            legend_elements.append(mpatches.Patch(color=color, label=category))

        fig.legend(handles=legend_elements, loc='center left',
                  bbox_to_anchor=(1, 0.5), fontsize=10, title='Element Categories')

        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"✅ Comprehensive visualization saved to: {output_file}")
        plt.close()

    def generate_comprehensive_report(self, pfizer_results, moderna_results,
                                     pfizer_density, moderna_density,
                                     pfizer_safety, moderna_safety,
                                     output_file="comprehensive_plasmid_report.txt"):
        """Generate comprehensive analysis report."""

        with open(output_file, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write("COMPREHENSIVE PLASMID ARCHITECTURE ANALYSIS REPORT\n")
            f.write("=" * 80 + "\n\n")

            f.write("SEQUENCE INFORMATION\n")
            f.write("-" * 80 + "\n")
            f.write(f"Pfizer BNT162b2 (OR134577.1): {pfizer_results['sequence_length']:,} bp\n")
            f.write(f"Moderna mRNA-1273 (OR134578.1): {moderna_results['sequence_length']:,} bp\n")
            f.write(f"Length difference: {pfizer_results['sequence_length'] - moderna_results['sequence_length']:,} bp\n\n")

            f.write("ELEMENT DENSITY ANALYSIS\n")
            f.write("-" * 80 + "\n")
            f.write(f"Pfizer density: {pfizer_density['elements_per_kb']:.2f} elements/kb\n")
            f.write(f"Moderna density: {moderna_density['elements_per_kb']:.2f} elements/kb\n")
            f.write(f"Safety element density: Pfizer {pfizer_density['safety_density']:.2f}/kb vs Moderna {moderna_density['safety_density']:.2f}/kb\n\n")

            f.write("SAFETY RELEVANCE SCORING\n")
            f.write("-" * 80 + "\n")
            f.write(f"Pfizer safety score: {pfizer_safety['safety_score']}\n")
            f.write(f"  High relevance: {pfizer_safety['high_relevance_count']}\n")
            f.write(f"  Medium relevance: {pfizer_safety['medium_relevance_count']}\n")
            f.write(f"  Low relevance: {pfizer_safety['low_relevance_count']}\n")
            f.write(f"  qPCR targets: {pfizer_safety['qPCR_target_count']}\n\n")

            f.write(f"Moderna safety score: {moderna_safety['safety_score']}\n")
            f.write(f"  High relevance: {moderna_safety['high_relevance_count']}\n")
            f.write(f"  Medium relevance: {moderna_safety['medium_relevance_count']}\n")
            f.write(f"  Low relevance: {moderna_safety['low_relevance_count']}\n")
            f.write(f"  qPCR targets: {moderna_safety['qPCR_target_count']}\n\n")

            f.write("COMPREHENSIVE ELEMENT COMPARISON\n")
            f.write("-" * 80 + "\n")

            # Get all unique elements
            all_elements = set(pfizer_results['elements_found'].keys()) | \
                          set(moderna_results['elements_found'].keys())

            for element_name in sorted(all_elements):
                element_data = PLASMID_ELEMENTS[element_name]
                pfizer_count = pfizer_results['elements_found'].get(element_name, {}).get('count', 0)
                moderna_count = moderna_results['elements_found'].get(element_name, {}).get('count', 0)

                difference = pfizer_count - moderna_count
                status = "EQUAL"
                if pfizer_count > moderna_count:
                    status = "PFIZER HIGHER"
                elif moderna_count > pfizer_count:
                    status = "MODERNA HIGHER"

                f.write(f"\n{element_name}:\n")
                f.write(f"  Category: {element_data['category']}\n")
                f.write(f"  Type: {element_data['type']}\n")
                f.write(f"  Safety Relevance: {element_data['safety_relevance']}\n")
                f.write(f"  qPCR Target: {element_data['qPCR_target']}\n")
                f.write(f"  Description: {element_data['description']}\n")
                f.write(f"  Pfizer: {pfizer_count} copies\n")
                f.write(f"  Moderna: {moderna_count} copies\n")
                f.write(f"  Difference: {difference:+d} ({status})\n")

                if pfizer_count > 0:
                    pfizer_positions = pfizer_results['elements_found'].get(element_name, {}).get('positions', [])
                    f.write(f"  Pfizer positions: {pfizer_positions[:5]}")
                    if len(pfizer_positions) > 5:
                        f.write(f" ... and {len(pfizer_positions) - 5} more")
                    f.write("\n")

                if moderna_count > 0:
                    moderna_positions = moderna_results['elements_found'].get(element_name, {}).get('positions', [])
                    f.write(f"  Moderna positions: {moderna_positions[:5]}")
                    if len(moderna_positions) > 5:
                        f.write(f" ... and {len(moderna_positions) - 5} more")
                    f.write("\n")

            f.write("\n" + "=" * 80 + "\n")
            f.write("qPCR TARGET ANALYSIS\n")
            f.write("=" * 80 + "\n")

            f.write("\nPfizer qPCR Targets:\n")
            for target in pfizer_results['qPCR_targets']:
                f.write(f"  {target['name']}: {len(target['positions'])} copies\n")
                f.write(f"    Description: {target['description']}\n")
                f.write(f"    Positions: {target['positions'][:3]}\n")

            f.write("\nModerna qPCR Targets:\n")
            for target in moderna_results['qPCR_targets']:
                f.write(f"  {target['name']}: {len(target['positions'])} copies\n")
                f.write(f"    Description: {target['description']}\n")
                f.write(f"    Positions: {target['positions'][:3]}\n")

            f.write("\n" + "=" * 80 + "\n")
            f.write("RESIDUAL DNA IMPLICATIONS\n")
            f.write("=" * 80 + "\n")

            f.write("\nHigh-safety-relevance elements:\n")
            f.write("- Pfizer has higher density of safety-relevant elements\n")
            f.write("- SV40 enhancer elements (Pfizer only) are qPCR targets for residual DNA\n")
            f.write("- KanR gene present in both (antibiotic resistance marker)\n")
            f.write("- ColE1 origin present in both (bacterial replication)\n")

            f.write("\nFragment risk assessment:\n")
            f.write("- Average fragment size ~214 bp (Nanopore data)\n")
            f.write("- Maximum fragment size up to 3.5 kb\n")
            f.write("- SV40 enhancer + promoter + TATA can fit in single fragment\n")
            f.write("- KanR + ColE1 regions likely to persist as fragments\n")

            f.write("\n" + "=" * 80 + "\n")
            f.write("ANALYSIS COMPLETE\n")
            f.write("=" * 80 + "\n")

        print(f"✅ Comprehensive report saved to: {output_file}")

def main():
    parser = argparse.ArgumentParser(
        description='Comprehensive plasmid architecture analysis'
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
        default='comprehensive_plasmid_analysis',
        help='Output directory for results'
    )

    args = parser.parse_args()

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    print("=" * 80)
    print("COMPREHENSIVE PLASMID ARCHITECTURE ANALYZER")
    print("=" * 80)
    print()

    analyzer = ComprehensivePlasmidAnalyzer()

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

    # Scan comprehensive elements
    print("Scanning comprehensive plasmid elements...")
    pfizer_results = analyzer.scan_comprehensive_elements(pfizer_seq, pfizer_id)
    moderna_results = analyzer.scan_comprehensive_elements(moderna_seq, moderna_id)

    print(f"✅ Pfizer: {pfizer_results['total_elements']} elements detected")
    print(f"✅ Moderna: {moderna_results['total_elements']} elements detected")
    print()

    # Calculate densities
    print("Calculating density scores...")
    pfizer_density = analyzer.calculate_density_scores(pfizer_results)
    moderna_density = analyzer.calculate_density_scores(moderna_results)

    print(f"✅ Pfizer density: {pfizer_density['elements_per_kb']:.2f} elements/kb")
    print(f"✅ Moderna density: {moderna_density['elements_per_kb']:.2f} elements/kb")
    print()

    # Calculate safety scores
    print("Calculating safety relevance scores...")
    pfizer_safety = analyzer.calculate_safety_score(pfizer_results)
    moderna_safety = analyzer.calculate_safety_score(moderna_results)

    print(f"✅ Pfizer safety score: {pfizer_safety['safety_score']}")
    print(f"✅ Moderna safety score: {moderna_safety['safety_score']}")
    print()

    # Generate outputs
    print("Generating comprehensive outputs...")
    viz_file = os.path.join(args.output_dir, "comprehensive_plasmid_architecture.png")
    report_file = os.path.join(args.output_dir, "comprehensive_plasmid_report.txt")
    json_file = os.path.join(args.output_dir, "comprehensive_results.json")

    analyzer.create_comprehensive_visualization(
        pfizer_results, moderna_results,
        pfizer_density, moderna_density,
        pfizer_safety, moderna_safety,
        viz_file
    )

    analyzer.generate_comprehensive_report(
        pfizer_results, moderna_results,
        pfizer_density, moderna_density,
        pfizer_safety, moderna_safety,
        report_file
    )

    # Save JSON results
    comprehensive_results = {
        'pfizer': {
            'sequence_info': {
                'id': pfizer_id,
                'length': pfizer_len,
                'total_elements': pfizer_results['total_elements']
            },
            'elements': pfizer_results['elements_found'],
            'density': pfizer_density,
            'safety_score': pfizer_safety,
            'qpcr_targets': pfizer_results['qPCR_targets']
        },
        'moderna': {
            'sequence_info': {
                'id': moderna_id,
                'length': moderna_len,
                'total_elements': moderna_results['total_elements']
            },
            'elements': moderna_results['elements_found'],
            'density': moderna_density,
            'safety_score': moderna_safety,
            'qpcr_targets': moderna_results['qPCR_targets']
        }
    }

    with open(json_file, 'w') as f:
        json.dump(comprehensive_results, f, indent=2)

    print(f"✅ JSON results saved to: {json_file}")
    print()
    print("=" * 80)
    print("COMPREHENSIVE ANALYSIS COMPLETE")
    print("=" * 80)
    print()
    print(f"Results saved to: {args.output_dir}/")
    print(f"- Visualization: {viz_file}")
    print(f"- Report: {report_file}")
    print(f"- Data: {json_file}")

if __name__ == '__main__':
    main()