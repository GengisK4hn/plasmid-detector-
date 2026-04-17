#!/usr/bin/env python3
"""
Comprehensive Batch Plasmid Analysis

Analyzes ALL sequences for plasmid content including:
- Vaccine sequences (Pfizer, Moderna)
- SARS-CoV-2 natural sequences (negative controls)
- SV40 positive controls
- Additional viral sequences
- Bacterial controls

Author: Enhanced Plasmid Finder Suite
Date: April 2026
"""

import sys
sys.path.insert(0, '.')

from pathlib import Path
from Bio import SeqIO
from plasmid_finder_integration import PlasmidEnhancedDetector
import json
from datetime import datetime

def analyze_all_sequences():
    """Analyze all FASTA sequences for plasmid content."""

    print("=" * 80)
    print("COMPREHENSIVE BATCH PLASMID ANALYSIS")
    print("=" * 80)
    print(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print()

    # Initialize detector
    detector = PlasmidEnhancedDetector()
    print(f"Abricate available: {detector.plasmid_finder.abricate_available}")
    print()

    # Find all FASTA files
    seq_files = (list(Path('data/sequences').glob("*.fasta")) +
                list(Path('data/sequences').glob("*.fa")) +
                list(Path('data/sequences/additional').glob("*.fasta")) +
                list(Path('data/sequences/additional').glob("*.fa")))

    print(f"Total files found: {len(seq_files)}")
    print()

    # Categorize sequences
    categories = {
        'vaccines': [],
        'sars_cov2': [],
        'sv40_controls': [],
        'viral_other': [],
        'bacterial': [],
        'other': []
    }

    # Categorization logic
    for seq_file in seq_files:
        name = seq_file.name.lower()

        if any(x in name for x in ['pfizer', 'moderna', 'bnt', 'mrna']):
            categories['vaccines'].append(seq_file)
        elif any(x in name for x in ['nc_045512', 'mt020880', 'sars_cov', 'wuhan', 'delta', 'omicron', 'nmadoh']):
            categories['sars_cov2'].append(seq_file)
        elif any(x in name for x in ['nc_001669', 'sv40']):
            categories['sv40_controls'].append(seq_file)
        elif any(x in name for x in ['staphylococcus', 'bacterial', 'coli']):
            categories['bacterial'].append(seq_file)
        elif any(x in name for x in ['influenza', 'hiv', 'mers', 'prion']):
            categories['viral_other'].append(seq_file)
        else:
            categories['other'].append(seq_file)

    # Analyze each category
    all_results = {}
    plasmids_found = 0

    for category, files in categories.items():
        if not files:
            continue

        print("=" * 80)
        print(f"CATEGORY: {category.upper()}")
        print("=" * 80)

        for seq_file in files:
            try:
                print(f"\n[{len(all_results)+1}] {seq_file.name}")
                print("-" * 80)

                result = detector.analyze_fasta_file_abricate(str(seq_file))

                if result:
                    # Store result
                    all_results[seq_file.name] = {
                        'category': category,
                        'plasmid_detected': result.get('plasmid_detected', False),
                        'method': result.get('method', 'unknown'),
                        'confidence': result.get('confidence', 'NONE'),
                        'hit_count': result.get('hit_count', 0),
                        'markers_found': result.get('markers_found', []),
                        'hits': result.get('hits', [])
                    }

                    # Display result
                    if result.get('plasmid_detected'):
                        plasmids_found += 1
                        print(f"  ✅ PLASMID DETECTED!")
                        print(f"  Confidence: {result.get('confidence', 'N/A')}")
                        print(f"  Hit Count: {result.get('hit_count', 0)}")
                        print(f"  Markers: {', '.join(result.get('markers_found', [])[:5])}")
                        if len(result.get('markers_found', [])) > 5:
                            print(f"    ... and {len(result.get('markers_found', [])) - 5} more")
                    else:
                        print(f"  ❌ No plasmid detected")
                        print(f"  Method: {result.get('method', 'unknown')}")
                else:
                    print(f"  ⚠️  Analysis failed")
                    all_results[seq_file.name] = {
                        'category': category,
                        'plasmid_detected': False,
                        'error': 'Analysis failed'
                    }

            except Exception as e:
                print(f"  ❌ Error: {e}")
                all_results[seq_file.name] = {
                    'category': category,
                    'plasmid_detected': False,
                    'error': str(e)
                }

    # Summary
    print()
    print("=" * 80)
    print("ANALYSIS SUMMARY")
    print("=" * 80)
    print(f"Total sequences analyzed: {len(all_results)}")
    print(f"Plasmids detected: {plasmids_found}")
    print()

    # Category breakdown
    print("BY CATEGORY:")
    print("-" * 80)
    for category, files in categories.items():
        if not files:
            continue
        category_plasmids = sum(1 for f in files
                               if all_results.get(f.name, {}).get('plasmid_detected', False))
        print(f"{category:20s}: {category_plasmids}/{len(files)} plasmids detected")

    print()
    print("=" * 80)
    print("PLASMID-POSITIVE SEQUENCES:")
    print("=" * 80)

    plasmid_positive = {name: data for name, data in all_results.items()
                       if data.get('plasmid_detected', False)}

    if plasmid_positive:
        for name, data in sorted(plasmid_positive.items(),
                                 key=lambda x: x[1].get('hit_count', 0),
                                 reverse=True):
            print(f"\n{name}")
            print(f"  Category: {data.get('category', 'unknown')}")
            print(f"  Confidence: {data.get('confidence', 'N/A')}")
            print(f"  Hit Count: {data.get('hit_count', 0)}")
            print(f"  Markers: {', '.join(data.get('markers_found', [])[:3])}")
    else:
        print("No plasmids detected in any sequences")

    print()
    print("=" * 80)
    print("EXPECTED RESULTS VALIDATION:")
    print("=" * 80)

    # Validation checks
    checks = []

    # Check 1: SARS-CoV-2 sequences should have NO plasmids
    sars_results = {name: data for name, data in all_results.items()
                   if data.get('category') == 'sars_cov2'}
    sars_plasmids = sum(1 for data in sars_results.values()
                       if data.get('plasmid_detected', False))
    checks.append(("SARS-CoV-2 negative control",
                  sars_plasmids == 0,
                  f"{sars_plasmids} plasmids found (expected 0)"))

    # Check 2: SV40 should have plasmids
    sv40_results = {name: data for name, data in all_results.items()
                   if data.get('category') == 'sv40_controls'}
    sv40_plasmids = sum(1 for data in sv40_results.values()
                       if data.get('plasmid_detected', False))
    checks.append(("SV40 positive control",
                  sv40_plasmids > 0,
                  f"{sv40_plasmids} plasmids found (expected >0)"))

    # Check 3: Bacterial should have plasmids
    bacterial_results = {name: data for name, data in all_results.items()
                        if data.get('category') == 'bacterial'}
    bacterial_plasmids = sum(1 for data in bacterial_results.values()
                            if data.get('plasmid_detected', False))
    checks.append(("Bacterial positive control",
                  bacterial_plasmids > 0,
                  f"{bacterial_plasmids} plasmids found (expected >0)"))

    # Check 4: Vaccines may have plasmids
    vaccine_results = {name: data for name, data in all_results.items()
                      if data.get('category') == 'vaccines'}
    vaccine_plasmids = sum(1 for data in vaccine_results.values()
                          if data.get('plasmid_detected', False))
    checks.append(("mRNA vaccines",
                  vaccine_plasmids >= 0,  # May or may not have plasmids
                  f"{vaccine_plasmids} plasmids found (variable)"))

    for check_name, passed, message in checks:
        status = "✅ PASS" if passed else "❌ FAIL"
        print(f"{status}: {check_name:30s} - {message}")

    # Save results
    output_file = f"plasmid_analysis_results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    with open(output_file, 'w') as f:
        json.dump({
            'timestamp': datetime.now().isoformat(),
            'total_sequences': len(all_results),
            'plasmids_detected': plasmids_found,
            'results': all_results,
            'validation_checks': [
                {'name': name, 'passed': passed, 'message': message}
                for name, passed, message in checks
            ]
        }, f, indent=2)

    print()
    print(f"Results saved to: {output_file}")
    print()
    print("=" * 80)
    print("ANALYSIS COMPLETE")
    print("=" * 80)

    return all_results


if __name__ == "__main__":
    analyze_all_sequences()
