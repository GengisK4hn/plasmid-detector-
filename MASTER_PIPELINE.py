#!/usr/bin/env python3
"""
Master Integrated Pipeline for Enhanced Plasmid Finder Suite
Orchestrates all tools for comprehensive analysis with cross-referencing
"""

import os
import sys
import json
import subprocess
from pathlib import Path
from datetime import datetime
from collections import defaultdict

class MasterPipeline:
    def __init__(self, input_dir, output_dir):
        self.input_dir = Path(input_dir)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.results = {
            'timestamp': datetime.now().isoformat(),
            'input_directory': str(self.input_dir),
            'phases': {},
            'summary': {},
            'cross_references': {}
        }
        self.sequence_files = []

    def print_header(self, title):
        print("\n" + "="*70)
        print(f"  {title}")
        print("="*70)

    def discover_sequences(self):
        """Phase 1: Discover and validate all sequence files"""
        self.print_header("PHASE 1: SEQUENCE DISCOVERY")

        extensions = ['.fasta', '.fa', '.fna']
        self.sequence_files = []

        for ext in extensions:
            self.sequence_files.extend(self.input_dir.glob(f"**/*{ext}"))
            self.sequence_files.extend(self.input_dir.glob(f"**/*{ext.upper()}"))

        # Remove duplicates
        self.sequence_files = list(set(self.sequence_files))

        valid_files = []
        for seq_file in self.sequence_files:
            try:
                # Quick validation - check if file is readable and has content
                if seq_file.stat().st_size > 0:
                    valid_files.append(seq_file)
                    print(f"✅ Found: {seq_file.name} ({seq_file.stat().st_size} bytes)")
            except Exception as e:
                print(f"❌ Error: {seq_file.name} - {e}")

        self.sequence_files = valid_files
        self.results['phases']['discovery'] = {
            'total_files': len(self.sequence_files),
            'files': [str(f) for f in self.sequence_files]
        }

        print(f"\n📊 Total valid sequence files: {len(self.sequence_files)}")
        return len(self.sequence_files) > 0

    def run_plasmid_detection(self):
        """Phase 2: Plasmid detection with ABRICATE"""
        self.print_header("PHASE 2: PLASMID DETECTION (ABRICATE)")

        results = {}
        plasmid_count = 0

        for seq_file in self.sequence_files:
            print(f"\n🔬 Analyzing: {seq_file.name}")

            try:
                # Use plasmid_finder.py
                result = subprocess.run(
                    ['python', 'plasmid_finder.py', str(seq_file)],
                    capture_output=True,
                    text=True,
                    timeout=60
                )

                if result.returncode == 0:
                    # Parse output for plasmid detection
                    output = result.stdout
                    has_plasmid = 'plasmid' in output.lower() or 'col' in output.lower()

                    results[seq_file.name] = {
                        'has_plasmid': has_plasmid,
                        'output': output,
                        'exit_code': result.returncode
                    }

                    if has_plasmid:
                        plasmid_count += 1
                        print(f"   ✅ Plasmid detected")
                    else:
                        print(f"   ❌ No plasmid detected")
                else:
                    print(f"   ⚠️  Error running analysis")
                    results[seq_file.name] = {'error': result.stderr}

            except subprocess.TimeoutExpired:
                print(f"   ⏱️  Timeout - file too large")
                results[seq_file.name] = {'error': 'timeout'}
            except Exception as e:
                print(f"   ❌ Error: {e}")
                results[seq_file.name] = {'error': str(e)}

        self.results['phases']['plasmid_detection'] = {
            'files_analyzed': len(results),
            'plasmids_found': plasmid_count,
            'results': results
        }

        print(f"\n📊 Plasmids detected: {plasmid_count}/{len(self.sequence_files)}")
        return results

    def run_sv40_validation(self):
        """Phase 3: SV40 enhancer validation"""
        self.print_header("PHASE 3: SV40 ENHANCER VALIDATION")

        results = {}
        sv40_count = 0

        for seq_file in self.sequence_files:
            print(f"\n🧬 Analyzing: {seq_file.name}")

            try:
                # Use batch_mrna_scanner.py for SV40 detection
                result = subprocess.run(
                    ['python', 'batch_mrna_scanner.py', str(seq_file)],
                    capture_output=True,
                    text=True,
                    timeout=60
                )

                output = result.stdout
                # Parse for SV40 detection - look for specific patterns
                has_sv40 = False
                sv40_copies = 0

                # Check for SV40 enhancer patterns
                if 'SV40 72bp enhancer' in output:
                    has_sv40 = True
                    # Try to extract copy number
                    if '2 copies' in output or 'copies: 2' in output:
                        sv40_copies = 2
                    elif '1 copy' in output or 'copies: 1' in output:
                        sv40_copies = 1
                    else:
                        # Count occurrences
                        sv40_copies = output.count('SV40')

                results[seq_file.name] = {
                    'has_sv40': has_sv40,
                    'sv40_copies': sv40_copies,
                    'output': output
                }

                if has_sv40:
                    sv40_count += 1
                    if sv40_copies > 0:
                        print(f"   ✅ SV40 detected: {sv40_copies} copies")
                    else:
                        print(f"   ✅ SV40 detected")
                else:
                    print(f"   ❌ No SV40 detected")

            except subprocess.TimeoutExpired:
                print(f"   ⏱️  Timeout")
                results[seq_file.name] = {'error': 'timeout'}
            except Exception as e:
                print(f"   ❌ Error: {e}")
                results[seq_file.name] = {'error': str(e)}

        self.results['phases']['sv40_validation'] = {
            'files_analyzed': len(results),
            'sv40_found': sv40_count,
            'results': results
        }

        print(f"\n📊 SV40 detected: {sv40_count}/{len(self.sequence_files)}")
        return results

    def run_promoter_analysis(self):
        """Phase 4: Promoter element analysis"""
        self.print_header("PHASE 4: PROMOTER ELEMENT ANALYSIS")

        results = {}

        for seq_file in self.sequence_files:
            print(f"\n🎯 Analyzing: {seq_file.name}")

            try:
                # Use enhanced_sv40_analyzer.py
                result = subprocess.run(
                    ['python', 'enhanced_sv40_analyzer.py', str(seq_file)],
                    capture_output=True,
                    text=True,
                    timeout=60
                )

                output = result.stdout
                results[seq_file.name] = {
                    'output': output
                }

                # Check for promoter elements
                if 'promoter' in output.lower():
                    print(f"   ✅ Promoter elements detected")
                else:
                    print(f"   ℹ️  Analysis complete")

            except Exception as e:
                print(f"   ⚠️  Skipped: {e}")
                results[seq_file.name] = {'error': str(e)}

        self.results['phases']['promoter_analysis'] = {
            'files_analyzed': len(results),
            'results': results
        }

        return results

    def cross_reference_results(self):
        """Phase 5: Cross-reference all findings"""
        self.print_header("PHASE 5: CROSS-REFERENCING RESULTS")

        cross_ref = defaultdict(lambda: {
            'plasmid': False,
            'sv40': False,
            'promoters': False,
            'tools_used': []
        })

        # Get plasmid results
        if 'plasmid_detection' in self.results['phases']:
            for file, result in self.results['phases']['plasmid_detection']['results'].items():
                if isinstance(result, dict) and result.get('has_plasmid'):
                    cross_ref[file]['plasmid'] = True
                    cross_ref[file]['tools_used'].append('ABRICATE')

        # Get SV40 results
        if 'sv40_validation' in self.results['phases']:
            for file, result in self.results['phases']['sv40_validation']['results'].items():
                if isinstance(result, dict) and result.get('has_sv40'):
                    cross_ref[file]['sv40'] = True
                    cross_ref[file]['tools_used'].append('SV40_Validator')

        # Get promoter results
        if 'promoter_analysis' in self.results['phases']:
            for file, result in self.results['phases']['promoter_analysis']['results'].items():
                if isinstance(result, dict) and 'promoter' in str(result.get('output', '')).lower():
                    cross_ref[file]['promoters'] = True
                    cross_ref[file]['tools_used'].append('Promoter_Analyzer')

        # Print cross-referenced results
        print("\n📊 CROSS-REFERENCED FINDINGS:")
        print("-" * 70)

        for file, findings in sorted(cross_ref.items()):
            status = []
            if findings['plasmid']:
                status.append("Plasmid✅")
            if findings['sv40']:
                status.append("SV40✅")
            if findings['promoters']:
                status.append("Promoters✅")

            tools = ", ".join(findings['tools_used']) if findings['tools_used'] else "None"
            print(f"{file:50} | {' | '.join(status):20} | Tools: {tools}")

        self.results['cross_references'] = dict(cross_ref)

        # Summary statistics
        total = len(cross_ref)
        plasmid_only = sum(1 for v in cross_ref.values() if v['plasmid'] and not v['sv40'])
        sv40_only = sum(1 for v in cross_ref.values() if v['sv40'] and not v['plasmid'])
        both = sum(1 for v in cross_ref.values() if v['plasmid'] and v['sv40'])

        print(f"\n📈 SUMMARY:")
        print(f"   Total files analyzed: {total}")
        print(f"   Plasmid only: {plasmid_only}")
        print(f"   SV40 only: {sv40_only}")
        print(f"   Both: {both}")

        self.results['summary'] = {
            'total_files': total,
            'plasmid_only': plasmid_only,
            'sv40_only': sv40_only,
            'both': both
        }

        return cross_ref

    def generate_master_report(self):
        """Phase 6: Generate comprehensive master report"""
        self.print_header("PHASE 6: MASTER REPORT GENERATION")

        # Create report filename
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        report_file = self.output_dir / f"MASTER_REPORT_{timestamp}.json"
        txt_report = self.output_dir / f"MASTER_REPORT_{timestamp}.txt"

        # Save JSON report
        with open(report_file, 'w') as f:
            json.dump(self.results, f, indent=2)

        # Create text report
        with open(txt_report, 'w') as f:
            f.write("="*70 + "\n")
            f.write("ENHANCED PLASMID FINDER - MASTER ANALYSIS REPORT\n")
            f.write("="*70 + "\n\n")
            f.write(f"Generated: {self.results['timestamp']}\n")
            f.write(f"Input Directory: {self.results['input_directory']}\n\n")

            # Summary
            f.write("SUMMARY\n")
            f.write("-"*70 + "\n")
            if 'summary' in self.results:
                for key, value in self.results['summary'].items():
                    f.write(f"{key.replace('_', ' ').title()}: {value}\n")

            # Cross-referenced findings
            f.write("\nCROSS-REFERENCED FINDINGS\n")
            f.write("-"*70 + "\n")
            for file, findings in self.results.get('cross_references', {}).items():
                f.write(f"\n{file}\n")
                f.write(f"  Plasmid: {'YES' if findings['plasmid'] else 'NO'}\n")
                f.write(f"  SV40: {'YES' if findings['sv40'] else 'NO'}\n")
                f.write(f"  Promoters: {'YES' if findings['promoters'] else 'NO'}\n")
                f.write(f"  Tools: {', '.join(findings['tools_used'])}\n")

        print(f"\n📄 Reports generated:")
        print(f"   JSON: {report_file}")
        print(f"   Text: {txt_report}")

        return report_file, txt_report

    def run_complete_pipeline(self):
        """Execute all phases of the master pipeline"""
        print("\n╔════════════════════════════════════════════════════════════════════╗")
        print("║     MASTER PIPELINE - ENHANCED PLASMID FINDER SUITE             ║")
        print("╚════════════════════════════════════════════════════════════════════╝")

        try:
            # Phase 1: Discovery
            if not self.discover_sequences():
                print("❌ No sequence files found!")
                return False

            # Phase 2: Plasmid detection
            self.run_plasmid_detection()

            # Phase 3: SV40 validation
            self.run_sv40_validation()

            # Phase 4: Promoter analysis
            self.run_promoter_analysis()

            # Phase 5: Cross-reference
            self.cross_reference_results()

            # Phase 6: Generate report
            self.generate_master_report()

            self.print_header("PIPELINE COMPLETE")
            print("✅ All phases executed successfully")
            print("📊 Master report generated with cross-referenced findings")

            return True

        except KeyboardInterrupt:
            print("\n\n⚠️  Pipeline interrupted by user")
            return False
        except Exception as e:
            print(f"\n\n❌ Pipeline error: {e}")
            import traceback
            traceback.print_exc()
            return False

def main():
    if len(sys.argv) < 2:
        print("Usage: python MASTER_PIPELINE.py <input_directory> [output_directory]")
        print("\nExample:")
        print("  python MASTER_PIPELINE.py data/sequences results")
        sys.exit(1)

    input_dir = sys.argv[1]
    output_dir = sys.argv[2] if len(sys.argv) > 2 else "master_pipeline_results"

    pipeline = MasterPipeline(input_dir, output_dir)
    success = pipeline.run_complete_pipeline()

    sys.exit(0 if success else 1)

if __name__ == '__main__':
    main()
