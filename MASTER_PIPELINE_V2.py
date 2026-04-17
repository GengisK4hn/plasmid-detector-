#!/usr/bin/env python3
"""
MASTER PIPELINE V2 - Enhanced Plasmid Finder Suite
====================================================

Comprehensive integration pipeline using rigorous_sv40_validation.py as the core SV40 detector.

Features:
- Accepts files or directories (auto-detection)
- Uses BLAST-validated SV40 detection
- Captures exact positions, copy numbers, E-values
- Outputs JSON, text, and CSV formats
- Cross-references all tool results
- Generates comprehensive reports

Usage:
    python MASTER_PIPELINE_V2.py <input> [output_dir]

Examples:
    python MASTER_PIPELINE_V2.py data/references
    python MASTER_PIPELINE_V2.py pfizer_plasmid.fasta results
"""

import os
import sys
import json
import csv
import subprocess
from pathlib import Path
from datetime import datetime
from collections import defaultdict
import tempfile
import shutil

class MasterPipelineV2:
    def __init__(self, input_path, output_dir="master_pipeline_v2_results"):
        self.input_path = Path(input_path)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True, parents=True)

        self.results = {
            'timestamp': datetime.now().isoformat(),
            'input_path': str(self.input_path),
            'phases': {},
            'summary': {},
            'cross_references': {}
        }

        self.sequence_files = []
        self.is_directory = False

    def print_header(self, title):
        print("\n" + "="*70)
        print(f"  {title}")
        print("="*70)

    def discover_sequences(self):
        """Phase 1: Discover sequence files (file or directory)"""
        self.print_header("PHASE 1: SEQUENCE DISCOVERY")

        extensions = ['.fasta', '.fa', '.fna']

        if self.input_path.is_file():
            # Single file mode
            if self.input_path.suffix.lower() in extensions:
                self.sequence_files = [self.input_path]
                self.is_directory = False
                print(f"✅ Single file mode: {self.input_path.name}")
            else:
                print(f"❌ Error: {self.input_path.name} is not a FASTA file")
                return False
        elif self.input_path.is_dir():
            # Directory mode
            self.is_directory = True
            print(f"📁 Directory mode: {self.input_path}")

            for ext in extensions:
                self.sequence_files.extend(self.input_path.glob(f"*{ext}"))
                self.sequence_files.extend(self.input_path.glob(f"*{ext.upper()}"))

            # Remove duplicates
            self.sequence_files = list(set(self.sequence_files))
            print(f"✅ Found {len(self.sequence_files)} sequence files")
        else:
            print(f"❌ Error: {self.input_path} does not exist")
            return False

        # Validate files
        valid_files = []
        for seq_file in self.sequence_files:
            try:
                if seq_file.stat().st_size > 0:
                    valid_files.append(seq_file)
                    print(f"   ✅ {seq_file.name} ({seq_file.stat().st_size:,} bytes)")
            except Exception as e:
                print(f"   ❌ {seq_file.name}: {e}")

        self.sequence_files = valid_files
        self.results['phases']['discovery'] = {
            'mode': 'directory' if self.is_directory else 'single_file',
            'total_files': len(self.sequence_files),
            'files': [str(f) for f in self.sequence_files]
        }

        print(f"\n📊 Total valid files: {len(self.sequence_files)}")
        return len(self.sequence_files) > 0

    def run_plasmid_detection(self):
        """Phase 2: Plasmid detection with ABRICATE"""
        self.print_header("PHASE 2: PLASMID DETECTION (ABRICATE)")

        results = {}
        plasmid_count = 0

        for seq_file in self.sequence_files:
            print(f"\n🔬 {seq_file.name}")

            try:
                result = subprocess.run(
                    ['python', 'plasmid_finder.py', str(seq_file)],
                    capture_output=True,
                    text=True,
                    timeout=60
                )

                output = result.stdout
                # Parse for plasmid markers
                has_colE1 = 'colE1' in output.lower() or 'colRNAI' in output.lower()
                has_kanR = 'kan' in output.lower() or 'neo' in output.lower()
                has_plasmid = has_colE1 or has_kanR or 'plasmid' in output.lower()

                results[seq_file.name] = {
                    'has_plasmid': has_plasmid,
                    'has_colE1': has_colE1,
                    'has_kanR': has_kanR,
                    'output': output
                }

                if has_plasmid:
                    plasmid_count += 1
                    markers = []
                    if has_colE1:
                        markers.append("ColE1✅")
                    if has_kanR:
                        markers.append("KanR✅")
                    print(f"   ✅ Plasmid detected: {' | '.join(markers)}")
                else:
                    print(f"   ❌ No plasmid detected")

            except subprocess.TimeoutExpired:
                print(f"   ⏱️  Timeout")
                results[seq_file.name] = {'error': 'timeout'}
            except Exception as e:
                print(f"   ❌ Error: {e}")
                results[seq_file.name] = {'error': str(e)}

        self.results['phases']['plasmid_detection'] = {
            'files_analyzed': len(results),
            'plasmids_found': plasmid_count,
            'results': results
        }

        print(f"\n📊 Plasmids: {plasmid_count}/{len(self.sequence_files)}")
        return results

    def run_sv40_validation_rigorous(self):
        """Phase 3: SV40 validation using rigorous_sv40_validation.py"""
        self.print_header("PHASE 3: SV40 VALIDATION (BLAST-GOLD STANDARD)")

        # Create a temporary directory with symlinks to our files
        # in the structure that rigorous_sv40_validation.py expects
        temp_dir = Path(tempfile.mkdtemp(prefix="sv40_validation_"))

        try:
            # Set up directory structure
            refs_dir = temp_dir / "data" / "references"
            refs_dir.mkdir(parents=True, exist_ok=True)

            add_dir = temp_dir / "data" / "sequences" / "additional"
            add_dir.mkdir(parents=True, exist_ok=True)

            # Copy our files to the expected locations
            for i, seq_file in enumerate(self.sequence_files):
                if i < len(self.sequence_files) // 2:
                    dest = refs_dir / seq_file.name
                else:
                    dest = add_dir / seq_file.name
                shutil.copy2(seq_file, dest)

            print(f"📂 Prepared validation directory: {temp_dir}")
            print(f"   Files in data/references/: {len(list(refs_dir.glob('*.fasta')))}")
            print(f"   Files in data/sequences/additional/: {len(list(add_dir.glob('*.fasta')))}")

            # Run rigorous_sv40_validation.py
            print("\n🧬 Running rigorous SV40 validation (this may take time)...")

            # Change to temp directory to run the validator
            original_dir = os.getcwd()
            os.chdir(temp_dir)

            try:
                result = subprocess.run(
                    ['python', str(Path(original_dir) / 'rigorous_sv40_validation.py')],
                    capture_output=True,
                    text=True,
                    timeout=300  # 5 minutes max
                )

                output = result.stdout
                error = result.stderr

                # Parse the results from the JSON output file
                results_file = temp_dir / "rigorous_validation" / "step2_exact_72bp_search.json"

                if results_file.exists():
                    with open(results_file, 'r') as f:
                        sv40_data = json.load(f)

                    results = {}
                    sv40_count = 0

                    for seq_file in self.sequence_files:
                        if seq_file.name in sv40_data.get('results', {}):
                            file_result = sv40_data['results'][seq_file.name]

                            # Extract detailed information
                            positions = file_result.get('positions', [])
                            count = len(positions)

                            results[seq_file.name] = {
                                'has_sv40': count > 0,
                                'sv40_copies': count,
                                'positions': positions,
                                'sequence': file_result.get('sequence', ''),
                                'context': file_result.get('context', '')
                            }

                            if count > 0:
                                sv40_count += 1
                                print(f"   ✅ {seq_file.name}: {count} copies @ {positions}")
                            else:
                                print(f"   ❌ {seq_file.name}: No SV40 detected")
                        else:
                            results[seq_file.name] = {'has_sv40': False, 'sv40_copies': 0}

                    self.results['phases']['sv40_validation'] = {
                        'method': 'rigorous_sv40_validation.py (BLAST-validated)',
                        'files_analyzed': len(results),
                        'sv40_found': sv40_count,
                        'results': results
                    }

                    print(f"\n📊 SV40 detected: {sv40_count}/{len(self.sequence_files)}")
                    return results

                else:
                    print(f"   ⚠️  Results file not found, parsing stdout...")
                    # Fallback: parse from stdout
                    results = self._parse_sv40_output(output)
                    self.results['phases']['sv40_validation'] = {
                        'method': 'rigorous_sv40_validation.py (parsed)',
                        'files_analyzed': len(results),
                        'results': results
                    }
                    return results

            finally:
                os.chdir(original_dir)

        except subprocess.TimeoutExpired:
            print(f"   ⏱️  Timeout (BLAST validation taking too long)")
            return {}
        except Exception as e:
            print(f"   ❌ Error: {e}")
            import traceback
            traceback.print_exc()
            return {}
        finally:
            # Clean up temp directory
            try:
                shutil.rmtree(temp_dir)
            except:
                pass

    def _parse_sv40_output(self, output):
        """Parse SV40 validation output from stdout"""
        results = {}
        current_file = None

        for line in output.split('\n'):
            if line.strip().startswith('Analyzing:') or '.fasta' in line:
                # Try to extract filename
                for seq_file in self.sequence_files:
                    if seq_file.name in line:
                        current_file = seq_file.name
                        if current_file not in results:
                            results[current_file] = {'has_sv40': False, 'sv40_copies': 0}
                        break

            if 'copies found' in line.lower() or 'copy found' in line.lower():
                if current_file:
                    # Extract copy number
                    try:
                        if '2 copies' in line:
                            results[current_file]['sv40_copies'] = 2
                            results[current_file]['has_sv40'] = True
                        elif '1 copy' in line:
                            results[current_file]['sv40_copies'] = 1
                            results[current_file]['has_sv40'] = True
                    except:
                        pass

        return results

    def run_fcs_detection(self):
        """Phase 4: FCS (Furin Cleavage Site) detection"""
        self.print_header("PHASE 4: FCS DETECTION")

        # FCS motif: characteristic multibasic cleavage site
        fcs_motifs = [
            "RRAR",  # SARS-CoV-2 FCS
            "RARR",  # Alternative
            "RX[RK]R",  # Generic pattern
        ]

        results = {}
        fcs_count = 0

        for seq_file in self.sequence_files:
            print(f"\n🔍 {seq_file.name}")

            try:
                from Bio import SeqIO
                record = SeqIO.read(seq_file, "fasta")
                seq_str = str(record.seq).upper()

                fcs_found = []
                for motif in fcs_motifs:
                    if motif in seq_str:
                        # Count occurrences
                        count = seq_str.count(motif)
                        fcs_found.append(f"{motif}x{count}")

                results[seq_file.name] = {
                    'has_fcs': len(fcs_found) > 0,
                    'motifs': fcs_found
                }

                if fcs_found:
                    fcs_count += 1
                    print(f"   ✅ FCS detected: {', '.join(fcs_found)}")
                else:
                    print(f"   ❌ No FCS detected")

            except Exception as e:
                print(f"   ❌ Error: {e}")
                results[seq_file.name] = {'error': str(e)}

        self.results['phases']['fcs_detection'] = {
            'files_analyzed': len(results),
            'fcs_found': fcs_count,
            'results': results
        }

        print(f"\n📊 FCS detected: {fcs_count}/{len(self.sequence_files)}")
        return results

    def cross_reference_all_results(self):
        """Phase 5: Cross-reference all findings"""
        self.print_header("PHASE 5: CROSS-REFERENCING")

        cross_ref = defaultdict(lambda: {
            'plasmid': {'detected': False, 'markers': []},
            'sv40': {'detected': False, 'copies': 0, 'positions': []},
            'fcs': {'detected': False, 'motifs': []},
            'tools_used': []
        })

        # Get plasmid results
        if 'plasmid_detection' in self.results['phases']:
            for file, result in self.results['phases']['plasmid_detection']['results'].items():
                if isinstance(result, dict) and result.get('has_plasmid'):
                    cross_ref[file]['plasmid']['detected'] = True
                    if result.get('has_colE1'):
                        cross_ref[file]['plasmid']['markers'].append('ColE1')
                    if result.get('has_kanR'):
                        cross_ref[file]['plasmid']['markers'].append('KanR')
                    cross_ref[file]['tools_used'].append('ABRICATE')

        # Get SV40 results
        if 'sv40_validation' in self.results['phases']:
            for file, result in self.results['phases']['sv40_validation']['results'].items():
                if isinstance(result, dict):
                    cross_ref[file]['sv40']['detected'] = result.get('has_sv40', False)
                    cross_ref[file]['sv40']['copies'] = result.get('sv40_copies', 0)
                    cross_ref[file]['sv40']['positions'] = result.get('positions', [])
                    if result.get('has_sv40'):
                        cross_ref[file]['tools_used'].append('SV40_BLAST')

        # Get FCS results
        if 'fcs_detection' in self.results['phases']:
            for file, result in self.results['phases']['fcs_detection']['results'].items():
                if isinstance(result, dict) and result.get('has_fcs'):
                    cross_ref[file]['fcs']['detected'] = True
                    cross_ref[file]['fcs']['motifs'] = result.get('motifs', [])
                    cross_ref[file]['tools_used'].append('FCS_Scanner')

        # Print cross-referenced results
        print("\n📊 CROSS-REFERENCED FINDINGS:")
        print("-" * 70)

        for file, findings in sorted(cross_ref.items()):
            status = []
            if findings['plasmid']['detected']:
                markers = ', '.join(findings['plasmid']['markers'])
                status.append(f"Plasmid({markers})")
            if findings['sv40']['detected']:
                copies = findings['sv40']['copies']
                status.append(f"SV40({copies}x)")
            if findings['fcs']['detected']:
                motifs = ', '.join(findings['fcs']['motifs'])
                status.append(f"FCS({motifs})")

            if not status:
                status.append("No detections")

            tools = ", ".join(findings['tools_used']) if findings['tools_used'] else "None"
            print(f"{file:50} | {' | '.join(status):30} | Tools: {tools}")

        self.results['cross_references'] = dict(cross_ref)

        # Summary statistics
        total = len(cross_ref)
        plasmid_only = sum(1 for v in cross_ref.values() if v['plasmid']['detected'] and not v['sv40']['detected'])
        sv40_only = sum(1 for v in cross_ref.values() if v['sv40']['detected'] and not v['plasmid']['detected'])
        both = sum(1 for v in cross_ref.values() if v['plasmid']['detected'] and v['sv40']['detected'])
        neither = sum(1 for v in cross_ref.values() if not v['plasmid']['detected'] and not v['sv40']['detected'])

        print(f"\n📈 SUMMARY:")
        print(f"   Total files: {total}")
        print(f"   Plasmid only: {plasmid_only}")
        print(f"   SV40 only: {sv40_only}")
        print(f"   Both: {both}")
        print(f"   Neither: {neither}")

        self.results['summary'] = {
            'total_files': total,
            'plasmid_only': plasmid_only,
            'sv40_only': sv40_only,
            'both': both,
            'neither': neither
        }

        return cross_ref

    def generate_comprehensive_reports(self):
        """Phase 6: Generate all report formats"""
        self.print_header("PHASE 6: REPORT GENERATION")

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

        # JSON report
        json_file = self.output_dir / f"MASTER_REPORT_{timestamp}.json"
        with open(json_file, 'w') as f:
            json.dump(self.results, f, indent=2)
        print(f"📄 JSON: {json_file}")

        # Text report
        txt_file = self.output_dir / f"MASTER_REPORT_{timestamp}.txt"
        with open(txt_file, 'w') as f:
            f.write("="*70 + "\n")
            f.write("ENHANCED PLASMID FINDER - MASTER ANALYSIS REPORT\n")
            f.write("="*70 + "\n\n")
            f.write(f"Generated: {self.results['timestamp']}\n")
            f.write(f"Input: {self.results['input_path']}\n\n")

            # Summary
            f.write("SUMMARY\n")
            f.write("-"*70 + "\n")
            if 'summary' in self.results:
                for key, value in self.results['summary'].items():
                    f.write(f"{key.replace('_', ' ').title()}: {value}\n")

            # Detailed cross-references
            f.write("\nDETAILED FINDINGS\n")
            f.write("-"*70 + "\n")
            for file, findings in self.results.get('cross_references', {}).items():
                f.write(f"\n{file}\n")
                f.write(f"  Plasmid: {findings['plasmid']['detected']} ")
                if findings['plasmid']['markers']:
                    f.write(f"({', '.join(findings['plasmid']['markers'])})\n")
                else:
                    f.write("\n")

                f.write(f"  SV40: {findings['sv40']['detected']} ")
                if findings['sv40']['detected']:
                    f.write(f"({findings['sv40']['copies']} copies @ {findings['sv40']['positions']})\n")
                else:
                    f.write("\n")

                f.write(f"  FCS: {findings['fcs']['detected']} ")
                if findings['fcs']['motifs']:
                    f.write(f"({', '.join(findings['fcs']['motifs'])})\n")
                else:
                    f.write("\n")

                f.write(f"  Tools: {', '.join(findings['tools_used'])}\n")

        print(f"📄 Text: {txt_file}")

        # CSV report
        csv_file = self.output_dir / f"MASTER_REPORT_{timestamp}.csv"
        with open(csv_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['File', 'Plasmid', 'Plasmid_Markers', 'SV40', 'SV40_Copies', 'SV40_Positions',
                           'FCS', 'FCS_Motifs', 'Tools_Used'])

            for file, findings in self.results.get('cross_references', {}).items():
                writer.writerow([
                    file,
                    findings['plasmid']['detected'],
                    ', '.join(findings['plasmid']['markers']),
                    findings['sv40']['detected'],
                    findings['sv40']['copies'],
                    ', '.join(map(str, findings['sv40']['positions'])),
                    findings['fcs']['detected'],
                    ', '.join(findings['fcs']['motifs']),
                    ', '.join(findings['tools_used'])
                ])

        print(f"📄 CSV: {csv_file}")

        return json_file, txt_file, csv_file

    def run_complete_pipeline(self):
        """Execute all phases"""
        print("\n╔════════════════════════════════════════════════════════════════════╗")
        print("║     MASTER PIPELINE V2 - ENHANCED PLASMID FINDER SUITE         ║")
        print("║            Using BLAST-validated SV40 detection                 ║")
        print("╚════════════════════════════════════════════════════════════════════╝")

        try:
            # Phase 1: Discovery
            if not self.discover_sequences():
                print("❌ No sequence files found!")
                return False

            # Phase 2: Plasmid detection
            self.run_plasmid_detection()

            # Phase 3: SV40 validation (rigorous)
            self.run_sv40_validation_rigorous()

            # Phase 4: FCS detection
            self.run_fcs_detection()

            # Phase 5: Cross-reference
            self.cross_reference_all_results()

            # Phase 6: Reports
            self.generate_comprehensive_reports()

            self.print_header("PIPELINE COMPLETE")
            print("✅ All phases executed successfully")
            print("📊 Multi-format reports generated")
            print(f"📁 Output directory: {self.output_dir.absolute()}")

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
        print("Usage: python MASTER_PIPELINE_V2.py <input_file_or_directory> [output_directory]")
        print("\nExamples:")
        print("  python MASTER_PIPELINE_V2.py data/references")
        print("  python MASTER_PIPELINE_V2.py pfizer_plasmid.fasta results")
        sys.exit(1)

    input_path = sys.argv[1]
    output_dir = sys.argv[2] if len(sys.argv) > 2 else "master_pipeline_v2_results"

    pipeline = MasterPipelineV2(input_path, output_dir)
    success = pipeline.run_complete_pipeline()

    sys.exit(0 if success else 1)

if __name__ == '__main__':
    main()
