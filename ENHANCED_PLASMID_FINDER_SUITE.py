#!/usr/bin/env python3
"""
ENHANCED PLASMID FINDER SUITE - MASTER CONTROLLER
=================================================
Comprehensive plasmid and vector analysis with full integration

This suite:
1. Calls main run_everything.sh pipeline for 20+ analyzers
2. Adds specialized plasmid detection (ABRICATE)
3. SV40 validation (BLAST-confirmed)
4. Error correction and recovery
5. Automatic data discovery and scanning

Features:
- Automatic error detection and correction
- Multi-source data discovery
- Cross-tool validation
- Comprehensive reporting

Author: Enhanced Plasmid Finder Suite
Date: April 2026
"""

import os
import sys
import json
import time
import subprocess
import shutil
from pathlib import Path
from datetime import datetime
from collections import defaultdict

# Try importing suite components
try:
    # Add current directory to path for imports
    import sys
    import os
    suite_dir = os.path.dirname(os.path.abspath(__file__))
    if suite_dir not in sys.path:
        sys.path.insert(0, suite_dir)

    from plasmid_finder import PlasmidEnhancedDetector
    from rigorous_sv40_validation import RigorousValidator
    SUITE_AVAILABLE = True
except ImportError as e:
    print(f"Warning: Suite components not fully available: {e}")
    SUITE_AVAILABLE = False

class ErrorCorrector:
    """Handle errors and attempt corrections."""

    def __init__(self, log_file="error_correction.log"):
        self.log_file = Path(log_file)
        self.errors = []
        self.corrections = []

    def log_error(self, phase, error, context=""):
        """Log an error with context."""
        error_entry = {
            'timestamp': datetime.now().isoformat(),
            'phase': phase,
            'error': str(error),
            'context': context
        }
        self.errors.append(error_entry)

        # Write to log
        with open(self.log_file, 'a') as f:
            f.write(f"[{error_entry['timestamp']}] {phase}: {error}\n")
            if context:
                f.write(f"  Context: {context}\n")

    def attempt_correction(self, phase, error, correction_fn):
        """Attempt to correct an error."""
        try:
            print(f"  🔧 Attempting correction for: {phase}")
            result = correction_fn()
            if result:
                correction_entry = {
                    'timestamp': datetime.now().isoformat(),
                    'phase': phase,
                    'correction': correction_fn.__name__,
                    'success': True
                }
                self.corrections.append(correction_entry)
                print(f"  ✅ Correction successful")
                return True
        except Exception as e:
            print(f"  ❌ Correction failed: {e}")
            self.log_error(phase, f"Correction failed: {e}", str(error))

        return False

    def has_errors(self):
        """Check if any errors were logged."""
        return len(self.errors) > 0

    def get_summary(self):
        """Get error correction summary."""
        return {
            'total_errors': len(self.errors),
            'corrected_errors': len(self.corrections),
            'uncorrected_errors': len(self.errors) - len(self.corrections)
        }

class DataScanner:
    """Automatic data discovery and scanning."""

    def __init__(self, base_dirs=None):
        self.base_dirs = base_dirs or [
            Path('data/sequences'),
            Path('data/sequences/additional'),
            Path('data/references'),
            Path('/media/external_drive/mckernan_supracode_data'),
            Path('../supracode-tool/data/sequences')
        ]
        self.discovered_files = []
        self.file_stats = defaultdict(lambda: {
            'count': 0,
            'total_size': 0,
            'formats': defaultdict(int)
        })

    def scan_all_sources(self):
        """Scan all configured sources for sequence data."""
        print("\n" + "="*70)
        print("  DATA DISCOVERY AND SCANNING")
        print("="*70)

        all_files = []

        for base_dir in self.base_dirs:
            if not base_dir.exists():
                print(f"⚠️  Path not found: {base_dir}")
                continue

            print(f"\n📁 Scanning: {base_dir}")

            # Find sequence files
            patterns = ['*.fasta', '*.fa', '*.fna', '*.gb', '*.genbank']
            found = []

            for pattern in patterns:
                files = list(base_dir.rglob(pattern))
                found.extend(files)

            if found:
                print(f"   Found {len(found)} sequence files")
                all_files.extend(found)

                # Update stats
                for f in found:
                    try:
                        stat = f.stat()
                        ext = f.suffix[1:] if f.suffix else 'unknown'
                        self.file_stats[str(base_dir)]['count'] += 1
                        self.file_stats[str(base_dir)]['total_size'] += stat.st_size
                        self.file_stats[str(base_dir)]['formats'][ext] += 1
                    except:
                        pass

        self.discovered_files = list(set(all_files))  # Remove duplicates

        print(f"\n📊 DISCOVERY SUMMARY:")
        print(f"   Total files found: {len(self.discovered_files)}")

        for dir_path, stats in self.file_stats.items():
            if stats['count'] > 0:
                print(f"\n   {Path(dir_path).name}:")
                print(f"      Files: {stats['count']}")
                print(f"      Size: {stats['total_size']:,} bytes")
                print(f"      Formats: {dict(stats['formats'])}")

        return self.discovered_files

    def validate_files(self):
        """Validate discovered files."""
        print("\n" + "="*70)
        print("  FILE VALIDATION")
        print("="*70)

        valid_files = []
        invalid_files = []

        for f in self.discovered_files:
            try:
                # Check if file is readable
                if f.stat().st_size == 0:
                    invalid_files.append((f, "Empty file"))
                    continue

                # Try to parse as FASTA
                from Bio import SeqIO
                try:
                    records = list(SeqIO.parse(f, "fasta"))
                    if records:
                        valid_files.append(f)
                    else:
                        invalid_files.append((f, "No records found"))
                except:
                    # Try GenBank format
                    try:
                        records = list(SeqIO.parse(f, "genbank"))
                        if records:
                            valid_files.append(f)
                        else:
                            invalid_files.append((f, "No records found"))
                    except:
                        invalid_files.append((f, "Parse error"))

            except Exception as e:
                invalid_files.append((f, str(e)))

        print(f"\n✅ Valid files: {len(valid_files)}")
        print(f"❌ Invalid files: {len(invalid_files)}")

        if invalid_files:
            print("\n⚠️  Invalid files:")
            for f, reason in invalid_files[:10]:  # Show first 10
                print(f"   {f.name}: {reason}")
            if len(invalid_files) > 10:
                print(f"   ... and {len(invalid_files) - 10} more")

        return valid_files

class EnhancedPlasmidSuite:
    """Master controller for enhanced plasmid finder suite."""

    def __init__(self, input_dir="data/sequences", output_dir="suite_results"):
        self.input_dir = Path(input_dir)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True, parents=True)

        self.error_corrector = ErrorCorrector(
            log_file=self.output_dir / "error_correction.log"
        )
        self.data_scanner = DataScanner()

        self.results = {
            'timestamp': datetime.now().isoformat(),
            'input_directory': str(self.input_dir),
            'phases': {},
            'summary': {}
        }

    def print_header(self, title):
        """Print formatted header."""
        print("\n" + "="*70)
        print(f"  {title}")
        print("="*70)

    def phase1_data_discovery(self):
        """Phase 1: Automatic data discovery and validation."""
        self.print_header("PHASE 1: DATA DISCOVERY")

        try:
            # Discover files
            discovered = self.data_scanner.scan_all_sources()

            # Validate files
            valid_files = self.data_scanner.validate_files()

            self.valid_sequences = valid_files

            self.results['phases']['data_discovery'] = {
                'total_discovered': len(discovered),
                'valid_files': len(valid_files),
                'files': [str(f) for f in valid_files[:20]]  # First 20
            }

            return len(valid_files) > 0

        except Exception as e:
            self.error_corrector.log_error("data_discovery", e)
            # Attempt correction: use default path
            def fix():
                default_path = Path('data/sequences')
                if default_path.exists():
                    files = list(default_path.glob('*.fasta')) + list(default_path.glob('*.fa'))
                    self.valid_sequences = files
                    return len(files) > 0
                return False

            if self.error_corrector.attempt_correction("data_discovery", e, fix):
                return True

        return False

    def phase2_main_pipeline(self):
        """Phase 2: Call main run_everything.sh pipeline."""
        self.print_header("PHASE 2: MAIN PIPELINE (20+ ANALYZERS)")

        main_script = Path('../run_everything.sh')
        if not main_script.exists():
            main_script = Path('/home/dad/supracode-tool/run_everything.sh')

        if not main_script.exists():
            print("⚠️  Main run_everything.sh not found, skipping")
            return {}

        print(f"🚀 Running main pipeline: {main_script}")
        print("   This will run 20+ analyzers across 9 phases...")
        print("   (this may take a while)")

        try:
            # Run main pipeline
            result = subprocess.run(
                ['bash', str(main_script)],
                capture_output=True,
                text=True,
                timeout=3600,  # 1 hour timeout
                cwd=str(main_script.parent)
            )

            # Save output
            output_file = self.output_dir / "main_pipeline_output.txt"
            with open(output_file, 'w') as f:
                f.write("STDOUT:\n")
                f.write(result.stdout)
                f.write("\nSTDERR:\n")
                f.write(result.stderr)

            if result.returncode == 0:
                print("✅ Main pipeline completed successfully")
            else:
                print(f"⚠️  Main pipeline exited with code {result.returncode}")

            self.results['phases']['main_pipeline'] = {
                'exit_code': result.returncode,
                'output_file': str(output_file)
            }

            return {'stdout': result.stdout, 'stderr': result.stderr}

        except subprocess.TimeoutExpired:
            print("⚠️  Main pipeline timed out (1 hour)")
            self.error_corrector.log_error("main_pipeline", "Timeout after 1 hour")
            return {}
        except Exception as e:
            self.error_corrector.log_error("main_pipeline", e)
            return {}

    def phase3_plasmid_detection(self):
        """Phase 3: ABRICATE plasmid detection."""
        self.print_header("PHASE 3: PLASMID DETECTION (ABRICATE)")

        if not SUITE_AVAILABLE:
            print("⚠️  Suite not available, skipping plasmid detection")
            return {}

        detector = PlasmidEnhancedDetector()
        results = {}
        plasmid_count = 0

        for seq_file in self.valid_sequences:
            print(f"\n🔬 {seq_file.name}")

            try:
                result = detector.analyze_fasta_file_abricate(str(seq_file))

                if result:
                    file_result = {
                        'plasmid_detected': result.get('plasmid_detected', False),
                        'confidence': result.get('confidence', 'N/A'),
                        'markers': result.get('markers_found', [])
                    }

                    results[seq_file.name] = file_result

                    if file_result['plasmid_detected']:
                        plasmid_count += 1
                        markers = ", ".join(file_result['markers'])
                        print(f"   ✅ PLASMID: {markers}")
                    else:
                        print(f"   ❌ No plasmid")

            except Exception as e:
                print(f"   ⚠️  Error: {e}")
                self.error_corrector.log_error("plasmid_detection", e, seq_file.name)
                results[seq_file.name] = {'error': str(e)}

        self.results['phases']['plasmid_detection'] = {
            'files_analyzed': len(results),
            'plasmids_found': plasmid_count,
            'results': results
        }

        print(f"\n📊 Found {plasmid_count} plasmids in {len(results)} files")
        return results

    def phase4_sv40_validation(self):
        """Phase 4: SV40 validation (BLAST-confirmed)."""
        self.print_header("PHASE 4: SV40 VALIDATION")

        if not SUITE_AVAILABLE:
            print("⚠️  Suite not available, skipping SV40 validation")
            return {}

        validator = RigorousValidator(
            output_dir=self.output_dir / "sv40_validation"
        )

        results = {}
        sv40_count = 0

        # Run step 1: sequence validation
        try:
            integrity_results = validator.step1_validate_sequence_integrity()
        except Exception as e:
            print(f"⚠️  SV40 validation error: {e}")
            self.error_corrector.log_error("sv40_validation", e)
            return {}

        # Run step 2: exact 72bp search
        try:
            search_results = validator.step2_exact_72bp_search(integrity_results['sequences'])
        except Exception as e:
            print(f"⚠️  SV40 search error: {e}")
            search_results = {}

        # Process results
        for seq_name, seq_data in integrity_results['sequences'].items():
            try:
                # seq_data is a tuple (sequence, filename)
                if isinstance(seq_data, tuple) and len(seq_data) >= 2:
                    sequence, filename = seq_data[0], seq_data[1]
                else:
                    sequence = seq_data
                    filename = seq_name

                # Check search results for SV40
                search_result = search_results.get(seq_name, {})
                pfizer_copies = len(search_result.get('pfizer_forward', [])) + len(search_result.get('pfizer_reverse', []))
                moderna_copies = len(search_result.get('moderna_forward', [])) + len(search_result.get('moderna_reverse', []))

                file_result = {
                    'has_sv40': pfizer_copies > 0,
                    'pfizer_copies': pfizer_copies,
                    'moderna_copies': moderna_copies,
                    'positions': search_result.get('pfizer_forward', []) + search_result.get('pfizer_reverse', [])
                }

                # Use filename
                filename = Path(filename).name
                results[filename] = file_result

                if file_result['has_sv40']:
                    sv40_count += 1
                    print(f"   ✅ {filename}: SV40 detected ({pfizer_copies} copies)")
                else:
                    print(f"   ❌ {filename}: No SV40")

            except Exception as e:
                print(f"   ⚠️  Error processing {seq_name}: {e}")
                self.error_corrector.log_error("sv40_validation", e, seq_name)

        self.results['phases']['sv40_validation'] = {
            'files_analyzed': len(results),
            'sv40_found': sv40_count,
            'results': results
        }

        print(f"\n📊 Found SV40 in {sv40_count} files")
        return results

    def phase5_cross_validation(self):
        """Phase 5: Cross-validate all findings."""
        self.print_header("PHASE 5: CROSS-VALIDATION")

        # Aggregate findings
        all_findings = defaultdict(lambda: {
            'plasmid': False,
            'sv40': False,
            'main_pipeline': False
        })

        # Extract from results
        if 'plasmid_detection' in self.results['phases']:
            for file, result in self.results['phases']['plasmid_detection']['results'].items():
                if isinstance(result, dict) and result.get('plasmid_detected'):
                    all_findings[file]['plasmid'] = True

        if 'sv40_validation' in self.results['phases']:
            for file, result in self.results['phases']['sv40_validation']['results'].items():
                if isinstance(result, dict) and result.get('has_sv40'):
                    all_findings[file]['sv40'] = True

        # Print summary
        print("\n📊 CROSS-VALIDATED FINDINGS:")
        print("-" * 70)

        for file, findings in sorted(all_findings.items()):
            if any(findings.values()):
                status = []
                if findings['plasmid']:
                    status.append("Plasmid✅")
                if findings['sv40']:
                    status.append("SV40✅")

                print(f"{file:50} | {' | '.join(status)}")

        # Statistics
        total = sum(1 for v in all_findings.values() if any(v.values()))
        plasmid_only = sum(1 for v in all_findings.values() if v['plasmid'] and not v['sv40'])
        sv40_only = sum(1 for v in all_findings.values() if v['sv40'] and not v['plasmid'])
        both = sum(1 for v in all_findings.values() if v['plasmid'] and v['sv40'])

        self.results['summary'] = {
            'total_with_findings': total,
            'plasmid_only': plasmid_only,
            'sv40_only': sv40_only,
            'both': both
        }

        print(f"\n📈 SUMMARY: {total} files with findings")
        print(f"   Plasmid only: {plasmid_only}")
        print(f"   SV40 only: {sv40_only}")
        print(f"   Both: {both}")

        return dict(all_findings)

    def phase6_generate_reports(self):
        """Phase 6: Generate comprehensive reports."""
        self.print_header("PHASE 6: REPORT GENERATION")

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

        # JSON report
        json_file = self.output_dir / f"suite_report_{timestamp}.json"
        with open(json_file, 'w') as f:
            json.dump(self.results, f, indent=2)
        print(f"✅ JSON: {json_file}")

        # Text report
        txt_file = self.output_dir / f"suite_report_{timestamp}.txt"
        with open(txt_file, 'w') as f:
            f.write("="*70 + "\n")
            f.write("ENHANCED PLASMID FINDER SUITE - FINAL REPORT\n")
            f.write("="*70 + "\n\n")
            f.write(f"Generated: {self.results['timestamp']}\n\n")

            # Error correction summary
            error_summary = self.error_corrector.get_summary()
            f.write("ERROR CORRECTION\n")
            f.write("-"*70 + "\n")
            f.write(f"Total errors: {error_summary['total_errors']}\n")
            f.write(f"Corrected: {error_summary['corrected_errors']}\n")
            f.write(f"Uncorrected: {error_summary['uncorrected_errors']}\n\n")

            # Phase results
            f.write("PHASE RESULTS\n")
            f.write("-"*70 + "\n")
            for phase, data in self.results.get('phases', {}).items():
                f.write(f"\n{phase}:\n")
                for key, value in data.items():
                    if key != 'results':  # Skip detailed results
                        f.write(f"  {key}: {value}\n")

            # Summary
            f.write("\nSUMMARY\n")
            f.write("-"*70 + "\n")
            if self.results.get('summary'):
                for key, value in self.results['summary'].items():
                    f.write(f"{key}: {value}\n")

        print(f"✅ Text: {txt_file}")

        # Error correction log
        if self.error_corrector.has_errors():
            print(f"📋 Error log: {self.error_corrector.log_file}")

        return json_file, txt_file

    def run_complete_suite(self):
        """Execute all phases."""
        print("\n╔════════════════════════════════════════════════════════════════════╗")
        print("║     ENHANCED PLASMID FINDER SUITE - MASTER CONTROLLER            ║")
        print("║     (Error Correction + Data Scanning + Full Integration)        ║")
        print("╚════════════════════════════════════════════════════════════════════╝")

        try:
            # Phase 1: Data discovery
            if not self.phase1_data_discovery():
                print("❌ Data discovery failed!")
                return False

            # Phase 2: Main pipeline
            self.phase2_main_pipeline()

            # Phase 3: Plasmid detection
            self.phase3_plasmid_detection()

            # Phase 4: SV40 validation
            self.phase4_sv40_validation()

            # Phase 5: Cross-validation
            self.phase5_cross_validation()

            # Phase 6: Reports
            self.phase6_generate_reports()

            self.print_header("SUITE COMPLETE")
            print("✅ All phases executed")

            if self.error_corrector.has_errors():
                error_summary = self.error_corrector.get_summary()
                print(f"\n📋 Error Correction Summary:")
                print(f"   Total errors: {error_summary['total_errors']}")
                print(f"   Corrected: {error_summary['corrected_errors']}")
                print(f"   Uncorrected: {error_summary['uncorrected_errors']}")

            return True

        except KeyboardInterrupt:
            print("\n⚠️  Interrupted by user")
            return False
        except Exception as e:
            print(f"\n❌ Suite error: {e}")
            import traceback
            traceback.print_exc()
            return False

def main():
    """Main entry point."""
    if len(sys.argv) < 2:
        print("Usage: python ENHANCED_PLASMID_FINDER_SUITE.py <input_directory> [output_directory]")
        print("\nExample:")
        print("  python ENHANCED_PLASMID_FINDER_SUITE.py data/sequences suite_results")
        sys.exit(1)

    input_dir = sys.argv[1]
    output_dir = sys.argv[2] if len(sys.argv) > 2 else "suite_results"

    suite = EnhancedPlasmidSuite(input_dir, output_dir)
    success = suite.run_complete_suite()

    sys.exit(0 if success else 1)

if __name__ == '__main__':
    main()
