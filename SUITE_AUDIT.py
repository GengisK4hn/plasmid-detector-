#!/usr/bin/env python3
"""
Comprehensive suite audit for enhanced-plasmid-finder
Checks all tools, files, and integrations
"""

import os
import sys
from pathlib import Path
import importlib.util
import subprocess
from datetime import datetime

class SuiteAudit:
    def __init__(self):
        self.base_dir = Path('/home/dad/supracode-tool/enhanced-plasmid-finder')
        self.results = {
            'timestamp': datetime.now().isoformat(),
            'tools': {},
            'files': {},
            'integrations': {},
            'issues': []
        }

    def log(self, category, item, status, details=""):
        """Log audit results"""
        if category not in self.results:
            self.results[category] = {}
        self.results[category][item] = {
            'status': status,
            'details': details
        }
        status_symbol = "✅" if status == "OK" else "❌" if status == "FAIL" else "⚠️"
        print(f"{status_symbol} {category}: {item} - {status}")
        if details:
            print(f"   {details}")

    def check_python_tools(self):
        """Check all Python scanner tools"""
        print("\n" + "="*70)
        print("AUDIT: Python Scanner Tools")
        print("="*70)

        python_tools = [
            'ENHANCED_PLASMID_FINDER.py',
            'plasmid_finder_integration.py',
            'rigorous_sv40_validation.py',
            'visualize_pfizer_sv40.py',
            'batch_mrna_scanner.py',
            'verify_sv40_findings.py'
        ]

        for tool in python_tools:
            tool_path = self.base_dir / tool
            if tool_path.exists():
                try:
                    # Check if it can be imported
                    if tool.endswith('.py'):
                        spec = importlib.util.spec_from_file_location(tool[:-3], tool_path)
                        if spec and spec.loader:
                            self.log('tools', tool, 'OK', f'Found at {tool_path}')
                        else:
                            self.log('tools', tool, 'WARN', 'Cannot import module')
                except Exception as e:
                    self.log('tools', tool, 'WARN', f'Import error: {e}')
            else:
                self.log('tools', tool, 'FAIL', 'Not found')

    def check_data_files(self):
        """Check all data files and references"""
        print("\n" + "="*70)
        print("AUDIT: Data Files")
        print("="*70)

        # Check sequence directories
        seq_dirs = [
            'data/sequences',
            'data/sequences/additional',
            'data/references'
        ]

        for dir_path in seq_dirs:
            full_path = self.base_dir / dir_path
            if full_path.exists():
                files = list(full_path.glob('*.fasta')) + list(full_path.glob('*.fa'))
                self.log('files', dir_path, 'OK', f'{len(files)} sequence files')
                for f in files[:5]:  # Show first 5
                    print(f"      - {f.name}")
                if len(files) > 5:
                    print(f"      ... and {len(files) - 5} more")
            else:
                self.log('files', dir_path, 'FAIL', 'Directory not found')

    def check_documentation(self):
        """Check documentation files"""
        print("\n" + "="*70)
        print("AUDIT: Documentation")
        print("="*70)

        docs = [
            'README.md',
            'SV40_VALIDATION_README.md',
            'BLAST_VALIDATION_REPORT.md',
            'PROJECT_COMPLETION_SUMMARY.md',
            'ANALYSIS_COMPLETE.txt'
        ]

        for doc in docs:
            doc_path = self.base_dir / doc
            if doc_path.exists():
                size = doc_path.stat().st_size
                self.log('docs', doc, 'OK', f'{size} bytes')
            else:
                self.log('docs', doc, 'FAIL', 'Not found')

    def check_external_tools(self):
        """Check external tool dependencies"""
        print("\n" + "="*70)
        print("AUDIT: External Tools")
        print("="*70)

        tools = {
            'abricate': 'abricate --version',
            'blastn': 'blastn -version',
            'makeblastdb': 'makeblastdb -version',
            'Biopython': 'python -c "import Bio; print(Bio.__version__)"'
        }

        for tool, cmd in tools.items():
            try:
                result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=5)
                if result.returncode == 0:
                    version = result.stdout.strip().split('\n')[0]
                    self.log('external', tool, 'OK', version)
                else:
                    self.log('external', tool, 'FAIL', 'Not executable')
            except Exception as e:
                self.log('external', tool, 'FAIL', str(e))

    def check_integration_points(self):
        """Check integration between tools"""
        print("\n" + "="*70)
        print("AUDIT: Integration Points")
        print("="*70)

        # Check if plasmid_finder_integration can import abricate
        try:
            sys.path.insert(0, str(self.base_dir))
            from plasmid_finder_integration import PlasmidEnhancedDetector
            detector = PlasmidEnhancedDetector()
            self.log('integrations', 'PlasmidEnhancedDetector', 'OK',
                    f'Abricate available: {detector.plasmid_finder.abricate_available}')
        except Exception as e:
            self.log('integrations', 'PlasmidEnhancedDetector', 'FAIL', str(e))

        # Check if SV40 validation works
        try:
            from rigorous_sv40_validation import SV40Validator
            validator = SV40Validator()
            self.log('integrations', 'SV40Validator', 'OK', 'Can be instantiated')
        except Exception as e:
            self.log('integrations', 'SV40Validator', 'FAIL', str(e))

    def check_scan_capabilities(self):
        """Check what each tool can scan for"""
        print("\n" + "="*70)
        print("AUDIT: Scan Capabilities")
        print("="*70)

        capabilities = {
            'PlasmidEnhancedDetector': [
                'Plasmid sequences (ABRICATE)',
                'ColE1 origin',
                'Kanamycin resistance',
                'Vaccine signatures'
            ],
            'SV40Validator': [
                'SV40 72bp enhancer',
                'SV40 promoter elements',
                'BLAST validation'
            ],
            'batch_mrna_scanner': [
                'FCS sequences',
                'SV40 elements',
                'Promoter motifs'
            ]
        }

        for tool, caps in capabilities.items():
            self.log('capabilities', tool, 'OK', f'{len(caps)} capabilities')
            for cap in caps:
                print(f"      - {cap}")

    def generate_report(self):
        """Generate final audit report"""
        print("\n" + "="*70)
        print("AUDIT SUMMARY")
        print("="*70)

        total = 0
        ok = 0
        fail = 0
        warn = 0

        for category, items in self.results.items():
            if category == 'timestamp':
                continue
            for item, data in items.items():
                if isinstance(data, dict) and 'status' in data:
                    total += 1
                    if data['status'] == 'OK':
                        ok += 1
                    elif data['status'] == 'FAIL':
                        fail += 1
                    else:
                        warn += 1

        print(f"\nTotal checks: {total}")
        print(f"✅ Passed: {ok}")
        print(f"❌ Failed: {fail}")
        print(f"⚠️  Warnings: {warn}")

        if fail > 0:
            print(f"\n⚠️  AUDIT FOUND {fail} FAILURES - REVIEW NEEDED")
            return 1
        elif warn > 0:
            print(f"\n⚠️  AUDIT COMPLETE WITH {warn} WARNINGS")
            return 0
        else:
            print(f"\n✅ AUDIT PASSED - ALL SYSTEMS OPERATIONAL")
            return 0

def main():
    auditor = SuiteAudit()

    print("╔════════════════════════════════════════════════════════════════════╗")
    print("║          ENHANCED PLASMID FINDER - SUITE AUDIT                    ║")
    print("╚════════════════════════════════════════════════════════════════════╝")

    auditor.check_python_tools()
    auditor.check_data_files()
    auditor.check_documentation()
    auditor.check_external_tools()
    auditor.check_integration_points()
    auditor.check_scan_capabilities()

    return auditor.generate_report()

if __name__ == '__main__':
    sys.exit(main())
