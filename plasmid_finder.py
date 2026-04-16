#!/usr/bin/env python3
"""
Enhanced PlasmidFinder Integration
Detects plasmid sequences, replication systems, and mobility genes

Features:
- Direct Abricate integration with detailed hit information
- Position, coverage, and identity tracking for each marker
- Enhanced error handling and timeouts
- Fallback to built-in detection when tools unavailable
- Rich reporting with confidence scoring
"""

import os
import logging
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import re

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class PlasmidFinderIntegration:
    """
    Integration with PlasmidFinder for plasmid detection

    PlasmidFinder detects:
    - Plasmid replication systems
    - Plasmid mobility genes
    - Plasmid typing
    """

    # Plasmid replication system markers (simplified database)
    PLASMID_MARKERS = {
        # Inc groups (incompatibility groups)
        'IncF': ['traJ', 'finc', 'pilX'],
        'IncI': ['shufflon', 'pilV'],
        'IncA/C': ['repA', 'parA', 'parB'],
        'IncN': 'korA',
        'IncP': 'trfA',
        'IncQ': 'repB',
        'IncW': 'repA',
        # Rolling circle
        'RepL': 'rep',
        'RepC': 'rep',
        # Conjugative
        'traI': 'traI',
        'traG': 'traG',
        'mob': 'mob',
        # Resistance plasmids
        'bla_TEM': 'bla',
        'CTX-M': 'bla',
        'NDM': 'bla',
    }

    def __init__(self, conda_env: Optional[str] = None, plasmidfinder_db_path: Optional[str] = None):
        """
        Initialize PlasmidFinder integration

        Args:
            conda_env: Name of conda environment (default: None, uses current env)
            plasmidfinder_db_path: Path to PlasmidFinder database
        """
        self.conda_env = conda_env
        self.plasmidfinder_available = self._check_plasmidfinder()
        self.abricate_available = self._check_abricate()
        self.db_path = plasmidfinder_db_path
        self.markers_found = []

        logger.info(f"PlasmidFinder available: {self.plasmidfinder_available}")
        logger.info(f"Abricate available: {self.abricate_available}")

    def _check_plasmidfinder(self) -> bool:
        """Check if PlasmidFinder is installed"""
        try:
            result = subprocess.run(['plasmidfinder', '--help'],
                                  capture_output=True,
                                  timeout=5)
            return result.returncode == 0
        except (FileNotFoundError, subprocess.TimeoutExpired):
            return False

    def _check_abricate(self) -> bool:
        """Check if Abricate with plasmidfinder database is available"""
        try:
            # First try direct abricate command (works if running in conda env)
            result = subprocess.run(
                ['abricate', '--version'],
                capture_output=True,
                timeout=5
            )
            if result.returncode == 0:
                logger.info("Abricate found directly (in conda environment)")
                return True

            # Fallback: try conda run if environment specified
            if self.conda_env:
                result = subprocess.run(
                    ['conda', 'run', '-n', self.conda_env, 'abricate', '--version'],
                    capture_output=True,
                    timeout=5
                )
                if result.returncode == 0:
                    logger.info(f"Abricate found via conda env '{self.conda_env}'")
                    return True

            logger.warning("Abricate not found")
            return False
        except (FileNotFoundError, subprocess.TimeoutExpired) as e:
            logger.warning(f"Abricate check failed: {e}")
            return False

    def install_plasmidfinder(self) -> bool:
        """
        Install PlasmidFinder via conda

        Returns:
            True if successful
        """
        try:
            logger.info("Installing PlasmidFinder...")
            subprocess.run([
                'conda', 'install', '-y',
                '-c', 'bioconda',
                '-c', 'conda-forge',
                'plasmidfinder'
            ], check=True, capture_output=True)
            logger.info("PlasmidFinder installed successfully")
            return True
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to install PlasmidFinder: {e}")
            return False

    def detect_plasmid_markers_builtin(self, sequence: Seq) -> Dict[str, List[int]]:
        """
        Detect plasmid markers using built-in simplified database

        This is a fallback when PlasmidFinder is not available

        Args:
            sequence: Input sequence

        Returns:
            Dictionary of marker names and their positions
        """
        seq_str = str(sequence).upper()
        markers_found = {}

        for marker_name, patterns in self.PLASMID_MARKERS.items():
            positions = []

            if isinstance(patterns, str):
                patterns = [patterns]

            for pattern in patterns:
                # Search for marker (allow some mismatches)
                for match in re.finditer(f'(?=({pattern}))', seq_str):
                    positions.append(match.start())

            if positions:
                markers_found[marker_name] = positions

        return markers_found

    def detect_plasmid_markers_plasmidfinder(self, fasta_file: str) -> Dict:
        """
        Run PlasmidFinder on FASTA file

        Args:
            fasta_file: Path to input FASTA file

        Returns:
            Dictionary with plasmid detection results
        """
        if not self.plasmidfinder_available:
            logger.warning("PlasmidFinder not available, using built-in detection")
            return {}

        try:
            # Run PlasmidFinder
            cmd = [
                'plasmidfinder',
                '-i', fasta_file,
                '-o', 'plasmidfinder_results',
                '-p', 'plasmidfinder_db'  # Path to database
            ]

            result = subprocess.run(cmd,
                                    capture_output=True,
                                    text=True,
                                    timeout=300)

            if result.returncode == 0:
                # Parse results
                return self._parse_plasmidfinder_results('plasmidfinder_results')
            else:
                logger.warning(f"PlasmidFinder failed: {result.stderr}")
                return {}

        except subprocess.TimeoutExpired:
            logger.error("PlasmidFinder timed out")
            return {}
        except Exception as e:
            logger.error(f"PlasmidFinder error: {e}")
            return {}

    def _parse_plasmidfinder_results(self, result_dir: str) -> Dict:
        """Parse PlasmidFinder output files"""
        results = {}

        # Parse fsa file (summary)
        fsa_file = Path(result_dir) / 'results.fsa'
        if fsa_file.exists():
            with open(fsa_file) as f:
                for line in f:
                    if line.startswith('>'):
                        parts = line[1:].strip().split()
                        if len(parts) > 1:
                            seq_id = parts[0]
                            plasmid_type = parts[1]
                            results[seq_id] = plasmid_type

        return results

    def detect_plasmid_markers_abricate(self, fasta_file: str) -> Dict:
        """
        Run Abricate with plasmidfinder database on FASTA file

        Args:
            fasta_file: Path to input FASTA file

        Returns:
            Dictionary with plasmid detection results including detailed hit information
        """
        if not self.abricate_available:
            logger.warning("Abricate not available")
            return {}

        try:
            # Try direct abricate command first (if in conda env)
            cmd = ['abricate', '--db', 'plasmidfinder', fasta_file]

            result = subprocess.run(cmd,
                                    capture_output=True,
                                    text=True,
                                    timeout=300)

            # If direct command failed and conda env specified, try conda run
            if result.returncode != 0 and self.conda_env:
                cmd = ['conda', 'run', '-n', self.conda_env, 'abricate', '--db', 'plasmidfinder', fasta_file]
                result = subprocess.run(cmd,
                                        capture_output=True,
                                        text=True,
                                        timeout=300)

            if result.returncode == 0:
                return self._parse_abricate_results(result.stdout)
            else:
                logger.warning(f"Abricate failed: {result.stderr}")
                return {}

        except subprocess.TimeoutExpired:
            logger.error("Abricate timed out")
            return {}
        except Exception as e:
            logger.error(f"Abricate error: {e}")
            return {}

    def _parse_abricate_results(self, abricate_output: str) -> Dict:
        """
        Parse Abricate TSV output with detailed hit information

        Returns:
            Dictionary with plasmid detection results including:
            - plasmid_detected: bool
            - genes_found: list of gene names
            - plasmid_types: list of unique plasmid types
            - hits: list of detailed hit dictionaries with position, coverage, identity
        """
        results = {
            'plasmid_detected': False,
            'genes_found': [],
            'plasmid_types': [],
            'hits': []
        }

        lines = abricate_output.strip().split('\n')

        # Skip header and empty lines
        for line in lines[1:]:
            if not line.strip() or line.startswith('#'):
                continue

            parts = line.split('\t')
            if len(parts) >= 10:
                sequence_id = parts[1]  # SEQUENCE column
                gene = parts[4]  # GENE column
                start = parts[2]  # START column
                end = parts[3]  # END column
                coverage = parts[8]  # %COVERAGE column (index 8)
                identity = parts[9]  # %IDENTITY column (index 9)

                # Only count high-quality hits (80%+ coverage and identity)
                try:
                    cov = float(coverage.rstrip('%'))
                    ident = float(identity.rstrip('%'))

                    if cov >= 80.0 and ident >= 80.0:
                        hit_info = {
                            'gene': gene,
                            'sequence': sequence_id,
                            'start': int(start),
                            'end': int(end),
                            'coverage': cov,
                            'identity': ident
                        }
                        results['hits'].append(hit_info)
                        results['genes_found'].append(gene)
                        if gene not in results['plasmid_types']:
                            results['plasmid_types'].append(gene)
                except (ValueError, IndexError):
                    continue

        results['plasmid_detected'] = len(results['genes_found']) > 0
        return results

    def analyze_plasmid_content(self, sequence: Seq, sequence_name: str = "") -> Dict:
        """
        Complete plasmid analysis with priority-based method selection

        Priority:
        1. Abricate (most accurate, detailed hit information)
        2. PlasmidFinder CLI tool
        3. Built-in detection (limited accuracy, fallback only)

        Args:
            sequence: Input sequence
            sequence_name: Name of the sequence

        Returns:
            Dictionary with plasmid detection results
        """
        logger.info(f"Running plasmid analysis on: {sequence_name}")

        # Create temporary FASTA file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as tmp_fasta:
            # Create SeqRecord object
            from Bio.SeqRecord import SeqRecord
            record = SeqRecord(sequence, id=sequence_name, description="")
            SeqIO.write(record, tmp_fasta, 'fasta')
            tmp_path = tmp_fasta.name

        try:
            # PRIORITY 1: Use Abricate with plasmidfinder database (most reliable)
            if self.abricate_available:
                ab_results = self.detect_plasmid_markers_abricate(tmp_path)
                if ab_results:
                    # Calculate confidence based on hits
                    hit_count = len(ab_results.get('hits', []))
                    if hit_count >= 3:
                        confidence = 'HIGH'
                    elif hit_count >= 1:
                        confidence = 'MODERATE'
                    else:
                        confidence = 'NONE'

                    return {
                        'plasmid_detected': ab_results.get('plasmid_detected', False),
                        'method': 'Abricate',
                        'results': ab_results,
                        'markers_found': ab_results.get('genes_found', []),
                        'hits': ab_results.get('hits', []),
                        'hit_count': hit_count,
                        'confidence': confidence,
                        'database': 'plasmidfinder'
                    }

            # PRIORITY 2: Try PlasmidFinder command-line tool
            if self.plasmidfinder_available:
                pf_results = self.detect_plasmid_markers_plasmidfinder(tmp_path)
                if pf_results:
                    return {
                        'plasmid_detected': True,
                        'method': 'PlasmidFinder',
                        'results': pf_results
                    }

            # PRIORITY 3: Fallback to built-in detection (limited accuracy)
            logger.warning("Using built-in detection - results may be inaccurate")
            markers = self.detect_plasmid_markers_builtin(sequence)

            # Calculate plasmid confidence
            plasmid_detected = len(markers) > 0
            marker_count = len(markers)
            unique_markers = len(set(markers.keys()))

            # High confidence = multiple different markers
            confidence = 'HIGH' if unique_markers >= 3 else 'MODERATE' if unique_markers >= 2 else 'LOW'

            return {
                'plasmid_detected': plasmid_detected,
                'method': 'Built-in',
                'markers_found': list(markers.keys()),
                'marker_count': marker_count,
                'unique_markers': unique_markers,
                'confidence': confidence,
                'marker_details': markers,
                'warning': 'Built-in detection is limited - install abricate for accurate results'
            }

        finally:
            # Clean up temp file
            try:
                os.unlink(tmp_path)
            except:
                pass

    def calculate_plasmid_score(self, plasmid_results: Dict) -> float:
        """
        Calculate plasmid engineering score (0-20)

        Args:
            plasmid_results: Results from analyze_plasmid_content

        Returns:
            Plasmid score contribution
        """
        if not plasmid_results.get('plasmid_detected', False):
            return 0.0

        confidence = plasmid_results.get('confidence', 'LOW')
        unique_markers = plasmid_results.get('unique_markers', 0)

        # Score based on confidence and marker count
        if confidence == 'HIGH':
            score = 20
        elif confidence == 'MODERATE':
            score = 15
        else:
            score = 10

        # Bonus for multiple markers
        if unique_markers >= 5:
            score = min(20, score + 5)

        return float(score)


class PlasmidEnhancedDetector:
    """
    Enhanced plasmid detector with comprehensive analysis and reporting

    Features:
    - Direct Abricate file analysis (fastest method)
    - Rich hit metadata (positions, coverage, identity)
    - Enhanced error handling with timeouts
    - Detailed confidence scoring
    - Comprehensive reporting
    """

    def __init__(self, conda_env: Optional[str] = None):
        """
        Initialize the enhanced detector

        Args:
            conda_env: Name of conda environment (default: None)
        """
        self.plasmid_finder = PlasmidFinderIntegration(conda_env=conda_env)

    def analyze_fasta_file_abricate(self, fasta_file: str) -> Dict:
        """
        Analyze a FASTA file directly using Abricate (fastest method)

        This runs abricate directly on the file without loading into memory.

        Args:
            fasta_file: Path to FASTA file

        Returns:
            Analysis results dict with detailed hit information
        """
        if not self.plasmid_finder.abricate_available:
            logger.warning("Abricate not available, falling back to slower method")
            return self.analyze_sequence(fasta_file)

        try:
            # Build command with conda environment if specified
            if self.plasmid_finder.conda_env:
                cmd = ['conda', 'run', '-n', self.plasmid_finder.conda_env, 'abricate', '--db', 'plasmidfinder', fasta_file]
            else:
                cmd = ['abricate', '--db', 'plasmidfinder', fasta_file]

            result = subprocess.run(cmd,
                                    capture_output=True,
                                    text=True,
                                    timeout=3600)  # 1 hour timeout for large files

            if result.returncode != 0:
                logger.warning(f"Abricate failed: {result.stderr}")
                return self.analyze_sequence(fasta_file)

            # Parse the results
            ab_results = self.plasmid_finder._parse_abricate_results(result.stdout)

            # Calculate confidence based on hits
            hit_count = len(ab_results.get('hits', []))
            if hit_count >= 3:
                confidence = 'HIGH'
            elif hit_count >= 1:
                confidence = 'MODERATE'
            else:
                confidence = 'NONE'

            return {
                'file': fasta_file,
                'plasmid_detected': ab_results.get('plasmid_detected', False),
                'method': 'Abricate',
                'markers_found': ab_results.get('genes_found', []),
                'hits': ab_results.get('hits', []),
                'hit_count': hit_count,
                'confidence': confidence,
                'database': 'plasmidfinder',
                'analysis_complete': True
            }

        except subprocess.TimeoutExpired:
            logger.error(f"Abricate timed out on {fasta_file}")
            return {
                'file': fasta_file,
                'plasmid_detected': False,
                'error': 'timeout',
                'analysis_complete': False
            }
        except Exception as e:
            logger.error(f"Error analyzing {fasta_file}: {e}")
            return {
                'file': fasta_file,
                'plasmid_detected': False,
                'error': str(e),
                'analysis_complete': False
            }

    def analyze_sequence(self, fasta_file: str) -> Dict:
        """
        Analyze a FASTA file for plasmid content

        This method loads the FASTA file into memory and analyzes it.
        For large files, use analyze_fasta_file_abricate instead.

        Args:
            fasta_file: Path to FASTA file

        Returns:
            Analysis results dict
        """
        try:
            # Parse the FASTA file
            records = list(SeqIO.parse(fasta_file, "fasta"))
            if not records:
                logger.warning(f"No sequences found in {fasta_file}")
                return None

            # Process first sequence (can be extended for multi-record files)
            record = records[0]
            result = self.analyze_sequence_with_plasmid(record.seq, record.id)

            # Add file path to results
            result['file'] = fasta_file
            return result

        except Exception as e:
            logger.error(f"Error analyzing {fasta_file}: {e}")
            return None

    def analyze_sequence_with_plasmid(self, sequence: Seq, sequence_name: str = "") -> Dict:
        """
        Combined analysis with detailed plasmid detection

        Args:
            sequence: Input sequence
            sequence_name: Name of the sequence

        Returns:
            Combined analysis results with detailed hit information
        """
        logger.info(f"Running plasmid analysis on: {sequence_name}")

        # Plasmid analysis
        plasmid_results = self.plasmid_finder.analyze_plasmid_content(sequence, sequence_name)
        plasmid_score = self.plasmid_finder.calculate_plasmid_score(plasmid_results)

        # Build comprehensive result with detailed hit information
        result = {
            'sequence_name': sequence_name,
            'plasmid_detected': plasmid_results.get('plasmid_detected', False),
            'plasmid_confidence': plasmid_results.get('confidence', 'LOW'),
            'plasmid_markers': plasmid_results.get('markers_found', []),
            'plasmid_score': plasmid_score,
            'detection_method': plasmid_results.get('method', 'Built-in'),
            'analysis_complete': True
        }

        # Add detailed hit information if available (from Abricate)
        if 'hits' in plasmid_results:
            result['hits'] = plasmid_results['hits']
            result['hit_count'] = plasmid_results.get('hit_count', len(plasmid_results['hits']))

        # Add database information if available
        if 'database' in plasmid_results:
            result['database'] = plasmid_results['database']

        # Add warning if using built-in method
        if 'warning' in plasmid_results:
            result['warning'] = plasmid_results['warning']

        return result

    def generate_plasmid_report(self, results: List[Dict]) -> pd.DataFrame:
        """
        Generate plasmid analysis report with detailed hit information

        Args:
            results: List of analysis result dictionaries

        Returns:
            DataFrame with analysis results
        """
        data = []
        for result in results:
            row = {
                'Sequence': result.get('sequence_name', result.get('file', 'Unknown')),
                'Plasmid_Detected': 'YES' if result['plasmid_detected'] else 'NO',
                'Confidence': result.get('plasmid_confidence', result.get('confidence', 'UNKNOWN')),
                'Markers_Found': len(result.get('plasmid_markers', result.get('markers_found', []))),
                'Marker_Names': ', '.join(result.get('plasmid_markers', result.get('markers_found', []))[:5]),
                'Method': result.get('detection_method', result.get('method', 'Unknown')),
                'Plasmid_Score': result.get('plasmid_score', 0)
            }

            # Add hit details if available
            if 'hits' in result and result['hits']:
                hit_details = []
                for hit in result['hits']:
                    hit_details.append(f"{hit['gene']} ({hit['coverage']:.1f}% cov, {hit['identity']:.1f}% id)")
                row['Hit_Details'] = '; '.join(hit_details)
            else:
                row['Hit_Details'] = ''

            data.append(row)

        df = pd.DataFrame(data)

        print("\n" + "=" * 100)
        print("PLASMID ANALYSIS REPORT (Abricate/PlasmidFinder)")
        print("=" * 100)
        print(f"Sequences analyzed: {len(results)}")
        print(f"Plasmids detected: {sum(1 for r in results if r['plasmid_detected'])}")
        print()

        # Print detailed results for sequences with hits
        for result in results:
            if result['plasmid_detected'] and 'hits' in result:
                confidence = result.get('plasmid_confidence', result.get('confidence', 'UNKNOWN'))
                markers = result.get('plasmid_markers', result.get('markers_found', []))
                print(f"Sequence: {result.get('sequence_name', result.get('file', 'Unknown'))}")
                print(f"   Confidence: {confidence}")
                print(f"   Markers: {', '.join(markers)}")
                print(f"   Hit Details:")
                for hit in result['hits']:
                    print(f"      - {hit['gene']}: positions {hit['start']}-{hit['end']} "
                          f"({hit['coverage']:.1f}% coverage, {hit['identity']:.1f}% identity)")
                print()

        return df


def analyze_fasta_files(fasta_files: List[str], conda_env: Optional[str] = None) -> pd.DataFrame:
    """
    Analyze multiple FASTA files for plasmid content

    Args:
        fasta_files: List of FASTA file paths
        conda_env: Name of conda environment (optional)

    Returns:
        DataFrame with results
    """
    detector = PlasmidEnhancedDetector(conda_env=conda_env)
    results = []

    print("Running Plasmid Analysis on FASTA files...")
    print("Using Abricate with plasmidfinder database")
    print("=" * 100)

    for i, fasta_file in enumerate(fasta_files, 1):
        print(f"\n[{i}/{len(fasta_files)}] Analyzing: {fasta_file}")
        result = detector.analyze_fasta_file_abricate(fasta_file)
        if result:
            results.append(result)

            # Print summary for this file
            if result['plasmid_detected']:
                print(f"   PLASMID DETECTED: {result['hit_count']} hits")
                for hit in result.get('hits', []):
                    print(f"      - {hit['gene']}: positions {hit['start']}-{hit['end']}")
            else:
                print(f"   No plasmids detected")

    # Generate final report
    if results:
        detector.generate_plasmid_report(results)

    return detector.generate_plasmid_report(results)


if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print("Usage: python plasmid_finder.py <fasta_file> [conda_env]")
        print("\nExample:")
        print("  python plasmid_finder.py sequence.fasta")
        print("  python plasmid_finder.py sequence.fasta myenv")
        sys.exit(1)

    fasta_file = sys.argv[1]
    conda_env = sys.argv[2] if len(sys.argv) > 2 else None

    detector = PlasmidEnhancedDetector(conda_env=conda_env)
    result = detector.analyze_fasta_file_abricate(fasta_file)

    if result:
        print("\n" + "=" * 100)
        print("ANALYSIS COMPLETE")
        print("=" * 100)
        print(f"File: {fasta_file}")
        print(f"Plasmid Detected: {result['plasmid_detected']}")
        print(f"Confidence: {result.get('confidence', 'N/A')}")
        print(f"Hit Count: {result.get('hit_count', 0)}")
        print(f"Method: {result.get('method', 'N/A')}")
