#!/usr/bin/env python3
"""
Barcode Format Diagnostic Tool

This script analyzes dual-indexed sequencing data to determine the optimal barcode format
for demultiplexing. It tests various combinations of index1 and index2 barcodes against
the mapping file and reports which format gives the highest match rate.

Usage:
    python barcode_diagnostic.py --index1 /path/to/index1.fastq.gz 
                                --index2 /path/to/index2.fastq.gz 
                                --mapping_file /path/to/mapping.txt 
                                --sample_size 10000
                                --verbose
"""

import os
import sys
import gzip
import argparse
from collections import defaultdict, Counter
from Bio import SeqIO
import pandas as pd
import time
from itertools import islice

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Barcode format diagnostic tool for dual-indexed sequencing data."
    )
    
    # Main required arguments
    parser.add_argument("--index1", required=True,
                       help="Path to index1 FASTQ file (gzipped)")
    parser.add_argument("--index2", required=True,
                       help="Path to index2 FASTQ file (gzipped)")
    parser.add_argument("--mapping_file", required=True,
                       help="Path to mapping file with sample and barcode information")
    
    # Sample size for testing
    parser.add_argument("--sample_size", type=int, default=10000,
                       help="Number of reads to sample for barcode format testing")
    parser.add_argument("--min_success", type=float, default=1.0,
                       help="Minimum success rate to consider a format successful (percent)")
    parser.add_argument("--no_rc", action="store_true",
                       help="Skip testing reverse complement formats (faster)")
    parser.add_argument("--verbose", action="store_true",
                       help="Show detailed diagnostic information")
    parser.add_argument("--output", 
                       help="Output file for diagnostic results (optional)")
    
    args = parser.parse_args()
    
    # Validate that input files exist
    for file_path in [args.index1, args.index2, args.mapping_file]:
        if not os.path.exists(file_path):
            print(f"Error: File not found: {file_path}", file=sys.stderr)
            sys.exit(1)
            
    return args

def read_mapping_file(mapping_path):
    """
    Read the mapping file containing sample information and barcodes.
    
    Args:
        mapping_path: Path to the tab-delimited mapping file
        
    Returns:
        DataFrame containing mapping information
    """
    try:
        # Skip comment lines starting with #, but keep the header row
        with open(mapping_path, 'r') as f:
            # Find the first line that starts with #SampleID
            header_line = None
            for line in f:
                if line.startswith('#SampleID'):
                    header_line = line
                    break
                    
            if header_line is None:
                # Try without the # prefix
                with open(mapping_path, 'r') as f2:
                    for line in f2:
                        if line.startswith('SampleID'):
                            header_line = line
                            break
                
                if header_line is None:
                    raise ValueError("Mapping file does not contain a header starting with #SampleID or SampleID")
        
        # Read the file using pandas
        mapping_df = pd.read_csv(
            mapping_path, 
            sep='\t',
            comment='#',
            header=0,
            names=['SampleID', 'BarcodeSequence', 'LinkerPrimerSequence', 'Plate', 'Well', 'Description']
        )
        
        # Rename SampleID column if it has a # prefix in the dataframe
        if mapping_df.columns[0].startswith('#'):
            mapping_df.rename(columns={mapping_df.columns[0]: 'SampleID'}, inplace=True)
            
        # Check that required columns exist
        required_cols = ['SampleID', 'BarcodeSequence']
        missing_cols = [col for col in required_cols if col not in mapping_df.columns]
        if missing_cols:
            print(f"Warning: Missing required columns in mapping file: {', '.join(missing_cols)}")
            print("Available columns:", ", ".join(mapping_df.columns))
            if 'SampleID' in missing_cols:
                raise ValueError("SampleID column is required in mapping file")
                
        # Handle missing BarcodeSequence column
        if 'BarcodeSequence' not in mapping_df.columns:
            closest_match = next((col for col in mapping_df.columns if 'barcode' in col.lower()), None)
            if closest_match:
                print(f"Using '{closest_match}' as barcode column instead of 'BarcodeSequence'")
                mapping_df['BarcodeSequence'] = mapping_df[closest_match]
            else:
                raise ValueError("Could not find a column for barcode sequences")
                
        # Check for and remove any rows with duplicate sample IDs
        if mapping_df['SampleID'].duplicated().any():
            duplicates = mapping_df['SampleID'][mapping_df['SampleID'].duplicated()].unique()
            print(f"Warning: Found duplicate sample IDs: {', '.join(duplicates)}")
            print("Keeping only the first occurrence of each sample ID")
            mapping_df = mapping_df.drop_duplicates(subset=['SampleID'])
            
        return mapping_df
        
    except Exception as e:
        print(f"Error reading mapping file: {str(e)}", file=sys.stderr)
        sys.exit(1)

def create_sample_barcode_map(mapping_df):
    """
    Create a mapping from barcodes to sample IDs.
    
    Args:
        mapping_df: DataFrame from mapping file
        
    Returns:
        Dictionary mapping barcodes to sample IDs
    """
    barcode_to_sample = {}
    
    for _, row in mapping_df.iterrows():
        sample_id = row['SampleID']
        barcode = row['BarcodeSequence']
        
        # Skip empty barcodes
        if pd.isna(barcode) or barcode == '':
            continue
            
        # Remove any whitespace
        barcode = barcode.strip()
        
        # Add mapping for this barcode
        if barcode in barcode_to_sample:
            print(f"Warning: Barcode {barcode} is assigned to multiple samples. Using first occurrence.")
        else:
            barcode_to_sample[barcode] = sample_id
            
    return barcode_to_sample

def extract_barcodes_from_index_files(index1_path, index2_path, sample_size, verbose=False):
    """
    Extract barcodes from index FASTQ files.
    
    Args:
        index1_path: Path to index1 FASTQ file
        index2_path: Path to index2 FASTQ file
        sample_size: Number of records to sample
        verbose: Whether to show verbose output
        
    Returns:
        Dictionary of read IDs to barcode information
    """
    try:
        if verbose:
            print(f"Extracting barcodes from:\n  Index1: {index1_path}\n  Index2: {index2_path}")
            print(f"Sampling {sample_size:,} reads...")
            
        barcode_dict = {}
        total_processed = 0
        
        # Open and parse both index files simultaneously
        with gzip.open(index1_path, 'rt') as f1, gzip.open(index2_path, 'rt') as f2:
            # Create iterators for both files
            iter1 = SeqIO.parse(f1, 'fastq')
            iter2 = SeqIO.parse(f2, 'fastq')
            
            for record1, record2 in zip(islice(iter1, sample_size), islice(iter2, sample_size)):
                seq_id1 = record1.id
                seq_id2 = record2.id
                
                # Ensure sequence IDs match (sanity check)
                if seq_id1 != seq_id2:
                    if verbose:
                        print(f"Warning: Sequence ID mismatch: {seq_id1} != {seq_id2}")
                    continue
                
                # Get raw sequence ID without /1 or /2 suffix if present
                seq_id = seq_id1
                if seq_id.endswith('/1') or seq_id.endswith('/2'):
                    seq_id = seq_id[:-2]
                
                # Extract barcode sequences
                barcode1 = str(record1.seq)
                barcode2 = str(record2.seq)
                
                # Store in dictionary
                barcode_dict[seq_id] = {
                    'barcode1': barcode1,
                    'barcode2': barcode2
                }
                
                total_processed += 1
                
                # Print progress every 1,000 records in verbose mode
                if verbose and total_processed % 1000 == 0:
                    print(f"  Processed {total_processed:,} index records...")
                    
        if verbose:
            print(f"Extracted barcode information for {total_processed:,} reads")
            
        # Return early if we didn't find any barcodes
        if not barcode_dict:
            print("Error: No barcodes found in index files", file=sys.stderr)
            sys.exit(1)
            
        return barcode_dict
        
    except Exception as e:
        print(f"Error processing index files: {str(e)}", file=sys.stderr)
        sys.exit(1)

def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 
                 'N': 'N', 
                 'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 
                 'n': 'n'}
    return ''.join([complement.get(base, base) for base in seq[::-1]])

def test_barcode_formats(barcode_info, barcode_to_sample, no_rc=False, verbose=False):
    """
    Test different barcode combination formats to determine which format works best.
    
    Args:
        barcode_info: Dictionary mapping sequence IDs to barcodes
        barcode_to_sample: Dictionary mapping barcodes to sample IDs
        no_rc: Whether to skip testing reverse complement formats
        verbose: Whether to show verbose output
        
    Returns:
        Dictionary with success rates for different barcode formats
    """
    # Define barcode format functions - each takes barcode1 and barcode2 and returns a formatted barcode
    barcode_formats = {
        "barcode1": lambda b1, b2: b1,
        "barcode2": lambda b1, b2: b2,
        "barcode1+barcode2": lambda b1, b2: b1 + b2,
        "barcode2+barcode1": lambda b1, b2: b2 + b1,
        "barcode1-barcode2": lambda b1, b2: f"{b1}-{b2}",
        "barcode2-barcode1": lambda b1, b2: f"{b2}-{b1}"
    }
    
    # Add reverse and reverse complement formats
    barcode_formats["barcode1_rev"] = lambda b1, b2: b1[::-1]
    barcode_formats["barcode2_rev"] = lambda b1, b2: b2[::-1]
    barcode_formats["barcode1+barcode2_rev"] = lambda b1, b2: b1 + b2[::-1]
    barcode_formats["barcode2_rev+barcode1"] = lambda b1, b2: b2[::-1] + b1
    barcode_formats["barcode1-barcode2_rev"] = lambda b1, b2: f"{b1}-{b2[::-1]}"
    barcode_formats["barcode2_rev-barcode1"] = lambda b1, b2: f"{b2[::-1]}-{b1}"
    
    # Add reverse complement formats if not disabled
    if not no_rc:
        barcode_formats["barcode1_rc"] = lambda b1, b2: reverse_complement(b1)
        barcode_formats["barcode2_rc"] = lambda b1, b2: reverse_complement(b2)
        barcode_formats["barcode1+barcode2_rc"] = lambda b1, b2: b1 + reverse_complement(b2)
        barcode_formats["barcode2_rc+barcode1"] = lambda b1, b2: reverse_complement(b2) + b1
        barcode_formats["barcode1-barcode2_rc"] = lambda b1, b2: f"{b1}-{reverse_complement(b2)}"
        barcode_formats["barcode2_rc-barcode1"] = lambda b1, b2: f"{reverse_complement(b2)}-{b1}"
        barcode_formats["barcode1_rc+barcode2"] = lambda b1, b2: reverse_complement(b1) + b2
        barcode_formats["barcode1+barcode2_rc"] = lambda b1, b2: b1 + reverse_complement(b2)
        barcode_formats["barcode1_rc-barcode2"] = lambda b1, b2: f"{reverse_complement(b1)}-{b2}"
        barcode_formats["barcode1-barcode2_rc"] = lambda b1, b2: f"{b1}-{reverse_complement(b2)}"
    
    # Initialize counters for each format
    format_counts = {fmt: {'matched': 0, 'total': 0} for fmt in barcode_formats}
    
    # Keep track of all found barcodes for diagnostics
    found_barcodes = {}
    
    # Keep track of successful matches for each format
    sample_matches = {fmt: defaultdict(int) for fmt in barcode_formats}
    
    # Keep a few examples of each format for diagnostic purposes
    format_examples = {fmt: [] for fmt in barcode_formats}
    
    # Check each barcode against all formats
    total_reads = len(barcode_info)
    if verbose:
        print(f"Testing {len(barcode_formats)} barcode formats on {total_reads:,} reads...")
    
    for i, (read_id, barcodes) in enumerate(barcode_info.items()):
        barcode1 = barcodes['barcode1']
        barcode2 = barcodes['barcode2']
        
        # Store this for diagnostic output
        found_barcodes[read_id] = {'barcode1': barcode1, 'barcode2': barcode2}
        
        # Try each barcode format
        for fmt_name, fmt_func in barcode_formats.items():
            format_counts[fmt_name]['total'] += 1
            
            try:
                test_barcode = fmt_func(barcode1, barcode2)
                
                # Store a few examples of each format
                if len(format_examples[fmt_name]) < 3:
                    format_examples[fmt_name].append(test_barcode)
                
                # Check if this barcode maps to a sample
                if test_barcode in barcode_to_sample:
                    format_counts[fmt_name]['matched'] += 1
                    sample_id = barcode_to_sample[test_barcode]
                    sample_matches[fmt_name][sample_id] += 1
            except Exception as e:
                if verbose:
                    print(f"Error with format {fmt_name} for read {read_id}: {str(e)}")
                # If any format function fails, just continue
                pass
        
        # Print progress in verbose mode
        if verbose and (i+1) % 1000 == 0:
            print(f"  Processed {i+1:,} of {total_reads:,} reads...")
    
    # Calculate success rates
    format_success = {}
    for fmt_name, counts in format_counts.items():
        if counts['total'] > 0:
            success_rate = counts['matched'] / counts['total'] * 100
            samples_matched = len(sample_matches[fmt_name])
            format_success[fmt_name] = {
                'success_rate': success_rate,
                'matched': counts['matched'],
                'total': counts['total'],
                'samples_matched': samples_matched,
                'examples': format_examples[fmt_name]
            }
        else:
            format_success[fmt_name] = {
                'success_rate': 0,
                'matched': 0,
                'total': 0,
                'samples_matched': 0,
                'examples': []
            }
    
    # Find best format
    best_format = max(format_success.items(), key=lambda x: x[1]['success_rate'])
    
    # Diagnostic information
    diagnostic = {
        'success_rates': format_success,
        'best_format': best_format[0],
        'best_rate': best_format[1]['success_rate'],
        'sample_barcodes': list(found_barcodes.values())[:10],  # First 10 barcodes for reference
        'barcode_to_sample_sizes': {
            'mapping_barcodes': len(barcode_to_sample),
            'barcode_info': len(barcode_info)
        },
        'barcode_examples': format_examples,
        'sample_matches': {fmt: dict(matches) for fmt, matches in sample_matches.items()}
    }
    
    return diagnostic

def analyze_barcode_lengths(barcode_info, barcode_to_sample):
    """
    Analyze barcode lengths from index files and mapping file.
    
    Args:
        barcode_info: Dictionary mapping sequence IDs to barcodes
        barcode_to_sample: Dictionary mapping barcodes to sample IDs
        
    Returns:
        Dictionary with barcode length statistics
    """
    # Get lengths from index files
    barcode1_lengths = Counter()
    barcode2_lengths = Counter()
    for read_info in barcode_info.values():
        barcode1_lengths[len(read_info['barcode1'])] += 1
        barcode2_lengths[len(read_info['barcode2'])] += 1
    
    # Get lengths from mapping file
    mapping_barcode_lengths = Counter()
    for barcode in barcode_to_sample.keys():
        mapping_barcode_lengths[len(barcode)] += 1
    
    # Calculate statistics
    stats = {
        'index1': {
            'lengths': dict(barcode1_lengths),
            'most_common': barcode1_lengths.most_common(1)[0][0] if barcode1_lengths else 0
        },
        'index2': {
            'lengths': dict(barcode2_lengths),
            'most_common': barcode2_lengths.most_common(1)[0][0] if barcode2_lengths else 0
        },
        'mapping': {
            'lengths': dict(mapping_barcode_lengths),
            'most_common': mapping_barcode_lengths.most_common(1)[0][0] if mapping_barcode_lengths else 0
        }
    }
    
    # Add combined length stats
    combined_length = stats['index1']['most_common'] + stats['index2']['most_common']
    stats['combined'] = {
        'length': combined_length,
        'matches_mapping': combined_length == stats['mapping']['most_common']
    }
    
    # Check if mapping barcode length equals sum of index barcode lengths (common case)
    stats['likely_combined'] = stats['mapping']['most_common'] == stats['index1']['most_common'] + stats['index2']['most_common']
    
    # Check if mapping barcode length equals either index length (single index case)
    stats['likely_single_index'] = stats['mapping']['most_common'] in [stats['index1']['most_common'], stats['index2']['most_common']]
    
    return stats

def print_diagnostic_report(diagnostic, length_stats, barcode_to_sample, args, output_file=None):
    """
    Print or write a detailed diagnostic report.
    
    Args:
        diagnostic: Dictionary with diagnostic information
        length_stats: Dictionary with barcode length statistics
        barcode_to_sample: Dictionary mapping barcodes to sample IDs
        args: Command line arguments
        output_file: Optional file to write report to
    """
    # Initialize output
    output = []
    
    def print_or_write(line):
        output.append(line)
        if not output_file:  # Print immediately if not writing to file
            print(line)
    
    # Format success rates as a table
    print_or_write("\n=== BARCODE FORMAT DIAGNOSTIC REPORT ===\n")
    print_or_write(f"Sample size: {args.sample_size:,} reads")
    print_or_write(f"Mapping file: {args.mapping_file}")
    print_or_write(f"Index1 file: {args.index1}")
    print_or_write(f"Index2 file: {args.index2}")
    print_or_write(f"Reverse complement formats: {'Tested' if not args.no_rc else 'Skipped'}")
    print_or_write("\n=== BARCODE FORMAT SUCCESS RATES ===\n")
    print_or_write(f"{'Format':<30} | {'Success Rate':<12} | {'Matched/Total':<15} | {'Samples':<10} | {'Example'}")
    print_or_write(f"{'-'*30}-+-{'-'*12}-+-{'-'*15}-+-{'-'*10}-+-{'-'*20}")
    
    # Sort formats by success rate
    sorted_formats = sorted(
        diagnostic['success_rates'].items(), 
        key=lambda x: x[1]['success_rate'], 
        reverse=True
    )
    
    # Find formats with success rates above the minimum threshold
    min_success = args.min_success
    viable_formats = []
    
    for fmt_name, stats in sorted_formats:
        success_rate = stats['success_rate']
        matched = stats['matched']
        total = stats['total']
        samples = stats['samples_matched']
        examples = stats['examples'][0] if stats['examples'] else ""
        
        # Format the ratio
        ratio = f"{matched:,}/{total:,}"
        
        # Mark the best format with a star
        marker = "★ " if fmt_name == diagnostic['best_format'] else "  "
        
        # Add to viable formats if above threshold
        if success_rate >= min_success:
            viable_formats.append((fmt_name, success_rate))
        
        # Add to the report
        print_or_write(f"{marker}{fmt_name:<28} | {success_rate:<12.2f} | {ratio:<15} | {samples:<10} | {examples}")
    
    # Print barcode length analysis
    print_or_write("\n=== BARCODE LENGTH ANALYSIS ===\n")
    print_or_write("Index1 barcode lengths:")
    for length, count in sorted(length_stats['index1']['lengths'].items()):
        print_or_write(f"  {length} bp: {count:,} reads")
        
    print_or_write("\nIndex2 barcode lengths:")
    for length, count in sorted(length_stats['index2']['lengths'].items()):
        print_or_write(f"  {length} bp: {count:,} reads")
        
    print_or_write("\nMapping file barcode lengths:")
    for length, count in sorted(length_stats['mapping']['lengths'].items()):
        print_or_write(f"  {length} bp: {count:,} samples")
    
    print_or_write("\nLength Analysis:")
    print_or_write(f"  Most common index1 length: {length_stats['index1']['most_common']} bp")
    print_or_write(f"  Most common index2 length: {length_stats['index2']['most_common']} bp")
    print_or_write(f"  Most common mapping length: {length_stats['mapping']['most_common']} bp")
    print_or_write(f"  Combined index length: {length_stats['combined']['length']} bp")
    
    if length_stats['likely_combined']:
        print_or_write("  ✓ Mapping barcode length equals sum of index lengths (likely dual-indexed)")
    elif length_stats['likely_single_index']:
        print_or_write("  ✓ Mapping barcode length equals one index length (likely single-indexed)")
    else:
        print_or_write("  ⚠ Mapping barcode length doesn't match index lengths (unusual pattern)")
    
    # Print a sample of barcodes from index files
    print_or_write("\n=== SAMPLE BARCODES FROM INDEX FILES ===\n")
    for i, barcode in enumerate(diagnostic['sample_barcodes'][:5]):
        print_or_write(f"Read {i+1}:")
        print_or_write(f"  Index1: {barcode['barcode1']}")
        print_or_write(f"  Index2: {barcode['barcode2']}")
        
        # For each sample, show what it would be with the best format
        best_format_func = None
        try:
            if diagnostic['best_format'] == "barcode1":
                best_result = barcode['barcode1']
            elif diagnostic['best_format'] == "barcode2":
                best_result = barcode['barcode2']
            elif diagnostic['best_format'] == "barcode1+barcode2":
                best_result = barcode['barcode1'] + barcode['barcode2']
            elif diagnostic['best_format'] == "barcode1+barcode2_rev":
                best_result = barcode['barcode1'] + barcode['barcode2'][::-1]
            elif diagnostic['best_format'] == "barcode1+barcode2_rc":
                best_result = barcode['barcode1'] + reverse_complement(barcode['barcode2'])
            else:
                best_result = f"(complex format: {diagnostic['best_format']})"
        except:
            best_result = "(could not generate)"
            
        print_or_write(f"  Best format ({diagnostic['best_format']}): {best_result}")
    
    # Print a sample of barcodes from mapping file
    print_or_write("\n=== SAMPLE BARCODES FROM MAPPING FILE ===\n")
    mapping_samples = list(diagnostic['success_rates'][diagnostic['best_format']].get('sample_matches', {}).items())[:5]
    for i, (barcode, sample_id) in enumerate(list(barcode_to_sample.items())[:5]):
        print_or_write(f"Sample {i+1}:")
        print_or_write(f"  ID: {sample_id}")
        print_or_write(f"  Barcode: {barcode}")
    
    # Print recommendations
    print_or_write("\n=== RECOMMENDATIONS ===\n")
    
    if viable_formats:
        print_or_write(f"✓ BEST FORMAT: {diagnostic['best_format']} (Success rate: {diagnostic['best_rate']:.2f}%)")
        if len(viable_formats) > 1:
            print_or_write("  Other viable formats:")
            for fmt, rate in viable_formats[1:5]:  # List next 4 viable formats
                print_or_write(f"  - {fmt} ({rate:.2f}%)")
    else:
        print_or_write("⚠ WARNING: No format achieved the minimum success rate threshold.")
        print_or_write(f"  Best performing format was {diagnostic['best_format']} with only {diagnostic['best_rate']:.2f}%")
        print_or_write("  Consider checking your barcode files and mapping file for compatibility.")
    
    if diagnostic['best_rate'] < 50:
        print_or_write("\n⚠ Low match rate detected. Possible issues:")
        print_or_write("  - Barcodes in mapping file may not match those in index files")
        print_or_write("  - Index files may be misaligned or from a different run")
        print_or_write("  - Mapping file format may be unexpected (check for hidden characters)")
        print_or_write("  - Try looking at the actual contents of your mapping file and index files")
    
    # Write to file if specified
    if output_file:
        with open(output_file, 'w') as f:
            for line in output:
                f.write(line + "\n")
        print(f"Diagnostic report written to: {output_file}")

def main():
    """Main function to execute the script."""
    # Parse command line arguments
    args = parse_arguments()
    
    # Time the entire process
    start_time = time.time()
    
    # Read mapping file
    print(f"Reading mapping file: {args.mapping_file}")
    mapping_df = read_mapping_file(args.mapping_file)
    print(f"Found {len(mapping_df)} samples in mapping file")
    
    # Create mapping from barcodes to samples
    barcode_to_sample = create_sample_barcode_map(mapping_df)
    print(f"Created mapping for {len(barcode_to_sample)} unique barcodes")
    
    # Extract barcodes from index files
    barcode_info = extract_barcodes_from_index_files(
        args.index1, 
        args.index2, 
        args.sample_size,
        args.verbose
    )
    
    # Analyze barcode lengths
    length_stats = analyze_barcode_lengths(barcode_info, barcode_to_sample)
    
    # Test different barcode formats
    print(f"Testing barcode formats...")
    diagnostic = test_barcode_formats(
        barcode_info,
        barcode_to_sample,
        args.no_rc,
        args.verbose
    )
    
    # Print diagnostic report
    print_diagnostic_report(diagnostic, length_stats, barcode_to_sample, args, args.output)
    
    # Report total execution time
    elapsed = time.time() - start_time
    print(f"\nDiagnostic completed in {elapsed:.2f} seconds")

if __name__ == "__main__":
    main()
