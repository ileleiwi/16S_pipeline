#!/usr/bin/env python3
import gzip
import sys
from Bio import SeqIO
from Bio.Seq import Seq

def examine_read(file_path, read_id):
    """Examine a specific read in a FASTQ file"""
    print(f"Examining file: {file_path} for read: {read_id}")
    
    with gzip.open(file_path, 'rt') if file_path.endswith('.gz') else open(file_path, 'r') as handle:
        for record in SeqIO.parse(handle, "fastq"):
            if record.id == read_id:
                print(f"\nFound read: {record.id}")
                seq = str(record.seq)
                qual = record.letter_annotations["phred_quality"]
                
                print(f"Sequence length: {len(seq)}")
                print(f"Quality scores length: {len(qual)}")
                
                print("\nSequence:")
                print(seq)
                
                print("\nQuality string (ASCII):")
                print(record.format("fastq").split('\n')[3])
                
                print("\nQuality values (numeric):")
                print(qual)
                
                # Calculate ASCII offsets to check for unusual quality encoding
                ascii_repr = record.format("fastq").split('\n')[3]
                print("\nASCII code points:")
                print([ord(c) for c in ascii_repr])
                
                # Check if quality scores have been properly parsed
                if len(seq) != len(qual):
                    print(f"\n⚠️ WARNING: Sequence and quality lengths differ: {len(seq)} vs {len(qual)}")
                    
                    # Debug the quality string character by character
                    print("\nCharacter-by-character analysis of quality string:")
                    for i, c in enumerate(ascii_repr):
                        print(f"Pos {i+1}: '{c}' (ASCII: {ord(c)})")
                
                return
    
    print(f"Read {read_id} not found in {file_path}")

def count_and_validate_fastq(file_path):
    """Count reads and check for format errors in a FASTQ file"""
    print(f"\nAnalyzing file: {file_path}")
    
    total_reads = 0
    format_errors = 0
    format_error_examples = []
    
    try:
        with gzip.open(file_path, 'rt') if file_path.endswith('.gz') else open(file_path, 'r') as handle:
            for record in SeqIO.parse(handle, "fastq"):
                total_reads += 1
                
                seq_len = len(record.seq)
                qual_len = len(record.letter_annotations["phred_quality"])
                
                if seq_len != qual_len:
                    format_errors += 1
                    if len(format_error_examples) < 3:
                        format_error_examples.append({
                            "id": record.id,
                            "seq_len": seq_len,
                            "qual_len": qual_len
                        })
    except Exception as e:
        print(f"Error parsing file: {str(e)}")
        
    print(f"Total reads: {total_reads}")
    print(f"Format errors: {format_errors}")
    
    if format_errors > 0:
        print("\nFormat error examples:")
        for example in format_error_examples:
            print(f"Read ID: {example['id']}")
            print(f"Sequence length: {example['seq_len']}")
            print(f"Quality length: {example['qual_len']}")
            print()

def check_output_file(file_path):
    """Check gzip integrity and FASTQ format validity"""
    print(f"\nChecking file integrity: {file_path}")
    
    # Check if the file is a valid gzip file
    try:
        with gzip.open(file_path, 'rb') as f:
            # Try to read a small chunk to verify gzip format
            f.read(100)
        print("✓ Valid gzip format")
    except Exception as e:
        print(f"✗ Invalid gzip format: {str(e)}")
        return
    
    # Check FASTQ format
    try:
        record_count = 0
        with gzip.open(file_path, 'rt') as f:
            for record in SeqIO.parse(f, "fastq"):
                record_count += 1
                if record_count > 1000:  # Just check the first 1000 records
                    break
        print(f"✓ Valid FASTQ format (checked {record_count} records)")
    except ValueError as e:
        print(f"✗ FASTQ format error: {str(e)}")
    except Exception as e:
        print(f"✗ Error checking FASTQ format: {str(e)}")

# Main script
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python diagnostic.py <fastq_file> [read_id]")
        sys.exit(1)
    
    fastq_file = sys.argv[1]
    
    if len(sys.argv) >= 3:
        read_id = sys.argv[2]
        examine_read(fastq_file, read_id)
    else:
        count_and_validate_fastq(fastq_file)
        check_output_file(fastq_file)
