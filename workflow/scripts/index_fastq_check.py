#!/usr/bin/env python3
"""
Sequencing Data Analysis Script

This script processes Illumina index FASTQ files to extract barcodes and sequencing run information.
It creates a barcode file and a summary report that checks congruency between index files.

Usage:
    python index_fastq_check.py --index1 /path/to/index1.fastq.gz \
        --index2 /path/to/index2.fastq.gz --outdir /path/to/output/directory

    Single-index runs (EMP-style barcodes, or runs submitted with inline
    barcodes) omit --index2: the 12 nt barcodes are read from I1 alone and
    there is no I2 file, e.g.
    data/raw/20240209_16S/Undetermined_S0_L001_I1_001.fastq.gz.
"""

from Bio import SeqIO
import gzip
import os
import hashlib
import argparse
import sys
from collections import Counter

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Process index FASTQ files to extract barcodes and sequencing information.")
    
    parser.add_argument("--index1", required=True, help="Path to index1 FASTQ file (gzipped)")
    # Index2 is optional: single-index libraries (EMP-style barcodes, or runs
    # submitted with inline barcodes) ship an I1 read only.
    parser.add_argument("--index2", default=None,
                        help="Path to index2 FASTQ file (gzipped). Omit for single-index runs.")
    parser.add_argument("--outdir", required=True, help="Path to output directory")
    
    # Optional arguments
    parser.add_argument("--index1_md5", default="",
                        help="Expected MD5 checksum for index1 (blank to skip the check)")
    parser.add_argument("--index2_md5", default="",
                        help="Expected MD5 checksum for index2 (blank to skip the check)")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose output")
    
    return parser.parse_args()

def calculate_md5(filename):
    """Calculate MD5 hash for a file."""
    hash_md5 = hashlib.md5()
    with open(filename, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

def get_instrument_name(instrument_code, flowcell_id):
    """Infer instrument name from instrument code and flowcell ID."""
    # Extract flowcell pattern (last 4 characters) if available
    flowcell_pattern = flowcell_id[-4:] if len(flowcell_id) >= 4 else ""
    
    # Check instrument code pattern first
    if instrument_code.startswith('HWI-M') or instrument_code.startswith('M'):
        return "MiSeq"
    elif instrument_code.startswith('HWUSI'):
        return "Genome Analyzer IIx (GAIIx)"
    elif instrument_code.startswith('HWI-D'):
        return "HiSeq 2000/2500"
    elif instrument_code.startswith('K'):
        return "HiSeq 3000/4000"
    elif instrument_code.startswith('N'):
        return "NextSeq 500/550"
    elif instrument_code.startswith('A'):
        return "NovaSeq"
    elif instrument_code.startswith('V'):
        return "NextSeq 2000"
    elif instrument_code.startswith('AA'):
        return "NextSeq 2000 P1/P2/P3"
    elif instrument_code.startswith('_2225'):
        return "NextSeq P4"
    elif instrument_code.startswith('H'):
        return "NovaSeq S1/S2/S4"
    
    # Then check flowcell pattern
    if flowcell_pattern == "AAXX":
        return "Genome Analyzer"
    elif flowcell_pattern == "BCXX":
        return "HiSeq v1.5"
    elif flowcell_pattern == "ACXX":
        return "HiSeq High-Output v3"
    elif flowcell_pattern == "ANXX":
        return "HiSeq High-Output v4"
    elif flowcell_pattern == "ADXX":
        return "HiSeq RR v1"
    elif flowcell_pattern in ["AMXX", "BCXX"]:
        return "HiSeq RR v2"
    elif flowcell_pattern == "ALXX":
        return "HiSeqX"
    elif flowcell_pattern in ["BGXX", "AGXX"]:
        return "High-Output NextSeq"
    elif flowcell_pattern == "AFXX":
        return "Mid-Output NextSeq"
    elif len(instrument_code) == 5 and any(c.isdigit() for c in instrument_code):
        return "MiSeq"
    else:
        return "Unknown Instrument"

def extract_fastq_header_info(record):
    """Extract sequencing information from a FASTQ header."""
    # Example Illumina header format:
    # @INSTRUMENT:RUN:FLOWCELL:LANE:TILE:X:Y READ:IS_FILTERED:CONTROL_NUMBER:BARCODE
    header_parts = record.description.split(" ")[0].split(":")
    
    if len(header_parts) >= 4:
        instrument_code = header_parts[0].lstrip('@')
        run = header_parts[1]
        flowcell = header_parts[2]
        lane = header_parts[3]
        
        # Infer instrument name
        instrument_name = get_instrument_name(instrument_code, flowcell)
        
        return {
            "instrument_code": instrument_code,
            "instrument_name": instrument_name,
            "run": run,
            "flowcell": flowcell,
            "lane": lane
        }
    else:
        return {
            "instrument_code": "unknown",
            "instrument_name": "Unknown Instrument",
            "run": "unknown",
            "flowcell": "unknown",
            "lane": "unknown"
        }

def _empty_file_info(filename, expected_md5):
    """Build the per-file metrics container used by the summary writers."""
    return {
        "filename": os.path.basename(filename),
        "md5": calculate_md5(filename),
        "expected_md5": (expected_md5 or "").strip().lower(),
        "read_count": 0,
        "instrument_codes": set(),
        "instrument_names": set(),
        "runs": set(),
        "flowcells": set(),
        "lanes": set(),
        "header_info": []
    }


def write_single_index_summary(summary_output_file, file_info, index_lengths, barcodes):
    """Write the summary for a single-index (I1 only) run.

    There is no second index file to check congruency against, so instead of the
    dual-index table this reports the I1 metrics plus the barcode composition
    (index length distribution and unique barcode count) that demultiplexing and
    the barcode-format diagnostic depend on.
    """
    index1 = file_info["index1"]
    expected = index1["expected_md5"]
    md5_label = "N/A"
    if expected:
        md5_label = "Yes" if index1["md5"].lower() == expected else "No"

    lengths = index_lengths["index1"]
    total_barcodes = sum(lengths.values())
    unique_barcodes = len(set(entry[0] for entry in barcodes.values()))
    most_common = sorted(lengths.items(), key=lambda kv: kv[1], reverse=True)

    with open(summary_output_file, "w") as out:
        out.write("# Sequencing Run Information (single-index run)\n\n")
        out.write("Index layout: single (I1 only; no I2 file to compare)\n\n")

        out.write("## File Metrics\n\n")
        out.write(f"{'Parameter':<24}{'Index1':<40}\n")
        out.write(f"{'-' * 24}{'-' * 40}\n")
        out.write(f"{'Filename':<24}{index1['filename']:<40}\n")
        out.write(f"{'MD5':<24}{index1['md5']:<40}\n")
        out.write(f"{'MD5 matches expected':<24}{md5_label:<40}\n")
        out.write(f"{'Read count':<24}{index1['read_count']:<40}\n")
        out.write(f"{'Instrument code':<24}{', '.join(sorted(index1['instrument_codes'])):<40}\n")
        out.write(f"{'Instrument name':<24}{', '.join(sorted(index1['instrument_names'])):<40}\n")
        out.write(f"{'Run ID':<24}{', '.join(sorted(index1['runs'])):<40}\n")
        out.write(f"{'Flowcell':<24}{', '.join(sorted(index1['flowcells'])):<40}\n")
        out.write(f"{'Lane':<24}{', '.join(sorted(index1['lanes'])):<40}\n")

        out.write("\n## Index (barcode) Composition\n\n")
        out.write(f"Unique barcodes: {unique_barcodes}\n")
        for length, count in sorted(lengths.items()):
            pct = (count / total_barcodes * 100) if total_barcodes else 0.0
            out.write(f"  {length} bp: {count:,} reads ({pct:.2f}%)\n")

        out.write("\n## Summary\n\n")
        if most_common:
            out.write(f"Most common index length: {most_common[0][0]} bp\n")
        if len(lengths) > 1:
            out.write("- Index reads have mixed lengths; barcode matching must be length aware\n")
        if md5_label == "No":
            out.write("- MD5 checksum mismatch (Index1)\n")
        verdict = "Issues detected" if (md5_label == "No" or len(lengths) > 1) else "Good"
        out.write(f"Overall congruency: {verdict}\n")


def process_fastq_files(index1_file, index2_file, output_dir, verbose=False,
                        index1_md5_expected="", index2_md5_expected=""):
    """Process FASTQ files and generate output files."""
    # Create output directory structure
    data_dir = os.path.join(output_dir, "data")
    os.makedirs(data_dir, exist_ok=True)
    
    if verbose:
        print(f"Processing index files:\n  Index1: {index1_file}\n  Index2: {index2_file}")
        print(f"Output directory: {output_dir}")
    
    # Dictionary to store barcode sequences
    barcodes = {}

    # File metrics and header information. A single-index run has no index2
    # entry, which every block below tests for.
    single_index = not index2_file
    file_info = {"index1": _empty_file_info(index1_file, index1_md5_expected)}
    if not single_index:
        file_info["index2"] = _empty_file_info(index2_file, index2_md5_expected)

    # Index length composition (dict per index file), filled while parsing.
    index_lengths = {"index1": Counter()}
    if not single_index:
        index_lengths["index2"] = Counter()

    # Track mismatches
    mismatches = []

    if verbose:
        print("Calculating MD5 checksums...")
        print(f"Index1 MD5: {file_info['index1']['md5']}")
        if not single_index:
            print(f"Index2 MD5: {file_info['index2']['md5']}")
        print("Parsing FASTQ files...")

    # Open and parse the FASTQ files
    try:
        # Single-index runs read index1 only; index2_handle stays None and the
        # loop below yields (record1, None) pairs.
        index2_handle = None if single_index else gzip.open(index2_file, "rt")
        with gzip.open(index1_file, "rt") as index1_handle:
            records1 = SeqIO.parse(index1_handle, "fastq")
            if index2_handle is None:
                record_pairs = ((record1, None) for record1 in records1)
            else:
                record_pairs = zip(records1, SeqIO.parse(index2_handle, "fastq"))
            for record1, record2 in record_pairs:
                seq_id1 = record1.id  # Extract sequence ID from index1
                seq_id2 = record2.id if record2 is not None else ""
                
                # Increment read count
                file_info["index1"]["read_count"] += 1
                if record2 is not None:
                    file_info["index2"]["read_count"] += 1
                
                # Display progress every 10,000 reads
                if file_info["index1"]["read_count"] % 10000 == 0:
                    print(f"Processed {file_info['index1']['read_count']:,} reads...", flush=True)
                
                # Extract header information
                header_info1 = extract_fastq_header_info(record1)
                header_info2 = (extract_fastq_header_info(record2)
                                if record2 is not None else None)
                
                file_info["index1"]["header_info"].append(header_info1)
                if header_info2 is not None:
                    file_info["index2"]["header_info"].append(header_info2)
                
                # Update sets with header information
                file_info["index1"]["instrument_codes"].add(header_info1["instrument_code"])
                file_info["index1"]["instrument_names"].add(header_info1["instrument_name"])
                file_info["index1"]["runs"].add(header_info1["run"])
                file_info["index1"]["flowcells"].add(header_info1["flowcell"])
                file_info["index1"]["lanes"].add(header_info1["lane"])
                
                if header_info2 is not None:
                    file_info["index2"]["instrument_codes"].add(header_info2["instrument_code"])
                    file_info["index2"]["instrument_names"].add(header_info2["instrument_name"])
                    file_info["index2"]["runs"].add(header_info2["run"])
                    file_info["index2"]["flowcells"].add(header_info2["flowcell"])
                    file_info["index2"]["lanes"].add(header_info2["lane"])
                
                # Check for mismatches in sequence IDs
                if record2 is not None and seq_id1 != seq_id2:
                    mismatches.append(f"Sequence ID mismatch: {seq_id1} != {seq_id2}")
                    continue  # Skip this read pair due to mismatch
                
                # Check for mismatches in header information
                if header_info2 is not None and header_info1 != header_info2:
                    mismatch_details = []
                    for key in header_info1:
                        if header_info1[key] != header_info2[key]:
                            mismatch_details.append(f"{key}: {header_info1[key]} != {header_info2[key]}")
                    
                    mismatches.append(f"Header mismatch for {seq_id1}: {', '.join(mismatch_details)}")
                
                barcode1 = str(record1.seq)  # Extract barcode sequence from index1
                barcode2 = str(record2.seq) if record2 is not None else ""
                
                # Create full barcode (barcode1 + reverse of barcode2)
                barcode2_reverse = barcode2[::-1]
                full_barcode = barcode1 + barcode2_reverse
                
                # Store barcode sequences with sequence ID as key
                barcodes[seq_id1] = (barcode1, barcode2, full_barcode)

                # Index length composition (e.g. 12 nt Golay barcode reads).
                index_lengths["index1"][len(barcode1)] += 1
                if record2 is not None:
                    index_lengths["index2"][len(barcode2)] += 1
    except FileNotFoundError as e:
        print(f"Error: File not found - {str(e)}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error processing FASTQ files: {str(e)}", file=sys.stderr)
        sys.exit(1)

    if verbose:
        print(f"Processed {file_info['index1']['read_count']} reads from index1")
        if not single_index:
            print(f"Processed {file_info['index2']['read_count']} reads from index2")
        print(f"Found {len(mismatches)} mismatches between files")

    # Save barcode sequences to a text file
    barcodes_output_file = os.path.join(data_dir, "barcodes.txt")
    with open(barcodes_output_file, "w") as out:
        out.write("Seq_ID\tBarcode1\tBarcode2\tFullBarcode\n")
        for seq_id, (barcode1, barcode2, full_barcode) in barcodes.items():
            out.write(f"{seq_id}\t{barcode1}\t{barcode2}\t{full_barcode}\n")
    
    if verbose:
        print(f"Saved barcode sequences to {barcodes_output_file}")

    # Create summary report
    summary_output_file = os.path.join(data_dir, "sequencing_summary.txt")

    # Single-index runs have no second index file: report the I1 metrics and the
    # barcode composition instead of the dual-index congruency table.
    if single_index:
        write_single_index_summary(summary_output_file, file_info, index_lengths, barcodes)
        if verbose:
            print(f"Saved sequencing summary to {summary_output_file}")
        return barcodes_output_file, summary_output_file

    with open(summary_output_file, "w") as out:
        out.write("# Sequencing Run Information and Congruency Check\n\n")
        
        # Write file metrics
        out.write("## File Metrics\n\n")
        out.write(f"{'Parameter':<20}{'Index1':<40}{'Index2':<40}{'Match?':<10}\n")
        out.write(f"{'-'*20}{'-'*40}{'-'*40}{'-'*10}\n")
        
        # Filename
        out.write(f"{'Filename':<20}{file_info['index1']['filename']:<40}{file_info['index2']['filename']:<40}{'N/A':<10}\n")
        
        # MD5 - verify each index file against its own expected checksum.
        # A file "matches itself" when its computed MD5 equals the expected
        # value configured for that same file. When no expected value is
        # configured the check is skipped and reported as N/A.
        index1_expected = file_info['index1']['expected_md5']
        index2_expected = file_info['index2']['expected_md5']
        index1_ok = (not index1_expected) or (file_info['index1']['md5'].lower() == index1_expected)
        index2_ok = (not index2_expected) or (file_info['index2']['md5'].lower() == index2_expected)
        md5_checked = bool(index1_expected) or bool(index2_expected)
        md5_match = index1_ok and index2_ok

        if not md5_checked:
            md5_label = "N/A"
        elif md5_match:
            md5_label = "Yes"
        else:
            failed = []
            if not index1_ok:
                failed.append("Index1")
            if not index2_ok:
                failed.append("Index2")
            md5_label = f"No ({', '.join(failed)})"

        out.write(f"{'MD5':<20}{file_info['index1']['md5']:<40}{file_info['index2']['md5']:<40}{md5_label:<10}\n")
        
        # Read count
        read_count_match = file_info['index1']['read_count'] == file_info['index2']['read_count']
        out.write(f"{'Read count':<20}{file_info['index1']['read_count']:<40}{file_info['index2']['read_count']:<40}{'Yes' if read_count_match else 'No':<10}\n")
        
        # Instrument Codes
        instrument_codes_match = file_info['index1']['instrument_codes'] == file_info['index2']['instrument_codes']
        out.write(f"{'Instrument Codes':<20}{', '.join(file_info['index1']['instrument_codes']):<40}{', '.join(file_info['index2']['instrument_codes']):<40}{'Yes' if instrument_codes_match else 'No':<10}\n")
        
        # Instrument Names
        instrument_names_match = file_info['index1']['instrument_names'] == file_info['index2']['instrument_names']
        out.write(f"{'Instrument Names':<20}{', '.join(file_info['index1']['instrument_names']):<40}{', '.join(file_info['index2']['instrument_names']):<40}{'Yes' if instrument_names_match else 'No':<10}\n")
        
        # Runs (Changed to Run_IDs to avoid confusion)
        runs_match = file_info['index1']['runs'] == file_info['index2']['runs']
        run_ids_str1 = ', '.join([f'"{run_id}"' for run_id in file_info['index1']['runs']])
        run_ids_str2 = ', '.join([f'"{run_id}"' for run_id in file_info['index2']['runs']])
        out.write(f"{'Run_IDs':<20}{run_ids_str1:<40}{run_ids_str2:<40}{'Yes' if runs_match else 'No':<10}\n")
        
        # Flowcells
        flowcells_match = file_info['index1']['flowcells'] == file_info['index2']['flowcells']
        out.write(f"{'Flowcells':<20}{', '.join(file_info['index1']['flowcells']):<40}{', '.join(file_info['index2']['flowcells']):<40}{'Yes' if flowcells_match else 'No':<10}\n")
        
        # Lanes
        lanes_match = file_info['index1']['lanes'] == file_info['index2']['lanes']
        out.write(f"{'Lanes':<20}{', '.join(file_info['index1']['lanes']):<40}{', '.join(file_info['index2']['lanes']):<40}{'Yes' if lanes_match else 'No':<10}\n")
        
        # Mismatches
        if mismatches:
            out.write("\n## Detected Mismatches\n\n")
            for i, mismatch in enumerate(mismatches[:100], 1):  # Limit to first 100 mismatches
                out.write(f"{i}. {mismatch}\n")
            
            if len(mismatches) > 100:
                out.write(f"\n... and {len(mismatches) - 100} more mismatches\n")
        else:
            out.write("\n## No mismatches detected\n")
        
        # Summary
        out.write("\n## Summary\n\n")
        
        all_match = all([read_count_match, instrument_codes_match, instrument_names_match, 
                         runs_match, flowcells_match, lanes_match])
        # Only fold the MD5 check into the overall verdict when it was actually
        # performed (i.e. at least one expected checksum was configured).
        if md5_checked and not md5_match:
            all_match = False
        out.write(f"Overall congruency: {'Good' if all_match and not mismatches else 'Issues detected'}\n")
        
        if not all_match or mismatches:
            out.write("Detected issues:\n")
            if md5_checked and not md5_match:
                out.write(f"- MD5 checksum mismatch ({md5_label})\n")
            if not read_count_match:
                out.write("- Read count mismatch between index files\n")
            if not instrument_codes_match:
                out.write("- Instrument code mismatch between index files\n")
            if not instrument_names_match:
                out.write("- Instrument name mismatch between index files\n")
            if not runs_match:
                out.write("- Run ID mismatch between index files\n")
            if not flowcells_match:
                out.write("- Flowcell ID mismatch between index files\n")
            if not lanes_match:
                out.write("- Lane mismatch between index files\n")
            if mismatches:
                out.write(f"- {len(mismatches)} sequence ID or header information mismatches detected\n")
    
    if verbose:
        print(f"Saved sequencing summary to {summary_output_file}")
    
    return barcodes_output_file, summary_output_file

def main():
    """Main function to execute the script."""
    # Parse command line arguments
    args = parse_arguments()
    
    # Process FASTQ files
    barcodes_file, summary_file = process_fastq_files(
        args.index1, 
        args.index2, 
        args.outdir,
        args.verbose,
        index1_md5_expected=args.index1_md5,
        index2_md5_expected=args.index2_md5,
    )
    
    print(f"Analysis complete!")
    print(f"Barcode file: {barcodes_file}")
    print(f"Summary file: {summary_file}")

if __name__ == "__main__":
    main()
