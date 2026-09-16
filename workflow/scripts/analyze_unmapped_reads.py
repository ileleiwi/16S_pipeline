#!/usr/bin/env python3

import gzip
import argparse
from collections import defaultdict, Counter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import BiopythonParserWarning
import warnings

warnings.simplefilter("ignore", BiopythonParserWarning)

def read_fastq_safe(file):
    reads = {}
    skipped = 0
    skip_reasons = {
        "no_header": 0,
        "no_plus": 0,
        "len_mismatch": 0,
        "exception": 0
    }
    
    # Determine if file is gzipped based on extension
    is_gzipped = file.endswith('.gz')
    
    # Use appropriate open function
    open_func = gzip.open if is_gzipped else open
    open_mode = "rt" if is_gzipped else "r"
    
    with open_func(file, open_mode) as f:
        while True:
            try:
                header = f.readline().strip()
                if not header:  # EOF check
                    break
                    
                seq = f.readline().strip()
                plus = f.readline().strip()
                qual = f.readline().strip()
                
                if not header.startswith("@"):
                    skip_reasons["no_header"] += 1
                    skipped += 1
                    continue
                    
                if not plus.startswith("+"):
                    skip_reasons["no_plus"] += 1
                    skipped += 1
                    continue
                    
                if len(seq) != len(qual):
                    skip_reasons["len_mismatch"] += 1
                    skipped += 1
                    continue
                    
                read_id = header.split()[0][1:]
                reads[read_id] = seq
                
            except Exception as e:
                print(f"Exception: {str(e)}")
                skip_reasons["exception"] += 1
                skipped += 1
                continue
                
    return reads, skipped, skip_reasons

def load_contaminants(fasta_file):
    contaminants = {}
    with open(fasta_file, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            contaminants[record.id] = str(record.seq)
    return contaminants

def match_contaminants(reads, contaminant_seqs):
    match_counts = Counter()
    for read_seq in reads.values():
        for cid, cseq in contaminant_seqs.items():
            if cseq in read_seq:
                match_counts[cid] += 1
    return match_counts

def main():
    parser = argparse.ArgumentParser(description="Analyze unmapped reads from demux step.")
    parser.add_argument("--forward", required=True, help="Undetermined_R1.fastq.gz")
    parser.add_argument("--reverse", required=True, help="Undetermined_R2.fastq.gz")
    parser.add_argument("--contaminants", required=True, help="FASTA file of contaminant sequences")
    args = parser.parse_args()

    print("🔍 Parsing forward and reverse reads safely...")
    forward_reads, skipped_f, skip_reasons_f = read_fastq_safe(args.forward)
    reverse_reads, skipped_r, skip_reasons_r = read_fastq_safe(args.reverse)

    forward_ids = set(forward_reads.keys())
    reverse_ids = set(reverse_reads.keys())

    total_forward = len(forward_ids)
    total_reverse = len(reverse_ids)
    paired = forward_ids & reverse_ids  
    forward_only = forward_ids - reverse_ids
    reverse_only = reverse_ids - forward_ids

    print(f"Total forward reads: {total_forward}")
    print(f"Total reverse reads: {total_reverse}")
    print(f"Skipped forward records: {skipped_f}")
    print(f"Skipped reverse records: {skipped_r}")
    print(f"Forward skip reasons: {skip_reasons_f}")
    print(f"Reverse skip reasons: {skip_reasons_r}")
    print(f"Paired (orphaned) reads: {len(paired)} ({len(paired)/(total_forward + total_reverse)*100:.2f}%)")
    print(f"Forward-only reads: {len(forward_only)} ({len(forward_only)/total_forward*100:.2f}%)")
    print(f"Reverse-only reads: {len(reverse_only)} ({len(reverse_only)/total_reverse*100:.2f}%)")

    print("📦 Loading contaminant sequences...")
    contaminants = load_contaminants(args.contaminants)

    print("🧪 Searching forward reads for contaminants...")
    forward_matches = match_contaminants(forward_reads, contaminants)
    print("🧪 Searching reverse reads for contaminants...")
    reverse_matches = match_contaminants(reverse_reads, contaminants)

    total_reads = total_forward + total_reverse
    print("\n🔬 Contaminant match summary:")
    combined_matches = forward_matches + reverse_matches
    for cid, count in combined_matches.most_common():
        pct = (count / total_reads) * 100
        print(f"{cid}: {count} reads ({pct:.2f}%)")

if __name__ == "__main__":
    main()
