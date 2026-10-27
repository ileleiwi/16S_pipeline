#!/usr/bin/env python3
"""
Normalise a QIIME 2 demultiplex export into the pipeline's FASTQ naming.

`qiime tools export` writes one FASTQ pair per sample in Casava 1.8
single-lane-per-sample format:

    <sample>_S<n>_L001_R1_001.fastq.gz
    <sample>_S<n>_L001_R2_001.fastq.gz

Downstream rules expect `results/demultiplex/<sample>_R1.fastq.gz` /
`_R2.fastq.gz` (see workflow/rules/common.smk), so this script links each sample's
exported read pair to that name, counts the reads per sample and writes a stats
table. Samples that received no reads are listed explicitly: with an incomplete
barcode list those reads are dropped by the demultiplexer, so the count belongs in
the report rather than being silently lost.

Usage:
    python normalize_demux_export.py \
        --export_dir results/demultiplex/exported \
        --outdir results/demultiplex \
        --stats results/demultiplex/demux_stats.tsv \
        --samples sampleA sampleB ...
"""

import argparse
import glob
import gzip
import os
import re
import sys

# QIIME 2 appends "<num>_L001" (or "S<num>_L001") to the sample id on export:
#   4072_0523_53.015_d071_rnalater_18_L001_R1_001.fastq.gz
# so the read suffix is parsed from the end and the sample id is resolved against
# the expected sample list instead of being reconstructed from the filename.
READ_SUFFIX_RE = re.compile(r"^(?P<stem>.+)_R(?P<num>[12])(?:_001)?\.fastq(?:\.gz)?$")


def parse_arguments():
    parser = argparse.ArgumentParser(description="Normalise a QIIME 2 demux export.")
    parser.add_argument("--export_dir", required=True,
                        help="Directory produced by `qiime tools export`")
    parser.add_argument("--outdir", required=True,
                        help="Directory for <sample>_R{1,2}.fastq.gz outputs")
    parser.add_argument("--stats", required=True,
                        help="TSV with per-sample read counts")
    parser.add_argument("--samples", nargs="+", required=True,
                        help="Expected sample ids (from the samplesheet)")
    parser.add_argument("--no_count", action="store_true",
                        help="Skip read counting (faster; counts written as NA)")
    return parser.parse_args()


def find_reads(export_dir, expected_samples):
    """Map (sample, read) -> path for the FASTQs in a QIIME 2 export directory.

    `qiime tools export` writes one FASTQ pair per sample, adding a sample number
    and lane to the sample id (e.g. `<sample>_18_L001_R1_001.fastq.gz`). Sample ids
    are therefore resolved by matching the *longest* expected sample id that the
    filename starts with, which also copes with ids that themselves contain
    underscores or dots.
    """
    found = {}
    ordered = sorted(expected_samples, key=len, reverse=True)
    pattern = os.path.join(export_dir, "**", "*.fastq*")
    for path in sorted(glob.glob(pattern, recursive=True)):
        name = os.path.basename(path)
        match = READ_SUFFIX_RE.match(name)
        if not match:
            continue
        read = "r1" if match.group("num") == "1" else "r2"
        for sample in ordered:
            if name.startswith(sample + "_"):
                found[(sample, read)] = path
                break
    return found


def count_reads(path):
    """Count FASTQ records by streaming the (optionally gzipped) file."""
    opener = gzip.open if path.endswith(".gz") else open
    mode = "rb" if path.endswith(".gz") else "rb"
    lines = 0
    with opener(path, mode) as handle:
        for chunk in iter(lambda: handle.read(1 << 22), b""):
            lines += chunk.count(b"\n")
    return lines // 4


def main():
    args = parse_arguments()
    if not os.path.isdir(args.export_dir):
        sys.exit(f"Error: export directory not found: {args.export_dir}")

    os.makedirs(args.outdir, exist_ok=True)
    stats_parent = os.path.dirname(args.stats)
    if stats_parent:
        os.makedirs(stats_parent, exist_ok=True)

    found = find_reads(args.export_dir, args.samples)
    print(f"Exported FASTQ files discovered: {len(found)}")

    rows = []
    missing = []
    zero_reads = []
    for sample in args.samples:
        row = {"sample": sample}
        for read in ("r1", "r2"):
            source = found.get((sample, read))
            target = os.path.join(args.outdir, f"{sample}_{read.upper()}.fastq.gz")
            row[f"{read}_reads"] = "NA"
            if os.path.islink(target) or os.path.exists(target):
                os.remove(target)

            if source is None:
                # No reads were exported for this sample: its barcode was not found
                # in the barcode read (an incomplete barcode list looks exactly like
                # this). The file is declared as a workflow output, so write an
                # empty gzip FASTQ placeholder and record the status in the stats
                # table instead of failing the job; downstream stages must skip
                # samples whose status is "missing".
                row[f"{read}_status"] = "missing"
                with gzip.open(target, "wb"):
                    pass
                continue

            os.symlink(os.path.abspath(source), target)
            row[f"{read}_status"] = "ok"
            # Preserve the gzip flag of the exported file.
            if not source.endswith(".gz"):
                os.remove(target)
                with open(source, "rb") as src, gzip.open(target, "wb") as dst:
                    dst.write(src.read())
            if not args.no_count:
                row[f"{read}_reads"] = count_reads(source)
        if row["r1_status"] == "missing" and row["r2_status"] == "missing":
            missing.append(sample)
        elif row["r1_reads"] not in ("NA", 0) and row["r1_reads"] == 0:
            zero_reads.append(sample)
        rows.append(row)

    with open(args.stats, "w") as handle:
        handle.write("sample\tr1_status\tr1_reads\tr2_status\tr2_reads\n")
        for row in rows:
            handle.write(f"{row['sample']}\t{row['r1_status']}\t{row['r1_reads']}\t"
                         f"{row['r2_status']}\t{row['r2_reads']}\n")

    total_r1 = sum(r["r1_reads"] for r in rows if isinstance(r["r1_reads"], int))
    print(f"Wrote stats for {len(rows)} samples: {args.stats}")
    print(f"Total demultiplexed R1 read pairs: {total_r1:,}")
    print(f"Samples with no exported reads (not demultiplexable): {len(missing)}")
    if missing:
        print("  " + ", ".join(missing[:10]) + (" ..." if len(missing) > 10 else ""))
    if zero_reads:
        print(f"Samples with zero reads: {len(zero_reads)}")


if __name__ == "__main__":
    main()
