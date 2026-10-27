#!/usr/bin/env python3
"""
Export a barcode table (xlsx/csv/tsv) into the two files the pipeline needs:

  1. a QIIME 2 metadata file (`#SampleID` + barcode column + optional extra
     columns) used by `qiime demux emp-paired`
  2. a pipeline samplesheet (`config/samples.<run>.tsv`) with one row per sample
     that has a barcode

Both outputs are plain text and tracked in git: the raw barcode table lives under
`data/raw/` (gitignored), so these derived files are the reproducible, auditable
record of which sample carries which index. The source file's md5 and sheet name
are reported so the derivation can be traced back.

Rows without a barcode are dropped (they cannot be demultiplexed); duplicate
barcodes or barcodes with unexpected characters/length are reported as errors.

Usage:
    python export_barcodes.py \
        --barcode_table data/raw/20240209_16S/barcodes.xlsx \
        --sheet all --sample_column sample-id --barcode_column barcode-sequence \
        --out_metadata config/metadata.20240209_16S.tsv \
        --out_samplesheet config/samples.20240209_16S.tsv \
        --barcode_length 12
"""

import argparse
import hashlib
import os
import sys

import pandas as pd

SAMPLESHEET_COLUMNS = [
    "sample_id",
    "batch",
    "group",
    "sample_type",
    "subject",
    "r1",
    "r2",
    "barcode",
    "forward_primer",
    "reverse_primer",
]


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Export a barcode table to QIIME 2 metadata + pipeline samplesheet."
    )
    parser.add_argument("--barcode_table", required=True,
                        help="Barcode table (xlsx/csv/tsv)")
    parser.add_argument("--sheet", default=0,
                        help="Sheet name for xlsx input (default: first sheet)")
    parser.add_argument("--sample_column", default="sample-id",
                        help="Column holding sample identifiers")
    parser.add_argument("--barcode_column", default="barcode-sequence",
                        help="Column holding barcode sequences")
    parser.add_argument("--out_metadata", required=True,
                        help="QIIME 2 metadata file to write")
    parser.add_argument("--out_samplesheet", required=True,
                        help="Pipeline samplesheet to write")
    parser.add_argument("--barcode_length", type=int, default=0,
                        help="Expected barcode length (0 = accept any length)")
    parser.add_argument("--extra_columns", default="",
                        help="Comma-separated extra table columns to carry into the "
                             "metadata file (e.g. plate,well)")
    parser.add_argument("--barcode_field", default="BarcodeSequence",
                        help="Barcode column name in the metadata file")
    return parser.parse_args()


def md5sum(path):
    """Return the md5 checksum of a file (recorded for provenance)."""
    digest = hashlib.md5()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_table(path, sheet, sample_column, barcode_column, extra_columns):
    """Read the barcode table and return (dataframe, source description)."""
    lowered = path.lower()
    if lowered.endswith((".xlsx", ".xlsm", ".xls")):
        sheets = pd.read_excel(path, sheet_name=None, dtype=str)
        if sheet in (None, "", 0):
            name = list(sheets)[0]
        else:
            name = sheet
        if name not in sheets:
            sys.exit(f"Error: sheet '{name}' not in {path} (available: {', '.join(sheets)})")
        frame = sheets[name]
        source = f"{path} (sheet '{name}')"
    elif lowered.endswith((".tsv", ".txt")):
        frame = pd.read_csv(path, sep="\t", dtype=str, comment="#")
        source = path
    else:
        frame = pd.read_csv(path, sep=",", dtype=str, comment="#")
        source = path

    wanted = [sample_column, barcode_column]
    wanted += [c.strip() for c in extra_columns.split(",") if c.strip()]
    missing = [c for c in wanted if c not in frame.columns]
    if missing:
        sys.exit(
            f"Error: column(s) not found: {', '.join(missing)}\n"
            f"Available columns: {', '.join(map(str, frame.columns))}"
        )

    table = frame[wanted].copy()
    table[sample_column] = table[sample_column].fillna("").str.strip()
    table[barcode_column] = table[barcode_column].fillna("").str.strip().str.upper()
    return table, source


def main():
    args = parse_arguments()
    if not os.path.exists(args.barcode_table):
        sys.exit(f"Error: barcode table not found: {args.barcode_table}")

    table, source = read_table(
        args.barcode_table, args.sheet, args.sample_column, args.barcode_column,
        args.extra_columns
    )
    print(f"Barcode table source: {source}")
    print(f"Source md5: {md5sum(args.barcode_table)}")

    total_rows = len(table)
    table = table[table[args.sample_column] != ""]
    nameless = total_rows - len(table)

    no_barcode = int((table[args.barcode_column] == "").sum())
    usable = table[table[args.barcode_column] != ""].copy()

    # Barcode sanity checks: allowed characters, expected length, uniqueness.
    bad_chars = usable[~usable[args.barcode_column].str.fullmatch(r"[ACGTN]+")]
    if len(bad_chars):
        sys.exit("Error: barcodes with unexpected characters: "
                 + ", ".join(bad_chars[args.barcode_column].head(5)))
    if args.barcode_length:
        wrong_len = usable[usable[args.barcode_column].str.len() != args.barcode_length]
        if len(wrong_len):
            sys.exit(f"Error: barcodes not {args.barcode_length} bp long: "
                     + ", ".join(wrong_len[args.barcode_column].head(5)))
    if usable[args.barcode_column].duplicated().any():
        dups = usable[args.barcode_column][usable[args.barcode_column].duplicated()].unique()
        sys.exit("Error: barcode(s) assigned to more than one sample: " + ", ".join(dups[:5]))
    duplicated_samples = int(usable[args.sample_column].duplicated().sum())
    if duplicated_samples:
        # Duplicated ids cannot be written to QIIME 2 metadata (ids must be
        # unique) or demultiplexed reliably: keep the first occurrence and report
        # the affected ids so the barcode table can be curated upstream.
        dup_ids = usable[args.sample_column][usable[args.sample_column].duplicated()].unique()
        print(f"WARNING: {duplicated_samples} duplicated sample id(s) - keeping the first "
              f"occurrence of: {', '.join(dup_ids[:10])}")
        usable = usable.drop_duplicates(subset=[args.sample_column], keep="first")
        print(f"  usable samples after dropping duplicates: {len(usable)}")

    print(f"Rows: {total_rows}")
    print(f"  without a sample id: {nameless}")
    print(f"  without a barcode (cannot be demultiplexed): {no_barcode}  <- dropped")
    print(f"  usable samples: {len(usable)}")

    for path in (args.out_metadata, args.out_samplesheet):
        parent = os.path.dirname(path)
        if parent:
            os.makedirs(parent, exist_ok=True)

    # ------------------------------------------------------------------
    # QIIME 2 metadata file (used by qiime demux emp-paired)
    # ------------------------------------------------------------------
    metadata = pd.DataFrame({
        "#SampleID": usable[args.sample_column].values,
        args.barcode_field: usable[args.barcode_column].values,
    })
    for column in [c.strip() for c in args.extra_columns.split(",") if c.strip()]:
        metadata[column] = usable[column].values
    metadata.to_csv(args.out_metadata, sep="\t", index=False)
    print(f"Wrote QIIME 2 metadata ({len(metadata)} samples): {args.out_metadata}")

    # ------------------------------------------------------------------
    # Pipeline samplesheet
    # ------------------------------------------------------------------
    samplesheet = pd.DataFrame({column: "" for column in SAMPLESHEET_COLUMNS},
                               index=range(len(usable)))
    samplesheet["sample_id"] = usable[args.sample_column].values
    samplesheet["barcode"] = usable[args.barcode_column].values
    # r1/r2 stay empty on purpose: workflow/rules/common.smk resolves them to the
    # demultiplexed FASTQs when demultiplexing is enabled, or to raw_dir otherwise.
    samplesheet.to_csv(args.out_samplesheet, sep="\t", index=False)
    print(f"Wrote pipeline samplesheet ({len(samplesheet)} samples): {args.out_samplesheet}")


if __name__ == "__main__":
    main()
