# =============================================================================
# common.smk - shared helpers, wildcards, and constants
# =============================================================================
# Included by the main Snakefile. Defines the sample list and helper functions
# used across all rule files. Contains NO rules that do work.
# =============================================================================

import os
import pandas as pd

# -----------------------------------------------------------------------------
# Load the samplesheet
# -----------------------------------------------------------------------------
samplesheet_path = os.path.join(ROOT_DIR, config["samplesheet"])
samples = pd.read_csv(samplesheet_path, sep="\t", dtype=str).fillna("")

# The canonical list of sample IDs (wildcard {sample}).
SAMPLES = samples["sample_id"].tolist()

# FASTQ filename suffixes.
r1_suffix = config["sequencing"]["read_suffix"]["r1"]
r2_suffix = config["sequencing"]["read_suffix"]["r2"]


# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------
def _sample_fastq(sample_id, column, suffix):
    """Resolve a sample's FASTQ path.

    Prefer the explicit path recorded in the samplesheet column (``r1``/``r2``);
    fall back to ``raw_dir/<sample_id><suffix>`` when the column is empty.
    """
    row = samples.loc[samples["sample_id"] == sample_id]
    if not row.empty:
        value = str(row.iloc[0].get(column, "") or "").strip()
        if value:
            return os.path.join(ROOT_DIR, value)
    return os.path.join(ROOT_DIR, config["raw_dir"], sample_id + suffix)


def sample_r1(wildcards):
    """Return the path to the raw R1 FASTQ for a sample."""
    return _sample_fastq(wildcards.sample, "r1", r1_suffix)


def sample_r2(wildcards):
    """Return the path to the raw R2 FASTQ for a sample."""
    return _sample_fastq(wildcards.sample, "r2", r2_suffix)


def get_sample_row(sample_id):
    """Return the metadata row for a given sample ID as a dict."""
    row = samples.loc[samples["sample_id"] == sample_id]
    if row.empty:
        raise ValueError(f"Sample '{sample_id}' not found in samplesheet.")
    return row.iloc[0].to_dict()


# -----------------------------------------------------------------------------
# Directory constants (relative to repository root)
# -----------------------------------------------------------------------------
RESULTS_DIR = os.path.join(ROOT_DIR, "results")
LOGS_DIR = os.path.join(ROOT_DIR, "logs")
RESOURCES_DIR = os.path.join(ROOT_DIR, "resources")