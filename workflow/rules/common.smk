# =============================================================================
# common.smk - shared helpers, wildcards, and constants
# =============================================================================
# Included by the main Snakefile. Defines the sample list and helper functions
# used across all rule files. Contains NO rules that do work.
# =============================================================================

import os
import pandas as pd

# -----------------------------------------------------------------------------
# Directory constants (relative to repository root)
# -----------------------------------------------------------------------------
RESULTS_DIR = os.path.join(ROOT_DIR, "results")
LOGS_DIR = os.path.join(ROOT_DIR, "logs")
RESOURCES_DIR = os.path.join(ROOT_DIR, "resources")

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
# Demultiplexing
# -----------------------------------------------------------------------------
# When `demultiplex.enabled` is true the raw pool (e.g. an Undetermined_* run) is
# split into per-sample FASTQs by the demultiplex stage, and every downstream
# stage reads those instead of the samplesheet/raw_dir paths.
demux_cfg = config.get("demultiplex", {})
DEMULTIPLEX_ENABLED = bool(demux_cfg.get("enabled", False))
DEMUX_DIR = os.path.join(RESULTS_DIR, "demultiplex")
DEMUX_SUFFIX = {"r1": "_R1.fastq.gz", "r2": "_R2.fastq.gz"}


def demux_sample_fastq(sample_id, read="r1"):
    """Return the demultiplexed FASTQ path for a sample and read (r1/r2)."""
    return os.path.join(DEMUX_DIR, sample_id + DEMUX_SUFFIX[read])


# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------
def _sample_fastq(sample_id, column, suffix, read="r1"):
    """Resolve a sample's FASTQ path.

    Precedence:
      1. demultiplexed output, when `demultiplex.enabled` is true
      2. the explicit path recorded in the samplesheet column (``r1``/``r2``)
      3. ``raw_dir/<sample_id><suffix>``
    """
    if DEMULTIPLEX_ENABLED:
        return demux_sample_fastq(sample_id, read)
    row = samples.loc[samples["sample_id"] == sample_id]
    if not row.empty:
        value = str(row.iloc[0].get(column, "") or "").strip()
        if value:
            return os.path.join(ROOT_DIR, value)
    return os.path.join(ROOT_DIR, config["raw_dir"], sample_id + suffix)


def sample_r1(wildcards):
    """Return the path to the R1 FASTQ for a sample (demultiplexed when enabled)."""
    return _sample_fastq(wildcards.sample, "r1", r1_suffix, read="r1")


def sample_r2(wildcards):
    """Return the path to the R2 FASTQ for a sample (demultiplexed when enabled)."""
    return _sample_fastq(wildcards.sample, "r2", r2_suffix, read="r2")


def get_sample_row(sample_id):
    """Return the metadata row for a given sample ID as a dict."""
    row = samples.loc[samples["sample_id"] == sample_id]
    if row.empty:
        raise ValueError(f"Sample '{sample_id}' not found in samplesheet.")
    return row.iloc[0].to_dict()