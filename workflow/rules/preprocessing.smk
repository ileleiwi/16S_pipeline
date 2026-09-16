# =============================================================================
# preprocessing.smk - pre-QC input and index diagnostics
# =============================================================================
# This stage runs AHEAD of quality control. It validates the raw sequencing
# inputs and the library/index configuration so that problems are caught before
# any filtering or trimming happens. Every rule here is diagnostic: it reports
# issues (barcode format, index congruency, contaminant carry-over, malformed
# FASTQ) rather than modifying the data.
#
# Scripts wrapped here live in workflow/scripts/:
#   - index_fastq_check.py     -> index congruency + barcode extraction
#   - barcode-diagnostic.py    -> optimal barcode format for demultiplexing
#   - analyze_unmapped_reads.py-> contaminant screen of undetermined reads
#   - diagnostic.py            -> per-file FASTQ read count / format validation
# =============================================================================

import os

# -----------------------------------------------------------------------------
# Paths and configuration
# -----------------------------------------------------------------------------
SCRIPTS_DIR = os.path.join(str(WORKFLOW_DIR), "scripts")
DIAGNOSTICS_ENV = os.path.join(str(WORKFLOW_DIR), "envs", "diagnostics.yaml")

# The diagnostics block is optional; fall back to disabled if it is missing.
diag_cfg = config.get("diagnostics", {})
DIAGNOSTICS_ENABLED = diag_cfg.get("enabled", False)
barcode_cfg = diag_cfg.get("barcode", {})

# Resolve run-level inputs relative to the repository root.
INDEX1 = os.path.join(ROOT_DIR, diag_cfg.get("index1", ""))
INDEX2 = os.path.join(ROOT_DIR, diag_cfg.get("index2", ""))
MAPPING_FILE = os.path.join(ROOT_DIR, diag_cfg.get("mapping_file", ""))
UNMAPPED_R1 = os.path.join(ROOT_DIR, diag_cfg.get("unmapped_r1", ""))
UNMAPPED_R2 = os.path.join(ROOT_DIR, diag_cfg.get("unmapped_r2", ""))
CONTAMINANTS = os.path.join(ROOT_DIR, diag_cfg.get("contaminants", ""))

# Output locations for this stage.
PREPROC_DIR = os.path.join(RESULTS_DIR, "preprocessing")
INDEX_CHECK_DIR = os.path.join(PREPROC_DIR, "index_check")
BARCODE_DIAG_DIR = os.path.join(PREPROC_DIR, "barcode_diagnostic")
UNMAPPED_DIR = os.path.join(PREPROC_DIR, "unmapped_reads")
FASTQ_CHECK_DIR = os.path.join(PREPROC_DIR, "fastq_check")

MEM_MB = config["resources"]["mem_mb_default"]


def preprocessing_targets():
    """Return the outputs produced by the pre-QC diagnostics stage.

    Returns an empty list when diagnostics are disabled so the stage can be
    switched off from config without editing the Snakefile.
    """
    if not DIAGNOSTICS_ENABLED:
        return []
    targets = [
        os.path.join(INDEX_CHECK_DIR, "data", "sequencing_summary.txt"),
        os.path.join(INDEX_CHECK_DIR, "data", "barcodes.txt"),
        os.path.join(BARCODE_DIAG_DIR, "barcode_format_report.txt"),
    ]
    # The contaminant screen needs a contaminant reference file; only request
    # its output when the reference (and undetermined reads) are available.
    if os.path.exists(CONTAMINANTS) and os.path.exists(UNMAPPED_R1):
        targets.append(os.path.join(UNMAPPED_DIR, "contaminant_report.txt"))
    targets += [
        os.path.join(FASTQ_CHECK_DIR, f"{sample}_R1_validation.txt")
        for sample in SAMPLES
    ]
    return targets


# -----------------------------------------------------------------------------
# Rules
# -----------------------------------------------------------------------------
rule check_index_fastq:
    """Check index-file congruency and extract per-read barcodes (I1/I2)."""
    input:
        index1=INDEX1,
        index2=INDEX2,
    output:
        summary=os.path.join(INDEX_CHECK_DIR, "data", "sequencing_summary.txt"),
        barcodes=os.path.join(INDEX_CHECK_DIR, "data", "barcodes.txt"),
    log:
        os.path.join(LOGS_DIR, "preprocessing", "index_check.log"),
    params:
        script=os.path.join(SCRIPTS_DIR, "index_fastq_check.py"),
        outdir=INDEX_CHECK_DIR,
    conda:
        DIAGNOSTICS_ENV
    threads: 1
    resources:
        mem_mb=MEM_MB,
    shell:
        r"""
        python {params.script} \
            --index1 {input.index1} \
            --index2 {input.index2} \
            --outdir {params.outdir} \
            --verbose > {log} 2>&1
        """


rule diagnose_barcode_formats:
    """Determine the optimal barcode format for demultiplexing."""
    input:
        index1=INDEX1,
        index2=INDEX2,
        mapping=MAPPING_FILE,
    output:
        report=os.path.join(BARCODE_DIAG_DIR, "barcode_format_report.txt"),
    log:
        os.path.join(LOGS_DIR, "preprocessing", "barcode_diagnostic.log"),
    params:
        script=os.path.join(SCRIPTS_DIR, "barcode-diagnostic.py"),
        sample_size=barcode_cfg.get("sample_size", 10000),
        min_success=barcode_cfg.get("min_success", 1.0),
        # Pass --no_rc only when reverse-complement testing is disabled.
        rc_flag="--no_rc" if not barcode_cfg.get("test_reverse_complement", True) else "",
    conda:
        DIAGNOSTICS_ENV
    threads: 1
    resources:
        mem_mb=MEM_MB,
    shell:
        r"""
        python {params.script} \
            --index1 {input.index1} \
            --index2 {input.index2} \
            --mapping_file {input.mapping} \
            --sample_size {params.sample_size} \
            --min_success {params.min_success} \
            {params.rc_flag} \
            --output {output.report} \
            --verbose > {log} 2>&1
        """


rule analyze_unmapped_reads:
    """Screen undetermined (unmapped) read pairs for contaminant sequences."""
    input:
        forward_reads=UNMAPPED_R1,
        reverse_reads=UNMAPPED_R2,
        contaminants=CONTAMINANTS,
    output:
        report=os.path.join(UNMAPPED_DIR, "contaminant_report.txt"),
    log:
        os.path.join(LOGS_DIR, "preprocessing", "unmapped_reads.log"),
    params:
        script=os.path.join(SCRIPTS_DIR, "analyze_unmapped_reads.py"),
    conda:
        DIAGNOSTICS_ENV
    threads: 1
    resources:
        mem_mb=MEM_MB,
    shell:
        r"""
        python {params.script} \
            --forward {input.forward_reads} \
            --reverse {input.reverse_reads} \
            --contaminants {input.contaminants} > {output.report} 2> {log}
        """


rule validate_raw_fastq:
    """Count reads and validate FASTQ format for each raw R1 file."""
    input:
        fastq=sample_r1,
    output:
        report=os.path.join(FASTQ_CHECK_DIR, "{sample}_R1_validation.txt"),
    log:
        os.path.join(LOGS_DIR, "preprocessing", "fastq_check", "{sample}.log"),
    params:
        script=os.path.join(SCRIPTS_DIR, "diagnostic.py"),
    conda:
        DIAGNOSTICS_ENV
    threads: 1
    resources:
        mem_mb=MEM_MB,
    shell:
        r"""
        python {params.script} {input.fastq} > {output.report} 2> {log}
        """