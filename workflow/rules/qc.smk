# =============================================================================
# qc.smk - quality control stage
# =============================================================================
# Runs after the pre-QC diagnostics stage (preprocessing.smk). Produces:
#   - FastQC per-file raw-read quality reports (R1 and R2)
#   - fastp adapter/quality-trimmed reads plus per-sample HTML/JSON reports
#   - MultiQC aggregated report over all FastQC and fastp outputs
#
# All tunable parameters come from the `qc:` block in config/config.yaml. Tool
# versions are pinned in workflow/envs/qc.yaml.
# =============================================================================

import os

# -----------------------------------------------------------------------------
# Paths and configuration
# -----------------------------------------------------------------------------
QC_ENV = os.path.join(str(WORKFLOW_DIR), "envs", "qc.yaml")

qc_cfg = config.get("qc", {})

QC_DIR = os.path.join(RESULTS_DIR, "qc")
FASTQC_DIR = os.path.join(QC_DIR, "fastqc")
FASTP_DIR = os.path.join(QC_DIR, "fastp")
MULTIQC_DIR = os.path.join(QC_DIR, "multiqc")
STAGING_DIR = os.path.join(QC_DIR, "staging")

THREADS = config["resources"]["threads_default"]
MEM_MB = config["resources"]["mem_mb_default"]

# fastp parameters. A value of 0 (or an empty string) means "no limit" and the
# corresponding flag is omitted from the command line.
MIN_QUALITY = qc_cfg.get("min_quality", 20)
MIN_LENGTH = qc_cfg.get("min_length", 50)
MAX_LENGTH = qc_cfg.get("max_length", 0)
MAX_N = qc_cfg.get("max_n", 0)
TRIM_FRONT = qc_cfg.get("trim_front", 0)
TRIM_TAIL = qc_cfg.get("trim_tail", 0)
ADAPTER_DETECTION = qc_cfg.get("adapter_detection", True)

# Pre-compute optional fastp flags so the shell command stays readable.
MAX_LENGTH_FLAG = f"--length_limit {MAX_LENGTH}" if MAX_LENGTH else ""
MAX_N_FLAG = f"--n_base_limit {MAX_N}" if MAX_N else ""
TRIM_FRONT_FLAG = (
    f"--trim_front1 {TRIM_FRONT} --trim_front2 {TRIM_FRONT}" if TRIM_FRONT else ""
)
TRIM_TAIL_FLAG = (
    f"--trim_tail1 {TRIM_TAIL} --trim_tail2 {TRIM_TAIL}" if TRIM_TAIL else ""
)
ADAPTER_FLAG = "--detect_adapter_for_pe" if ADAPTER_DETECTION else ""


def qc_targets():
    """Return the outputs produced by the QC stage.

    Used by the main Snakefile to build the default target. Returns an empty
    list when there are no samples so the workflow still parses.
    """
    targets = [os.path.join(MULTIQC_DIR, "multiqc_report.html")]
    for sample in SAMPLES:
        targets.append(os.path.join(FASTQC_DIR, f"{sample}_R1_fastqc.html"))
        targets.append(os.path.join(FASTQC_DIR, f"{sample}_R2_fastqc.html"))
        targets.append(os.path.join(FASTP_DIR, f"{sample}_R1.fastq.gz"))
        targets.append(os.path.join(FASTP_DIR, f"{sample}_R2.fastq.gz"))
    return targets


# -----------------------------------------------------------------------------
# Rules
# -----------------------------------------------------------------------------
rule run_fastqc:
    """Run FastQC on raw reads for each sample (R1 and R2)."""
    input:
        r1=sample_r1,
        r2=sample_r2,
    output:
        html1=os.path.join(FASTQC_DIR, "{sample}_R1_fastqc.html"),
        html2=os.path.join(FASTQC_DIR, "{sample}_R2_fastqc.html"),
        zip1=os.path.join(FASTQC_DIR, "{sample}_R1_fastqc.zip"),
        zip2=os.path.join(FASTQC_DIR, "{sample}_R2_fastqc.zip"),
    log:
        os.path.join(LOGS_DIR, "qc", "fastqc", "{sample}.log"),
    params:
        outdir=FASTQC_DIR,
        staging=STAGING_DIR,
    conda:
        QC_ENV
    threads: THREADS
    resources:
        mem_mb=MEM_MB,
    shell:
        r"""
        mkdir -p {params.outdir} {params.staging}
        # Stage inputs under sample-based names so FastQC output is predictable
        # regardless of the original raw filenames.
        ln -sf {input.r1} {params.staging}/{wildcards.sample}_R1.fastq.gz
        ln -sf {input.r2} {params.staging}/{wildcards.sample}_R2.fastq.gz
        fastqc --threads {threads} --outdir {params.outdir} \
            {params.staging}/{wildcards.sample}_R1.fastq.gz \
            {params.staging}/{wildcards.sample}_R2.fastq.gz > {log} 2>&1
        rm -f {params.staging}/{wildcards.sample}_R1.fastq.gz \
              {params.staging}/{wildcards.sample}_R2.fastq.gz
        """


rule run_fastp:
    """Trim adapters and low-quality bases with fastp."""
    input:
        r1=sample_r1,
        r2=sample_r2,
    output:
        r1=os.path.join(FASTP_DIR, "{sample}_R1.fastq.gz"),
        r2=os.path.join(FASTP_DIR, "{sample}_R2.fastq.gz"),
        html=os.path.join(FASTP_DIR, "{sample}_fastp.html"),
        json=os.path.join(FASTP_DIR, "{sample}_fastp.json"),
    log:
        os.path.join(LOGS_DIR, "qc", "fastp", "{sample}.log"),
    params:
        min_quality=MIN_QUALITY,
        min_length=MIN_LENGTH,
        max_length_flag=MAX_LENGTH_FLAG,
        max_n_flag=MAX_N_FLAG,
        trim_front_flag=TRIM_FRONT_FLAG,
        trim_tail_flag=TRIM_TAIL_FLAG,
        adapter_flag=ADAPTER_FLAG,
    conda:
        QC_ENV
    threads: THREADS
    resources:
        mem_mb=MEM_MB,
    shell:
        r"""
        fastp \
            --in1 {input.r1} --in2 {input.r2} \
            --out1 {output.r1} --out2 {output.r2} \
            --qualified_quality_phred {params.min_quality} \
            --length_required {params.min_length} \
            {params.max_length_flag} \
            {params.max_n_flag} \
            {params.trim_front_flag} \
            {params.trim_tail_flag} \
            {params.adapter_flag} \
            --thread {threads} \
            --html {output.html} --json {output.json} > {log} 2>&1
        """


rule run_multiqc:
    """Aggregate FastQC and fastp reports into a single MultiQC report."""
    input:
        fastqc=expand(
            os.path.join(FASTQC_DIR, "{sample}_R1_fastqc.zip"), sample=SAMPLES
        )
        + expand(
            os.path.join(FASTQC_DIR, "{sample}_R2_fastqc.zip"), sample=SAMPLES
        ),
        fastp=expand(
            os.path.join(FASTP_DIR, "{sample}_fastp.json"), sample=SAMPLES
        ),
    output:
        html=os.path.join(MULTIQC_DIR, "multiqc_report.html"),
        data=os.path.join(MULTIQC_DIR, "multiqc_data", "multiqc_data.json"),
    log:
        os.path.join(LOGS_DIR, "qc", "multiqc.log"),
    params:
        outdir=MULTIQC_DIR,
        search_dir=QC_DIR,
    conda:
        QC_ENV
    threads: 1
    resources:
        mem_mb=MEM_MB,
    shell:
        r"""
        multiqc --force --outdir {params.outdir} {params.search_dir} > {log} 2>&1
        """