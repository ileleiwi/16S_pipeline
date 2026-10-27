# =============================================================================
# demultiplex.smk - split a multiplexed / Undetermined run into per-sample FASTQs
# =============================================================================
# Runs BEFORE quality control whenever `demultiplex.enabled: true`.
#
# Reads of a run submitted without per-sample indexes (e.g. with a dummy index)
# all land in the Undetermined_* pool; the real barcodes then live either in a
# separate index read (I1 - the `emp-paired` layout used by
# data/raw/20240209_16S) or inside the reads. The QIIME 2 EMP action consumes the
# barcode read plus the forward/reverse reads, assigns each pair to its sample
# using the barcode column of a QIIME 2 metadata file, and writes demultiplexed
# FASTQs.
#
# Why there is no `conda:` directive here: QIIME 2 publishes no osx-arm64 builds,
# so Snakemake cannot create that environment on Apple Silicon (it would target
# osx-arm64). The environment is pinned instead in
#   workflow/envs/qiime2.yaml       (verbatim QIIME 2 2026.7 distro file)
#   workflow/envs/16s_pipeline.yaml (QIIME 2 + pipeline tools, when compatible)
# and reached through `demultiplex.qiime_cmd`:
#   "qiime"                                  the driver env holds QIIME 2
#   "conda run -n 16s_pipeline-qiime2 qiime" QIIME 2 in its own environment
# On Linux/HPC use `demultiplex.qiime_container` with `--sdm container` instead.
# =============================================================================

import os

dmx_cfg = config.get("demultiplex", {})
DMX_ENABLED = bool(dmx_cfg.get("enabled", False))
DMX_BARCODE_COLUMN = dmx_cfg.get("barcode_column", "BarcodeSequence")
DMX_SAMPLE_SIZE = int(dmx_cfg.get("sample_size", 0) or 0)
DMX_IGNORE_DESC = bool(dmx_cfg.get("ignore_description_mismatch", False))

# -----------------------------------------------------------------------------
# How QIIME 2 is reached
# -----------------------------------------------------------------------------
# `demultiplex.qiime_cmd` is an explicit command prefix and wins when set.
# Otherwise, when `demultiplex.qiime_container` is configured, QIIME 2 is run
# inside that image. The repository is mounted at its own absolute path so every
# path Snakemake passes (results/, logs/, config/) is valid inside the container.
QIIME_IMAGE = (dmx_cfg.get("qiime_container") or "").replace("docker://", "")
QIIME_ENGINE = dmx_cfg.get("qiime_container_engine", "docker")
QIIME_CMD = dmx_cfg.get("qiime_cmd") or ""
if not QIIME_CMD and QIIME_IMAGE:
    if QIIME_ENGINE in ("apptainer", "singularity"):
        # Linux/HPC: `snakemake --sdm container` can also be used instead.
        QIIME_CMD = f"{QIIME_ENGINE} exec --bind {ROOT_DIR}:{ROOT_DIR} {QIIME_IMAGE} qiime"
    else:
        QIIME_CMD = (
            f"{QIIME_ENGINE} run --rm "
            f"{dmx_cfg.get('qiime_container_user') or ''} "
            f"-v {ROOT_DIR}:{ROOT_DIR} -w {ROOT_DIR} {QIIME_IMAGE} qiime"
        ).replace("  ", " ")

DMX_DIR = os.path.join(RESULTS_DIR, "demultiplex")
DMX_STAGING = os.path.join(DMX_DIR, "staging")
DMX_QZA = os.path.join(DMX_DIR, "emp-sequences.qza")
DMX_RESULT = os.path.join(DMX_DIR, "demux.qza")
DMX_EC = os.path.join(DMX_DIR, "demux-error-correction.qza")
DMX_QZV = os.path.join(DMX_DIR, "demux-summary.qzv")
DMX_EXPORT = os.path.join(DMX_DIR, "exported")
DMX_STATS = os.path.join(DMX_DIR, "demux_stats.tsv")

DMX_FORWARD = os.path.join(ROOT_DIR, dmx_cfg.get("forward") or "")
DMX_REVERSE = os.path.join(ROOT_DIR, dmx_cfg.get("reverse") or "")
DMX_BARCODES = os.path.join(ROOT_DIR, dmx_cfg.get("barcodes") or "")
DMX_METADATA = os.path.join(ROOT_DIR, dmx_cfg.get("barcodes_metadata") or "")

DMX_MEM_MB = config["resources"]["mem_mb_default"]


def demux_export_outputs():
    """Per-sample FASTQs produced by the export rule."""
    return [demux_sample_fastq(sample, read) for sample in SAMPLES for read in ("r1", "r2")]


def demultiplex_targets():
    """Return the outputs of the demultiplex stage (empty list when disabled)."""
    if not DMX_ENABLED:
        return []
    return [DMX_QZV, DMX_STATS] + demux_export_outputs()


# -----------------------------------------------------------------------------
# Rules
# -----------------------------------------------------------------------------
rule import_emp_sequences:
    """Stage the multiplexed reads and import them as EMPPairedEndSequences.

    The staging directory is created inside this job (`params.staging`) instead of
    being a declared `directory()` output: Snakemake writes a
    `.snakemake_timestamp` marker into declared directory outputs, and QIIME 2
    imports a directory whose only contents may be `forward.fastq.gz`,
    `reverse.fastq.gz` and `barcodes.fastq.gz`.
    """
    input:
        # NOTE: "forward"/"reverse" are reserved names in Snakemake, hence r1/r2.
        r1=DMX_FORWARD,
        r2=DMX_REVERSE,
        barcodes=DMX_BARCODES,
    output:
        qza=DMX_QZA,
    log:
        os.path.join(LOGS_DIR, "demultiplex", "import.log"),
    params:
        staging=DMX_STAGING,
        sample_size=DMX_SAMPLE_SIZE,
        qiime=QIIME_CMD,
    threads: 1
    resources:
        mem_mb=DMX_MEM_MB,
    shell:
        r"""
        rm -rf {params.staging}
        mkdir -p {params.staging}
        if [ {params.sample_size} -gt 0 ]; then
            # Smoke-test mode: keep only the first N read pairs of each file.
            # `set +o pipefail` is required because Snakemake runs job shells with
            # `set -euo pipefail` and gzip exits 141 (SIGPIPE) once head has
            # collected enough lines.
            for pair in "forward {input.r1}" "reverse {input.r2}" "barcodes {input.barcodes}"; do
                set -- $pair;
                ( set +o pipefail; gzip -dc "$2" | head -n $(( {params.sample_size} * 4 )) | gzip > "{params.staging}/$1".fastq.gz );
            done
        else
            ln -sf {input.r1} {params.staging}/forward.fastq.gz
            ln -sf {input.r2} {params.staging}/reverse.fastq.gz
            ln -sf {input.barcodes} {params.staging}/barcodes.fastq.gz
        fi
        {params.qiime} tools import \
            --type EMPPairedEndSequences \
            --input-path {params.staging} \
            --output-path {output.qza} > {log} 2>&1
        """


rule demultiplex_stage:
    """Convenience target: run only the demultiplex stage.

    Use `snakemake demultiplex_stage` to demultiplex without starting QC, e.g. for
    a smoke test on a slice of the run:
      snakemake --use-conda --cores 4 --config dataset=20240209_16S \
          demux_sample_size=100000 demultiplex_stage
    """
    input:
        demultiplex_targets()


rule demultiplex_emp:
    """Demultiplex the barcode read against the metadata barcode column."""
    input:
        qza=DMX_QZA,
        metadata=DMX_METADATA,
    output:
        demux=DMX_RESULT,
        error_correction=DMX_EC,
    log:
        os.path.join(LOGS_DIR, "demultiplex", "demux.log"),
    params:
        qiime=QIIME_CMD,
        barcode_column=DMX_BARCODE_COLUMN,
        golay_ec=("--p-golay-error-correction" if dmx_cfg.get("golay_error_correction", True)
                  else "--p-no-golay-error-correction"),
        rev_comp_barcodes=("--p-rev-comp-barcodes" if dmx_cfg.get("rev_comp_barcodes", False)
                           else "--p-no-rev-comp-barcodes"),
        rev_comp_mapping=("--p-rev-comp-mapping-barcodes"
                          if dmx_cfg.get("rev_comp_mapping_barcodes", False)
                          else "--p-no-rev-comp-mapping-barcodes"),
        ignore_desc=("--p-ignore-description-mismatch" if DMX_IGNORE_DESC
                     else "--p-no-ignore-description-mismatch"),
    threads: 1
    resources:
        mem_mb=DMX_MEM_MB,
    shell:
        r"""
        {params.qiime} demux emp-paired \
            --i-seqs {input.qza} \
            --m-barcodes-file {input.metadata} \
            --m-barcodes-column {params.barcode_column} \
            {params.golay_ec} \
            {params.rev_comp_barcodes} \
            {params.rev_comp_mapping} \
            {params.ignore_desc} \
            --o-per-sample-sequences {output.demux} \
            --o-error-correction-details {output.error_correction} > {log} 2>&1
        """


rule summarize_demultiplex:
    """Per-sample read counts and quality summary for the demux artifact."""
    input:
        demux=DMX_RESULT,
    output:
        qzv=DMX_QZV,
    log:
        os.path.join(LOGS_DIR, "demultiplex", "summarize.log"),
    params:
        qiime=QIIME_CMD,
    threads: 1
    resources:
        mem_mb=DMX_MEM_MB,
    shell:
        r"""
        {params.qiime} demux summarize \
            --i-data {input.demux} \
            --o-visualization {output.qzv} > {log} 2>&1
        """


rule export_demultiplexed_reads:
    """Export per-sample FASTQs and normalise them to <sample>_R{1,2}.fastq.gz."""
    input:
        demux=DMX_RESULT,
        # Declared so that edits to the normalisation logic re-trigger this rule.
        script=os.path.join(SCRIPTS_DIR, "normalize_demux_export.py"),
    output:
        stats=DMX_STATS,
        reads=demux_export_outputs(),
    log:
        os.path.join(LOGS_DIR, "demultiplex", "export.log"),
    params:
        script=os.path.join(SCRIPTS_DIR, "normalize_demux_export.py"),
        export_dir=DMX_EXPORT,
        outdir=DMX_DIR,
        qiime=QIIME_CMD,
        samples=" ".join(SAMPLES),
    threads: 1
    resources:
        mem_mb=DMX_MEM_MB,
    shell:
        r"""
        rm -rf {params.export_dir}
        {params.qiime} tools export --input-path {input.demux} \
            --output-path {params.export_dir} >> {log} 2>&1
        python {params.script} \
            --export_dir {params.export_dir} \
            --outdir {params.outdir} \
            --stats {output.stats} \
            --samples {params.samples} >> {log} 2>&1
        """
