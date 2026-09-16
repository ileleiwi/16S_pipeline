# 16S Amplicon Sequencing Pipeline

A reproducible workflow for processing 16S rRNA gene amplicon sequencing data.

## Workflow

```
Raw FASTQ files
→ Pre-QC input & index diagnostics
→ Quality control
→ Demultiplexing
→ Primer removal
→ Quality filtering
→ Denoising
→ ASV generation
→ Taxonomic assignment
→ Sequencing-depth assessment
→ Downstream analysis
```

The pipeline is orchestrated with Snakemake. Parameters live in
`config/config.yaml`, samples in `config/samples.tsv`, rules in
`workflow/rules/`, tool scripts in `workflow/scripts/`, and pinned conda
environments in `workflow/envs/`.

## Pre-QC input & index diagnostics

This stage runs **before quality control** and validates the raw inputs and the
library/index configuration so that problems are caught early. It is purely
diagnostic: it reports issues rather than modifying data.

| Rule | Script | Purpose | Output |
|---|---|---|---|
| `check_index_fastq` | `index_fastq_check.py` | Index (I1/I2) congruency check and barcode extraction | `results/preprocessing/index_check/data/` |
| `diagnose_barcode_formats` | `barcode-diagnostic.py` | Determine the optimal barcode format for demultiplexing | `results/preprocessing/barcode_diagnostic/barcode_format_report.txt` |
| `analyze_unmapped_reads` | `analyze_unmapped_reads.py` | Screen undetermined reads for contaminants | `results/preprocessing/unmapped_reads/contaminant_report.txt` |
| `validate_raw_fastq` | `diagnostic.py` | Read count / FASTQ format validation per raw R1 file | `results/preprocessing/fastq_check/` |

The stage is configured under the `diagnostics:` key in `config/config.yaml` and
can be switched off with `diagnostics.enabled: false`.

## Quality control

The QC stage runs after the pre-QC diagnostics and produces per-sample quality
reports plus adapter/quality-trimmed reads.

| Rule | Tool | Purpose | Output |
|---|---|---|---|
| `run_fastqc` | FastQC | Per-file raw-read quality reports (R1 and R2) | `results/qc/fastqc/` |
| `run_fastp` | fastp | Adapter/quality trimming + per-sample HTML/JSON reports | `results/qc/fastp/` |
| `run_multiqc` | MultiQC | Aggregate all FastQC and fastp reports | `results/qc/multiqc/multiqc_report.html` |

Parameters live under the `qc:` key in `config/config.yaml`; tool versions are
pinned in `workflow/envs/qc.yaml`.

## Test dataset

A small EMP-style dual-indexed test run is bundled under `data/raw/`:

| File | Description |
|---|---|
| `read1.fastq.gz` / `read2.fastq.gz` | Paired-end reads |
| `index1.fastq.gz` / `index2.fastq.gz` | Dual index (I1/I2) reads |
| `mapping_file.txt` | EMP barcode-to-sample mapping |

`config/samples.tsv` points the `test_sample` row at these reads, and the
`diagnostics:` block in `config/config.yaml` references the index and mapping
files. The contaminant screen is skipped automatically when
`resources/contaminants.fasta` is absent.

## Running

From the `workflow/` directory:

```bash
snakemake --use-conda --cores <N>
```

Dry run (plan the DAG without executing):

```bash
snakemake -n --use-conda
```

## Repository layout

| Path | Responsibility |
|---|---|
| `config/` | Samplesheet, parameters, reference DB definitions |
| `data/` | Input data (raw FASTQ is gitignored) |
| `docs/` | Methods, workflow notes, analysis write-ups |
| `workflow/Snakefile` | Entry point; includes rule files |
| `workflow/rules/` | Modular Snakemake rules (one concern per file) |
| `workflow/scripts/` | R/Python scripts invoked by rules |
| `workflow/envs/` | Pinned conda environments |
| `workflow/notebooks/` | Interactive downstream analysis |
| `results/`, `logs/`, `resources/` | Generated outputs (gitignored) |