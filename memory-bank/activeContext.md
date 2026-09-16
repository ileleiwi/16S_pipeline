# Active Context — 16S Amplicon Sequencing Pipeline

## Current Work Focus
**Quality-control stage + test-data wiring.** The pipeline skeleton (Snakefile,
config, samplesheet, `common.smk`) is in place, the pre-QC diagnostics stage is
implemented, and the QC stage (`workflow/rules/qc.smk`) now exists. A bundled
EMP-style test dataset under `data/raw/` is wired into the config and
samplesheet so the workflow can be exercised end-to-end.

## Recent Changes
- Created `.clinerules` with project conventions, coding standards, and
  Memory Bank instructions.
- Created the Memory Bank (`memory-bank/`) with foundational context files.
- Added `workflow/rules/preprocessing.smk` wrapping four diagnostic scripts
  ahead of QC:
  - `index_fastq_check.py` -> `check_index_fastq` (index congruency + barcodes)
  - `barcode-diagnostic.py` -> `diagnose_barcode_formats` (barcode format)
  - `analyze_unmapped_reads.py` -> `analyze_unmapped_reads` (contaminant screen)
  - `diagnostic.py` -> `validate_raw_fastq` (per-file FASTQ validation)
- Added `workflow/envs/diagnostics.yaml` (pinned Python/Biopython/pandas).
- Added a `diagnostics:` block to `config/config.yaml` (inputs + barcode
  parameters, with an `enabled` switch).
- Updated `workflow/Snakefile` to include `preprocessing.smk` ahead of QC and to
  add the diagnostic outputs to the default target. The QC include is now
  conditional so the workflow parses while `rules/qc.smk` is still missing.
- Expanded `README.md` with the pre-QC stage and repository layout.
- Added `workflow/rules/qc.smk` (FastQC + fastp + MultiQC) and
  `workflow/envs/qc.yaml` (pinned fastqc=0.12.1, fastp=0.23.4, multiqc=1.21).
- Wired the QC stage into `workflow/Snakefile` via a `qc_targets()` helper and
  removed the conditional/include guard now that the rule file exists.
- Wired the bundled test dataset (`data/raw/read1.fastq.gz`, `read2.fastq.gz`,
  `index1.fastq.gz`, `index2.fastq.gz`, `mapping_file.txt`) into
  `config/config.yaml` (diagnostics inputs) and `config/samples.tsv`
  (`test_sample` row; `tissue` column renamed to `sample_type`).
- Updated `common.smk` so `sample_r1`/`sample_r2` honour the explicit `r1`/`r2`
  samplesheet columns (falling back to `raw_dir/<sample>_R<1|2>.fastq.gz`).
- Made the contaminant screen (`analyze_unmapped_reads`) optional: its output is
  only requested when the contaminant reference and undetermined reads exist.
- Verified the full default target resolves with `snakemake -n` in the
  `16s_pipeline` conda env (7 jobs).

## Current State of the Repository
```
16S_pipeline/
|-- .clinerules              # Cline rules
|-- .gitignore               # Ignores raw data, results, logs, resources
|-- README.md                # Workflow overview + pre-QC stage
|-- config/
|   |-- config.yaml          # Parameters incl. `diagnostics:` block
|   `-- samples.tsv          # Samplesheet template
|-- data/                    # raw/ is gitignored
|-- docs/                    # EMPTY - needs methods docs
|-- memory-bank/             # Context files
`-- workflow/
    |-- Snakefile            # Entry point; includes preprocessing + QC
    |-- envs/
    |   |-- diagnostics.yaml # Pinned env for diagnostic scripts
    |   `-- qc.yaml          # Pinned env: FastQC/fastp/MultiQC (NEW)
    |-- notebooks/           # EMPTY - downstream analysis
    |-- rules/
    |   |-- common.smk       # Shared helpers/wildcards
    |   |-- preprocessing.smk# Pre-QC diagnostics
    |   `-- qc.smk           # Quality control (NEW)
    `-- scripts/
        |-- index_fastq_check.py
        |-- barcode-diagnostic.py
        |-- analyze_unmapped_reads.py
        `-- diagnostic.py
```

## Next Steps (Suggested)
1. **Provide `resources/contaminants.fasta`** - the contaminant screen is
   skipped until this reference exists. Alternatively leave it out and rely on
   the automatic skip.
2. **Add `docs/methods.md`** - methods write-up for publication.
3. **Continue the pipeline** - demultiplexing, primer removal, filtering,
   denoising, taxonomy, depth, diversity.
4. **Smoke-test / CI target** - run the QC stage on the bundled `data/raw/`
   test data (a full `--use-conda` run builds the FastQC/fastp/MultiQC envs).

## Open Questions / Decisions Pending
- Which denoiser? (DADA2 vs Deblur vs UNOISE3)
- Which reference database and version? (SILVA 138.1, Greengenes2, GTDB)
- Which 16S region / primer set? (V3-V4, V4, etc.)
- Single-end or paired-end reads?
- Target execution environment: local only, or HPC (Slurm)?
- Is QIIME 2 part of the stack, or a pure Snakemake + R/Python approach?

## Important Notes for Future Sessions
- Read ALL Memory Bank files at the start of every task.
- The pipeline order is fixed and documented - do not reorder without updating
  `docs/` and `README.md`.
- Keep raw data, results, logs, and reference DBs out of git.