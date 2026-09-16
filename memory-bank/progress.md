# Progress — 16S Amplicon Sequencing Pipeline

## Status Legend
- [x] Complete
- [~] In progress
- [ ] Not started

## What Works
- [x] Repository skeleton (`config/`, `data/`, `docs/`, `workflow/`)
- [x] `.gitignore` covering raw data, results, logs, resources, OS/IDE files
- [x] `README.md` with workflow overview, pre-QC stage, and repo layout
- [x] `.clinerules` with conventions and coding standards
- [x] Memory Bank with foundational context files
- [x] `workflow/Snakefile` entry point (includes `common.smk` + `preprocessing.smk`)
- [x] `config/config.yaml` (parameters incl. `diagnostics:` block)
- [x] `config/samples.tsv` samplesheet template
- [x] `workflow/rules/common.smk` (shared helpers, wildcards, constants)
- [x] **Pre-QC diagnostics stage** (`workflow/rules/preprocessing.smk`):
  - [x] `check_index_fastq` -> `index_fastq_check.py`
  - [x] `diagnose_barcode_formats` -> `barcode-diagnostic.py`
  - [x] `analyze_unmapped_reads` -> `analyze_unmapped_reads.py`
  - [x] `validate_raw_fastq` -> `diagnostic.py`
- [x] `workflow/envs/diagnostics.yaml` (pinned Python/Biopython/pandas)
- [x] **Quality-control stage** (`workflow/rules/qc.smk`):
  - [x] `run_fastqc` -> FastQC per-file reports (R1/R2)
  - [x] `run_fastp` -> adapter/quality trimming + HTML/JSON reports
  - [x] `run_multiqc` -> aggregated MultiQC report
- [x] `workflow/envs/qc.yaml` (pinned fastqc/fastp/multiqc)
- [x] Bundled EMP-style test dataset wired in (`data/raw/` + `config/`)
- [x] `common.smk` honours explicit `r1`/`r2` samplesheet columns
- [x] Full default target resolves with `snakemake -n` (7 jobs)

## What's Left (Roadmap)

### Phase 1 - Foundations
- [x] `workflow/Snakefile` entry point
- [x] `config/config.yaml` (parameters, primers, reference DB)
- [x] `config/samples.tsv` template
- [x] `workflow/rules/common.smk` (shared helpers, wildcards)
- [x] Python env pinned in `workflow/envs/` for base tooling

### Phase 2 - Preprocessing Rules
- [x] `workflow/rules/preprocessing.smk` + `envs/diagnostics.yaml`
      (pre-QC input/index diagnostics)
- [x] `workflow/rules/qc.smk` + `envs/qc.yaml` (FastQC/MultiQC/fastp)
- [ ] `workflow/rules/demultiplex.smk` + env (if needed)
- [ ] `workflow/rules/primer_removal.smk` + `envs/cutadapt.yaml`
- [ ] `workflow/rules/filtering.smk` + env

### Phase 3 - Denoising & ASVs
- [ ] `workflow/rules/denoising.smk` + `envs/dada2.yaml`
- [ ] `workflow/scripts/denoise.R`
- [ ] Chimera removal step
- [ ] ASV table + representative sequences outputs

### Phase 4 - Taxonomy
- [ ] `workflow/rules/taxonomy.smk` + env
- [ ] Reference DB download/record step (version + checksum)
- [ ] `workflow/scripts/assign_taxonomy.R`
- [ ] Taxonomy table output

### Phase 5 - Depth & Downstream
- [ ] `workflow/rules/depth.smk` (rarefaction curves, depth assessment)
- [ ] `workflow/rules/diversity.smk` (alpha/beta diversity)
- [ ] `workflow/notebooks/` R Markdown analyses
- [ ] Differential abundance analysis
- [ ] Figures/tables for publication

### Phase 6 - Documentation & CI
- [ ] Expand `README.md` (install, usage, outputs)
- [ ] `docs/methods.md` (methods write-up)
- [ ] `docs/` parameter reference
- [ ] Test dataset under `data/test/`
- [ ] Smoke-test / CI target (e.g., `snakemake -n` dry-run)

## Known Issues
- **Contaminant reference missing** - `resources/contaminants.fasta` is not
  present, so the `analyze_unmapped_reads` rule is skipped automatically. Add
  the reference to enable the contaminant screen.
- **No documentation yet** - `docs/` is empty.
- **Open technical decisions** - denoiser, reference DB, region/primers,
  read layout, and target execution environment are all undecided (see
  `activeContext.md`).
- **QC not yet executed** - the DAG resolves, but a full `--use-conda` run
  (which builds the FastQC/fastp/MultiQC envs) has not been performed.

## Next Milestone
**Milestone 1: Runnable skeleton.** A `Snakefile` + `config.yaml` + QC rules
that run end-to-end on the bundled `data/raw/` test dataset with
`snakemake --use-conda`. The pre-QC diagnostics and QC stages are both wired in
and the DAG resolves; the remaining step is an actual `--use-conda` execution.

## Notes
- Update this file at the end of every work session.
- Move completed items up to "What Works" and remove them from the roadmap.