# System Patterns — 16S Amplicon Sequencing Pipeline

## High-Level Architecture
The pipeline is a **directed acyclic graph (DAG)** of steps managed by Snakemake.
Data flows linearly through preprocessing stages, then branches into taxonomy and
diversity analyses. Each node in the DAG is a Snakemake rule that declares its
inputs, outputs, logs, resources, and conda environment.

```
Raw FASTQ
   |
   v
Pre-QC diagnostics [preprocessing.smk]
   |
   v
QC (fastp/FastQC) --> Demultiplexing --> Primer removal (cutadapt)
   |
   v
Quality filtering --> Denoising (DADA2/Deblur) --> ASV table + seqs
   |
   v
Taxonomic assignment (SILVA/GTDB) --> Sequencing-depth assessment
   |
   v
Downstream analysis (diversity, ordination, differential abundance)
```

## Directory Layout & Responsibilities
| Path | Responsibility |
|---|---|
| `config/` | Samplesheets, parameter YAML, reference DB definitions |
| `data/` | Input data; `data/raw/` is gitignored |
| `docs/` | Methods, workflow notes, analysis write-ups |
| `workflow/Snakefile` | Entry point; includes all rule files |
| `workflow/rules/` | Modular Snakemake rules, one concern per file |
| `workflow/envs/` | Pinned conda env YAMLs, one per rule group |
| `workflow/scripts/` | R/Python scripts invoked by rules |
| `workflow/notebooks/` | Rmd/Jupyter for interactive downstream analysis |
| `results/` | Generated outputs (gitignored) |
| `logs/` | Per-rule logs (gitignored) |
| `resources/` | Reference databases (gitignored) |

## Key Technical Decisions
1. **Snakemake as the orchestrator** — provides DAG resolution, partial
   re-runs, cluster/HPC support, and per-rule conda envs.
2. **One concern per rule file** — e.g., `preprocessing.smk`, `qc.smk`,
   `primer_removal.smk`, `denoising.smk`, `taxonomy.smk`, `diversity.smk`.
3. **One conda env per tool/rule group** — isolates incompatible tool
   dependencies and pins versions.
4. **Config-driven design** — primers, parameters, paths, and DB versions live
   in `config/`, never hardcoded in rules or scripts.
5. **Per-rule `log:`** — every rule writes a log; errors are diagnosable.
6. **Per-rule `threads:`/`resources:`** — explicit resource requests for HPC.
7. **Scripts are thin wrappers** — R/Python scripts do the science; Snakemake
   handles wiring, environment, and provenance.
8. **Relative paths everywhere** — snapshots and containers work unchanged.

## Design Patterns
### Rule Pattern
```python
rule denoise:
    input:
        r1="results/filtered/{sample}_R1.fastq.gz",
        r2="results/filtered/{sample}_R2.fastq.gz"
    output:
        table="results/denoising/asv_table.tsv",
        seqs="results/denoising/asv_seqs.fasta"
    log:
        "logs/denoising/{sample}.log"
    threads: 8
    conda:
        "../envs/dada2.yaml"
    script:
        "../scripts/denoise.R"
```

### Config Pattern
```yaml
# config/config.yaml
samples: config/samples.tsv
primers:
  forward: GTGCCAGCMGCCGCGGTAA
  reverse: GGACTACHVGGGTWTCTAAT
denoising:
  trunc_len_f: 240
  trunc_len_r: 160
reference:
  database: silva
  version: "138.1"
  path: resources/silva_138.1.fasta
```

### Naming Conventions
- Rules: `snake_case`, verb-first (`run_qc`, `remove_primers`, `assign_taxonomy`).
- Outputs: `results/<stage>/<sample>_<suffix>.<ext>`.
- Wildcards: `{sample}`, `{region}`, `{db}`.
- Env files: `<toolgroup>.yaml` (e.g., `dada2.yaml`, `qiime2.yaml`).

## Data Flow Contracts
- FASTQ files: `*_R1.fastq.gz` / `*_R2.fastq.gz` naming.
- ASV table: TSV with ASVs as rows, samples as columns (or vice versa — document).
- Representative sequences: FASTA keyed by ASV ID.
- Taxonomy: TSV with `ASV`, `Kingdom`...`Species`, `Confidence`.
- Metadata: samplesheet with `sample`, `group`, `batch`, and covariates.

## Reproducibility Patterns
- Record tool versions (`conda list`, `--version`) into `results/versions/`.
- Record reference DB name + version + checksum in config and logs.
- Set and record random seeds for any stochastic step.
- Use Snakemake `--use-conda` and, where possible, containerize.