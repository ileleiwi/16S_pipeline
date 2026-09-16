# Project Brief — 16S Amplicon Sequencing Pipeline

## Purpose
Build a **reproducible, documented bioinformatics workflow** for processing 16S
rRNA gene amplicon sequencing data — from raw FASTQ reads through taxonomic
assignment and downstream ecological analysis.

## Core Goals
1. **Reproducibility** — every step deterministic, versioned, and re-runnable
   on any machine with the same inputs.
2. **Modularity** — each processing step is an independent, testable unit.
3. **Transparency** — the full path from raw reads to biological conclusions is
   documented and auditable.
4. **Portability** — environments are self-contained and pinned; no reliance on
   the host system's software.
5. **Scalability** — capable of running on a laptop for test data and on an HPC
   cluster for full datasets.

## Scope

### In Scope
- Raw FASTQ quality control and demultiplexing
- Primer removal (e.g., 16S V3–V4, V4 primers)
- Quality filtering and denoising (ASV generation)
- Taxonomic assignment against a reference database
- Sequencing-depth assessment / rarefaction
- Downstream ecological/diversity analysis

### Out of Scope (for now)
- Metagenomic (shotgun) data
- Long-read (PacBio/Nanopore) full-length 16S
- Clinical diagnostic reporting
- Real-time/live sequencing analysis

## Workflow Stages
1. Quality control
2. Demultiplexing
3. Primer removal
4. Quality filtering
5. Denoising
6. ASV generation
7. Taxonomic assignment
8. Sequencing-depth assessment
9. Downstream analysis

## Success Criteria
- A single command runs the entire pipeline end-to-end.
- All tool and reference database versions are recorded.
- Intermediate and final outputs are reproducible bit-for-bit (or close enough
  for the denoising step to be deterministic given fixed seeds).
- Documentation allows a new lab member to run the pipeline unaided.

## Key Constraints
- Do not commit raw sequencing data, large outputs, or reference databases.
- Pin all tool versions in `workflow/envs/*.yaml`.
- Follow the documented workflow order.