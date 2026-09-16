# Product Context — 16S Amplicon Sequencing Pipeline

## Why This Project Exists
16S rRNA gene amplicon sequencing is a foundational method in microbial ecology,
microbiome research, and environmental microbiology. However, analysis pipelines
are notoriously difficult to reproduce: tool versions drift, parameters are
undocumented, and reference databases are unversioned. This project exists to
provide a **single, self-contained, reproducible workflow** that removes those
pain points.

## Problems It Solves
| Problem | How this pipeline addresses it |
|---|---|
| Irreproducible analyses | Pinned conda envs + Snakemake DAG + recorded versions |
| Ad-hoc, undocumented steps | One rule per workflow stage, documented in `docs/` |
| Environment drift | `workflow/envs/*.yaml` with pinned tool versions |
| "Works on my machine" | Relative paths, container-ready, no `~`/`$HOME` |
| Silent data loss | Explicit QC/denoising reports at each stage |
| Inconsistent taxonomy | Versioned reference DB recorded in config |

## Target Users
- **Microbial ecologists / microbiome researchers** running 16S studies.
- **Lab members** who need to re-run or extend the pipeline.
- **Bioinformaticians** adapting the workflow to new markers or tools.
- **Reviewers / collaborators** auditing the analysis for publication.

## Expected Outcomes
- Cleaned, denoised ASV tables and representative sequences.
- Taxonomic classifications with confidence scores.
- Diversity and community-composition results suitable for publication.
- Full provenance: every output traceable to inputs, parameters, and versions.

## User Experience Goals
- A newcomer can run the pipeline on demo data following `README.md` alone.
- Errors are actionable and point to the failing rule and log file.
- Configuration lives in one obvious place (`config/`), not scattered in code.
- Adding a new tool means adding one rule + one env file.

## Domain Background (16S)
- **Amplicon:** 16S rRNA gene, typically hypervariable regions (V3–V4, V4).
- **Primers:** Region-specific; must be removed before denoising to avoid
  spurious ASVs.
- **ASVs vs OTUs:** ASVs (amplicon sequence variants) give single-nucleotide
  resolution and are preferred over 97% OTU clustering.
- **Denoising tools:** DADA2 (R), Deblur (Python), UNOISE3 (usearch/vsearch).
- **Reference DBs:** SILVA, Greengenes (incl. Greengenes2), GTDB, RDP.
- **Common confounders:** sequencing depth, primer bias, chimeras, contamination.