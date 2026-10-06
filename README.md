# Carrot early callus transcriptome reproducibility archive

Authors: Tae-Yoon Kim and Heesung Woo.

Target manuscript: “Genotype-specific transcriptional reprogramming and a shared response core during early carrot callus induction”.

Publication: [MANUSCRIPT UNDER PREPARATION]. Contact: [TO BE FINALIZED BY PI]. This is a local public-release candidate awaiting author and license approval.

This study asks how shared transcriptional responses coexist with genotype-specific reprogramming across the Day0–Day20 callus-induction interval in C815, C819, C824 and Ws. The analysis supports the claim: “A stable shared transcriptional response core is embedded within broad genotype-specific reprogramming during the early carrot callus-induction interval.” It does not establish genotype-level regeneration competence or causal regeneration mechanisms.

Raw RNA-seq data are not redistributed in this repository.
The original sequencing data are publicly available under
NCBI BioProject PRJNA1398431.

See [the BioProject](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1398431), metadata/sample_manifest.tsv and metadata/accession_crosswalk.tsv. The selected design is 24 libraries: four genotype labels, two times, three dish replicates per cell. Replicate suffixes were not assumed to represent cross-time pairing.

## Contents

- metadata/: sample/accession and retrieval/checksum records.
- scripts/: acquisition, QC, quantification, aggregation, statistics, robustness, annotation and figure code.
- results/source_data/: derived primary/architecture tables and two numerical checkpoints.
- results/supplementary_tables/: S1–S9.
- figures/source_data/: retained figure inputs; future renders go to ignored figures/reproduced/.
- environment/: recorded Python package pins and external tool versions.
- provenance/: source-to-candidate hashes, edits, references, ontology and candidate checksums.
- docs/: execution order and limitations.

Scripts remain flat to preserve their original project-root resolution and cross-script calls. Stage mapping appears in docs/analysis_overview.md. Runtime work/ and outputs/ are ignored and contain no retained publication files.

## Reproduction

Use Python 3.12.14 and the recorded packages in environment/requirements.txt. Read environment/README_environment.md before installing external tools. No software, reference or raw sequencing binaries are bundled.

To render Figures 1–5 from retained inputs, run `python scripts/render_figures.py`. To reproduce primary statistics from the retained numeric checkpoint, first run `python scripts/prepare_workspace.py primary`, then `python scripts/dedif_primary_analysis.py`. For architecture and omission sensitivity, run `python scripts/prepare_workspace.py architecture`, then `python scripts/interaction_architecture.py`. Functional annotation additionally requires the pinned reference GAF and GO ontology; see docs/analysis_overview.md.

DH1 v3.0 assembly GCF_001625215.2 with RefSeq annotation GCF_001625215.2-RS_2024_03 was used. Exact retrieval URLs and historical checksums are in provenance/reference_resources.tsv. Reference files are not redistributed.

Expected outputs include response similarity and dispersion tables, the fixed 24,961-gene universe, interaction tests/effect sizes, omission sensitivity, the 968-gene shared core, stability and predefined GO summaries, and five figures. Retained values are unchanged; this archive preparation executed no analysis.

Citation: use CITATION.cff for the archive and separately credit the original PRJNA1398431 dataset/publication. No publication DOI or release DOI has been assigned here. LICENSE is an MIT candidate for project code/documentation only; LICENSE_REVIEW.md explains pending approval and separate data/resource rights. Do not interpret the candidate as blanket relicensing of all files.
