# Carrot early callus transcriptome reproducibility release v1.2.0

This release retains the original 24-library Day0/Day20 analysis and adds final augmented homology-supported annotation and fixed-prior PyDESeq2 coefficient sensitivity. No new biological experiment, primary FDR or gene-set reselection occurred.

- v1.0.0: original primary analysis; commit 98cf5a3551a007b82ee1ed42a47c526c5568f70d; DOI https://doi.org/10.5281/zenodo.23180748.
- v1.1.0: earlier sensitivity, annotation and QC additions in secondary_analysis_r1_r2/; commit de216d9ef5b53316b9f2fbfc76037ad823ddea3b; DOI https://doi.org/10.5281/zenodo.24259913.
- v1.2.0: final complementary annotation and standard fixed-prior coefficient sensitivity in secondary_analysis_r4/. See that directory for scope, complete outputs, commands, limitations and copy-only integrity audit.

Concept DOI remains https://doi.org/10.5281/zenodo.23180747. The actual new version DOI will be recorded on the published GitHub release after Zenodo assigns it; no DOI is guessed or embedded as a self-referential release commit.

Primary results/, metadata/, figures/source_data/, original scripts/ and environment/ remain unchanged. The fixed universe is 24,961 genes, with 15,053 interactions and a 968-gene core (800 internally stable). New annotation has usable GO for 9,508 universe genes and 380 core genes. New coefficient sensitivity has 14,986 evaluable interactions, 67 all-zero exclusions and 69.4982% LARGE. A clearly identified reporting-field erratum preserves historic source bytes and does not alter scientific values.

Raw FASTQ/SRA files are not redistributed. Original data: Wang et al., BioProject PRJNA1398431. Read docs/analysis_overview.md for the historical primary workflow; read the two additional-analysis README files for their respective scopes. Historical directory/provenance names remain identifiable and do not imply prior journal review.

## Citation and roles

CITATION.cff describes software/archive creators Tae-Yoon Kim and Heesung Woo. Manuscript authors are Tae-Yoon Kim, Seung-Yop Lee and Heesung Woo. Seung-Yop Lee is already a Researcher contributor in Zenodo; confirmed manuscript roles are Conceptualization and Writing – review & editing. No software-production role or automatic creator promotion is inferred. See LICENSE_SCOPE.md for code-only MIT and third-party/data/documentation scope.
