# Carrot early callus transcriptome reproducibility release v1.2.0

This release retains the original 24-library Day0/Day20 analysis and adds final augmented homology-supported annotation and fixed-prior PyDESeq2 coefficient sensitivity. No new biological experiment, primary FDR or gene-set reselection occurred.

- v1.0.0: original primary analysis; commit 98cf5a3551a007b82ee1ed42a47c526c5568f70d; DOI https://doi.org/10.5281/zenodo.23180748.
- v1.1.0: earlier sensitivity, annotation and QC additions in secondary_analysis_r1_r2/; commit de216d9ef5b53316b9f2fbfc76037ad823ddea3b; DOI https://doi.org/10.5281/zenodo.24259913.
- v1.2.0: final complementary annotation and standard fixed-prior coefficient sensitivity in secondary_analysis_r4/. See that directory for scope, complete outputs, commands, limitations and copy-only integrity audit.

The final v1.2.0 scientific archive is available at https://doi.org/10.5281/zenodo.24293704. Concept DOI https://doi.org/10.5281/zenodo.23180747 links archive versions. The frozen v1.2.0 tag and archived file snapshot are unchanged by current main-branch metadata corrections.

Primary results/, metadata/, figures/source_data/, original scripts/ and environment/ remain unchanged. The fixed universe is 24,961 genes, with 15,053 interactions and a 968-gene core (800 internally stable). New annotation has usable GO for 9,508 universe genes and 380 core genes. New coefficient sensitivity has 14,986 evaluable interactions, 67 all-zero exclusions and 69.4982% LARGE. A clearly identified reporting-field erratum preserves historic source bytes and does not alter scientific values.

Raw FASTQ/SRA files are not redistributed. Original data: Wang et al., BioProject PRJNA1398431. Read docs/analysis_overview.md for the historical primary workflow; read the two additional-analysis README files for their respective scopes. Historical directory/provenance names remain identifiable and do not imply prior journal review.

## Citation and roles

Manuscript authors, in order: Tae-Yoon Kim, Taejin Kim and Seung-Yop Lee. Seung-Yop Lee is the corresponding author (sylee@sogang.ac.kr).

The archive/software creator is Tae-Yoon Kim. CITATION.cff cites this analysis resource rather than automatically treating every manuscript author as a software author. Archive contributors, in order, are Taejin Kim (Methodology; DataCite Other) and Seung-Yop Lee (Conceptualization, Writing – review & editing, Supervision; DataCite Supervisor). The controlled DataCite roles are recorded in .zenodo.json, with the full contribution basis in its notes.

Tae-Yoon Kim's confirmed manuscript contributions are Conceptualization, Formal analysis, Validation, Biological interpretation, Writing – review & editing, and Supervision. These manuscript roles are distinct from the archive creator designation; no superseded manuscript roles are restored.

Historical citation/provenance files and frozen releases retain the metadata recorded at release time; current metadata does not rewrite those snapshots. Adding .zenodo.json on main does not by itself change an existing Zenodo record. See LICENSE_SCOPE.md for code-only MIT and third-party/data/documentation scope.
