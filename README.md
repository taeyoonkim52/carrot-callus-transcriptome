# Carrot early callus transcriptome reproducibility release v1.1.0

Version v1.1.0. Adds expert-review secondary sensitivity, functional-annotation, mapping-QC and provenance analyses supporting the revised manuscript. The original primary analysis is unchanged. This release contains no new biological experiment.

## Primary analysis

The retained results/, metadata/, figures/source_data/, original scripts/ and primary environment/ describe the original v1.0.0 scientific analysis: 24 Day0/Day20 libraries, fixed 24,961-gene universe, genotype-dependent responses and a smaller stable shared core. Raw sequencing files are not redistributed. Original dataset: Wang et al., BioProject PRJNA1398431.
The frozen primary archive identifies commit 98cf5a3551a007b82ee1ed42a47c526c5568f70d and version DOI https://doi.org/10.5281/zenodo.23180748. It is not rewritten or replaced. The concept DOI is https://doi.org/10.5281/zenodo.23180747.

## Secondary expert review analyses

secondary_analysis_r1_r2/ adds mapping QC, complete fixed-core annotation, canonical regulators, cross-version orthology and exploratory GO, approximate expected-Fisher and expression-support sensitivities, numerical warnings and criteria provenance. These support the R2 manuscript and were not part of v1.0.0. See its README for execution order, inputs, outputs and limitations. The original primary numerical results remain unchanged; the arithmetic display 14,486/24,961 is 58.03% rather than 58.04%. A new Figure5 copy changes only the operational-category heading; historical figures/code remain retained.

## Citation and contributorship

CITATION.cff describes the software/analysis archive, retaining historical archive creators Tae-Yoon Kim and Heesung Woo. Manuscript authors are Tae-Yoon Kim, Seung-Yop Lee and Heesung Woo; manuscript authorship and archive production are distinct. Seung-Yop Lee's documented contribution is Conceptualization; no software-production role is inferred. Cite the frozen v1.0.0 version for historical primary results and v1.1.0 for these secondary analyses. The actual version DOI is reported on the GitHub v1.1.0 release and under the existing Zenodo concept record after archival publication; no version DOI is guessed here.

Original project code remains within the MIT scope in LICENSE_SCOPE.md. External reference/GO/Ensembl/publication resources retain their own terms; see secondary_analysis_r1_r2/resources/RESOURCE_TERMS.md. Corporate affiliation and competing-interest wording require factual author confirmation before submission; this release does not invent a company-role declaration.

Read docs/analysis_overview.md for PRIMARY reproduction only; it preserves historical operational labels. Read secondary_analysis_r1_r2/README.md for additions. Audit checksums offline with python secondary_analysis_r1_r2/scripts/audit_release.py. No statistical analysis was run during R3-A preparation.
