# Final additional annotation and coefficient sensitivity

This directory contains the retained final homology-supported annotation and fixed-prior PyDESeq2 sensitivity supporting the submitted Day0-to-Day20 carrot manuscript. No new primary analysis or biological experiment was performed during release preparation.

## Version scope

v1.0.0 is the original primary analysis. v1.1.0 adds the earlier sensitivity, annotation and quality-control materials retained under secondary_analysis_r1_r2/. v1.2.0 adds the final complementary annotation and coefficient sensitivity here. Historical commits/tags/DOIs are retained; primary results, FDR, gene sets and core membership are unchanged.

## Retained evidence

- annotation/direct_sequence_mapping: complete bidirectional DH1 v3.0/TAIR10.1 search hits, qualifying hits, gene/protein/identical-sequence groups and reference gene space. No BLAST database cache or technical probe is included.
- annotation/augmented_annotation and coverage_comparison: old/direct/augmented mappings, discordances, complete coverage. Direct-only annotation did not improve usable GO coverage; augmented annotation combines the prior high-confidence layer and strict reciprocal homology only for previously unmapped genes. This is not a new gene-tree orthology framework.
- functional_enrichment: complete five-fixed-set family, 13,925 tests with one global BH correction, nulls, contributing genes, full old/new comparisons. Core-down ribosome is 29/212 versus 188/9,508, 6.918205539943798-fold, FDR 1.1325509472926264e-13. Core-up ubiquitin-ligase complex is 10/168 versus 150/9,508, 3.773015873015873-fold, FDR 0.018993775667778955. Persistence is annotation sensitivity, not independent biological replication.
- canonical_regulators: complete 23-regulator reassessment, one existing LBD29-associated candidate with strict reciprocal homology, 17 ambiguous and five unresolved outcomes. No named carrot ortholog or gene-absence claim.
- lfc_shrinkage: original 15,053 interaction records, four genotype-specific outputs, 67 all-zero-block exclusions and unchanged full-precision summary. 14,986 evaluable; 601 SMALL, 3,970 MODERATE, 10,415 LARGE (69.49819831843054%), compared with 70.34565594554918% unshrunk on the same subset. Primary 70.48% remains primary. This is fixed-prior PyDESeq2 coefficient sensitivity conditional on saved size factors/dispersions, not adaptive R apeglm, a new interaction FDR or a joint posterior for range.
- numerical_qc: explicit reporting erratum. The historical summary field nonzero_optimizer_failures counted nonzero blocks owing to object-dtype Boolean inversion; all actual nonzero convergence flags are True, failures zero. The original summary is preserved. No scientific coefficients/statistics changed.
- resources/provenance/checksums: pinned source resources, exact normalized commands, software hashes, script adaptation ledger and source-matching manifest. NCBI BLAST+ binaries are retrieved separately using the pinned manifest, not redistributed.

Usable GO is 9,508/24,961 (38.09%) and 380/968 core (39.26%); 394 core genes map. Noncoding genes lacking proteins and ambiguous homology remain unresolved.

## Future independent reproduction

These instructions are not analyses executed during release packaging. Use a fresh checkout, Python 3.12.14 and requirements.txt. Install NCBI BLAST+ 2.17.0 on PATH. Versioned input hashes govern replay, rather than mutable current resource downloads.

1. Run python secondary_analysis_r4/scripts/audit_release.py for hashes only.
2. Run python secondary_analysis_r1_r2/scripts/prepare_secondary_workspace.py to copy existing frozen primary checkpoints to ignored work/dedifferentiation/. No fit is performed by this preparation step.
3. Run python secondary_analysis_r4/scripts/run_protein_search.py only for a separately authorized full BLAST replay. It reads pinned resources and writes work/r4_reproduction/protein_search/. To reuse retained searches instead, copy the five retained search/mapping files and PROTEIN_SUMMARY.json to that directory and use the command provenance as the completed-search sentinel EXECUTION.json. The preparation helper below performs this byte-copy operation without running a search.
4. Run python secondary_analysis_r4/scripts/prepare_retained_search.py for the copy-only route, then compare_annotation.py and refine_coverage_reporting.py in that order for an authorized derived-only replay. They use exact existing baseline GO/Compara snapshots, fixed S4 genes and the predeclared panel; no primary model is rerun.
5. run_standard_shrinkage.py can reproduce the four conditional coefficient fits from the saved primary checkpoint, writing only work/r4_reproduction/. It does not call deseq2(), summary(), refit dispersions or replace primary FDR. Existing generated per-genotype files are checkpoints; validate them before reusing.
6. report_optimizer_status.py creates a separate reporting-only corrected summary from retained flags; it preserves the historical source summary. It does not compute any new coefficient, statistic or association.

Portable adapters were syntax checked, not re-executed end to end for release preparation. Scientific byte copies match authoritative sources; a future computational replay can differ in row ordering/floating-point last digits and should compare keyed outputs with justified tolerance. All primary/data provenance remains in the unchanged baseline repository. Do not interpret unresolved homology as absence or annotation enrichment as causal activity.
