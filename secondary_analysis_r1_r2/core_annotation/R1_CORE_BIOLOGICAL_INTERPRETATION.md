# Complete fixed core and canonical panel — MODERATE interpretability

All 968 original GeneIDs are retained: 441 up and 527 down; 800 retain the original robustness designation. All 968 have a local NCBI feature-table name (some remain uncharacterized LOC names); 747 have local GO annotation. The complete table carries all original statistics and all four effects, bootstrap direction support, omission retention, RefSeq identifiers and orthology statuses. No genes were reselected and no core thresholds were relaxed.

The canonical panel was fixed in R1_PANEL_PREDECLARED.tsv before joining to memberships. Stable Arabidopsis loci were cross-checked against curated BioMart records and primary accession lists (WIND: doi:10.1111/nph.17594; PLT: https://pmc.ncbi.nlm.nih.gov/articles/PMC5240741/; WOX: doi:10.1105/tpc.114.122887). Hecht's initially mistyped DOI was corrected transparently in R1_PANEL_BIBLIOGRAPHIC_CORRECTION.md; panel loci and selection did not change. SERK2–5 are family comparators, not experimentally equivalent embryogenesis regulators.

Resolved predicted homologs: LBD16-like GeneID108213528 and LBD29-like GeneID108196780 are INTERACTION ONLY; WIND4-associated GeneID108210944 is NEITHER. Each has MODERATE confidence because the Compara relationship is one-to-many, despite confidence=1. Low-confidence LBD16/WOX11/BBM/WIND4 records are explicitly UNRESOLVED, not positive matches. Most other requested loci have no qualifying bridge in this expressed universe. This is not evidence of their biological absence, no expression, or non-core status. No HIGH-confidence one-to-one named canonical core member was established by this method. Missing panel mappings must not be converted to NEITHER.

Selected display of the down-core ribosome annotation group (first eight GeneIDs numerically, not chosen by effect or a desired story; all contributors in full enrichment table):

|   gene_id | carrot_symbol   | carrot_name                                                       | core_direction   |
|----------:|:----------------|:------------------------------------------------------------------|:-----------------|
| 108192452 | LOC108192452    | large ribosomal subunit protein eL32z                             | DOWN             |
| 108192518 | LOC108192518    | small ribosomal subunit protein eS26y                             | DOWN             |
| 108193285 | LOC108193285    | large ribosomal subunit protein eL28y                             | DOWN             |
| 108193745 | LOC108193745    | uncharacterized LOC108193745;uncharacterized protein LOC108193745 | DOWN             |
| 108194160 | LOC108194160    | large ribosomal subunit protein eL42                              | DOWN             |
| 108195143 | LOC108195143    | large ribosomal subunit protein eL39x                             | DOWN             |
| 108196740 | LOC108196740    | small ribosomal subunit protein eS7                               | DOWN             |
| 108197444 | LOC108197444    | large ribosomal subunit protein eL14                              | DOWN             |

The shared core contains common directional bulk responses with a substantial housekeeping/protein-homeostasis annotation component. It cannot be relabeled 'regeneration genes' or a proven embryogenesis program. Ribosome-associated down-core enrichment and up-core SCF/ubiquitin-ligase annotations give a more concrete biological description than an anonymous core count. Their functional implications remain hypotheses, and restricted orthology coverage prevents a complete program-level interpretation. The original local-GAF chromatin association is retained separately; chromatin activity and causality were not measured.

R2 should expose the complete core and canonical negative/unresolved results, rather than select attractive individual genes. An optional future sequence-based homology audit could resolve missing mappings, but was not silently substituted for the declared Compara method.
