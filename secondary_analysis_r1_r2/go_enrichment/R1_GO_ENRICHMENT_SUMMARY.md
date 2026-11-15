# Secondary ortholog-transferred GO layer — PARTIAL

This supplements, and does not replace, the frozen 12-category/36-test analysis. The original chromatin result remains unchanged and does not establish chromatin activity or causation.

Method and exact resource versions: R1_ANNOTATION_PROVENANCE.json; raw queries and hashes: sources/resource_manifest.json; complete bridge table: R1_ORTHOLOG_MAPPING.tsv. Retrieval 2026-11-15. Ensembl Plants 63 uses older ASM162521v1, not project DH1 v3.0. A unique explicit NCBI GeneID↔Ensembl bridge is required in BOTH directions; high-confidence Compara ortholog records only. No loose gene-name mapping. Multiple qualifying Arabidopsis orthologs are unioned; genes are counted once per term. All mappings and rejected statuses are retained.

Coverage: 6,192/24,961 genes have a usable high-confidence bridge; 6,075/24,961 (24.34%) have non-root GO. Core coverage is 246/968 mapped and 235/968 (24.28%) GO-annotated. There are 15,089 no-explicit-xref genes, 259 ambiguous-bridge genes and 3,421 uniquely bridged genes without a high-confidence ortholog. Thus annotation **depth** improves for a restricted subset; whole-universe coverage does not improve over the original 18,409-gene local-GAF background. This restriction is substantial and non-random.

All BP/MF/CC terms with 5–5,000 testable background genes were tested (2,213 terms × five sets = 11,065 tests), using one-sided hypergeometric enrichment and global BH. is_a/part_of ancestors were propagated. Obsolete/unknown terms and roots were excluded. All tested null terms, raw p, global FDR and every contributing GeneID are in R1_GO_ENRICHMENT_FULL.tsv. These five sets overlap; GO ancestry makes many significant terms redundant. Counts of significant terms are not counts of independent biological programs.

| set                   |   tested_terms |   testable_set_size |   original_set_size |   FDR_le_0_05 |
|:----------------------|---------------:|--------------------:|--------------------:|--------------:|
| CORE_ALL              |           2213 |                 235 |                 968 |            16 |
| CORE_DOWN             |           2213 |                 140 |                 527 |            28 |
| CORE_UP               |           2213 |                  95 |                 441 |             4 |
| LARGE_INTERACTION     |           2213 |                2587 |               10609 |            93 |
| PRACTICAL_INTERACTION |           2213 |                3805 |               14486 |            47 |

## Shared core, all

| GO         | term                               |   gene_count |   testable_set_size |   background_count |   fold_enrichment |   global_FDR |
|:-----------|:-----------------------------------|-------------:|--------------------:|-------------------:|------------------:|-------------:|
| GO:0005840 | ribosome                           |           26 |                 235 |                148 |           4.5414  |  2.6864e-08  |
| GO:0003735 | structural constituent of ribosome |           24 |                 235 |                127 |           4.88524 |  2.84123e-08 |
| GO:0022625 | cytosolic large ribosomal subunit  |           12 |                 235 |                 36 |           8.61702 |  1.39463e-06 |
| GO:0022626 | cytosolic ribosome                 |           17 |                 235 |                 84 |           5.23176 |  3.91772e-06 |
| GO:0032991 | protein-containing complex         |           56 |                 235 |                733 |           1.97498 |  5.39371e-05 |
| GO:0006412 | translation                        |           25 |                 235 |                211 |           3.06292 |  9.13106e-05 |

## Core down

| GO         | term                               |   gene_count |   testable_set_size |   background_count |   fold_enrichment |   global_FDR |
|:-----------|:-----------------------------------|-------------:|--------------------:|-------------------:|------------------:|-------------:|
| GO:0005840 | ribosome                           |           26 |                 140 |                148 |           7.62307 |  8.84713e-13 |
| GO:0003735 | structural constituent of ribosome |           24 |                 140 |                127 |           8.20022 |  1.71431e-12 |
| GO:0032991 | protein-containing complex         |           51 |                 140 |                733 |           3.01915 |  7.78962e-11 |
| GO:0022626 | cytosolic ribosome                 |           17 |                 140 |                 84 |           8.78189 |  3.47625e-09 |
| GO:0022625 | cytosolic large ribosomal subunit  |           12 |                 140 |                 36 |          14.4643  |  8.61447e-09 |
| GO:0006412 | translation                        |           24 |                 140 |                211 |           4.93568 |  2.6864e-08  |

## Core up

| GO         | term                                                                    |   gene_count |   testable_set_size |   background_count |   fold_enrichment |   global_FDR |
|:-----------|:------------------------------------------------------------------------|-------------:|--------------------:|-------------------:|------------------:|-------------:|
| GO:0000151 | ubiquitin ligase complex                                                |            8 |                  95 |                 99 |           5.16746 |    0.0125336 |
| GO:0031461 | cullin-RING ubiquitin ligase complex                                    |            7 |                  95 |                 82 |           5.45892 |    0.0201486 |
| GO:0031146 | SCF-dependent proteasomal ubiquitin-dependent protein catabolic process |            4 |                  95 |                 27 |           9.47368 |    0.0452009 |
| GO:0009416 | response to light stimulus                                              |           13 |                  95 |                302 |           2.7527  |    0.0461184 |

The mapped down-core includes ribosomal/translation-associated genes; the mapped up-core includes ubiquitin-ligase/SCF-related annotations. These describe annotation associations in the mapped subset. They do not establish changes in translation rate, protein degradation flux, cell identity or regeneration competence. No claim that all other pathways are absent is justified.

Limitations: cross-version bridges discard many genes; paralog transfer and one-to-many union can broaden inferred functions; annotation ascertainment bias remains; BioMart does not expose annotation dates/NOT qualifiers, so this is not a fully qualifier-audited TAIR GAF; all evidence codes including electronic evidence were retained. Contemporary refers to the retrieved version, not newly validated carrot function. GO master ontology content is pinned by SHA256. KEGG was not forced.
