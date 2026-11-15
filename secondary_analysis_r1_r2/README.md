# Expert review secondary analyses supporting R2

These are retained R1 outputs integrated into R2, prepared for an unreleased v1.1.0. They supplement the frozen primary analysis; they do not refit it, replace FDR values or select a new hypothesis.

| Directory | Claim supported | Inputs and output | Limitation |
|---|---|---|---|
| mapping_qc | No gross Ws pseudoalignment deficit | Retained 24-library QC → genotype/time summary | Aggregate efficiency cannot exclude allelic/gene-specific bias |
| core_annotation | Fixed 968 genes, 441 up/527 down | S7 + RefSeq + mapped GO → complete core table | Names can be uncharacterised; membership unchanged |
| canonical_regulators | Qualified regulator matches and unresolved cases | Predeclared 23 Arabidopsis loci + Compara + S4/S7 → all panel records | Predicted orthology is not functional equivalence; unresolved is not absence |
| ortholog_mapping | 6,192 high-confidence mappings; 6,075 GO genes | Frozen BioMart snapshots and explicit NCBI bridge | Older ASM162521v1 vs DH1v3; ambiguous/one-to-many retained |
| go_enrichment | Ribosome/core-down and ubiquitin-ligase/core-up annotation | Restricted mapped background + frozen ontology → all 11,065 tests | 235/968 core GO coverage; global BH; ancestry/qualifier limits; no activity inference |
| effect_size_sensitivity | LARGE percentage support dependence | Saved NB and CPM checkpoints + S4 → final Fisher/CPM tables | Approximate normal/Fisher, not refit/new FDR; CPM can exclude on/off genes |
| numerical_qc | Gene108213928 precision warning | Saved fitted information → two warning records | No targeted refit; exact extreme range/Wald precision not persuasive |
| criteria_provenance | Operational rather than preregistered criteria | Frozen plan files + original Git commit chronology | First committed with results; pre-inspection fixation unresolved |
| figures | Figure5 label consistency | Unchanged S7/S8 and original plot workflow → new PNG/PDF/SVG | Historical plot files remain unchanged |

## Offline independent reproduction instructions

Commands below are instructions for a future authorised independent reproduction, NOT actions executed in R3-A. Run from the candidate root in a fresh copy. Install Python3.12.14 and requirements-secondary.txt in an independent environment. Full primary plotting needs environment/requirements.txt. Retained resources and checkpoints are already present; no FASTQ/SRA or current BioMart download is needed for these derived-only secondary results.

1. python secondary_analysis_r1_r2/scripts/audit_release.py — integrity checks only.
2. python secondary_analysis_r1_r2/scripts/prepare_secondary_workspace.py — seed two unchanged numerical checkpoints into ignored work/dedifferentiation/; no fitting.
3. python secondary_analysis_r1_r2/scripts/analyze_r1_portable.py — future secondary mapping/GO/panel/support and initial HMH computation; writes only work/r1_reproduction/. Reads original S4/S7, sample metadata, RefSeq feature/GAF, frozen resources and saved checkpoints. The initial HMH output is a retained diagnostic, not the final sensitivity.
4. python secondary_analysis_r1_r2/scripts/fisher_sensitivity_r1_portable.py — final expected-Fisher sensitivity/numerical warning, keeping original fit/FDR fixed. Read effect_size_sensitivity/R1_SHRINKAGE_METHOD_ADDENDUM.md first.
5. python secondary_analysis_r1_r2/scripts/finalize_secondary_tables.py — future offline metadata finalization from the original finalizer: correct only the recorded SERK1 reference DOI and combine FastQC module outcomes across mates exactly as the authoritative final R1 tables. No model/hypothesis test.
6. Compare numerical tables against retained outputs. Set/hash iteration order or platform floating point can change byte order/low digits; compare keyed records and documented numerical tolerance rather than claiming all independent executions produce identical TSV hashes. The recorded source-to-candidate byte copies themselves MUST match exactly.

The portability edits were AST syntax checked only, not executed end to end. Historical Git criteria evidence is carried forward, not reconstructed from this candidate's Git history. No repository needed for the audit itself. No network script is invoked in this order. The finalizer's two table-changing metadata steps are included; its narrative-report generation and live DOI checks are not required to generate biological outputs. The retained environment manifest lists the historical run order, while the commands above govern the portable candidate.

## Versions and provenance

See provenance/R1_ENVIRONMENT_PUBLIC.json, requirements-secondary.txt, original environment/requirements.txt, resources/resource_manifest.json and provenance/R1_ANNOTATION_PROVENANCE.json. Ensembl Plants63 (June2026 archive metadata), Arabidopsis TAIR10; project DH1v3/RefSeqRS_2024_03; GO snapshot SHA256 and file-changing commit are pinned. Resource snapshot hashes govern replay, not mutable web endpoints.
The canonical-panel predeclaration's incorrect Hecht DOI is retained as history and corrected in the final panel; see its correction note. Annotation qualifiers/dates are incomplete in BioMart, three roots/obsolete terms are excluded, is_a/part_of ancestry only, global BH covers 2,213 terms × five queries.

supplementary_tables/ preserves exact R2 S10–S17 bytes and all sub-tables; provenance/SUPPLEMENT_ALIGNMENT.tsv links them to R1 and verifies the R2 ZIP. S1–S9 remain under results/supplementary_tables/ unchanged. Identical named R1 and S-table copies are deliberate traceability views, not separate analyses.
