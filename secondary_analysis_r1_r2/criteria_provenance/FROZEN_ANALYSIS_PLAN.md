# Frozen dedifferentiation analysis plan

Frozen before new expression results. Scope: C815, C819, C824, Ws; days 0 and 20; three independent dish replicates per cell. No late libraries or differentiation phenotype values enter this analysis.

Reference: DH1 v3.0, GCF_001625215.2, RefSeq GCF_001625215.2-RS_2024_03 (2024-03-27). Quantify the full annotated RNA transcriptome, aggregate by stable GeneID with transcript effective-length normalization. Kallisto 0.50.1, paired reads, unstranded default because strand metadata are not established; retain raw reads without trimming. Record reference MD5 and SHA256, input MD5/SHA256 and output SHA256. A strand diagnostic and read-quality inspection must precede acceptance of production quantifications.

Technical gate: all 24 libraries and checksum-verified mates; complete quantifications; no catastrophic read/library failure. Pseudoalignment <50% triggers investigation rather than automatic exclusion. Replicate cohesion and potential label problems assessed without changing metadata based only on expression. Unresolved major failure stops biological testing.

Count inference: use a validated RNA-seq negative-binomial implementation, full genotype * time model with three interaction degrees of freedom; balanced average time contrast and four within-genotype contrasts; BH FDR 0.05. Freeze the exact implementation and length-offset method before inference; do not substitute ordinary least squares for count inference.

Universe: CPM >=1 in at least three libraries, without DEG selection. Transform normalized gene abundance using log2(CPM+1) for descriptive vectors and distances. All response similarities use this entire fixed universe. Report Pearson correlation, cosine, direction agreement, RMS response and 2,000 independent dish bootstrap resamples within each genotype/time cell (seed 20261004). Time samples are not paired.

Conserved core: all four estimated within-genotype log2 fold changes have the same sign, minimum absolute effect >=0.5, maximum minus minimum fold change <=1 log2 unit, and balanced time contrast FDR <=0.05. Report interaction significance separately; non-significance does not establish equivalence. Report sensitivity summaries at effect thresholds 0.25 and 1 without changing the primary core.

Convergence: mean squared Euclidean separation of genotype centroids per gene, averaged over six pairs, at each time. Also report sample-level distances and centroid distances corrected for within-cell sampling variance. Bootstrap biological dishes within eight cells for dispersion ratio and difference, 2,000 resamples; convergence supported only when the 95% ratio interval is below 1 and primary reduction >=20%. No PCA-only conclusion. Small n and bootstrap uncertainty remain limitations.

Models: A supported if all six response correlations have positive bootstrap lower bounds, median correlation >=0.5, and core >=10% of testable genes; partial if directional concordance exists but these thresholds are not all met. B supported by the frozen quantitative convergence criterion; partial if point reduction exists without full uncertainty support. C supported if interaction FDR genes >=10% of universe and >=50% have a between-genotype response range >=1 log2 unit; partial if significant interaction exists below these thresholds. Otherwise not supported. These are operational thresholds, not universal biological laws.

Program analysis only after primary results are saved: ten requested categories, defined by GO terms/descendants before enrichment; fixed universe, BH across tested categories. Photosynthesis suppression requires negative core effects. No post hoc attractive categories or candidate selection.

Publication gate depends jointly on technical validity, interpretable replication, clear predefined model support, substantive novelty beyond published PCA/selected contrasts, and independence from unresolved late phenotype. Original regeneration NO_GO remains unchanged.
