# Public-safety audit

PASS — no prohibited private material or credential detected in selected candidate files.

SAFE: accession metadata, project-authored reproducibility scripts, byte-identical derived tables/checkpoints, figure inputs, S1–S9 and sanitized provenance. Scanned 85 pre-audit files, including full decompressed gzip content, NPZ headers and complete string-array contents; numeric payloads were verified by CRC and source hash rather than interpreted as text. Python syntax/matrix multiplication and public URL syntax were manually distinguished from email/path scan false positives. All file names and selected scripts were reviewed; no contact email, telephone/home address, journal login, private URL, private correspondence or conversation/prompt content was retained.

REDACTED: local executable locations and administrative command strings in selected provenance; unnecessary machine-specific Kallisto source path. Script changes appear in provenance/script_edits.tsv.

REMOVED: raw FASTQ/SRA, vendor binaries/packages, caches, source Git history, unpublished manuscript/cover letter, correspondence, author declarations, internal PI/gate/audit reports, rejected/superseded helpers and unrelated project materials. The source inventory lists relative excluded paths without reproducing their content.

PI_REVIEW_REQUIRED: license/data rights, final authors/contact and release approval. This does not represent a detected private-file or credential exception. No material is published by this audit.

Coverage includes drive-qualified and user-home paths, personal email strings, credential assignments, key/token patterns, private-key headers, archive member content and manual file-role inspection. This is a bounded static review, not a guarantee against every possible undisclosed secret. Final staged-file scan must repeat before commit.
