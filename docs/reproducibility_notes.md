# Static audit scope and limitations

No biological analysis, quantification, figure rendering, package installation or network download was executed to prepare this archive. Static checks validate Python syntax, manifest design/accession consistency, table structure, retained source hashes, gzip/ZIP integrity, documented runtime input contracts and package records. They do not establish successful end-to-end execution on a new computer.

The acquisition/quantification scripts target the original native Windows execution layout. Other operating systems require an independently reviewed launcher adaptation. Install Kallisto, Java, FastQC and SRA Toolkit externally; no binary is bundled. Exact pins are recorded, not asserted to be available indefinitely. Numeric libraries, parallel computation and future public data changes can affect bitwise reproducibility. Model checkpoints allow reproduction without downloading raw sequencing data. Reference, index and ontology hashes identify the original versions.

The public script edits are enumerated in provenance/script_edits.tsv. Removed blocks were administrative checks, superseded figure/report generation and unreported threshold-sensitivity outputs; statistical function bodies, filtering/model parameters and seeds were preserved. A new standard-library workspace seeder handles public layout inputs without executing analyses. Full-file gene IDs are stored as Unicode rather than pickled object arrays, preserving the already retained Unicode checkpoint.

All retained derived tables and numerical checkpoints have byte-identical source provenance. Sanitized software/provenance JSON removes personal executable locations and administrative commands; reference/tool hashes remain. Numerical results were not recomputed or altered.

provenance/file_checksums.tsv covers every archive file except itself and .git/. Excluded raw/private files are not checksummed there. Source repository history was not copied into this archive.

Original project analysis code is licensed under MIT; see LICENSE_SCOPE.md for the distinct scope of source data, derived outputs and external resources.
