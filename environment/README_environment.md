# Recorded environment

software_versions.txt and requirements.txt are transcribed from frozen project records, not inferred versions. Python 3.12.14 was used. Create an isolated environment and install requirements.txt; external Kallisto 0.50.1, FastQC 0.12.1, Temurin Java 21.0.8+9 and SRA Toolkit 3.2.1 must be obtained independently. No installations or availability checks were performed here.

Native Windows was the production platform. Kallisto is expected at work/dedifferentiation/tools/kallisto.exe. FastQC/Java and SRA conversion layouts are explicit in their scripts. The original scripts prepend ignored runtime python_packages directories; absent directories permit normal installed-package imports. No virtual environment or vendor code is included.

Use the pinned reference/ontology checksums before executing statistical reproduction. Exact environment reconstruction and end-to-end numerical agreement remain untested in this static-only preparation. Model inputs/checkpoints, algorithms, seeds and thread limits are retained to support such validation.
