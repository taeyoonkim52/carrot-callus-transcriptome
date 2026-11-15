# R1 reference / mapping-bias QC

Decision: **NO_GROSS_SIGNAL**. Read-only review of existing Kallisto 0.50.1 unstranded quantification against DH1 v3.0, with full-file FastQC. No new alignment, quantification or raw-data download.

| genotype   |   libraries |   pseudoalignment_percent_mean |   pseudoalignment_percent_min |   pseudoalignment_percent_max |   estimated_fragments_mean |   estimated_fragments_min |   estimated_fragments_max |   detected_genes_min |   detected_genes_max |
|:-----------|------------:|-------------------------------:|------------------------------:|------------------------------:|---------------------------:|--------------------------:|--------------------------:|---------------------:|---------------------:|
| C815       |           6 |                        92.2333 |                          91.5 |                          92.9 |                2.78621e+07 |               2.23696e+07 |               3.33984e+07 |                28854 |                30105 |
| C819       |           6 |                        92.3333 |                          91.7 |                          92.9 |                2.61887e+07 |               2.15556e+07 |               2.94319e+07 |                26944 |                29156 |
| C824       |           6 |                        91.8667 |                          90.9 |                          92.8 |                2.52725e+07 |               1.95016e+07 |               3.01619e+07 |                25254 |                28735 |
| Ws         |           6 |                        91.85   |                          91.5 |                          92.1 |                2.74461e+07 |               2.48599e+07 |               3.24819e+07 |                29156 |                29701 |

Individual Ws libraries span 91.5–92.1%; the three cultivated accessions together span 90.9–92.9%. Ws is not grossly depleted in pseudoalignment. The review's Ws 91.7–92.0% and cultivated 91.6–92.8% correspond to rounded genotype×time averages, not the full individual-library ranges. Both aggregation levels are provided in the TSV.

Full-file FastQC module totals across 48 mates:

```json
{
  "Basic Statistics": {
    "PASS": 48
  },
  "Per base sequence quality": {
    "PASS": 48
  },
  "Per sequence quality scores": {
    "PASS": 48
  },
  "Per base sequence content": {
    "WARN": 48
  },
  "Per sequence GC content": {
    "PASS": 48
  },
  "Per base N content": {
    "PASS": 48
  },
  "Sequence Length Distribution": {
    "PASS": 48
  },
  "Sequence Duplication Levels": {
    "FAIL": 24,
    "WARN": 24
  },
  "Overrepresented sequences": {
    "PASS": 48
  },
  "Adapter Content": {
    "PASS": 48
  }
}
```

Duplicated sequences and base-composition warnings must be interpreted in RNA-seq context and are not proof of a Ws-specific failure. The shared warnings are retained, not discarded. No read trimming or deduplication was introduced.

Gross pseudoalignment efficiency was comparable among genotypes, providing no evidence that the Ws response pattern was driven simply by markedly poorer mapping efficiency. This does not eliminate gene-specific allele/reference bias, paralog ambiguity, ancestry confounding or tissue-composition differences. No gene-specific reference-bias test or alternate-genome experiment was performed.
