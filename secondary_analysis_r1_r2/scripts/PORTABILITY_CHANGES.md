# Secondary script portability changes

The original calculation statements were retained; R3-A did not execute these scripts.
Original analyze_r1.py SHA256: 588ec87a6e93cf05bbd230338efb7769117d0fe6f9226fbe90be74e57b578945.
Original fisher_sensitivity_r1.py SHA256: de2fa8f588bda818c1a72f226d56a29f0c86c6695863d48fae40830aa5a14658.
Paths now resolve relative to the candidate; retained resources replace network acquisition;
outputs go to ignored work/r1_reproduction; a completed-output guard prevents accidental overwrite.
The original desktop-manuscript guard and local Git executable/history calls are removed.
Frozen historical criteria TSV is copied rather than inventing chronology in a new repository.
Primary-model/Unicode checkpoints are seeded from the unchanged primary archive.
The original initial HMH approximation is still run first during a FUTURE full reproduction,
then Fisher sensitivity replaces only its designated summaries, exactly as documented in the addendum.
The reference-check/final-report scripts are not biological result generators and are not needed here.
Only AST syntax validation was performed in R3-A; no end-to-end numerical reproduction is claimed.
