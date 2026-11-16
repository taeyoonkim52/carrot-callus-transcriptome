"""Validate existing convergence flags and correct reporting metadata only. No fit or statistics."""
from pathlib import Path
import csv,json
ROOT=Path(__file__).resolve().parents[2]
SRC=ROOT/'secondary_analysis_r4/lfc_shrinkage'
OUT=ROOT/'work/r4_reproduction';OUT.mkdir(parents=True,exist_ok=True)
summary=json.loads((SRC/'STANDARD_SHRINKAGE_SUMMARY.json').read_text(encoding='utf-8'))
corrected={}
for g in ['C815','C819','C824','Ws']:
 with (SRC/f'shrinkage_{g}.tsv').open(encoding='utf-8') as f:rr=list(csv.DictReader(f,delimiter='\t'))
 corrected[g]=sum(r['nonzero_genotype'].lower()=='true' and r['converged'].lower()!='true' for r in rr)
summary['historical_nonzero_optimizer_failures']=summary['nonzero_optimizer_failures']
summary['nonzero_optimizer_failures']=corrected
summary['reporting_note']='Explicit Boolean-string parser fixes object-dtype inversion. All scientific values retained; historical summary and erratum preserved.'
(OUT/'STANDARD_SHRINKAGE_SUMMARY_REPORTING_CORRECTED.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
print(corrected)
