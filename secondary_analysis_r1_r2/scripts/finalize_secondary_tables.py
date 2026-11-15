"""Offline final table metadata corrections retained from finalize_r1.py. No fitting."""
from pathlib import Path
import json,collections
import pandas as pd
R=Path(__file__).resolve().parents[2];O=R/'work/r1_reproduction'
p=O/'R1_CANONICAL_REGULATOR_PANEL.tsv'
p.write_text(p.read_text(encoding='utf-8').replace('10.1105/tpc.13.12.2907','10.1104/pp.010324'),encoding='utf-8')
qc=pd.read_csv(O/'R1_MAPPING_BIAS_SAMPLE_QC.tsv',sep='\t');qsum=pd.read_csv(O/'R1_MAPPING_BIAS_QC.tsv',sep='\t')
for i,row in qsum.iterrows():
 f=qc[qc.genotype==row.genotype]
 if row.summary_level=='genotype_time':f=f[f.day==int(row.day)]
 counts=collections.defaultdict(collections.Counter)
 for s in f.fastqc_modules:
  for mate,mods in json.loads(s).items():
   for name,status in mods.items():counts[name][status]+=1
 qsum.at[i,'fastqc_summary']=json.dumps(counts,sort_keys=True)
qsum.to_csv(O/'R1_MAPPING_BIAS_QC.tsv',sep='\t',index=False)
print('Retained bibliographic correction and FastQC module grouping applied; no model or hypothesis test.')
