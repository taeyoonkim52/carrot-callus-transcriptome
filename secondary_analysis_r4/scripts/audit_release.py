"""Offline hash audit only; no scientific computation."""
from pathlib import Path
import hashlib,csv
R=Path(__file__).resolve().parents[2]
with (R/'secondary_analysis_r4/checksums/V1_2_0_MANIFEST.tsv').open(encoding='utf-8') as f:rows=list(csv.DictReader(f,delimiter='\t'))
for r in rows:
 p=R/r['path'];assert p.stat().st_size==int(r['size']) and hashlib.sha256(p.read_bytes()).hexdigest()==r['SHA256'],r['path']
print('PASS',len(rows),'scientific source copies; no analysis performed.')
