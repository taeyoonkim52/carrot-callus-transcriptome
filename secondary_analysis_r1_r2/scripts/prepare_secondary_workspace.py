"""Seed retained inputs only. No analysis and no network."""
from pathlib import Path
import hashlib,shutil
R=Path(__file__).resolve().parents[2]
for n in ['primary_model.npz','technical_expression_unicode.npz']:
 a=R/'results/source_data/checkpoints'/n;b=R/'work/dedifferentiation'/n;b.parent.mkdir(parents=True,exist_ok=True)
 if b.exists():
  assert hashlib.sha256(a.read_bytes()).digest()==hashlib.sha256(b.read_bytes()).digest(),'Conflicting checkpoint'
 else:shutil.copyfile(a,b)
print('Retained checkpoints seeded; no analysis executed.')
