"""Offline checksum audit only. Does not run analyses or download resources."""
from pathlib import Path
import csv,hashlib
R=Path(__file__).resolve().parents[2]
with (R/'provenance/file_checksums.tsv').open(encoding='utf-8') as f:rows=list(csv.DictReader(f,delimiter='\t'))
for row in rows:
 p=R/row['relative_path'];assert p.is_file(),row['relative_path'];assert p.stat().st_size==int(row['size_bytes']),row['relative_path'];assert hashlib.sha256(p.read_bytes()).hexdigest()==row['sha256'],row['relative_path']
print('PASS:',len(rows),'candidate files; no statistical analysis executed.')
