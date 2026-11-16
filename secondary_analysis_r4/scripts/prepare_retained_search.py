"""Copy retained search evidence into ignored replay workspace. No search or analysis."""
from pathlib import Path
import shutil,hashlib,json
R=Path(__file__).resolve().parents[2];S=R/'secondary_analysis_r4/annotation/direct_sequence_mapping';D=R/'work/r4_reproduction/protein_search';D.mkdir(parents=True,exist_ok=True)
for n in ['DH1_gene_protein_map.tsv','TAIR_gene_protein_map.tsv','DH1_vs_TAIR.tsv','TAIR_vs_DH1.tsv','PROTEIN_SUMMARY.json']:
 p=D/n
 if p.exists():assert hashlib.sha256(p.read_bytes()).digest()==hashlib.sha256((S/n).read_bytes()).digest(),n
 else:shutil.copyfile(S/n,p)
(D/'EXECUTION.json').write_bytes((R/'secondary_analysis_r4/provenance/BLAST_COMMANDS.json').read_bytes())
print('Retained search checkpoint copied; no BLAST or statistical analysis executed.')
