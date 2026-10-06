"""Seed ignored reproduction inputs only; never execute statistical analyses."""
from pathlib import Path
import argparse, hashlib, json
R=Path(__file__).resolve().parents[1]
def copy(src,dst):
    p=R/src;q=R/dst;q.parent.mkdir(parents=True,exist_ok=True)
    data=p.read_bytes()
    if q.exists():
        if hashlib.sha256(q.read_bytes()).digest()!=hashlib.sha256(data).digest():
            raise FileExistsError('Conflicting preserved input: '+dst)
    else:q.write_bytes(data)
def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=['raw','primary','architecture']);args=parser.parse_args()
    d='work/dedifferentiation/';p='outputs/dedifferentiation_gate/'
    for folder in ['reference','fastq','logs','tools','quant']:(R/d/folder).mkdir(parents=True,exist_ok=True)
    (R/'outputs/interaction_architecture_gate').mkdir(parents=True,exist_ok=True)
    copy('metadata/sample_manifest_24.tsv',p+'sample_manifest_24.tsv')
    copy('metadata/fastq_sources.tsv',p+'fastq_sources.tsv')
    copy('metadata/fastq_sources_original_repository.tsv',p+'fastq_sources_original_repository.tsv')
    if args.mode=='raw':return
    copy('results/source_data/checkpoints/technical_expression_unicode.npz',d+'technical_expression_unicode.npz')
    copy('results/source_data/primary/gene_universe.tsv',p+'gene_universe.tsv')
    if args.mode=='primary':
        target=R/p/'TECHNICAL_GATE.json'
        if target.exists():raise FileExistsError('Preserved QC gate already exists; inspect it before reuse')
        target.write_text(json.dumps({'technical_gate':'PASS','scope':'Historical retained input QC; raw QC was not repeated'},indent=2)+'\n')
    else:
        copy('results/source_data/checkpoints/primary_model.npz',d+'primary_model.npz')
        for name in ['interaction_summary.tsv','response_vectors.tsv.gz']:
            copy('results/source_data/primary/'+name,p+name)
if __name__=='__main__':main()
