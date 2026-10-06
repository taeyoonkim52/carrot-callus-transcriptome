"""Read/quantification provenance checks and technical QC; no biological tests."""
import sys,json,csv,gzip,zipfile,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];W=ROOT/'work/dedifferentiation';O=ROOT/'outputs/dedifferentiation_gate'
sys.path.insert(0,str(W/'python_packages'))
import numpy as np,pandas as pd
from scipy.spatial.distance import pdist,squareform
samples=pd.read_csv(O/'sample_manifest_24.tsv',sep='\t',dtype={'day':str}).set_index('sample')
assert len(samples)==24 and set(samples.day)=={'0','20'}
def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for c in iter(lambda:f.read(8*1024*1024),b''):h.update(c)
    return h.hexdigest()
if any(not (W/'quant'/s/'COMPLETE.json').exists() for s in samples.index):raise SystemExit('PENDING: quantifications incomplete; no gate assigned')
files=pd.read_csv(O/'fastq_sources.tsv',sep='\t')
for name in files.filename:
    base=name.removesuffix('.gz').removesuffix('.fq').removesuffix('.fastq')
    if not (O/'fastqc'/(base+'_fastqc.zip')).exists():raise SystemExit('PENDING: full-file FastQC incomplete; no gate assigned')
mapping=pd.read_csv(O/'transcript_to_gene.tsv',sep='\t',dtype=str).set_index('transcript_id').gene_id
genes=sorted(mapping.unique());gid={g:i for i,g in enumerate(genes)};indices=np.array([gid[g] for g in mapping])
counts=[];tpms=[];lengths=[];qc=[];flags=[]
for s,row in samples.iterrows():
    target=W/'quant'/s;marker=json.loads((target/'COMPLETE.json').read_text())
    for filename,sha in marker['outputs'].items():
        if digest(target/filename)!=sha:raise RuntimeError('Completed output hash mismatch '+s)
    a=pd.read_csv(target/'abundance.tsv',sep='\t').set_index('target_id').reindex(mapping.index)
    if a.isna().any().any():raise RuntimeError('Missing transcript target '+s)
    ct=np.bincount(indices,weights=a.est_counts,minlength=len(genes));tpm=np.bincount(indices,weights=a.tpm,minlength=len(genes))
    weighted=np.bincount(indices,weights=a.tpm*a.eff_length,minlength=len(genes))
    fallback=np.bincount(indices,weights=a.eff_length,minlength=len(genes))/np.bincount(indices,minlength=len(genes))
    effective=np.divide(weighted,tpm,out=fallback.copy(),where=tpm>0)
    counts.append(ct);tpms.append(tpm);lengths.append(effective)
    raw=(target/'run_info.json').read_text()
    try:info=json.loads(raw)
    except json.JSONDecodeError:info=json.loads(raw.replace('\\','\\\\'))
    prefix=json.loads((W/'read_qc'/s/'read_qc.json').read_text())
    modules={};read_counts=[]
    for name in files.loc[files['sample']==s,'filename']:
        base=name.removesuffix('.gz').removesuffix('.fq').removesuffix('.fastq');archive=O/'fastqc'/(base+'_fastqc.zip')
        with zipfile.ZipFile(archive) as z:
            if z.testzip():raise RuntimeError('FastQC CRC failure '+name)
            summary=z.read(next(n for n in z.namelist() if n.endswith('/summary.txt'))).decode()
            modules[name]={line.split('\t')[1]:line.split('\t')[0] for line in summary.splitlines() if len(line.split('\t'))>=2}
            data=z.read(next(n for n in z.namelist() if n.endswith('/fastqc_data.txt'))).decode()
            read_counts.append(int(next(line.split('\t')[1] for line in data.splitlines() if line.startswith('Total Sequences\t'))))
    if read_counts[0]!=read_counts[1] or read_counts[0]!=info['n_processed']:flags.append(s+': full-file read counts inconsistent')
    if s=='C815_20d_3' and (O/'sra_fallback_provenance.json').exists() and info['n_processed']!=30104493:flags.append(s+': derived FASTQ spot count does not match repository 30104493 spots')
    if info['p_pseudoaligned']<50:flags.append(s+': pseudoalignment <50%')
    if min(prefix['mean_phred'])<20:flags.append(s+': mean sampled quality <20')
    for name,m in modules.items():
        if m.get('Per base sequence quality')=='FAIL':flags.append(s+': FastQC base quality failure '+name)
    qc.append(dict(sample=s,genotype=row.genotype,day=row.day,replicate=row.biological_replicate,processed_pairs=info['n_processed'],pseudoaligned_pairs=info['n_pseudoaligned'],pseudoalignment_percent=info['p_pseudoaligned'],estimated_fragments=ct.sum(),detected_genes=int((ct>0).sum()),orientation='unstranded',prefix_FR_percent=prefix['pilot_pseudoalignment_percent']['FR'],prefix_RF_percent=prefix['pilot_pseudoalignment_percent']['RF'],prefix_unstranded_percent=prefix['pilot_pseudoalignment_percent']['unstranded'],fastqc_modules=json.dumps(modules)))
counts=np.array(counts).T;tpms=np.array(tpms).T;lengths=np.array(lengths).T
scaled=tpms*lengths.mean(axis=1,keepdims=True);scaled=scaled/scaled.sum(axis=0,keepdims=True)*counts.sum(axis=0,keepdims=True)
frame=pd.DataFrame(scaled,index=genes,columns=samples.index);frame.index.name='gene_id'
frame.to_csv(O/'gene_expression_matrix.tsv.gz',sep='\t',compression={'method':'gzip','mtime':0})
pd.DataFrame(counts,index=genes,columns=samples.index).to_csv(O/'gene_estimated_counts.tsv.gz',sep='\t',compression={'method':'gzip','mtime':0})
pd.DataFrame(tpms,index=genes,columns=samples.index).to_csv(O/'gene_tpm.tsv.gz',sep='\t',compression={'method':'gzip','mtime':0})
pd.DataFrame(lengths,index=genes,columns=samples.index).to_csv(O/'gene_effective_lengths.tsv.gz',sep='\t',compression={'method':'gzip','mtime':0})
identical=pd.read_csv(O/'identical_sequence_groups.tsv',sep='\t')
ambiguous=set(g for ids in identical.loc[identical.cross_gene_ambiguity.astype(str).str.lower()=='true','gene_ids'] for g in ids.split(';'))
cpm=scaled/scaled.sum(axis=0,keepdims=True)*1e6
keep=(cpm>=1).sum(axis=1)>=3;keep&=~np.isin(genes,list(ambiguous));Y=np.log2(cpm[keep].T+1)
if Y.shape[1]<1000:flags.append('Too few testable unambiguous genes')
pd.DataFrame({'gene_id':genes,'cross_gene_identical':np.isin(genes,list(ambiguous)),'primary_testable':keep}).to_csv(O/'gene_universe.tsv',sep='\t',index=False)
dist=squareform(pdist(Y))/np.sqrt(Y.shape[1]);corr=np.corrcoef(Y)
pd.DataFrame(dist,index=samples.index,columns=samples.index).to_csv(O/'sample_distances.tsv',sep='\t')
pd.DataFrame(corr,index=samples.index,columns=samples.index).to_csv(O/'sample_correlations.tsv',sep='\t')
groups=samples.genotype+'_'+samples.day;cohesion=[]
for i,s in enumerate(samples.index):
    own=np.flatnonzero((groups==groups.iloc[i]).to_numpy());own=own[own!=i]
    own_dist=np.sqrt(np.mean((Y[i]-Y[own].mean(axis=0))**2));alternatives={g:np.sqrt(np.mean((Y[i]-Y[(groups==g).to_numpy()].mean(axis=0))**2)) for g in groups.unique() if g!=groups.iloc[i]}
    nearest=min(alternatives,key=alternatives.get);mincorr=float(corr[i,own].min())
    potential_swap=alternatives[nearest]<0.5*own_dist
    if mincorr<0.8:flags.append(s+': poor replicate correlation (<0.8)')
    if potential_swap:flags.append(s+': alternative group substantially closer; label review required')
    cohesion.append(dict(sample=s,own_leave_one_out_rms_distance=own_dist,nearest_alternative=nearest,nearest_alternative_distance=alternatives[nearest],minimum_replicate_correlation=mincorr,potential_label_problem=potential_swap))
pd.DataFrame(cohesion).to_csv(O/'replicate_cohesion.tsv',sep='\t',index=False)
pd.DataFrame(qc).to_csv(O/'sample_qc.tsv',sep='\t',index=False)
np.savez_compressed(W/'technical_expression_unicode.npz',Y=Y,genes=np.asarray(genes,dtype=str)[keep],scaled_counts=scaled[keep],sample_ids=np.asarray(samples.index,dtype=str),cpm=cpm[keep])
gate={'technical_gate':'FAIL' if flags else 'PASS','flags':flags,'samples':24,'fastqs':48,'primary_testable_genes':int(keep.sum()),'cross_gene_identical_genes':len(ambiguous),'scope':'Day0/Day20 only','read_qc':'FastQC full files plus 200000-pair orientation/quality pilots','limitations':['Expression-based cohesion cannot prove physical sample identity','Historical dish replication is reported, not independently observed','Unstranded quantification limits resolution of overlapping antisense transcripts','Single DH1 reference may create genotype-dependent mapping bias']}
gate['input_integrity_scope']='46 original FASTQs: repository MD5; 2 derived mates: verified full-quality SRA Normalized parent archive and local MD5/SHA256. Full read counts compared with quantification; derived run checked against repository total spots.'
(O/'TECHNICAL_GATE.json').write_text(json.dumps(gate,indent=2));print(json.dumps(gate,indent=2))
