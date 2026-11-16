from pathlib import Path
import gzip,csv,json,hashlib,subprocess,concurrent.futures,time,shutil
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
RETAINED=ROOT/'secondary_analysis_r4'
O=ROOT/'work/r4_reproduction'
O.mkdir(parents=True,exist_ok=True)
R=ROOT
S=RETAINED/'resources';D=O/'protein_search';D.mkdir(exist_ok=True)
def fasta(p):
 with gzip.open(p,'rt') as f:
  name=None;seq=[]
  for l in f:
   if l.startswith('>'):
    if name:yield name,''.join(seq)
    name=l[1:].split()[0];seq=[]
   else:seq.append(l.strip().rstrip('*'))
  if name:yield name,''.join(seq)
maps={};summary={}
for species in ['DH1','TAIR']:
 with gzip.open(S/(species+'_feature_table.txt.gz'),'rt') as f:
  rows=list(csv.DictReader(f,delimiter='\t'))
 acc={r['product_accession']:r for r in rows if r['# feature']=='CDS' and r['product_accession'] and r['GeneID']}
 by={}
 for a,seq in fasta(S/(species+'_protein.faa.gz')):
  if a not in acc:continue
  row=acc[a];g=row['GeneID'];cur=by.get(g)
  if cur is None or (-len(seq),a)<(-len(cur['sequence']),cur['protein']):by[g]={'gene_id':g,'protein':a,'AGI':row['locus_tag'],'sequence':seq,'symbol':row['symbol']}
 groups={}
 for g,r in by.items():groups.setdefault(r['sequence'],[]).append(g)
 rep={};out=[]
 for seq,gs in groups.items():
  gs=sorted(gs,key=int);r=gs[0];rep[r]=gs
  out.append('>'+r+'\n'+seq+'\n')
  for g in gs:by[g]['representative']=r;by[g]['identical_gene_count']=len(gs)
 (D/(species+'.faa')).write_text(''.join(out),encoding='ascii')
 maps[species]={g:{k:v for k,v in r.items() if k!='sequence'} for g,r in by.items()}
 summary[species]={'protein_coding_genes':len(by),'unique_longest_proteins':len(groups),'identical_sequence_gene_groups':sum(len(v)>1 for v in rep.values()),'identical_sequence_genes':sum(len(v) for v in rep.values() if len(v)>1)}
 with (D/(species+'_gene_protein_map.tsv')).open('w',encoding='utf-8',newline='') as f:
  wr=csv.DictWriter(f,fieldnames=list(next(iter(maps[species].values()))),delimiter='\t');wr.writeheader();wr.writerows(maps[species].values())
(D/'PROTEIN_SUMMARY.json').write_text(json.dumps(summary,indent=2));print(summary,flush=True)
import os
blastp=shutil.which(os.environ.get('BLASTP','blastp')) or shutil.which('blastp.exe')
assert blastp, 'Install NCBI BLAST+ 2.17.0 and add its bin directory to PATH'
bin=Path(blastp).parent
def run(cmd,stem):
 print('START',stem,flush=True);start=time.time()
 cmd=[str(x.relative_to(D)) if isinstance(x,Path) and x.is_relative_to(D) else str(x) for x in cmd]
 with (D/(stem+'.stdout.log')).open('w') as out,(D/(stem+'.stderr.log')).open('w') as err:p=subprocess.run(cmd,stdout=out,stderr=err,cwd=D)
 assert p.returncode==0,(stem,p.returncode)
 print('DONE',stem,round(time.time()-start,2),flush=True)
 return {'step':stem,'command':[str(x) for x in cmd],'seconds':time.time()-start,'returncode':p.returncode}
records=[]
records.append(run([blastp,'-version'],'version'))
for sp in ['DH1','TAIR']:records.append(run([shutil.which('makeblastdb') or shutil.which('makeblastdb.exe'),'-in',D/(sp+'.faa'),'-dbtype','prot','-out',D/sp],sp+'_makeblastdb'))
commands=[]
for query,target in [('DH1','TAIR'),('TAIR','DH1')]:
 stem=query+'_vs_'+target;out=D/(stem+'.tsv')
 commands.append(([blastp,'-query',D/(query+'.faa'),'-db',D/target,'-out',out,'-evalue','1e-10','-seg','yes','-soft_masking','true','-max_hsps','1','-max_target_seqs','10','-num_threads','6','-outfmt','6 qseqid sseqid pident length qstart qend sstart send qlen slen evalue bitscore'],stem))
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:records.extend(ex.map(lambda x:run(*x),commands))
(D/'EXECUTION.json').write_text(json.dumps(records,indent=2));print('SEARCH_COMPLETE',flush=True)
