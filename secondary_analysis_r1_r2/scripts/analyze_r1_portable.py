"""Secondary R1 evidence audit. Read-only inputs; writes exclusively beside this script."""
from pathlib import Path
import os,sys,json,hashlib,gzip,collections,functools,subprocess,xml.etree.ElementTree as ET,zipfile,datetime
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='1'
# Use an independent environment with requirements-secondary.txt.
import numpy as np,pandas as pd
from scipy.stats import hypergeom
R=Path(__file__).resolve().parents[2]; O=R/'work/r1_reproduction'; O.mkdir(parents=True,exist_ok=True); S=R/'secondary_analysis_r1_r2/resources'; REF=S
if (O/'R1_CORE_968_COMPLETE.tsv').exists(): raise FileExistsError('Completed secondary output exists; use a fresh reproduction workspace')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(df,name):df.to_csv(O/name,sep='\t',index=False)
def md(name,s):(O/name).write_text(s,encoding='utf-8')
def bh(p):
 p=np.asarray(p);order=np.argsort(p);v=p[order]*len(p)/np.arange(1,len(p)+1);q=np.empty(len(p));q[order]=np.minimum(1,np.minimum.accumulate(v[::-1])[::-1]);return q
# Protected release inputs are checked by audit_release.py; all outputs go to work/r1_reproduction.
guard={}
primary=pd.read_csv(R/'results/supplementary_tables/S4_interaction_statistics.tsv',sep='\t',dtype={'gene_id':str}).set_index('gene_id')
core=pd.read_csv(R/'results/supplementary_tables/S7_conserved_core.tsv',sep='\t',dtype={'gene_id':str}).set_index('gene_id')
assert len(primary)==24961 and len(core)==968 and (primary.interaction_FDR<=.05).sum()==15053
G=['C815','C819','C824','Ws'];universe=set(primary.index);sig=set(primary.index[primary.interaction_FDR<=.05])
practical=set(primary.index[(primary.interaction_FDR<=.05)&(primary.response_effect_range_log2>=.5)])
large=set(primary.index[(primary.interaction_FDR<=.05)&(primary.response_effect_range_log2>=1)])
assert len(practical)==14486 and len(large)==10609
print('Input assertions passed',flush=True)
# Exact biological context from retained primary article, independently re-read.
tree=ET.parse(R/'secondary_analysis_r1_r2/resources/Wang_2026_article.xml');sections=[]
for sec in tree.findall('.//sec'):
 title=sec.findtext('title','')
 if title in ['Sample collection','Data Description']:sections.append(title+'\n'+' '.join(sec.itertext()))
md('R1_ORIGINAL_ARTICLE_EVIDENCE.txt','\n\n'.join(sections))
rows=[]
for g,pct in zip(G,[25,12.5,2,0]):
 rows.append(dict(genotype=g,source_classification='wild Daucus carota subsp. carota accession' if g=='Ws' else 'local cultivated Daucus carota var. sativus accession; named cultivar/breeding-line status not specified',reported_differentiation_percent=pct,phenotype_evidence_status='Reported value VERIFIED; scoring endpoint, denominator and phenotype replication UNRESOLVED',Day0_material='hypocotyl explant from 10-day sterile seedling',Day20_material='whole cultured tissues after callus induction; not demonstrated pure cell type',early_medium='MS + 1 mg/L 2,4-D + 1 mg/L 6-BA; 26 C; 16h light/8h dark',later_stage_context='Day20-70 transferred induced calli to hormone-free MS for shoot regeneration; Day120 sampled; Day70-120 transfer details unspecified',RNA_replication='3 entire tissues pooled per Petri dish; 3 independent dishes per genotype/time; not evidence of phenotype n=3',source='Wang et al. 2026 doi:10.1016/j.dib.2026.112728; sections 3 and 4.1; retained work/article.xml',Day0_Day20_rationale='SUPPORTED for culture interval; does not establish identical cell composition or intervening developmental path'))
write(pd.DataFrame(rows),'R1_DATASET_CONTEXT.tsv')
# Mapping audit: existing whole-file QC, no new quantification.
qc=pd.read_csv(R/'secondary_analysis_r1_r2/mapping_qc/R1_MAPPING_BIAS_SAMPLE_QC.tsv',sep='\t');write(qc,'R1_MAPPING_BIAS_SAMPLE_QC.tsv')
numeric=['pseudoalignment_percent','processed_pairs','pseudoaligned_pairs','estimated_fragments','detected_genes']
qrows=[]
for level,groups in [('genotype_time',qc.groupby(['genotype','day'])),('genotype',qc.groupby('genotype'))]:
 for key,f in groups:
  key=key if isinstance(key,tuple) else (key,)
  row={'summary_level':level,'genotype':key[0],'day':key[1] if len(key)>1 else 'ALL','libraries':len(f)}
  for n in numeric:
   for stat in ['mean','min','max']:row[n+'_'+stat]=getattr(f[n],stat)()
  modules=collections.defaultdict(collections.Counter)
  for x in f.fastqc_modules:
   j=json.loads(x)
   for mate,vals in j.items():
    if isinstance(vals,dict):
     for name,status in vals.items():modules[mate+':'+name][status]+=1
  row['fastqc_summary']=json.dumps(modules,sort_keys=True);qrows.append(row)
qsum=pd.DataFrame(qrows);write(qsum,'R1_MAPPING_BIAS_QC.tsv')
md('R1_MAPPING_BIAS_QC.md','# Reference / pseudoalignment bias audit\n\nExisting 24-library Kallisto and full-file FastQC results only. DH1 v3.0 reference unchanged.\n\n'+qsum[qsum.summary_level=='genotype'][['genotype','libraries','pseudoalignment_percent_mean','pseudoalignment_percent_min','pseudoalignment_percent_max']].to_markdown(index=False)+'\n\nComparable aggregate pseudoalignment efficiency can exclude only a simple gross mapping-loss explanation. It cannot exclude gene-specific allelic reference bias, paralog ambiguity, wild/cultivated ancestry or composition effects. No genotype-level biological inference or mapping re-run was performed. Full module outcomes are retained in sample and genotype-time TSVs.\n')
# Local annotation, complete fixed core.
ft=pd.read_csv(REF/'GCF_001625215.2_DH1_v3.0_feature_table.txt.gz',sep='\t',dtype=str).fillna('');ft=ft.rename(columns={'# feature':'feature'})
annotation={}
for gid,f in ft.groupby('GeneID'):
 annotation[gid]={k:';'.join(sorted(set(f[c])-{''})) for k,c in [('carrot_symbol','symbol'),('carrot_name','name'),('RefSeq_product_accessions','product_accession'),('RefSeq_related_accessions','related_accession')]}
localgo=collections.defaultdict(set)
with gzip.open(REF/'GCF_001625215.2-RS_2024_03_gene_ontology.gaf.gz','rt') as f:
 for line in f:
  if line.startswith('!'):continue
  a=line.rstrip('\n').split('\t')
  if len(a)>4 and 'NOT' not in a[3].split('|'):localgo[a[1]].add(a[4])
# Parse contemporary hierarchy; is_a and part_of only. Ignore regulates relations.
obo=S/'go-edit.obo';terms={};aliases={}
for block in obo.read_text(encoding='utf-8').split('\n[Term]\n')[1:]:
 lines=block.splitlines();d={};parents=[];alts=[]
 for line in lines:
  if line.startswith('id: '):d['id']=line[4:]
  elif line.startswith('name: '):d['name']=line[6:]
  elif line.startswith('namespace: '):d['namespace']=line[11:]
  elif line.startswith('is_a: '):parents.append(line[6:].split()[0])
  elif line.startswith('relationship: part_of '):parents.append(line.split()[2])
  elif line.startswith('alt_id: '):alts.append(line[8:])
  elif line=='is_obsolete: true':d['obsolete']=True
 if 'id' in d:
  d['parents']=parents;terms[d['id']]=d
  for alt in alts:aliases[alt]=d['id']
@functools.lru_cache(None)
def ancestors(t):
 t=aliases.get(t,t)
 if t not in terms or terms[t].get('obsolete'):return frozenset()
 return frozenset([t]).union(*(ancestors(p) for p in terms[t]['parents']))
xref=pd.read_csv(S/'carrot_xref.tsv',sep='\t',names=['ensembl_gene','gene_id'],dtype=str).dropna().drop_duplicates()
ncbi2ens=xref.groupby('gene_id').ensembl_gene.agg(set).to_dict();ens2ncbi=xref.groupby('ensembl_gene').gene_id.agg(set).to_dict()
orth=pd.read_csv(S/'carrot_orthology.tsv',sep='\t',names=['ensembl_gene','arabidopsis_id','arabidopsis_symbol','orthology_type','target_percent_identity','query_percent_identity','confidence'],dtype=str).fillna('').drop_duplicates()
orthgroups={k:f for k,f in orth.groupby('ensembl_gene')}
ago=pd.read_csv(S/'arabidopsis_go.tsv',sep='\t',names=['arabidopsis_id','symbol','go_id','go_name','namespace','evidence'],dtype=str).fillna('')
agomap=collections.defaultdict(set)
for a in ago.itertuples(index=False):
 if a.go_id.startswith('GO:'):
  agomap[a.arabidopsis_id].update(ancestors(a.go_id))
roots={'GO:0008150','GO:0003674','GO:0005575'}
mapping=[];usable={};alllinks=[]
for gid in sorted(universe):
 es=ncbi2ens.get(gid,set());status='UNMAPPED_NO_EXPLICIT_NCBI_XREF';ats=set();good=[]
 if es:
  if len(es)!=1 or any(len(ens2ncbi[e])!=1 for e in es):status='AMBIGUOUS_CROSS_VERSION_BRIDGE'
  else:
   e=next(iter(es));f=orthgroups.get(e)
   if f is not None:
    for a in f.itertuples(index=False):
     if a.arabidopsis_id:alllinks.append(dict(gene_id=gid,**a._asdict()))
    good=f[(f.confidence=='1')&(f.arabidopsis_id!='')].to_dict('records')
   ats={a['arabidopsis_id'] for a in good};status='USABLE_HIGH_CONFIDENCE' if ats else 'NO_HIGH_CONFIDENCE_ORTHOLOG'
   if ats:usable[gid]=ats
 gos=set().union(*(agomap[a] for a in ats))-roots if ats else set()
 mapping.append(dict(gene_id=gid,ensembl_gene_ids=';'.join(sorted(es)),bridge_status=status,arabidopsis_ids=';'.join(sorted(ats)),orthology_types=';'.join(sorted({a['orthology_type'] for a in good})),one_to_many=len(ats)>1,transferred_GO_count=len(gos),transferred_GO=';'.join(sorted(gos))))
mapping=pd.DataFrame(mapping);write(mapping,'R1_ORTHOLOG_MAPPING.tsv');write(pd.DataFrame(alllinks),'R1_ORTHOLOG_ALL_BRIDGED_RECORDS.tsv')
gene2go={gid:set().union(*(agomap[a] for a in ats))-roots for gid,ats in usable.items()};background={g for g,ts in gene2go.items() if ts};term2genes=collections.defaultdict(set)
for g in background:
 for t in gene2go[g]:term2genes[t].add(g)
tested={t:gs for t,gs in term2genes.items() if 5<=len(gs)<=5000}
sets={'CORE_ALL':set(core.index),'CORE_UP':set(core.index[core.balanced_time_log2FC>0]),'CORE_DOWN':set(core.index[core.balanced_time_log2FC<0]),'PRACTICAL_INTERACTION':practical,'LARGE_INTERACTION':large}
enr=[]
for name,selected in sets.items():
 selected=selected&background;n=len(selected);N=len(background)
 for t,gs in sorted(tested.items()):
  hit=selected&gs;k=len(hit);K=len(gs);p=hypergeom.sf(k-1,N,K,n) if n else 1.
  enr.append(dict(set=name,GO=t,term=terms[t].get('name',''),namespace=terms[t].get('namespace',''),gene_count=k,testable_set_size=n,original_set_size=len(sets[name]),background_count=K,background_size=N,fold_enrichment=(k/n)/(K/N) if n else np.nan,raw_p=p,contributing_genes=';'.join(sorted(hit))))
enr=pd.DataFrame(enr);enr['global_FDR']=bh(enr.raw_p);write(enr,'R1_GO_ENRICHMENT_FULL.tsv')
coverage={'conceptual_universe':24961,'bridge_status_counts':mapping.bridge_status.value_counts().to_dict(),'high_confidence_mapped':len(usable),'GO_background':len(background),'tested_terms_per_set':len(tested),'tests_total':len(enr),'core_mapped':len(set(core.index)&set(usable)),'core_GO_annotated':len(set(core.index)&background),'significant_counts':enr[enr.global_FDR<=.05].groupby('set').size().to_dict()}
(O/'R1_ANNOTATION_SUMMARY.json').write_text(json.dumps(coverage,indent=2))
fullcore=core.copy()
for c in ['carrot_symbol','carrot_name','RefSeq_product_accessions','RefSeq_related_accessions']:fullcore[c]=[annotation.get(g,{}).get(c,'') for g in core.index]
fullcore['existing_GO']=[';'.join(sorted(localgo[g])) for g in core.index];fullcore['core_direction']=np.where(core.balanced_time_log2FC>0,'UP','DOWN')
fullcore=fullcore.join(mapping.set_index('gene_id'));write(fullcore.reset_index(),'R1_CORE_968_COMPLETE.tsv')
# Canonical panel was fixed on disk before any joins to memberships.
panel=pd.read_csv(R/'secondary_analysis_r1_r2/canonical_regulators/R1_PANEL_PREDECLARED.tsv',sep='\t');links=pd.DataFrame(alllinks);prows=[]
for a in panel.itertuples(index=False):
 matches=links[links.arabidopsis_id==a.arabidopsis_id] if len(links) else pd.DataFrame()
 if matches.empty:
  prows.append(dict(regulator=a.regulator,arabidopsis_id=a.arabidopsis_id,primary_DOI=a.primary_DOI,carrot_gene_id='',mapping_confidence='UNRESOLVED',membership='UNRESOLVED',reason='No unique explicit NCBI bridge to a Compara ortholog in expressed universe; not evidence of absence from carrot'));continue
 for b in matches.itertuples(index=False):
  g=b.gene_id;high=b.confidence=='1';incore=g in core.index;insig=g in sig
  membership=('BOTH CORE AND INTERACTION' if incore and insig else 'SHARED CORE' if incore else 'INTERACTION ONLY' if insig else 'NEITHER') if high else 'UNRESOLVED'
  conf='HIGH' if high and b.orthology_type=='ortholog_one2one' else 'MODERATE' if high else 'LOW'
  row=dict(regulator=a.regulator,arabidopsis_id=a.arabidopsis_id,primary_DOI=a.primary_DOI,carrot_gene_id=g,carrot_symbol=annotation.get(g,{}).get('carrot_symbol',''),carrot_name=annotation.get(g,{}).get('carrot_name',''),mapping_confidence=conf,Compara_confidence=b.confidence,orthology_type=b.orthology_type,query_percent_identity=b.query_percent_identity,target_percent_identity=b.target_percent_identity,membership=membership,reason='Predicted orthology, not demonstrated regulatory functional equivalence')
  row.update(primary.loc[g,[k+'_log2FC' for k in G]+['interaction_FDR','balanced_time_FDR']].to_dict());prows.append(row)
panelout=pd.DataFrame(prows);write(panelout,'R1_CANONICAL_REGULATOR_PANEL.tsv')
print('Annotation complete',coverage,flush=True)
# Sensitivity: saved NB point fit/covariance; no primary refit or replacement tests.
d=np.load(R/'work/dedifferentiation/primary_model.npz');tech=np.load(R/'work/dedifferentiation/technical_expression_unicode.npz')
assert d['genes'].tolist()==primary.index.tolist()==tech['genes'].tolist()
X=d['design'];names=d['design_columns'].tolist();C=np.zeros((4,8));C[:,names.index('time[T.D20]')]=1
for i,g in enumerate(G[1:],1):C[i,names.index('genotype[T.'+g+']:time[T.D20]')]=1
b=d['log_coefficients']@C.T/np.log(2);orig=primary[[g+'_log2FC' for g in G]].to_numpy();assert np.allclose(b,orig,atol=1e-8)
idx=np.flatnonzero(primary.interaction_FDR.to_numpy()<=.05);rawrange=np.ptp(b[idx],axis=1);srows=[];detail=[]
def summarize(name,e,indices):
 rg=np.ptp(e,axis=1);small=int((rg<.5).sum());moderate=int(((rg>=.5)&(rg<1)).sum());large_=int((rg>=1).sum())
 row=dict(analysis=name,retained_interaction_genes=len(rg),SMALL=small,MODERATE=moderate,LARGE=large_,LARGE_percent=100*large_/len(rg),practical_genes=moderate+large_,practical_percent_original_universe=100*(moderate+large_)/24961,median_response_range=np.median(rg),Q25=np.quantile(rg,.25),Q75=np.quantile(rg,.75),P90=np.quantile(rg,.90),P95=np.quantile(rg,.95),P99=np.quantile(rg,.99),maximum=rg.max(),above5=int((rg>5).sum()),above10=int((rg>10).sum()),above15=int((rg>15).sum()),above20=int((rg>20).sum()));srows.append(row)
 for j,v,eff in zip(indices,rg,e):detail.append(dict(analysis=name,gene_id=d['genes'][j],response_range=v,**dict(zip([g+'_log2FC' for g in G],eff.tolist()))))
summarize('ORIGINAL_UNSHRUNK',b[idx],idx)
shrunk={1.:[],2.:[]}
for j in idx:
 beta=d['log_coefficients'][j];mu=d['size_factors']*np.exp(X@beta);w=mu/(1+mu*d['dispersions'][j]);M=(X.T*w)@X;H=np.linalg.inv(M+np.eye(8)*1e-6);cov=C@H@M@H@C.T/np.log(2)**2
 for tau in shrunk:shrunk[tau].append(tau*tau*np.linalg.solve(cov+np.eye(4)*tau*tau,b[j]))
for tau,e in shrunk.items():summarize('NORMAL_APPROX_SHRINK_TAU_'+str(tau),np.asarray(e),idx)
meta=pd.read_csv(R/'metadata/sample_manifest_24.tsv',sep='\t',dtype={'day':str}).set_index('sample').loc[tech['sample_ids']]
for threshold in [1,5]:
 supported=np.ones(24961,dtype=bool)
 for g in G:
  for t in ['0','20']:
   inds=np.flatnonzero(((meta.genotype==g)&(meta.day==t)).to_numpy());assert len(inds)==3;supported&=(tech['cpm'][:,inds]>=threshold).sum(axis=1)>=2
 retain=idx[supported[idx]];summarize('CPM_'+str(threshold)+'_2of3_EVERY_CELL',b[retain],retain)
sensitivity=pd.DataFrame(srows);write(sensitivity,'R1_LFC_SENSITIVITY.tsv');write(pd.DataFrame(detail),'R1_LFC_GENE_DETAILS.tsv')
robust='ROBUST' if ((sensitivity.iloc[1:].practical_percent_original_universe>=10)&(sensitivity.iloc[1:].LARGE_percent>=50)).all() else 'PARTIALLY_ROBUST' if (sensitivity.iloc[1:].practical_percent_original_universe>=10).all() else 'NOT_ROBUST'
(O/'R1_SENSITIVITY_GATE.json').write_text(json.dumps({'classification':robust,'extreme_above10':int((rawrange>10).sum()),'max':float(rawrange.max())},indent=2))
md('R1_LFC_SENSITIVITY_REPORT.md','# LFC sensitivity — '+robust+'\n\nOriginal significant interaction membership and FDR remain fixed; no new inferential tests.\n\n'+sensitivity.to_markdown(index=False)+'\n\nA: normal-likelihood multivariate shrinkage of four responses, using saved NB HMH covariance. Independent N(0,tau²) priors at tau=1 and 2 log2 units. This is an approximate sensitivity, not apeglm/PyDESeq2 shrinkage, and does not re-estimate dispersion or produce replacement FDR. Multi-genotype max-minus-min is nonlinear; single-coefficient shrinkage was not misapplied. Normal approximation is weakest for sparse/on-off genes; prior dependence is explicit.\n\nB: CPM>=1 and companion CPM>=5 in at least two of three dishes in EVERY genotype/time cell. Original effects retained. This excludes genuine on/off genes and prioritizes well-supported effects; the subset is not a new universe or a phenotype-defined subset.\n\nThe classification rule was fixed before these secondary calculations in R1_SECONDARY_PROTOCOL.md. Primary descriptive response correlations/dispersion and omission fits are preserved and are not replaced by these checks.\n')
# Historical Git chronology is retained evidence, not recreated from a new checkout.
import shutil
shutil.copyfile(R/'secondary_analysis_r1_r2/criteria_provenance/R1_PREDEFINED_CRITERIA_PROVENANCE.tsv',O/'R1_PREDEFINED_CRITERIA_PROVENANCE.tsv')
# Verify protected original files stayed unchanged.
bad=[p for p,h in guard.items() if sha(Path(p))!=h];assert not bad,bad
(O/'R1_INTEGRITY_CHECK.json').write_text(json.dumps({'protected_files':len(guard),'changed_original_files':bad,'primary_assertions':'PASS','core_rows':len(fullcore),'gene_universe':len(primary),'timestamp':datetime.datetime.now().astimezone().isoformat()},indent=2))
print('R1 numeric work complete',robust,flush=True);print(sensitivity.to_string(index=False),flush=True)
