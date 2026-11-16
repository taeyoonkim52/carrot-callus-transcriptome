from pathlib import Path
import os,csv,json,collections,functools
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
RETAINED=ROOT/'secondary_analysis_r4'
O=ROOT/'work/r4_reproduction'
O.mkdir(parents=True,exist_ok=True)
R=ROOT
A=ROOT/'secondary_analysis_r1_r2';D=O/'protein_search'
import numpy as np,pandas as pd
from scipy.stats import hypergeom
assert (D/'EXECUTION.json').exists(),'Search not complete'
maps={s:pd.read_csv(D/(s+'_gene_protein_map.tsv'),sep='\t',dtype=str).fillna('').set_index('gene_id') for s in ['DH1','TAIR']}
groups={s:{k:list(x.index) for k,x in m.groupby('representative')} for s,m in maps.items()}
cols=['query','target','identity','length','qstart','qend','sstart','send','qlen','slen','evalue','bitscore']
hits={};bests={};decisive={}
for q,t in [('DH1','TAIR'),('TAIR','DH1')]:
 f=pd.read_csv(D/f'{q}_vs_{t}.tsv',sep='\t',names=cols,dtype={'query':str,'target':str})
 raw_limits={gid:(len(ff),float(ff.bitscore.min())) for gid,ff in f.groupby('query')}
 f['query_coverage']=(f.qend-f.qstart+1)/f.qlen;f['target_coverage']=(f.send-f.sstart+1)/f.slen
 f=f[(f.identity>=30)&(f.query_coverage>=.7)&(f.target_coverage>=.7)&(f.evalue<=1e-10)].sort_values(['query','bitscore','target'],ascending=[True,False,True])
 f.to_csv(O/f'R4A_{q}_QUALIFYING_HITS.tsv',sep='\t',index=False)
 hits[q]=f;bests[q]={};decisive[q]={}
 for gid,rows in f.groupby('query'):
  b=rows.iloc[0];bests[q][gid]=str(b.target)
  truncated_near_tie=raw_limits[gid][0]>=10 and raw_limits[gid][1]>=float(b.bitscore)/1.1
  decisive[q][gid]=(len(rows)==1 or float(b.bitscore)>=1.1*float(rows.iloc[1].bitscore)) and not truncated_near_tie and len(groups[q][gid])==1 and len(groups[t][str(b.target)])==1
strict={};status={}
for g,row in maps['DH1'].iterrows():
 rep=row.representative;tar=bests['DH1'].get(rep);s='NO_QUALIFYING_HIT'
 if int(row.identical_gene_count)>1:s='IDENTICAL_SEQUENCE_GENE_AMBIGUITY'
 elif tar:
  if not decisive['DH1'][rep]:s='NEAR_TIE_OR_TARGET_SEQUENCE_AMBIGUITY'
  elif bests['TAIR'].get(tar)!=rep:s='NONRECIPROCAL_HOMOLOGY'
  elif not decisive['TAIR'].get(tar,False):s='REVERSE_NEAR_TIE_OR_SEQUENCE_AMBIGUITY'
  else:
   agi=maps['TAIR'].loc[tar,'AGI']
   if agi:strict[g]={agi};s='STRICT_RECIPROCAL_HOMOLOGY'
   else:s='NO_EXPLICIT_AGI'
 status[g]=s
# Reuse the exact old GO layer/hierarchy to isolate the impact of the protein bridge.
terms={};aliases={}
for block in (A/'resources/go-edit.obo').read_text(encoding='utf-8').split('\n[Term]\n')[1:]:
 d={};parents=[];alts=[]
 for line in block.splitlines():
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
ago=pd.read_csv(A/'resources/arabidopsis_go.tsv',sep='\t',names=['AGI','symbol','GO','term','namespace','evidence'],dtype=str).fillna('');agogo=collections.defaultdict(set)
for agi,go in ago[['AGI','GO']].itertuples(index=False,name=None):
 if go.startswith('GO:'):agogo[agi].update(ancestors(go))
roots={'GO:0008150','GO:0003674','GO:0005575'}
for agi in agogo:agogo[agi]-=roots
original=pd.read_csv(ROOT/'results/supplementary_tables/S4_interaction_statistics.tsv',sep='\t',dtype={'gene_id':str}).set_index('gene_id');universe=set(original.index)
old=pd.read_csv(A/'ortholog_mapping/R1_ORTHOLOG_MAPPING.tsv',sep='\t',dtype=str).fillna('').set_index('gene_id')
oldmap={g:set(row.arabidopsis_ids.split(';')) for g,row in old.iterrows() if row.bridge_status=='USABLE_HIGH_CONFIDENCE'}
layers={'OLD':oldmap,'DIRECT_RBH':{g:a for g,a in strict.items() if g in universe}}
layers['AUGMENTED_CANDIDATE']=dict(oldmap)
layers['AUGMENTED_CANDIDATE'].update({g:a for g,a in layers['DIRECT_RBH'].items() if g not in oldmap})
mapping=[]
for g in original.index:
 direct=layers['DIRECT_RBH'].get(g,set());previous=oldmap.get(g,set());new=layers['AUGMENTED_CANDIDATE'].get(g,set())
 mapping.append({'gene_id':g,'protein_status':status.get(g,'NO_DH1_PROTEIN'),'DH1_protein':maps['DH1'].loc[g,'protein'] if g in maps['DH1'].index else '', 'OLD_AGI':';'.join(sorted(previous)),'DIRECT_RBH_AGI':';'.join(sorted(direct)),'direct_old_discordance':bool(direct and previous and not direct.issubset(previous)),'candidate_AGI':';'.join(sorted(new)),'candidate_GO':';'.join(sorted(set().union(*(agogo[a] for a in new)) if new else set()))})
pd.DataFrame(mapping).to_csv(O/'R4A_GENE_ANNOTATION_CANDIDATE.tsv',sep='\t',index=False)
core=set(original.index[original.conserved_core]);assert len(core)==968
sets={'UNIVERSE':universe,'CORE_ALL':core,'CORE_UP':set(original.index[original.conserved_core&(original.balanced_time_log2FC>0)]),'CORE_DOWN':set(original.index[original.conserved_core&(original.balanced_time_log2FC<0)]),'PRACTICAL_INTERACTION':set(original.index[(original.interaction_FDR<=.05)&(original.response_effect_range_log2>=.5)]),'LARGE_INTERACTION':set(original.index[(original.interaction_FDR<=.05)&(original.response_effect_range_log2>=1)])}
assert [len(sets[n]) for n in ['CORE_UP','CORE_DOWN','PRACTICAL_INTERACTION','LARGE_INTERACTION']]==[441,527,14486,10609]
panel=pd.read_csv(A/'canonical_regulators/R1_PANEL_PREDECLARED.tsv',sep='\t',dtype=str).fillna('');oldpanel=pd.read_csv(A/'canonical_regulators/R1_CANONICAL_REGULATOR_PANEL.tsv',sep='\t',dtype=str).fillna('');pan=[]
for row in panel.itertuples(index=False):
 agi=row.arabidopsis_id;found=[g for g,ids in layers['DIRECT_RBH'].items() if agi in ids]
 targets=set(maps['TAIR'].index[maps['TAIR'].AGI==agi]);reps=set(maps['TAIR'].loc[list(targets),'representative']) if targets else set()
 homolog=hits['DH1'][hits['DH1'].target.isin(reps)]
 family=set(g for rep in homolog['query'] for g in groups['DH1'][rep])&universe
 prev=oldpanel[oldpanel.regulator==row.regulator]
 pan.append({'regulator':row.regulator,'Arabidopsis_AGI':agi,'OLD_confidence':';'.join(prev.mapping_confidence.unique()),'OLD_carrot_gene_ids':';'.join(sorted(set(prev.carrot_gene_id)-{''})),'direct_RBH_candidates':';'.join(sorted(found)),'qualifying_family_homology_candidates':';'.join(sorted(family)),'resolution':'RECIPROCAL_HOMOLOGY_CANDIDATE_NOT_VALIDATED_ORTHOLOG' if found else 'AMBIGUOUS_FAMILY_HOMOLOGY' if family else 'UNRESOLVED','core_candidates':';'.join(sorted(set(found)&core)),'interaction_candidates':';'.join(sorted(set(found)&set(original.index[original.interaction_FDR<=.05]))),'phylogeny':'NOT_PERFORMED; no specific paralog orthology claimed','limitation':'RBH is sequence evidence, not gene-tree/species-tree orthology or functional validation; incomplete hits do not establish absence.'})
pd.DataFrame(pan).to_csv(O/'R4A_CANONICAL_REGULATOR_REANNOTATION.tsv',sep='\t',index=False)
canonical=set(g for g,row in oldpanel.set_index('carrot_gene_id').iterrows() if g) # Report this fixed old carrot panel denominator separately.
sets['CANONICAL_OLD_CARROT_CANDIDATES']=canonical
coverage=[];go_layers={}
for layer,mp in layers.items():
 go={g:set().union(*(agogo[a] for a in ids)) for g,ids in mp.items()};go_layers[layer]=go;bg={g for g,v in go.items() if v}
 for name,genes in sets.items():
   coverage.append({'layer':layer,'set':name,'genes':len(genes),'mapped':len(genes&set(mp)),'mapped_percent':100*len(genes&set(mp))/len(genes) if genes else 0,'usable_GO':len(genes&bg),'usable_GO_percent':100*len(genes&bg)/len(genes) if genes else 0,'ambiguous_direct':sum(status.get(g,'').endswith('AMBIGUITY') or 'NEAR_TIE' in status.get(g,'') for g in genes),'unmapped':len(genes-set(mp))})
 # Canonical panel denominator is the fixed 23 Arabidopsis regulators, not the
 # three older carrot candidates. This measures candidate coverage, not confirmed orthology.
 mapped_panel=0;go_panel=0;ambiguous_panel=0
 for pr in panel.itertuples(index=False):
  gg={g for g,aa in mp.items() if pr.arabidopsis_id in aa}
  mapped_panel+=bool(gg);go_panel+=bool(gg&bg);ambiguous_panel+=len(gg)>1
 coverage.append({'layer':layer,'set':'CANONICAL_PANEL_23_REGULATORS','genes':23,'mapped':mapped_panel,'mapped_percent':100*mapped_panel/23,'usable_GO':go_panel,'usable_GO_percent':100*go_panel/23,'ambiguous_direct':ambiguous_panel,'unmapped':23-mapped_panel})
pd.DataFrame(coverage).to_csv(O/'R4A_ANNOTATION_COVERAGE_COMPARISON.tsv',sep='\t',index=False)
def classify(layer):
 gains=[]
 for s in ['UNIVERSE','CORE_ALL']:
  v={x['layer']:x['usable_GO_percent'] for x in coverage if x['set']==s};gains.append(v[layer]-v['OLD'])
 return ('MAJOR_IMPROVEMENT' if min(gains)>=15 else 'MODERATE_IMPROVEMENT' if min(gains)>=5 else 'MINOR_IMPROVEMENT' if max(gains)>0 else 'NO_IMPROVEMENT'),gains
summary={'coverage':coverage,'direct_classification':classify('DIRECT_RBH'),'candidate_classification':classify('AUGMENTED_CANDIDATE'),'protein_status_counts':dict(collections.Counter(x['protein_status'] for x in mapping)),'discordant_direct_old':sum(x['direct_old_discordance'] for x in mapping),'canonical_strict_candidate_regulators':sum(bool(x['direct_RBH_candidates']) for x in pan),'canonical_family_only':sum(x['resolution']=='AMBIGUOUS_FAMILY_HOMOLOGY' for x in pan)}
(O/'R4A_ANNOTATION_SUMMARY.json').write_text(json.dumps(summary,indent=2));print(json.dumps({k:v for k,v in summary.items() if k!='coverage'},indent=2),flush=True)
material=summary['candidate_classification'][0] in ['MAJOR_IMPROVEMENT','MODERATE_IMPROVEMENT']
def bh(p):
 p=np.asarray(p);ix=np.argsort(p);v=p[ix]*len(p)/np.arange(1,len(p)+1);q=np.empty(len(p));q[ix]=np.minimum(1,np.minimum.accumulate(v[::-1])[::-1]);return q
if material:
 go=go_layers['AUGMENTED_CANDIDATE'];bg={g for g,v in go.items() if v};termgenes=collections.defaultdict(set)
 for g in bg:
  for t in go[g]:termgenes[t].add(g)
 tested={t:gs for t,gs in termgenes.items() if 5<=len(gs)<=5000};enr=[]
 for name in ['CORE_ALL','CORE_UP','CORE_DOWN','PRACTICAL_INTERACTION','LARGE_INTERACTION']:
  selected=sets[name]&bg;n=len(selected);N=len(bg)
  for t,gs in sorted(tested.items()):
   hit=selected&gs;k=len(hit);K=len(gs);p=hypergeom.sf(k-1,N,K,n) if n else 1
   enr.append({'set':name,'GO':t,'term':terms[t].get('name',''),'namespace':terms[t].get('namespace',''),'gene_count':k,'testable_set_size':n,'original_set_size':len(sets[name]),'background_count':K,'background_size':N,'fold_enrichment':(k/n)/(K/N) if n else np.nan,'raw_p':p,'contributing_genes':';'.join(sorted(hit))})
 f=pd.DataFrame(enr);f['global_FDR']=bh(f.raw_p);f.to_csv(O/'R4A_REANNOTATION_GO_FULL.tsv',sep='\t',index=False)
 oldgo=pd.read_csv(A/'go_enrichment/R1_GO_ENRICHMENT_FULL.tsv',sep='\t');comparison=oldgo.merge(f,on=['set','GO'],how='outer',suffixes=('_old','_new'));comparison.to_csv(O/'R4A_GO_OLD_NEW_ALL.tsv',sep='\t',index=False)
 target=f[f.GO.isin(['GO:0005840','GO:0000151'])];target.to_csv(O/'R4A_FROZEN_FUNCTIONAL_PATTERNS.tsv',sep='\t',index=False)
 summary['GO_tests']=len(f);summary['GO_background']=len(bg);summary['new_significant_counts']=f[f.global_FDR<=.05].groupby('set').size().to_dict();summary['frozen_terms']=target.to_dict('records')
 (O/'R4A_ANNOTATION_SUMMARY.json').write_text(json.dumps(summary,indent=2))
else:
 pd.DataFrame(columns=['set','GO','term','raw_p','global_FDR','status']).to_csv(O/'R4A_REANNOTATION_GO_FULL.tsv',sep='\t',index=False)
 print('NO_ENRICHMENT: coverage gate not material')
print('ANNOTATION_COMPARISON_COMPLETE',flush=True)
