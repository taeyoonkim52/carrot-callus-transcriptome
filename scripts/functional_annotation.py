import os,sys,json,gzip,collections,hashlib
from pathlib import Path
R=Path(__file__).resolve().parents[1];W=R/'work/interaction_architecture';P=R/'outputs/dedifferentiation_gate';D=R/'work/dedifferentiation';O=R/'outputs/interaction_architecture_gate'
import numpy as np,pandas as pd
from scipy.stats import hypergeom
loo_saved=pd.read_csv(O/'leave_one_genotype_out.tsv',sep='\t')
for i,row in loo_saved.iterrows():
    with np.load(W/('exclude_'+row.excluded_genotype+'.npz'),allow_pickle=False) as fit:
        ranges=np.ptp(fit['effects'],axis=1);sig=fit['interaction_q']<=.05
        loo_saved.loc[i,'valid_inference_genes']=int(fit['valid'].sum())
        for name,hit in [('SMALL',ranges<.5),('MODERATE',(ranges>=.5)&(ranges<1)),('LARGE',ranges>=1)]:
            loo_saved.loc[i,name+'_significant_genes']=int((sig&hit).sum())
            loo_saved.loc[i,name+'_percent_significant']=100*(sig&hit).sum()/sig.sum() if sig.any() else np.nan
        for name,p in [('q10',.1),('q25',.25),('median',.5),('q75',.75),('q90',.9)]:loo_saved.loc[i,'significant_response_range_'+name]=float(np.quantile(ranges[sig],p))
loo_saved.to_csv(O/'leave_one_genotype_out.tsv',sep='\t',index=False)
(O/'model_sensitivity_manifest.json').write_text(json.dumps(dict(count_model='PyDESeq2 0.5.4 ~genotype*time, original integer lengthScaledTPM counts',fixed_genes=24961,genotype_omissions=4,dish_omissions=24,interaction_tests='joint 2-df genotype omission; 3-df dish omission Wald using same covariance as original',fdr=.05,gene_filter='original universe unchanged; Cook/finite/convergence validity, missing convergence means ineligible',parallel_workers=4,bootstrap_resamples=2000,seed=20261005,primary_reference='GCF_001625215.2-RS_2024_03',upstream_reruns=0),indent=2))
cats={'cell_cycle':'GO:0007049','cell_division':'GO:0051301','chromatin_organization':'GO:0006325','auxin_response':'GO:0009733','cytokinin_response':'GO:0009735','wound_response':'GO:0009611','stress_response':'GO:0006950','cell_wall_organization':'GO:0071555','meristem_maintenance':'GO:0010073','developmental_process':'GO:0032502','photosynthesis':'GO:0015979','primary_metabolic_process':'GO:0044238'}
terms={};cur={};kind=None
for line in (D/'reference/go-edit.obo').read_text(encoding='utf8').splitlines()+['[END]']:
    if line.startswith('['):
        if kind=='[Term]' and 'id' in cur:terms[cur['id']]=cur
        cur={'parents':[],'alt_ids':[]};kind=line
    elif kind=='[Term]':
        if line.startswith('id: '):cur['id']=line[4:]
        elif line.startswith('name: '):cur['name']=line[6:]
        elif line.startswith('is_a: '):cur['parents'].append(line[6:].split()[0])
        elif line.startswith('relationship: part_of '):cur['parents'].append(line.split()[2])
        elif line.startswith('is_obsolete: true'):cur['obsolete']=True
        elif line.startswith('alt_id: '):cur['alt_ids'].append(line.split()[1])
active={i:v for i,v in terms.items() if not v.get('obsolete')};aliases={a:i for i,v in active.items() for a in v['alt_ids']};children=collections.defaultdict(set)
for i,t in active.items():
    for p in t['parents']:children[p].add(i)
def descendants(root):
    seen={root};stack=[root]
    while stack:
        for c in children[stack.pop()]:
            if c not in seen:seen.add(c);stack.append(c)
    return seen
cat_terms={k:descendants(v) for k,v in cats.items()};sets={k:set() for k in cats};old=pd.read_csv(P/'interaction_summary.tsv',sep='\t',dtype={'gene_id':str});universe=set(old.gene_id);annotated=set();unresolved=0
with gzip.open(D/'reference/GCF_001625215.2-RS_2024_03_gene_ontology.gaf.gz','rt',encoding='utf8') as handle:
    for line in handle:
        if line.startswith('!'):continue
        a=line.rstrip().split('\t')
        if len(a)<9 or 'NOT' in a[3].split('|') or a[1] not in universe:continue
        go=aliases.get(a[4],a[4])
        if go not in active:unresolved+=1;continue
        annotated.add(a[1])
        for k,t in cat_terms.items():
            if go in t:sets[k].add(a[1])
core=old[old.conserved_core.astype(str).str.lower()=='true'];inter=old[(old.interaction_FDR<=.05)&(old.response_effect_range_log2>=.5)]
queries={'core_up':set(core.loc[core.balanced_time_log2FC>0,'gene_id']),'core_down':set(core.loc[core.balanced_time_log2FC<0,'gene_id']),'practical_interaction':set(inter.gene_id)};rows=[]
for name,q in queries.items():
    test=q&annotated
    for cat,m in sets.items():
        hit=len(test&m);p=hypergeom.sf(hit-1,len(annotated),len(m),len(test)) if test and m else 1
        rows.append(dict(query=name,category=cat,term=cats[cat],root_name=active[cats[cat]]['name'],selection='PREDEFINED_CATEGORY',direction='Day20 up in all four' if name=='core_up' else 'Day20 down in all four' if name=='core_down' else 'heterogeneous; no pooled direction',query_all_genes=len(q),query_annotated_genes=len(test),expressed_universe=len(universe),annotated_background=len(annotated),category_background=len(m),gene_count=hit,pvalue=p,fold_enrichment=(hit/len(test))/(len(m)/len(annotated)) if test and m else np.nan,zero_annotated_background=len(m)==0))
en=pd.DataFrame(rows);idx=np.argsort(en.pvalue);q=np.empty(len(en));q[idx]=np.minimum(1,np.minimum.accumulate((en.pvalue.to_numpy()[idx]*len(en)/np.arange(1,len(en)+1))[::-1])[::-1]);en['FDR_36_predefined_tests']=q
en[en['query'].str.startswith('core')].to_csv(O/'conserved_core_enrichment.tsv',sep='\t',index=False);en[en['query']=='practical_interaction'].to_csv(O/'interaction_functional_structure.tsv',sep='\t',index=False)
(O/'annotation_provenance.json').write_text(json.dumps(dict(annotated_background=len(annotated),unannotated_testable=len(universe)-len(annotated),unresolved_annotation_records=unresolved,reference_GAF_sha256=hashlib.sha256((D/'reference/GCF_001625215.2-RS_2024_03_gene_ontology.gaf.gz').read_bytes()).hexdigest(),ontology_sha256=hashlib.sha256((D/'reference/go-edit.obo').read_bytes()).hexdigest(),tests=36,background='GO-annotated genes within frozen expressed universe; no annotation treated as missing, not negative',limitations='Electronic/homology annotation; zero-category genes cannot establish absence of biology; broad categories are not mechanistic pathways'),indent=2))
