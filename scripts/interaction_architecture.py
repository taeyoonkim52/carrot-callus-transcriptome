"""Frozen interaction architecture audit, reusing completed upstream results."""
import os,sys,json,time,itertools,hashlib,subprocess,gc
from pathlib import Path
from datetime import datetime,timezone
R=Path(__file__).resolve().parents[1];W=R/'work/interaction_architecture';W.mkdir(exist_ok=True,parents=True)
D=R/'work/dedifferentiation';P=R/'outputs/dedifferentiation_gate';O=R/'outputs/interaction_architecture_gate'
os.environ['PYTHONPATH']=str(W/'python_packages');sys.path.insert(0,str(W/'python_packages'))
os.environ['MPLCONFIGDIR']=str(W/'mpl');os.environ['PYTHONUTF8']='1'
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import numpy as np,pandas as pd
from scipy.stats import chi2,spearmanr,hypergeom
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
start=time.monotonic()
def status(stage,**kwargs):
    (W/'STATUS.json').write_text(json.dumps(dict(stage=stage,pid=os.getpid(),utc=datetime.now(timezone.utc).isoformat(),elapsed_seconds=time.monotonic()-start,**kwargs),indent=2))
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
if (O/'ARCHITECTURE_FROZEN.json').exists():raise SystemExit('Completed architecture exists; preserved')
old=pd.read_csv(P/'interaction_summary.tsv',sep='\t',dtype={'gene_id':str});genes=old.gene_id.to_numpy();N=len(genes)
original_sig=(old.interaction_FDR<=.05).to_numpy();core=old.conserved_core.astype(str).str.lower().eq('true').to_numpy()
assert N==24961 and original_sig.sum()==15053 and core.sum()==968
G=['C815','C819','C824','Ws'];pairs=list(itertools.combinations(range(4),2));effects=old[[g+'_log2FC' for g in G]].to_numpy()
model=np.load(D/'primary_model.npz',allow_pickle=False);assert model['genes'].tolist()==genes.tolist()
counts=model['imported_integer_counts'];ids=model['sample_ids'];meta=pd.read_csv(P/'sample_manifest_24.tsv',sep='\t',dtype={'day':str}).set_index('sample').loc[ids];meta['time']='D'+meta.day
data=np.load(D/'technical_expression_unicode.npz',allow_pickle=False);Y=data['Y'];assert data['genes'].tolist()==genes.tolist()
Rv=pd.read_csv(P/'response_vectors.tsv.gz',sep='\t',dtype={'gene_id':str}).set_index('gene_id').loc[genes,G].to_numpy().T
def bh(p):
    q=np.full(len(p),np.nan);i=np.flatnonzero(np.isfinite(p));i=i[np.argsort(p[i])];q[i]=np.minimum(1,np.minimum.accumulate((p[i]*len(i)/np.arange(1,len(i)+1))[::-1])[::-1]);return q
def magnitudes(e):return np.ptp(e,axis=1),np.median(np.array([np.abs(e[:,a]-e[:,b]) for a,b in itertools.combinations(range(e.shape[1]),2)]),axis=0)
def correlations(rv):return np.array([np.corrcoef(rv[a],rv[b])[0,1] for a,b in pairs])
def topology(rs):
    mat=np.eye(4)
    for k,(a,b) in enumerate(pairs):mat[a,b]=mat[b,a]=rs[k]
    singles=[]
    for a in range(4):
        rem=[b for b in range(4) if b!=a];within=np.median([mat[b,c] for b,c in itertools.combinations(rem,2)]);cross=np.median(mat[a,rem])
        if within>=.5 and cross<=.2 and within-cross>=.3:singles.append(G[a])
    if len(singles)==1:return 'one_outlier:'+singles[0]
    groups=[]
    for b in range(1,4):
        first=[0,b];second=[i for i in range(4) if i not in first];within=np.median([mat[first[0],first[1]],mat[second[0],second[1]]]);cross=np.median([mat[a,c] for a in first for c in second])
        if within>=.5 and cross<=.2 and within-cross>=.3:groups.append('two_group:'+','.join(G[i] for i in first)+'|'+','.join(G[i] for i in second))
    return groups[0] if len(groups)==1 else 'continuous_or_unstructured'
def fit(keep,label):
    out=W/(label+'.npz')
    if out.exists():
        with np.load(out,allow_pickle=False) as z:return {k:z[k] for k in z.files}
    m=meta.iloc[np.flatnonzero(keep)].copy();gg=sorted(m.genotype.unique());n=len(gg)
    cf=pd.DataFrame(counts[keep],index=m.index,columns=genes)
    dds=DeseqDataSet(counts=cf,metadata=m,design='~genotype * time',n_cpus=4,quiet=True);dds.deseq2()
    X=np.asarray(dds.obsm['design_matrix']);names=list(dds.obsm['design_matrix'].columns);ti=names.index('time[T.D20]');ix=[i for i,k in enumerate(names) if ':' in k]
    assert len(ix)==n-1 and np.linalg.matrix_rank(X)==2*n
    c=np.zeros(len(names));c[ti]=1;c[ix]=1/n
    st=DeseqStats(dds,contrast=c,independent_filter=False,cooks_filter=True,n_cpus=4,quiet=True);st.summary()
    valid=np.isfinite(st.results_df.pvalue.to_numpy())
    if '_LFC_converged' in dds.var:valid&=dds.var['_LFC_converged'].fillna(False).to_numpy().astype(bool)
    beta=np.asarray(dds.varm['LFC']);disp=dds.var.dispersions.to_numpy();sf=dds.obs.size_factors.to_numpy();ps=[]
    for j,b in enumerate(beta):
        mu=sf*np.exp(X@b);weights=mu/(1+mu*disp[j]);M=(X.T*weights)@X;H=np.linalg.inv(M+np.eye(X.shape[1])*1e-6);cov=(H@M@H)[np.ix_(ix,ix)];v=b[ix]
        ps.append(chi2.sf(v@np.linalg.solve(cov,v),n-1))
    ps=np.array(ps);ps[~valid]=np.nan;eff=[]
    for g in gg:
        ec=np.zeros(len(names));ec[ti]=1
        if g!=gg[0]:ec[names.index(f'genotype[T.{g}]:time[T.D20]')]=1
        eff.append(beta@ec/np.log(2))
    eff=np.array(eff).T;result=dict(effects=eff,interaction_q=bh(ps),time_q=bh(np.where(valid,st.results_df.pvalue.to_numpy(),np.nan)),valid=valid,genotypes=np.asarray(gg,dtype=str))
    np.savez_compressed(out,**result);del dds,st;gc.collect();return result
ran,med=magnitudes(effects);bins=np.where(ran<.5,'SMALL',np.where(ran<1,'MODERATE','LARGE'));practical=original_sig&(ran>=.5)
pd.DataFrame(dict(gene_id=genes[original_sig],maximum_response_difference_log2=ran[original_sig],median_pairwise_response_difference_log2=med[original_sig],response_range_log2=ran[original_sig],effect_bin=bins[original_sig])).to_csv(O/'interaction_effect_sizes.tsv',sep='\t',index=False)
summ=[]
for b in ['SMALL','MODERATE','LARGE']:
    hit=original_sig&(bins==b);summ.append(dict(effect_bin=b,genes=int(hit.sum()),percent_interaction_genes=100*hit.sum()/original_sig.sum(),median_range=float(np.median(ran[hit])) if hit.any() else np.nan))
pd.DataFrame(summ).to_csv(O/'interaction_effect_size_summary.tsv',sep='\t',index=False)
pairrows=[];rs=correlations(Rv);basepattern=topology(rs);mat=np.eye(4)
for k,(a,b) in enumerate(pairs):
    r1,r2=Rv[a],Rv[b];strong=(effects[:,a]*effects[:,b]<0)&(np.abs(effects[:,a])>=.5)&(np.abs(effects[:,b])>=.5)
    pairrows.append(dict(genotype_a=G[a],genotype_b=G[b],pearson_r=rs[k],spearman_r=spearmanr(r1,r2).statistic,cosine_similarity=np.dot(r1,r2)/np.linalg.norm(r1)/np.linalg.norm(r2),rms_a=np.sqrt(np.mean(r1*r1)),rms_b=np.sqrt(np.mean(r2*r2)),magnitude_ratio_b_a=np.linalg.norm(r2)/np.linalg.norm(r1),concordant_direction_fraction=np.mean(np.sign(r1)==np.sign(r2)),strongly_discordant_fraction_universe=strong.mean(),strongly_discordant_fraction_interaction=strong[original_sig].mean()))
    mat[a,b]=mat[b,a]=rs[k]
pd.DataFrame(pairrows).to_csv(O/'pairwise_response_architecture.tsv',sep='\t',index=False);pd.DataFrame(mat,index=G,columns=G).to_csv(O/'response_similarity_matrix.tsv',sep='\t',index_label='genotype')
loo=[]
for omitted in G:
    status('LEAVE_ONE_GENOTYPE_OUT',current=omitted,completed=len(loo));keep=(meta.genotype!=omitted).to_numpy();f=fit(keep,'exclude_'+omitted);rr,mm=magnitudes(f['effects']);sig=f['interaction_q']<=.05;pr=sig&(rr>=.5);remaining=[G.index(g) for g in G if g!=omitted];concordance=np.median([np.corrcoef(Rv[a],Rv[b])[0,1] for a,b in itertools.combinations(remaining,2)])
    retained=np.sum(pr&practical)/practical.sum()
    loo.append(dict(excluded_genotype=omitted,samples=18,interaction_df=2,significant_genes=int(sig.sum()),significant_fraction_universe=sig.mean(),practical_heterogeneous_genes=int(pr.sum()),practical_fraction_universe=pr.mean(),original_practical_signal_retention=retained,original_interaction_genes_retained=int((sig&original_sig).sum()),median_range_original_interaction=np.median(rr[original_sig]),median_pairwise_difference_original_interaction=np.median(mm[original_sig]),median_response_correlation=concordance,substantial=pr.mean()>=.1))
    pd.DataFrame(loo).to_csv(O/'leave_one_genotype_out.tsv',sep='\t',index=False)
core_sign=np.sign(old.balanced_time_log2FC.to_numpy()[core]);cells=[[np.flatnonzero(((meta.genotype==g)&(meta.day==t)).to_numpy()) for t in ['0','20']] for g in G]
rep=[];core_ret=[];core_same=[];core_ranges=[];pattern_hits=[]
for i,sample in enumerate(ids):
    status('LEAVE_ONE_DISH_OUT',current=str(sample),completed=i,total=24);keep=np.arange(24)!=i;f=fit(keep,'omit_'+str(sample));rr,mm=magnitudes(f['effects']);pr=(f['interaction_q']<=.05)&(rr>=.5)
    cs=(f['effects'][core]>0).all(axis=1)|(f['effects'][core]<0).all(axis=1);ret=cs&(np.min(np.abs(f['effects'][core]),axis=1)>=.5)&(rr[core]<=1)&(f['time_q'][core]<=.05)
    same=(np.sign(f['effects'][core])==core_sign[:,None]).all(axis=1)
    core_ret.append(ret);core_same.append(same);core_ranges.append(rr[core])
    rv=np.array([Y[[j for j in cells[g][1] if j!=i]].mean(axis=0)-Y[[j for j in cells[g][0] if j!=i]].mean(axis=0) for g in range(4)])
    cr=correlations(rv);pat=topology(cr);pattern_hits.append(pat==basepattern)
    rep.append(dict(omitted_sample=str(sample),genotype=meta.iloc[i].genotype,day=meta.iloc[i].day,significant_interaction_genes=int((f['interaction_q']<=.05).sum()),practical_interaction_genes=int(pr.sum()),practical_fraction_universe=pr.mean(),core_original_criteria_retained=int(ret.sum()),core_direction_retained=int(same.sum()),maximum_pairwise_correlation_change=float(np.max(np.abs(cr-rs))),pattern=pat,substantial_and_response_stable=pr.mean()>=.1 and np.max(np.abs(cr-rs))<=.2,**{f'pearson_{G[a]}_{G[b]}':cr[k] for k,(a,b) in enumerate(pairs)}))
    pd.DataFrame(rep).to_csv(O/'replicate_sensitivity.tsv',sep='\t',index=False)
status('DISH_BOOTSTRAP',resamples=2000);rng=np.random.default_rng(20261005);sign_count=np.zeros(core.sum());boot_patterns=[]
for b in range(2000):
    means=np.array([[Y[rng.choice(idx,3,replace=True)].mean(axis=0) for idx in cell] for cell in cells]);rv=means[:,1]-means[:,0]
    sign_count+=(np.sign(rv[:,core].T)==core_sign[:,None]).all(axis=1);boot_patterns.append(topology(correlations(rv))==basepattern)
core_ret=np.array(core_ret);core_same=np.array(core_same);core_ranges=np.array(core_ranges);boot_sign=sign_count/2000
gene_robust=(boot_sign>=.9)&(core_ret.mean(axis=0)>=.8)
cr=pd.DataFrame(dict(gene_id=genes[core],**{g+'_original_log2FC':effects[core,k] for k,g in enumerate(G)},original_direction_consistent=True,original_response_range=ran[core],original_median_pairwise_difference=med[core],original_interaction_FDR=old.interaction_FDR.to_numpy()[core],bootstrap_all_four_original_direction_probability=boot_sign,dish_omission_original_criteria_retention_fraction=core_ret.mean(axis=0),dish_omission_original_direction_retention_fraction=core_same.mean(axis=0),dish_omission_max_response_range=core_ranges.max(axis=0),robust=gene_robust))
cr.to_csv(O/'conserved_core_robustness.tsv',sep='\t',index=False)

(O/'topology_support.json').write_text(json.dumps(dict(pattern=basepattern,bootstrap_support=np.mean(boot_patterns),dish_omission_support=np.mean(pattern_hits)),indent=2))
