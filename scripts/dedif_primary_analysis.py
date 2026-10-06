"""Predefined primary analyses. Requires a completed PASS technical gate."""
import sys,os,json,itertools,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];W=ROOT/'work/dedifferentiation';O=ROOT/'outputs/dedifferentiation_gate'
sys.path.insert(0,str(W/'python_packages'));os.environ['PYTHONPATH']=str(W/'python_packages');os.environ['MPLCONFIGDIR']=str(W/'mpl_config')
for variable in ['OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS']:os.environ[variable]='1'
import numpy as np,pandas as pd
from scipy.stats import chi2
from scipy.spatial.distance import pdist,squareform
from sklearn.decomposition import PCA
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
gate=json.loads((O/'TECHNICAL_GATE.json').read_text())
if gate['technical_gate']!='PASS':raise SystemExit('Biological analysis blocked: technical gate not PASS')
for file in ['interaction_summary.tsv','conserved_core.tsv','model_decision.tsv']:
    if (O/file).exists():raise SystemExit('Preserved primary output exists; refusing to overwrite '+file)
data=np.load(W/'technical_expression_unicode.npz',allow_pickle=False);Y=data['Y'];genes=data['genes'];scaled=data['scaled_counts'];ids=data['sample_ids']
meta=pd.read_csv(O/'sample_manifest_24.tsv',sep='\t',dtype={'day':str}).set_index('sample').loc[ids];meta['time']='D'+meta.day
G=['C815','C819','C824','Ws'];pairs=list(itertools.combinations(range(4),2));cells=[[np.flatnonzero(((meta.genotype==g)&(meta.day==t)).to_numpy()) for t in ['0','20']] for g in G]
assert all(len(i)==3 for cell in cells for i in cell)
prior_files=['variance_partition.tsv','pca_coordinates.tsv','response_vector_similarity.tsv','response_vectors.tsv.gz','genotype_dispersion.tsv']
if all((O/f).exists() for f in prior_files):
    # Resume from preserved descriptive results; no PCA/bootstrap recomputation.
    vectors=pd.read_csv(O/'response_vectors.tsv.gz',sep='\t',dtype={'gene_id':str}).set_index('gene_id')
    assert vectors.index.tolist()==genes.tolist()
    R=vectors[G].to_numpy().T
    previous=pd.read_csv(O/'response_vector_similarity.tsv',sep='\t')
    assert list(zip(previous.genotype_a,previous.genotype_b))==[(G[a],G[b]) for a,b in pairs]
    sim=previous[['pearson_r','cosine_similarity','shared_direction_fraction']].to_numpy()
    response_ci=previous[['pearson_r_ci_low','pearson_r_ci_high']].to_numpy()
    ds=pd.read_csv(O/'genotype_dispersion.tsv',sep='\t');dr=ds.loc[ds.record=='primary_summary'].iloc[0]
    disp=np.array([dr.dispersion_day0,dr.dispersion_day20]);ratio=dr.ratio_D20_D0;low=dr.ratio_ci_low;high=dr.ratio_ci_high
    pcs=pd.read_csv(O/'pca_coordinates.tsv',sep='\t').set_index('sample').loc[ids]
    coordinates=pcs[[f'PC{i+1}' for i in range(8)]].to_numpy()
    from types import SimpleNamespace
    total=float(pd.read_csv(O/'variance_partition.tsv',sep='\t').sum_squares.sum())
    pca=SimpleNamespace(explained_variance_ratio_=np.sum(coordinates**2,axis=0)/total)
else:
    if any((O/f).exists() for f in prior_files):raise SystemExit('Incomplete descriptive checkpoint; preserve for review')
    grand=Y.mean(axis=0);means=np.array([[Y[idx].mean(axis=0) for idx in cell] for cell in cells]);gmean=means.mean(axis=1);tmean=means.mean(axis=0)
    ss_g=6*np.sum((gmean-grand)**2);ss_t=12*np.sum((tmean-grand)**2);ss_i=3*np.sum((means-gmean[:,None,:]-tmean[None,:,:]+grand)**2)
    residual=0
    for g in range(4):
        for t in range(2):residual+=np.sum((Y[cells[g][t]]-means[g,t])**2)
    total=np.sum((Y-grand)**2)
    assert np.isclose(ss_g+ss_t+ss_i+residual,total)
    pd.DataFrame({'component':['genotype','time','genotype_x_time','within_cell_residual'],'sum_squares':[ss_g,ss_t,ss_i,residual],'fraction_total':[x/total for x in [ss_g,ss_t,ss_i,residual]],'degrees_freedom':[3,1,3,16]}).to_csv(O/'variance_partition.tsv',sep='\t',index=False)
    pca=PCA(n_components=8,svd_solver='full');coordinates=pca.fit_transform(Y)
    pd.DataFrame(coordinates,index=ids,columns=[f'PC{i+1}' for i in range(8)]).assign(genotype=meta.genotype,day=meta.day).to_csv(O/'pca_coordinates.tsv',sep='\t',index_label='sample')
    R=means[:,1]-means[:,0]
    def similarity(r):
        return np.array([[np.corrcoef(r[a],r[b])[0,1],np.dot(r[a],r[b])/(np.linalg.norm(r[a])*np.linalg.norm(r[b])),np.mean(np.sign(r[a])==np.sign(r[b]))] for a,b in pairs])
    def dispersion(m):return np.array([np.mean([np.mean((m[a,t]-m[b,t])**2) for a,b in pairs]) for t in range(2)])
    sim=similarity(R);disp=dispersion(means);rng=np.random.default_rng(20261004);boot_sim=[];boot_disp=[]
    for b in range(2000):
        m=np.array([[Y[rng.choice(idx,3,replace=True)].mean(axis=0) for idx in cell] for cell in cells]);boot_sim.append(similarity(m[:,1]-m[:,0]));boot_disp.append(dispersion(m))
    boot_sim=np.array(boot_sim);boot_disp=np.array(boot_disp)
    rows=[]
    for k,(a,b) in enumerate(pairs):
        row={'genotype_a':G[a],'genotype_b':G[b]}
        for j,key in enumerate(['pearson_r','cosine_similarity','shared_direction_fraction']):row[key]=sim[k,j];row[key+'_ci_low'],row[key+'_ci_high']=np.quantile(boot_sim[:,k,j],[.025,.975])
        row['rms_response_a']=np.sqrt(np.mean(R[a]**2));row['rms_response_b']=np.sqrt(np.mean(R[b]**2));rows.append(row)
    pd.DataFrame(rows).to_csv(O/'response_vector_similarity.tsv',sep='\t',index=False)
    pd.DataFrame(R.T,index=genes,columns=G).to_csv(O/'response_vectors.tsv.gz',sep='\t',compression={'method':'gzip','mtime':0},index_label='gene_id')
    ratios=boot_disp[:,1]/boot_disp[:,0];ratio=disp[1]/disp[0];low,high=np.quantile(ratios,[.025,.975]);drows=[]
    for a,b in pairs:
        for t in range(2):
            squared=np.mean((means[a,t]-means[b,t])**2)
            noise=(np.var(Y[cells[a][t]],axis=0,ddof=1)/3+np.var(Y[cells[b][t]],axis=0,ddof=1)/3).mean()
            drows.append({'record':'centroid_pair','genotype_a':G[a],'genotype_b':G[b],'day':[0,20][t],'squared_rms_distance':squared,'noise_corrected_squared_rms':squared-noise})
    drows.append({'record':'primary_summary','dispersion_day0':disp[0],'dispersion_day20':disp[1],'ratio_D20_D0':ratio,'ratio_ci_low':low,'ratio_ci_high':high,'bootstrap_replicates':2000,'sampling_unit':'independent dish within genotype/time; unpaired'})
    pd.DataFrame(drows).to_csv(O/'genotype_dispersion.tsv',sep='\t',index=False)
    response_ci=np.quantile(boot_sim[:,:,0],[.025,.975],axis=0).T
countframe=pd.DataFrame(np.rint(scaled.T).astype(np.int64),index=ids,columns=genes)
dds=DeseqDataSet(counts=countframe,metadata=meta,design='~genotype * time',n_cpus=4,quiet=False)
dds.deseq2();X=np.asarray(dds.obsm['design_matrix']);names=list(dds.obsm['design_matrix'].columns)
assert X.shape==(24,8) and np.linalg.matrix_rank(X)==8
timeidx=names.index('time[T.D20]');interactionidx=[i for i,n in enumerate(names) if ':' in n];assert len(interactionidx)==3
contrast=np.zeros(8);contrast[timeidx]=1;contrast[interactionidx]=.25
stats=DeseqStats(dds,contrast=contrast,independent_filter=False,cooks_filter=True,n_cpus=4);stats.summary();common=stats.results_df.copy()
lfcs=[];gstats=[]
for g in G:
    c=np.zeros(8);c[timeidx]=1
    if g!='C815':c[names.index(f'genotype[T.{g}]:time[T.D20]')]=1
    s=DeseqStats(dds,contrast=c,independent_filter=False,cooks_filter=True,n_cpus=4);s.summary();lfcs.append(s.results_df.log2FoldChange.to_numpy());gstats.append(s.results_df.pvalue.to_numpy())
lfcs=np.array(lfcs).T;betas=np.asarray(dds.varm['LFC']);dispersions=dds.var.dispersions.to_numpy();factors=dds.obs.size_factors.to_numpy();pvalues=[]
for j,beta in enumerate(betas):
    mu=factors*np.exp(X@beta);weights=mu/(1+mu*dispersions[j]);M=(X.T*weights)@X;H=np.linalg.inv(M+np.eye(8)*1e-6);cov=(H@M@H)[np.ix_(interactionidx,interactionidx)];v=beta[interactionidx]
    pvalues.append(chi2.sf(v@np.linalg.solve(cov,v),3))
pvalues=np.array(pvalues);valid=np.isfinite(common.pvalue.to_numpy())
if '_LFC_converged' in dds.var:valid&=dds.var['_LFC_converged'].to_numpy().astype(bool)
pvalues[~valid]=np.nan
def bh(p):
    result=np.full(len(p),np.nan);idx=np.flatnonzero(np.isfinite(p));order=idx[np.argsort(p[idx])];q=p[order]*len(order)/np.arange(1,len(order)+1);result[order]=np.minimum(1,np.minimum.accumulate(q[::-1])[::-1]);return result
q=bh(pvalues);timeq=bh(np.where(valid,common.pvalue.to_numpy(),np.nan));effectrange=lfcs.max(axis=1)-lfcs.min(axis=1)
same=(lfcs>0).all(axis=1)|(lfcs<0).all(axis=1);core=same&(np.abs(lfcs).min(axis=1)>=.5)&(effectrange<=1)&(timeq<=.05)
result=pd.DataFrame({'gene_id':genes,'balanced_time_log2FC':common.log2FoldChange.to_numpy(),'balanced_time_p':common.pvalue.to_numpy(),'balanced_time_FDR':timeq,'interaction_3df_p':pvalues,'interaction_FDR':q,'response_effect_range_log2':effectrange,'conserved_core':core,'valid_inference':valid})
for k,g in enumerate(G):result[g+'_log2FC']=lfcs[:,k];result[g+'_p']=gstats[k]
result.to_csv(O/'interaction_summary.tsv',sep='\t',index=False);result.loc[core].to_csv(O/'conserved_core.tsv',sep='\t',index=False)
np.savez_compressed(W/'primary_model.npz',design=X,design_columns=np.array(names),log_coefficients=betas,dispersions=dispersions,size_factors=factors,genes=genes,sample_ids=ids,imported_integer_counts=countframe.to_numpy(),valid_inference=valid)
summary={'testable_genes':len(genes),'valid_inference_genes':int(valid.sum()),'common_time_FDR_genes':int((timeq<=.05).sum()),'interaction_FDR_genes':int((q<=.05).sum()),'conserved_core_genes':int(core.sum()),'core_and_interaction_overlap_genes':int((core&(q<=.05)).sum()),'dispersion_trend_fit':dds.uns['disp_function_type'],'median_response_correlation':float(np.median(sim[:,0])),'dispersion_ratio':float(ratio),'dispersion_ratio_ci':[float(low),float(high)],'gene_count_approximation':'lengthScaledTPM, nearest-integer imported estimate counts','scope':'Day0/Day20; no late or phenotype inputs'}
strong_a=(response_ci[:,0]>0).all() and np.median(sim[:,0])>=.5 and core.mean()>=.1
partial_a=np.median(sim[:,0])>0 and same.mean()>.5
strong_b=ratio<=.8 and high<1
strong_c=(q<=.05).mean()>=.1 and np.mean(effectrange[q<=.05]>=1)>=.5
decisions={'A':'SUPPORTED' if strong_a else 'PARTIALLY_SUPPORTED' if partial_a else 'NOT_SUPPORTED','B':'SUPPORTED' if strong_b else 'PARTIALLY_SUPPORTED' if ratio<1 else 'NOT_SUPPORTED','C':'SUPPORTED' if strong_c else 'PARTIALLY_SUPPORTED' if (q<=.05).any() else 'NOT_SUPPORTED'}
summary['model_decisions']=decisions
pd.DataFrame([{'model':k,'decision':v,'criteria':'Frozen plan; primary quantitative metrics, not PCA'} for k,v in decisions.items()]).to_csv(O/'model_decision.tsv',sep='\t',index=False)
(O/'PRIMARY_RESULTS_FROZEN.json').write_text(json.dumps(summary,indent=2))
