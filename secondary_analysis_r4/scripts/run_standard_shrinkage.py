from pathlib import Path
import os,json,time,hashlib,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
RETAINED=ROOT/'secondary_analysis_r4'
O=ROOT/'work/r4_reproduction'
O.mkdir(parents=True,exist_ok=True)
R=ROOT
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
os.environ['PYTHONDONTWRITEBYTECODE']='1'
import numpy as np,pandas as pd,pydeseq2
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference
def main():
 d=np.load(R/'work/dedifferentiation/primary_model.npz');t=np.load(R/'work/dedifferentiation/technical_expression_unicode.npz')
 assert np.array_equal(d['imported_integer_counts'],np.rint(t['scaled_counts'].T).astype(np.int64))
 assert np.array_equal(d['genes'],t['genes']) and np.array_equal(d['sample_ids'],t['sample_ids'])
 original=pd.read_csv(ROOT/'results/supplementary_tables/S4_interaction_statistics.tsv',sep='\t',dtype={'gene_id':str});assert original.gene_id.tolist()==d['genes'].tolist()
 ix=np.flatnonzero(original.interaction_FDR.to_numpy()<=.05);assert len(ix)==15053
 m=pd.read_csv(ROOT/'metadata/sample_manifest_24.tsv',sep='\t',dtype={'day':str}).set_index('sample').loc[d['sample_ids']]
 eta=d['design']@d['log_coefficients'].T;Y=d['imported_integer_counts'];sf=d['size_factors'];results=[];records=[]
 for g in ['C815','C819','C824','Ws']:
  dest=O/f'shrinkage_{g}.tsv'
  if dest.exists():results.append(pd.read_csv(dest,sep='\t',dtype={'gene_id':str}).set_index('gene_id'));continue
  start=time.time();samples=np.flatnonzero((m.genotype==g).to_numpy());mm=m.iloc[samples]
  design=pd.DataFrame({'Intercept':np.ones(6),'time_D20':(mm.day=='20').astype(float)},index=mm.index)
  baseline=eta[samples[np.flatnonzero(mm.day.to_numpy()=='0')[0]],ix]
  effect=original.loc[ix,g+'_log2FC'].to_numpy()*np.log(2)
  frozen=np.column_stack([baseline,effect]); predicted=design.to_numpy()@frozen.T
  assert np.allclose(predicted,eta[np.ix_(samples,ix)],atol=1e-10,rtol=1e-10)
  assert np.linalg.matrix_rank(design)==2
  ds=DeseqDataSet(counts=pd.DataFrame(Y[np.ix_(samples,ix)],index=mm.index,columns=d['genes'][ix]),metadata=mm,design=design,refit_cooks=False,quiet=True,n_cpus=2)
  ds.obs['size_factors']=sf[samples];ds.var['dispersions']=d['dispersions'][ix];ds.var['_normed_means']=(ds.X/sf[samples,None]).mean(axis=0)
  ds.varm['LFC']=pd.DataFrame(frozen,index=ds.var_names,columns=design.columns)
  nz=np.flatnonzero((ds.X!=0).any(axis=0)&np.isfinite(ds.var.dispersions.to_numpy())&(ds.var.dispersions.to_numpy()>0))
  ds.non_zero_idx=nz;ds.non_zero_genes=ds.var_names[nz]
  stats=DeseqStats(ds,contrast=np.array([0.,1.]),independent_filter=False,cooks_filter=False,n_cpus=2,quiet=True)
  stats.SE=pd.Series(np.nan,index=ds.var_names)
  # Supported public API. No summary(), deseq2(), dispersion or FDR estimation.
  stats.lfc_shrink(coeff='time_D20',adapt=False)
  conv=stats._LFC_shrink_converged.astype('boolean')
  out=pd.DataFrame({'gene_id':ds.var_names,'original_log2FC':effect/np.log(2),'shrunk_log2FC':stats.LFC.time_D20.to_numpy()/np.log(2),'MAP_SE_log2':stats.SE.to_numpy()/np.log(2),'converged':conv.to_numpy(),'nonzero_genotype':np.isin(np.arange(len(ix)),nz)})
  out.to_csv(dest,sep='\t',index=False);results.append(out.set_index('gene_id'))
  rec={'genotype':g,'seconds':time.time()-start,'genes':len(ix),'nonzero_genes':len(nz),'converged':int(conv.fillna(False).sum()),'design_prediction_max_absolute_difference':float(np.max(np.abs(predicted-eta[np.ix_(samples,ix)])))};records.append(rec)
  (O/'STANDARD_SHRINKAGE_PROGRESS.json').write_text(json.dumps(records,indent=2));print(rec,flush=True)
 out=pd.DataFrame({'gene_id':d['genes'][ix]});valid=np.ones(len(ix),dtype=bool)
 for g,res in zip(['C815','C819','C824','Ws'],results):
  res=res.loc[out.gene_id];out[g+'_shrunk_log2FC']=res.shrunk_log2FC.to_numpy();valid&=res.converged.fillna(False).to_numpy().astype(bool)&res.nonzero_genotype.to_numpy()&np.isfinite(res.shrunk_log2FC.to_numpy())&np.isfinite(res.MAP_SE_log2.to_numpy())
 out['all_four_valid']=valid;arr=out[[g+'_shrunk_log2FC' for g in ['C815','C819','C824','Ws']]].to_numpy();ranges=np.ptp(arr,axis=1);out['response_range']=np.where(valid,ranges,np.nan);out.to_csv(O/'R4A_STANDARD_SHRINKAGE_GENE_DETAILS.tsv',sep='\t',index=False)
 rg=ranges[valid];un=original.loc[ix,'response_effect_range_log2'].to_numpy()[valid]
 summary={'method':'PyDESeq2 0.5.4 DeseqStats.lfc_shrink; adapt=False; Cauchy scale 1 natural-log units; nuisance scale 15','classification':'VALID_STANDARD_SHRINKAGE_AVAILABLE' if valid.all() else 'PARTIALLY_VALID','original_fixed_interaction_genes':len(ix),'retained_all_four_converged':int(valid.sum()),'invalid_or_nonconverged':int((~valid).sum()),'SMALL':int((rg<.5).sum()),'MODERATE':int(((rg>=.5)&(rg<1)).sum()),'LARGE':int((rg>=1).sum()),'LARGE_percent':float(100*np.mean(rg>=1)),'range_quantiles':dict(zip(['min','Q25','median','Q75','P90','P95','P99','max'],np.quantile(rg,[0,.25,.5,.75,.9,.95,.99,1]).tolist())),'original_vs_shrunk_range_Pearson':float(np.corrcoef(un,rg)[0,1]),'genotype_effect_correlations':{g:float(np.corrcoef(original.loc[ix,g+'_log2FC'].to_numpy()[valid],arr[valid,k])[0,1]) for k,g in enumerate(['C815','C819','C824','Ws'])},'primary_p_FDR_changed':False,'core_membership_changed':False,'conditional_likelihood_equivalence':'saturated 8-cell NB likelihood factors into four independent two-endpoint genotype blocks with fixed primary size factors/dispersions','limitation':'Four marginal MAP effects yield a descriptive range, not a joint posterior for range; fixed prior is not adaptive apeglm in R.'}
 (O/'R4A_STANDARD_SHRINKAGE_SUMMARY.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2),flush=True)
if __name__=='__main__':main()
