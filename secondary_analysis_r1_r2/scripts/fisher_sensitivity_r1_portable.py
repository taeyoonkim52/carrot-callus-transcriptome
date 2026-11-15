from pathlib import Path
import os,sys,json
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
R=Path(__file__).resolve().parents[2];O=R/'work/r1_reproduction'
import numpy as np,pandas as pd
d=np.load(R/'work/dedifferentiation/primary_model.npz');p=pd.read_csv(R/'results/supplementary_tables/S4_interaction_statistics.tsv',sep='\t',dtype={'gene_id':str})
m=pd.read_csv(R/'metadata/sample_manifest_24.tsv',sep='\t',dtype={'day':str}).set_index('sample').loc[d['sample_ids']]
X=d['design'];mu=d['size_factors'][:,None]*np.exp(X@d['log_coefficients'].T);w=mu/(1+mu*d['dispersions'][None,:]);G=['C815','C819','C824','Ws'];cells=[];variances=[]
for g in G:
 infos=[]
 for t in ['0','20']:
  mask=((m.genotype==g)&(m.day==t)).to_numpy();infos.append(w[mask].sum(axis=0));cells.append(infos[-1])
 variances.append((1/infos[0]+1/infos[1])/np.log(2)**2)
v=np.array(variances).T;idx=np.flatnonzero(p.interaction_FDR.to_numpy()<=.05);b=p[[g+'_log2FC' for g in G]].to_numpy()
original=pd.read_csv(O/'R1_LFC_SENSITIVITY.tsv',sep='\t');original.to_csv(O/'R1_LFC_INITIAL_HMH_DIAGNOSTIC.tsv',sep='\t',index=False)
rows=[original.iloc[0].to_dict()];details=[]
for tau in [1.,2.]:
 e=b[idx]*(tau*tau/(v[idx]+tau*tau));rg=np.ptp(e,axis=1);sm=int((rg<.5).sum());mo=int(((rg>=.5)&(rg<1)).sum());la=int((rg>=1).sum())
 row=dict(analysis='EXPECTED_FISHER_SHRINK_TAU_'+str(tau),retained_interaction_genes=len(rg),SMALL=sm,MODERATE=mo,LARGE=la,LARGE_percent=100*la/len(rg),practical_genes=mo+la,practical_percent_original_universe=100*(mo+la)/24961,median_response_range=np.median(rg),Q25=np.quantile(rg,.25),Q75=np.quantile(rg,.75),P90=np.quantile(rg,.9),P95=np.quantile(rg,.95),P99=np.quantile(rg,.99),maximum=rg.max(),above5=int((rg>5).sum()),above10=int((rg>10).sum()),above15=int((rg>15).sum()),above20=int((rg>20).sum()));rows.append(row)
 for j,eff,ran in zip(idx,e,rg):details.append(dict(gene_id=d['genes'][j],analysis=row['analysis'],response_range=ran,**dict(zip([g+'_log2FC' for g in G],eff))))
rows+=original[original.analysis.str.startswith('CPM_')].to_dict('records');out=pd.DataFrame(rows);out.to_csv(O/'R1_LFC_SENSITIVITY.tsv',sep='\t',index=False);pd.DataFrame(details).to_csv(O/'R1_LFC_FISHER_GENE_DETAILS.tsv',sep='\t',index=False)
# Batched condition numbers use only saved fitted quantities, not refitting.
M=np.einsum('si,sj,sg->gij',X,X,w);cond=np.linalg.cond(M);mininfo=np.min(cells,axis=0);warning=(cond>1e12)|(mininfo<1e-6)
warn=p.loc[warning].copy();warn['information_condition_number']=cond[warning];warn['minimum_cell_information']=mininfo[warning];warn['significant_original_interaction']=warn.interaction_FDR<=.05;warn.to_csv(O/'R1_LFC_NUMERICAL_WARNINGS.tsv',sep='\t',index=False)
decision='ROBUST' if ((out.iloc[1:].practical_percent_original_universe>=10)&(out.iloc[1:].LARGE_percent>=50)).all() else 'PARTIALLY_ROBUST' if (out.iloc[1:].practical_percent_original_universe>=10).all() else 'NOT_ROBUST'
summary={'classification':decision,'extreme_above10':int(original.iloc[0].above10),'original_max':float(original.iloc[0].maximum),'numerical_warning_genes':int(warning.sum()),'numerical_warning_significant':int((warning&(p.interaction_FDR.to_numpy()<=.05)).sum()),'numerical_warning_core':int((warning&p.conserved_core.to_numpy()).sum()),'condition_threshold':1e12,'minimum_cell_information_threshold':1e-6}
(O/'R1_SENSITIVITY_GATE.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2));print(out.to_string(index=False))
