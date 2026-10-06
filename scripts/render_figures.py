from pathlib import Path
import os,sys,json,itertools
R=Path(__file__).resolve().parents[1];P=R/'results';os.environ['MPLCONFIGDIR']=str(R/'work/publication/matplotlib');sys.path.insert(0,str(R/'work/interaction_architecture/python_packages'))
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
(R/'figures/reproduced').mkdir(exist_ok=True,parents=True)
G=['C815','C819','C824','Ws'];C=dict(zip(G,['#0072B2','#009E73','#D55E00','#CC79A7']))
def read(p):return pd.read_csv(p,sep='\t',dtype={'gene_id':str},float_precision='round_trip')
def table(n):return read(next((P/'supplementary_tables').glob(f'S{n}_*')))
def source(n):return read(R/'figures/source_data'/n)
def save(fig,n):
 fig.set_facecolor('white');fig.savefig(R/f'figures/reproduced/Figure{n}.pdf',bbox_inches='tight');fig.savefig(R/f'figures/reproduced/Figure{n}.svg',bbox_inches='tight');fig.savefig(R/f'figures/reproduced/Figure{n}.png',dpi=350,bbox_inches='tight');plt.close(fig)
def title(ax,s):ax.set_title(s,loc='left',fontweight='bold',pad=12)
pca=source('pca_coordinates.tsv');var=source('variance_partition.tsv');dis=source('genotype_dispersion.tsv').query("record=='primary_summary'").iloc[0]
fig,ax=plt.subplots(2,2,figsize=(11,8),layout='constrained');a=ax[0,0];a.axis('off');title(a,'A  Experimental design')
a.text(0,.85,'4 accession labels × 2 collection times × 3 dishes',fontsize=12);a.text(0,.64,'C815    C819    C824    Ws\nDay0 → Day20 callus-induction interval',linespacing=2);a.text(0,.30,'24 paired-end libraries\nThree whole cultured tissues pooled per dish\nDish is the biological replicate (n = 3 per cell)',linespacing=1.8)
a=ax[0,1]
for g in G:
 for d,marker in [(0,'o'),(20,'^')]:
  t=pca[(pca.genotype==g)&(pca.day==d)];a.scatter(t.PC1,t.PC2,c=C[g],marker=marker,label=f'{g} D{d}',s=45)
title(a,'B  Descriptive global structure');a.set(xlabel='PC1 (saved coordinates)',ylabel='PC2 (saved coordinates)');a.legend(ncol=2,fontsize=8,frameon=False)
a=ax[1,0];a.barh(['Genotype','Time','Interaction','Within-cell residual'],100*var.fraction_total,color=['#0072B2','#999999','#D55E00','#cccccc']);a.set(xlim=(0,100),xlabel='Descriptive transformed-expression SS (%)');title(a,'C  Global variance partition')
a=ax[1,1];a.bar(['Day0','Day20'],[dis.dispersion_day0,dis.dispersion_day20],color=['#999999','#0072B2']);a.set(ylim=(0,1.4),ylabel='Mean squared centroid distance per gene');title(a,'D  No established global convergence');a.text(.02,.96,f'Ratio = {dis.ratio_D20_D0:.4f}\nDish-bootstrap 95% CI\n{dis.ratio_ci_low:.4f}–{dis.ratio_ci_high:.4f}',transform=a.transAxes,va='top',fontsize=10);save(fig,1)
s=table(3);vectors=source('response_vectors.tsv.gz');matrix=np.eye(4)
for t in s.itertuples():i,j=G.index(t.genotype_a),G.index(t.genotype_b);matrix[i,j]=matrix[j,i]=t.pearson_r
fig=plt.figure(figsize=(12,11),layout='constrained');gs=fig.add_gridspec(3,3);a=fig.add_subplot(gs[0,0]);im=a.imshow(matrix,cmap='RdBu_r',vmin=-1,vmax=1);a.set(xticks=range(4),yticks=range(4),xticklabels=G,yticklabels=G);title(a,'A  Pearson response similarity')
for i in range(4):
 for j in range(4):a.text(j,i,f'{matrix[i,j]:.2f}',ha='center',va='center',color='white' if abs(matrix[i,j])>.6 else 'black')
fig.colorbar(im,ax=a,shrink=.7,label='Pearson r')
a=fig.add_subplot(gs[0,1:]);y=np.arange(6);a.hlines(y,s.pearson_r_ci_low,s.pearson_r_ci_high,color='black');a.scatter(s.pearson_r,y,c='#0072B2');a.axvline(0,c='#999999',ls='--');a.set(xlim=(-1,1),yticks=y,yticklabels=s.genotype_a+'–'+s.genotype_b,xlabel='Pearson r; 95% within-cell dish-bootstrap interval');title(a,'B  All six pairs; 24,961 genes')
lim=float(np.abs(vectors[G].to_numpy()).max());h=[]
for k,t in enumerate(s.itertuples()):
 a=fig.add_subplot(gs[1+k//3,k%3]);v=a.hexbin(vectors[t.genotype_a],vectors[t.genotype_b],gridsize=50,extent=(-lim,lim,-lim,lim),mincnt=1,cmap='viridis',norm=LogNorm(vmin=1,vmax=24961));h.append(v);a.plot([-lim,lim],[-lim,lim],c='#999999',lw=.6);a.set(xlim=(-lim,lim),ylim=(-lim,lim),xlabel=f'{t.genotype_a} response',ylabel=f'{t.genotype_b} response');title(a,f'{chr(67+k)}  r = {t.pearson_r:.3f}; ρ = {t.spearman_r:.3f}')
fig.colorbar(h[-1],ax=fig.axes[2:],shrink=.6,label='Genes per hexagon (log scale)');save(fig,2)
ef=table(5);f=table(4);fig,ax=plt.subplots(1,3,figsize=(12,4),layout='constrained');a=ax[0];a.bar(['Interaction','Other / ineligible'],[15053,24961-15053],color=['#D55E00','#bbbbbb']);a.set(ylabel='Genes',ylim=(0,24961));a.tick_params(axis='x',rotation=15);title(a,'A  15,053 / 24,961 (60.31%)')
a=ax[1];a.hist(ef.response_range_log2,bins=40,color='#0072B2');a.axvline(.5,c='black',ls='--');a.axvline(1,c='black',ls=':');a.set(xlabel='Genotype response range (log2FC)',ylabel='Significant interaction genes');title(a,'B  Full effect-size distribution')
a=ax[2];counts=[567,3877,10609];a.bar(['Small\n<0.5','Moderate\n0.5–<1','Large\n≥1'],counts,color=['#bbbbbb','#E69F00','#D55E00']);a.set(ylim=(0,12500),ylabel='Significant interaction genes');title(a,'C  Predominantly nontrivial effects')
for i,v in enumerate(counts):a.text(i,v+220,f'{v:,}\n{v/15053:.2%}',ha='center',fontsize=9)
save(fig,3)
loo=table(6);rep=table(9);fig,ax=plt.subplots(1,3,figsize=(13,4.6),layout='constrained');a=ax[0];a.bar(loo.excluded_genotype,100*loo.practical_fraction_universe,color=[C[g] for g in loo.excluded_genotype]);a.set(ylim=(0,100),xlabel='Excluded genotype (n = 18 per fit)',ylabel='Practical interaction genes / 24,961 (%)');title(a,'A  Signal persists without Ws')
for i,v in enumerate(loo.practical_heterogeneous_genes):a.text(i,v/24961*100+2,f'{v:,}',ha='center',fontsize=9)
a=ax[1];a.scatter(s.pearson_r,range(6),color=['#CC79A7' if 'Ws' in [t.genotype_a,t.genotype_b] else '#0072B2' for t in s.itertuples()]);a.set(xlim=(-1,1),yticks=range(6),yticklabels=s.genotype_a+'–'+s.genotype_b,xlabel='Pearson response r');a.axvline(0,c='#999999',ls='--');a.set_title('B  Ws outlier topology',loc='left',fontweight='bold',pad=30);a.text(.02,1.02,'Bootstrap support 98.6%',transform=a.transAxes,fontsize=9)
a=ax[2];a.bar(range(24),100*rep.practical_fraction_universe,color=[C[g] for g in rep.genotype]);a.set(ylim=(0,100),xticks=range(24),xticklabels=rep.omitted_sample,ylabel='Practical interaction genes / 24,961 (%)',xlabel='Omitted dish (n = 23 per fit)');a.tick_params(axis='x',rotation=90,labelsize=7);title(a,'C  24 / 24 omissions stable');save(fig,4)
core=table(7);en=table(8);fig,ax=plt.subplots(2,2,figsize=(12,9),layout='constrained');a=ax[0,0];a.bar(['Fixed core','Secondary stability'],[968,800],color=['#0072B2','#009E73']);a.set(ylim=(0,1100),ylabel='Genes');title(a,'A  Shared subset within 24,961 genes');a.text(0,990,'968 (3.88%)',ha='center');a.text(1,820,'800 / 968 (82.64%)',ha='center')
a=ax[0,1]
for k,g in enumerate(G):
 vals=core[g+'_log2FC'];v=a.violinplot([vals[vals>0],vals[vals<0]],positions=[k-.16,k+.16],widths=.28,showmedians=True,showextrema=False)
 for b in v['bodies']:b.set_facecolor(C[g]);b.set_alpha(.65)
 v['cmedians'].set_color('black')
a.axhline(0,c='#999999',lw=.7);a.set(xticks=range(4),xticklabels=G,ylabel='Original NB response (log2FC)');title(a,'B  Direction and effect consistency; n = 968')
a=ax[1,0];a.scatter(core.bootstrap_all_four_original_direction_probability,core.dish_omission_original_criteria_retention_fraction,c=np.where(core.robust,'#009E73','#999999'),s=12,alpha=.5);a.axvline(.9,c='black',ls='--');a.axhline(.8,c='black',ls='--');a.set(xlim=(0,1.02),ylim=(0,1.02),xlabel='Original-direction bootstrap probability',ylabel='Original-criteria dish-omission retention');title(a,'C  Internal stability; fixed membership')
a=ax[1,1];cats=list(en.category.unique())
for qi,(query,col,offset) in enumerate([('core_up','#0072B2',-.12),('core_down','#D55E00',.12)]):
 t=en[en['query']==query]
 for j,row in enumerate(t.itertuples()):
  if row.zero_annotated_background:a.text(0,j+offset,'N/A',fontsize=7,va='center');continue
  a.scatter(row.fold_enrichment,j+offset,s=15+row.gene_count*2,c=col,edgecolors='black' if row.FDR_36_predefined_tests<=.05 else 'none',linewidth=1.4)
a.axvline(1,c='#999999',ls='--');a.set(yticks=range(12),yticklabels=[x.replace('_',' ') for x in cats],xlabel='Annotated fold enrichment (up blue; down orange)');a.tick_params(axis='y',labelsize=8);a.set_title('D  Twelve predefined categories',loc='left',fontweight='bold',pad=30);a.text(.02,1.015,'Black edge: joint 36-test FDR ≤ 0.05; size: overlap',transform=a.transAxes,fontsize=8);save(fig,5)
