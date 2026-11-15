from pathlib import Path
import os,sys,json,itertools
R=Path(__file__).resolve().parents[2];P=R/'results';os.environ['MPLCONFIGDIR']=str(R/'work/figure5_matplotlib')
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
G=['C815','C819','C824','Ws'];C=dict(zip(G,['#0072B2','#009E73','#D55E00','#CC79A7']))
def read(p):return pd.read_csv(p,sep='\t',dtype={'gene_id':str},float_precision='round_trip')
def table(n):return read(next((P/'supplementary_tables').glob(f'S{n}_*')))
def source(n):return read(P/'figure_data'/n)
def save(fig,n):
 import hashlib,json
 from matplotlib.collections import Collection
 from matplotlib.lines import Line2D
 from matplotlib.text import Text
 out=R/'secondary_analysis_r1_r2/figures';out.mkdir(parents=True,exist_ok=True)
 qa=R.parent/'qa';qa.mkdir(exist_ok=True)
 def snapshot():
  rows=[]
  for ai,a in enumerate(fig.axes):
   rows.append(('axes',ai,list(a.get_position().bounds),list(a.get_xlim()),list(a.get_ylim())))
   for art in a.get_children():
    if isinstance(art,Collection):
     rows.append(('collection',ai,[p.vertices.tolist() for p in art.get_paths()],art.get_offsets().tolist(),art.get_sizes().tolist() if hasattr(art,'get_sizes') else [],art.get_facecolors().tolist(),art.get_edgecolors().tolist(),art.get_linewidths().tolist()))
    elif isinstance(art,Line2D):rows.append(('line',ai,np.asarray(art.get_xdata()).tolist(),np.asarray(art.get_ydata()).tolist(),str(art.get_color()),art.get_linewidth(),str(art.get_linestyle())))
    elif isinstance(art,Text) and art is not a._left_title:rows.append(('text',ai,art.get_text(),list(art.get_position())))
  return rows
 fig.set_facecolor('white');fig.canvas.draw();fig.set_layout_engine(None)
 before=snapshot()
 fig.savefig(qa/'Figure5_original_label_regenerated.png',dpi=350,bbox_inches='tight')
 target=fig.axes[3]._left_title
 assert target.get_text()=='D  Twelve predefined categories'
 target.set_text('D  Twelve operational categories')
 fig.canvas.draw();after=snapshot();assert before==after,'Scientific artists changed'
 for ext in ['pdf','svg','png']:
  kw={'dpi':350} if ext=='png' else {}
  fig.savefig(out/('Figure5_R3_operational.'+ext),bbox_inches='tight',**kw)
 digest=hashlib.sha256(json.dumps(before,sort_keys=True).encode()).hexdigest()
 (out/'Figure5_artist_integrity.json').write_text(json.dumps({'scientific_artists_before_after_identical':True,'scientific_artist_SHA256':digest,'axes_position_locked_after_original_layout':True,'old_label':'D  Twelve predefined categories','new_label':'D  Twelve operational categories'},indent=2))
 plt.close(fig)
def title(ax,s):ax.set_title(s,loc='left',fontweight='bold',pad=9,fontsize=8)
core=table(7);en=table(8);fig,ax=plt.subplots(2,2,figsize=(6.85,8.5),layout='constrained');a=ax[0,0];a.bar(['Fixed\ncore','Secondary\nstability'],[968,800],color=['#0072B2','#009E73']);a.set(ylim=(0,1100),ylabel='Genes');title(a,'A  Shared subset\nwithin 24,961 genes');a.text(0,990,'968 (3.88%)',ha='center');a.text(1,820,'800 / 968 (82.64%)',ha='center')
a=ax[0,1]
for k,g in enumerate(G):
 vals=core[g+'_log2FC'];v=a.violinplot([vals[vals>0],vals[vals<0]],positions=[k-.16,k+.16],widths=.28,showmedians=True,showextrema=False)
 for b in v['bodies']:b.set_facecolor(C[g]);b.set_alpha(.65)
 v['cmedians'].set_color('black')
a.axhline(0,c='#999999',lw=.7);a.set(xticks=range(4),xticklabels=G,ylabel='Original NB response (log2FC)');title(a,'B  Direction and effect\nconsistency; n = 968')
a=ax[1,0];a.scatter(core.bootstrap_all_four_original_direction_probability,core.dish_omission_original_criteria_retention_fraction,c=np.where(core.robust,'#009E73','#999999'),s=12,alpha=.5);a.axvline(.9,c='black',ls='--');a.axhline(.8,c='black',ls='--');a.set(xlim=(0,1.02),ylim=(0,1.02),xlabel='Original-direction\nbootstrap probability',ylabel='Original-criteria\ndish-omission retention');title(a,'C  Internal stability;\nfixed membership')
a=ax[1,1];cats=list(en.category.unique())
for qi,(query,col,offset) in enumerate([('core_up','#0072B2',-.12),('core_down','#D55E00',.12)]):
 t=en[en['query']==query]
 for j,row in enumerate(t.itertuples()):
  if row.zero_annotated_background:a.text(0,j+offset,'N/A',fontsize=8,va='center');continue
  a.scatter(row.fold_enrichment,j+offset,s=15+row.gene_count*2,c=col,edgecolors='black' if row.FDR_36_predefined_tests<=.05 else 'none',linewidth=1.4)
a.axvline(1,c='#999999',ls='--');a.set(yticks=range(12),yticklabels=[x.replace('_',' ').replace('primary metabolic process','primary metabolic\nprocess').replace('chromatin organization','chromatin\norganization').replace('cell wall organization','cell wall\norganization').replace('meristem maintenance','meristem\nmaintenance') for x in cats],xlabel='Annotated fold enrichment\n(up blue; down orange)');a.tick_params(axis='y',labelsize=8);a.set_title('D  Twelve predefined categories',loc='left',fontweight='bold',pad=40,fontsize=8);a.text(.02,1.015,'Black edge: joint 36-test FDR ≤ 0.05\nSize: overlap',transform=a.transAxes,fontsize=8);save(fig,5)
