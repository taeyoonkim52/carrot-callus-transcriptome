"""Typography/layout-only reuse of frozen figure renderer and source data."""
from pathlib import Path
import re
P=Path(__file__).resolve().parent;R=P.parents[1]
source=R/'publication/scripts/render_publication.py'
s=source.read_text(encoding='utf8').split("(P/'FIGURE_CAPTIONS.md').write_text")[0]
s=s.replace("P/f'figures/", "R/f'publication/submission_plant_cell_reports/final_submission/Figures/")
s=s.replace("'font.size':10", "'font.size':8")
s=s.replace('figsize=(11,8)','figsize=(6.85,6.5)').replace('figsize=(12,11)','figsize=(6.85,9.0)').replace('figsize=(12,4)','figsize=(6.85,7.5)').replace('figsize=(13,4.6)','figsize=(6.85,8.5)').replace('figsize=(12,9)','figsize=(6.85,8.5)')
s=s.replace('plt.subplots(1,3,','plt.subplots(3,1,')
s=s.replace('fontweight=\'bold\',pad=12','fontweight=\'bold\',pad=9,fontsize=9')
s=s.replace('fontweight=\'bold\',pad=30','fontweight=\'bold\',pad=18,fontsize=9')
s=re.sub(r'fontsize=(?:7|9|10|12)(?![0-9])','fontsize=8',s)
s=s.replace('labelsize=7','labelsize=8').replace('labelsize=8','labelsize=8')
s=s.replace('gs=fig.add_gridspec(3,3)','gs=fig.add_gridspec(4,2)')
s=s.replace('gs[1+k//3,k%3]','gs[1+k//2,k%2]')
s=s.replace("f'{chr(67+k)}  r = {t.pearson_r:.3f}; ρ = {t.spearman_r:.3f}'", "f'{chr(67+k)}  r = {t.pearson_r:.3f}\\nρ = {t.spearman_r:.3f}'")
s=s.replace("fig.colorbar(h[-1],ax=fig.axes[2:],shrink=.6,label=", "fig.colorbar(h[-1],ax=fig.axes[2:],orientation='horizontal',fraction=.025,pad=.025,shrink=.8,label=")
s=s.replace("'A  Pearson response similarity'", "'A  Pearson response\\nsimilarity'")
s=s.replace("['Fixed core','Secondary stability']", "['Fixed\\ncore','Secondary\\nstability']")
s=s.replace("'Original-direction bootstrap probability'", "'Original-direction\\nbootstrap probability'")
s=s.replace("'Original-criteria dish-omission retention'", "'Original-criteria\\ndish-omission retention'")
s=s.replace("[x.replace('_',' ') for x in cats]", "[x.replace('_',' ').replace('primary metabolic process','primary metabolic\\nprocess').replace('chromatin organization','chromatin\\norganization').replace('cell wall organization','cell wall\\norganization').replace('meristem maintenance','meristem\\nmaintenance') for x in cats]")
s=s.replace("'Annotated fold enrichment (up blue; down orange)'", "'Annotated fold enrichment\\n(up blue; down orange)'")
s=s.replace("'Black edge: joint 36-test FDR ≤ 0.05; size: overlap'", "'Black edge: joint 36-test FDR ≤ 0.05\\nSize: overlap'")
s=s.replace("'Mean squared centroid distance per gene'","'Mean squared centroid distance\nper gene'") if False else s
s=s.replace("'C  15,053 / 24,961 (60.31%)'","'C  15,053 / 24,961 (60.31%)'")
# Reflow only long presentation titles; statistical quantities remain unmodified.
for a,b in {
 'B  Direction and effect consistency; n = 968':'B  Direction and effect\nconsistency; n = 968',
 'C  Internal stability; fixed membership':'C  Internal stability;\nfixed membership',
 'C  Practical interactions after omission':'C  Practical interactions\nafter omission',
 'B  All six pairs; 24,961 genes':'B  Six pairs; 24,961 genes',
 'A  Shared subset within 24,961 genes':'A  Shared subset\nwithin 24,961 genes',
 'D  No established global convergence':'D  No established\nglobal convergence',
 'C  Global variance partition':'C  Global variance\npartition'
}.items():s=s.replace(a,b.replace('\n','\\n'))
s=s.replace("a.set_title('D  Twelve predefined categories',loc='left',fontweight='bold',pad=18,fontsize=8)", "a.set_title('D  Twelve predefined categories',loc='left',fontweight='bold',pad=40,fontsize=8)")
(P/'final_submission/Figures').mkdir(parents=True,exist_ok=True)
exec(compile(s,str(source),'exec'),{'__file__':str(source),'__name__':'__main__'})
