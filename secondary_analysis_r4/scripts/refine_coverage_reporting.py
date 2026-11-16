from pathlib import Path
import pandas as pd,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
RETAINED=ROOT/'secondary_analysis_r4'
O=ROOT/'work/r4_reproduction'
O.mkdir(parents=True,exist_ok=True)
R=ROOT
c=pd.read_csv(O/'R4A_ANNOTATION_COVERAGE_COMPARISON.tsv',sep='\t');g=pd.read_csv(O/'R4A_GENE_ANNOTATION_CANDIDATE.tsv',sep='\t',dtype=str).fillna('').set_index('gene_id');p=pd.read_csv(ROOT/'results/supplementary_tables/S4_interaction_statistics.tsv',sep='\t',dtype={'gene_id':str}).set_index('gene_id');old=pd.read_csv(ROOT/'secondary_analysis_r1_r2/ortholog_mapping/R1_ORTHOLOG_MAPPING.tsv',sep='\t',dtype=str).fillna('').set_index('gene_id');panel=pd.read_csv(ROOT/'secondary_analysis_r1_r2/canonical_regulators/R1_CANONICAL_REGULATOR_PANEL.tsv',sep='\t',dtype=str).fillna('');newpanel=pd.read_csv(O/'R4A_CANONICAL_REGULATOR_REANNOTATION.tsv',sep='\t',dtype=str).fillna('')
sets={'UNIVERSE':set(p.index),'CORE_ALL':set(p.index[p.conserved_core]),'CORE_UP':set(p.index[p.conserved_core&(p.balanced_time_log2FC>0)]),'CORE_DOWN':set(p.index[p.conserved_core&(p.balanced_time_log2FC<0)]),'PRACTICAL_INTERACTION':set(p.index[(p.interaction_FDR<=.05)&(p.response_effect_range_log2>=.5)]),'LARGE_INTERACTION':set(p.index[(p.interaction_FDR<=.05)&(p.response_effect_range_log2>=1)]),'CANONICAL_OLD_CARROT_CANDIDATES':set(panel.carrot_gene_id)-{''}}
c['ambiguous_mappings']=0;c['one_to_many_accepted']=0
for ix,row in c.iterrows():
 layer=row['layer'];name=row['set']
 if name=='CANONICAL_PANEL_23_REGULATORS':
  # This row counts regulators; all other rows count carrot genes.
  c.loc[ix,'ambiguous_mappings']=int((panel.mapping_confidence!='HIGH').sum()-(panel.mapping_confidence=='UNRESOLVED').sum()) if layer=='OLD' else int((newpanel.resolution=='AMBIGUOUS_FAMILY_HOMOLOGY').sum())
  continue
 genes=sets[name];amb_old=set(old.index[old.bridge_status=='AMBIGUOUS_CROSS_VERSION_BRIDGE'])&genes
 amb_new={x for x in genes if 'AMBIGUITY' in g.loc[x,'protein_status'] or 'NEAR_TIE' in g.loc[x,'protein_status']}
 mapped_candidate={x for x in genes if g.loc[x,'candidate_AGI']}
 c.loc[ix,'ambiguous_mappings']=len(amb_old) if layer=='OLD' else len(amb_new) if layer=='DIRECT_RBH' else len((amb_old|amb_new)-mapped_candidate)
 col='OLD_AGI' if layer=='OLD' else 'DIRECT_RBH_AGI' if layer=='DIRECT_RBH' else 'candidate_AGI'
 c.loc[ix,'one_to_many_accepted']=sum(';' in g.loc[x,col] for x in genes)
c.to_csv(O/'R4A_ANNOTATION_COVERAGE_COMPARISON.tsv',sep='\t',index=False)
s=json.loads((O/'R4A_ANNOTATION_SUMMARY.json').read_text());s['coverage']=c.to_dict('records');s['coverage_ambiguity_note']='ambiguous_mappings is layer-specific; ambiguous_direct is the same diagnostic direct-search uncertainty count for comparison, not the old bridge ambiguity count. Canonical rows count 23 regulators rather than carrot genes; accepted one-to-many GO unions are separately reported.'
(O/'R4A_ANNOTATION_SUMMARY.json').write_text(json.dumps(s,indent=2))
print('Layer-specific ambiguity counts clarified; no annotation, enrichment or genotype join rerun.')
