import csv,gzip,hashlib,json,collections
from pathlib import Path
P=Path(__file__).resolve().parents[1];W=P/'work/dedifferentiation/reference';O=P/'outputs/dedifferentiation_gate'
mapping=collections.defaultdict(set)
with gzip.open(W/'GCF_001625215.2_DH1_v3.0_feature_table.txt.gz','rt') as f:
    reader=csv.DictReader(f,delimiter='\t')
    for r in reader:
        if r['product_accession'] and r['GeneID'] and r['# feature']!='CDS':mapping[r['product_accession']].add(r['GeneID'])
seqs={};name=None;pieces=[]
with gzip.open(W/'GCF_001625215.2_DH1_v3.0_rna.fna.gz','rt') as f:
    for line in f:
        if line.startswith('>'):
            if name:seqs[name]=''.join(pieces)
            name=line[1:].split()[0];assert name not in seqs;pieces=[]
        else:pieces.append(line.strip().upper())
    if name:seqs[name]=''.join(pieces)
groups=collections.defaultdict(list)
for name,seq in seqs.items():groups[hashlib.sha256(seq.encode()).hexdigest()].append(name)
with (O/'transcript_to_gene.tsv').open('w',newline='') as f:
    w=csv.writer(f,delimiter='\t');w.writerow(['transcript_id','gene_id'])
    for name in seqs:
        ids=mapping[name]
        if len(ids)!=1:raise RuntimeError('Ambiguous/unmapped transcript '+name+' '+str(ids))
        w.writerow([name,next(iter(ids))])
with (O/'identical_sequence_groups.tsv').open('w',newline='') as f:
    w=csv.writer(f,delimiter='\t');w.writerow(['sequence_sha256','transcripts','gene_ids','cross_gene_ambiguity'])
    for h,names in groups.items():
        if len(names)>1:
            ids=sorted(set().union(*(mapping[n] for n in names)));w.writerow([h,';'.join(names),';'.join(ids),len(ids)>1])
audit={'actual_rna_fasta_records':len(seqs),'actual_mapped_genes':len(set().union(*(mapping[n] for n in seqs))),'report_all_transcripts':60543,'shorter_than_k31':sum(len(s)<31 for s in seqs.values()),'identical_sequence_groups':sum(len(n)>1 for n in groups.values()),'cross_gene_identical_groups':sum(len(n)>1 and len(set().union(*(mapping[x] for x in n)))>1 for n in groups.values()),'policy':'Retain original transcript IDs; flag cross-gene identical groups. Remove cross-gene-identical genes from primary gene-specific inference; report sensitivity separately. No arbitrary reassignment.'}
(O/'reference_content_audit.json').write_text(json.dumps(audit,indent=2));print(json.dumps(audit,indent=2))
rows=list(csv.DictReader((O/'reference_manifest.tsv').open(),delimiter='\t'))
for row in rows:row['actual_rna_fasta_records']=len(seqs);row['actual_mapped_genes']=audit['actual_mapped_genes']
with (O/'reference_manifest.tsv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
