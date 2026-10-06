"""Resumable acquisition only. No expression inference or biological testing."""
import csv, hashlib, json, os, time, urllib.request, msvcrt,subprocess,sys
from pathlib import Path
from datetime import datetime, timezone
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/dedifferentiation_gate'
WORK=ROOT/'work/dedifferentiation'
for folder in ('reference','fastq','logs'): (WORK/folder).mkdir(parents=True,exist_ok=True)
lock=(WORK/'acquisition.lock').open('a+b');lock.write(b'0');lock.flush();lock.seek(0)
try: msvcrt.locking(lock.fileno(),msvcrt.LK_NBLCK,1)
except OSError: raise SystemExit('Another acquisition owns the lock')
state={'pid':os.getpid(),'stage':'REFERENCE','verified':[]}
def now(): return datetime.now(timezone.utc).isoformat()
def log(message):
    with (WORK/'logs/acquisition.log').open('a',encoding='utf8') as f:f.write(now()+' '+message+'\n')
def status(**changes):
    state.update(changes);state['updated_utc']=now()
    temp=WORK/'EXECUTION_STATUS.tmp';temp.write_text(json.dumps(state,indent=2));os.replace(temp,WORK/'EXECUTION_STATUS.json')
def hashes(path):
    a=hashlib.md5();b=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(8*1024*1024),b''):a.update(chunk);b.update(chunk)
    return a.hexdigest(),b.hexdigest()
def download(url,path,expected=None,md5=None):
    if path.exists():
        if expected is not None and path.stat().st_size!=expected:raise RuntimeError('Existing size mismatch; preserved '+str(path))
        a,b=hashes(path)
        if md5 and a!=md5:raise RuntimeError('Existing checksum mismatch; preserved '+str(path))
        log('REUSE VERIFIED '+path.name);return a,b
    partial=Path(str(path)+'.partial')
    for attempt in range(1,9):
        try:
            offset=partial.stat().st_size if partial.exists() else 0
            if expected is None or offset<expected:
                headers={'User-Agent':'CarrotDataAudit/1.0'}
                if offset:headers['Range']=f'bytes={offset}-'
                with urllib.request.urlopen(urllib.request.Request(url,headers=headers),timeout=90) as response:
                    if offset and (response.status!=206 or not response.headers.get('Content-Range','').startswith(f'bytes {offset}-')):raise RuntimeError('Range not honored; preserved partial')
                    with partial.open('ab' if offset else 'wb') as f:
                        last=time.monotonic()
                        while True:
                            chunk=response.read(1024*1024)
                            if not chunk:break
                            f.write(chunk)
                            if time.monotonic()-last>=15:
                                f.flush();status(current_download=path.name,download_bytes=f.tell());last=time.monotonic()
            if expected is not None and partial.stat().st_size!=expected:raise RuntimeError('Incomplete download size')
            a,b=hashes(partial)
            if md5 and a!=md5:raise RuntimeError('Checksum mismatch; partial preserved for inspection')
            os.replace(partial,path);log('DOWNLOAD VERIFIED '+path.name);return a,b
        except Exception as e:
            log(f'RETRY {attempt} {path.name}: {type(e).__name__}: {e}')
            if 'checksum' in str(e).lower() or attempt==8:raise
            time.sleep(min(60,10*attempt))
def record(kind,path,a,b):
    with (OUT/'file_hashes.jsonl').open('a') as f:f.write(json.dumps(dict(kind=kind,path=str(path.relative_to(ROOT)),bytes=path.stat().st_size,md5=a,sha256=b,verified_utc=now()))+'\n')
try:
    status()
    base='https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/79200/GCF_001625215.2-RS_2024_03/'
    prefix='GCF_001625215.2_DH1_v3.0'
    path=WORK/'reference/md5checksums.txt';a,b=download(base+path.name,path);record('reference',path,a,b)
    checks={line.split()[1].lstrip('./'):line.split()[0] for line in path.read_text().splitlines() if len(line.split())==2}
    names=[prefix+'_rna.fna.gz',prefix+'_genomic.gff.gz',prefix+'_feature_table.txt.gz','GCF_001625215.2-RS_2024_03_annotation_report.xml','GCF_001625215.2-RS_2024_03_gene_ontology.gaf.gz']
    refs=[]
    for name in names:
        if name not in checks:raise RuntimeError('Official checksum not found: '+name)
        path=WORK/'reference'/name;a,b=download(base+name,path,md5=checks[name]);record('reference',path,a,b)
        refs.append(['GCF_001625215.2 DH1 v3.0','GCF_001625215.2-RS_2024_03','2024-03-27',33079,39086,60543,name,base+name,a,b])
    with (OUT/'reference_manifest.tsv').open('w',newline='') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['assembly','annotation','release_date','protein_coding_genes','genes_and_pseudogenes','transcripts','file','url','md5','sha256']);w.writerows(refs)
    subprocess.run([sys.executable,str(ROOT/'scripts/dedif_reference_audit.py')],check=True)
    rows=list(csv.DictReader((OUT/'fastq_sources.tsv').open(),delimiter='\t'))
    samples=list(csv.DictReader((OUT/'sample_manifest_24.tsv').open(),delimiter='\t'))
    assert len(rows)==48 and len(samples)==24 and all(r['day'] in ('0','20') for r in samples)
    import shutil
    status(stage='FASTQ_ACQUISITION')
    for row in rows:
        path=WORK/'fastq'/row['filename'];status(current_sample=row['sample'],current_download=path.name,download_bytes=0)
        if row['url'].startswith('local-derived://') and not path.exists():raise RuntimeError('Derived mate missing; rerun SRA recovery rather than downloading archive as FASTQ')
        partial=Path(str(path)+'.partial');remaining=int(row['bytes'])-(partial.stat().st_size if partial.exists() else 0)
        if not path.exists() and shutil.disk_usage(WORK).free<remaining+25_000_000_000:raise RuntimeError('25 GB disk reserve would be breached')
        a,b=download(row['url'],path,int(row['bytes']),row['md5']);record('fastq',path,a,b)
        state['verified'].append(path.name);status(current_download=None)
    status(stage='ACQUISITION_COMPLETE_AWAITING_READ_QC',current_sample=None,current_download=None)
except Exception as e:
    log('FAILED '+repr(e));status(stage='FAILED',error=repr(e));raise
finally:
    lock.seek(0);msvcrt.locking(lock.fileno(),msvcrt.LK_UNLCK,1);lock.close()
