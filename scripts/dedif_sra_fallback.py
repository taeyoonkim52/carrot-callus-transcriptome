"""Recover one approved run from checksum-verified SRA Normalized, not SRA Lite."""
import csv,gzip,json,hashlib,os,subprocess,time,urllib.request,shutil,concurrent.futures
from pathlib import Path
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parents[1];W=ROOT/'work/dedifferentiation';O=ROOT/'outputs/dedifferentiation_gate'
RUN='SRR36683244';SAMPLE='C815_20d_3';SIZE=5721440007;MD5='f3bb9714fdf2d446c7459ce26a3e184b';URL=f'https://sra-pub-run-odp.s3.amazonaws.com/sra/{RUN}/{RUN}'
for d in ['sra','sra_conversion','sra_scratch']:(W/d).mkdir(exist_ok=True)
def now():return datetime.now(timezone.utc).isoformat()
def status(stage,**kw):
    tmp=W/'SRA_FALLBACK_STATUS.tmp';tmp.write_text(json.dumps(dict(pid=os.getpid(),stage=stage,updated_utc=now(),**kw),indent=2));os.replace(tmp,W/'SRA_FALLBACK_STATUS.json')
def hashes(path):
    a=hashlib.md5();b=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(8*1024*1024),b''):a.update(chunk);b.update(chunk)
    return a.hexdigest(),b.hexdigest()
def download():
    path=W/'sra'/RUN;partial=Path(str(path)+'.partial')
    if path.exists():
        if path.stat().st_size!=SIZE or hashes(path)[0]!=MD5:raise RuntimeError('Existing SRA archive invalid; preserved')
        return path
    for attempt in range(1,9):
        try:
            offset=partial.stat().st_size if partial.exists() else 0;headers={'User-Agent':'CarrotDataAudit/1.0'}
            if offset:headers['Range']=f'bytes={offset}-'
            if offset<SIZE:
                with urllib.request.urlopen(urllib.request.Request(URL,headers=headers),timeout=90) as r:
                    if offset and (r.status!=206 or not r.headers.get('Content-Range','').startswith(f'bytes {offset}-')):raise RuntimeError('Range rejected; preserved partial')
                    with partial.open('ab' if offset else 'wb') as f:
                        last=time.monotonic()
                        while True:
                            chunk=r.read(1024*1024)
                            if not chunk:break
                            f.write(chunk)
                            if time.monotonic()-last>15:f.flush();status('SRA_DOWNLOAD',bytes=f.tell(),attempt=attempt);last=time.monotonic()
            if partial.stat().st_size!=SIZE:raise RuntimeError('Incomplete SRA size')
            if hashes(partial)[0]!=MD5:raise RuntimeError('SRA checksum mismatch; preserved')
            os.replace(partial,path);return path
        except Exception as e:
            with (W/'logs/sra_fallback.log').open('a') as f:f.write(now()+f' RETRY {attempt}: {e}\n')
            if 'checksum' in str(e).lower() or attempt==8:raise
            time.sleep(min(60,10*attempt))
try:
    status('SRA_DOWNLOAD',bytes=0)
    assert SAMPLE in {r['sample'] for r in csv.DictReader((O/'sample_manifest_24.tsv').open(),delimiter='\t')}
    if shutil.disk_usage(W).free<50_000_000_000:raise RuntimeError('Insufficient 50 GB reserve for streaming SRA conversion; no conversion started')
    archive=download();archive_hashes=hashes(archive)
    toolkit=W/'tools/sratoolkit';exe=next(toolkit.rglob('fastq-dump.exe'));validator=next(toolkit.rglob('vdb-validate.exe'))
    env=os.environ.copy();env['NCBI_SETTINGS']=str(W/'sra_user_settings.mkfg')
    version=subprocess.check_output([str(exe),'--version'],env=env,text=True)
    total_spots=30104493
    commands=[]
    for part in range(4):
        start=part*total_spots//4+1;end=(part+1)*total_spots//4
        folder=W/'sra_conversion'/f'part_{part+1}';folder.mkdir(exist_ok=True)
        commands.append([str(exe),str(archive),'--split-files','--gzip','-N',str(start),'-X',str(end),'-O',str(folder)])
    marker=W/'sra_conversion/CONVERSION_COMPLETE.json'
    if not marker.exists():
        status('SRA_VALIDATION')
        existing_log=W/'logs/sra_validate.log'
        if not existing_log.exists() or "Database 'SRR36683244' is consistent" not in existing_log.read_text():
            with existing_log.open('a') as f:subprocess.run([str(validator),str(archive)],env=env,stdout=f,stderr=f,check=True)
        status('FASTQ_CONVERSION')
        def convert(part):
            folder=W/'sra_conversion'/f'part_{part+1}';done=folder/'DONE.json'
            if done.exists():return
            if any(folder.glob('*.gz')):raise RuntimeError('Incomplete partition preserved; inspect before retry')
            log=W/'logs'/f'fastq_dump_part_{part+1}.log'
            with log.open('a') as f:subprocess.run(commands[part],env=env,stdout=f,stderr=f,check=True)
            import re
            written=[int(x.replace(',','')) for x in re.findall(r'Written ([\d,]+) spots',log.read_text())]
            expected=(part+1)*total_spots//4-part*total_spots//4
            if not written or written[-1]!=expected:raise RuntimeError('Partition spot count mismatch')
            if len(list(folder.glob('*.gz')))!=2:raise RuntimeError('Partition does not contain two biological mates')
            done.write_text(json.dumps({'part':part+1,'first_spot':part*total_spots//4+1,'last_spot':(part+1)*total_spots//4,'written_spots':expected,'command':commands[part]}))
        with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:list(pool.map(convert,range(4)))
        marker.write_text(json.dumps({'run':RUN,'sample':SAMPLE,'expected_spots':total_spots,'archive_md5':archive_hashes[0],'archive_sha256':archive_hashes[1],'software':version,'commands':commands,'partition_policy':'Four contiguous non-overlapping exhaustive spot ranges; concatenate gzip members in range order','completed_utc':now()},indent=2))
    status('FASTQ_AGGREGATION_AND_CRC_CHECK')
    def compress(mate):
        dest=W/'fastq'/f'{SAMPLE}_R{mate}.fq.gz'
        if dest.exists():return dest,hashes(dest)
        temp=Path(str(dest)+'.compression_partial')
        with temp.open('wb') as raw:
            for part in range(4):
                source=W/'sra_conversion'/f'part_{part+1}'/f'{RUN}_{mate}.fastq.gz'
                if not source.exists():raise RuntimeError('Missing converted mate partition '+str(source))
                with source.open('rb') as f:shutil.copyfileobj(f,raw,8*1024*1024)
        with gzip.open(temp,'rb') as f:
            for chunk in iter(lambda:f.read(8*1024*1024),b''):pass
        os.replace(temp,dest);return dest,hashes(dest)
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:results=list(pool.map(compress,[1,2]))
    original=O/'fastq_sources_original_repository.tsv'
    if not original.exists():shutil.copy2(O/'fastq_sources.tsv',original)
    rows=list(csv.DictReader((O/'fastq_sources.tsv').open(),delimiter='\t'))
    for row in rows:
        row['source_type']='ORIGINAL_FASTQ' if row['sample']!=SAMPLE else 'SRA_NORMALIZED_DERIVED_FASTQ'
        if row['sample']==SAMPLE:
            path,h=next((p,h) for p,h in results if p.name==row['filename']);row['bytes']=str(path.stat().st_size);row['md5']=h[0];row['url']='local-derived://'+RUN+'/'+path.name
            with (O/'file_hashes.jsonl').open('a') as f:f.write(json.dumps({'kind':'fastq','path':str(path.relative_to(ROOT)),'bytes':path.stat().st_size,'md5':h[0],'sha256':h[1],'verified_utc':now(),'source_type':row['source_type'],'parent_archive_md5':MD5,'parent_archive_sha256':archive_hashes[1]})+'\n')
    with (O/'fastq_sources.tsv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');writer.writeheader();writer.writerows(rows)
    provenance=json.loads(marker.read_text());provenance.update({'archive_url':URL,'archive_bytes':SIZE,'original_fastq_access':'Original nested S3 source returns HTTP403; ENA filereport HTTP500','original_fastq_cloud_locations':[{'filename':f'{SAMPLE}_R{mate}.fq.gz','url':f's3://sra-pub-src-14/{RUN}/{SAMPLE}_R{mate}.fq.gz.1','repository_md5':md5,'repository_bytes':size} for mate,md5,size in [(1,'9223b2eb7d123d3256cf5248043e6538',3371668357),(2,'01ab0b46a7b08b4e40d42f7f2af3f081',3857300954)]],'normalization':'Full-quality SRA Normalized archive, not quality-binned SRA Lite','fastq_outputs':[{'filename':p.name,'bytes':p.stat().st_size,'md5':h[0],'sha256':h[1]} for p,h in results],'checksum_scope':'46 original FASTQs use repository MD5; two derived FASTQs are locally hashed from a repository-MD5-verified and VDB-validated archive','tool_archive_sha256':hashes(W/'tools/sratoolkit.3.2.1-win64.zip')[1]})
    (O/'sra_fallback_provenance.json').write_text(json.dumps(provenance,indent=2))
    # Keep partition outputs and full-quality archive as retained recovery provenance.
    status('SRA_FALLBACK_COMPLETE',run=RUN,sample=SAMPLE)
except Exception as e:status('FAILED',error=repr(e));raise
