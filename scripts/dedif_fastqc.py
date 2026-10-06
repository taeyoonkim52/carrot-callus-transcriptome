"""Full-file read QC worker, restricted to the 48 approved mates."""
import csv,json,os,subprocess,time,msvcrt,shutil
from pathlib import Path
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parents[1];WORK=ROOT/'work/dedifferentiation';OUT=ROOT/'outputs/dedifferentiation_gate';DEST=OUT/'fastqc';DEST.mkdir(exist_ok=True)
lock=(WORK/'fastqc.lock').open('a+b');lock.write(b'0');lock.flush();lock.seek(0)
try:msvcrt.locking(lock.fileno(),msvcrt.LK_NBLCK,1)
except OSError:raise SystemExit('Another FastQC worker owns lock')
state={'pid':os.getpid(),'completed':[],'stage':'WAITING_SOFTWARE'}
def update(**kw):
    state.update(kw);state['updated_utc']=datetime.now(timezone.utc).isoformat();tmp=WORK/'FASTQC_STATUS.tmp';tmp.write_text(json.dumps(state,indent=2));os.replace(tmp,WORK/'FASTQC_STATUS.json')
try:
    update()
    while not (OUT/'read_qc_software.json').exists():time.sleep(15)
    software=json.loads((OUT/'read_qc_software.json').read_text())
    rows=list(csv.DictReader((OUT/'fastq_sources.tsv').open(),delimiter='\t'));assert len(rows)==48
    for row in rows:
        name=row['filename'];path=WORK/'fastq'/name;update(stage='WAITING_VERIFIED_FASTQ',current_file=name)
        while True:
            ledger=OUT/'file_hashes.jsonl';records=[]
            if ledger.exists():
                for line in ledger.read_text().splitlines():
                    try:records.append(json.loads(line))
                    except json.JSONDecodeError:pass
            if path.exists() and any(r['path']==str(path.relative_to(ROOT)) and r['md5']==row['md5'] for r in records):break
            status=WORK/'EXECUTION_STATUS.json'
            if status.exists() and json.loads(status.read_text()).get('stage')=='FAILED':raise RuntimeError('Acquisition failed')
            time.sleep(15)
        basename=name.removesuffix('.gz').removesuffix('.fq').removesuffix('.fastq')
        expected=DEST/(basename+'_fastqc.zip')
        legacy=WORK/'fastq'/(basename+'_fastqc.zip')
        if not expected.exists() and legacy.exists():
            import zipfile
            with zipfile.ZipFile(legacy) as z:
                if z.testzip() is not None:raise RuntimeError('Corrupt prior QC archive')
            shutil.copy2(legacy,expected)
            html=legacy.with_suffix('.html')
            if html.exists():shutil.copy2(html,DEST/html.name)
        if not expected.exists():
            update(stage='FASTQC_RUNNING',current_file=name)
            cmd=[software['java_path'],'-Xmx1024m','-XX:+ExitOnOutOfMemoryError','-Djava.awt.headless=true','-Dfastqc.output_dir='+str(DEST),'-Dfastqc.nogroup=true','-Dfastqc.threads=1','-cp',software['fastqc_path']+';'+software['fastqc_path']+'/*','uk.ac.babraham.FastQC.FastQCApplication',str(path)]
            with (OUT/'commands.jsonl').open('a') as f:f.write(json.dumps({'worker':'FastQC','command':cmd})+'\n')
            with (WORK/'logs'/f'{name}_fastqc.log').open('a') as f:subprocess.run(cmd,stdout=f,stderr=f,check=True)
            if not expected.exists():raise RuntimeError('FastQC did not produce expected file '+name)
        import zipfile
        with zipfile.ZipFile(expected) as z:
            if z.testzip() is not None:raise RuntimeError('Corrupt FastQC archive '+name)
        state['completed'].append(name);update()
    update(stage='FULL_READ_QC_COMPLETE',current_file=None)
except Exception as e:update(stage='FAILED',error=repr(e));raise
finally:lock.seek(0);msvcrt.locking(lock.fileno(),msvcrt.LK_UNLCK,1);lock.close()
