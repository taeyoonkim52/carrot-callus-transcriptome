"""Uniform unstranded quantification with read and orientation diagnostics.
Never performs biological inference. Requires checksum ledger entries.
"""
import csv,gzip,hashlib,json,os,shutil,subprocess,time,msvcrt
from pathlib import Path
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/dedifferentiation_gate';WORK=ROOT/'work/dedifferentiation'
for d in ['quant','tools','read_qc','logs']:(WORK/d).mkdir(parents=True,exist_ok=True)
lock=(WORK/'quantification.lock').open('a+b');lock.write(b'0');lock.flush();lock.seek(0)
try:msvcrt.locking(lock.fileno(),msvcrt.LK_NBLCK,1)
except OSError:raise SystemExit('Another quantification worker owns lock')
state={'pid':os.getpid(),'completed':[],'stage':'WAITING_REFERENCE'}
def now():return datetime.now(timezone.utc).isoformat()
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for chunk in iter(lambda:f.read(8*1024*1024),b''):h.update(chunk)
    return h.hexdigest()
def status(**kw):
    state.update(kw);state['updated_utc']=now();tmp=WORK/'QUANTIFICATION_STATUS.tmp';tmp.write_text(json.dumps(state,indent=2));os.replace(tmp,WORK/'QUANTIFICATION_STATUS.json')
def command(cmd,log):
    cmd=[arg.replace('\\','/') for arg in cmd]
    with (OUT/'commands.jsonl').open('a') as f:f.write(json.dumps({'utc':now(),'command':cmd})+'\n')
    with log.open('a') as f:subprocess.run(cmd,stdout=f,stderr=f,check=True)
def run_info(path):
    raw=path.read_text()
    try:return json.loads(raw)
    except json.JSONDecodeError:
        # Native 0.50.1 writes Windows call paths without JSON escaping.
        # Preserve original file; recover literal backslashes for metadata only.
        return json.loads(raw.replace('\\','\\\\'))
def wait_for_ledger(path):
    while True:
        ledger=OUT/'file_hashes.jsonl'
        if ledger.exists():
            entries=[]
            for line in ledger.read_text().splitlines():
                try:entries.append(json.loads(line))
                except json.JSONDecodeError:pass
            hits=[x for x in entries if x['path']==str(path.relative_to(ROOT))]
            if hits and path.exists():
                if sha(path)!=hits[-1]['sha256']:raise RuntimeError('Input hash mismatch '+str(path))
                return hits[-1]
        s=WORK/'EXECUTION_STATUS.json'
        if s.exists() and json.loads(s.read_text()).get('stage')=='FAILED':raise RuntimeError('Acquisition failed; waiting worker stopped')
        time.sleep(15)
try:
    status()
    ref=WORK/'reference/GCF_001625215.2_DH1_v3.0_rna.fna.gz';wait_for_ledger(ref)
    exe=WORK/'tools/kallisto.exe'
    if not exe.exists():raise FileNotFoundError('Install the recorded Kallisto 0.50.1 executable at work/dedifferentiation/tools/kallisto.exe')
    version=subprocess.check_output([str(exe),'version'],text=True).strip()
    if '0.50.1' not in version:raise RuntimeError('Wrong executable version')
    index=WORK/'reference/dh1_v3_rna.idx'
    if not index.exists():command([str(exe),'index','-i',str(index),str(ref)],WORK/'logs/index.log')
    (OUT/'quantification_software.json').write_text(json.dumps({'kallisto':version,'executable_sha256':sha(exe),'index_sha256':sha(index),'reference_sha256':sha(ref),'production_orientation':'unstranded','threads':8,'python':os.sys.version},indent=2))
    files=list(csv.DictReader((OUT/'fastq_sources.tsv').open(),delimiter='\t'))
    samples=list(csv.DictReader((OUT/'sample_manifest_24.tsv').open(),delimiter='\t'))
    assert len(samples)==24 and all(s['day'] in ('0','20') for s in samples)
    # Execution order follows acquisition; analytical matrix order stays in the manifest.
    priority=list(dict.fromkeys(r['sample'] for r in files));samples.sort(key=lambda s:priority.index(s['sample']))
    for sample in samples:
        sid=sample['sample'];status(stage='WAITING_VERIFIED_FASTQS',current_sample=sid)
        mates=sorted([r for r in files if r['sample']==sid],key=lambda x:x['filename']);assert len(mates)==2
        paths=[WORK/'fastq'/r['filename'] for r in mates];inputs=[wait_for_ledger(p) for p in paths]
        target=WORK/'quant'/sid;marker=target/'COMPLETE.json'
        if marker.exists():
            old=json.loads(marker.read_text())
            if any(sha(target/n)!=h for n,h in old['outputs'].items()):raise RuntimeError('Completed quantification mismatch '+sid)
        else:
            if target.exists():raise RuntimeError('Incomplete output preserved; manual inspection required '+sid)
            status(stage='READ_QC',current_sample=sid)
            qcdir=WORK/'read_qc'/sid;qcdir.mkdir(exist_ok=True)
            pilot=[qcdir/'R1.fq.gz',qcdir/'R2.fq.gz'];n=0;quality=[0,0];bases=[0,0];ns=[0,0];q30=[0,0];lengths=[{},{}]
            with gzip.open(paths[0],'rb') as a,gzip.open(paths[1],'rb') as b,gzip.open(pilot[0],'wb') as x,gzip.open(pilot[1],'wb') as y:
                for _ in range(200000):
                    records=[[f.readline() for j in range(4)] for f in (a,b)]
                    if not records[0][0] and not records[1][0]:break
                    if any(not r[3] or not r[0].startswith(b'@') or not r[2].startswith(b'+') or len(r[1].strip())!=len(r[3].strip()) for r in records):raise RuntimeError('Malformed FASTQ '+sid)
                    names=[r[0].split()[0].removesuffix(b'/1').removesuffix(b'/2') for r in records]
                    if names[0]!=names[1]:raise RuntimeError('Mate identifier mismatch '+sid)
                    for k,(r,f) in enumerate(zip(records,(x,y))):
                        seq=r[1].strip();qs=r[3].strip();f.write(b''.join(r));bases[k]+=len(seq);ns[k]+=seq.upper().count(b'N');quality[k]+=sum(q-33 for q in qs);q30[k]+=sum(q>=63 for q in qs);lengths[k][str(len(seq))]=lengths[k].get(str(len(seq)),0)+1
                    n+=1
            if n==0:raise RuntimeError('Empty read QC '+sid)
            orientation={}
            for flag,label in [([], 'unstranded'),(['--fr-stranded'],'FR'),(['--rf-stranded'],'RF')]:
                dest=qcdir/label
                if not (dest/'run_info.json').exists():command([str(exe),'quant','-i',str(index),'-o',str(dest),'-t','4']+flag+[str(p) for p in pilot],WORK/'logs'/f'{sid}_pilot_{label}.log')
                orientation[label]=run_info(dest/'run_info.json')['p_pseudoaligned']
            qc={'sample':sid,'sampled_pairs':n,'sampling':'first 200000 pairs; descriptive prefix QC, not random whole-library QC','mean_phred':[quality[k]/bases[k] for k in range(2)],'q30_fraction':[q30[k]/bases[k] for k in range(2)],'N_fraction':[ns[k]/bases[k] for k in range(2)],'length_histograms':lengths,'pilot_pseudoalignment_percent':orientation,'production_orientation':'unstranded'}
            (qcdir/'read_qc.json').write_text(json.dumps(qc,indent=2))
            if min(qc['mean_phred'])<20:raise RuntimeError('Low prefix read quality; review required '+sid)
            status(stage='QUANTIFYING',current_sample=sid)
            cmd=[str(exe),'quant','-i',str(index),'-o',str(target),'-t','8']+[str(p) for p in paths]
            command(cmd,WORK/'logs'/f'{sid}.log')
            info=run_info(target/'run_info.json')
            if info['n_processed']<=0:raise RuntimeError('Empty quantification '+sid)
            marker.write_text(json.dumps({'sample':sid,'completed_utc':now(),'command':cmd,'inputs':inputs,'outputs':{n:sha(target/n) for n in ['abundance.tsv','abundance.h5','run_info.json']}},indent=2))
        state['completed'].append(sid);status()
    status(stage='QUANTIFICATION_COMPLETE_AWAITING_TECHNICAL_GATE',current_sample=None)
except Exception as e:
    status(stage='FAILED',error=repr(e));raise
finally:
    lock.seek(0);msvcrt.locking(lock.fileno(),msvcrt.LK_UNLCK,1);lock.close()
