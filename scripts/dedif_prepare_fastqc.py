import urllib.request,json,hashlib,zipfile,subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];W=ROOT/'work/dedifferentiation/tools';W.mkdir(exist_ok=True)
url='https://api.github.com/repos/adoptium/temurin21-binaries/releases/tags/jdk-21.0.8%2B9'
metadata=json.load(urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'CarrotDataAudit/1.0'}),timeout=60))
asset=next(a for a in metadata['assets'] if a['name']=='OpenJDK21U-jre_x64_windows_hotspot_21.0.8_9.zip')
checksum=urllib.request.urlopen(asset['browser_download_url']+'.sha256.txt',timeout=60).read().decode().split()[0]
binary={'name':asset['name'],'link':asset['browser_download_url'],'checksum':checksum}
archive=W/binary['name']
if not archive.exists():urllib.request.urlretrieve(binary['link'],archive)
assert hashlib.sha256(archive.read_bytes()).hexdigest()==binary['checksum']
javafolder=W/'java'
if not javafolder.exists():
    with zipfile.ZipFile(archive) as z:z.extractall(javafolder)
java=next(javafolder.glob('*/bin/java.exe'))
q=W/'fastqc_v0.12.1.zip'
if not q.exists():urllib.request.urlretrieve('https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip',q)
if not (W/'FastQC').exists():
    with zipfile.ZipFile(q) as z:z.extractall(W)
record={'java_download_metadata_url':url,'java_archive_url':binary['link'],'java_archive_sha256':binary['checksum'],'java_version':subprocess.check_output([str(java),'-version'],stderr=subprocess.STDOUT,text=True),'fastqc_url':'https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip','fastqc_archive_sha256':hashlib.sha256(q.read_bytes()).hexdigest(),'fastqc_version':subprocess.check_output([str(java),'-Dfastqc.show_version=true','-cp',str(W/'FastQC'),'uk.ac.babraham.FastQC.FastQCApplication'],text=True),'java_path':str(java),'fastqc_path':str(W/'FastQC')}
(ROOT/'outputs/dedifferentiation_gate/read_qc_software.json').write_text(json.dumps(record,indent=2));print(record['java_version'],record['fastqc_version'])
