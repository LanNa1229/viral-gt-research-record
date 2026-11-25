"""Pin sources and inspect existing MIBiG records; no models or clustering."""
import ast,collections,csv,datetime,hashlib,importlib.util,io,json,os,shutil,struct,sys,tarfile
from pathlib import Path
ROOT=Path(r'/redacted-home/260908_MiBG') if os.name=='nt' else Path('/redacted-home/260908_MiBG')
OUT=ROOT/'analysis/viral_gt_reference_chemotype_curation_20261027'
PRIVATE=Path(r'/redacted-home/projects\REDACTED_PRIVATE_REPOSITORY') if os.name=='nt' else Path('/redacted-home/projects/REDACTED_PRIVATE_REPOSITORY')
SCOPE='results/WP07_viral_gt_figure3C_clustering_20261024/'
def digest(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def read(p,sep='\t'):
 with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter=sep))
def write(p,rows):
 with p.open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
ev=json.loads((OUT/'PRIVATE_ARCHIVE_IDENTITY.json').read_text())
checks=[]
for r in ev['files']:
 p=PRIVATE/r['path'];b=p.read_bytes()
 assert hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()==r['sha']
 checks.append(dict(path=str(p),sha256=digest(p),size_bytes=len(b),git_blob_sha=r['sha']))
metadata=read(PRIVATE/SCOPE/'inputs/workspace/analysis/viral_glycan_benchmark_20261016/mibig_matched_reference_metadata.tsv')
assert len(metadata)==2502 and len({r['BGC'] for r in metadata})==2502
assert [r['embedding_row'] for r in metadata]==[str(i) for i in range(2502)]
assert all(r['n_biological_reference_rows']=='1' and r['n_source_augmentations']=='5' for r in metadata)
global_rows=read(PRIVATE/SCOPE/'assignments/global_joint_leiden_grid.csv',',')
global_map={r['sample_id']:r['k20_r1_s42'] for r in global_rows}
assert len(global_rows)==len(global_map)==2789
for r in metadata:r['GC']=global_map[r['BGC']]
sac={r['BGC'] for r in metadata if 'Saccharide' in r['label'].split()}
g8={r['BGC'] for r in metadata if r['GC']=='GC08'};g5={r['BGC'] for r in metadata if r['GC']=='GC05'}
priority=sac|g8|g5
assert (len(sac),len(g8),len(g8&sac),len(g5),len(g5&sac),len(sac-g8),len(priority))==(195,105,100,53,24,95,229)
emb=PRIVATE/SCOPE/'inputs/workspace/analysis/viral_glycan_benchmark_20261016/mibig_CORE_ONLY_embeddings.npy'
with emb.open('rb') as f:
 assert f.read(6)==b'\x93NUMPY';major,minor=f.read(2)
 n=struct.unpack('<H' if major==1 else '<I',f.read(2 if major==1 else 4))[0]
 header=ast.literal_eval(f.read(n).decode('latin1'))
assert header['shape']==(2502,320)
source=ROOT/'sources/step1b/downloads/mibig_json_3.1.tar.gz'
download_manifest=ROOT/'sources/step1b/additional_download_manifest.json'
dm=json.loads(download_manifest.read_text())
proof=next(r for r in dm if r['name']==source.name)
assert digest(source)==proof['sha256']=='87221d7ead87709fbcfe81176d4f1684520017e1d30d5f1f24ef7a39310d1e94'
rawdir=OUT/'raw_sources';rawdir.mkdir()
shutil.copyfile(source,rawdir/source.name);shutil.copyfile(download_manifest,rawdir/'original_download_manifest.json')
json_records={};record_bytes={}
with tarfile.open(source,'r:gz') as tar:
 for member in tar.getmembers():
  if member.isfile() and member.name.endswith('.json'):
   b=tar.extractfile(member).read();j=json.loads(b)
   accession=Path(member.name).stem
   assert accession not in json_records
   json_records[accession]=j;record_bytes[accession]=b
assert set(r['BGC'] for r in metadata)<=set(json_records)
jdir=rawdir/'mibig_json_3.1';jdir.mkdir()
for r in metadata:(jdir/(r['BGC']+'.json')).write_bytes(record_bytes[r['BGC']])
gbk=ROOT/'sources/mibig_gbk_3.1.tar.gz'
shutil.copyfile(gbk,rawdir/gbk.name)
source_info=dict(status='PASS',snapshot='MIBiG 3.1',mibig_json_url=proof['url'],json_original_retrieval_utc=proof['download_utc'],
 json_sha256=proof['sha256'],gbk_url='https://dl.secondarymetabolites.org/mibig/mibig_gbk_3.1.tar.gz',
 gbk_sha256=digest(gbk),gbk_original_retrieval_utc='NOT_RECORDED_IN_INSPECTED_MANIFEST',
 json_records_in_snapshot=len(json_records),reference_matches=2502,priority_records=229,
 snapshot_reused_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
 source_files=checks,embedding_header=header,biopython_available=importlib.util.find_spec('Bio') is not None)
(OUT/'SOURCE_PREFLIGHT.json').write_text(json.dumps(source_info,indent=2)+'\n',encoding='utf-8')
cards=[]
for r in metadata:
 if r['BGC'] in priority:
  j=json_records[r['BGC']]
  cards.append(dict(accession=r['BGC'],GC=r['GC'],original_label=r['label'],raw_record=j))
(OUT/'PRIORITY_SOURCE_RECORDS.json').write_text(json.dumps(cards,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
summary=[]
for r in metadata:
 if r['BGC'] in priority:
  j=json_records[r['BGC']];c=j.get('cluster',j)
  products=c.get('compounds',[])
  summary.append(dict(mibig_accession=r['BGC'],community=r['GC'],original_label=r['label'],
   products=json.dumps(products,ensure_ascii=False),biosyn_class=json.dumps(c.get('biosyn_class'),ensure_ascii=False),
   status=c.get('mibig_accession',c.get('accession','')),record_keys=';'.join(c.keys())))
write(OUT/'PRIORITY_PRODUCT_REVIEW.tsv',summary)
print(json.dumps(source_info,indent=2))
for accession in ['BGC0000727',metadata[0]['BGC']]:
 print(accession,json.dumps(json_records[accession],ensure_ascii=False)[:18000])
print('RECORD KEY COUNTS',dict(collections.Counter(k for j in json_records.values() for k in j)))
