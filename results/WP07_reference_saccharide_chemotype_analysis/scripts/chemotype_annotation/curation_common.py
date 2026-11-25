"""Standard-library I/O shared by the curation/validation scripts."""
import collections,csv,datetime,hashlib,json,os,re,sys,unicodedata
from pathlib import Path
OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[1]
PRIVATE=Path(r'/redacted-home/projects\REDACTED_PRIVATE_REPOSITORY') if os.name=='nt' else Path('/redacted-home/projects/REDACTED_PRIVATE_REPOSITORY')
SCOPE='results/WP07_viral_gt_figure3C_clustering_20261024'
ARCHIVE=PRIVATE/SCOPE
VIRAL=ROOT/'analysis/viral_gt_community_data_integration_20261027/viral_287_integrated_community_mapping.tsv'
COMMIT='284d4f1655cf036c1890bdffcb4aa2dc93f6e41e'
NR='NOT_REPORTED'
def jtext(x):return json.dumps(x,ensure_ascii=False,separators=(',',':'))
def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def read_table(p,sep=None):
 p=Path(p)
 with p.open(encoding='utf-8-sig',newline='') as f:
  reader=csv.DictReader(f,delimiter=sep or ('\t' if p.suffix=='.tsv' else ','))
  result=list(reader)
  assert reader.fieldnames and len(reader.fieldnames)==len(set(reader.fieldnames))
  assert all(None not in r and None not in r.values() for r in result)
  return result
def write_table(name,rows,fields=None):
 p=OUT/name
 with p.open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields or list(rows[0]),delimiter='\t',lineterminator='\n')
  w.writeheader();w.writerows(rows)
def read_json(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def write_json(name,obj):(OUT/name).write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
def as_cell(v):
 if v is None or v=='' or v==[] or v=={}:return NR
 if isinstance(v,(dict,list)):return jtext(v)
 return str(v)
def unique(rows,key):
 result={}
 for r in rows:
  assert r[key] and r[key] not in result,(key,r[key])
  result[r[key]]=r
 return result
def load_mibig():
 import tarfile
 result={}
 with tarfile.open(OUT/'raw_sources/mibig_json_3.1.tar.gz','r:gz') as tar:
  for m in tar:
   if m.isfile() and m.name.endswith('.json'):
    b=tar.extractfile(m).read();a=Path(m.name).stem
    assert a not in result
    result[a]=dict(raw=json.loads(b),sha256=hashlib.sha256(b).hexdigest(),member=m.name)
 return result
def recursive_evidence(x):
 result=[]
 if isinstance(x,dict):
  for k,v in x.items():
   if k=='evidence' and isinstance(v,list):result+=v
   else:result+=recursive_evidence(v)
 elif isinstance(x,list):
  for v in x:result+=recursive_evidence(v)
 return sorted(set(str(v) for v in result))
