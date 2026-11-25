"""Shared I/O and explicit runtime paths. No fitting at import."""
import csv,hashlib,json,os,sys
from pathlib import Path
OUT=Path(__file__).resolve().parents[1]
CONFIG=json.loads((OUT/'config.json').read_text(encoding='utf-8-sig'))
ROOT=Path(CONFIG['workspace_root'])
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):os.environ[k]='1'
os.environ['MPLCONFIGDIR']=str(OUT/'runtime/matplotlib_cache')
sys.path.insert(0,str(ROOT/CONFIG['analysis_dependencies']));sys.path.insert(0,CONFIG['graph_runtime'])
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def read(p,sep=','):
 with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter=sep))
def write(name,rows,fields=None):
 p=OUT/name;p.parent.mkdir(parents=True,exist_ok=True)
 with p.open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields or list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
def save(name,value):
 p=OUT/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8')
