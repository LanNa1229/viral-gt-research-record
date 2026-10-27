"""Python-level file-access audit around the authorized representative regression."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import os,json,runpy,time
D=Path(__file__).resolve().parents[1];V=D/'validation';G=D/'github_package/workflows/viral_gt_bgc_prophet_inference';C=V/'runtime';ROOT=D.parents[1];OUTPUT=V/'cold_cache_regression'
allowed=[C,G,Path(sys.prefix),Path(sys.base_prefix)]+[Path(p) for p in ['/lib','/lib64','/usr/lib','/usr/local/lib','/usr/share','/etc','/proc','/sys','/dev']]
accesses=[];denied=[];network=[];canary=False
os.environ['TORCH_HOME']=str(V/'torch_home');os.environ['XDG_CACHE_HOME']=str(V/'xdg_cache');os.environ['PYTHONDONTWRITEBYTECODE']='1'
def under(p,base):return p==base or p.is_relative_to(base)
def guard(event,args):
    global canary
    if event in ['socket.connect','subprocess.Popen','os.system']:
        network.append(event);raise PermissionError('No network/subprocesses during inference')
    if event=='open' and args and isinstance(args[0],(str,bytes,os.PathLike)):
        p=Path(os.fsdecode(args[0])).resolve();mode=args[1];flags=args[2] if len(args)>2 and isinstance(args[2],int) else 0
        write=bool(flags&(os.O_WRONLY|os.O_RDWR|os.O_CREAT|os.O_TRUNC|os.O_APPEND)) or (isinstance(mode,str) and any(x in mode for x in 'wax+'))
        good=(under(p,V) and not under(p,C)) if write else any(under(p,b) for b in allowed)
        accesses.append(dict(path=str(p),write=write,allowed=good))
        if not good:
            if p==ROOT/'AGENTS.md' and canary:raise PermissionError('Canary blocked')
            denied.append(dict(event=event,path=str(p),write=write));raise PermissionError('Undocumented file access: '+str(p))
    if event=='os.listdir' and args and isinstance(args[0],(str,bytes,os.PathLike)):
        p=Path(os.fsdecode(args[0])).resolve()
        if not (any(under(p,b) for b in allowed) or under(p,V)):
            denied.append(dict(event=event,path=str(p)));raise PermissionError('Undocumented directory access: '+str(p))
sys.addaudithook(guard)
canary=True
try:
    (ROOT/'AGENTS.md').read_text();raise RuntimeError('Read guard did not block original workspace')
except PermissionError:pass
canary=False
started=time.monotonic();error=None
try:
    sys.argv=[str(G/'portable/step2.py'),'regress','--cache',str(C),'--output',str(OUTPUT)]
    runpy.run_path(str(G/'portable/step2.py'),run_name='__main__')
except Exception as e:
    error=repr(e)
finally:
    receipt=dict(passed=not denied and not network and error is None,canary_original_workspace_read_blocked=True,denied=denied,network_or_subprocess_events=network,error=error,scope='Python audit hooks; not kernel sandboxing',elapsed_seconds=time.monotonic()-started,accesses=accesses,declared_runtime_roots=[str(x) for x in allowed])
    (V/'FILESYSTEM_GUARD.json').write_text(json.dumps(receipt,indent=2)+'\n')
if error:raise RuntimeError(error)
