"""External audit harness: enforce Python file/network/import boundaries, leave package intact."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import argparse,json,os,runpy,time,collections,traceback,platform
p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--code',type=Path,required=True);p.add_argument('--data',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--canary',type=Path,required=True);a=p.parse_args()
root=a.root.resolve();code=a.code.resolve();data=a.data.resolve();out=a.output.resolve();log=root/'filesystem_audit.jsonl';handle=log.open('w');counts=collections.Counter();opened=collections.Counter();imports=set();denials=[];phase='self_test';report_fd=os.open(root/'guard_summary.json',os.O_WRONLY|os.O_CREAT|os.O_EXCL,0o600)
runtime_roots=sorted({Path(sys.prefix).resolve(),Path(sys.base_prefix).resolve(),Path(sys.exec_prefix).resolve(),Path(sys.executable).resolve().parents[1]},key=str)
os_roots=[Path(x) for x in ['/usr/lib','/usr/local/lib','/usr/share/locale','/usr/share/zoneinfo','/lib','/lib64','/etc','/proc','/dev']]
def inside(path,parent):return path==parent or path.is_relative_to(parent)
def record(event,path,category,**extra):
    counts[event]+=1;handle.write(json.dumps(dict(phase=phase,event=event,path=str(path),category=category,**extra))+'\n');handle.flush()
def blocked(event,path,why):
    item=dict(phase=phase,event=event,path=str(path),reason=why);denials.append(item);handle.write(json.dumps(dict(category='DENIED',**item))+'\n');handle.flush();raise PermissionError(why+': '+str(path))
def audit(event,args):
    if event=='import':
        name=str(args[0]);imports.add(name)
        if name.split('.')[0] in ['torch','esm']:blocked(event,name,'Model import prohibited')
    if event in ['socket.__new__','socket.connect','socket.getaddrinfo','subprocess.Popen','os.system']:blocked(event,event,'Network/subprocess operation prohibited')
    if event in ['open','sqlite3.connect','os.listdir','os.scandir']:
        raw=args[0]
        if not isinstance(raw,(str,bytes,os.PathLike)):return
        path=Path(os.fsdecode(raw));path=(Path.cwd()/path).resolve() if not path.is_absolute() else path.resolve()
        write=False
        if event=='sqlite3.connect':write=True
        if event=='open':
            mode=args[1];flags=args[2] if len(args)>2 and isinstance(args[2],int) else 0
            write=bool(isinstance(mode,str) and any(c in mode for c in 'wax+') or flags&(os.O_WRONLY|os.O_RDWR|os.O_CREAT|os.O_TRUNC|os.O_APPEND))
        if write:
            if not inside(path,out):blocked(event,path,'Writes allowed only in fresh regenerated output directory')
            category='regenerated_output'
        elif inside(path,code):category='packaged_code'
        elif inside(path,data):category='packaged_data'
        elif inside(path,out):category='regenerated_output'
        elif path==Path(__file__).resolve():category='external_audit_harness'
        elif any(inside(path,q) for q in runtime_roots):category='documented_python_runtime'
        elif any(inside(path,q) for q in os_roots):category='operating_system_runtime'
        else:blocked(event,path,'Read outside unpacked packages and declared runtime roots')
        record(event,path,category,write=write);opened[(str(path),category,write)]+=1
sys.addaudithook(audit)
# A real pre-existing project file must be blocked; its contents are never read.
canary_ok=False
try:a.canary.open('rb')
except PermissionError:canary_ok=True
assert canary_ok,'File-boundary self-test failed'
phase='extraction';start=time.perf_counter();failure=None
try:
    sys.argv=[str(code/'portable/step1.py'),'extract','--data-root',str(data),'--output',str(out),'--expected-anchors','594']
    runpy.run_path(str(code/'portable/step1.py'),run_name='__main__')
except BaseException as e:
    failure=dict(type=type(e).__name__,message=str(e),traceback=traceback.format_exc());raise
finally:
    handle.flush();handle.close()
    summary=dict(succeeded=failure is None,elapsed_seconds=time.perf_counter()-start,python=sys.version,isolated_interpreter=sys.flags.isolated,user_site_disabled=sys.flags.no_user_site,root=str(root),code_root=str(code),data_root=str(data),output_root=str(out),runtime_roots=list(map(str,runtime_roots)),operating_system_read_roots=list(map(str,os_roots)),canary_blocked=canary_ok,denied_events=denials,event_counts=dict(counts),imports=sorted(imports),file_accesses=[dict(path=p,category=c,write=w,count=n) for (p,c,w),n in sorted(opened.items())],failure=failure,scope_note='Python audit-hook enforcement, not kernel syscall tracing or a container. Scientific file opens are restricted; OS/Python shared-library runtime reads are permitted.')
    # The harness report is not an algorithm output and is written after disabling sampling of events.
    phase='harness_report'
    # Existing handle bypass avoids asking the algorithm to write outside its output.
    os.write(report_fd,json.dumps(summary,indent=2).encode()) if 'report_fd' in globals() else None

