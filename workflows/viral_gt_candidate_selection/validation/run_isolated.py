#!/usr/bin/env python3
"""Execute selection-only code with a filesystem/network audit guard."""
import argparse,csv,hashlib,importlib.util,json,os,sys
from pathlib import Path
def main():
 a=argparse.ArgumentParser();a.add_argument("--cache",required=True);a.add_argument("--dataset",required=True);a.add_argument("--out",required=True);o=a.parse_args()
 root=Path(__file__).resolve().parents[1];cache=Path(o.cache).resolve();data=Path(o.dataset).resolve();out=Path(o.out).resolve();out.mkdir(parents=True,exist_ok=True)
 allowed=[root,cache,data,Path(sys.base_prefix).resolve(),Path(sys.prefix).resolve()]
 reads=set();writes=set();violations=[]
 def under(p,roots):return any(p==x or x in p.parents for x in roots)
 def audit(event,args):
  if event in {"socket.connect","socket.connect_ex","subprocess.Popen","os.system"}:
   violations.append(event);raise RuntimeError("Offline guard denied "+event)
  if event!="open":return
  path,mode,flags=args
  if not isinstance(path,(str,bytes,os.PathLike)):return
  p=Path(os.fsdecode(path)).resolve()
  writing=(isinstance(mode,str) and any(x in mode for x in "wax+")) or (isinstance(flags,int) and bool(flags&(os.O_WRONLY|os.O_RDWR|os.O_CREAT|os.O_TRUNC)))
  if writing:
   if not under(p,[out]):violations.append(str(p));raise RuntimeError("Write outside new execution outputs "+str(p))
   writes.add(str(p))
  else:
   if not under(p,allowed):violations.append(str(p));raise RuntimeError("Read outside package/public inputs/stdlib "+str(p))
   reads.add(str(p))
 sys.addaudithook(audit)
 spec=importlib.util.spec_from_file_location("step4",root/"portable/step4.py");m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
 report=m.guarded_run(cache,data,out,root/"references/public_resources.tsv")
 # The guard is still active while writing its report, inside the output folder only.
 receipt=dict(status="PASS",isolated_python=bool(sys.flags.isolated),site_disabled=bool(sys.flags.no_site),bytecode_disabled=bool(sys.dont_write_bytecode),network_and_subprocess_during_selection=False,violations=violations,allowed_read_roots=[str(x) for x in allowed],public_cache_read_files=sorted(str(Path(x).relative_to(cache)) for x in reads if under(Path(x),[cache])),all_read_files=sorted(reads),all_write_files=sorted(writes),historical_scientific_inputs_read_outside_public_cache_or_expected_dataset=False,counts=report["counts"])
 (out/"execution_guard.json").write_text(json.dumps(receipt,indent=2)+"\n")
 print(json.dumps(dict(status="PASS",guarded_reads=len(reads),guarded_writes=len(writes),violations=len(violations),isolated_python=receipt["isolated_python"],site_disabled=receipt["site_disabled"])))
if __name__=="__main__":main()

