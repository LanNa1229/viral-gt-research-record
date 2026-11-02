"""Read historical Ubuntu sources; write only to the explicitly supplied staging root."""
import argparse, csv, hashlib, json, os, platform, shutil, subprocess, sys
from datetime import datetime, timezone
from pathlib import Path

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def dump(p, obj):
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(obj, indent=2, default=str) + '\n')

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--stage',type=Path); ap.add_argument('--probe',action='store_true'); ap.add_argument('--validate',action='store_true'); ap.add_argument('--supplement',action='store_true'); ap.add_argument('--verify',action='store_true'); ap.add_argument('--correct',action='store_true'); ap.add_argument('--validate-package',action='store_true')
    a=ap.parse_args()
    if a.probe:
        home=Path('/home/rahnn')
        print('HOME DIRECTORIES', [p.name for p in home.iterdir() if p.is_dir()])
        for pattern in ['*/bin/python*','*/envs/*/bin/python*','.local/lib/python*/site-packages/matplotlib-*.dist-info','*/lib/python*/site-packages/matplotlib-*.dist-info','*/envs/*/lib/python*/site-packages/matplotlib-*.dist-info','2606*Foldseek*/**/*manifest*','2606*Foldseek*/**/*metadata*']:
            print('PATTERN',pattern)
            print('\n'.join(str(p) for p in home.glob(pattern)))
        return
    out=a.stage.resolve(); src=Path('/home/rahnn/260630_Temp')
    if not str(out).startswith('/mnt/d/Projects/Viral GT cluster 8/'): raise ValueError(out)
    if a.correct or a.validate_package:
        root=out.parent/'WP08_release_candidate/workflows/WP08_cluster8_gene_neighborhood_heatmap' if a.validate_package else out
        dest=out.parent/'WP08_package_final_validation' if a.validate_package else out/'corrected_interpretation_linux'
        env=dict(os.environ,MPLBACKEND='Agg',PYTHONDONTWRITEBYTECODE='1',MPLCONFIGDIR=str(out/'correction_mplconfig'))
        script='reproduce.py' if a.validate_package else 'correct_interpretation.py'
        cmd=['/home/rahnn/miniforge3/bin/python','-B',str(root/script),'--root',str(root),'--out',str(dest)]
        if a.validate_package: cmd+=['--require-exact-render']
        r=subprocess.run(cmd,env=env); sys.exit(r.returncode)
    if a.verify:
        provenance=json.loads((out/'audit/provenance.json').read_text())
        results=[dict(source=r['source'],sha256_unchanged=sha(Path(r['source']))==r['sha256'],copy_matches=sha(out/r['archive_path'])==r['sha256'],mtime_unchanged=Path(r['source']).stat().st_mtime_ns==r['source_mtime_ns']) for r in provenance['files']]
        dump(out/'audit/source_preservation_final.json',results)
        assert all(r['sha256_unchanged'] and r['copy_matches'] and r['mtime_unchanged'] for r in results)
        print('Verified original source bytes, staged bytes and source mtimes:',len(results))
        return
    if a.validate:
        env=dict(os.environ,MPLBACKEND='Agg',PYTHONDONTWRITEBYTECODE='1',MPLCONFIGDIR=str(out/'validation_linux/mplconfig'))
        cmd=['/home/rahnn/miniforge3/bin/python','-B',str(out/'validate_recovery.py'),'--root',str(out),'--out',str(out/'validation_linux_final')]
        r=subprocess.run(cmd,env=env); sys.exit(r.returncode)
    out.mkdir(parents=True,exist_ok=a.supplement)
    records=json.loads((out/'audit/provenance.json').read_text())['files'] if a.supplement else []
    def preserve(p, rel):
        st=p.stat(); h=sha(p); dest=out/'originals'/rel
        dest.parent.mkdir(parents=True,exist_ok=True); shutil.copy2(p,dest)
        assert sha(dest)==h and sha(p)==h
        records.append(dict(source=str(p),archive_path=str(dest.relative_to(out)),sha256=h,
                            bytes=st.st_size,source_mtime_ns=st.st_mtime_ns,
                            source_mtime_utc=datetime.fromtimestamp(st.st_mtime,timezone.utc).isoformat(),
                            source_resolved=str(p.resolve()),symlink=p.is_symlink()))
    for p in ([] if a.supplement else sorted(src.rglob('*'))):
        if p.is_file() and ':Zone.Identifier' not in p.name:
            preserve(p,Path('260630_Temp')/p.relative_to(src))
    scripts=[]
    patterns=['foldseek_manifest*.tsv','foldseek_manifest*.tsv.bak'] if a.supplement else ['foldseek_batch_api_v2.py','analyze_foldseek_results_v2.py']
    for pattern in patterns:
        for p in sorted(Path('/home/rahnn').glob('*/'+pattern)):
            preserve(p,Path('related_scripts')/p.parent.name/p.name)
            scripts.append(dict(path=str(p),sha256=sha(p),bytes=p.stat().st_size))
    dump(out/'audit/provenance.json',dict(recovered_utc=datetime.now(timezone.utc).isoformat(),
        source_root=str(src),platform=platform.platform(),python=sys.version,files=records,
        note='Content and mtime preserved; source atime may change during reading. Zone.Identifier sidecar excluded; no scientific files excluded.'))
    dump(out/('audit/supplemental_manifests.json' if a.supplement else 'audit/script_variants.json'),scripts)
    (out/'SHA256SUMS-originals.txt').write_text(''.join(f"{r['sha256']}  {r['archive_path']}\n" for r in records))
    repos=[]
    for root in [Path('/home/rahnn/projects/viral-gt-research-record'),Path('/home/rahnn/projects/viral-gt-internal-analyses')]:
        if not root.exists(): continue
        repo=dict(path=str(root))
        for key,cmd in [('remote',['git','remote','-v']),('status',['git','status','--short']),('head',['git','rev-parse','HEAD'])]:
            r=subprocess.run(cmd,cwd=root,text=True,capture_output=True); repo[key]=r.stdout
        for n in ['AGENTS.md','README.md','CITATION.cff','LICENSE']:
            p=root/n
            if p.exists():
                dest=out/'audit/repository_context'/root.name/n; dest.parent.mkdir(parents=True,exist_ok=True); shutil.copy2(p,dest)
        repo['workflows']=[str(p.relative_to(root)) for p in (root/'workflows').glob('*')] if (root/'workflows').exists() else []
        repos.append(repo)
    dump(out/'audit/repository_discovery.json',repos)
    print(json.dumps(dict(files=len(records),bytes=sum(r['bytes'] for r in records),scripts=scripts,repositories=repos),indent=2))
    for name in ['numpy','pandas','matplotlib','PIL','fitz']:
        try:
            mod=__import__(name); print(name,getattr(mod,'__version__','available'))
        except ImportError: print(name,'UNAVAILABLE')

if __name__=='__main__': main()
