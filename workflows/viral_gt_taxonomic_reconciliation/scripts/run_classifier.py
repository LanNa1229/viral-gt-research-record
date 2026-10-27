import sys,pathlib,subprocess,os,json,datetime,resource,shutil,time
OUT=pathlib.Path(__file__).resolve().parents[1]; BIN=OUT/'env/taxonomy/bin'
task=sys.argv[1];mode=sys.argv[2] if len(sys.argv)>2 else 'full'
if task=='genomad':
    cmd=[str(BIN/'genomad'),'annotate','--threads','2','--splits','32','--full-ictv-lineage',str(OUT/'inputs/unresolved_genomes.fna'),str(OUT/'runs/genomad'),str(OUT/'db/genomad_db')];name='genomad_annotate'
else:
    assert task in ['prokaryotes','eukaryotes'] and mode in ['full','no_gt']
    key=task+('_crosscheck' if len(sys.argv)>3 and sys.argv[3]=='crosscheck' else '')+'_'+mode;name='vcontact3_'+key
    cmd=[str(BIN/'vcontact3'),'run','--proteins',str(OUT/'inputs'/(key+'.faa')),'--gene2genome',str(OUT/'inputs'/(key+'_gene2genome.parquet')),'--len-nucleotide',str(OUT/'inputs'/(key+'_lengths.parquet')),'--db-domain',task,'--db-version','230','--db-path',str(OUT/'db/vcontact3'),'--output',str(OUT/'runs'/name),'--threads','4','--reduce-memory','--keep-temp','--no-progress','--mmseqs-bin',str(BIN/'mmseqs'),'--exports','completeness']
assert shutil.disk_usage('/mnt/c').free>8*1024**3,'C: headroom below 8 GiB; refusing disk-intensive run'
env=os.environ.copy();env.update(PATH=str(BIN)+os.pathsep+env.get('PATH',''),PYTHONDONTWRITEBYTECODE='1',CUDA_VISIBLE_DEVICES='-1',OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2',MPLCONFIGDIR=str(OUT/'cache/matplotlib'),XDG_CACHE_HOME=str(OUT/'cache'))
if task!='genomad':env['PYTHONPATH']=os.pathsep.join([str(OUT/'compat/vendor'),str(OUT/'compat/pandas2'),env.get('PYTHONPATH','')])
record=dict(task=name,command=cmd,start_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),status='running')
status=OUT/'logs'/(name+'.json')
if status.exists():
    stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S');shutil.copyfile(status,OUT/'logs'/(name+'_'+stamp+'_previous.json'))
    old=OUT/'logs'/(name+'.log')
    if old.exists():shutil.copyfile(old,OUT/'logs'/(name+'_'+stamp+'_previous.log'))
status.write_text(json.dumps(record,indent=2))
with (OUT/'logs'/(name+'.log')).open('w') as log:
    process=subprocess.Popen(cmd,cwd=OUT,env=env,stdout=log,stderr=subprocess.STDOUT)
    while process.poll() is None:
        time.sleep(25)
        free=shutil.disk_usage('/mnt/c').free/1024**3
        lines=(OUT/'logs'/(name+'.log')).read_text(errors='replace').splitlines()
        print(name,'C_free_GiB',round(free,2),'last:',lines[-1] if lines else 'starting',flush=True)
    rc=process.returncode
record.update(status='complete' if rc==0 else 'failed',returncode=rc,end_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),child_max_rss_kib=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
status.write_text(json.dumps(record,indent=2));print(json.dumps(record));sys.exit(rc)
