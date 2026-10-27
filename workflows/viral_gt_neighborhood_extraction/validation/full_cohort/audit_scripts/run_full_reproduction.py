"""Unpack only the sealed handoff archives into a new /tmp tree and run the unmodified CLI."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import csv,json,hashlib,zipfile,tempfile,subprocess,datetime,shutil,os
D=Path(__file__).resolve().parents[1];ROOT=D.parents[1];RELEASE=ROOT/'analysis/viral_gt_step1_release_20261026'
def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):
    with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def main():
    start=datetime.datetime.now(datetime.timezone.utc).isoformat();receipt=json.loads((RELEASE/'review/archive_receipt.json').read_text());archives=[(RELEASE/receipt['github']['handoff_zip'],receipt['github']['handoff_zip_sha256'],'code'),(RELEASE/receipt['zenodo']['archive'],receipt['zenodo']['archive_sha256'],'dataset')]
    temp=Path(tempfile.mkdtemp(prefix='viral_gt_step1_fullcohort_20261026_',dir='/tmp'));(D/'temporary_workspace.txt').write_text(str(temp)+'\n')
    archive_records=[]
    for p,expected,label in archives:
        assert sha(p)==expected,p
        dest=temp/label;dest.mkdir()
        with zipfile.ZipFile(p) as z:
            assert z.testzip() is None
            assert all(not Path(n).is_absolute() and '..' not in Path(n).parts for n in z.namelist())
            z.extractall(dest)
        archive_records.append(dict(path=str(p),sha256=expected,unpacked_to=str(dest)))
    code=temp/'code/workflows/viral_gt_neighborhood_extraction';dataset=temp/'dataset';data=dataset/'step1';out=temp/'regenerated';input_files={}
    for base in [code,dataset]:
        for r in read(base/'SHA256SUMS.tsv'):
            p=base/r['path'];assert sha(p)==r['sha256'];input_files[str(p)]=r['sha256']
        input_files[str(base/'SHA256SUMS.tsv')]=sha(base/'SHA256SUMS.tsv')
    shutil.copyfile(D/'scripts/guarded_execution.py',temp/'guarded_execution.py')
    py=ROOT/'analysis/bgcprophet_sequence_recovery_20260912/.venv/bin/python';env=dict(os.environ);env['PYTHONPATH']='';env['PYTHONNOUSERSITE']='1';env['PYTHONDONTWRITEBYTECODE']='1'
    command=[str(py),'-I','-B',str(temp/'guarded_execution.py'),'--root',str(temp),'--code',str(code),'--data',str(data),'--output',str(out),'--canary',str(ROOT/'AGENTS.md')]
    setup=dict(started_utc=start,temporary_workspace=str(temp),archives=archive_records,command=command,cwd=str(temp),environment_overrides={k:env[k] for k in ['PYTHONPATH','PYTHONNOUSERSITE','PYTHONDONTWRITEBYTECODE']},runtime='documented existing Python3.12.11 / Biopython1.85 environment; no packages installed',code_and_dataset_manifest_files_verified=len(input_files),fresh_empty_output=True,full_cohort_expected=594)
    (D/'qc/execution_setup.json').write_text(json.dumps(setup,indent=2)+'\n')
    (D/'qc/unpacked_inputs_before.json').write_text(json.dumps(input_files,indent=2)+'\n')
    with (D/'logs/full_extraction.stdout.log').open('w') as stdout,(D/'logs/full_extraction.stderr.log').open('w') as stderr:
        p=subprocess.run(command,cwd=temp,env=env,stdout=stdout,stderr=stderr)
    result=dict(completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),returncode=p.returncode,temporary_workspace=str(temp),regenerated_output=str(out),guard_report_exists=(temp/'guard_summary.json').exists())
    if (temp/'guard_summary.json').exists():shutil.copyfile(temp/'guard_summary.json',D/'qc/filesystem_guard_summary.json')
    if (temp/'filesystem_audit.jsonl').exists():shutil.copyfile(temp/'filesystem_audit.jsonl',D/'logs/filesystem_audit.jsonl')
    if (out/'portable_execution.json').exists():shutil.copyfile(out/'portable_execution.json',D/'qc/portable_execution.json')
    for path,digest in input_files.items():assert sha(path)==digest,('Package input modified',path)
    for path,digest,label in [(p,e,l) for p,e,l in archives]:assert sha(path)==digest
    result['unpacked_code_and_data_unchanged']=True;result['release_archives_unchanged']=True
    (D/'qc/execution_result.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
    if p.returncode:print((D/'logs/full_extraction.stderr.log').read_text()[-6000:]);sys.exit(p.returncode)
    print((D/'logs/full_extraction.stdout.log').read_text()[-3500:])
if __name__=='__main__':main()
