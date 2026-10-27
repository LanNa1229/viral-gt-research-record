import sys
sys.dont_write_bytecode=True
import pathlib,hashlib,tarfile,json,datetime,shutil
OUT=pathlib.Path(__file__).resolve().parents[1];p=OUT/'db/vcontact3/v230.tar.gz'
assert shutil.disk_usage('/mnt/c').free>15*1024**3,'Insufficient C: headroom before extraction'
md5=hashlib.file_digest(p.open('rb'),'md5').hexdigest()
assert md5=='dbf915c9da182e0c3dbec60d2156f8b9',md5
print('Verified official Zenodo MD5',md5,flush=True)
with tarfile.open(p,'r:gz') as archive:
    selected=[];excluded=[]
    for m in archive.getmembers():
        parts=pathlib.PurePosixPath(m.name).parts
        reason='publisher_scratch_or_Mac_metadata' if '_tmp' in m.name or any(x.startswith('._') for x in parts) else 'absolute_link' if m.issym() and pathlib.PurePosixPath(m.linkname).is_absolute() else ''
        if reason:excluded.append(dict(path=m.name,reason=reason))
        else:selected.append(m)
    print('Selected unpacked bytes:',sum(m.size for m in selected),flush=True)
    archive.extractall(p.parent,members=selected,filter='data')
    (OUT/'reports/reference_archive_excluded_entries.json').write_text(json.dumps(excluded,indent=2))
(OUT/'reports/vcontact_reference_provenance.json').write_text(json.dumps(dict(version='230',ncbi_refseq_release='230',genbank_release='266',source_url='https://zenodo.org/records/15598380/files/v230.tar.gz',source_doi='10.5281/zenodo.15598380',archive_md5=md5,verified_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()),indent=2))
print('Reference extracted',flush=True)
