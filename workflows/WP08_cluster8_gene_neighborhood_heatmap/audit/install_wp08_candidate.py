"""Install only a validated new workflow; never commit, stage or push."""
import hashlib, json, os, shutil, subprocess, zipfile
from datetime import datetime, timezone
from pathlib import Path

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    base=Path('/mnt/d/Projects/Viral GT cluster 8')
    release=base/'WP08_release_candidate'
    src=release/'workflows/WP08_cluster8_gene_neighborhood_heatmap'
    repo=Path('/home/rahnn/projects/viral-gt-research-record')
    dest=repo/'workflows/WP08_cluster8_gene_neighborhood_heatmap'
    v=json.loads((base/'WP08_package_final_validation/validation.json').read_text())
    assert v['matrix_equal'] and v['evidence_multiset_equal'] and v['originals_preserved']
    assert v['plots'][0]['byte_equal'] and v['plots'][0]['pixel_equal'] and v['pdf']['equal_except_creation_date']
    assert not dest.exists(), 'Refuse to overwrite an existing workflow'
    assert sha(repo/'AGENTS.md')==sha(src/'audit/repository_context/viral-gt-research-record/AGENTS.md'), 'Repository instructions changed'
    for p in [repo.parent/'AGENTS.md',repo.parent.parent/'AGENTS.md',repo/'workflows/AGENTS.md']:
        if p.exists(): raise RuntimeError('Additional instructions need review: '+str(p)+'\n'+p.read_text())
    # Clarify that repository non-modification statements describe candidate assembly.
    updates={
      'README.md': [('The existing repository was inspected read-only and remains unchanged.', 'At candidate assembly the existing repository was inspected read-only and remained unchanged; local installation may subsequently add this validated directory, with a separate installation receipt.')],
      'docs/PROVENANCE.md': [('The existing checkout was not edited;', 'During recovery and candidate assembly the existing checkout was not edited;')],
      'release/PUBLICATION_CHECKLIST.md': [('The candidate was intentionally prepared outside the existing checkout, as requested.', 'The candidate was intentionally prepared and validated outside the existing checkout, as requested. A later local installation may add only this validated directory; that operation is recorded separately and does not publish it.')]
    }
    for n,pairs in updates.items():
        p=src/n; text=p.read_text()
        for before,after in pairs:
            assert before in text; text=text.replace(before,after)
        p.write_text(text)
    shutil.copy2(base/'WP08_package_final_validation/validation.json',src/'validation/assembled_package_validation.json')
    shutil.copy2(__file__,src/'audit/install_wp08_candidate.py')
    checks=src/'SHA256SUMS-package.txt'
    checks.write_text(''.join(f'{sha(p)}  {p.relative_to(src).as_posix()}\n' for p in sorted(src.rglob('*')) if p.is_file() and p!=checks))
    archive=release/'WP08_cluster8_gene_neighborhood_heatmap_v1.0.0_candidate.zip'
    with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in sorted(src.rglob('*')):
            if p.is_file(): z.write(p,p.relative_to(release).as_posix())
    archive.with_suffix('.zip.sha256').write_text(f'{sha(archive)}  {archive.name}\n')
    env=dict(os.environ,GIT_OPTIONAL_LOCKS='0')
    def git(*args): return subprocess.run(['git',*args],cwd=repo,env=env,capture_output=True,text=True,check=True).stdout
    before=git('status','--porcelain=v1','-z')
    changed=git('ls-files','-m','-z').split('\0')
    existing={n:sha(repo/n) for n in changed if n and (repo/n).is_file()}
    head=git('rev-parse','HEAD').strip()
    shutil.copytree(src,dest,copy_function=shutil.copy2)
    files=sorted(p for p in src.rglob('*') if p.is_file())
    assert all(sha(p)==sha(dest/p.relative_to(src)) for p in files)
    assert all(sha(repo/n)==h for n,h in existing.items())
    assert head==git('rev-parse','HEAD').strip()
    after=git('status','--porcelain=v1','-z')
    before_entries=set(before.split('\0'))-{''}; after_entries=set(after.split('\0'))-{''}
    assert before_entries<=after_entries
    additions=after_entries-before_entries
    assert additions=={'?? workflows/WP08_cluster8_gene_neighborhood_heatmap/'}
    with zipfile.ZipFile(archive) as z:
        assert all(hashlib.sha256(z.read(p.relative_to(release).as_posix())).hexdigest()==sha(p) for p in files)
    receipt=dict(timestamp_utc=datetime.now(timezone.utc).isoformat(),repository=str(repo),destination=str(dest),
       local_head=head,files=len(files),all_installed_bytes_match_candidate=True,
       existing_modified_file_hashes_preserved=existing,git_status_additions=sorted(additions),
       candidate_sha256=sha(archive),candidate_bytes=archive.stat().st_size,
       validation='Strict assembled-package run passed before local repository installation; subsequent package changes only clarify documentation and include validation/installation utilities',
       commits_created=False,index_staging_performed=False,pushed=False,zenodo_changed=False)
    (release/'LOCAL_INSTALLATION_RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))

if __name__=='__main__': main()
