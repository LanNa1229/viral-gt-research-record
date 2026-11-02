"""Scoped archival publication. No scientific scripts are executed by this utility."""
import argparse,csv,hashlib,json,os,shutil,subprocess,sys,time
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path

BASE=Path('/mnt/d/Projects/Viral GT cluster 8')
RUN=BASE/'WP08_github_publication_20261102'
REPO=Path('/home/rahnn/projects/viral-gt-research-record')
REL=Path('workflows/WP08_cluster8_gene_neighborhood_heatmap')
DEST=REPO/REL
SOURCE=BASE/'WP08_release_candidate'/REL
ENV=dict(os.environ,GIT_OPTIONAL_LOCKS='0')
def digest(p):
 h=hashlib.sha256();g=hashlib.sha1();n=p.stat().st_size;g.update(f'blob {n}\0'.encode())
 with p.open('rb') as f:
  for block in iter(lambda:f.read(4*1024*1024),b''):h.update(block);g.update(block)
 return h.hexdigest(),g.hexdigest()
def sha(p):return digest(p)[0]
def git(*args,check=True,input=None):
 r=subprocess.run(['git',*map(str,args)],cwd=REPO,env=ENV,input=input,capture_output=True)
 if check and r.returncode:raise RuntimeError(f'git {args}: '+r.stderr.decode(errors='replace'))
 return r
def textgit(*a):return git(*a).stdout.decode()
def dump(p,d):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(d,indent=2)+'\n')
def write(rel,text):p=DEST/rel;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(text.strip()+'\n');return p
def snapshot():
 status=textgit('status','--porcelain=v1','-z')
 modified=textgit('ls-files','-m','-z').split('\0')
 return dict(head=textgit('rev-parse','HEAD').strip(),branch=textgit('branch','--show-current').strip(),status=status,
    modified_hashes={p:sha(REPO/p) for p in modified if p and (REPO/p).is_file()},
    index_empty=not bool(textgit('diff','--cached','--name-only')))
def verify_unrelated(before):
 assert all(sha(REPO/p)==h for p,h in before['modified_hashes'].items()),'Unrelated user file changed'
def copy_exact(src,dest):
 h,g=digest(src);st=src.stat();dest.parent.mkdir(parents=True,exist_ok=True)
 if dest.exists():assert sha(dest)==h,'Existing source copy changed: '+str(dest)
 else:shutil.copy2(src,dest)
 assert sha(dest)==h
 assert src.stat().st_mtime_ns==st.st_mtime_ns
 return h,g,st

def prepare():
 RUN.mkdir(exist_ok=True)
 before=snapshot();assert before['index_empty'] and before['branch']=='main'
 remote=textgit('ls-remote','origin','refs/heads/main').split()[0]
 assert remote==before['head'],'Local and remote main differ; reconcile before publication'
 assert not textgit('ls-tree','HEAD',str(REL)).strip(),'WP08 already tracked; inspect before changing'
 dump(RUN/'working_tree_before.json',before)
 # The full local candidate is preserved separately. Publish only historical science.
 excluded=[];selected=set()
 for p in SOURCE.rglob('*'):
  if not p.is_file():continue
  rel=p.relative_to(SOURCE)
  if rel.parts[0] in ['corrected_interpretation','release'] or rel.as_posix()=='correct_interpretation.py':
   excluded.append(dict(path=rel.as_posix(),reason='Proposed correction or superseded release/Zenodo draft; preserved locally, outside this historical-only publication'));continue
  if rel.as_posix()=='SHA256SUMS-package.txt':continue
  if rel.as_posix()=='README.md':continue
  copy_exact(p,DEST/rel);selected.add(rel.as_posix())
 prov=json.loads((SOURCE/'audit/provenance.json').read_text())
 assert len(prov['files'])==346
 for r in prov['files']:assert sha(DEST/r['archive_path'])==r['sha256']
 print('Original validated files: 346/346 byte-identical',flush=True)
 plan=json.loads((RUN/'source_supplement_plan.json').read_text())
 assert not plan['missing'] and not plan['over_100_MiB']
 assert all(t['matches_historical'] for l in plan['loci'] for t in l['run_raw_tables'])
 supplements=[]
 assert shutil.disk_usage(REPO).free>plan['bytes']*2+1024**3,'Insufficient free disk for preservation and Git objects'
 for i,r in enumerate(plan['files'],1):
  src=Path(r['source']);assert src.is_file() and str(src).startswith('/home/rahnn/2606')
  rel=Path('originals/source_runs')/src.relative_to('/home/rahnn')
  h,g,st=copy_exact(src,DEST/rel);selected.add(rel.as_posix())
  supplements.append(dict(**r,archive_path=rel.as_posix(),sha256=h,git_blob_sha1=g,source_mtime_ns=st.st_mtime_ns))
  if i%100==0:print(f'Preserved source supplement {i}/{len(plan["files"])}',flush=True)
 dump(DEST/'audit/source_supplement_provenance.json',dict(recovered_utc=datetime.now(timezone.utc).isoformat(),files=supplements,loci=plan['loci'],missing=[]))
 selected.add('audit/source_supplement_provenance.json')
 (DEST/'SHA256SUMS-source-supplement.txt').write_text(''.join(f"{r['sha256']}  {r['archive_path']}\n" for r in supplements));selected.add('SHA256SUMS-source-supplement.txt')
 # Include the targeted manuscript audit; exclude contact sheets of other manuscript figures.
 audit=BASE/'WP08_manuscript_comparison_20261102'
 audit_exclusions=[]
 for p in audit.rglob('*'):
  if not p.is_file():continue
  rel=p.relative_to(audit)
  if p.name.startswith('contact_sheet_') or p.name=='SHA256SUMS-audit.txt':
   audit_exclusions.append(dict(path=rel.as_posix(),reason='Contact sheet of other manuscript figures, or superseded full-directory checksum'));continue
  target=Path('validation/manuscript_comparison')/rel
  copy_exact(p,DEST/target);selected.add(target.as_posix())
 # Existing licenses and notices are copied without changing their scoped grants.
 for name in ['LICENSE-MIT.txt','LICENSE-CC-BY-4.0.md','THIRD_PARTY_NOTICES.md']:
  target=Path('licenses/repository_existing')/name;copy_exact(REPO/name,DEST/target);selected.add(target.as_posix())
 def add(rel,text):write(rel,text);selected.add(rel)
 add('README.md','''
# Main Figure 4 — Cluster 8 genomic-neighborhood heatmap

This historical research archive preserves the computational record used to compare the genomic context of 11 viral Cluster 8 glycosyltransferase loci across 33 functional gene classes. The figure summarizes the strongest named structural evidence per locus/class; it is not a count of genes or proof of biochemical activity.

**Lan Na — research conception, methodology, computational analysis, data curation, and visualization.**

All 346 files from the validated October 29 recovery are published individually and unchanged. Additional original query models, saved raw API responses, per-query tables and run records are preserved under `originals/source_runs/`, traced to the exact run-level hit tables used by the heatmap. No Foldseek searches or scientific analyses were rerun for publication. Proposed corrected outputs remain local and are not included in this historical-only version.

## Four-stage computational workflow

1. **Structural search:** `foldseek_batch_api_v2.py` reads query structure manifests, submits structures to the Foldseek API and saves response JSON, run ledgers and hit tables. Original script variants, query models and responses are archived; current online databases cannot recreate historical server state exactly.
2. **Annotation summaries:** `analyze_foldseek_results_v2.py` processes saved `foldseek_all_hits.tsv`, derives readable names and keyword tags, and writes annotation/high-confidence/ORF summary tables including `10_database_top_hit_names.tsv`.
3. **Neighborhood matrix:** `build_cluster8_neighborhood_matrix.py` reads `cluster8_analysis_manifest.fixed.tsv` and the 11 top-hit tables, applies the original 0–3 evidence and class rules, then writes the 11 × 33 matrix and long evidence table. The original focal-ID mismatches and row classifications are retained.
4. **Plotting:** `plot_cluster8_panelB_heatmap.py` and `plot_cluster8_panelB_heatmap_clean.py` render the original and clean heatmaps. No biological curation or figure redesign was applied to this archive.

## Inputs, code and results

- [Original fixed manifest](originals/260630_Temp/cluster8_analysis_manifest.fixed.tsv).
- [Matrix builder](originals/260630_Temp/build_cluster8_neighborhood_matrix.py), [original plotter](originals/260630_Temp/plot_cluster8_panelB_heatmap.py), [clean plotter](originals/260630_Temp/plot_cluster8_panelB_heatmap_clean.py).
- [Batch and annotation script variants / manifests](originals/related_scripts/) and [source query models / complete saved runs](originals/source_runs/).
- `originals/260630_Temp/<candidate_id>/final/`: all 11 original `foldseek_all_hits.tsv` files and `analysis_v2/` tables, including all 11 `10_database_top_hit_names.tsv` files.
- [Figure source matrix](originals/260630_Temp/cluster8_global_neighborhood_panelB/cluster8_neighborhood_matrix.tsv) and [4,919-record evidence table](originals/260630_Temp/cluster8_global_neighborhood_panelB/cluster8_neighborhood_evidence_long.tsv).
- [Original and clean PDF/PNG figures](originals/260630_Temp/cluster8_global_neighborhood_panelB/).
- [Main Figure 4 source mapping](docs/FIGURE_SOURCE_MAP.tsv), [Methods and confidence scoring](docs/METHODS_AND_SCORING.md), [all 33 gene-class definitions](docs/GENE_CLASS_ASSIGNMENTS.md).

## Dependencies and reproduction

Offline reconstruction from the preserved top-hit tables requires Python, pandas, NumPy, Matplotlib and Pillow; see [environment.yml](environment.yml). The tested matching environment used Python 3.12.11, pandas 2.3.3, NumPy 2.4.6, Matplotlib 3.10.8, Pillow 12.0.0 and FreeType 2.6.1. The API script additionally requires `requests`, but it is not invoked by offline reproduction.

From this workflow directory, choose a new output directory:

```bash
python -B reproduce.py --out /absolute/path/to/new_validation
```

For the exact rendering check in the matching Linux environment:

```bash
python -B reproduce.py --out /absolute/path/to/new_validation --require-exact-render
```

To plot the immutable historical matrix directly:

```bash
python -B originals/260630_Temp/plot_cluster8_panelB_heatmap_clean.py \\
  --matrix originals/260630_Temp/cluster8_global_neighborhood_panelB/cluster8_neighborhood_matrix.tsv \\
  --outdir /absolute/path/to/new_figure_output
```

Expected offline validation outputs include rebuilt matrix/evidence TSVs, validation JSON, and the clean heatmap PDF/PNG in `plot_with_title/`. The wrapper resolves packaged input tables without modifying preserved absolute paths or focal accession IDs. [Detailed reproducibility instructions](docs/REPRODUCIBILITY.md) also describe the original builder CLI. The original API stage is preserved for provenance, not automatically rerun.

## Completed validation

- 11 loci, 33 classes, **363/363 score cells reproduced**.
- **4,919/4,919 evidence records reproduced**; evidence row order is compared as a multiset because the original code iterates sets.
- Original clean **PNG reproduced byte-for-byte**; the reproduced PDF differs only in creation timestamp.
- Original-file hashes and source mtimes were checked during recovery; publication verifies every archived byte against its recorded checksum and its Git blob.
- [Manuscript comparison](validation/manuscript_comparison/AUDIT_REPORT.md): **363/363 cells match** in actual v1.4 and final available v1.5; zero score discrepancies or uncertain cells. Their embedded PNGs are byte-identical. Word crops the bottom 14.262%, hiding the group legend but not the score grid; this is a presentation difference, not a scientific correction.

## Limitations, provenance and rights

Retain the [known historical limitations](docs/KNOWN_LIMITATIONS.md), including A0A6J5M498 row-group fall-through, broad nucleotide-metabolism categories, ambiguous ORF assignments, and five manifest-query identity mismatches. No IDs, scores, annotations, classifications or original figures have been corrected. The earlier recovery reports may describe local proposed corrections; those outputs are outside this historical publication and any future curated version must be separate.

[Recovery report](docs/RECOVERY_REPORT.md), [provenance](docs/PROVENANCE.md), [publication scope and source supplement](docs/PUBLICATION_SCOPE.md), and source checksums document the record. `SHA256SUMS-originals.txt` covers the original 346; `SHA256SUMS-source-supplement.txt` covers additional recovered source files; `SHA256SUMS-package.txt` covers the exact historical publication payload except itself. Historical local paths are execution provenance, not portable defaults.

The repository's existing [MIT notice](licenses/repository_existing/LICENSE-MIT.txt), [CC BY notice](licenses/repository_existing/LICENSE-CC-BY-4.0.md) and [third-party attribution document](licenses/repository_existing/THIRD_PARTY_NOTICES.md) are retained verbatim. Their original Section 1 scope is not expanded to WP08 by this archive. [WP08 rights and attribution](THIRD_PARTY_NOTICES.md) distinguishes existing grants from this authorized historical publication. No Zenodo deposition or release is part of this task.
''')
 # Keep previous reports intact, while explicitly identifying their earlier scope.
 add('docs/PUBLICATION_SCOPE.md',f'''
# Historical publication scope — November 2, 2026

This publication is authorized by Lan Na's explicit request to commit and push the complete historical Main Figure 4 record to `LanNa1229/viral-gt-research-record` on main. It introduces only this workflow directory. No repository-root files or unrelated work are included.

All 346 validated recovered originals are present at the same relative paths and match `SHA256SUMS-originals.txt`. The additional source inventory contains {len(supplements)} unchanged original files ({sum(r['bytes'] for r in supplements):,} bytes), located by saved raw-JSON and query-structure references in the historical hit tables. All 11 source-run `foldseek_all_hits.tsv` tables match their corresponding historical copies. No scientific transformation was performed during this curation.

The B2ZYD6 run output was relocated from `260617_Foldseek_batch/260618_foldseek_api_BAG41682_neighbors` to the same run-folder name under `260622_Foldseek_batch`. Its original input structures remain under `260617_Foldseek_batch/structures/`. Both actual source locations are recorded; original path strings are unchanged. The source run for A0A159B7X0 is historically spelled `A0A169B7X0_foldseek`; that spelling is retained, not silently corrected. Additional no-hit query models and saved empty-query responses are retained with their run folders.

`audit/source_supplement_provenance.json` records actual source paths, sizes, mtimes and SHA-256. The largest supplemental file is {max(r['bytes'] for r in supplements):,} bytes, below GitHub's 100 MiB file threshold; ordinary Git blobs preserve all files individually. No LFS pointer or ZIP-only substitute is used.

The completed manuscript audit is included under `validation/manuscript_comparison/`. Contact sheets of unrelated manuscript figures are excluded; the two actual heatmap images, source document hashes/XML, original-figure copies, complete 363-cell comparisons, decoded figure matrix and audit scripts/reports remain. A package-wide checksum covers the selected audit files; its former full-directory checksum remains in the untouched local audit.

Local proposed correction outputs and superseded release/Zenodo drafts are not selected for this historical-only commit. Their exact paths and reasons appear in `audit/publication_scope.json`. These are generated artifacts, not any of the 346 recovered original files. All remain available in the untouched Windows release candidate or as untracked local files. Prior recovery/validation reports describe their own execution-time state and are retained as historical records; the current README and this scope document govern the public payload.

Attribution: Lan Na — research conception, methodology, computational analysis, data curation, and visualization. The figure identification as Main Figure 4 is supplied by the user's publication request. Original code, input data and figures remain byte-preserved. No new biology, annotation correction, figure redesign, Foldseek search, scientific rerun or Zenodo operation is authorized or performed by this publication utility.
''')
 add('THIRD_PARTY_NOTICES.md','''
# WP08 attribution and existing license scope

Lan Na — research conception, methodology, computational analysis, data curation, and visualization. Publication of this complete historical record is explicitly authorized by the user on November 2, 2026.

The repository's existing license and attribution documents are copied byte-for-byte in `licenses/repository_existing/`. Their original MIT and CC BY grants apply to their stated Section 1 scope; this archival publication does not silently extend those grants to unrelated code or data. No additional WP08 redistribution license is invented here. Original copyright and third-party rights remain applicable.

The historical inputs and exports retain the Foldseek target identifiers, database names and query model provenance. Relevant original sources include Foldseek, PDB, AlphaFold Database, BFVD, UniProt and CATH. The existing attribution document provides their original source citations and notices. Preserve those credits and database identifiers when using this record; a public copy does not change third-party terms.

Saved query models and API results are archived as original research inputs/outputs, without claiming new ownership of database material. No Foldseek executable or live reference database is bundled. Historical source database builds and server state are not reconstructed. Reproduction of the historical heatmap uses the included saved results and requires no live API query.
''')
 mapping=[
  ('Main Figure 4','historical matrix','originals/260630_Temp/cluster8_global_neighborhood_panelB/cluster8_neighborhood_matrix.tsv','11 loci x 33 functional classes; unchanged'),
  ('Main Figure 4','historical evidence','originals/260630_Temp/cluster8_global_neighborhood_panelB/cluster8_neighborhood_evidence_long.tsv','4919 original evidence records; unchanged'),
  ('Main Figure 4','historical clean PNG','originals/260630_Temp/cluster8_global_neighborhood_panelB/figures_clean/PanelB_cluster8_global_neighborhood_heatmap_clean.png','Original 4144 x 1958 PNG; unchanged'),
  ('Main Figure 4','historical clean PDF','originals/260630_Temp/cluster8_global_neighborhood_panelB/figures_clean/PanelB_cluster8_global_neighborhood_heatmap_clean.pdf','Original PDF; unchanged'),
  ('Main Figure 4','clean plotting script','originals/260630_Temp/plot_cluster8_panelB_heatmap_clean.py','Recovered original script'),
  ('Main Figure 4','fixed locus manifest','originals/260630_Temp/cluster8_analysis_manifest.fixed.tsv','Retains five focal query mismatches'),
  ('Main Figure 4','manuscript audit','validation/manuscript_comparison/AUDIT_REPORT.md','v1.4 and final v1.5; all 363 cells match'),
  ('Main Figure 4','cell-by-cell source mapping','validation/manuscript_comparison/all_363_score_cell_comparisons.csv','One record per locus/class cell'),
  ('Main Figure 4','manuscript embedded image','validation/manuscript_comparison/manuscript_images/Viral_GT_draftv1.4_261022_LanNa_references_added/word/media/image19.png','SHA256 8a938331cff92bd986ec142b63c18cf316ad6557cf6f2172701ee9178764e5bb')]
 p=DEST/'docs/FIGURE_SOURCE_MAP.tsv'
 with p.open('w',newline='') as f:
  w=csv.writer(f,delimiter='\t');w.writerow(['figure','role','path','notes']);w.writerows(mapping)
 selected.add('docs/FIGURE_SOURCE_MAP.tsv')
 dump(DEST/'audit/publication_scope.json',dict(date='2026-11-02',figure='Main Figure 4 — Cluster 8 Genomic Neighborhood Heatmap',original_files_preserved=346,supplemental_original_files=len(supplements),supplemental_bytes=sum(r['bytes'] for r in supplements),local_generated_artifacts_not_selected=excluded,manuscript_audit_files_not_selected=audit_exclusions,original_science_modified=False,scientific_analysis_executed=False,zenodo_operation_performed=False))
 selected.add('audit/publication_scope.json')
 for name in ['inspect_wp08_publication.py','plan_wp08_source_supplement.py','publish_wp08_github.py']:
  rel=Path('audit/publication_tools')/name;copy_exact(BASE/name,DEST/rel);selected.add(rel.as_posix())
 for name in ['git_history_before.json','raw_input_references.json','related_directory_inventory.json','source_supplement_plan.json']:
  rel=Path('audit/publication_preflight')/name;copy_exact(RUN/name,DEST/rel);selected.add(rel.as_posix())
 # Generated documents replacing candidate-only files are included; preserved candidate files stay unchanged.
 manifest=DEST/'SHA256SUMS-package.txt'
 manifest.write_text(''.join(f'{sha(DEST/p)}  {p}\n' for p in sorted(selected)))
 selected.add('SHA256SUMS-package.txt')
 inventory=[]
 for p in sorted(selected):
  path=DEST/p;h,g=digest(path)
  assert path.stat().st_size<100*1024*1024,'GitHub oversized file: '+p
  inventory.append(dict(path=(REL/p).as_posix(),workflow_path=p,bytes=path.stat().st_size,sha256=h,git_blob_sha1=g))
 dump(RUN/'publication_inventory.json',inventory)
 (RUN/'publish_pathspec.nul').write_bytes(b''.join(r['path'].encode()+b'\0' for r in inventory))
 verify_unrelated(before)
 print(json.dumps(dict(selected_files=len(inventory),original_files=346,supplemental_files=len(supplements),total_bytes=sum(r['bytes'] for r in inventory),largest_bytes=max(r['bytes'] for r in inventory),scope=str(REL)),indent=2),flush=True)

def commit():
 before=json.loads((RUN/'working_tree_before.json').read_text());verify_unrelated(before)
 expected={r['path'] for r in json.loads((RUN/'publication_inventory.json').read_text())}
 current_staged=set(textgit('diff','--cached','--name-only','-z').split('\0'))-{''}
 assert current_staged<=expected,'Index contains unrelated staged changes'
 assert textgit('rev-parse','HEAD').strip()==before['head']
 assert textgit('ls-remote','origin','refs/heads/main').split()[0]==before['head'],'Remote changed since preparation'
 inventory=json.loads((RUN/'publication_inventory.json').read_text())
 for r in inventory:assert sha(REPO/r['path'])==r['sha256']
 print('All publication files verified before staging',flush=True)
 git('add','-f','--pathspec-from-file='+str(RUN/'publish_pathspec.nul'),'--pathspec-file-nul')
 staged=textgit('diff','--cached','--name-only','-z').split('\0');staged=set(staged)-{''}
 assert staged=={r['path'] for r in inventory}
 assert all(p.startswith(REL.as_posix()+'/') for p in staged)
 index={}
 for line in textgit('ls-files','--stage','-z',str(REL)).split('\0'):
  if line:
   left,path=line.split('\t',1);index[path]=left.split()[1]
 assert all(index[r['path']]==r['git_blob_sha1'] for r in inventory),'Staged Git blob differs from source bytes'
 status=textgit('diff','--cached','--stat');(RUN/'commit_diff_stat.txt').write_text(status)
 git('commit','-m','Archive original Cluster 8 genomic-neighborhood heatmap analysis and source data')
 new=textgit('rev-parse','HEAD').strip()
 changed=set(textgit('diff-tree','--no-commit-id','--name-only','-r',new).splitlines())
 assert changed==staged
 assert textgit('rev-parse',new+'^').strip()==before['head']
 verify_unrelated(before)
 dump(RUN/'commit_receipt.json',dict(commit_sha=new,parent=before['head'],files=len(inventory),original_files=346,supplemental_original_files=len(json.loads((DEST/'audit/source_supplement_provenance.json').read_text())['files']),only_requested_directory_modified=True,all_staged_blobs_match_source_bytes=True,unrelated_user_changes_preserved=True,pushed=False))
 print(json.dumps(json.loads((RUN/'commit_receipt.json').read_text()),indent=2),flush=True)

def repair_staging():
    """Resolve generic ignore patterns without resetting the index or excluding originals."""
    before=json.loads((RUN/'working_tree_before.json').read_text());verify_unrelated(before)
    inventory=json.loads((RUN/'publication_inventory.json').read_text())
    expected={r['path'] for r in inventory}
    staged=set(textgit('diff','--cached','--name-only','-z').split('\0'))-{''}
    assert staged<=expected and textgit('rev-parse','HEAD').strip()==before['head']
    target='audit/publication_tools/publish_wp08_github.py'
    shutil.copy2(Path(__file__),DEST/target)
    scope=DEST/'docs/PUBLICATION_SCOPE.md'
    scope.write_text(scope.read_text()+'\nThe repository has generic ignore rules for CIF models and log files. This authorized historical publication stages only the explicit checksum-verified WP08 path list with `git add -f --pathspec-from-file`; it does not change repository-wide ignore rules, omit source models/logs, or stage any unrelated path.\n')
    changed={target,'docs/PUBLICATION_SCOPE.md'}
    for r in inventory:
        if r['workflow_path'] in changed:
            r['sha256'],r['git_blob_sha1']=digest(REPO/r['path']);r['bytes']=(REPO/r['path']).stat().st_size
    manifest=DEST/'SHA256SUMS-package.txt'
    manifest.write_text(''.join(f"{r['sha256']}  {r['workflow_path']}\n" for r in inventory if r['workflow_path']!='SHA256SUMS-package.txt'))
    for r in inventory:
        if r['workflow_path']=='SHA256SUMS-package.txt':
            r['sha256'],r['git_blob_sha1']=digest(manifest);r['bytes']=manifest.stat().st_size
    dump(RUN/'publication_inventory.json',inventory)
    print('Scoped staging repaired; original scientific bytes and unrelated index entries are unchanged',flush=True)

def push():
 before=json.loads((RUN/'working_tree_before.json').read_text());receipt=json.loads((RUN/'commit_receipt.json').read_text())
 assert textgit('rev-parse','HEAD').strip()==receipt['commit_sha']
 remote=textgit('ls-remote','origin','refs/heads/main').split()[0]
 assert remote in [receipt['parent'],receipt['commit_sha']],'Remote main changed; no overwrite permitted'
 verify_unrelated(before)
 if remote!=receipt['commit_sha']:
  # Inherit the current user's Git authentication. Never print credential configuration.
  r=subprocess.run(['git','push','origin','HEAD:refs/heads/main'],cwd=REPO,env=ENV)
  if r.returncode:sys.exit(r.returncode)
 remote=textgit('ls-remote','origin','refs/heads/main').split()[0]
 assert remote==receipt['commit_sha'];verify_unrelated(before)
 receipt.update(pushed=True,remote_main_sha=remote,directory_url='https://github.com/LanNa1229/viral-gt-research-record/tree/main/'+REL.as_posix())
 dump(RUN/'push_receipt.json',receipt);print(json.dumps(receipt,indent=2),flush=True)

def verify():
 receipt=json.loads((RUN/'push_receipt.json').read_text());before=json.loads((RUN/'working_tree_before.json').read_text())
 remote=textgit('ls-remote','origin','refs/heads/main').split()[0];assert remote==receipt['commit_sha']
 inventory=json.loads((RUN/'publication_inventory.json').read_text())
 tree={}
 for line in textgit('ls-tree','-r','-z',remote,str(REL)).split('\0'):
  if line:
   left,path=line.split('\t',1);tree[path]=left.split()[2]
 assert set(tree)=={r['path'] for r in inventory}
 assert all(tree[r['path']]==r['git_blob_sha1'] for r in inventory)
 # Verify every committed byte from Git object contents independently of the working files.
 proc=subprocess.Popen(['git','cat-file','--batch'],cwd=REPO,env=ENV,stdin=subprocess.PIPE,stdout=subprocess.PIPE)
 for i,r in enumerate(inventory,1):
  proc.stdin.write((tree[r['path']]+'\n').encode());proc.stdin.flush()
  header=proc.stdout.readline().decode().split();assert header[1]=='blob';remain=int(header[2]);h=hashlib.sha256()
  assert remain==r['bytes']
  while remain:
   block=proc.stdout.read(min(remain,4*1024*1024));assert block;h.update(block);remain-=len(block)
  assert proc.stdout.read(1)==b'\n';assert h.hexdigest()==r['sha256']
  if i%200==0:print(f'Verified committed SHA-256 {i}/{len(inventory)}',flush=True)
 proc.stdin.close();assert proc.wait()==0
 verify_unrelated(before)
 receipt.update(all_committed_file_sha256_verified=True,remote_main_confirmed=True,tracked_originals=346,
  uploaded_file_count=len(inventory),total_source_bytes=sum(r['bytes'] for r in inventory),
  unresolved_upload_problems=[],zenodo_created=False,git_status_after=textgit('status','--short'))
 dump(RUN/'final_verification.json',receipt)
 print(json.dumps(receipt,indent=2),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('phase',choices=['prepare','commit','push','verify','repair-staging']);a=p.parse_args()
    {'prepare':prepare,'commit':commit,'push':push,'verify':verify,'repair-staging':repair_staging}[a.phase]()
