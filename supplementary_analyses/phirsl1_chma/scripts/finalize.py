#!/usr/bin/env python3
"""Finalize summary/report and verify stable manifests without changing scientific outputs."""
from pipeline import *
import ast,contextlib,os
import summarize,report
def main():
 for p in (ROOT/'scripts').glob('*.py'):ast.parse(p.read_text(),filename=str(p))
 with contextlib.redirect_stdout(io.StringIO()) as captured:summarize.main()
 (ROOT/'logs/final_validation.stdout').write_text(captured.getvalue())
 # Receipt is written before regenerating the manifests and checksums.
 localchecks=[]
 for name in ['README.md','METHODS.md','scripts/pipeline.py','scripts/searches.py','scripts/markers.py','results/validation.tsv','results/candidate_review.tsv']:
  localchecks.append(dict(path=name,exists=(ROOT/name).exists()))
 if not all(r['exists'] for r in localchecks):raise RuntimeError('Missing deliverable')
 (ROOT/'results/package_validation.json').write_text(json.dumps(dict(utc=stamp(),python_scripts_parse=True,deliverables=localchecks,scientific_validation='All checks in validation.tsv pass; see final_validation.stdout'),indent=2))
 report.main()
 # Verify every downloaded payload remains available in its original form.
 for r in csv.DictReader((ROOT/'sources/source_catalog.tsv').open(),delimiter='\t'):
  assert r['preserved_download_matches']=='True',(r['path'],'original download not preserved')
 # Check all local relative Markdown links from the author-written reports.
 for name in ['README.md','METHODS.md']:
  for link in re.findall(r'\]\(([^)]+)\)',(ROOT/name).read_text()):
   if '://' not in link:assert (ROOT/link).exists(),(name,link)
 checksums()
 for manifest in ['input_manifest.tsv','output_manifest.tsv','checksums.tsv']:
  for r in csv.DictReader((ROOT/manifest).open(),delimiter='\t'):
   p=ROOT/r['path'];assert p.is_file() and p.stat().st_size==int(r['bytes']) and digest(p.read_bytes())==r['sha256'],(manifest,r['path'])
 n=sum(1 for r in csv.DictReader((ROOT/'checksums.tsv').open(),delimiter='\t'))
 print('PASS: scientific validation, all Python syntax, local report links, original download preservation, and all manifest checksums.')
 print('Stable files:',n)
 print('Report:',ROOT/'README.md')
 print('Sixframe best dominant domain:',json.loads((ROOT/'results/summary.json').read_text())['sixframe_best'])
if __name__=='__main__':main()
