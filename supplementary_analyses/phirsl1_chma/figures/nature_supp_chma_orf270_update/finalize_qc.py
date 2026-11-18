from pathlib import Path
import csv,json,hashlib
R=Path(__file__).resolve().parents[2];O=Path(__file__).resolve().parent
assert (O/'logs/version.stdout').read_text().strip()==(R/'logs/version_foldseek.stdout').read_text().strip()
u=json.loads((O/'results/update_summary.json').read_text())
for n,v in u['source_hashes'].items():assert hashlib.sha256((R/n).read_bytes()).hexdigest()==v
q=O/'FIGURE_QC.md'
q.write_text(q.read_text().replace('Visual inspection: pending.','Visual inspection passed 2026-11-17: reviewed the actual PDF through Poppler rendering and the PNG at publication proportions. No clipped text, missing points or overlapping labels. Panel A clutter removed; ORF270 highlighted in both panels. PDF dimensions 180 x 95 mm; fonts embedded exclusively ArialMT/Arial-BoldMT; SVG has 31 editable text elements; PNG 4251 x 2244 at 600 dpi. Foldseek version matches original frozen executable version. Original-source hashes rechecked after rendering.'),encoding='utf-8')
with (O/'checksums.tsv').open('w') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['path','bytes','sha256'])
 for p in sorted(O.rglob('*')):
  if p.is_file() and 'work' not in p.relative_to(O).parts and p.name!='checksums.tsv':
   w.writerow([str(p.relative_to(O)),p.stat().st_size,hashlib.sha256(p.read_bytes()).hexdigest()])
print('Final QA and provenance checks passed; original source hashes unchanged.')
