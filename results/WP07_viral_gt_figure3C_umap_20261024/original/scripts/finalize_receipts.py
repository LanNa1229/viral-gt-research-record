"""Record final static-export QA and refresh the local artifact inventory."""
import csv,hashlib,json
from pathlib import Path
OUT=Path(__file__).resolve().parents[1]
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
validation=json.loads((OUT/'independent_validation.json').read_text())
assert validation['status']=='PASS' and all(validation['checks'].values())
with (OUT/'output_checksums.csv').open(encoding='utf-8',newline='') as f: old={r['path'].replace('\\','/'):r['sha256'] for r in csv.DictReader(f)}
for name in ('joint_umap_samples.csv','tables/joint_umap_coordinates.npy','tables/representative_15gene_embeddings.npy'):
    assert sha(OUT/name)==old[name], 'Saved coordinates or representative inputs changed'
provenance=json.loads((OUT/'execution_provenance.json').read_text())
provenance.setdefault('fit_script_sha256',provenance['script_sha256'])
provenance['script_sha256']=sha(OUT/'scripts/exploratory_umap.py')
provenance['effective_optimization_epochs']=500
provenance['effective_epoch_source']='Installed umap/umap_.py: default 500 epochs for <=10,000 samples, n_epochs=None, densmap=False'
provenance['postfit_rendering']={'script_sha256':sha(OUT/'scripts/render_saved_coordinates.py'),'reason':'Set Matplotlib colorbar solids rasterized=False for fully vector SVG. Same saved coordinates; no additional fit.','coordinates_verified_unchanged':True}
(OUT/'execution_provenance.json').write_text(json.dumps(provenance,indent=2)+'\n',encoding='utf-8')
qa=dict(status='PASS',inspected=['figures/preview_ABC.png','figures/Figure_A_reference_classes.png','figures/Figure_B_candidate_status.png','figures/Figure_C_saccharide_P10.png'],
    checks=['Readable legends and labels','No clipping or annotation overlap','Same coordinate system and axis limits','Viral triangles visible above references','Seven existing classes retained, including split-color hybrid markers','Figure B uses original candidate labels','Figure C color bar spans 0 to 1','PNG 5100x3420 at 600 dpi','SVG contains no raster image elements','Saved coordinates unchanged after SVG correction'],
    known_limitations=['Dense/exact-duplicate loci overplot','Reference class color is needed to distinguish seven annotations','No color-vision-deficiency simulations performed'],
    no_refit_for_export_correction=True)
(OUT/'visual_qa.json').write_text(json.dumps(qa,indent=2)+'\n',encoding='utf-8')
rows=[dict(path=p.relative_to(OUT).as_posix(),sha256=sha(p)) for p in sorted(OUT.rglob('*')) if p.is_file() and 'cache' not in p.relative_to(OUT).parts and p.name!='output_checksums.csv']
with (OUT/'output_checksums.csv').open('w',encoding='utf-8',newline='') as f:
    w=csv.DictWriter(f,['path','sha256'],lineterminator='\n'); w.writeheader(); w.writerows(rows)
assert all(sha(OUT/r['path'])==r['sha256'] for r in rows)
print(json.dumps(dict(status='PASS',artifacts_in_inventory=len(rows),coordinate_files_unchanged=True,fit_count=provenance['fit_count'])))
