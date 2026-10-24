"""Independent arithmetic and export QA for the initial UMAP; never refits the map."""
import csv, hashlib, importlib.util, json, sys, xml.etree.ElementTree as ET
from pathlib import Path
OUT=Path(__file__).resolve().parents[1]; ROOT=OUT.parents[1]
if sys.platform=='win32': sys.path.insert(0,str(ROOT/'analysis/bgcprophet_reference/.deps'))
import numpy as np
from scipy.spatial.distance import cdist
from PIL import Image
BASE=ROOT/'analysis/viral_gt_mining_20261017/runs/baseline_594_v2'; REF=ROOT/'analysis/viral_glycan_benchmark_20261016'
def read(p,sep=','):
    with p.open(encoding='utf-8-sig',newline='') as f: return list(csv.DictReader(f,delimiter=sep))
def write(name,rows):
    with (OUT/name).open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,list(rows[0]),lineterminator='\n'); w.writeheader(); w.writerows(rows)
checks={}
def require(ok,name):
    if not ok: raise ValueError(name)
    checks[name]=True
ref=np.load(REF/'mibig_CORE_ONLY_embeddings.npy'); viral=np.load(BASE/'viral_bgc_embeddings.npy')
rm=read(REF/'mibig_matched_reference_metadata.tsv','\t'); vm=read(BASE/'embedding_manifest.tsv','\t')
meta={r['window_id']:r for r in vm}; selected=read(OUT/'tables/included_locus_representatives.csv'); registry=read(OUT/'tables/locus_registry.csv')
reps=np.load(OUT/'tables/representative_15gene_embeddings.npy'); xy=np.load(OUT/'tables/joint_umap_coordinates.npy'); samples=read(OUT/'joint_umap_samples.csv')
require(len(registry)==565 and len(selected)==len(reps)==556,'565_baseline_556_included')
require(len(samples)==len(xy)==len(ref)+len(reps),'sample_coordinate_dimensions')
require(len({r['sample_id'] for r in samples})==len(samples),'unique_sample_ids')
require(np.array_equal(np.array([[float(r['UMAP1']),float(r['UMAP2'])] for r in samples],dtype=xy.dtype),xy),'coordinate_CSV_roundtrip')
require([r['sample_id'] for r in samples[:len(ref)]]==[r['BGC'] for r in rm],'reference_order')
require([r['sample_id'] for r in samples[len(ref):]]==[r['viral_locus_id'] for r in selected],'viral_order')
require(np.array_equal(reps,np.stack([viral[int(meta[r['representative_window_id']]['embedding_row'])] for r in selected])),'observed_vectors_only')
medoid_ok=True
for r in registry:
    ids=sorted(w for w in r['all_window_ids'].split(';') if w in meta)
    if not ids:
        medoid_ok &= r['included']=='False' and not r['representative_window_id']; continue
    vectors=viral[[int(meta[w]['embedding_row']) for w in ids]].astype(float)
    if len(ids)==1: best=ids[0]
    else:
        sim=1-cdist(vectors,vectors,metric='cosine'); np.fill_diagonal(sim,0); means=sim.sum(axis=1)/(len(ids)-1)
        best=ids[np.flatnonzero(np.isclose(means,means.max(),atol=1e-12,rtol=0))[0]]
    medoid_ok &= best==r['representative_window_id']
require(medoid_ok,'independent_scipy_cosine_medoids')
cos=1-cdist(reps.astype(float),ref.astype(float),metric='cosine'); ids=np.array([r['BGC'] for r in rm]); pos=np.array(['Saccharide' in r['label'].split() for r in rm])
p10=[]; top=[]
for i,r in enumerate(selected):
    order=np.lexsort((ids,-cos[i])); p10.append(pos[order[:10]].mean()); top.append(ids[order[0]])
require(np.array_equal(p10,[float(r['P_at_10']) for r in selected]),'independent_P10')
require(top==[r['nearest_reference'] for r in selected],'independent_nearest_reference')
require(all(float(s['P_at_10'])==float(r['P_at_10']) and s['candidate_status']==r['candidate_status'] and s['representative_window_id']==r['representative_window_id'] for s,r in zip(samples[len(ref):],selected)),'plotted_scores_labels_representatives')
# Exact duplicate vectors are a QC equivalence relation, not biological clusters.
unique,inverse,counts=np.unique(reps,axis=0,return_inverse=True,return_counts=True); duplicate_rows=[]
for g in np.flatnonzero(counts>1):
    ix=np.flatnonzero(inverse==g); contexts={meta[selected[i]['representative_window_id']]['ordered_sequence_sha256'] for i in ix}
    for i in ix: duplicate_rows.append(dict(exact_vector_group=int(g),n_loci=int(counts[g]),n_unique_ordered_contexts=len(contexts),viral_locus_id=selected[i]['viral_locus_id'],representative_window_id=selected[i]['representative_window_id'],n_genes=selected[i]['representative_n_genes']))
write('tables/exact_duplicate_vector_QC.csv',duplicate_rows)
exports=[]
for name in ('Figure_A_reference_classes','Figure_B_candidate_status','Figure_C_saccharide_P10'):
    png=OUT/'figures'/(name+'.png'); svg=OUT/'figures'/(name+'.svg')
    with Image.open(png) as im:
        require(im.size==(5100,3420),name+'_PNG_resolution'); resolution=list(im.size); dpi=im.info.get('dpi')
    root=ET.parse(svg).getroot(); require(root.tag.endswith('svg'),name+'_SVG_valid_XML')
    require(not root.findall('.//{http://www.w3.org/2000/svg}image'),name+'_SVG_vector_marks')
    exports.append(dict(name=name,png_pixels=resolution,png_dpi=dpi,svg_valid=True,svg_embedded_raster_images=0))
summary=dict(status='PASS',checks=checks,reference_BGCs=len(ref),viral_loci=len(reps),total_points=len(samples),
    unique_representative_vectors=len(unique),exact_duplicate_vector_sets=int(np.count_nonzero(counts>1)),loci_in_duplicate_sets=int(counts[counts>1].sum()),
    duplicate_sets_with_different_ordered_contexts=len({r['exact_vector_group'] for r in duplicate_rows if r['n_unique_ordered_contexts']>1}),exports=exports)
(OUT/'independent_validation.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
print(json.dumps(summary,indent=2))
