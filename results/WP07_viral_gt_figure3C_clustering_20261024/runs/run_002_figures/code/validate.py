"""Independent alignment, graph/objective arithmetic, membership and export checks."""
import ast,collections,csv,json,math,xml.etree.ElementTree as ET
import common as c
import numpy as np
from scipy.spatial.distance import cdist
from PIL import Image
checks={}
def check(v,name):
 if not v:raise ValueError('VALIDATION: '+name)
 checks[name]=True
old=c.ROOT/c.CONFIG['previous_analysis'];meta=c.read(old/'tables/included_locus_representatives.csv');x=np.load(old/'tables/representative_15gene_embeddings.npy');selected=[i for i,r in enumerate(meta) if r['candidate_status']=='candidate'];source={r['viral_locus_id']:r for r in meta}
rows=c.read(c.OUT/'tables/per_locus_selected287_clustering.csv');allrows=c.read(c.OUT/'tables/per_locus_all556_clustering.csv');ids=[r['viral_locus_id'] for r in rows];labels=np.array([r['Leiden_cluster'] for r in rows]);lookup={sid:i for i,sid in enumerate(ids)}
check(len(rows)==287 and len(allrows)==556 and ids==[meta[i]['viral_locus_id'] for i in selected],'primary_and_all_cohort_alignment')
required=['viral_locus_id','genome_unit_id','representative_window_id','representative_n_genes','candidate_status','status11','status15','P_at_10','Leiden_cluster','HDBSCAN_cluster','HDBSCAN_probability','nearest_reference','nearest_reference_class','nearest_reference_cosine','duplicate_context_group','UMAP1','UMAP2'];check(set(required)<=set(rows[0]),'required_per_locus_fields')
for r in allrows:
 o=source[r['viral_locus_id']];check(all(r[k]==v for k,v in o.items()),'original_metadata_preserved')
um={r['sample_id']:r for r in c.read(old/'joint_umap_samples.csv')}
check(all(float(r['UMAP1'])==float(um[r['viral_locus_id']]['UMAP1']) and float(r['UMAP2'])==float(um[r['viral_locus_id']]['UMAP2']) for r in allrows),'coordinates_exactly_reused')
check(all(0<=float(r['HDBSCAN_probability'])<=1 for r in rows),'HDBSCAN_probability_range')
check(all(float(r['HDBSCAN_probability'])==0 for r in rows if r['HDBSCAN_cluster']=='Noise'),'HDBSCAN_noise_probability_zero')
# SciPy distances independently verify the positive-cosine union kNN baseline graph.
v=x[selected].astype(float);d=cdist(v,v,metric='cosine');keys=[row.tobytes() for row in x[selected]]
for i in range(len(v)):
 for j in range(len(v)):
  if keys[i]==keys[j]:d[i,j]=0
np.fill_diagonal(d,0);d=np.maximum(d,0)
edges=c.read(c.OUT/'graphs/selected287_k20_edges.csv');actual={(min(lookup[r['source']],lookup[r['target']]),max(lookup[r['source']],lookup[r['target']])) for r in edges};expected=set();idarray=np.array(ids)
for i in range(len(v)):
 order=np.lexsort((idarray,d[i]));order=order[order!=i][:20]
 for j in order:
  if 1-d[i,j]>0:expected.add((min(i,int(j)),max(i,int(j))))
# Differences from independent accumulation can only be floating ties; require every
# edge disagreement to be within 1e-12 of an independently tied k boundary.
extra=actual^expected;ambiguous=[]
for i,j in extra:
 legitimate=False
 for a,b in [(i,j),(j,i)]:
  values=np.delete(d[a],a);boundary=np.partition(values,19)[19]
  if abs(d[a,b]-boundary)<=1e-12:legitimate=True
 if legitimate:ambiguous.append((i,j))
check(len(ambiguous)==len(extra),'independent_cosine_kNN_union_or_documented_roundoff_boundary')
check(all(abs(float(r['weight'])-(1-d[lookup[r['source']],lookup[r['target']]))<1e-12 for r in edges),'independent_edge_cosine_weights')
check(len(actual)==len(edges) and all(i!=j for i,j in actual),'simple_undirected_no_self_edges')
# Independent connected components and RB quality from saved edge weights.
adj=[set() for _ in ids];strength=np.zeros(len(ids));within=0.;total=0.
for r in edges:
 i,j=lookup[r['source']],lookup[r['target']];w=float(r['weight']);adj[i].add(j);adj[j].add(i);strength[i]+=w;strength[j]+=w;total+=w
 if labels[i]==labels[j]:within+=w
unseen=set(range(len(ids)));sizes=[]
while unseen:
 start=unseen.pop();visited={start};stack=[start]
 while stack:
  for j in adj[stack.pop()]:
   if j not in visited:visited.add(j);unseen.discard(j);stack.append(j)
 sizes.append(len(visited))
q=2*within-sum(float(strength[labels==label].sum())**2/(2*total) for label in set(labels))
grid=c.read(c.OUT/'stability/selected287_grid_summary.csv');baseline=next(r for r in grid if r['run_id']=='k20_r1_s42')
check(abs(q-float(baseline['rb_objective']))<1e-8,'independent_RB_configuration_objective')
check(abs(q/(2*total)-float(baseline['weighted_modularity_gamma1']))<1e-12,'independent_weighted_modularity')
diag=next(r for r in c.read(c.OUT/'graphs/selected287_diagnostics.csv') if r['k']=='20')
check(sorted(sizes,reverse=True)==list(map(int,diag['component_sizes'].split(';'))),'independent_connected_components')
check(abs(float(diag['density'])-2*len(actual)/(len(ids)*(len(ids)-1)))<1e-12,'independent_graph_density')
# Every prespecified grid run is present for all six cohorts, with one extra baseline.
for name in c.CONFIG['cohorts']:
 table=c.read(c.OUT/f'assignments/{name}_leiden_grid.csv');summ=c.read(c.OUT/f'stability/{name}_grid_summary.csv')
 check(len(summ)==81 and len(table[0])-1==81,name+'_complete_81_run_grid')
 check(len({r['sample_id'] for r in table})==len(table),name+'_unique_ids')
 check(all(r['converged']=='True' for r in summ),name+'_optimizer_convergence')
full=c.read(c.OUT/'assignments/selected_full15_leiden_grid.csv');check(len(full)==243 and all(source[r['sample_id']]['representative_n_genes']=='15' for r in full),'full15_selected_subset')
short=c.read(c.OUT/'tables/short_representatives_excluded_from_full15.csv');check(len(short)==44 and all(int(r['representative_n_genes'])<15 for r in short),'all_44_short_exclusions_recorded')
# Manual ARI using unordered pair contingency counts verifies shared-locus comparisons.
def choose2(n):return n*(n-1)/2
def ari(a,b):
 cells=collections.Counter(zip(a,b));ac=collections.Counter(a);bc=collections.Counter(b);n=len(a);tot=choose2(n);s=sum(choose2(v) for v in cells.values());aa=sum(choose2(v) for v in ac.values());bb=sum(choose2(v) for v in bc.values());expected=aa*bb/tot;den=(aa+bb)/2-expected
 return 1.0 if den==0 else (s-expected)/den
cross=c.read(c.OUT/'tables/cross_cohort_locus_memberships.csv');stats=c.read(c.OUT/'tables/cross_cohort_partition_agreement.csv')
for stat in stats:
 r=[row for row in cross if row['comparison']==stat['comparison']];value=ari([row['Leiden_cluster'] for row in r],[row['comparison_cluster'] for row in r]);check(abs(value-float(stat['ARI']))<1e-12,'manual_ARI_'+stat['comparison'])
check(all(float(r['nearest_reference_cosine'])==float(source[r['viral_locus_id']]['nearest_reference_cosine']) for r in rows),'prior_nearest_reference_values_preserved')
# Figure formats and resolution are checked without any re-fit or re-render.
export=json.loads((c.OUT/'provenance/figure_exports.json').read_text())
for entry in export['plots']:
 name=entry['name']
 with Image.open(c.OUT/'figures'/(name+'.png')) as im:check(im.size==tuple(int(round(v*600)) for v in entry['figure_inches']),name+'_600dpi_pixels')
 root=ET.parse(c.OUT/'figures'/(name+'.svg')).getroot();check(not root.findall('.//{http://www.w3.org/2000/svg}image'),name+'_vector_SVG')
for p in (c.OUT/'scripts').glob('*.py'):ast.parse(p.read_text(encoding='utf-8-sig'))
check(True,'source_syntax')
science=json.loads((c.OUT/'provenance/scientific_execution.json').read_text());check(science['umap_refitted']==False and science['embedding_inference_run']==False,'no_old_model_or_UMAP_refitting')
c.save('VALIDATION.json',dict(status='PASS',checks=checks,n_checks=len(checks),independent_graph_edge_differences_due_to_1e_12_roundoff_ties=len(extra),edge_count=len(actual),component_sizes=sorted(sizes,reverse=True),primary_RB_objective=q,scientific_clustering_refits_during_validation=False))
print('VALIDATION PASS',len(checks),'checks; independent graph, objective, alignment and ARI verified.')
