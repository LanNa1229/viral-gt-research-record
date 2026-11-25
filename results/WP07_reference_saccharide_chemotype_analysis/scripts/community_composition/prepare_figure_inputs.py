"""Verify frozen curated chemotypes and reconstruct three disjoint figure populations."""
import collections,csv,hashlib,json,sys
from pathlib import Path
from fractions import Fraction
csv.field_size_limit(50*1024*1024)
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[1];CUR=ROOT/'analysis/viral_gt_reference_chemotype_curation_20261027';ARCH=ROOT/'analysis/viral_gt_clustering_20261024_packaging/archive';UMAP=ROOT/'analysis/viral_gt_suppfig_global_leiden_redesign_20261027'
GROUPS=['GC08','GC05','Other GCs']
CATS=[
 ('Capsular polysaccharide','Capsule',['Capsular polysaccharide'],'#6F4E0B'),
 ('Extracellular/matrix polysaccharide','EPS/matrix polysaccharide',['Extracellular/matrix polysaccharide'],'#C8A44C'),
 ('LPS/LOS/O-antigen','LPS/LOS/O-antigen',['Lipopolysaccharide (LPS)','Lipooligosaccharide (LOS)','O-antigen'],'#D55E00'),
 ('S-layer/other surface glycan','S-layer/other surface glycan',['S-layer glycan','Other cell-surface glycan'],'#F3D7A5'),
 ('Aminoglycoside/aminocyclitol','Aminoglycoside/aminocyclitol',['Aminoglycoside / aminocyclitol-related metabolite'],'#0072B2'),
 ('Glycosylated specialized metabolite','Glycosylated specialized metabolite',['Glycosylated specialized metabolite'],'#B68BC4'),
 ('Other/unresolved','Other/unresolved',['Glycopeptidolipid','Teichuronic acid','Phosphonoglycan','Other carbohydrate-related product','Non-carbohydrate-associated product','Unresolved'],'#A8A8A8')]
def read(p):
 with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t' if p.suffix=='.tsv' else ','))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(n,x):(OUT/n).write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
def write(n,rows,fields=None):
 with (OUT/n).open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields or list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def main():
 sources={'chemotype_master':CUR/'reference_2502_chemotype_master.tsv','chemotype_dictionary':CUR/'CHEMOTYPE_DATA_DICTIONARY.md','curation_methods':CUR/'CURATION_METHODS.md','reference_QC':CUR/'REFERENCE_QC.md','prior_validation':CUR/'validation_results.json','curation_decisions':CUR/'CURATION_DECISIONS.json','curation_completion':CUR/'COMPLETION_RECEIPT.json','curation_output_manifest':CUR/'OUTPUT_MANIFEST_SHA256.tsv','curation_output_manifest_hash':CUR/'OUTPUT_MANIFEST_SHA256.tsv.sha256','curation_source_preflight':CUR/'SOURCE_PREFLIGHT.json','archive_manifest':ARCH/'ARCHIVE_SHA256.tsv','archive_pins':ROOT/'analysis/viral_gt_community_data_integration_20261027/INPUT_SHA256.tsv','global_assignments':ARCH/'assignments/global_joint_leiden_grid.csv','original_reference_labels':ARCH/'inputs/workspace/analysis/viral_gt_umap_20261024/joint_umap_samples.csv','UMAP_palette':UMAP/'COLOR_PALETTE.tsv','UMAP_plot_script':UMAP/'plot_integrated_umap.py','UMAP_config':UMAP/'config.json','Arial_font':Path('C:/Windows/Fonts/arial.ttf')}
 write('INPUT_MANIFEST_SHA256.tsv',[dict(role=k,path=str(p),size_bytes=p.stat().st_size,sha256=sha(p)) for k,p in sources.items()])
 assert sha(sources['curation_output_manifest'])==sources['curation_output_manifest_hash'].read_text().split()[0]
 oldhashes={r['path']:r['sha256'] for r in read(sources['curation_output_manifest'])}
 for k,p in sources.items():
  if p.is_relative_to(CUR) and k not in ['curation_output_manifest','curation_output_manifest_hash']:assert sha(p)==oldhashes[p.relative_to(CUR).as_posix()]
 completion=json.loads(sources['curation_completion'].read_text());assert sha(sources['chemotype_master'])==completion['master_sha256']=='5f86b03ff3e071eec4fcafa02eac36c981f1561777a1d1acbe6f69e51315cbc9'
 pins={r['sha256'] for r in read(sources['archive_pins']) if r['archive_path'].endswith('/ARCHIVE_SHA256.tsv')};assert pins=={sha(sources['archive_manifest'])}
 archived={r['path']:r['sha256'] for r in read(sources['archive_manifest'])}
 for key in ['global_assignments','original_reference_labels']:assert sha(sources[key])==archived[sources[key].relative_to(ARCH).as_posix()]
 prior=json.loads(sources['prior_validation'].read_text());assert prior['status']=='PASS' and prior['failed_checks']==0 and all(r['status']=='PASS' for r in prior['checks'])
 save('PRIOR_VALIDATION_REVIEW.json',dict(status=prior['status'],all_checks_parsed=len(prior['checks']),status_counts=dict(collections.Counter(r['status'] for r in prior['checks'])),check_family_counts=dict(collections.Counter(r['check'].split(':')[0] for r in prior['checks'])),subgroup_counts=prior['subgroup_counts'],dataset_accounting=prior['dataset_accounting'],historical_GC08=prior['historical_GC08'],new_chemical_review=False))
 master=read(sources['chemotype_master']);assert len(master)==2502 and len(master[0])==91 and len({r['mibig_accession'] for r in master})==2502
 grid={r['sample_id']:r['k20_r1_s42'] for r in read(sources['global_assignments'])}
 labels={r['sample_id']:r['reference_class'] for r in read(sources['original_reference_labels']) if r['sample_type']=='reference_BGC'}
 decisions={r['mibig_accession']:r for r in json.loads(sources['curation_decisions'].read_text())['decisions']}
 selected=[]
 for r in master:
  a=r['mibig_accession'];tokens=json.loads(r['original_class_multilabel'])
  assert r['global_leiden_community']==grid[a] and r['original_bgcprophet_class']==labels[a]
  assert tokens==labels[a].split() and (r['is_saccharide_reference']=='True')==('Saccharide' in tokens)
  if 'Saccharide' in tokens:selected.append(r)
 lookup={detail:(i,c[0]) for i,c in enumerate(CATS) for detail in c[2]}
 unknown=sorted({r['primary_chemotype'] for r in selected}-set(lookup))
 if unknown:
  save('DISCREPANCY_REPORT.json',dict(status='STOP',unknown_primary_chemotypes=unknown));raise SystemExit(2)
 mapping=[]
 for r in selected:
  a=r['mibig_accession'];assert r['source_version']=='MIBiG 3.1'
  for f in ['primary_chemotype','chemotype_annotation_confidence','curated_chemotype_status','curation_notes']:assert r[f]==decisions[a][f]
  i,category=lookup[r['primary_chemotype']];group=r['global_leiden_community'] if r['global_leiden_community'] in ['GC08','GC05'] else 'Other GCs'
  membership='IN_DEFINED_SET' if i<4 else 'UNRESOLVED_PRIMARY_CHEMOTYPE' if r['primary_chemotype']=='Unresolved' else 'OUTSIDE_DEFINED_SET'
  mapping.append({**r,'figure_population':group,'display_chemotype_group':category,'display_category_order':str(i+1),'four_category_surface_membership':membership})
 counts=collections.Counter((r['figure_population'],r['display_chemotype_group']) for r in mapping);sizes=collections.Counter(r['figure_population'] for r in mapping)
 # Historical numbers are validation targets only; plotting data are reconstructed above.
 expected={'GC08':[31,19,18,4,20,4,4],'GC05':[4,2,3,0,7,5,3],'Other GCs':[1,1,1,0,7,51,10]}
 actual={g:[counts[g,c[0]] for c in CATS] for g in GROUPS}
 bad={g:dict(expected=expected[g],observed=actual[g]) for g in GROUPS if actual[g]!=expected[g]}
 if len(mapping)!=195 or dict(sizes)!={'GC08':100,'GC05':24,'Other GCs':71}:bad['population_sizes']=dict(sizes)
 if bad:save('DISCREPANCY_REPORT.json',dict(status='STOP',discrepancies=bad));raise SystemExit(2)
 write('ACCESSION_LEVEL_CHEMOTYPE_MAPPING.tsv',mapping)
 rows=[]
 for g in GROUPS:
  left=Fraction(0)
  for i,(name,short,details,h) in enumerate(CATS):
   n=counts[g,name];pct=Fraction(100*n,sizes[g]);right=left+pct
   rows.append(dict(population=g,population_order=GROUPS.index(g)+1,chemotype_group=name,legend_label=short,category_order=i+1,n_BGCs=n,denominator_saccharide_BGCs=sizes[g],fraction=float(Fraction(n,sizes[g])),percentage=float(pct),left_percentage=float(left),right_percentage=float(right),accessions_json=json.dumps([r['mibig_accession'] for r in mapping if (r['figure_population'],r['display_chemotype_group'])==(g,name)],separators=(',',':'))))
   left=right
  assert left==100
 write('FIGURE_SOURCE_DATA.tsv',rows)
 write('CHEMOTYPE_GROUP_MAPPING.tsv',[dict(source_primary_chemotype=d,display_chemotype_group=c[0],legend_label=c[1],category_order=i+1,in_four_category_surface_definition=str(i<4)) for i,c in enumerate(CATS) for d in c[2]])
 write('COLOR_PALETTE.tsv',[dict(chemotype_group=c[0],legend_label=c[1],category_order=i+1,HEX=c[3],design_role='surface/extracellular operational set' if i<4 else 'other chemotype grouping',existing_UMAP_Saccharide_HEX='#C8A44C') for i,c in enumerate(CATS)])
 surface=[dict(population=g,n_four_category_surface=sum(actual[g][:4]),denominator=sizes[g],percentage=100*sum(actual[g][:4])/sizes[g],definition='display categories 1-4 only; excludes phosphonoglycan and broader source flags') for g in GROUPS]
 write('SURFACE_EXTRACELLULAR_SUMMARY.tsv',surface)
 unresolved=[r for r in mapping if r['primary_chemotype']=='Unresolved' or r['chemotype_annotation_confidence'].startswith('PROVISIONAL')]
 write('ANNOTATION_REVIEW_FLAGS.tsv',[{k:r[k] for k in ['mibig_accession','figure_population','global_leiden_community','primary_chemotype','display_chemotype_group','chemotype_annotation_confidence','curated_chemotype_status','manual_review_required','curation_notes']} for r in unresolved])
 config=dict(inputs={k:str(p) for k,p in sources.items()},populations=GROUPS,categories=[dict(name=c[0],short_label=c[1],source_primaries=c[2],HEX=c[3],order=i+1) for i,c in enumerate(CATS)],width_mm=183,height_mm=95,dpi=600,font_path=str(sources['Arial_font']),axis_title_pt=8.5,tick_pt=8,bar_label_pt=8.5,legend_pt=8,bar_height=.56,separator_linewidth_pt=.35,axes_bounds_mm=[32,42,145,43],x_limits=[0,100],x_ticks=[0,25,50,75,100],x_title='Saccharide reference BGCs (%)',surface_operational_definition='categories 1-4 only; phosphonoglycans remain category 7',reference_release='MIBiG 3.1',Leiden_baseline=dict(k=20,resolution=1.0,seed=42),Leiden_archival_commit=completion['private_archive_commit'],no_new_clustering=True,no_new_annotation=True,no_enrichment_testing=True)
 save('config.json',config)
 save('INPUT_VALIDATION.json',dict(status='PASS',master_rows=2502,master_columns=91,unique_saccharide_BGCs=195,population_sizes=dict(sizes),seven_group_counts=actual,counts_match_brief=True,hybrid_saccharide_references=sum(r['original_is_hybrid']=='True' for r in mapping),nonoverlapping_populations=True,surface_summary=surface,annotation_confidence_counts=dict(collections.Counter(r['chemotype_annotation_confidence'] for r in mapping)),unresolved_or_provisional_accessions=[r['mibig_accession'] for r in unresolved],new_models=False,new_annotations=False))
 print((OUT/'INPUT_VALIDATION.json').read_text(),flush=True)
if __name__=='__main__':main()
