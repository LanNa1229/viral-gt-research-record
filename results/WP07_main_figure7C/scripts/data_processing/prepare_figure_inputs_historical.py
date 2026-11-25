"""Verify frozen sources and prepare two presentation-only UMAP designs."""
import collections,csv,hashlib,json,sys
from pathlib import Path
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[1]
sys.path.append(str(ROOT/'analysis/bgcprophet_reference/.deps'))
import numpy as np
A=ROOT/'analysis/viral_gt_clustering_20261024_packaging/archive'
U=ROOT/'analysis/viral_gt_umap_20261024';I=ROOT/'analysis/viral_gt_community_data_integration_20261027'
C=ROOT/'analysis/viral_gt_reference_chemotype_curation_20261027'
S=ROOT/'analysis/viral_gt_suppfig_chemotype_20261027';R=ROOT/'analysis/viral_gt_suppfig_global_leiden_redesign_20261027'
VCS={'VC01':68,'VC02':55,'VC03':74,'VC04':35,'VC05':55}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):
 with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t' if p.suffix=='.tsv' else ','))
def save(n,x):(OUT/n).write_text(json.dumps(x,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
def write(n,rows):
 with (OUT/n).open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def unique(rows,k):
 d={}
 for r in rows:assert r[k] and r[k] not in d;d[r[k]]=r
 return d
def main():
 named={
 'original_UMAP':U/'joint_umap_samples.csv','original_saved_coordinates':U/'tables/joint_umap_coordinates.npy',
 'original_plot_script':U/'scripts/exploratory_umap.py','original_UMAP_provenance':U/'execution_provenance.json','original_UMAP_checksums':U/'output_checksums.csv','original_UMAP_methods':U/'README.md',
 'viral_master':I/'viral_287_integrated_community_mapping.tsv','viral_dictionary':I/'DATA_DICTIONARY.md','viral_input_pins':I/'INPUT_SHA256.tsv','viral_output_manifest':I/'OUTPUT_SHA256.tsv','private_archive_identity':I/'PRIVATE_ARCHIVE_IDENTITY.json',
 'chemotype_master':C/'reference_2502_chemotype_master.tsv','chemotype_dictionary':C/'CHEMOTYPE_DATA_DICTIONARY.md','curation_methods':C/'CURATION_METHODS.md','curation_completion':C/'COMPLETION_RECEIPT.json','curation_manifest':C/'OUTPUT_MANIFEST_SHA256.tsv','curation_manifest_pin':C/'OUTPUT_MANIFEST_SHA256.tsv.sha256',
 'archive_manifest':A/'ARCHIVE_SHA256.tsv','global_assignments':A/'assignments/global_joint_leiden_grid.csv','viral_assignments':A/'assignments/selected287_leiden_grid.csv','clustering_methods':A/'METHODS_RECORD.md','clustering_config':A/'config.json','archived_UMAP':A/'inputs/workspace/analysis/viral_gt_umap_20261024/joint_umap_samples.csv',
 'SupFig2_mapping':S/'CHEMOTYPE_GROUP_MAPPING.tsv','SupFig2_accession_data':S/'ACCESSION_LEVEL_CHEMOTYPE_MAPPING.tsv','SupFig2_palette':S/'COLOR_PALETTE.tsv','SupFig2_data':S/'FIGURE_SOURCE_DATA.tsv','SupFig2_dictionary':S/'DATA_DICTIONARY.md','SupFig2_methods':S/'METHODS_RECORD.md','SupFig2_config':S/'config.json','SupFig2_manifest':S/'OUTPUT_MANIFEST_SHA256.tsv','SupFig2_seal':S/'CHECKSUM_VERIFICATION.json',
 'SupUMAP_dictionary':R/'DATA_DICTIONARY.md','SupUMAP_methods':R/'METHODS_RECORD.md','SupUMAP_config':R/'config.json','SupUMAP_palette':R/'COLOR_PALETTE.tsv','SupUMAP_manifest':R/'OUTPUT_MANIFEST_SHA256.tsv','SupUMAP_seal':R/'CHECKSUM_VERIFICATION.json','SupUMAP_style_PNG':R/'qc/preview_183mm_150dpi.png','publisher_style_crop':R/'qc/Flamholz_2024_Figure3a_style_crop.png','publisher_style_provenance':R/'REFERENCE_STYLE_PROVENANCE.json',
 'first_SupUMAP_dictionary':ROOT/'analysis/viral_gt_suppfig_global_leiden_20261027/DATA_DICTIONARY.md',
 'Arial_font':Path('C:/Windows/Fonts/arial.ttf'),'Arial_bold_font':Path('C:/Windows/Fonts/arialbd.ttf')}
 write('INPUT_MANIFEST_SHA256.tsv',[dict(role=k,path=str(p),size_bytes=p.stat().st_size,sha256=sha(p)) for k,p in named.items()])
 assert sha(named['original_UMAP'])=='edcbfaeadbc07540e6e605530b2d51cdfa91d1b502dad8a1f5051874376560ba'
 assert sha(named['viral_master'])=='c3b4d6c8967051d5a1ec488137f777096e3d7c49779154fa9f97c59208060e4d'
 assert sha(named['chemotype_master'])=='5f86b03ff3e071eec4fcafa02eac36c981f1561777a1d1acbe6f69e51315cbc9'
 assert sha(named['curation_manifest'])=='54a32a9e523faee4d08473a8d023635ad2e077739a5b9de18fc187171f09d4b3'
 pins=read(named['viral_input_pins']);assert {r['sha256'] for r in pins if r['archive_path'].endswith('/ARCHIVE_SHA256.tsv')}=={sha(named['archive_manifest'])}
 ah={r['path']:r['sha256'] for r in read(named['archive_manifest'])}
 for p in named.values():
  if p.is_relative_to(A) and p.name!='ARCHIVE_SHA256.tsv':assert sha(p)==ah[p.relative_to(A).as_posix()]
 assert sha(named['original_UMAP'])==sha(named['archived_UMAP'])
 uh={r['path']:r['sha256'] for r in read(named['original_UMAP_checksums'])}
 for role in ['original_UMAP','original_saved_coordinates','original_UMAP_provenance']:
  p=named[role];assert sha(p)==uh[p.relative_to(U).as_posix()]
 ch={r['path']:r['sha256'] for r in read(named['curation_manifest'])}
 for p in named.values():
  if p.is_relative_to(C) and p.name not in ['OUTPUT_MANIFEST_SHA256.tsv','OUTPUT_MANIFEST_SHA256.tsv.sha256']:assert sha(p)==ch[p.relative_to(C).as_posix()]
 vh={r['path']:r['sha256'] for r in read(named['viral_output_manifest'])}
 for p in named.values():
  if p.is_relative_to(I) and p.name!='OUTPUT_SHA256.tsv':assert sha(p)==vh[p.relative_to(I).as_posix()]
 preservation=[]
 for base,expected in [(S,'7bd488b4bdaf9e0fcc2942668f3840b2f1a76ad5bed5aec8778e2c07ca8f159d'),(R,'a7dd32df3b54c2f0ebde67b4272c1f86f7428cd97163bd21bcd53a0a648ecdc0'),(ROOT/'analysis/viral_gt_suppfig_global_leiden_20261027','a269c7cd59e6291f33b2c23cd8dab2008bc75e261ccb102699760462fac9c39d')]:
  assert sha(base/'OUTPUT_MANIFEST_SHA256.tsv')==expected
  for r in read(base/'OUTPUT_MANIFEST_SHA256.tsv'):
   p=base/r['relative_path'];assert sha(p)==r['sha256'];assert p.stat().st_size==int(r['size_bytes'])
   preservation.append(dict(path=str(p),sha256=r['sha256'],size_bytes=r['size_bytes'],before='PASS'))
 write('PREVIOUS_FIGURE_PRESERVATION.tsv',preservation)
 urows=read(named['original_UMAP']);u=unique(urows,'sample_id');assert len(u)==3058
 xy=np.load(named['original_saved_coordinates'],allow_pickle=False)
 assert xy.shape==(3058,2) and np.array_equal(xy,np.array([[float(r['UMAP1']),float(r['UMAP2'])] for r in urows]))
 vi=unique(read(named['viral_master']),'viral_locus_id');cm=unique(read(named['chemotype_master']),'mibig_accession')
 gg=unique(read(named['global_assignments']),'sample_id');vg=unique(read(named['viral_assignments']),'sample_id')
 ss=unique(read(named['SupFig2_accession_data']),'mibig_accession');mapping=unique(read(named['SupFig2_mapping']),'source_primary_chemotype')
 palette=read(named['SupFig2_palette']);pc={r['chemotype_group']:r for r in palette}
 assert len(vi)==287 and len(cm)==2502 and len(ss)==195
 assert dict(collections.Counter(r['viral_only_leiden_cluster'] for r in vi.values()))==VCS
 assert {s for s,r in u.items() if r['sample_type']=='reference_BGC'}==set(cm)
 assert set(vi)<=set(u) and set(gg)==set(cm)|set(vi)
 vc_colors={'VC01':'#003F5C','VC02':'#00866A','VC03':'#9E2F50','VC04':'#514478','VC05':'#566B30'}
 config=dict(width_mm=183,A_height_mm=150,B_height_mm=111,dpi=600,font_path=str(named['Arial_font']),bold_font_path=str(named['Arial_bold_font']),axis_titles_pt=8.5,tick_pt=7.5,legend_pt=7.5,legend_title_pt=8,panel_title_pt=8.5,axis_limits={'UMAP1':[-6.5,16.5],'UMAP2':[-1.5,12]},x_ticks=[-5,0,5,10,15],y_ticks=[0,5,10],frame_linewidth_pt=.5,A_axes_mm=[12,46,165,165*13.5/23],B_axes_mm=[[12,56,76,76*13.5/23],[104,56,76,76*13.5/23]],reference_categories=[dict(name=r['chemotype_group'],legend_label=r['legend_label'],HEX=r['HEX'],order=int(r['category_order'])) for r in palette],viral_categories=[dict(name=vc,HEX=hc,n=VCS[vc],fill='none' if vc in ['VC04','VC05'] else 'solid') for vc,hc in vc_colors.items()],non_saccharide_HEX='#DEDEDE',context_saccharide_HEX='#D1BD86',context_viral_HEX='#717171',original_projection_n=3058,original_reference_n=2502,original_viral_n=556,displayed_n=2789,omitted_unselected_viral_n=269,original_space_dimensions=320,viral_Leiden_baseline='k20_r1_s42',archive_commit='284d4f1655cf036c1890bdffcb4aa2dc93f6e41e',selected_design=None,no_fitting=True,local_only=True)
 save('config.json',config)
 rows=[];keep=['chemotype_annotation_confidence','curated_chemotype_status','manual_review_required','secondary_chemotype_tags','curation_notes','source_record_URL','source_record_SHA256','source_version']
 # Cached identity receipts avoid repeatedly hashing the large source master per row.
 provenance_entries={k:dict(path=str(p),sha256=sha(p)) for k,p in named.items()}
 for sid in sorted(set(cm)|set(vi)):
  o=u[sid];isref=o['sample_type']=='reference_BGC';cr=cm.get(sid,{});v=vi.get(sid,{});labs=o['reference_class'].split() if isref else []
  sac='Saccharide' in labs;group=mapping[cr['primary_chemotype']]['display_chemotype_group'] if sac else ''
  assert all(o[k] for k in ['UMAP1','UMAP2']) and np.isfinite([float(o['UMAP1']),float(o['UMAP2'])]).all()
  gc=gg[sid]['k20_r1_s42']
  if isref:
   assert cr['original_bgcprophet_class']==o['reference_class'] and json.loads(cr['original_class_multilabel'])==labs
   assert cr['global_leiden_community']==gc and (cr['is_saccharide_reference']=='True')==sac
   if sac:
    assert group==ss[sid]['display_chemotype_group']
    for k,val in cr.items():assert ss[sid][k]==val
  else:
   assert sid==o['viral_locus_id']==v['viral_locus_id']
   assert v['UMAP1']==o['UMAP1'] and v['UMAP2']==o['UMAP2'] and v['global_joint_leiden_cluster']==gc
   assert v['viral_only_leiden_cluster']==vg[sid]['k20_r1_s42']
  hc=pc[group]['HEX'] if sac else (config['non_saccharide_HEX'] if isref else vc_colors[v['viral_only_leiden_cluster']])
  provenance={k:provenance_entries[k] for k in (['original_UMAP','chemotype_master','global_assignments'] if isref else ['original_UMAP','viral_master','viral_assignments','global_assignments'])}
  r=dict(sample_id=sid,sample_type=o['sample_type'],UMAP1=o['UMAP1'],UMAP2=o['UMAP2'],reference_class_original=o['reference_class'],reference_class_labels_json=json.dumps(labs,separators=(',',':')),is_saccharide_reference=str(sac),is_hybrid_reference=str(len(labs)>1),reference_primary_chemotype=cr.get('primary_chemotype',''),reference_chemotype_visualization_group=group,viral_locus_id=o['viral_locus_id'],viral_only_leiden_community=v.get('viral_only_leiden_cluster',''),global_leiden_community=gc,saccharide_joint_leiden_community=v.get('saccharide_joint_leiden_cluster',''),plotted_marker_type='circle' if isref else 'triangle',plotted_HEX_color=hc,A_fill_style='solid' if isref else next(x['fill'] for x in config['viral_categories'] if x['name']==v['viral_only_leiden_cluster']),B_left_HEX_color=hc if isref else config['context_viral_HEX'],B_right_HEX_color=config['context_saccharide_HEX'] if sac else hc,data_provenance=json.dumps(provenance,separators=(',',':')))
  r.update({k:cr.get(k,'') for k in keep});rows.append(r)
 write('FIGURE_SOURCE_DATA.tsv',rows)
 counts=collections.Counter(r['reference_chemotype_visualization_group'] for r in rows if r['is_saccharide_reference']=='True');priorcounts=collections.Counter()
 for r in read(named['SupFig2_data']):priorcounts[r['chemotype_group']]+=int(r['n_BGCs'])
 assert counts==priorcounts and sum(counts.values())==195
 assert len(rows)==2789 and sum(r['sample_type']=='reference_BGC' for r in rows)==2502
 write('CHEMOTYPE_COUNTS.tsv',[dict(chemotype_group=r['chemotype_group'],HEX=r['HEX'],n_references=counts[r['chemotype_group']]) for r in palette])
 write('ANNOTATION_REVIEW_FLAGS.tsv',[r for r in rows if r['is_saccharide_reference']=='True' and r['chemotype_annotation_confidence'] in ['UNRESOLVED','PROVISIONAL_RECORD_LEVEL']])
 write('FINAL_PALETTE.tsv',[dict(section='Reference chemotype',category=r['chemotype_group'],HEX=r['HEX'],marker='circle',fill='solid') for r in palette]+[dict(section='Viral Leiden community',category=r['name'],HEX=r['HEX'],marker='triangle',fill=r['fill']) for r in config['viral_categories']]+[dict(section='Context',category='Other reference BGCs',HEX=config['non_saccharide_HEX'],marker='circle',fill='solid'),dict(section='Context',category='Saccharide landmarks in B right',HEX=config['context_saccharide_HEX'],marker='circle',fill='solid'),dict(section='Context',category='Selected viral context in B left',HEX=config['context_viral_HEX'],marker='triangle',fill='none')])
 save('INPUT_VALIDATION.json',dict(status='PASS',inputs_hashed_before_join=len(named),prior_figure_files_preserved=len(preservation),original_coordinate_array_exact=True,source_coordinate_strings_exact=True,original_projection_counts=dict(collections.Counter(r['sample_type'] for r in urows)),displayed_counts={'reference_BGC':2502,'viral_locus':287},saccharide=195,other_reference=2307,saccharide_hybrids=sum(r['is_saccharide_reference']=='True' and r['is_hybrid_reference']=='True' for r in rows),chemotype_counts=dict(counts),VC_counts=VCS,all_joins='Exact unique sample_id / mibig_accession / viral_locus_id; never row order',annotations_changed=0,new_assignments=0))
 print(json.dumps(json.loads((OUT/'INPUT_VALIDATION.json').read_text()),indent=2))
if __name__=='__main__':main()
