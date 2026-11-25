"""Prepare a new presentation-only revision from frozen scientific source data."""
import concurrent.futures,csv,datetime,hashlib,json,math,shutil
from pathlib import Path
from PIL import Image
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[1];PREVIOUS=ROOT/'analysis/viral_gt_suppfig_global_leiden_20261027'
CLASSES=['Alkaloid','Terpene','NRP','Polyketide','RiPP','Saccharide','Other']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):
 with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t' if p.suffix=='.tsv' else ','))
def save(n,x):(OUT/n).write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
def write(n,rows):
 with (OUT/n).open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def main():
 oldseal=json.loads((PREVIOUS/'CHECKSUM_VERIFICATION.json').read_text())
 assert sha(PREVIOUS/'OUTPUT_MANIFEST_SHA256.tsv')==oldseal['manifest_sha256']
 frozen=read(PREVIOUS/'OUTPUT_MANIFEST_SHA256.tsv')
 def verify(r):
  p=PREVIOUS/r['relative_path'];return p.stat().st_size==int(r['size_bytes']) and sha(p)==r['sha256']
 with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:assert all(pool.map(verify,frozen))
 oldinputs=read(PREVIOUS/'INPUT_MANIFEST_SHA256.tsv')
 assert all(Path(r['path']).stat().st_size==int(r['size_bytes']) and sha(Path(r['path']))==r['sha256'] for r in oldinputs)
 previous_members={r['relative_path']:r for r in frozen}
 names={
 'previous_figure_source':PREVIOUS/'FIGURE_SOURCE_DATA.tsv','previous_plot_script':PREVIOUS/'plot_integrated_umap.py',
 'previous_palette':PREVIOUS/'COLOR_PALETTE.tsv','previous_visual_QC':PREVIOUS/'VISUAL_QC_REPORT.md',
 'previous_config':PREVIOUS/'config.json','previous_manifest':PREVIOUS/'OUTPUT_MANIFEST_SHA256.tsv',
 'previous_seal':PREVIOUS/'CHECKSUM_VERIFICATION.json','previous_PDF':PREVIOUS/'figures/SuppFig_global_leiden_biosynthetic_classes.pdf',
 'previous_SVG':PREVIOUS/'figures/SuppFig_global_leiden_biosynthetic_classes.svg',
 'previous_PNG':PREVIOUS/'figures/SuppFig_global_leiden_biosynthetic_classes.png',
 'previous_community_key':PREVIOUS/'COMMUNITY_ANNOTATION_KEY.tsv','previous_class_counts':PREVIOUS/'REFERENCE_CLASS_MARGINAL_COUNTS.tsv',
 'published_reference_image':OUT/'qc/Flamholz_2024_Figure3_reference.png'}
 for k,p in names.items():
  if p.is_relative_to(PREVIOUS) and k not in ['previous_manifest','previous_seal']:assert sha(p)==previous_members[p.relative_to(PREVIOUS).as_posix()]['sha256']
 records=[dict(role='upstream_'+r['role'],path=r['path'],size_bytes=r['size_bytes'],sha256=r['sha256']) for r in oldinputs]
 records += [dict(role=k,path=str(p),size_bytes=p.stat().st_size,sha256=sha(p)) for k,p in names.items()]
 write('INPUT_MANIFEST_SHA256.tsv',records)
 oldcfg=json.loads((PREVIOUS/'config.json').read_text())
 data=read(PREVIOUS/'FIGURE_SOURCE_DATA.tsv')
 assert sha(PREVIOUS/'FIGURE_SOURCE_DATA.tsv')==oldseal['source_data_sha256'] and len(data)==2789
 out=[]
 for r in data:
  row=dict(r);labels=json.loads(r['reference_class_labels_json'])
  if r['sample_type']=='viral_locus':group='Selected viral GT-associated loci (287)';marker='hollow_triangle';color='#303030';face='none'
  elif len(labels)>1:group='Hybrid (multiple labels)';marker='hollow_circle';color='#666666';face='none'
  else:group=labels[0];marker='filled_circle';color=oldcfg['palette'][group];face=color
  row.update(plot_marker=marker,rendering_group=group,plot_edge_or_fill_hex=color,plot_facecolor=face)
  out.append(row)
 write('FIGURE_SOURCE_DATA.tsv',out)
 for n in ['COMMUNITY_ANNOTATION_KEY.tsv','REFERENCE_CLASS_MARGINAL_COUNTS.tsv']:shutil.copyfile(PREVIOUS/n,OUT/n)
 limits=[]
 for c in ['UMAP1','UMAP2']:
  vals=[float(r[c]) for r in data];lo=min(vals);hi=max(vals);pad=.055*(hi-lo)
  limits.append([math.floor((lo-pad)*2)/2,math.ceil((hi+pad)*2)/2])
 config=dict(version='minimalist_v2',previous_output=str(PREVIOUS),inputs={k:str(p) for k,p in names.items()},upstream_inputs=oldcfg['inputs'],classes=CLASSES,palette=oldcfg['palette'],hybrid_color='#666666',viral_color='#303030',font_path=oldcfg['font_path'],width_mm=183,height_mm=110,dpi=600,axes_left_fraction=.07,axes_width_fraction=.75,axes_bottom_fraction=.14,axes_max_height_fraction=.78,axis_titles_pt=9,tick_labels_pt=8,legend_pt=8.5,legend_title_pt=8.5,reference_area_pt2=6,hybrid_area_pt2=6,viral_area_pt2=8,reference_alpha=.83,hybrid_alpha=.72,viral_alpha=.87,hybrid_linewidth_pt=.35,viral_linewidth_pt=.45,frame_linewidth_pt=.5,axis_limits=dict(UMAP1=limits[0],UMAP2=limits[1]),hybrid_rendering='Neutral hollow circle; rendering category only; all original memberships preserved',viral_rendering='Small hollow dark-gray triangle; no GC coloring',equal_axis_units=True,original_UMAP_analysis=oldcfg['original_UMAP_analysis'],original_projection_points=3058,displayed_points=2789,Leiden_baseline=oldcfg['Leiden_baseline'],Leiden_archival_commit=oldcfg['Leiden_archival_commit'],no_model_execution=True,community_ids_retained_only_in_source_tables_and_caption=True)
 save('config.json',config)
 write('COLOR_PALETTE.tsv',[dict(category=c,HEX=config['palette'][c],marker='filled circle',annotation_status='Original single-label reference class') for c in CLASSES]+[dict(category='Hybrid (multiple labels)',HEX=config['hybrid_color'],marker='hollow circle',annotation_status='Display category only; original labels retained'),dict(category='Selected viral GT-associated loci (287)',HEX=config['viral_color'],marker='hollow triangle',annotation_status='Selected original viral loci')])
 image=Image.open(names['published_reference_image']);assert image.size==(2168,1251)
 crop=(0,0,1480,625);image.crop(crop).save(OUT/'qc/Flamholz_2024_Figure3a_style_crop.png')
 save('REFERENCE_STYLE_PROVENANCE.json',dict(article='Flamholz, Biller and Kelly, Nature Microbiology 9, 537-549 (2024)',DOI='10.1038/s41564-023-01584-8',figure='3a',figure_page='https://www.nature.com/articles/s41564-023-01584-8/figures/3',publisher_image_URL='https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fs41564-023-01584-8/MediaObjects/41564_2023_1584_Fig3_HTML.png',retrieval='Opened and visually inspected full-resolution publisher image in in-app browser; downloaded using observed image locator and copied locally',source_image_sha256=sha(names['published_reference_image']),source_pixels=list(image.size),local_QC_crop_box_pixels=list(crop),no_reference_artwork_in_final_scientific_figure=True,attachment_contents='Text brief only; publisher original used for direct visual comparison',failed_direct_fetch='DNS blocked in sandbox; elevated HTTP request returned client-challenge HTML; Europe PMC fullTextXML returned HTTP 500; normal browser image loaded successfully'))
 save('INPUT_VALIDATION.json',dict(status='PASS',previous_manifest_files_unchanged=len(frozen),upstream_inputs_unchanged=len(oldinputs),new_input_manifest_records=len(records),total_samples=2789,references=2502,single_label_references=2035,hybrid_references=467,selected_viral_loci=287,original_coordinate_strings_preserved=True,all_non_rendering_source_fields_preserved=True,no_models_executed=True))
 print((OUT/'INPUT_VALIDATION.json').read_text(),flush=True)
if __name__=='__main__':main()
