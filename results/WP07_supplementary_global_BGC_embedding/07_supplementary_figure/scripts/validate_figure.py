"""Validate the presentation-only redesign against frozen and original sources."""
import concurrent.futures,csv,datetime,hashlib,json,re,subprocess,sys,xml.etree.ElementTree as ET
from pathlib import Path
from PIL import Image
from pypdf import PdfReader
OUT=Path(__file__).resolve().parent;CONFIG=json.loads((OUT/'config.json').read_text());PREVIOUS=Path(CONFIG['previous_output'])
CHECKS=[]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):
 with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t' if p.suffix=='.tsv' else ','))
def unique(rr,key):
 d={r[key]:r for r in rr};assert len(d)==len(rr);return d
def save(n,x):(OUT/n).write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
def check(name,value,evidence):
 CHECKS.append(dict(check=name,status='PASS' if value else 'FAIL',evidence=evidence))
 if not value:raise AssertionError(name+': '+str(evidence))
def md(n,s):(OUT/n).write_text(s.strip()+'\n',encoding='utf-8')
def table(headers,rows):
 return '| '+' | '.join(headers)+' |\n| '+' | '.join(['---']*len(headers))+' |\n'+''.join('| '+' | '.join(str(v).replace('|',' / ') for v in r)+' |\n' for r in rows)
def main():
 data=read(OUT/'FIGURE_SOURCE_DATA.tsv');new=unique(data,'sample_id')
 old=unique(read(Path(CONFIG['inputs']['previous_figure_source'])),'sample_id')
 joint=unique(read(Path(CONFIG['upstream_inputs']['saved_joint_UMAP'])),'sample_id')
 grid=unique(read(Path(CONFIG['upstream_inputs']['global_assignments'])),'sample_id')
 master=unique(read(Path(CONFIG['upstream_inputs']['integrated_viral_master'])),'viral_locus_id')
 check('population_and_exact_ID_sets',len(new)==2789 and len(joint)==3058 and set(new)==set(old)==set(grid) and {sid for sid,r in new.items() if r['sample_type']=='viral_locus'}==set(master) and len(master)==287,'All 2502 references and all 287 selected viral IDs retained; original projection has 3058 samples; 269 unselected loci remain omitted')
 science=[c for c in next(iter(old.values())) if c!='plot_marker']
 check('all_original_scientific_source_fields',all(all(r[c]==old[sid][c] for c in science) for sid,r in new.items()),dict(unchanged_columns=science,rendering_only_changes=['plot_marker','rendering_group','plot_edge_or_fill_hex','plot_facecolor']))
 check('saved_coordinates_unchanged',all(all(r[c]==joint[sid][c] for c in ['UMAP1','UMAP2']) for sid,r in new.items()),'5578 coordinate strings exactly match the original saved joint UMAP; no refit or jitter')
 check('original_Leiden_assignments',all(r['global_leiden_community']==grid[sid]['k20_r1_s42'] for sid,r in new.items()) and all(r['global_leiden_community']==master[sid]['global_joint_leiden_cluster'] for sid,r in new.items() if r['sample_type']=='viral_locus'),'All ten baseline GC IDs unchanged; original-space k=20, resolution=1.0, seed=42')
 refs={sid:r for sid,r in new.items() if r['sample_type']=='reference_BGC'}
 hybrid={sid:r for sid,r in refs.items() if len(json.loads(r['reference_class_labels_json']))>1}
 single={sid:r for sid,r in refs.items() if len(json.loads(r['reference_class_labels_json']))==1}
 check('reference_classes_and_multilabel_memberships',len(refs)==2502 and len(single)==2035 and len(hybrid)==467 and all(r['reference_class_original']==joint[sid]['reference_class'] and json.loads(r['reference_class_labels_json'])==joint[sid]['reference_class'].split() for sid,r in refs.items()),'All original class strings and JSON memberships retained; neutral hybrid display does not replace scientific annotations')
 check('hollow_rendering_categories',all(r['plot_marker']=='hollow_circle' and r['plot_facecolor']=='none' and r['rendering_group']=='Hybrid (multiple labels)' for r in hybrid.values()) and all(r['plot_marker']=='hollow_triangle' and r['plot_facecolor']=='none' for r in new.values() if r['sample_type']=='viral_locus') and all(r['rendering_group']==json.loads(r['reference_class_labels_json'])[0] and r['plot_marker']=='filled_circle' for r in single.values()),'467 neutral open circles and 287 neutral open triangles; single-label references use their original class colors')
 actual=unique(read(OUT/'PLOTTED_SAMPLE_RECEIPTS.tsv'),'sample_id')
 check('rendered_sample_receipts',set(actual)==set(new) and all(all(actual[sid][c]==r[c] for c in ['UMAP1','UMAP2','global_leiden_community','rendering_group']) and actual[sid]['original_class_labels_json']==r['reference_class_labels_json'] for sid,r in new.items()),'One plotted glyph per biological sample; 2789 unique render receipts')
 stem=OUT/'figures/SuppFig_global_leiden_biosynthetic_classes_minimalist_v2'
 ns={'s':'http://www.w3.org/2000/svg'};root=ET.parse(stem.with_suffix('.svg')).getroot()
 expected=read(OUT/'RENDER_GROUP_RECEIPTS.tsv');render_counts={}
 for r in expected:
  gs=[g for g in root.findall('.//s:g',ns) if g.attrib.get('id')==r['group_id']];assert len(gs)==1
  g=gs[0];n=len(g.findall('.//s:use',ns))+len(g.findall('s:path',ns))
  assert n==int(r['n_samples']);render_counts[r['group_id']]=n
  if r['group_id'] in ['samples_Hybrid','samples_SelectedViral']:assert re.search(r'fill:\s*none',ET.tostring(g,encoding='unicode'))
 check('actual_SVG_point_counts_and_open_fills',sum(render_counts.values())==2789 and render_counts['samples_Hybrid']==467 and render_counts['samples_SelectedViral']==287,render_counts)
 texts=[''.join(t.itertext()) for t in root.findall('.//s:text',ns)]
 check('no_community_labels_statistics_or_notes',not any(re.search(r'\bGC\d{2}\b|%|n/N|Global Leiden|Dominant reference|heterogeneous|fractions',t) for t in texts),'Artwork contains only numeric axes, UMAP titles, reference legend and selected-viral entry; no GC labels, leaders, table or detailed statistics')
 normalized=' '.join(' '.join(t.split()) for t in texts)
 check('simple_legend_labels',all(c in texts for c in CONFIG['classes']) and 'Hybrid' in texts and '(multiple labels)' in texts and 'Selected viral' in texts and 'GT-associated' in texts and 'loci (287)' in texts and 'Reference BGCs' in texts and '(2,502)' in texts,'Seven original single-label color entries plus Hybrid (multiple labels) and Selected viral GT-associated loci (287), wrapped for readability')
 layout=json.loads((OUT/'LAYOUT_RECEIPTS.json').read_text())
 check('full_cartesian_axes_and_numerical_ticks',layout['all_four_spines_visible'] and not layout['grid_visible'] and len(layout['visible_x_ticks'])>=4 and len(layout['visible_y_ticks'])>=4 and 'UMAP 1' in texts and 'UMAP 2' in texts,'Full rectangular black frame, numerical x/y ticks and visible UMAP 1 / UMAP 2 titles; no grid')
 check('plot_dominates_width',layout['scientific_axes']==1 and abs(layout['axes_width_fraction']-.75)<1e-9,dict(plot_width_percent=100*layout['axes_width_fraction'],plot_size_mm=[layout['axes_width_mm'],layout['axes_height_mm']],equal_units=layout['equal_axis_units']))
 for r in new.values():
  for c in ['UMAP1','UMAP2']:assert layout['axis_limits'][c][0]<float(r[c])<layout['axis_limits'][c][1]
 check('no_sample_clipped_by_axis_limits',True,layout['axis_limits'])
 check('readable_requested_typography',all(r['font_pt']==(9 if r['kind']=='axis_title' else 8 if r['kind'].endswith('tick') else 8.5) for r in layout['text']),'Arial: axis titles 9 pt, ticks 8 pt, legend/title 8.5 pt; current user-specified targets control this revision')
 w,h=layout['canvas_pixels']
 boundsok=all(r['x0']>=0 and r['x1']<=w and r['y0']>=0 and r['y1']<=h for r in layout['text'])
 overlap=[]
 for i,a in enumerate(layout['text']):
  for b in layout['text'][i+1:]:
   dx=min(a['x1'],b['x1'])-max(a['x0'],b['x0']);dy=min(a['y1'],b['y1'])-max(a['y0'],b['y0'])
   if dx>.5 and dy>.5:overlap.append([a['text'],b['text'],dx,dy])
 check('text_bounds_and_overlap',boundsok and not overlap,overlap or 'All final-size text bounds fit the canvas with no nontrivial overlaps')
 reader=PdfReader(stem.with_suffix('.pdf'));assert len(reader.pages)==1
 page=reader.pages[0];mm=[float(page.mediabox.width)*25.4/72,float(page.mediabox.height)*25.4/72]
 fonts=[]
 for name,ref in page['/Resources']['/Font'].items():
  f=ref.get_object()
  for child in f.get('/DescendantFonts',[f]):
   d=child.get_object();descriptor=d.get('/FontDescriptor');descriptor=descriptor.get_object() if descriptor else {}
   fonts.append(dict(base_font=str(f.get('/BaseFont')),descendant_type=str(d.get('/Subtype')),TrueType_embedded='/FontFile2' in descriptor))
 check('fixed_size_editable_vector_PDF',abs(mm[0]-183)<1e-6 and abs(mm[1]-110)<1e-6 and len(page.images)==0 and all(f['TrueType_embedded'] and f['descendant_type']=='/CIDFontType2' and 'Arial' in f['base_font'] for f in fonts),dict(dimensions_mm=mm,fonts=fonts,raster_images=len(page.images)))
 check('editable_vector_SVG',not root.findall('.//s:image',ns) and len(texts)>=20,'SVG text remains text; all point and frame marks are vector; no embedded raster images or publisher artwork')
 im=Image.open(stem.with_suffix('.png'));dpi=im.info.get('dpi')
 check('600dpi_PNG',list(im.size)==[int(183/25.4*600),int(110/25.4*600)] and dpi and abs(dpi[0]-600)<.1,dict(pixels=list(im.size),dpi=dpi))
 upstreamcomp={(r['joint_community'],r['reference_class']):r for r in read(Path(CONFIG['upstream_inputs']['reference_class_composition'])) if r['analysis']=='global_joint'}
 localcomp=read(OUT/'REFERENCE_CLASS_MARGINAL_COUNTS.tsv');key=unique(read(OUT/'COMMUNITY_ANNOTATION_KEY.tsv'),'community')
 for r in localcomp:
  n=sum(x['global_leiden_community']==r['community'] and r['reference_class'] in json.loads(x['reference_class_labels_json']) for x in refs.values())
  assert n==int(r['n_carrying_label'])==int(upstreamcomp[r['community'],r['reference_class']]['n_reference_with_label'])
  assert abs(float(r['fraction'])-n/int(r['n_references']))<1e-14
 check('community_statistics_preserved_outside_artwork',len(localcomp)==70 and len(key)==10 and int(key['GC08']['n_dominant_label'])==100 and int(key['GC08']['n_references'])==105 and int(key['GC05']['n_dominant_label'])==24 and int(key['GC05']['n_references'])==53,'All ten assignment records and all 70 marginal class counts/fractions preserved; GC08 100/105 Saccharide, GC05 24/53')
 inp=read(OUT/'INPUT_MANIFEST_SHA256.tsv')
 def verify(r):
  p=Path(r['path']);return p.stat().st_size==int(r['size_bytes']) and sha(p)==r['sha256']
 with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:inputpass=all(pool.map(verify,inp))
 prior=read(PREVIOUS/'OUTPUT_MANIFEST_SHA256.tsv')
 with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:priorpass=all(pool.map(lambda r:sha(PREVIOUS/r['relative_path'])==r['sha256'],prior))
 check('inputs_and_previous_figure_preserved',inputpass and priorpass,dict(frozen_input_records=len(inp),previous_package_files_unchanged=len(prior),previous_PDF_unchanged=sha(Path(CONFIG['inputs']['previous_PDF']))==json.loads((PREVIOUS/'CHECKSUM_VERIFICATION.json').read_text())['PDF_sha256']))
 ref=Image.open(OUT/'qc/Flamholz_2024_Figure3_reference.png')
 check('published_style_reference_and_accessibility_views',ref.size==(2168,1251) and all((OUT/('qc/'+n+'_preview.png')).exists() for n in ['protanopia','deuteranopia','tritanopia','grayscale']),'Full-resolution publisher Figure 3 retrieved and visually inspected; Figure 3a crop compared directly; normal, three CVD and grayscale views inspected')
 command=['C:/Users/archival-user/.cache/codex-runtimes/codex-primary-runtime/dependencies/native/poppler/Library/bin/pdftoppm.exe','-r','300','-singlefile','-png',str(stem.with_suffix('.pdf')),str(OUT/'qc/PDF_render_300dpi')]
 record=dict(command=command,started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),PDF_sha256=sha(stem.with_suffix('.pdf')))
 save('PDF_RENDER_EXECUTION.json',record)
 result=subprocess.run(command,capture_output=True,text=True)
 record.update(exit_code=result.returncode,stdout=result.stdout,stderr=result.stderr,finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),PNG_sha256=sha(OUT/'qc/PDF_render_300dpi.png'))
 save('PDF_RENDER_EXECUTION.json',record)
 check('independent_Poppler_PDF_render',result.returncode==0 and Image.open(OUT/'qc/PDF_render_300dpi.png').width>2100,'Final PDF independently rendered at 300 dpi for visual review')
 save('VALIDATION_RESULTS.json',dict(status='PASS',finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),checks=CHECKS,models_run=False,clustering_run=False,coordinates_changed=False,previous_figure_overwritten=False,committed=False,pushed=False))
 make_reports()
 print(json.dumps(dict(status='PASS',checks=len(CHECKS),displayed_samples=len(new),physical_dimensions_mm=mm,previous_package_preserved=True),indent=2),flush=True)

def make_reports():
 md('SUPPLEMENTARY_FIGURE_LEGEND.md',"""# Supplementary figure legend

**Biosynthetic class organization of reference BGCs and selected viral GT-associated loci in learned representation space.** The figure displays 2,502 reference biosynthetic gene clusters (BGCs) and 287 selected viral glycosyltransferase-associated loci using the previously fitted joint UMAP of original 320-dimensional BGC-Prophet embeddings. The original projection contained 2,502 references and 556 viral loci; this view retains the same coordinates for all references and the 287 selected loci, omitting 269 unselected viral loci. Filled reference circles are colored by their original single biosynthetic class annotations. The 467 reference BGCs with multiple class annotations are shown as neutral hollow circles for visual clarity; this rendering does not constitute an additional biological class, and all original memberships are retained in the source data. Selected viral loci are hollow dark-gray triangles. Separate global Leiden clustering in the original 320-dimensional space (k = 20, resolution = 1.0, seed = 42) yielded ten baseline communities whose reference annotations showed substantial descriptive correspondence with broad biosynthetic classes. GC08 contained 100/105 Saccharide-labeled references, whereas the mixed GC05 community contained 24/53. Community assignments and marginal reference-class composition are supplied in the supplementary source-data tables; multilabel fractions may sum above 100%. Other is a heterogeneous broad annotation. UMAP is a visualization, not the clustering procedure. Correspondence within the supervised BGC-Prophet representation does not independently validate viral functional predictions.
""")
 executions=[]
 for p in sorted((OUT/'executions').glob('*/before.json')):
  b=json.loads(p.read_text());executions.append([str(p.parent.relative_to(OUT)),b['started_utc'],json.dumps(b['command'])])
 palette=read(OUT/'COLOR_PALETTE.tsv')
 palettetable=table(['Display category','HEX','Marker'],[[r['category'],r['HEX'],r['marker']] for r in palette])
 md('METHODS_RECORD.md',"""# Methods record: minimalist UMAP redesign

## Presentation-only revision

This new sibling directory preserves the previous figure package at analysis/viral_gt_suppfig_global_leiden_20261027/. All 81 files in that package's sealed output manifest were checksum verified before and after this task, including the original PDF/SVG/PNG, plotting code, source data and QC records. No previous file was overwritten. No scientific model, UMAP, Leiden clustering or embeddings were rerun. No Git commit or push was performed.

The input manifest records 34 sources: 21 prior frozen scientific/font inputs plus the prior figure/source/palette/QC/configuration/manifests and the observed publisher reference image. Their SHA-256 values were recorded before rendering. The updated TSV preserves all previous scientific fields exactly. Only rendering metadata change: plot_marker and new rendering_group, plot_edge_or_fill_hex and plot_facecolor fields. All stable IDs, original coordinate strings, class annotations, hybrid memberships and baseline GC IDs remain unchanged.

## Data and lineage

Original UMAP analysis: viral_gt_umap_20261024. The fitted projection contains 3,058 samples: 2,502 references and 556 viral loci. The current visualization uses the same saved coordinates for all 2,502 references and exactly the 287 selected viral locus IDs, giving 2,789 displayed glyphs. The 269 unselected viral loci remain omitted. Filtering does not recompute or transform coordinates.

Global Leiden assignments are read from the existing original-320D baseline column k20_r1_s42 (k=20, resolution=1.0, seed=42). Private archival commit: 284d4f1655cf036c1890bdffcb4aa2dc93f6e41e in LanNa1229/viral-gt-internal-analyses/results/WP07_viral_gt_figure3C_clustering_20261024/. All assignments were independently compared with the archived grid and, for selected viral loci, the integrated viral master. UMAP geometry does not define these communities.

The original seven broad reference labels remain Alkaloid, Terpene, NRP, Polyketide, RiPP, Saccharide and Other. Exactly 2,035 references have a single label and use that label's filled-circle color. All 467 multilabel references use the same neutral hollow-circle rendering. Their original class strings and JSON lists remain unchanged. Hybrid is a rendering category, not a new biological class or a substituted single annotation. The community key and all 70 marginal reference-class composition rows were copied unchanged and independently cross-checked against the original annotation table. They are separate supporting data, absent from the artwork.

## Graphical specifications

Canvas: 183 x 110 mm. One data axes occupies exactly 75% of total figure width (137.25 mm), with an actual height of 80.56 mm after equal-unit aspect adjustment. The plot has a full rectangular black frame with 0.5 pt lines, bottom/left numeric ticks and UMAP 1 / UMAP 2 titles. No background grid. Axes are equal in physical units. Limits are calculated from the displayed saved coordinate extrema with 5.5% padding, expanded to half-unit boundaries: x [-6.5, 16.5], y [-1.5, 12.0]. All displayed sample centers lie strictly inside those limits. No jitter or coordinate movement is used.

Arial typography follows the user's revised readability targets: axis titles 9 pt, tick labels 8 pt, legend/header 8.5 pt. These explicit revised targets take precedence over the earlier figure's 5-7 pt styling. There is no in-artwork title or panel letter.

Filled single-label reference circles have area 6 pt squared, opacity 0.83 and no outlines. Hollow hybrid circles have area 6 pt squared, neutral #666666 outlines 0.35 pt wide and opacity 0.72. Hollow viral triangles have area 8 pt squared, dark-gray #303030 outlines 0.45 pt wide and opacity 0.87. These point sizes encode sample positions only. Hollow markers reduce occlusion without removing overlapping samples.

The vertical right-hand legend has the reference heading and count, seven single-label colors, Hybrid (multiple labels), and one Selected viral GT-associated loci (287) entry. Long labels are wrapped rather than shrunk. The artwork contains no GC01-GC10 labels, leader lines, community annotation table, percentages, n/N values, bottom notes or other headings. Community results are retained in the separate caption and source tables.

"""+palettetable+"""

## Reference-style comparison

The current attachment contained the text brief only. The full-resolution published [Flamholz et al. (2024), Figure 3](https://www.nature.com/articles/s41564-023-01584-8/figures/3) was opened in the in-app browser and visually inspected. The observed image was downloaded through its image locator, copied to qc/Flamholz_2024_Figure3_reference.png and hashed. A local QC crop isolates panel 3a for direct comparison. REFERENCE_STYLE_PROVENANCE.json records the actual retrieval and crop.

The revised graphical organization follows panel 3a's full rectangular numeric axes, dominant landscape scatter, small points, white background and simple vertical class legend. The original paper's data, exact colors and other panels are not used in the new scientific figure. Equal-unit axes preserve our saved projection geometry; the data shape and exact frame proportions therefore remain specific to this dataset. The publisher image is a local style-reference asset only and is not embedded in the delivered PDF/SVG/PNG.

An initial direct HTTP image fetch returned a client-challenge page and the public Europe PMC XML endpoint returned an error. These did not supply the comparison image. The normal browser subsequently displayed the full 2,168 x 1,251-pixel publisher image successfully. No CAPTCHA was solved, authentication altered or paywall bypassed.

## Export and accessibility review

Matplotlib exports a fully vector PDF using font type 42, an SVG retaining text nodes (svg.fonttype=none), and a 600-dpi PNG. Export uses the fixed canvas without tight cropping. The PDF page size and embedded Arial TrueType text were inspected independently; no raster images appear in either vector output. Actual SVG point instances were counted per explicit sample group, including all hybrid circles and viral triangles.

Protanopia, deuteranopia and tritanopia previews use the severity-100 Machado et al. matrices already recorded in the prior plotting work, applied in linear sRGB with gamut clipping. A grayscale preview is also retained. The seven class colors remain the previous final HEX map; filled versus hollow circles and triangle shape distinguish single-label references, hybrids and viral loci. Some class colors remain less distinct under severe CVD or grayscale; the display is not an accessibility certification. No class membership is inferred from neutral hybrid color.

The reference image, normal 150/96-dpi previews and all accessibility views were inspected. The final vector PDF is independently rendered at 300 dpi using Poppler, with its actual command and input/output hashes recorded in PDF_RENDER_EXECUTION.json. The final PDF visual review is recorded in VISUAL_QA.json.

## Executed commands and reproducibility

"""+table(['Frozen execution directory','Start UTC','Actual child argument array'],executions)+"""

Each run preserves the exact code/configuration, pre-execution command array, Python/platform version, input-manifest hash, stdout, stderr and final exit status under executions/. SOFTWARE_VERSIONS.json records the actual plotting runtime and resolved Arial font. The final config contains all display choices; rendering is deterministic.

VALIDATION_RESULTS.json records the source, glyph, typography, layout and editable-vector checks. INPUT_MANIFEST_SHA256.tsv and the final OUTPUT_MANIFEST_SHA256.tsv preserve source and output hashes. The checksum manifest and its own verification receipt exclude themselves to avoid recursive hashes. Any later GitHub upload requires explicit approval; all current files remain local.
""")
 md('DATA_DICTIONARY.md',"""# Source-data and display dictionary

FIGURE_SOURCE_DATA.tsv has one row for each of 2,789 unique displayed samples. The twelve non-rendering fields from the preceding figure source are preserved exactly: sample_id, sample_type, global_leiden_community, UMAP1, UMAP2, reference_class_original, reference_class_labels_json, n_reference_class_labels, is_hybrid_reference, viral_locus_id, genome_unit_id and representative_window_id. Stable identifiers are join keys; original UMAP decimal strings and reference label ordering are unmodified.

Rendering fields:
- plot_marker: filled_circle for a single-label reference, hollow_circle for a hybrid reference, hollow_triangle for a selected viral locus.
- rendering_group: original single reference class, Hybrid (multiple labels), or Selected viral GT-associated loci (287). The hybrid label describes presentation only and never replaces original scientific annotations.
- plot_edge_or_fill_hex: fill HEX for single-label references or outline HEX for hollow glyphs.
- plot_facecolor: class HEX or none.

The seven original functional labels and all 467 multilabel memberships remain in reference_class_original and reference_class_labels_json. A hollow hybrid glyph does not indicate missing labels.

COMMUNITY_ANNOTATION_KEY.tsv and REFERENCE_CLASS_MARGINAL_COUNTS.tsv are byte-identical copies of the preceding verified tables. They retain ten original GC labels, per-community reference/viral counts, dominant marginal labels and all 70 class memberships/fractions. Multilabel marginal fractions can exceed 100% in sum. These tables support the separate figure legend and are not rendered as an in-artwork table.

PLOTTED_SAMPLE_RECEIPTS.tsv records each actual plotted ID, display group, marker, unchanged coordinate strings and original class memberships. RENDER_GROUP_RECEIPTS.tsv records the nine plotted groups and their marker parameters. LAYOUT_RECEIPTS.json records actual axes geometry, ticks, font sizes and measured text bounds. COLOR_PALETTE.tsv includes seven annotation colors plus neutral hybrid and viral display colors.
""")
 md('VISUAL_QC_REPORT.md',"""# Visual QC: minimalist supplementary UMAP v2

The source, layout and vector-export checks passed. Final visual-review hashes and any material limitations are recorded in VISUAL_QA.json.

## Direct comparison with Flamholz Figure 3a

The full publisher image was successfully viewed in the browser and saved for local comparison. The panel 3a crop was inspected beside the revised render. Both use a dominant landscape scatter, complete rectangular numerical axes, small functional-category points, a white background and a vertical right-hand legend. The new figure uses this dataset's unchanged saved coordinates, original class vocabulary and existing palette.

The redesigned plot occupies 75% of figure width, with a compact legend instead of a numerical annotation column. Our axes retain equal physical units, preserving the existing point-cloud geometry; its exact aspect and data shape differ from the reference. The prior figure remains intact in its original directory.

## Final-size assessment

- Full X/Y frame and numerical labels are visible. UMAP 1 and UMAP 2 titles are clear; no grid is present.
- Arial 9 pt axis titles, 8 pt tick labels and 8.5 pt legend text are readable in the fixed 183 x 110 mm layout.
- The legend uses seven filled class circles, one hollow hybrid circle and one hollow viral triangle. Long labels wrap within the reserved margin.
- No GC labels, leader lines, community table, percentages, n/N values or explanatory bottom notes remain in the artwork.
- Single-label reference color groups are immediately visible. Hollow hybrid circles reduce the former pie-sector complexity while preserving the location of all 467 multilabel references.
- All 287 selected viral loci are hollow triangles, preserving visibility of reference colors beneath them. Dense exact/near overlaps remain; no samples were moved or jittered.
- Normal, protanopia, deuteranopia, tritanopia and grayscale previews were inspected. Sample types remain distinguishable by shape/fill. Color discrimination among all seven classes is reduced in severe deficiency simulations and grayscale; hybrid memberships must be read from source data.
- The PDF is independently rendered at 300 dpi for final review; 150-dpi and 96-dpi previews retain the same physical layout. Screen previews do not constitute a calibrated print.

## Integrity and interpretation

All 2,502 reference BGCs and 287 selected viral loci are retained, with unchanged saved UMAP coordinates, original reference labels and ten original-space GC assignments. Hybrid is a neutral rendering choice, not a new biological class. All community statistics remain in the separate caption and source tables.

The layout communicates the broad biosynthetic organization of the learned representation more directly than the earlier annotated version. This is descriptive correspondence within a supervised representation; it does not independently validate viral functions. UMAP was not used as a clustering algorithm, and no scientific model was rerun.
""")
 md('README.md',"""# Minimalist supplementary UMAP, version 2

A new presentation-only revision with full numerical Cartesian axes, larger Arial text, a dominant scatter panel and a compact vertical legend. The previous figure package is preserved.

## Revised figure

- [Editable SVG](figures/SuppFig_global_leiden_biosynthetic_classes_minimalist_v2.svg)
- [Vector PDF](figures/SuppFig_global_leiden_biosynthetic_classes_minimalist_v2.pdf)
- [600-dpi PNG](figures/SuppFig_global_leiden_biosynthetic_classes_minimalist_v2.png)
- [150-dpi preview](qc/preview_183mm_150dpi.png)
- [96-dpi preview](qc/preview_183mm_96dpi.png)
- [Independent PDF render](qc/PDF_render_300dpi.png)

Canvas 183 x 110 mm; plot width 75%. Arial axis titles 9 pt, ticks 8 pt, legend 8.5 pt. All 2,502 references and 287 selected viral loci remain. The 467 hybrids are neutral hollow circles; viral loci are hollow triangles. Original multilabel annotations and all GC memberships remain in the source data.

## Supporting records

- [Publication-style figure legend](SUPPLEMENTARY_FIGURE_LEGEND.md)
- [Figure source data](FIGURE_SOURCE_DATA.tsv)
- [Original community summaries, outside artwork](COMMUNITY_ANNOTATION_KEY.tsv)
- [All reference-class marginal counts](REFERENCE_CLASS_MARGINAL_COUNTS.tsv)
- [Final HEX palette](COLOR_PALETTE.tsv)
- [Data dictionary](DATA_DICTIONARY.md)
- [Methods record](METHODS_RECORD.md)
- [Visual QC report](VISUAL_QC_REPORT.md)
- [Machine-readable validation](VALIDATION_RESULTS.json)

prepare_figure_inputs.py verifies the prior sealed package and freezes sources. plot_integrated_umap.py renders only the saved data using config.json. validate_figure.py independently checks scientific fields, original coordinates/assignments, actual SVG glyph counts, axes, fonts and export structure. run_stage.py stores contemporaneous code/config snapshots, actual commands, software/platform, stdout/stderr and exit status under executions/.

INPUT_MANIFEST_SHA256.tsv lists 34 input records. PLOTTED_SAMPLE_RECEIPTS.tsv and RENDER_GROUP_RECEIPTS.tsv record every point and display group. LAYOUT_RECEIPTS.json, SOFTWARE_VERSIONS.json, PDF_RENDER_EXECUTION.json and VISUAL_QA.json retain layout and visual evidence. finalize_package.py seals OUTPUT_MANIFEST_SHA256.tsv and CHECKSUM_VERIFICATION.json after all execution files exist.

qc/ contains final-size, CVD and grayscale previews plus the actual publisher image/crop used for direct style comparison. The publisher image is a local reference asset and is not embedded in the scientific exports.

No UMAP/Leiden/model analysis was rerun. No previous figure was overwritten. Nothing was committed or pushed. Any later upload requires explicit scoped approval.
""")
if __name__=='__main__':main()
