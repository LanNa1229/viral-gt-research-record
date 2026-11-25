"""Independent accession, denominator, category, geometry and export validation."""
import collections,concurrent.futures,csv,datetime,hashlib,json,re,subprocess,sys,xml.etree.ElementTree as ET
from pathlib import Path
from fractions import Fraction
from PIL import Image
from pypdf import PdfReader
csv.field_size_limit(50*1024*1024)
OUT=Path(__file__).resolve().parent;CONFIG=json.loads((OUT/'config.json').read_text());CHECKS=[]
def read(p):
 with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t' if p.suffix=='.tsv' else ','))
def unique(rr,k):
 d={r[k]:r for r in rr};assert len(d)==len(rr);return d
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(n,x):(OUT/n).write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
def md(n,s):(OUT/n).write_text(s.strip()+'\n',encoding='utf-8')
def table(h,rows):return '| '+' | '.join(h)+' |\n| '+' | '.join(['---']*len(h))+' |\n'+''.join('| '+' | '.join(str(x).replace('|',' / ') for x in r)+' |\n' for r in rows)
def check(n,ok,evidence):
 CHECKS.append(dict(check=n,status='PASS' if ok else 'FAIL',evidence=evidence))
 if not ok:raise AssertionError(n+': '+str(evidence))
def main():
 master=read(Path(CONFIG['inputs']['chemotype_master']));mm=unique(master,'mibig_accession')
 selected={a:r for a,r in mm.items() if 'Saccharide' in json.loads(r['original_class_multilabel'])}
 mapping=unique(read(OUT/'ACCESSION_LEVEL_CHEMOTYPE_MAPPING.tsv'),'mibig_accession')
 check('exact_195_unique_original_Saccharide_BGCs',len(master)==2502 and len(master[0])==91 and len(selected)==195 and set(mapping)==set(selected),'All and only original Saccharide members, including hybrids; no duplicate accessions')
 check('all_91_original_fields_preserved',all(all(r[k]==selected[a][k] for k in selected[a]) for a,r in mapping.items()),'Every copied field remains exact, including detailed chemistry, confidence, source URLs/hashes and review notes')
 populations={g:{a for a,r in mapping.items() if r['figure_population']==g} for g in CONFIG['populations']}
 check('mutually_exclusive_exhaustive_populations',[len(populations[g]) for g in CONFIG['populations']]==[100,24,71] and set.union(*populations.values())==set(selected) and all(not populations[a]&populations[b] for i,a in enumerate(CONFIG['populations']) for b in CONFIG['populations'][i+1:]),'100 GC08 + 24 GC05 + 71 other GCs = 195; GC05 is excluded from Other GCs')
 grid={r['sample_id']:r['k20_r1_s42'] for r in read(Path(CONFIG['inputs']['global_assignments']))}
 labels={r['sample_id']:r['reference_class'] for r in read(Path(CONFIG['inputs']['original_reference_labels'])) if r['sample_type']=='reference_BGC'}
 check('original_communities_and_labels_unchanged',all(r['global_leiden_community']==grid[a] and r['original_bgcprophet_class']==labels[a] and json.loads(r['original_class_multilabel'])==labels[a].split() for a,r in mapping.items()),'Exact archived k20_r1_s42 communities and original broad class strings/memberships')
 check('hybrids_and_unresolved_retained',sum(r['original_is_hybrid']=='True' for r in mapping.values())==65 and {a for a,r in mapping.items() if r['primary_chemotype']=='Unresolved'}=={'BGC0000262','BGC0000727'} and {a for a,r in mapping.items() if r['chemotype_annotation_confidence'].startswith('PROVISIONAL')}=={'BGC0001765'},'65 hybrid records, two unresolved and one provisional source assignment retained; unresolved is not a negative')
 decisions=unique(json.loads(Path(CONFIG['inputs']['curation_decisions']).read_text())['decisions'],'mibig_accession')
 check('primary_chemotypes_not_reassigned',all(r['primary_chemotype']==decisions[a]['primary_chemotype'] and r['chemotype_annotation_confidence']==decisions[a]['chemotype_annotation_confidence'] for a,r in mapping.items()),'Primary chemistry and confidence match frozen curation decisions; no new annotation')
 detailmap=read(OUT/'CHEMOTYPE_GROUP_MAPPING.tsv');lookup=unique(detailmap,'source_primary_chemotype')
 check('exactly_one_display_group_per_BGC',len(lookup)==15 and all(r['display_chemotype_group']==lookup[r['primary_chemotype']]['display_chemotype_group'] for r in mapping.values()),'Explicit total mapping of 15 detailed source categories to seven mutually exclusive display groups; no fuzzy classification or missing-value imputation')
 counts=collections.Counter((r['figure_population'],r['display_chemotype_group']) for r in mapping.values());src=read(OUT/'FIGURE_SOURCE_DATA.tsv')
 expected={'GC08':[31,19,18,4,20,4,4],'GC05':[4,2,3,0,7,5,3],'Other GCs':[1,1,1,0,7,51,10]}
 actual={g:[counts[g,c['name']] for c in CONFIG['categories']] for g in CONFIG['populations']}
 check('historical_expectations_independently_reproduced',actual==expected,actual)
 for g in CONFIG['populations']:
  rows=[r for r in src if r['population']==g];left=Fraction(0)
  assert len(rows)==7 and len({r['chemotype_group'] for r in rows})==7
  for i,r in enumerate(rows):
   n=counts[g,CONFIG['categories'][i]['name']];N=len(populations[g]);value=Fraction(100*n,N)
   assert int(r['n_BGCs'])==n and int(r['denominator_saccharide_BGCs'])==N
   assert set(json.loads(r['accessions_json']))=={a for a in populations[g] if mapping[a]['display_chemotype_group']==r['chemotype_group']}
   assert abs(float(r['percentage'])-float(value))<1e-12 and abs(float(r['left_percentage'])-float(left))<1e-12
   left+=value;assert abs(float(r['right_percentage'])-float(left))<1e-12
  assert left==100
 check('all_21_source_cells_and_100_percent_totals',len(src)==21,'Accession-derived numerators, Saccharide-only denominators and exact rational cumulative widths checked; structural zeros retained')
 surface=read(OUT/'SURFACE_EXTRACELLULAR_SUMMARY.tsv')
 for g,n in [('GC08',72),('GC05',9),('Other GCs',3)]:
  ids={a for a in populations[g] if int(mapping[a]['display_category_order'])<=4}
  assert len(ids)==n
  assert not any(mapping[a]['primary_chemotype']=='Phosphonoglycan' for a in ids)
  r=next(r for r in surface if r['population']==g);assert int(r['n_four_category_surface'])==n and abs(float(r['percentage'])-100*n/len(populations[g]))<1e-12
 check('surface_operational_definition_consistent',True,'Categories 1-4: 72/100, 9/24, 3/71. GC08 phosphonoglycans remain in Other/unresolved; broader source flags are preserved but not counted')
 receipts=read(OUT/'PLOTTED_SEGMENT_RECEIPTS.tsv');assert len(receipts)==21
 for r,s in zip(receipts,src):
  assert (r['population'],r['chemotype_group'])==(s['population'],s['chemotype_group'])
  assert abs(float(r['width_percentage'])-float(s['percentage']))<1e-12 and abs(float(r['left_percentage'])-float(s['left_percentage']))<1e-12
 check('plotted_segment_values_match_source',True,'All 21 segment records, including two zeros, match exported source data in stable category order')
 stem=OUT/'figures/SuppFig2_saccharide_reference_chemotype_distribution'
 ns={'s':'http://www.w3.org/2000/svg'};root=ET.parse(stem.with_suffix('.svg')).getroot()
 for r in receipts:
  gs=[g for g in root.findall('.//s:g',ns) if g.attrib.get('id')==r['gid']];assert len(gs)==1
  path=gs[0].find('s:path',ns);assert path is not None
  vals=[float(x) for x in re.findall(r'[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?',path.attrib['d'])]
  assert len(vals)%2==0
  xs=vals[::2];actual_width=max(xs)-min(xs);expected_width=float(r['width_percentage'])/100*CONFIG['axes_bounds_mm'][2]*72/25.4
  assert abs(actual_width-expected_width)<2e-5,(r['gid'],actual_width,expected_width)
 check('actual_SVG_segment_geometry',True,'Each SVG rectangle width independently matches its percentage; zero segments stay zero; smallest positive segment is 1/71 = 1.40845% (~2.04 mm)')
 layout=json.loads((OUT/'LAYOUT_RECEIPTS.json').read_text())
 texts=[''.join(t.itertext()) for t in root.findall('.//s:text',ns)]
 check('axis_labels_scale_and_order',layout['x_limits']==[0,100] and layout['x_ticks']==[0,25,50,75,100] and layout['y_order_top_to_bottom']==['GC08','GC05','Other GCs'] and CONFIG['x_title'] in texts and all(str(t)+'%' in texts for t in [0,25,50,75,100]),'GC08, GC05, Other GCs top-to-bottom; x axis 0-100% at 25-point intervals')
 w,h=layout['canvas_pixels'];overlap=[]
 for i,a in enumerate(layout['text']):
  assert a['x0']>=0 and a['x1']<=w and a['y0']>=0 and a['y1']<=h
  for b in layout['text'][i+1:]:
   if min(a['x1'],b['x1'])-max(a['x0'],b['x0'])>.5 and min(a['y1'],b['y1'])-max(a['y0'],b['y0'])>.5:overlap.append([a['text'],b['text']])
 check('typography_and_text_layout',not overlap and all(7.5<=r['font_pt']<=9 for r in layout['text']) and layout['in_bar_text_count']==0,overlap or 'Arial 8-8.5 pt; no clipping/overlapping labels and no crowded in-bar text')
 check('seven_category_legend_and_fixed_colors',all(c['short_label'] in texts for c in CONFIG['categories']) and len({c['HEX'] for c in CONFIG['categories']})==7 and CONFIG['categories'][1]['HEX']=='#C8A44C','Seven solid, distinct HEX colors; four warm surface groups remain consecutive; UMAP Saccharide gold retained for EPS')
 reader=PdfReader(stem.with_suffix('.pdf'));assert len(reader.pages)==1;p=reader.pages[0]
 dims=[float(p.mediabox.width)*25.4/72,float(p.mediabox.height)*25.4/72];fonts=[]
 for ref in p['/Resources']['/Font'].values():
  f=ref.get_object()
  for child in f.get('/DescendantFonts',[f]):
   d=child.get_object();desc=d.get('/FontDescriptor');desc=desc.get_object() if desc else {}
   fonts.append(dict(base_font=str(f.get('/BaseFont')),type=str(d.get('/Subtype')),TrueType_embedded='/FontFile2' in desc))
 check('editable_vector_exports_and_dimensions',abs(dims[0]-183)<1e-6 and abs(dims[1]-95)<1e-6 and not p.images and not root.findall('.//s:image',ns) and len(texts)>=15 and all(f['TrueType_embedded'] and f['type']=='/CIDFontType2' and 'Arial' in f['base_font'] for f in fonts),dict(dimensions_mm=dims,fonts=fonts,SVG_text_nodes=len(texts)))
 im=Image.open(stem.with_suffix('.png'));dpi=im.info['dpi']
 check('600dpi_PNG',list(im.size)==[int(183/25.4*600),int(95/25.4*600)] and abs(dpi[0]-600)<.1,dict(pixels=list(im.size),dpi=dpi))
 inp=read(OUT/'INPUT_MANIFEST_SHA256.tsv')
 with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:good=list(pool.map(lambda r:sha(Path(r['path']))==r['sha256'] and Path(r['path']).stat().st_size==int(r['size_bytes']),inp))
 check('all_frozen_inputs_unchanged',all(good),dict(inputs_rehashed=len(inp),original_master_SHA256=sha(Path(CONFIG['inputs']['chemotype_master']))))
 check('color_vision_and_grayscale_previews',all((OUT/('qc/'+m+'_preview.png')).exists() for m in ['protanopia','deuteranopia','tritanopia','grayscale']),'Three CVD and grayscale views generated and inspected; exact display simulation matrices and pairwise diagnostics preserved')
 command=['/redacted-user/.cache/codex-runtimes/codex-primary-runtime/dependencies/native/poppler/Library/bin/pdftoppm.exe','-r','300','-singlefile','-png',str(stem.with_suffix('.pdf')),str(OUT/'qc/PDF_render_300dpi')]
 record=dict(command=command,started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),PDF_SHA256=sha(stem.with_suffix('.pdf')));save('PDF_RENDER_EXECUTION.json',record)
 r=subprocess.run(command,capture_output=True,text=True);record.update(exit_code=r.returncode,stdout=r.stdout,stderr=r.stderr,finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),PNG_SHA256=sha(OUT/'qc/PDF_render_300dpi.png'));save('PDF_RENDER_EXECUTION.json',record)
 check('independent_PDF_render',r.returncode==0 and Image.open(OUT/'qc/PDF_render_300dpi.png').width>2100,'Poppler independently rendered the final vector PDF for final-size visual review')
 save('VALIDATION_RESULTS.json',dict(status='PASS',checks=CHECKS,finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),n_BGCs=195,counts_match_historical=True,chemical_annotations_changed=False,community_assignments_changed=False,models_rerun=False,enrichment_tests=False,committed=False,pushed=False))
 make_reports(src,surface)
 print(json.dumps(dict(status='PASS',checks=len(CHECKS),unique_BGCs=195,population_sizes=[100,24,71],all_bars_total_100_percent=True),indent=2),flush=True)

def make_reports(src,surface):
 counts=table(['Chemotype group','GC08 (100)','GC05 (24)','Other GCs (71)'],[[c['short_label'],*[next(r['n_BGCs'] for r in src if r['population']==g and r['chemotype_group']==c['name']) for g in CONFIG['populations']]] for c in CONFIG['categories']])
 palette=table(['Chemotype group','Final HEX'],[[c['short_label'],c['HEX']] for c in CONFIG['categories']])
 commands=[]
 for p in sorted((OUT/'executions').glob('*/before.json')):
  b=json.loads(p.read_text());commands.append([str(p.parent.relative_to(OUT)),b['started_utc'],json.dumps(b['command'])])
 md('SUPPLEMENTARY_FIGURE_LEGEND.md',"""# Supplementary figure legend

**Supplementary Fig. 2 | Chemotype distribution of Saccharide reference BGCs across global embedding communities.** The 195 reference biosynthetic gene clusters (BGCs) carrying the original BGC-Prophet Saccharide label, including 65 multilabel references, were partitioned into three nonoverlapping populations using existing global Leiden assignments: GC08 (n = 100), GC05 (n = 24) and all other communities (n = 71). Leiden clustering was previously performed in the original 320-dimensional embedding space (k = 20, resolution = 1.0, seed = 42). Bars show within-population percentages of seven mutually exclusive groups consolidated from the existing primary chemotype annotations based on frozen MIBiG 3.1 records. Original detailed chemotypes, hybrid class memberships and evidence/review states are retained in the accession-level source data. The operational surface/extracellular glycan grouping comprises capsule, extracellular/matrix polysaccharide, LPS/LOS/O-antigen and S-layer/other surface glycan; it excludes phosphonoglycans and is distinct from broader source flags. This grouping accounts for 72/100 GC08 references (72.0%), 9/24 GC05 references (37.5%) and 3/71 references in other communities (4.2%). Other/unresolved combines specified miscellaneous chemotypes and unresolved annotations; it is not a biological-negative category. Assignments are source-based proposed annotations, not uniformly independently human-verified: the displayed cohort includes two unresolved and one provisional record-level assignment. These distributions are descriptive, without enrichment testing or independent validation of the supervised representation, and do not assign the corresponding products to viral loci.
""")
 md('METHODS_RECORD.md',"""# Methods record: Saccharide reference chemotype figure

## Scope and controlling sources

This is source-data verification and visualization only. No new chemical annotation, clustering, embeddings, model inference, nearest-neighbor calculation or enrichment test was run. All original analysis/curation files were read without modification. No commit or push was performed.

The controlling source is reference_2502_chemotype_master.tsv, 2,502 biological BGCs x 91 fields, from analysis/viral_gt_reference_chemotype_curation_20261027/. Its frozen SHA-256 is 5f86b03ff3e071eec4fcafa02eac36c981f1561777a1d1acbe6f69e51315cbc9. The dictionary, Methods, reference QC and complete prior validation JSON were inspected. All 10,036 prior checks were programmatically parsed and PASS status confirmed; repeated per-accession check families are summarized in PRIOR_VALIDATION_REVIEW.json rather than treating truncated console output as a complete review.

Every consumed curation artifact was checked against its sealed output manifest, whose own SHA-256 was verified. Existing source decisions were used only to verify unchanged primary chemotypes and confidence. The original global assignment grid and original reference class labels were checked against the pinned archival checksum file, tied to private archival commit 284d4f1655cf036c1890bdffcb4aa2dc93f6e41e. Baseline: k20_r1_s42, original-320D k=20, resolution=1.0, seed=42. Existing MIBiG release identity is MIBiG 3.1; no current live release or new literature source was substituted.

## Populations, categories and denominators

Saccharide selection is membership in the original broad class-token list, including hybrids. Exactly 195 unique accessions qualify: 100 GC08, 24 GC05, 71 in neither GC08 nor GC05. These sets are pairwise disjoint and their union is the selected 195. The earlier 95 non-GC08 references include the 24 GC05 references and are not used as a separate population. Non-Saccharide references, including the five non-Saccharide GC08 members and 29 non-Saccharide GC05 members, are excluded from these denominators.

Existing mutually exclusive primary_chemotype values are mapped through an explicit 15-to-7 lookup, without fuzzy matching or reinterpretation. The exact source spellings, including LPS/LOS abbreviations, are retained in CHEMOTYPE_GROUP_MAPPING.tsv. All 91 original master fields are preserved in ACCESSION_LEVEL_CHEMOTYPE_MAPPING.tsv; four additional fields record the display population, consolidated chemotype, category order and explicit membership in the operational four-category set. OUTSIDE_DEFINED_SET means outside those named bins, not absence of surface chemistry. Unresolved primary chemotypes have an explicit UNRESOLVED_PRIMARY_CHEMOTYPE membership state.

The operational surface/extracellular definition is categories 1-4 only. The two GC08 phosphonoglycans remain in Other/unresolved; the source's broader 74/100 flag is neither substituted nor overwritten. Glycopeptidolipid and other miscellaneous classes also stay in the requested final bin. Unresolved records remain in denominators and in that bin, without conversion to biological negatives.

The 21 population-category cells were reconstructed from individual accessions. Historical numbers in the request were used solely as validation targets; any discrepancy would stop preparation before plotting. All expected values matched. Bar widths are n/N x 100, using each population's Saccharide BGC count as N. Exact rational accumulation verifies 100% for every row before conversion to display floats. FIGURE_SOURCE_DATA.tsv retains counts, fractions, percentages, cumulative boundaries and exact accession lists for every cell, including zeros.

"""+counts+"""

## Graphical design

One 100% stacked horizontal bar chart, 183 x 95 mm, with bars ordered GC08, GC05, Other GCs. Category order is identical across rows. The first four consecutive segments form the operational surface/extracellular set, with dark gold, the previous UMAP Saccharide gold, orange and light amber. Blue, purple and gray complete the stack. No gradients, decorative borders, explanatory in-artwork paragraphs or numerical table are used.

Arial: x-axis title and group labels 8.5 pt, percentage ticks 8 pt, legend 8 pt. The axis is fixed at 0-100%, with labeled 0%, 25%, 50%, 75%, 100% ticks and the title Saccharide reference BGCs (%). White 0.35 pt separators preserve boundaries without inflating segment widths. Bar height is 0.56 category units, leaving clear gaps. A two-column legend below the chart keeps the four warm categories together and uses unwrapped short labels.

Zero-count S-layer/other-surface segments in GC05 and Other GCs remain zero-width, with explicit zero source rows; no minimum-width inflation is applied. The smallest positive segment, 1/71, is 1.40845% or approximately 2.04 mm across the 145 mm plotting width. All exact values remain in the source table rather than labels inside narrow bars.

"""+palette+"""

## Palette and export verification

The initial palette was inspected under severity-100 Machado protanopia, deuteranopia and tritanopia matrices in linear sRGB plus grayscale. Blue/purple separation was weak under protanopia, so purple was lightened to #B68BC4 and capsule gold darkened to #6F4E0B. The final palette retains #C8A44C from the previous UMAP. Exact matrices and pairwise CIE76 diagnostics are preserved. The final simulations were visually checked. Diagnostic differences do not certify accessibility for every observer, and grayscale cannot uniquely encode all seven categories; fixed order and segment boundaries remain helpful.

Exports use a fixed canvas, 600-dpi PNG, SVG text nodes and PDF font type 42 with embedded Arial TrueType text. No raster images occur in the vector files. The final PDF is independently rendered at 300 dpi through Poppler, with the actual command and input/output hashes recorded. Both 150-dpi and 96-dpi previews preserve the intended layout. Actual SVG rectangle widths are independently checked against source percentages, including zeros.

## Annotation limitations

The selected cohort contains 168 SOURCE_ANNOTATION_SUPPORTED, 24 LITERATURE_SUPPORTED, one PROVISIONAL_RECORD_LEVEL and two UNRESOLVED annotations. None of these status names establishes uniform independent human manual verification.

- BGC0000262 (Other GCs): aglycone products and hybrid/tailoring context leave product scope unresolved.
- BGC0000727 (GC08): MIBiG indigoidine product versus an auricin/deoxysugar-related sequence description remains an unresolved source conflict.
- BGC0001765 (GC08): the pinned salecan record is strain S33 while cited experimental chemistry concerns ZX09; the EPS category is retained provisionally without transferring experimental claims.

The original detailed notes, statuses, source URLs and record hashes remain intact. ANNOTATION_REVIEW_FLAGS.tsv lists these three records. Broader curation limitations remain in the original curation reports. This task did not attempt to resolve them.

## Actual execution evidence

"""+table(['Execution folder','Start UTC','Actual child argument array'],commands)+"""

run_stage.py freezes code/configuration and records pre-execution command arrays, Python/platform, input-manifest hash, stdout, stderr and final exit status. The initial plotting attempt produced the artwork but failed to serialize NumPy tick integers into its layout receipt; that implementation error was corrected by converting ticks to native floats, and its unsuccessful log/code are retained. The final palette adjustment has its own preparation and plotting records. No counts or chemistry were changed during these display iterations.

INPUT_MANIFEST_SHA256.tsv freezes all consumed sources before calculations; every source is rehashed after verification. SOFTWARE_VERSIONS.json records actual plotting libraries and resolved Arial. VALIDATION_RESULTS.json contains the independent source, mathematical, geometry, typography and vector-export checks. Output hashes are sealed after all execution files and final visual review records exist. The manifest and its own verification receipt are excluded from recursive hashing.

All files remain local. Any later private archival push requires explicit scoped approval and PRIVATE visibility verification; the public repository is unchanged.
""")
 md('DATA_DICTIONARY.md',"""# Figure source-data dictionary

ACCESSION_LEVEL_CHEMOTYPE_MAPPING.tsv contains one row per selected Saccharide reference (195 unique MIBiG accessions), preserving all 91 source-master fields exactly. Source detailed chemistry, original multilabel model classes, evidence states, curation notes, source release, URLs and checksums remain unchanged.

Four display fields are appended:
- figure_population: GC08, GC05 or Other GCs. The last group excludes both GC08 and GC05.
- display_chemotype_group: the requested consolidated mutually exclusive seven-category display value, not a new primary annotation.
- display_category_order: 1-7, fixed across bars.
- four_category_surface_membership: IN_DEFINED_SET for categories 1-4, OUTSIDE_DEFINED_SET for other named categories, UNRESOLVED_PRIMARY_CHEMOTYPE for Unresolved. These are operational set memberships, not new biological presence/absence claims.

FIGURE_SOURCE_DATA.tsv contains 21 population x display-category rows, including zeros. n_BGCs is the accession count; denominator_saccharide_BGCs is 100, 24 or 71. fraction and percentage use that denominator. left_percentage and right_percentage give the cumulative boundaries for the fixed stack order. accessions_json lists exactly which biological BGCs contribute to that cell. A hybrid reference remains one BGC, not one observation per original class label.

CHEMOTYPE_GROUP_MAPPING.tsv gives all 15 original detailed category spellings and their seven display bins. It includes Teichuronic acid even though no selected Saccharide reference has that primary category. No unknown or missing primary category was silently mapped.

SURFACE_EXTRACELLULAR_SUMMARY.tsv uses only the first four display categories. Original broader source flags, including phosphonoglycan-inclusive ones, are preserved separately in the accession mapping. ANNOTATION_REVIEW_FLAGS.tsv retains the two unresolved and one provisional selected records and exact source notes.

PLOTTED_SEGMENT_RECEIPTS.tsv records actual Matplotlib rectangle geometry for all 21 cells. COLOR_PALETTE.tsv records fixed category order, short legend labels and final HEX colors. PALETTE_DISTINGUISHABILITY.tsv records display-only CIE76 pair diagnostics for normal, three CVD simulations and grayscale, without an accessibility certification threshold.
""")
 md('DATA_VALIDATION_REPORT.md',"""# Data validation report

PASS. The complete machine-readable results are in VALIDATION_RESULTS.json.

"""+table(['Check','Status','Evidence'],[[r['check'],r['status'],json.dumps(r['evidence'],ensure_ascii=False)] for r in CHECKS])+"""

## Verified composition

"""+counts+"""

No discrepancy from the supplied historical expectations was found. The three populations are disjoint and exhaustive, with Saccharide-only denominators. All 65 hybrid records and all unresolved/provisional source states are retained. All 91 source fields are unchanged at accession level. Every BGC has exactly one display category. All three stacks total exactly 100% under rational accounting.

The final figure uses existing annotations as proposed/source-based categories. This validation establishes provenance and display consistency, not experimental confirmation of those chemical assignments. No new annotation, statistical test or model was performed.
""")
 md('VISUAL_QC_REPORT.md',"""# Visual QC: Supplementary Figure 2

The final normal preview, three CVD simulations and grayscale were inspected at the intended 183 x 95 mm layout, with a 96-dpi screen-size preview as an additional readability check. The final vector PDF is independently rendered at 300 dpi; its reviewed file hash is recorded in VISUAL_QA.json.

- GC08, GC05 and Other GCs are immediately comparable on the same 0-100% scale, with clear nonoverlapping population labels.
- The first four consecutive warm-colored segments occupy 72% of GC08, 37.5% of GC05 and 4.225% of Other GCs. Their grouping is visually apparent without in-bar text.
- The Other-GC bar is dominated by the purple glycosylated-specialized-metabolite bin (51/71), contrasting with GC08's warm segments.
- The dark, medium, orange and light-amber segments are separable in the final normal and simulated views. Purple was lightened after the initial protanopia view showed weak blue/purple separation. Final CIE76 diagnostics are retained without claiming universal perceptual thresholds.
- All positive narrow segments remain visible at true widths. The minimum is approximately 2.04 mm (1/71). Zero-width categories remain zero, rather than being drawn as false positives.
- The two-column legend is readable without line wrapping: the four operational surface/extracellular bins form the first column, with the remaining three bins in the second.
- Arial 8-8.5 pt, a white background, thin axis and restrained solid colors match the minimalist style of the previous supplementary UMAP. No decorative frame, unnecessary grid, title block, GC-specific prose or numerical table appears in the artwork.
- Normal, protanopia, deuteranopia and tritanopia views preserve useful category separation. Grayscale loses some seven-way category distinction; the stable order and white boundaries remain visible. Color reproduction is preferable for category identification.
- PDF text remains embedded editable Arial TrueType; SVG text stays as text. Both vector outputs contain no raster images. Measured text boxes fit the canvas without overlaps.

The chart supports a descriptive difference in the proposed chemotype distributions, with GC08 dominated by the defined surface/extracellular set. It does not demonstrate independently validated enrichment, experimentally confirmed chemical families, shared viral products or evolutionary mechanisms. Two unresolved and one provisional source assignment remain in the plotted cohort.
""")
 md('README.md',"""# Supplementary Figure 2: Saccharide reference chemotypes

A 100% stacked horizontal bar chart using 195 original Saccharide-annotated references, including hybrids, in three disjoint populations: GC08 (100), GC05 (24), Other GCs (71). Counts are reconstructed from the frozen 2,502 x 91 chemotype master, not from hard-coded plotting values.

## Figure and caption

- [Vector PDF](figures/SuppFig2_saccharide_reference_chemotype_distribution.pdf)
- [Editable SVG](figures/SuppFig2_saccharide_reference_chemotype_distribution.svg)
- [600-dpi PNG](figures/SuppFig2_saccharide_reference_chemotype_distribution.png)
- [Review preview](qc/preview_183mm_150dpi.png)
- [Independent PDF render](qc/PDF_render_300dpi.png)
- [Supplementary figure legend](SUPPLEMENTARY_FIGURE_LEGEND.md)

Canvas 183 x 95 mm. Arial: axis title/group labels 8.5 pt, ticks/legend 8 pt. Seven fixed colors; EPS uses the prior UMAP Saccharide gold. The four-category surface/extracellular proportions are 72/100, 9/24 and 3/71. This definition excludes phosphonoglycans.

## Data, code and validation

- [Figure source table](FIGURE_SOURCE_DATA.tsv)
- [Full accession-level mapping with all original 91 fields](ACCESSION_LEVEL_CHEMOTYPE_MAPPING.tsv)
- [Explicit category consolidation](CHEMOTYPE_GROUP_MAPPING.tsv)
- [Surface/extracellular summary](SURFACE_EXTRACELLULAR_SUMMARY.tsv)
- [Final HEX palette](COLOR_PALETTE.tsv)
- [Unresolved/provisional annotation records](ANNOTATION_REVIEW_FLAGS.tsv)
- [Data dictionary](DATA_DICTIONARY.md)
- [Methods](METHODS_RECORD.md)
- [Data validation report](DATA_VALIDATION_REPORT.md)
- [Visual QC report](VISUAL_QC_REPORT.md)
- [Machine-readable validation](VALIDATION_RESULTS.json)

prepare_figure_inputs.py verifies frozen source hashes and existing assignments, constructs disjoint populations and consolidates existing primary chemotypes. Any count discrepancy stops preparation. plot_chemotype_bars.py renders FIGURE_SOURCE_DATA.tsv using config.json. validate_figure.py independently checks accession membership, original annotations, category counts, percentages, actual SVG rectangle geometry, dimensions, fonts and final exports.

run_stage.py preserves actual commands, code/configuration snapshots, Python/platform, stdout/stderr and exit codes under executions/. INPUT_MANIFEST_SHA256.tsv freezes consumed source identities. SOFTWARE_VERSIONS.json, PLOTTED_SEGMENT_RECEIPTS.tsv, LAYOUT_RECEIPTS.json, PDF_RENDER_EXECUTION.json and VISUAL_QA.json preserve execution/rendering evidence. The source validation JSON's 10,036 checks are parsed and summarized in PRIOR_VALIDATION_REVIEW.json.

qc/ includes normal, CVD and grayscale previews. Palette iterations affect display only; original chemical classifications never change. finalize_package.py seals OUTPUT_MANIFEST_SHA256.tsv and CHECKSUM_VERIFICATION.json after all final execution records exist; the manifest and its own seal exclude themselves to avoid recursion. Exact commands and source paths are available in the recorded execution folders; another environment must adapt local paths while preserving source hashes.

No prior analysis or curation database was overwritten. No clustering, inference, chemical annotation or enrichment test was run. Nothing was committed or pushed. All outputs remain local pending review and explicit archival approval.
""")
if __name__=='__main__':main()
