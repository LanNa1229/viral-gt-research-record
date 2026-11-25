"""Build source-backed MIBiG chemical/biochemical tables; no scientific models."""
from curation_common import *
import re
DEFS={}
def setfield(row,name,value,description,source):
 row[name]=as_cell(value)
 DEFS.setdefault(name,dict(definition=description,source=source))
def gene_roles(product,ann,gts,sugar_context=False):
 # Source-function indexing, not sequence/domain prediction.
 text=product.lower();roles=[]
 if gts or re.search(r'(glycosyl|glucosyl|mannosyl|rhamnosyl|galactosyl|xylosyl|fucosyl|arabinosyl|hexosyl)transferase',text):
  roles.append('glycosyltransferase')
 if re.search(r'\b(?:[dD]?[tT][dD][pP]|[uU][dD][pP]|[gG][dD][pP]|[cC][dD][pP]|NDP|nucleotide)[ -]?(?:sugar|gluc|rham|mann|galact|fuc|hex|amin)',product) or 'nucleotide sugar' in text:
  roles.append('nucleotide_sugar_biosynthesis')
 sugar= sugar_context or any(w in text for w in ['sugar','glucose','rhamnose','mannose','galactose','fucose','hexose','nucleotide','deoxysugar'])
 if sugar and any(w in text for w in ['epimerase','dehydratase','reductase','aminotransferase','methyltransferase','acetyltransferase','deacetylase','sulfotransferase','kinase']):
  roles.append('sugar_modification')
 if any(w in text for w in ['polysaccharide polymerase','polysaccharide biosynthesis','glycan assembly','polysaccharide assembly','repeat unit polymerase','o-antigen polymerase']):
  roles.append('glycan_polymerization_or_assembly')
 if any(w in text for w in ['polysaccharide export','polysaccharide transport','glycan export','glycan transport','flippase','o-antigen transport','capsule export']):
  roles.append('glycan_export_or_transport')
 if not roles and (ann or product and product.lower() not in ['hypothetical protein','uncharacterized protein','putative protein']):
  roles.append('other_documented_gene_annotation')
 return roles
def extract_genes(acc,c,gbk,gbksha,jsonsha):
 cds=[f for f in gbk['features'] if f['type']=='CDS']
 annotations=c.get('genes',{}).get('annotations',[])
 gtrows=c.get('saccharide',{}).get('glycosyltransferases',[])
 contexts={str(v) for group in c.get('saccharide',{}).get('sugar_subclusters',[]) for v in group}
 entities=[];aliasmap=collections.defaultdict(list)
 for i,f in enumerate(cds,1):
  q=f['qualifiers'];aliases=set(q.get('protein_id',[])+q.get('locus_tag',[])+q.get('gene',[]))
  entity=dict(q=q,location=f['location'],annotations=[],gts=[],key=f'GBK_CDS_{i:05d}',mapping=[])
  entities.append(entity)
  for alias in aliases:aliasmap[alias].append(entity)
 unresolved={}
 def entity_for(identifier):
  matches=aliasmap.get(identifier,[])
  if len(matches)==1:return matches[0],'exact_source_identifier'
  if identifier not in unresolved:
   unresolved[identifier]=dict(q={},location=NR,annotations=[],gts=[],key='MIBIG_ID:'+identifier,mapping=[])
  return unresolved[identifier],('ambiguous_GBK_alias_kept_separate' if len(matches)>1 else 'MIBiG_gene_not_matched_to_GBK_CDS')
 for a in annotations:
  ent,status=entity_for(str(a.get('id','')));ent['annotations'].append(a);ent['mapping'].append(status)
 for a in gtrows:
  ent,status=entity_for(str(a.get('gene_id','')));ent['gts'].append(a);ent['mapping'].append(status)
 entities+=list(unresolved.values())
 rows=[]
 for e in entities:
  q=e['q'];anns=e['annotations'];gts=e['gts']
  ids=list(dict.fromkeys([str(a.get('id','')) for a in anns]+[str(a.get('gene_id','')) for a in gts]))
  pgbk='; '.join(q.get('product',[]));pjson='; '.join(str(a.get('product','')) for a in anns if a.get('product'))
  ev=recursive_evidence(anns+gts)
  level='DIRECT_ACTIVITY_ASSAY_RECORDED' if any('activity assay' in a.lower() or 'enzymatic assay' in a.lower() for a in ev) else 'GENETIC_OR_EXPRESSION_EVIDENCE_RECORDED' if any(w in ' '.join(ev).lower() for w in ['knock-out','heterologous expression','complementation']) else 'SOURCE_PREDICTED_FUNCTION' if any('prediction' in a.lower() or 'inference' in a.lower() for a in ev) else 'SOURCE_FUNCTION_ANNOTATION_EVIDENCE_UNSPECIFIED'
  prod=pgbk+' '+pjson
  roles=gene_roles(prod,anns,gts,any(k in contexts for k in ids+q.get('gene',[])+q.get('protein_id',[])))
  qual={k:v for k,v in q.items() if k!='translation'}
  rows.append(dict(mibig_accession=acc,gene_record_id=acc+':'+e['key'],gene_id=as_cell(ids or q.get('protein_id') or q.get('locus_tag') or q.get('gene') or [e['key']]),
   protein_id=as_cell(q.get('protein_id',[])),locus_tag=as_cell(q.get('locus_tag',[])),gene_name=as_cell(q.get('gene',[]) or [a['name'] for a in anns if a.get('name')]),
   location=as_cell(e['location']),documented_product_GBK=as_cell(pgbk),documented_product_MIBiG=as_cell(pjson),
   indexed_function_categories=as_cell(roles),MIBiG_function_annotations=as_cell(anns),MIBiG_glycosyltransferase_specificities=as_cell(gts),
   function_evidence_types=as_cell(ev),function_evidence_level=level,experimental_enzyme_activity_recorded=str(level=='DIRECT_ACTIVITY_ASSAY_RECORDED'),
   source_gene_mapping_status=as_cell(sorted(set(e['mapping'])) or ['GBK_CDS_no_explicit_MIBiG_gene_annotation']),
   raw_GBK_qualifiers_except_translation=jtext(qual),source_GBK_record_SHA256=gbksha,source_JSON_record_SHA256=jsonsha))
 return rows,len(cds)
def main():
 data=load_mibig();gbks=read_json(OUT/'GBK_PARSED_ANNOTATIONS.json')
 gbkhash={r['mibig_accession']:r['sha256'] for r in read_table(OUT/'GBK_RECORD_MANIFEST.tsv')}
 meta=read_table(ARCHIVE/'inputs/workspace/analysis/viral_glycan_benchmark_20261016/mibig_matched_reference_metadata.tsv')
 global_map=unique(read_table(ARCHIVE/'assignments/global_joint_leiden_grid.csv'),'sample_id')
 decisions={r['mibig_accession']:r for r in read_json(OUT/'CURATION_DECISIONS.json')['decisions']}
 features=read_json(OUT/'CURATED_CHEMICAL_FEATURES.json');features_by=collections.defaultdict(list)
 for f in features:features_by[f['mibig_accession']].append(f)
 pubs=read_json(OUT/'PUBMED_PRIMARY_RECORDS.json')
 master=[];compounds=[];genes=[];audit=[];chemical_features=[];source_manifest=[];gene_observations=[]
 archived=read_json(OUT/'SOURCE_PREFLIGHT.json')['source_files']
 for r in archived:source_manifest.append(dict(role='verified_private_archive_input',path=r['path'],sha256=r['sha256'],size_bytes=r['size_bytes'],source_url='REDACTED_PRIVATE_REPOSITORY'+COMMIT,source_version=COMMIT,retrieval_date='2026-10-27',mibig_accession='',record_status='',annotation_completeness='',supporting_publications=''))
 source_manifest.append(dict(role='existing_viral_master_input',path=str(VIRAL),sha256=sha(VIRAL),size_bytes=VIRAL.stat().st_size,source_url='LOCAL_VERIFIED_INTEGRATION',source_version='2026-10-27',retrieval_date='NOT_APPLICABLE',mibig_accession='',record_status='',annotation_completeness='',supporting_publications=''))
 for m in meta:
  acc=m['BGC'];record=data.get(acc);row={};p=OUT/'raw_sources/mibig_json_3.1'/(acc+'.json')
  c=record['raw']['cluster'] if record else {}
  bg=gbks.get(acc,{});clist=c.get('compounds',[]);decision=decisions.get(acc);gc=global_map[acc]['k20_r1_s42'];is_sac='Saccharide' in m['label'].split()
  def put(name,value,description,source):setfield(row,name,value,description,source)
  put('mibig_accession',acc,'Original biological reference accession; unique primary key.','frozen CORE_ONLY metadata:BGC')
  put('original_embedding_row',m['embedding_row'],'Original zero-based CORE_ONLY row, retained exactly; augmentation rows are not observations.','frozen metadata:embedding_row')
  put('original_bgcprophet_class',m['label'],'Original label string used in clustering; never replaced by source chemotypes.','frozen metadata:label')
  put('original_class_multilabel',m['label'].split(),'Original class tokens in original order as a JSON list; overlapping hybrid classes retained.','frozen metadata:label')
  put('original_is_hybrid',len(m['label'].split())>1,'More than one original class token.','original label tokens')
  put('global_leiden_community',gc,'Existing global joint baseline, k20_r1_s42; no recomputation.','global_joint_leiden_grid.csv:k20_r1_s42')
  for name,value,desc in [('is_saccharide_reference',is_sac,'Any Saccharide token in original labels, including hybrids.'),('is_GC08',gc=='GC08','GC08 membership of the reference.'),('is_GC05',gc=='GC05','GC05 membership of the reference.'),('is_saccharide_outside_GC08',is_sac and gc!='GC08','Original Saccharide reference outside GC08.')]:
   put(name,value,desc,'original labels and saved global assignment')
  put('reference_core_gene_count',m['core_gene_count'],'Frozen CORE_ONLY reference gene count; distinct from all annotated CDS count.','frozen metadata:core_gene_count')
  put('is_priority_curation_reference',acc in decisions,'Union of all Saccharide references and all GC08/GC05 references; 229 unique BGCs.','original labels and saved assignments')
  source_features=[f for f in bg.get('features',[]) if f['type']=='source'];sq=source_features[0]['qualifiers'] if source_features else {}
  gorg=(sq.get('organism') or [NR])[0];organism=c.get('organism_name',NR)
  taxids=[s.split(':',1)[1] for s in sq.get('db_xref',[]) if s.startswith('taxon:')]
  headerorg=bg.get('headers',{}).get('ORGANISM','')
  lineage=headerorg[len(gorg):].strip() if gorg!=NR and headerorg.startswith(gorg) else headerorg or NR
  tokens=organism.split();genus=NR;species=NR
  if tokens and re.fullmatch(r'\[?[A-Z][a-z]+\]?',tokens[0]) and tokens[0] not in ['Candidatus','Uncultured','Unknown']:
   genus=tokens[0]
   if len(tokens)>1 and re.fullmatch('[a-z][a-z-]+',tokens[1]) and tokens[1] not in ['sp','spp','bacterium','symbiont','strain']:species=' '.join(tokens[:2])
  for name,value,desc,source in [
   ('organism_name',organism,'Original MIBiG organism string; no current taxonomy substitution.','JSON:/cluster/organism_name'),
   ('ncbi_taxid',c.get('ncbi_tax_id'),'Recorded taxon ID kept as text.','JSON:/cluster/ncbi_tax_id'),
   ('genus',genus,'Lexical genus component of recorded organism name where unambiguous; not a new taxonomic determination.','recorded organism string; conservative lexical parsing'),
   ('species',species,'Recorded binomial component where unambiguous; strain text omitted without taxonomy lookup.','recorded organism string; conservative lexical parsing'),
   ('broader_taxonomy',lineage,'Documented GenBank lineage, retained separately from JSON organism.','frozen GBK:ORGANISM continuation'),
   ('taxonomy_provenance','FROZEN_JSON_AND_GBK_3.1; name components lexically parsed; no taxonomy retrieval','Distinguishes recorded taxonomy and derived name components from new retrievals.','curation method'),
   ('genbank_organism_name',gorg,'Frozen GenBank organism qualifier, including any disagreement with JSON.','GBK:source/organism'),
   ('genbank_ncbi_taxid',taxids,'Taxon cross-references in GenBank source qualifier.','GBK:source/db_xref'),
   ('source_record_accession',c.get('loci',{}).get('accession'),'Original nucleotide record accession/version, not the BGC accession.','JSON:/cluster/loci/accession'),
   ('bgc_start',c.get('loci',{}).get('start_coord'),'Original JSON BGC start coordinate; no coordinate-system conversion or inferred default.','JSON:/cluster/loci/start_coord'),
   ('bgc_end',c.get('loci',{}).get('end_coord'),'Original JSON BGC end coordinate; no inferred default.','JSON:/cluster/loci/end_coord'),
   ('bgc_completeness',c.get('loci',{}).get('completeness'),'Original locus completeness, retaining Unknown separately from missing.','JSON:/cluster/loci/completeness')]:put(name,value,desc,source)
  put('genbank_definition',bg.get('headers',{}).get('DEFINITION'),'Frozen sequence description, retained for source conflict audit.','GBK:DEFINITION')
  put('genbank_comment',bg.get('headers',{}).get('COMMENT'),'Original sequence warnings/comments; sequence replacement warnings are retained.','GBK:COMMENT')
  original_names=[x.get('compound','') for x in clist]
  normalized=[' '.join(unicodedata.normalize('NFC',x).split()) for x in original_names]
  put('compound_names_original',original_names,'All original compound names, ordered JSON list; no synonym merging.','JSON:/cluster/compounds/*/compound')
  put('normalized_product_names',normalized,'Unicode NFC plus whitespace normalization only; no chemical name reassignment.','original compound names')
  put('number_of_annotated_compounds',len(clist),'Number of source BGC-compound relationship records, including all congeners.','JSON:/cluster/compounds')
  for col,key,desc in [('chemical_formula','molecular_formula','Original formula per source compound.'),('molecular_mass','mol_mass','Original source mol_mass values; no mass computation or unit reinterpretation.'),('SMILES','chem_struct','Original source chemical structure strings; no structure generation.')]:
   put(col,[dict(compound_index=i,value=x[key]) for i,x in enumerate(clist,1) if x.get(key) is not None],desc,'JSON:/cluster/compounds/*/'+key)
  nstruct=sum(bool(x.get('chem_struct')) for x in clist)
  put('chemical_structure_availability','AVAILABLE_ALL' if clist and nstruct==len(clist) else 'AVAILABLE_PARTIAL' if nstruct else NR,'Availability of recorded structures, not new predictions.','source chem_struct presence')
  put('source_mibig_biosynthetic_classes',c.get('biosyn_class'),'MIBiG 3.1 classes retained independently of frozen training labels.','JSON:/cluster/biosyn_class')
  put('source_saccharide_subclass',c.get('saccharide',{}).get('subclass'),'Exact structured source Saccharide subclass; not replaced by curated primary category.','JSON:/cluster/saccharide/subclass')
  cat=decision['primary_chemotype'] if decision else 'NOT_REVIEWED'
  glycan_cats=['Capsular polysaccharide','Extracellular/matrix polysaccharide','Lipopolysaccharide (LPS)','Lipooligosaccharide (LOS)','O-antigen','S-layer glycan','Other cell-surface glycan','Phosphonoglycan','Teichuronic acid']
  architecture='GLYCAN_POLYMER_OR_ASSEMBLY' if cat in glycan_cats else 'GLYCOSYLATED_SPECIALIZED_METABOLITE' if cat=='Glycosylated specialized metabolite' else 'AMINOGLYCOSIDE_OR_AMINOCYCLITOL_METABOLITE' if cat.startswith('Aminoglycoside') else 'GLYCOPEPTIDOLIPID' if cat=='Glycopeptidolipid' else 'OTHER_CARBOHYDRATE_PRODUCT' if cat=='Other carbohydrate-related product' else 'NON_GLYCAN_FINAL_PRODUCT' if cat=='Non-carbohydrate-associated product' else cat
  glycan_kind='GLYCAN_OR_GLYCAN_ASSEMBLY' if cat in glycan_cats else 'GLYCOSYLATED_METABOLITE' if cat=='Glycosylated specialized metabolite' else 'AMINOGLYCOSIDE_OR_AMINOCYCLITOL' if cat.startswith('Aminoglycoside') else 'GLYCOLIPID' if cat=='Glycopeptidolipid' else 'CARBOHYDRATE_RELATED_PRODUCT_OR_PRECURSOR' if cat=='Other carbohydrate-related product' else 'NON_GLYCOSYLATED_FINAL_PRODUCT' if cat=='Non-carbohydrate-associated product' else cat
  polymer='POLYMERIC_OR_ASSEMBLED_GLYCAN' if cat in glycan_cats else 'NONPOLYMERIC' if cat in ['Glycosylated specialized metabolite','Non-carbohydrate-associated product','Glycopeptidolipid'] or cat.startswith('Aminoglycoside') else cat if cat in ['Unresolved','NOT_REVIEWED'] else 'NOT_DETERMINED'
  if acc=='BGC0001747':glycan_kind='SUGAR_PRECURSOR_FOR_SURFACE_GLYCAN';polymer='NONPOLYMERIC'
  if acc=='BGC0000899':glycan_kind='AMINO_SUGAR_PRECURSOR';polymer='NONPOLYMERIC'
  if acc=='BGC0001510':glycan_kind='NON_GLYCOSYLATED_FINAL_PRODUCT';architecture='PYRROLIDINE_WITH_CRYPTIC_INTERMEDIATE_GLYCOSYLATION'
  if acc=='BGC0000780':polymer='MIXED_O_AND_K_ANTIGEN_SCOPE'
  put('product_architecture',architecture,'Source-supported gross architecture; NOT_REVIEWED is not a negative call.','curation decisions')
  put('primary_chemotype',cat,'Mutually exclusive curated category; independent from original model classes.','CURATION_DECISIONS.json')
  put('secondary_chemotype_tags',decision['secondary_chemotype_tags'] if decision else 'NOT_REVIEWED','Overlapping descriptors and ambiguity flags; do not add category counts as independent observations.','curation decisions')
  put('product_glycan_or_glycosylated_metabolite',glycan_kind,'Distinguishes assembled glycan, glycometabolite, sugar precursor and non-glycosylated final product.','source-supported curation')
  put('polymeric_or_nonpolymeric',polymer,'Gross source-supported architecture, not a computed molecular descriptor.','source-supported curation')
  association='CRYPTIC_INTERMEDIATE_GLYCOSYLATION' if acc=='BGC0001510' else 'SURFACE_GLYCAN_PRECURSOR' if acc=='BGC0001747' else 'SOURCE_DOCUMENTED_GLYCAN_OR_GLYCOSYLATION' if cat in glycan_cats+['Glycosylated specialized metabolite','Glycopeptidolipid'] else 'AMINO_SUGAR_OR_AMINOCYCLITOL_ASSOCIATION' if cat.startswith('Aminoglycoside') else cat
  put('glycan_association_type',association,'Scope-specific relation; final-product chemistry and biosynthetic glycosylation are separate.','curation decisions and cited sources')
  surf=decision['bacterial_surface_glycan_related'] if decision else 'NOT_REVIEWED';ext=decision['bacterial_extracellular_glycan_related'] if decision else 'NOT_REVIEWED'
  target='BACTERIAL_CELL_SURFACE_AND_EXTRACELLULAR_POLYMER' if surf==ext=='True' else 'BACTERIAL_CELL_SURFACE' if surf=='True' else 'BACTERIAL_EXTRACELLULAR_GLYCAN' if ext=='True' else 'NOT_ESTABLISHED_FOR_SURFACE_GLYCANS' if decision else 'NOT_REVIEWED'
  put('biological_localization_or_target',target,'Only documented glycan location scope; not inferred from a generic biosynthetic label.','curation decisions')
  for col in ['bacterial_surface_glycan_related','bacterial_extracellular_glycan_related','surface_or_extracellular_glycan_combined','LPS_LOS_O_antigen_combined','chemotype_annotation_confidence','chemotype_evidence_source']:
   put(col,decision[col] if decision else 'NOT_REVIEWED','Explicit source-supported curation descriptor. True/False apply to documented product scope; NOT_DETERMINED and NOT_REVIEWED are not biological negatives.','CURATION_DECISIONS.json')
  chem_fields=['known_monosaccharide_identities','unusual_or_rare_monosaccharides','amino_sugars','deoxy_sugars','sugar_nucleotide_precursors','reported_glycosidic_linkages','branching_or_polymeric_architecture','methylation','acetylation','sulfation','phosphorylation','other_glycan_tailoring_modifications']
  for field in chem_fields:
   matched=[f for f in features_by[acc] if f['field']==field]
   put(field,[v for f in matched for v in f['value']],'Only explicit curated source claims; missing means not reported/reviewed in inspected sources, not chemical absence. See long evidence table for scope, level and source.','reference_chemical_features.tsv')
  put('chemical_composition_evidence_level',sorted({f['evidence_level'] for f in features_by[acc]}),'Chemical evidence levels distinguish primary experimental chemistry, literature annotation and review claims.','curated chemical features')
  put('sugar_specificity_annotations',c.get('saccharide',{}).get('glycosyltransferases'),'Recorded enzyme substrate specificities and evidence; not silently treated as demonstrated product composition.','JSON:/cluster/saccharide/glycosyltransferases')
  grows,ncds=extract_genes(acc,c,bg,gbkhash.get(acc,NR),record['sha256'] if record else NR);genes+=grows
  put('number_of_annotated_genes',ncds,'Number of CDS source annotations in the frozen GenBank BGC record; not the count of experimentally characterized enzymes.','frozen GBK:CDS features')
  put('number_of_explicit_MIBiG_gene_annotations',len(c.get('genes',{}).get('annotations',[])),'Number of explicitly curated JSON gene annotations; distinguish from all CDS.','JSON:/cluster/genes/annotations')
  categories={'glycosyltransferase_genes':'glycosyltransferase','nucleotide_sugar_biosynthesis_genes':'nucleotide_sugar_biosynthesis','sugar_modification_genes':'sugar_modification','glycan_polymerization_or_assembly_genes':'glycan_polymerization_or_assembly','glycan_export_or_transport_genes':'glycan_export_or_transport','other_biosynthetic_genes':'other_documented_gene_annotation'}
  for col,tag in categories.items():
   values=[dict(gene_record_id=g['gene_record_id'],gene_id=g['gene_id'],documented_product=g['documented_product_MIBiG'] if g['documented_product_MIBiG']!=NR else g['documented_product_GBK'],evidence_level=g['function_evidence_level']) for g in grows if g['indexed_function_categories']!=NR and tag in json.loads(g['indexed_function_categories'])]
   put(col,values,'Index of explicit source function descriptions. Roles may overlap; no sequence-based annotation and no absence claim.','reference_biosynthetic_genes.tsv')
  put('experimentally_characterized_enzymes',[dict(gene_record_id=g['gene_record_id'],gene_id=g['gene_id'],evidence=g['function_evidence_types']) for g in grows if g['experimental_enzyme_activity_recorded']=='True'],'Only genes with explicit source Activity assay evidence; knock-outs alone are not direct catalytic characterization.','MIBiG gene-level evidence')
  put('pathway_description',dict(recorded_products=original_names,source_saccharide_subclass=c.get('saccharide',{}).get('subclass'),source_description=bg.get('headers',{}).get('DEFINITION')),'Source descriptions only; no new pathway reconstruction.','frozen JSON and GBK')
  put('pathway_evidence_level',c.get('loci',{}).get('evidence'),'Exact locus/BGC evidence list; not transferred to every enzyme or product claim.','JSON:/cluster/loci/evidence')
  put('pathway_source_annotations_json',{k:c[k] for k in ['saccharide','polyketide','nrp','terpene','ripp','other'] if k in c},'Lossless source pathway-class annotations including predictions; not newly calculated chemistry.','MIBiG 3.1 class-specific sections')
  put('mibig_record_status',c.get('status'),'Original record status.','JSON:/cluster/status')
  put('mibig_minimal_record',c.get('minimal'),'Original minimal-record flag.','JSON:/cluster/minimal')
  all_evidence=recursive_evidence(c)
  experimental=[e for e in all_evidence if e not in ['Sequence-based prediction','Structure-based inference','Homology','Other']]
  put('experimental_evidence_types',experimental,'Explicit source experimental/genetic/chemical evidence labels, excluding prediction, homology and inference. Aggregate labels do not imply every enzyme is characterized.','explicit JSON evidence arrays, experimental terms only')
  put('source_evidence_types_all',all_evidence,'Complete original source evidence vocabulary, retaining predictions/inferences separately from experimental_evidence_types.','all explicit JSON evidence arrays')
  publications=c.get('publications',[]);pmids=[p[7:] for p in publications if p.startswith('pubmed:')];dois=[p[4:] for p in publications if p.startswith('doi:')]
  put('publication_PMIDs',pmids,'Original source PubMed IDs; no replacement by retrieved identifiers.','JSON:/cluster/publications')
  put('publication_DOIs',dois,'Original source DOI strings.','JSON:/cluster/publications')
  put('publication_identifiers_original',publications,'All original publication identifier strings, unchanged.','JSON:/cluster/publications')
  put('retrieved_publication_DOIs',sorted({d for pm in pmids if pm in pubs for d in pubs[pm]['doi']}),'Supplemental PubMed metadata, kept separate from original DOI annotations.','473-record PubMed retrieval')
  put('source_record_URL','https://dl.secondarymetabolites.org/mibig/mibig_json_3.1.tar.gz#'+record['member'] if record else 'NOT_AVAILABLE','Pinned official snapshot URL with archive member; current live records are not substituted.','frozen official MIBiG 3.1')
  put('source_record_SHA256',record['sha256'] if record else 'NOT_AVAILABLE','SHA-256 of exact raw JSON record bytes.','raw snapshot member')
  put('source_version','MIBiG 3.1','Single release used across all reference annotations.','frozen source manifest')
  put('curated_chemotype_status',decision['curated_chemotype_status'] if decision else 'NOT_REVIEWED','Rule-based or evidence-guided status is explicit; no human manual verification is claimed.','curation decision method')
  put('manual_review_required',decision['manual_review_required'] if decision else 'True','Records needing further human/scientific review, including unreviewed background entries.','curation decisions')
  put('curation_notes',decision['curation_notes'] if decision else 'Automatic source extraction only; detailed chemotype curation NOT_REVIEWED.','Accession-specific reasoning and unresolved evidence limits.','curation decisions')
  class_conflict=set(c.get('biosyn_class',[]))!=set(m['label'].split())
  named_conflict=acc in ['BGC0000727','BGC0001412','BGC0000763','BGC0001765','BGC0001846','BGC0000953']
  organism_conflict=organism!=gorg and gorg!=NR
  put('original_vs_mibig_class_disagreement',class_conflict,'Exact original/source class token-set disagreement; both remain unchanged.','frozen metadata versus MIBiG 3.1')
  put('JSON_GBK_organism_disagreement',organism_conflict,'Exact source organism string disagreement; may reflect nomenclature or strain differences.','JSON organism versus GBK source organism')
  put('annotation_conflict_flag',class_conflict or named_conflict or organism_conflict,'Recorded source or chemical-scope conflict; inspect notes and retained source values.','source cross-checks')
  master.append(row)
  audit.append(dict(mibig_accession=acc,original_product_designation=jtext(original_names),proposed_curated_category=cat,
   exact_supporting_record_or_publication=row['source_record_URL']+'; '+(decision['chemotype_evidence_source'] if decision else 'NOT_REVIEWED'),
   source_record_SHA256=row['source_record_SHA256'],curation_confidence=row['chemotype_annotation_confidence'],
   reasoning_for_ambiguous_cases=row['curation_notes'],review_status=row['curated_chemotype_status'],
   curation_method=decision['curation_method'] if decision else 'AUTOMATED_EXTRACTION_ONLY',
   human_manual_verification='NOT_PERFORMED',manual_review_required=row['manual_review_required'],is_priority_reference=str(acc in decisions)))
  for i,x in enumerate(clist,1):
   compounds.append(dict(mibig_accession=acc,compound_relationship_id=f'{acc}:compound:{i}',compound_index=str(i),
    compound_name_original=x.get('compound',NR),normalized_product_name=normalized[i-1],chemical_formula=as_cell(x.get('molecular_formula')),
    molecular_mass=as_cell(x.get('mol_mass')),SMILES=as_cell(x.get('chem_struct')),source_compound_annotations_json=jtext(x),
    original_bgcprophet_class=m['label'],parent_curated_chemotype=cat,compound_specific_chemotype_status='NOT_INDIVIDUALLY_CLASSIFIED; parent BGC category is not asserted for every congener',
    source_record_URL=row['source_record_URL']+f'::/cluster/compounds/{i-1}',source_record_SHA256=row['source_record_SHA256'],supporting_publications=jtext(publications)))
  def collect_gene_observations(node,pointer):
   if isinstance(node,dict):
    explicit=[]
    for key in ['genes','gene_id']:
     v=node.get(key)
     if isinstance(v,str):explicit.append(v)
     elif isinstance(v,list):explicit.extend(x for x in v if isinstance(x,str))
    for gid in dict.fromkeys(explicit):
     gene_observations.append(dict(mibig_accession=acc,source_gene_identifier=gid,source_JSON_pointer=pointer,source_annotation_json=jtext(node),
      source_evidence_types=as_cell(recursive_evidence(node)),evidence_scope='Exact source node; multi-gene/module/domain evidence is not automatically assigned as whole-gene catalytic characterization',
      source_record_SHA256=row['source_record_SHA256']))
    for k,v in node.items():collect_gene_observations(v,pointer+'/'+str(k))
   elif isinstance(node,list):
    for i,v in enumerate(node):collect_gene_observations(v,pointer+'/'+str(i))
  for section in ['saccharide','polyketide','nrp','terpene','ripp','other']:
   if section in c:collect_gene_observations(c[section],'/cluster/'+section)
  source_manifest.append(dict(role='frozen_MIBiG_JSON_record',path=str(p),sha256=row['source_record_SHA256'],size_bytes=p.stat().st_size,
   source_url=row['source_record_URL'],source_version='MIBiG 3.1',retrieval_date='2026-09-08T05:00:37.926106+00:00',mibig_accession=acc,record_status=row['mibig_record_status'],
   annotation_completeness='minimal='+row['mibig_minimal_record']+'; loci='+row['bgc_completeness'],supporting_publications=jtext(publications)))
  for f in features_by[acc]:chemical_features.append({k:as_cell(v) for k,v in f.items()})
 # Join original Top-20 relationships to the separate reference annotation table.
 byacc=unique(master,'mibig_accession');viral=unique(read_table(VIRAL),'viral_locus_id')
 neighbors=read_table(ARCHIVE/'tables/selected_reference_neighbors_top20.csv');links=[]
 for n in neighbors:
  v=viral[n['viral_locus_id']];ref=byacc[n['reference_id']]
  assert n['Leiden_cluster']==v['viral_only_leiden_cluster']
  assert n['reference_class']==ref['original_bgcprophet_class']
  links.append(dict(viral_locus_id=n['viral_locus_id'],viral_only_community=v['viral_only_leiden_cluster'],global_joint_community=v['global_joint_leiden_cluster'],saccharide_focused_community=v['saccharide_joint_leiden_cluster'],
   reference_rank=n['rank'],reference_mibig_accession=n['reference_id'],original_320D_cosine_similarity=n['cosine_320D'],
   original_reference_bgcprophet_class=n['reference_class'],curated_reference_chemotype=ref['primary_chemotype'],reference_compound_names=ref['compound_names_original'],
   chemotype_evidence_level=ref['chemotype_annotation_confidence'],chemotype_review_status=ref['curated_chemotype_status'],
   original_representative_P_at_10=v['P_at_10'],association_interpretation='REPRESENTATION_ASSOCIATION_ONLY; no viral product/function assignment'))
 write_table('reference_2502_chemotype_master.tsv',master)
 write_table('reference_compounds.tsv',compounds);write_table('reference_biosynthetic_genes.tsv',genes)
 write_table('chemotype_curation_audit.tsv',audit);write_table('viral_287_reference_chemotype_links.tsv',links)
 write_table('reference_chemical_features.tsv',chemical_features)
 write_table('reference_biosynthetic_gene_observations.tsv',gene_observations)
 write_table('SOURCE_MANIFEST_SHA256.tsv',source_manifest)
 write_json('FIELD_DEFINITIONS.json',DEFS)
 print('Built',len(master),'references,',len(compounds),'compound relationships,',len(genes),'gene records,',len(links),'viral-reference links,',len(DEFS),'master columns.',flush=True)
if __name__=='__main__':
 sys.stdout.reconfigure(encoding='utf-8');main()
