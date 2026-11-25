"""Validate exact source joins and generate scientific curation documentation."""
from curation_common import *
import math,re,tarfile
def main():
 checks=[]
 def check(name,ok,details=None):
  checks.append(dict(check=name,status='PASS' if ok else 'FAIL',details=details))
  if not ok:
   write_json('validation_results.json',dict(status='FAIL',checks=checks))
   raise ValueError(name+': '+str(details))
 master=read_table(OUT/'reference_2502_chemotype_master.tsv');byacc=unique(master,'mibig_accession')
 meta=read_table(ARCHIVE/'inputs/workspace/analysis/viral_glycan_benchmark_20261016/mibig_matched_reference_metadata.tsv');rm=unique(meta,'BGC')
 check('2502_unique_original_biological_references',len(master)==len(byacc)==len(rm)==2502 and set(byacc)==set(rm))
 check('original_embedding_row_correspondence', [r['mibig_accession'] for r in master]==[r['BGC'] for r in meta] and all(r['original_embedding_row']==rm[r['mibig_accession']]['embedding_row'] for r in master))
 check('2502_unique_embedding_rows',len({r['original_embedding_row'] for r in master})==2502)
 check('original_functional_labels_unchanged',all(r['original_bgcprophet_class']==rm[r['mibig_accession']]['label'] and json.loads(r['original_class_multilabel'])==rm[r['mibig_accession']]['label'].split() for r in master))
 check('reference_core_gene_counts_unchanged',all(r['reference_core_gene_count']==rm[r['mibig_accession']]['core_gene_count'] for r in master))
 gm=unique(read_table(ARCHIVE/'assignments/global_joint_leiden_grid.csv'),'sample_id')
 check('single_existing_global_baseline_assignment',all(r['global_leiden_community']==gm[r['mibig_accession']]['k20_r1_s42'] for r in master))
 sac=[r for r in master if r['is_saccharide_reference']=='True'];g8=[r for r in master if r['is_GC08']=='True'];g5=[r for r in master if r['is_GC05']=='True']
 g8s=[r for r in g8 if r['is_saccharide_reference']=='True'];g5s=[r for r in g5 if r['is_saccharide_reference']=='True'];outside=[r for r in sac if r['is_GC08']=='False']
 priority=[r for r in master if r['is_priority_curation_reference']=='True']
 counts=dict(all_references=len(master),all_saccharide=len(sac),GC08_all=len(g8),GC08_saccharide=len(g8s),GC05_all=len(g5),GC05_saccharide=len(g5s),saccharide_outside_GC08=len(outside),priority=len(priority))
 check('expected_reference_subsets',counts==dict(all_references=2502,all_saccharide=195,GC08_all=105,GC08_saccharide=100,GC05_all=53,GC05_saccharide=24,saccharide_outside_GC08=95,priority=229),counts)
 check('overlap_GC05_saccharide_subset_outside_GC08',{r['mibig_accession'] for r in g5s}<={r['mibig_accession'] for r in outside})
 check('all_priority_have_explicit_curation_status',all(r['curated_chemotype_status'] not in ['',NR,'NOT_REVIEWED'] for r in priority))
 check('background_not_reviewed_not_negative',all(r['primary_chemotype']=='NOT_REVIEWED' and r['surface_or_extracellular_glycan_combined']=='NOT_REVIEWED' for r in master if r['is_priority_curation_reference']=='False'))
 check('all_priority_categories_traceable',all(r['chemotype_evidence_source'] not in ['',NR,'NOT_REVIEWED'] for r in priority))
 check('predictions_excluded_from_experimental_evidence_types',all(not any(x in ['Sequence-based prediction','Structure-based inference','Homology','Other'] for x in (json.loads(r['experimental_evidence_types']) if r['experimental_evidence_types']!=NR else [])) for r in master))
 check('EPS_surface_localization_not_assumed',all(r['bacterial_surface_glycan_related']=='NOT_DETERMINED' for r in priority if r['primary_chemotype']=='Extracellular/matrix polysaccharide' and r['mibig_accession']!='BGC0001894'))

 mibig=load_mibig();compounds=read_table(OUT/'reference_compounds.tsv');genes=read_table(OUT/'reference_biosynthetic_genes.tsv')
 check('compound_relationship_keys_unique',len(unique(compounds,'compound_relationship_id'))==len(compounds))
 check('gene_record_keys_unique',len(unique(genes,'gene_record_id'))==len(genes))
 observations=read_table(OUT/'reference_biosynthetic_gene_observations.tsv')
 check('class_specific_gene_observation_parent_integrity',all(r['mibig_accession'] in byacc and r['source_gene_identifier'] and r['source_record_SHA256']==byacc[r['mibig_accession']]['source_record_SHA256'] for r in observations))
 cg=collections.defaultdict(list);gg=collections.defaultdict(list)
 for r in compounds:cg[r['mibig_accession']].append(r)
 for r in genes:gg[r['mibig_accession']].append(r)
 for acc,row in byacc.items():
  raw=mibig[acc]['raw']['cluster'];source_compounds=raw.get('compounds',[])
  check('compound_names_and_payload:'+acc,len(cg[acc])==len(source_compounds) and [json.loads(r['source_compound_annotations_json']) for r in cg[acc]]==source_compounds)
  check('master_compound_names:'+acc,json.loads(row['compound_names_original'])==[c.get('compound','') for c in source_compounds])
  source_ids={str(g['id']) for g in raw.get('genes',{}).get('annotations',[])}
  represented=set()
  for g in gg[acc]:
   if g['MIBiG_function_annotations']!=NR:represented.update(str(a['id']) for a in json.loads(g['MIBiG_function_annotations']))
  check('all_explicit_JSON_gene_ids_preserved:'+acc,source_ids<=represented)
  check('raw_JSON_SHA256_match:'+acc,row['source_record_SHA256']==mibig[acc]['sha256'])
  # Direct assay flags must derive from gene-level evidence, never BGC product identity.
  for g in gg[acc]:
   if g['experimental_enzyme_activity_recorded']=='True':
    ev=json.loads(g['function_evidence_types'])
    assert any('activity assay' in e.lower() or 'enzymatic assay' in e.lower() for e in ev)
 # Independent count directly from raw GBK tar, not from the parser's normalized output.
 raw_cds=0;raw_member_hashes={}
 with tarfile.open(OUT/'raw_sources/mibig_gbk_3.1.tar.gz','r:gz') as tar:
  for m in tar:
   if m.isfile() and m.name.endswith('.gbk'):
    acc=Path(m.name).stem;b=tar.extractfile(m).read()
    raw_cds+=len(re.findall(rb'^     CDS\s',b,re.M))
    raw_member_hashes[acc]=hashlib.sha256(b).hexdigest()
 check('raw_GBK_CDS_count_preserved',sum(int(r['number_of_annotated_genes']) for r in master)==raw_cds,dict(raw_GBK_CDS=raw_cds,normalized_gene_records=len(genes)))
 check('all_gene_source_record_hashes_preserved',all(g['source_GBK_record_SHA256']==raw_member_hashes[g['mibig_accession']] for g in genes))
 v=unique(read_table(VIRAL),'viral_locus_id');links=read_table(OUT/'viral_287_reference_chemotype_links.tsv');nn=read_table(ARCHIVE/'tables/selected_reference_neighbors_top20.csv')
 check('exact_287_viral_loci_and_5740_Top20_links',len(v)==287 and len(links)==len(nn)==5740 and {r['viral_locus_id'] for r in links}==set(v))
 lp=collections.Counter((r['viral_locus_id'],r['reference_rank']) for r in links)
 check('one_each_saved_Top20_rank_per_locus',len(lp)==5740 and all(n==1 for n in lp.values()) and all({r['reference_rank'] for r in links if r['viral_locus_id']==k}=={str(i) for i in range(1,21)} for k in v))
 check('all_neighbor_ids_ranks_and_cosines_unchanged',all(a['viral_locus_id']==b['viral_locus_id'] and a['reference_mibig_accession']==b['reference_id'] and a['reference_rank']==b['rank'] and a['original_320D_cosine_similarity']==b['cosine_320D'] and a['original_reference_bgcprophet_class']==b['reference_class'] for a,b in zip(links,nn)))
 check('viral_GC_VC_SC_and_P10_unchanged',all(r['viral_only_community']==v[r['viral_locus_id']]['viral_only_leiden_cluster'] and r['global_joint_community']==v[r['viral_locus_id']]['global_joint_leiden_cluster'] and r['saccharide_focused_community']==v[r['viral_locus_id']]['saccharide_joint_leiden_cluster'] and r['original_representative_P_at_10']==v[r['viral_locus_id']]['P_at_10'] for r in links))
 check('viral_master_exact_prior_integrated_hash',sha(VIRAL)=='c3b4d6c8967051d5a1ec488137f777096e3d7c49779154fa9f97c59208060e4d')
 check('all_reference_links_resolve_once',all(r['reference_mibig_accession'] in byacc for r in links))
 preflight=read_json(OUT/'SOURCE_PREFLIGHT.json')
 check('private_archival_inputs_unchanged_after_curation',all(sha(Path(r['path']))==r['sha256'] for r in preflight['source_files']))
 check('frozen_json_snapshot_unchanged',sha(OUT/'raw_sources/mibig_json_3.1.tar.gz')==preflight['json_sha256'])
 check('frozen_gbk_snapshot_unchanged',sha(OUT/'raw_sources/mibig_gbk_3.1.tar.gz')==preflight['gbk_sha256'])
 # Appendix accounting and tables include all denominators, including unresolved.
 cats=read_json(OUT/'CURATION_DECISIONS.json')['categories'];groups={'GC08_Saccharide':g8s,'GC05_Saccharide':g5s,'Non_GC08_Saccharide':outside,'GC08_all_references':g8,'GC05_all_references':g5}
 composition=[]
 for cat in cats:
  row={'primary_chemotype':cat}
  for group,rs in groups.items():
   n=sum(r['primary_chemotype']==cat for r in rs);row[group+'_n']=n;row[group+'_percent']=format(100*n/len(rs),'.6f')
  composition.append(row)
 check('mutually_exclusive_primary_category_counts',all(sum(int(r[group+'_n']) for r in composition)==len(rs) for group,rs in groups.items()))
 write_table('chemotype_group_comparison.tsv',composition)
 legacy_classes=['Capsular polysaccharide','Extracellular/matrix polysaccharide','S-layer glycan']
 legacy=sum(r['primary_chemotype'] in legacy_classes or r['LPS_LOS_O_antigen_combined']=='True' for r in g8s)
 broad=sum(r['surface_or_extracellular_glycan_combined']=='True' for r in g8s)
 historical=dict(target_capsular=29,target_EPS=21,target_LPS_LOS_O=19,target_S_layer=3,target_combined=72,
  source_reconstruction_capsular=sum(r['primary_chemotype']=='Capsular polysaccharide' for r in g8s),
  source_reconstruction_EPS=sum(r['primary_chemotype']=='Extracellular/matrix polysaccharide' for r in g8s),
  source_reconstruction_LPS_LOS_O=sum(r['LPS_LOS_O_antigen_combined']=='True' for r in g8s),
  source_reconstruction_S_layer=sum(r['primary_chemotype']=='S-layer glycan' for r in g8s),
  legacy_four_group_union=legacy,inclusive_surface_extracellular_union=broad,
  original_accession_level_manual_ledger='NOT_AVAILABLE_IN_SUPPLIED_SOURCES; do not claim recovery of original individual decisions')
 write_json('GC08_HISTORICAL_COUNT_AUDIT.json',historical)
 selected=[dict(mibig_accession=r['mibig_accession'],recorded_product=r['compound_names_original'],primary_chemotype=r['primary_chemotype'],
  legacy_four_group_included=str(r['primary_chemotype'] in legacy_classes or r['LPS_LOS_O_antigen_combined']=='True'),
  inclusive_surface_extracellular=r['surface_or_extracellular_glycan_combined'],evidence=r['chemotype_evidence_source'],notes=r['curation_notes']) for r in g8s]
 write_table('GC08_accession_level_audit.tsv',selected)
 accounting=dict(mibig_matches=sum(r['source_record_SHA256']!='NOT_AVAILABLE' for r in master),priority=len(priority),
  priority_assigned=sum(r['primary_chemotype']!='Unresolved' for r in priority),
  saccharide_assigned=sum(r['primary_chemotype']!='Unresolved' for r in sac),
  priority_unresolved=sum(r['primary_chemotype']=='Unresolved' for r in priority),
  background_NOT_REVIEWED=len(master)-len(priority),
  priority_provisional=sum(r['chemotype_annotation_confidence']=='PROVISIONAL_RECORD_LEVEL' for r in priority),
  hybrid_reference_records=sum(r['original_is_hybrid']=='True' for r in master),
  BGCs_with_recorded_structures=sum(r['chemical_structure_availability']!=NR for r in master),
  BGCs_with_explicit_MIBiG_gene_annotations=sum(int(r['number_of_explicit_MIBiG_gene_annotations'])>0 for r in master),
  raw_GBK_CDS=raw_cds,normalized_gene_records=len(genes),
  genes_with_explicit_activity_assay=sum(r['experimental_enzyme_activity_recorded']=='True' for r in genes),
  gene_mapping_ambiguous=sum('ambiguous_GBK_alias' in r['source_gene_mapping_status'] for r in genes),
  compound_relationships=len(compounds),viral_links=len(links),
  original_vs_mibig_class_disagreements=sum(r['original_vs_mibig_class_disagreement']=='True' for r in master))
 result=dict(status='PASS',checks=checks,failed_checks=0,subgroup_counts=counts,dataset_accounting=accounting,historical_GC08=historical,scientific_models_executed=False,enrichment_tests=False,commits=False,pushes=False)
 write_json('validation_results.json',result)
 # Data dictionary: each master field has a definition and source.
 defs=read_json(OUT/'FIELD_DEFINITIONS.json')
 lines=['# Chemotype data dictionary','','The master has '+str(len(master[0]))+' fields and exactly 2,502 original biological BGC rows. UTF-8 TSV; source lists/structured annotations use JSON cells. Original labels and original embedding-row strings are retained. Missingness and review states are explicit.','',
  '## Missingness and review codes','','- NOT_REPORTED: the inspected source/curated evidence does not supply this fact; not an absence claim.',
  '- NOT_REVIEWED: detailed chemical curation has not been performed; never a biological negative.',
  '- NOT_DETERMINED / Unresolved: evidence does not support a precise conclusion or sources conflict.',
  '- True/False surface/extracellular flags describe the documented product scope only, not the absence of genes or all possible functions.',
  '- Human manual verification was not performed. Methods distinguish rule-based source annotation from evidence-guided record/literature curation.','',
  '## Master fields','','| Field | Definition | Source |','| --- | --- | --- |']
 for field in master[0]:lines.append('| '+field+' | '+defs[field]['definition'].replace('|','/')+' | '+defs[field]['source'].replace('|','/')+' |')
 lines+=['','## Supporting tables','',
  '- reference_compounds.tsv: one BGC–source-compound relationship per row, keyed by compound_relationship_id. Original compound JSON and all names/structures are retained. The parent BGC chemotype does not assign a chemical category to every congener.',
  '- reference_biosynthetic_genes.tsv: one source-documented gene entity per row, keyed by gene_record_id. GBK CDS rows use source ordinal; explicit unmatched MIBiG gene IDs remain separate. Mapping ambiguity is explicit. Protein, locus and gene aliases are preserved. JSON function arrays and GBK qualifiers retain exact documented functions.',
  '- Gene indexed_function_categories are transparent lexical indexes of existing annotations, not new sequence/domain predictions. Empty categories mean unclassified/not reported. Genes may appear in several master function lists.',
  '- Only explicit gene-level activity-assay evidence sets experimental_enzyme_activity_recorded=True. Genetic evidence, sequence predictions and unspecified source-function annotations retain separate levels.',
  '- reference_biosynthetic_gene_observations.tsv: normalized gene references in class-specific source sections, with exact JSON pointers and lossless annotation fragments. Module/domain or multi-gene evidence is not silently treated as whole-gene catalytic validation.',
  '- chemotype_curation_audit.tsv: all 2,502 records, including 229 priority records with explicit source/evidence/rationale/status and 2,273 NOT_REVIEWED background records.',
  '- reference_chemical_features.tsv: one accession–field–source claim. Values are JSON lists with evidence_level, source, claim_scope and notes. No chemistry is inferred from broad BGC classes.',
  '- viral_287_reference_chemotype_links.tsv: 5,740 original Top-20 associations. original_320D_cosine_similarity is copied verbatim; community labels and P@10 come unchanged from the existing viral master. These associations are not viral product predictions.',
  '- chemotype_group_comparison.tsv: descriptive category counts and percentages; GC05 Saccharide is a subset of non-GC08 Saccharide, not an independent comparison sample.',
  '- GC08_accession_level_audit.tsv: source-based inclusion and notes for every GC08 Saccharide reference.','']
 (OUT/'CHEMOTYPE_DATA_DICTIONARY.md').write_text('\n'.join(lines),encoding='utf-8')
 methods=['# Curation methods and evidence limits','',
  '## Source snapshot and immutable identity','',
  'MIBiG 3.1 JSON and GenBank archives already existed in the workspace and were reused without combining releases. JSON official URL: https://dl.secondarymetabolites.org/mibig/mibig_json_3.1.tar.gz . Original retrieval: 2026-09-08T05:00:37.926106+00:00; SHA-256 '+preflight['json_sha256']+'. GenBank official URL: https://dl.secondarymetabolites.org/mibig/mibig_gbk_3.1.tar.gz ; SHA-256 '+preflight['gbk_sha256']+'. Its original download time was not recovered from the inspected manifest; this limitation is not replaced by filesystem timestamps. The reuse timestamp is recorded in SOURCE_PREFLIGHT.json. Raw tarballs and all raw JSON records are retained.','',
  'The exact private archive commit is '+COMMIT+'. Remote Git blob identities were checked against local consumed clustering inputs. The 2,502 original reference accessions, zero-based embedding rows, labels, and core gene counts were preserved; the NPY header/whole-file hash were inspected for identity only. Existing k20_r1_s42 global assignments were selected by exact column name. The 287-locus master and 5,740 saved neighbor identities/ranks/cosines were joined by exact identifiers. No models, nearest-neighbor searches, clustering or UMAP fitting were run.','',
  '## Scope and annotation rules','',
  'The priority population is the 229-record unique union of all 195 original Saccharide references with all GC08 and GC05 references. Detailed chemotype decisions are independent of community labels and use the same rules throughout. All other 2,273 records receive automatic source extraction and explicit NOT_REVIEWED chemical status. Primary categories are mutually exclusive; secondary tags may overlap. Hybrid original labels are preserved as strings and JSON token lists.','',
  'The exact rules, accession-level overrides and explanations are in curate_priority_records.py and CURATION_DECISIONS.json. Obvious product terms and structured source subclasses receive RULE_BASED_SOURCE_SCREEN status; record/structure/literature interpretation receives explicit evidence-guided status. No human manual verification is claimed. No category is chosen to satisfy the 72/100 historical target. A capsule takes primary precedence over a broad extracellular-polysaccharide label when capsule is specifically supported; both descriptions may remain in secondary scope flags. Mixed O/K annotation is retained as Other cell-surface glycan rather than forced into a single O-antigen/LOS/CPS identity.','',
  'Gross glycosylated architecture can be curated from a recorded chemical structure. No structure was generated, no molecular mass recalculated, and no sugar stereochemistry, linkage, precursor or tailoring chemistry was invented from broad biosynthetic labels. Normalized names use Unicode NFC and whitespace only, preserving original names independently. Compound relationships are kept separate so congeners/aglycones do not multiply the BGC master.','',
  '## Evidence hierarchy','',
  '1. Explicit biochemical or chemical-structure experiments in the primary study support only the tested enzyme/compound claim.',
  '2. Genetic or heterologous-expression evidence supports the tested BGC/gene involvement, not every catalytic function.',
  '3. MIBiG product/structure/subclass annotations support record-level chemistry and are retained separately from direct experiment claims.',
  '4. Recorded source genomic predictions, GBK product annotations and MIBiG sequence-based predictions remain predicted/unspecified evidence.',
  '5. Missing, conflicting or unreviewed information is explicit. Product identity never makes all enzymes experimentally characterized.','',
  'PubMed metadata/abstracts for all 473 cited PMID records of the priority set were retrieved separately on 2026-10-27. Original source publication IDs are not replaced; supplemental DOI metadata remains a separate field. Exact request URLs, raw XML and hashes are retained. This is record/abstract-guided curation, not exhaustive full-text experimental replication. Supplemental publisher/author source notes are in SUPPLEMENTARY_SOURCE_NOTES.json.','',
  'JSON gene annotations and GT specificities are mapped to GBK CDS by unique exact protein/gene/locus aliases. Unmatched or ambiguous aliases remain separately documented, never merged by a fuzzy guess. The standard-library GenBank parser preserves qualifier text/locations and raw archive hashes; the normalized table omits protein translation text because lossless source archives retain it. Existing antiSMASH-related qualifiers in the frozen GBK are source predictions, not newly run analyses.','',
  'Chemical claims are accession-specific. Experimental salecan data for ZX09 are not copied to the S33 entry. Phosphonate-containing polysaccharides are not called phosphorylated glycans. A sugar precursor or cryptic glycosylation step is distinguished from a final glycan. Genus/species columns are explicitly lexical components of the recorded source organism; no new taxonomy retrieval was performed.','',
  '## Historical GC08 audit','',
  'The original accession-level preliminary manual ledger was not supplied or found in the inspected sources. Therefore this is an independent accession-level source reconstruction, not a claimed recovery of every old manual decision. The legacy four-group union is compared with the broader surface/extracellular union separately; overlapping definitions are not independent observations. No enrichment p-values or viral chemistry predictions were calculated.','',
  '## Source rights and archive destination','',
  'MIBiG identifies CC-BY-4.0 attribution licensing: https://mibig.secondarymetabolites.org/ . Cite the MIBiG contributors and source release. Retrieved publication abstracts may have separate third-party rights; raw literature files are local research-source records and redistribution requires checking applicable terms. No Git commit/push is authorized by completion of this task. Proposed private destination: REDACTED_PRIVATE_REPOSITORY/results/WP07_reference_chemotype_curation_20261027/ .','']
 (OUT/'CURATION_METHODS.md').write_text('\n'.join(methods),encoding='utf-8')
 # Source manifest also covers archives, retrieval receipts, curation inputs and raw XML.
 sm=read_table(OUT/'SOURCE_MANIFEST_SHA256.tsv')
 existing={r['path'] for r in sm}
 pmmanifest=read_json(OUT/'PUBMED_RETRIEVAL_MANIFEST.json')
 raw_url={str(OUT/b['path']):b for b in pmmanifest['batches']}
 extra=[OUT/'raw_sources/mibig_json_3.1.tar.gz',OUT/'raw_sources/mibig_gbk_3.1.tar.gz',OUT/'raw_sources/original_download_manifest.json',
  OUT/'PRIVATE_ARCHIVE_IDENTITY.json',OUT/'CURATED_CHEMICAL_FEATURES.json',OUT/'CURATION_DECISIONS.json',
  OUT/'SUPPLEMENTARY_SOURCE_NOTES.json',OUT/'PUBMED_RETRIEVAL_MANIFEST.json',OUT/'PUBMED_PRIMARY_RECORDS.json']+list((OUT/'raw_sources/pubmed').glob('*.xml'))
 for p in extra:
  b=raw_url.get(str(p),{})
  sm.append(dict(role='snapshot_or_evidence_input',path=str(p),sha256=sha(p),size_bytes=p.stat().st_size,source_url=b.get('url','OFFICIAL_FROZEN_SOURCE_OR_LOCAL_CURATED_RECORD; see provenance'),source_version='MIBiG 3.1' if p.name.startswith('mibig_') else 'curation_2026-10-27',retrieval_date=b.get('retrieval_utc','SEE_SOURCE_PREFLIGHT_AND_ORIGINAL_MANIFEST'),mibig_accession='',record_status='',annotation_completeness='',supporting_publications=''))
 write_table('SOURCE_MANIFEST_SHA256.tsv',sm)
 # Concise scientific report with exact denominators and evidence qualifications.
 lines=['# Reference chemotype QC and scientific summary','','**PASS**: identity, joins, original labels/rows, priority statuses, normalized tables and saved viral crosswalk verified.','',
  '## Dataset accounting','']+['- '+k+': '+str(v) for k,v in accounting.items()]+['',
  'All 229 priority entries have an explicit curation status; '+str(accounting['priority_assigned'])+' have a proposed category and '+str(accounting['priority_unresolved'])+' remain unresolved. '+str(accounting['saccharide_assigned'])+'/195 Saccharide references have a proposed category. NOT_REVIEWED background entries are not biological negatives. Provisional and unresolved entries require further review.','',
  '## Chemical composition','',
  'GC05 Saccharide references (n=24) are a subset of non-GC08 Saccharide references (n=95). These columns overlap and must not be treated as independent biological samples.','',
  '| Primary category | GC08 Saccharide (100) | GC05 Saccharide (24) | Non-GC08 Saccharide (95) |','| --- | ---: | ---: | ---: |']
 for r in composition:
  lines.append('| '+r['primary_chemotype']+' | '+' | '.join(str(r[k+'_n'])+' ('+format(float(r[k+'_percent']),'.1f')+'%)' for k in ['GC08_Saccharide','GC05_Saccharide','Non_GC08_Saccharide'])+' |')
 lines+=['','### All reference members, including non-Saccharide exceptions','','| Primary category | GC08 all (105) | GC05 all (53) |','| --- | ---: | ---: |']
 for r in composition:lines.append('| '+r['primary_chemotype']+' | '+str(r['GC08_all_references_n'])+' | '+str(r['GC05_all_references_n'])+' |')
 lines+=['','## GC08 historical audit','',
  'The four legacy product groups yield '+str(legacy)+'/100. Source-specific capsule precedence yields '+str(historical['source_reconstruction_capsular'])+' CPS and '+str(historical['source_reconstruction_EPS'])+' EPS rather than the earlier 29/21 split. BGC0000768 (gellan) and BGC0000797 (S88 sphingan) are explicitly capsular in their cited papers and can also be described as secreted EPS; moving those two to an EPS-first convention would give 29/21 without changing the 72-record legacy union. This numerical reconciliation does not establish the unavailable original accession-level decisions.','',
  'The broader source-supported surface/extracellular union is '+str(broad)+'/100. The additional BGC0000806 and BGC0000807 are chemically characterized extracellular phosphonoglycans (PMID 24584498). They were outside the earlier four named bins. BGC0000769 glycopeptidolipid has an explicit chemical category but the inspected source did not establish a specific surface location; its surface flag remains NOT_DETERMINED.','',
  '## Ambiguities and annotation limitations','',
  '- BGC0000727: indigoidine in JSON versus an auricin deoxysugar-tailoring sequence description. Retained as Unresolved; original Saccharide label unchanged.',
  '- BGC0001846: original Saccharide/NRP annotation versus documented valinomycin/montanastatin cyclodepsipeptides; chemistry and model label remain separate.',
  '- BGC0001510: cryptic glycosylation in anisomycin biosynthesis; non-glycosylated final product retained.',
  '- BGC0001747: pseudaminic-acid precursor for a cell-surface protein glycan, not a complete surface polymer.',
  '- BGC0000780: mixed O/K antigen scope; BGC0000793: LPS core oligosaccharide, not automatically a complete O-antigen or LOS pathway.',
  '- BGC0000763 cites an amylovoran paper while the pinned record describes stewartan; BGC0001412 retains historical PS B/C naming ambiguity.',
  '- BGC0001765: cited salecan chemistry concerns ZX09 rather than the S33 sequence. Only provisional record-level chemotype is retained.',
  '- BGC0002450: sarasinoside/glycosidation candidate locus lacks demonstrated BGC-product linkage in the inspected source; provisional.',
  '- Product identity does not demonstrate every enzyme. Gene assay flags are restricted to explicit gene-level source assay evidence. Unknown chemistry/enzymes are not coded as absent.','',
  '## Reproducibility and scope','',
  'See CHEMOTYPE_DATA_DICTIONARY.md, CURATION_METHODS.md, validation_results.json, SOURCE_MANIFEST_SHA256.tsv, OUTPUT_MANIFEST_SHA256.tsv, execution logs and executed code snapshots. All models/clustering/nearest-neighbor calculations and enrichment tests were excluded. No previous archive or viral master was modified, and no commits/pushes were performed. Proposed private destination is recorded in CURATION_METHODS.md; explicit approval is required before archival push.','']
 (OUT/'REFERENCE_QC.md').write_text('\n'.join(lines),encoding='utf-8')
 print(json.dumps(dict(status='PASS',checks=len(checks),counts=counts,accounting=accounting,historical=historical),indent=2),flush=True)
if __name__=='__main__':
 sys.stdout.reconfigure(encoding='utf-8');main()
