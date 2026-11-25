"""Transparent record/primary-literature curation rules; no model predictions."""
import collections,csv,json,re,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent
CATEGORIES=['Capsular polysaccharide','Extracellular/matrix polysaccharide','Lipopolysaccharide (LPS)','Lipooligosaccharide (LOS)','O-antigen','S-layer glycan','Other cell-surface glycan','Glycopeptidolipid','Teichuronic acid','Phosphonoglycan','Aminoglycoside / aminocyclitol-related metabolite','Glycosylated specialized metabolite','Other carbohydrate-related product','Non-carbohydrate-associated product','Unresolved']
# These accession-level decisions are explicit and inspectable; they do not use GC labels.
# Each maps to the original pinned record and/or the retrieved publication.
OVERRIDES={
'BGC0000727':('Unresolved','MIBiG JSON product indigoidine conflicts with GenBank auricin deoxysugar-tailoring cluster; source identity/product attribution cannot be resolved here.','MIBiG3.1 JSON /cluster/compounds/0/compound; GBK DEFINITION; DOI:10.1007/s12223-010-0018-5'),
'BGC0000758':('Capsular polysaccharide','Product, GenBank definition and capsule genetics support capsule; retain source saccharide subclass exopolysaccharide as an overlapping broader annotation, not a replacement product label.','PMID:7961465; PMID:12057957'),
'BGC0000768':('Capsular polysaccharide','The cited study explicitly calls gellan a capsular polysaccharide; a second describes exopolysaccharide released to medium. Capsule takes primary precedence; extracellular tag is retained.','PMID:14767675; PMID:15741740'),
'BGC0000797':('Capsular polysaccharide','The cited study explicitly identifies a secreted sphingan capsule. Capsule takes primary precedence over extracellular polymer.','PMID:8626338'),
'BGC0000780':('Other cell-surface glycan','O&K product and cited paper cover distinct O-core and K antigen regions. Do not collapse this mixed annotation into pure O-antigen, LPS or capsule.','PMID:18557895'),
'BGC0000793':('Lipopolysaccharide (LPS)','Source product and GenBank definition identify core oligosaccharide; retain LPS-core secondary tag, without asserting a complete O-antigen or a LOS pathway.','MIBiG3.1 /cluster/compounds; GBK DEFINITION; DOI:10.1046/j.1365-2958.2002.03243.x'),
'BGC0000792':('Other carbohydrate-related product','The cited study establishes an antigenic polysaccharide, but does not establish the precise capsular, wall or extracellular product architecture.','PMID:9712783'),
'BGC0001410':('Capsular polysaccharide','Named PS A is explicitly capsular in GenBank and the cited study.','PMID:11401972'),
'BGC0001411':('Capsular polysaccharide','Named PS B is explicitly capsular. Paper documents historical reassignment of B/C loci; preserve this pinned record name.','PMID:11035722'),
'BGC0001412':('Capsular polysaccharide','GenBank names PS C; the cited study historically called its deletion product PS B. Retain frozen name, flag historical B/C naming ambiguity; both are capsular.','PMID:10377135; PMID:11035722'),
'BGC0001586':('Extracellular/matrix polysaccharide','The cited study identifies EPS273 as an exopolysaccharide and tests its biosynthesis locus by deletions.','PMID:28698510'),
'BGC0001380':('Extracellular/matrix polysaccharide','Salecan is explicitly a soluble glucan; GenBank identifies an exopolysaccharide locus.','PMID:27678113; GBK DEFINITION'),
'BGC0001765':('Extracellular/matrix polysaccharide','Frozen product is salecan but the cited study concerns ZX09 while this sequence is S33. Keep record-level product category provisionally; do not transfer experimentally established ZX09 chemistry or enzyme evidence to S33.','MIBiG3.1 /cluster/compounds; PMID:27678113; GBK DEFINITION'),
'BGC0001894':('Extracellular/matrix polysaccharide','The cited study explicitly describes a surface-associated EPS with structural analysis. Retain surface-associated secondary tag.','PMID:30286997'),
'BGC0002302':('Extracellular/matrix polysaccharide','The cited study establishes a secreted glucorhamnan. Chemical structure is experimentally supported, but BGC involvement is described as putative.','PMID:31182571'),
'BGC0001747':('Other carbohydrate-related product','Recorded product is a nonulosonic acid precursor for a cell-surface protein glycan, not the complete surface glycan or polymer. Preserve precursor-versus-assembled-glycan distinction.','PMID:27986835'),
'BGC0001846':('Non-carbohydrate-associated product','Valinomycin/montanastatin are cyclodepsipeptides in the cited study, despite the frozen Saccharide/NRP and hybrid-tailoring annotations. Keep both original annotations and flag this chemotype discordance.','DOI:10.1002/slct.201900922; Wiley ChemistrySelect 4(27) contents'),
'BGC0001510':('Non-carbohydrate-associated product','Final anisomycin is a pyrrolidine alkaloid; cryptic glycosylation during biosynthesis does not mean a glycosylated final product.','PMID:28373542'),
'BGC0000262':('Unresolved','Record lists aglycone products and a hybrid/tailoring saccharide annotation. Do not assume each listed product is glycosylated; product scope requires review.','MIBiG3.1 /cluster/compounds; /cluster/saccharide/subclass'),
'BGC0000953':('Glycosylated specialized metabolite','Primary study identifies amicetin as a disaccharide nucleoside antibiotic, not a canonical aminoglycoside. Retain original MIBiG aminoglycoside subclass as a source annotation.','PMID:22267658'),
'BGC0002020':('Glycosylated specialized metabolite','Photolose is a characterized glycoamino acid containing an anhydro-N-acetylglucosamine moiety, not a surface polysaccharide.','PMID:30694662'),
'BGC0000899':('Other carbohydrate-related product','Frozen product is standalone desosamine. Do not equate the sugar-biosynthesis segment with a complete glycosylated macrolide pathway.','MIBiG3.1 /cluster/compounds; PMID:9770448'),
'BGC0002405':('Other carbohydrate-related product','Plant medium-chain acylsugars are defensive specialized products; not bacterial surface or extracellular polymers.','PMID:32613943'),
'BGC0000855':('Non-carbohydrate-associated product','Ectoine is the recorded compatible-solute product; no glycan or glycosylated final product is documented.','MIBiG3.1 /cluster/compounds; DOI:10.1007/s00792-011-0396-x'),
'BGC0000916':('Unresolved','Molybdenum cofactor is recorded but its precise sugar-derived cofactor architecture is outside glycan subclass definitions; retain unresolved rather than calling it a biological carbohydrate negative.','PMID:9675851'),
'BGC0002656':('Unresolved','Named aromatic polyketides are retained; the inspected abstract does not resolve carbohydrate association of every congener. Broad Polyketide class is not sufficient for a negative call.','PMID:35384223'),
'BGC0002450':('Glycosylated specialized metabolite','Frozen sarasinoside product and GenBank glycosidation locus support a provisional glycosylated-metabolite annotation. The cited multi-omics study proposes candidate BGCs; experimental BGC-product validation is not established.','MIBiG3.1 /cluster/compounds; GBK DEFINITION; PMID:32092934')}
# Product-specific chemistry assignments supported by explicit retrieved records.
GLYCO_PMIDS={
'BGC0000025':'10449723','BGC0000193':'11016846','BGC0000265':'MIBiG3.1 structure and GT annotations',
'BGC0000277':'10658661','BGC0000278':'11094336','BGC0000632':'18503195',
'BGC0000813':'17396254','BGC0000814':'19756308','BGC0000825':'12617516',
'BGC0000826':'MIBiG3.1 structure; staurosporine identity','BGC0000827':'MIBiG3.1 structure; staurosporine identity',
'BGC0001076':'20039253','BGC0001525':'26235728','BGC0001526':'27680198',
'BGC0002173':'30977375','BGC0002295':'31816530','BGC0002407':'34404922',
'BGC0002462':'32350138','BGC0002573':'31057776','BGC0002652':'35263115',
'BGC0002722':'23788733','BGC0002733':'MIBiG3.1 SMILES; product sinapigladioside',
'BGC0002721':'MIBiG3.1 product explicitly hydroxynitrile glucosides',
'BGC0000798':'MIBiG3.1 chemical structure; dhurrin'}
NONCARB_PMIDS={'BGC0000316':'MIBiG3.1 capreomycin chemical structures','BGC0000897':'20416511',
'BGC0001925':'30643149','BGC0001929':'30398325','BGC0002004':'29567715','BGC0002099':'33321099',
'BGC0002163':'28802024','BGC0002410':'33789052','BGC0002423':'34378964',
'BGC0002583':'33727351','BGC0002631':'34946566','BGC0002698':'35027450'}
AMINO=['2-deoxystreptamine','hydroxystreptomycin','acarbose','apramycin','butirosin','fortimicin','gentamicin','hygromycin a','istamycin','jinggangmycin','kanamycin','kasugamycin','lividomycin','neomycin','paromomycin','ribostamycin','sisomicin','actinospectacin','spectinomycin','tobramycin','validamycin','streptomycin','acarviostatin']
EPS_NAMES=['diutan polysaccharide','emulsan','stewartan','gellan polysaccharide','sphingan polysaccharide','colanic acid','xanthan','galactoglucan','succinoglycan','salecan','alginate']

# Follow-up source review resolved additional products using retrieved primary
# literature or the explicitly recorded chemical structure; no new structure
# or sugar stereochemistry was generated.
OVERRIDES.update({
'BGC0000419':('Glycosylated specialized metabolite','Ristomycin is explicitly described as a type III glycopeptide in the cited primary study.','PMID:25114137'),
'BGC0001072':('Glycosylated specialized metabolite','The frozen chemical structure shows a glycosidic sugar-bearing specialized scaffold; the record documents D-olivose GT specificity. No new stereochemical assignment is made.','MIBiG3.1 /cluster/compounds/0/chem_struct; /cluster/saccharide/glycosyltransferases'),
'BGC0001074':('Glycosylated specialized metabolite','The frozen cosmomycin D structure explicitly contains glycosidic sugar chains. This is gross structure curation, not a new structure prediction.','MIBiG3.1 /cluster/compounds/0/chem_struct'),
'BGC0002054':('Glycosylated specialized metabolite','Amipurimycin is explicitly a C9-sugar nucleoside antibiotic in its primary biosynthetic study.','PMID:30990701'),
'BGC0002349':('Glycosylated specialized metabolite','The primary study explicitly characterizes glycosylated epemicin macrolides and links the BGC by heterologous expression and inactivation.','PMID:34279911'),
'BGC0002350':('Glycosylated specialized metabolite','The primary study explicitly discusses bleomycin disaccharide modification and deoxy-bleomycin production.','PMID:28256838'),
'BGC0002361':('Glycosylated specialized metabolite','The primary structure study identifies amicetose/aculose sugars on baikalomycin angucyclines.','PMID:32392775'),
'BGC0002430':('Glycosylated specialized metabolite','Gausemycins are explicitly lipoglycopeptides with tyrosine arabinose glycosylation in the primary study.','PMID:34009717'),
'BGC0002431':('Glycosylated specialized metabolite','The frozen maduramicin structure explicitly shows a sugar attached to its polyether scaffold; retained as source-structure curation without assigning new stereochemistry.','MIBiG3.1 /cluster/compounds/0/chem_struct'),
'BGC0002479':('Glycosylated specialized metabolite','The primary study identifies peptidyl nucleoside antibiotics with a branched C9 sugar core.','PMID:32167762')
})

def main():
 cards=json.loads((ROOT/'PRIORITY_EVIDENCE_CARDS.json').read_text(encoding='utf-8'));pubs=json.loads((ROOT/'PUBMED_PRIMARY_RECORDS.json').read_text(encoding='utf-8'));out=[]
 for card in cards:
  acc=card['accession'];j=json.loads((ROOT/'raw_sources/mibig_json_3.1'/(acc+'.json')).read_text(encoding='utf-8'))['cluster']
  products=card['products'];text='; '.join(products).lower();definition=card['source_definition'].lower()
  subclass=j.get('saccharide',{}).get('subclass','')
  category='Unresolved';method='RULE_BASED_SOURCE_SCREEN';confidence='UNRESOLVED';reason='Inspected source does not resolve the requested chemical subclass; no class-only imputation.'
  source='MIBiG3.1 /cluster/compounds; /cluster/saccharide'
  if 'capsular polysaccharide' in text or 'capsular polysaccharide' in definition:category=CATEGORIES[0]
  elif 'lipo-oligosaccharide' in text or 'lipooligosaccharide' in text:category=CATEGORIES[3]
  elif 'lipopolysaccharide' in text:category=CATEGORIES[2]
  elif 'o&k-antigen' in text:category=CATEGORIES[6]
  elif 'o-antigen' in text:category=CATEGORIES[4]
  elif 's-layer glycan' in text:category=CATEGORIES[5]
  elif 'glycopeptidolipid' in text:category=CATEGORIES[7]
  elif 'teichuronic acid' in text:category=CATEGORIES[8]
  elif 'phosphonoglycan' in text:category=CATEGORIES[9]
  elif 'exopolysaccharide' in text or 'exopolysacch' in definition or any(n==text for n in EPS_NAMES):category=CATEGORIES[1]
  elif any(n in text for n in AMINO):category=CATEGORIES[10]
  elif subclass=='hybrid/tailoring':category=CATEGORIES[11]
  elif subclass=='oligosaccharide':category=CATEGORIES[12]
  elif acc in GLYCO_PMIDS:category=CATEGORIES[11];source='PMID:'+GLYCO_PMIDS[acc] if GLYCO_PMIDS[acc].isdigit() else GLYCO_PMIDS[acc]
  elif acc in NONCARB_PMIDS:category=CATEGORIES[13];source='PMID:'+NONCARB_PMIDS[acc] if NONCARB_PMIDS[acc].isdigit() else NONCARB_PMIDS[acc]
  if category!='Unresolved':
   confidence='SOURCE_ANNOTATION_SUPPORTED';reason='Exact source product/architecture or documented saccharide subclass supports this category; no viral-function interpretation.'
  if acc in GLYCO_PMIDS or acc in NONCARB_PMIDS:method='EVIDENCE_GUIDED_SOURCE_CURATION'
  if acc in OVERRIDES:
   category,reason,source=OVERRIDES[acc];method='EVIDENCE_GUIDED_AMBIGUITY_REVIEW'
   confidence='LITERATURE_SUPPORTED' if 'PMID:' in source or 'Wiley' in source else 'SOURCE_ANNOTATION_SUPPORTED'
   if category=='Unresolved':confidence='UNRESOLVED'
  tags=[];surface='NOT_DETERMINED';extracellular='NOT_DETERMINED'
  if category in [CATEGORIES[i] for i in (0,2,3,4,5,6,8)]:surface='True'
  if category==CATEGORIES[1] or category==CATEGORIES[9]:extracellular='True'
  if category in [CATEGORIES[10],CATEGORIES[11],CATEGORIES[13]]:
   surface='False';extracellular='False'
  if category==CATEGORIES[0]:tags.append('capsular_polymer')
  if acc in ['BGC0000768','BGC0000797']:extracellular='True';tags.append('capsular_and_secreted_polymer')
  if acc=='BGC0001894':surface='True';tags.append('surface_associated_EPS')
  if acc=='BGC0001747':surface='True';tags+=['surface_protein_glycan_precursor','nonulosonic_acid']
  if acc=='BGC0000769':
   # Product chemistry alone is insufficient for a specific cellular location.
   surface='NOT_DETERMINED';tags.append('mycobacterial_glycopeptidolipid')
  if acc=='BGC0000793':tags.append('LPS_core_oligosaccharide')
  if acc=='BGC0000780':tags+=['O_antigen_annotation','K_antigen_annotation','mixed_or_ambiguous_antigen_scope']
  if acc=='BGC0001510':tags+=['cryptic_intermediate_glycosylation','non_glycosylated_final_product']
  if acc=='BGC0001846':tags.append('original_Saccharide_label_chemistry_discordance')
  if acc in ['BGC0001765','BGC0002450']:confidence='PROVISIONAL_RECORD_LEVEL';tags.append('BGC_product_experimental_link_requires_review')
  if acc=='BGC0001412':tags.append('historical_PS_B_PS_C_name_conflict')
  if acc=='BGC0000763':tags.append('cited_DOI_names_amylovoran_not_recorded_stewartan')
  if acc=='BGC0000727':tags.append('JSON_GBK_product_conflict')
  if acc=='BGC0000953':tags.append('source_aminoglycoside_subclass_retained_but_nucleoside_primary_category')
  if j.get('saccharide',{}).get('subclass'):tags.append('source_saccharide_subclass:'+subclass)
  if any('aglycon' in name.lower() for name in products) or acc in ['BGC0000285','BGC0002173','BGC0000898']:tags.append('mixed_glycosylated_aglycone_or_precursor_products')
  combined='True' if 'True' in (surface,extracellular) else 'False' if surface==extracellular=='False' else 'NOT_DETERMINED'
  broader='True' if category in [CATEGORIES[2],CATEGORIES[3],CATEGORIES[4]] or acc=='BGC0000780' else 'False'
  review='True' if category=='Unresolved' or confidence=='PROVISIONAL_RECORD_LEVEL' or acc in ['BGC0000780','BGC0001412','BGC0000763','BGC0000953'] else 'False'
  out.append(dict(mibig_accession=acc,original_product_designation=products,primary_chemotype=category,secondary_chemotype_tags=tags,
   bacterial_surface_glycan_related=surface,bacterial_extracellular_glycan_related=extracellular,
   surface_or_extracellular_glycan_combined=combined,LPS_LOS_O_antigen_combined=broader,
   chemotype_annotation_confidence=confidence,chemotype_evidence_source=source,
   curation_method=method,curated_chemotype_status='REVIEWED_UNRESOLVED' if category=='Unresolved' else method,
   manual_review_required=review,curation_notes=reason))
 assert len(out)==len({r['mibig_accession'] for r in out})==229
 (ROOT/'CURATION_DECISIONS.json').write_text(json.dumps(dict(categories=CATEGORIES,decisions=out,rule_version='2026-10-27.3',manual_human_review=False),ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
 print('Priority decisions:',len(out),'unresolved:',sum(r['primary_chemotype']=='Unresolved' for r in out))
 for card,r in zip(cards,out):
  if r['primary_chemotype']=='Unresolved':print('UNRESOLVED',r['mibig_accession'],card['products'],r['curation_notes'])
 for group in ['GC08','GC05','NON_GC08']:
  chosen=[r for card,r in zip(cards,out) if ('Saccharide' in card['original_label'].split()) and (card['GC']==group if group.startswith('GC') else card['GC']!='GC08')]
  print(group,len(chosen),dict(collections.Counter(r['primary_chemotype'] for r in chosen)), 'surface/extracellular',sum(r['surface_or_extracellular_glycan_combined']=='True' for r in chosen))
if __name__=='__main__':
 sys.stdout.reconfigure(encoding='utf-8');main()
