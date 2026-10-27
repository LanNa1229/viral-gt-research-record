"""Accession-based reconciliation with conservative, explicit evidence decisions."""
import sys
sys.dont_write_bytecode=True
import re,json,collections,xml.etree.ElementTree as ET,hashlib
from fetch_metadata import OUT,ROOT,read,write
GENERIC=re.compile(r'uncultured|unclassified|environmental|metagenom|\bsp\.',re.I)
RANKS=['realm','subrealm','kingdom','subkingdom','phylum','subphylum','class','subclass','order','suborder','family','subfamily','genus','subgenus','species']
def main():
    original=read(ROOT/'analysis/viral_gt_figureB_taxonomy_20261018/tables/master_genome_dataset.tsv')
    summaries=json.loads((OUT/'sources/ncbi_accession_summaries.json').read_text()); summary={r['accessionversion']:r for r in summaries}
    assert len(summary)==515
    vmr=read(OUT/'sources/ictv/VMR_MSL41.v1.20260729__VMR_MSL41.tsv');msl=read(OUT/'sources/ictv/ICTV_Master_Species_List_2025_MSL41.v1__MSL.tsv')
    species={r['Species']:r for r in msl if r.get('Species')}; index=collections.defaultdict(list);genera=collections.defaultdict(list)
    for r in msl:
        if r.get('Genus'):genera[r['Genus']].append(r)
    for r in vmr:
        if not r.get('Species'):continue
        assert r['Species'] in species,r['Species']
        for rank in RANKS:
            assert r.get(rank.capitalize(),'')==species[r['Species']].get(rank.capitalize(),''),(r['Species'],rank)
        for a in re.findall(r'\b[A-Z]{1,6}_?\d{5,12}(?:\.\d+)?\b',r.get('Virus GENBANK accession','')):index[a.split('.')[0]].append(r)
    tax={};nodes={}
    for p in (OUT/'sources/http').iterdir():
        if p.suffix=='.json' or not p.read_bytes().lstrip().startswith(b'<?xml'):continue
        doc=ET.parse(p)
        for t in doc.getroot().findall('Taxon'):
            path=t.findall('LineageEx/Taxon')+[t];lineage={}
            for n in path:
                rank=n.findtext('Rank');name=n.findtext('ScientificName');tid=n.findtext('TaxId');nodes[tid]=dict(name=name,rank=rank)
                if rank in RANKS:lineage[rank]=name;lineage[rank+'_taxid']=tid
            row=dict(scientific_name=t.findtext('ScientificName'),taxid=t.findtext('TaxId'),lineage=lineage,source=str(p.relative_to(OUT)),ncbi_path=[n.findtext('TaxId') for n in path])
            tax[t.findtext('TaxId')]=row
            for a in t.findall('AkaTaxIds/TaxId'):tax[a.text]=row
    accession_rows=[];bygen=[];decisions=[]
    for g in original:
        evidence=[]
        for a in g['source_genome_accessions'].split(';'):
            assert a in summary,('no exact-version summary',a)
            s=summary[a];t=tax[str(s['taxid'])];lin=t['lineage'];matches=index[a.split('.')[0]]
            # ICTV accession tables often omit versions; accession entity is matched, not an ANI claim.
            exact=[v for v in matches if a in re.findall(r'\b[A-Z]{1,6}_?\d{5,12}(?:\.\d+)?\b',v['Virus GENBANK accession'])]
            candidates=exact or matches
            unique={v['Species'] for v in candidates}
            # Here "sp." qualifies the host Bathycoccus, not an unspecified virus.
            named_host_virus=bool(re.fullmatch(r'Bathycoccus sp\. RCC\d+ virus (?:BpV\d+|\d+)',s.get('organism','')) and s.get('organism')==t['scientific_name'])
            specific=named_host_virus or not bool(GENERIC.search(s.get('organism','')+' '+t['scientific_name']))
            chosen={};method='';status='Unresolved';conflict='';note=''
            if len(unique)==1:
                v=candidates[0];chosen={rank:v.get(rank.capitalize(),'') for rank in RANKS};method='ICTV_VMR_accession_exact' if exact else 'ICTV_VMR_accession_base';status='Authoritative'
                if lin.get('family') and chosen.get('family') and lin['family']!=chosen['family']:conflict='NCBI_vs_ICTV_family_difference_ICTV_MSL41_preferred'
                note=';'.join(sorted({v['Isolate ID'] for v in candidates}))
            elif len(unique)>1:status='Tentative';method='ICTV_accession_conflict';conflict='accession_maps_to_multiple_ICTV_species'
            elif specific:
                if lin.get('species') in species:
                    v=species[lin['species']];chosen={rank:v.get(rank.capitalize(),'') for rank in RANKS};method='NCBI_specific_species_to_ICTV_MSL41';status='Authoritative'
                    if lin.get('family') and chosen.get('family') and lin['family']!=chosen['family']:conflict='NCBI_vs_ICTV_family_difference_ICTV_MSL41_preferred'
                elif lin.get('genus') in genera and len({v['Family'] for v in genera[lin['genus']]})==1:
                    v=genera[lin['genus']][0];chosen={rank:v.get(rank.capitalize(),'') if rank not in ['species','subgenus'] else lin.get(rank,'') for rank in RANKS};method='NCBI_specific_accession_genus_to_ICTV_MSL41';status='Authoritative'
                    if lin.get('family') and chosen.get('family') and lin['family']!=chosen['family']:conflict='NCBI_vs_ICTV_family_difference_ICTV_MSL41_preferred'
                elif lin.get('family') or lin.get('genus'):
                    chosen={rank:lin.get(rank,'') for rank in RANKS};method='NCBI_specific_accession_taxid';status='Authoritative'
                else:chosen={rank:lin.get(rank,'') for rank in RANKS};method='NCBI_specific_accession_higher_rank_only';status='Authoritative'
            else:
                chosen={rank:lin.get(rank,'') for rank in RANKS};method='generic_taxid_inherited_only';status='Tentative' if lin.get('family') or lin.get('genus') else 'Unresolved'
            row=dict(genome_unit_id=g['genome_unit_id'],accession=a,updated_ncbi_taxid=str(s['taxid']),updated_organism=s.get('organism',''),updated_taxon_name=t['scientific_name'],updated_title=s['title'],length=s['slen'],topology=s.get('topology',''),ncbi_update_date=s.get('updatedate',''),source_taxid_changed=str(s['taxid'])!=g['source_ncbi_taxid'],specific_ncbi_source=specific,vmr_isolate_ids=note,metadata_status=status,metadata_method=method,conflict_flag=conflict,ncbi_family=lin.get('family',''),ncbi_genus=lin.get('genus',''),ncbi_species=lin.get('species',''),family=chosen.get('family',''),genus=chosen.get('genus',''),species=chosen.get('species',''),lineage_json=json.dumps(chosen,sort_keys=True),ncbi_summary_source=s['_source_file'],ncbi_taxonomy_source=t['source'])
            accession_rows.append(row);evidence.append(row)
        authoritative=[r for r in evidence if r['metadata_status']=='Authoritative']
        families={r['family'] for r in authoritative if r['family']};gens={r['genus'] for r in authoritative if r['genus']}
        conflict=len(families)>1 or len(gens)>1
        chosen=next((r for r in authoritative if r['family']),authoritative[0] if authoritative else evidence[0])
        status='Tentative' if conflict else 'Authoritative' if authoritative else chosen['metadata_status']
        family=chosen['family'] if not conflict else '';genus=chosen['genus'] if not conflict else ''
        new=dict(g)
        for rank in ['family','genus','species']:new['original_'+rank]=new.pop(rank);new['original_'+rank+'_taxid']=new.pop(rank+'_taxid')
        new.update(dict(metadata_status=status,metadata_family=family,metadata_genus=genus,metadata_species=chosen['species'] if not conflict else '',metadata_method=chosen['metadata_method'],metadata_lineage_json=chosen['lineage_json'] if not conflict else '{}',metadata_supporting_accessions=';'.join(r['accession'] for r in authoritative),metadata_conflict='alias_accessions_conflict' if conflict else ';'.join(sorted({r['conflict_flag'] for r in evidence if r['conflict_flag']})),updated_ncbi_taxids=';'.join(sorted({r['updated_ncbi_taxid'] for r in evidence})),updated_organisms=';'.join(sorted({r['updated_organism'] for r in evidence})),metadata_family_resolved=status=='Authoritative' and bool(family),metadata_genus_resolved=status=='Authoritative' and bool(genus),newly_family_resolved=status=='Authoritative' and bool(family) and g['conservative_family_included']!='True',metadata_reference='ICTV MSL41.v1; VMR_MSL41.v1.20260729; fresh NCBI nuccore ESummary/taxonomy EFetch',needs_genome_classification=not(status=='Authoritative' and bool(family) and bool(genus))))
        bygen.append(new)
        decisions.append(dict(genome_unit_id=g['genome_unit_id'],status=status,family=family,genus=genus,method=chosen['metadata_method'],supporting_accessions=new['metadata_supporting_accessions'],conflict=new['metadata_conflict'],rationale='No specific-rank rescue from generic taxid alone; ICTV accession matches or specific NCBI records required.'))
    write(OUT/'tables/accession_metadata_evidence.tsv',accession_rows);write(OUT/'tables/phase1_genome_taxonomy_master.tsv',bygen);write(OUT/'tables/metadata_decisions.tsv',decisions)
    recovered=[r for r in bygen if r['newly_family_resolved']];positive=lambda rs:sum(r['combined_status']=='candidate' for r in rs)
    brief=dict(total_genome_units=len(bygen),unique_source_taxids=len({r['source_ncbi_taxid'] for r in original}),unique_source_organism_labels=len({r['source_organism_labels'] for r in original}),unique_accessions=len(summary),authoritative_family_units=sum(r['metadata_family_resolved'] for r in bygen),authoritative_family_candidates=positive([r for r in bygen if r['metadata_family_resolved']]),additional_family_units=len(recovered),additional_family_candidates=positive(recovered),remaining_family_unresolved_units=sum(not r['metadata_family_resolved'] for r in bygen),remaining_family_unresolved_candidates=positive([r for r in bygen if not r['metadata_family_resolved']]),authoritative_genus_units=sum(r['metadata_genus_resolved'] for r in bygen),authoritative_genus_candidates=positive([r for r in bygen if r['metadata_genus_resolved']]),ictv_direct_accession_units=sum(any(x['metadata_method'].startswith('ICTV_VMR_accession') for x in accession_rows if x['genome_unit_id']==r['genome_unit_id']) for r in bygen),changed_accession_taxids=sum(r['source_taxid_changed'] for r in accession_rows),original_conservative_lost=[r['genome_unit_id'] for r in bygen if r['conservative_family_included']=='True' and not r['metadata_family_resolved']],metadata_method_counts=dict(collections.Counter(r['metadata_method'] for r in bygen)),taxid2100421=dict(units=sum(r['source_ncbi_taxid']=='2100421' for r in bygen),new_authoritative_family=sum(r['source_ncbi_taxid']=='2100421' and r['metadata_family_resolved'] for r in bygen)))
    (OUT/'reports/phase1_recovery.json').write_text(json.dumps(brief,indent=2));print(json.dumps(brief,indent=2))
if __name__=='__main__':main()
