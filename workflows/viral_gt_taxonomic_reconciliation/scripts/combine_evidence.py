"""Combine authoritative metadata and independent genome evidence, rank by rank."""
import sys
sys.dont_write_bytecode=True
import csv,json,re,collections,pathlib
from fetch_metadata import OUT,read,write
from taxonomy_mapping import TaxonomyMap,RANKS
GENERIC=re.compile(r'uncultured|unclassified|environmental|metagenom|\bsp\.',re.I)
THRESHOLDS=dict(non_gt_marker_count=5,non_gt_distinct_marker_profiles=3,non_gt_weight_fraction=.8,non_gt_marker_span_fraction=.25,fragment_marker_count=10,fragment_marker_span_fraction=.5,min_nonself_specific_references=2,min_reference_purity=.8,max_ambiguous_fraction=.05)
def f(r,k):return float(r.get(k,0) or 0)
def b(r,k):return str(r.get(k,''))=='True'
def main():
    master=read(OUT/'tables/phase1_genome_taxonomy_master.tsv');byid={r['genome_unit_id']:r for r in master};mapping=TaxonomyMap()
    marker={r['genome_unit_id']:r for r in read(OUT/'tables/genomad_genome_marker_evidence.tsv')};quality={r['genome_unit_id']:r for r in read(OUT/'tables/full_genome_quality.tsv')}
    exact=collections.defaultdict(set)
    reference_authority_options=collections.defaultdict(list)
    for row in read(OUT/'tables/exact_reference_sequence_matches.tsv'):
        rid=row['reference_genome_id'].split('.')[0];exact[row['genome_unit_id']].add(rid);g=byid[row['genome_unit_id']]
        if b(g,'metadata_family_resolved'):reference_authority_options[rid].append(json.loads(g['metadata_lineage_json']))
    reference_authority={rid:paths[0] for rid,paths in reference_authority_options.items() if len({p.get('family') for p in paths})==1}
    ledger=[];pred=collections.defaultdict(list);completeness=collections.defaultdict(list);method_coverage=[]
    specs=[(d,m,'vcontact3_'+d+'_'+m,None) for d in ['prokaryotes','eukaryotes'] for m in ['full','no_gt']]
    cross=json.loads((OUT/'reports/additional_reference_domain_checks.json').read_text())
    specs.extend((r['domain'],r['mode'],r['run'],r['genome_unit_ids']) for r in cross)
    for domain,mode,run,explicit in specs:
        p=OUT/'runs'/run/'exports/final_assignments.csv'
        assert p.exists(),str(p)
        with p.open(newline='') as stream:rows=list(csv.DictReader(stream))
        refs=[r for r in rows if r.get('Reference','').lower()=='true'];queries=[r for r in rows if r.get('Genome') in byid]
        assert len({r['Genome'] for r in refs})==len(refs),'Duplicate reference IDs'
        expected=set(explicit) if explicit is not None else {gid for gid,q in quality.items() if q['classification_pending']=='True' and domain in q['vcontact_reference_domains'].split(';')}
        observed={r['Genome'] for r in queries};missing=expected-observed
        method_coverage.append(dict(domain=domain,mode=mode,expected_queries=len(expected),reported_queries=len(observed),no_assignment_returned=len(missing),missing_genome_ids=';'.join(sorted(missing))))
        cp=p.parent/'completeness/completeness.csv'
        if cp.exists():
            with cp.open(newline='') as stream:
                for cr in csv.DictReader(stream):
                    if cr.get('') in byid:completeness[cr['']].append(dict(domain=domain,mode=mode,values=cr))
        ref_groups={rank:collections.defaultdict(list) for rank in ['family','genus']}
        for r in refs:
            for rank in ref_groups:ref_groups[rank][r.get(rank+'_prediction','')].append(r)
        for r in queries:
            gid=r['Genome'];g=byid[gid];aliases={a.split('.')[0] for a in g['source_genome_accessions'].split(';')}|exact[gid]
            lineage,reason=mapping.from_lineage({rank:r.get(rank+'_prediction','') for rank in RANKS})
            item=dict(genome_unit_id=gid,method='vConTACT3 3.1.6',reference='v230 RefSeq230/GenBank266',domain=domain,mode=mode,raw_family=r.get('family_prediction',''),raw_genus=r.get('genus_prediction',''),current_family=lineage.get('family',''),current_genus=lineage.get('genus',''),mapping_reason=reason,source_file=str(p.relative_to(OUT)),all_original_fields_json=json.dumps(r,sort_keys=True))
            for rank in ['family','genus']:
                raw=r.get(rank+'_prediction','');options={mapping.current(rank,n) for n in raw.split('|')}-{''}
                item['ambiguous_'+rank+'_label']='|' in raw
                item['ambiguous_'+rank+'_options']=';'.join(sorted(options)) if '|' in raw else ''
            for rank in ['family','genus']:
                counts=collections.Counter();ids=[];support_ids=collections.defaultdict(set);excluded_self=[];excluded_generic=[]
                for ref in ref_groups[rank][r.get(rank+'_prediction','')]:
                    rid=ref.get('Genome','')
                    if rid.split('.')[0] in aliases:excluded_self.append(rid);continue
                    name=ref.get('GenomeName','')
                    curated=reference_authority.get(rid.split('.')[0])
                    if GENERIC.search(name) and not curated and not re.search(r'Bathycoccus sp\. RCC\d+ virus',name):excluded_generic.append(rid);continue
                    current=curated or mapping.from_lineage({k:ref.get(k+'_reference','') for k in RANKS})[0]
                    if current.get(rank):counts[current[rank]]+=1;ids.append(rid);support_ids[current[rank]].add(rid)
                target=item['current_'+rank]
                item[rank+'_nonself_specific_reference_count']=counts[target] if target else 0
                item[rank+'_reference_purity']=counts[target]/sum(counts.values()) if target and counts else 0
                item[rank+'_reference_counts_json']=json.dumps(counts,sort_keys=True)
                item[rank+'_supporting_reference_ids']=';'.join(sorted(support_ids[target]))
                item[rank+'_excluded_self_refs']=';'.join(excluded_self);item[rank+'_excluded_generic_refs']=';'.join(excluded_generic)
            ledger.append(item);pred[gid].append(item)
        for gid in sorted(missing):
            item={k:'' for k in ledger[-1]};item.update(genome_unit_id=gid,method='vConTACT3 3.1.6',reference='v230 RefSeq230/GenBank266',domain=domain,mode=mode,mapping_reason='method_returned_no_assignment',source_file=str(p.relative_to(OUT)),all_original_fields_json='{}')
            for rank in ['family','genus']:
                item[rank+'_nonself_specific_reference_count']=0;item[rank+'_reference_purity']=0;item[rank+'_reference_counts_json']='{}'
            ledger.append(item);pred[gid].append(item)
    write(OUT/'tables/vcontact_genome_evidence.tsv',ledger)
    write(OUT/'tables/vcontact_query_coverage.tsv',method_coverage)
    decisions=[];final=[];sensitivity=[]
    for original in master:
        g=dict(original);gid=g['genome_unit_id'];m=marker.get(gid,{});q=quality[gid];prs=pred[gid]
        authoritative_family=g['metadata_family'] if b(g,'metadata_family_resolved') else ''
        authoritative_genus=g['metadata_genus'] if b(g,'metadata_genus_resolved') else ''
        eligible=[r for r in prs if r['mode']=='full' and r['current_family'] and r['family_nonself_specific_reference_count']>=1]
        vfamilies={r['current_family'] for r in eligible}
        ambiguous_family=any(r.get('ambiguous_family_label') for r in prs if r['mode']=='full')
        for r in prs:
            if r['mode']=='full' and r.get('ambiguous_family_options'):vfamilies.update(r['ambiguous_family_options'].split(';'))
        vconflict=len(vfamilies)>1 or ambiguous_family
        vf=next(iter(vfamilies)) if len(vfamilies)==1 else ''
        mf=m.get('no_gt_current_family','') if f(m,'no_gt_current_family_weight_fraction')>=.5 and f(m,'no_gt_current_family_marker_count')>=2 else ''
        vgenera={r['current_genus'] for r in eligible if r['current_genus'] and r['genus_nonself_specific_reference_count']>=1};vg=next(iter(vgenera)) if len(vgenera)==1 else ''
        ambiguous_genus=any(r.get('ambiguous_genus_label') for r in prs if r['mode']=='full')
        mg=m.get('no_gt_current_genus','') if f(m,'no_gt_current_genus_weight_fraction')>=.5 and f(m,'no_gt_current_genus_marker_count')>=2 else ''
        family_conflict=vconflict or bool(vf and mf and vf!=mf)
        genus_conflict=len(vgenera)>1 or ambiguous_genus or bool(vg and mg and vg!=mg)
        best=None;best_genus=None
        label_stable_family=False;label_stable_genus=False
        for full in [r for r in prs if r['mode']=='full']:
            for ablated in [r for r in prs if r['mode']=='no_gt' and r['domain']==full['domain']]:
                label_stable_family |= bool(full['current_family'] and full['current_family']==ablated['current_family'] and not full.get('ambiguous_family_label') and not ablated.get('ambiguous_family_label'))
                label_stable_genus |= bool(full['current_genus'] and full['current_genus']==ablated['current_genus'] and not full.get('ambiguous_genus_label') and not ablated.get('ambiguous_genus_label'))
        for r in eligible:
            ab=[x for x in prs if x['mode']=='no_gt' and x['domain']==r['domain']]
            if not r.get('ambiguous_family_label') and any(not x.get('ambiguous_family_label') and x['current_family']==r['current_family'] and x['family_nonself_specific_reference_count']>=2 and x['family_reference_purity']>=.8 for x in ab) and r['family_nonself_specific_reference_count']>=2 and r['family_reference_purity']>=.8:
                best=r
            if r['current_genus'] and not r.get('ambiguous_genus_label') and any(not x.get('ambiguous_genus_label') and x['current_genus']==r['current_genus'] and x['genus_nonself_specific_reference_count']>=2 and x['genus_reference_purity']>=.8 for x in ab) and r['genus_nonself_specific_reference_count']>=2 and r['genus_reference_purity']>=.8:
                best_genus=r
        fragment=q['source_completeness']=='fragment_or_scaffold';count=10 if fragment else 5;span=.5 if fragment else .25
        marker_strong=(f(m,'no_gt_current_family_marker_count')>=count and f(m,'no_gt_current_family_unique_marker_count')>=3 and f(m,'no_gt_current_family_weight_fraction')>=.8 and f(m,'no_gt_current_family_marker_span_fraction')>=span)
        quality_ok=f(q,'n_ambiguous')/f(q,'genome_length')<.05 and b(q,'gt_ablation_complete')
        high_family=bool(not family_conflict and best and mf==best['current_family'] and marker_strong and quality_ok)
        predicted_family=best['current_family'] if high_family else ''
        family=authoritative_family or predicted_family
        family_status='Authoritative' if authoritative_family else 'High-confidence predicted' if high_family else 'Tentative' if family_conflict or vf or mf or g['metadata_status']=='Tentative' else 'Unresolved'
        tentative_family=next(iter(({vf,mf}-{''}))) if not family_conflict and len({vf,mf}-{''})==1 else ''
        authority_family_conflict=bool(authoritative_family and high_family and predicted_family!=authoritative_family)
        if authority_family_conflict:family='';family_status='Tentative';tentative_family=''
        high_genus=bool(family and not genus_conflict and best_genus and best_genus['current_genus']==mg and best_genus['current_family']==family and f(m,'no_gt_current_genus_marker_count')>=count and f(m,'no_gt_current_genus_unique_marker_count')>=3 and f(m,'no_gt_current_genus_weight_fraction')>=.8 and f(m,'no_gt_current_genus_marker_span_fraction')>=span and quality_ok)
        genus=authoritative_genus or (best_genus['current_genus'] if high_genus else '')
        genus_status='Authoritative' if authoritative_genus else 'High-confidence predicted' if high_genus else 'Tentative' if genus_conflict or vg or mg or g['metadata_status']=='Tentative' else 'Unresolved'
        tentative_genus=next(iter(({vg,mg}-{''}))) if not genus_conflict and len({vg,mg}-{''})==1 else ''
        authority_genus_conflict=bool(authoritative_genus and high_genus and best_genus['current_genus']!=authoritative_genus)
        if authority_genus_conflict or authority_family_conflict:genus='';genus_status='Tentative';tentative_genus=''
        # A named genus gives a current family only if the MSL41 parent is explicit.
        exploratory_family=family or tentative_family;exploratory_genus=genus or tentative_genus
        if exploratory_genus in mapping.index['genus']:
            gp=mapping.index['genus'][exploratory_genus].get('family','')
            if exploratory_family and gp and gp!=exploratory_family:exploratory_genus='';genus_conflict=True
        reasons=[]
        if not family:
            if family_conflict:reasons.append('conflicting_genome_based_families')
            if not best:reasons.append('insufficient_or_GT_sensitive_vcontact_reference_support')
            if not marker_strong:reasons.append('insufficient_non_GT_marker_strength_or_distribution')
            if not quality_ok:reasons.append('sequence_quality_or_GT_ablation_uncertain')
            if not vf and not mf:reasons.append('no_independent_named_family_consensus')
        if authoritative_family and (vf and vf!=authoritative_family or mf and mf!=authoritative_family):reasons.append('prediction_metadata_difference_review_no_overwrite')
        if authoritative_family and family_conflict:reasons.append('model_ambiguity_recorded_authoritative_taxonomy_retained')
        if authority_family_conflict:reasons.append('authoritative_metadata_vs_corroborated_genome_family_conflict_withheld')
        if authority_genus_conflict:reasons.append('authoritative_metadata_vs_corroborated_genome_genus_conflict_withheld')
        higher_rank='';higher_taxon='';higher_agreement=0;higher_markers=0
        for rank in ['order','class','phylum','realm']:
            if f(m,'no_gt_'+rank+'_weight_fraction')>=.5 and f(m,'no_gt_'+rank+'_marker_count')>=3:
                name=mapping.current(rank,m.get('no_gt_'+rank,''))
                if name:
                    higher_rank=rank;higher_taxon=name;higher_agreement=f(m,'no_gt_'+rank+'_weight_fraction');higher_markers=f(m,'no_gt_'+rank+'_marker_count');break
        g.update(dict(authoritative_family=authoritative_family,authoritative_genus=authoritative_genus,reconciled_family=family,reconciled_genus=genus,family_assignment_status=family_status,genus_assignment_status=genus_status,taxonomy_evidence_status=family_status,
            exploratory_family=exploratory_family,exploratory_genus=exploratory_genus,tentative_family=tentative_family,tentative_genus=tentative_genus,genome_based_family_conflict=family_conflict,genome_based_genus_conflict=genus_conflict,decision_reasons=';'.join(reasons),
            vcontact_current_family_options=';'.join(sorted(vfamilies)),vcontact_current_genus_options=';'.join(sorted(vgenera)),vcontact_family_stable_without_selected_GT=bool(best),vcontact_genus_stable_without_selected_GT=bool(best_genus),marker_supported_higher_rank=higher_rank,marker_supported_higher_taxon=higher_taxon,marker_supported_higher_agreement=higher_agreement,marker_supported_higher_marker_count=higher_markers,higher_rank_evidence='geNomad_non_GT_markers_single_method' if higher_rank else '',
            genomad_non_GT_current_family=m.get('no_gt_current_family',''),genomad_non_GT_family_marker_count=f(m,'no_gt_current_family_marker_count'),genomad_non_GT_family_weight_fraction=f(m,'no_gt_current_family_weight_fraction'),genomad_non_GT_current_genus=m.get('no_gt_current_genus',''),genomad_non_GT_genus_marker_count=f(m,'no_gt_current_genus_marker_count'),genomad_non_GT_genus_weight_fraction=f(m,'no_gt_current_genus_weight_fraction'),
            genomad_taxonomic_gene_coverage=f(m,'taxonomic_gene_coverage'),genomad_non_GT_marker_span_fraction=f(m,'no_gt_marker_span_fraction'),GT_taxonomic_bitscore_fraction=f(m,'gt_taxonomic_bitscore_fraction'),genome_length=q['genome_length'],source_genome_completeness=q['source_completeness'],n_cached_proteins=q['n_cached_translations'],n_ambiguous_nucleotides=q['n_ambiguous'],
            new_authoritative_family_relative_to_original=b(g,'newly_family_resolved'),new_predicted_family=bool(not authoritative_family and predicted_family),reference_versions='vConTACT3 3.1.6/v230; geNomad 1.13.0/db1.9(MSL39); current mapping MSL41.v1',vcontact_runtime_overlay='pandas2.3.3;pytz2026.5',vcontact_evidence_file='tables/vcontact_genome_evidence.tsv',genomad_evidence_file='tables/genomad_genome_marker_evidence.tsv'))
        g.update(vcontact_experimental_completeness_json=json.dumps(completeness[gid],sort_keys=True),vcontact_no_assignment_runs=';'.join(r['domain']+'_'+r['mode'] for r in prs if r['mapping_reason']=='method_returned_no_assignment'))
        g.update(genomad_non_GT_family_marker_span_fraction=f(m,'no_gt_current_family_marker_span_fraction'),genomad_non_GT_family_unique_marker_profiles=f(m,'no_gt_current_family_unique_marker_count'),genomad_non_GT_genus_marker_span_fraction=f(m,'no_gt_current_genus_marker_span_fraction'),genomad_non_GT_genus_unique_marker_profiles=f(m,'no_gt_current_genus_unique_marker_count'))
        g.update(vcontact_family_stable_without_selected_GT=label_stable_family,vcontact_genus_stable_without_selected_GT=label_stable_genus,vcontact_family_reference_and_ablation_rule_passed=bool(best),vcontact_genus_reference_and_ablation_rule_passed=bool(best_genus))
        g['vcontact_reference_domains_used']=';'.join(sorted({r['domain'] for r in prs}))
        g['exact_reference_sequence_ids']=';'.join(sorted(exact[gid]))
        g.update(authoritative_vs_corroborated_family_conflict=authority_family_conflict,authoritative_vs_corroborated_genus_conflict=authority_genus_conflict)
        final.append(g)
        decisions.append(dict(genome_unit_id=gid,family=family,family_status=family_status,genus=genus,genus_status=genus_status,metadata_method=g['metadata_method'],reason=g['decision_reasons'],non_GT_marker_rule_passed=marker_strong,vcontact_reference_and_ablation_rule_passed=bool(best),sequence_quality_rule_passed=quality_ok,family_method_conflict=family_conflict,fragment_more_stringent_rule=fragment))
        for threshold in [.5,.7,.8,.9]:
            sensitivity.append(dict(genome_unit_id=gid,weight_threshold=threshold,candidate=g['combined_status']=='candidate',new_predicted_family=bool(not authoritative_family and not family_conflict and best and m.get('no_gt_current_family')==best['current_family'] and f(m,'no_gt_current_family_marker_count')>=count and f(m,'no_gt_current_family_unique_marker_count')>=3 and f(m,'no_gt_current_family_weight_fraction')>=threshold and f(m,'no_gt_current_family_marker_span_fraction')>=span and quality_ok)))
    assert len(final)==466 and collections.Counter(r['combined_status'] for r in final)=={'candidate':245,'not_selected':215,'unresolved':6}
    write(OUT/'tables/final_genome_taxonomy_master.tsv',final);write(OUT/'tables/final_classification_decisions.tsv',decisions);write(OUT/'tables/confidence_threshold_sensitivity.tsv',sensitivity)
    (OUT/'reports/confidence_rules.json').write_text(json.dumps(THRESHOLDS,indent=2))
    write(OUT/'tables/unresolved_family_genomes.tsv',[r for r in final if not r['reconciled_family']]);write(OUT/'tables/unresolved_genus_genomes.tsv',[r for r in final if not r['reconciled_genus']])
    brief=dict(genomes=466,candidates=245,negative=215,screening_unresolved=6,family_status_counts=dict(collections.Counter(r['family_assignment_status'] for r in final)),genus_status_counts=dict(collections.Counter(r['genus_assignment_status'] for r in final)),authoritative_family_units=sum(bool(r['authoritative_family']) for r in final),authoritative_family_candidates=sum(bool(r['authoritative_family']) and r['combined_status']=='candidate' for r in final),reconciled_family_units=sum(bool(r['reconciled_family']) for r in final),reconciled_family_candidates=sum(bool(r['reconciled_family']) and r['combined_status']=='candidate' for r in final),reconciled_genus_units=sum(bool(r['reconciled_genus']) for r in final),reconciled_genus_candidates=sum(bool(r['reconciled_genus']) and r['combined_status']=='candidate' for r in final),new_predicted_family_units=sum(r['new_predicted_family'] for r in final),new_predicted_family_candidates=sum(r['new_predicted_family'] and r['combined_status']=='candidate' for r in final),taxid2100421=dict(units=sum(r['source_ncbi_taxid']=='2100421' for r in final),reconciled_family=sum(r['source_ncbi_taxid']=='2100421' and bool(r['reconciled_family']) for r in final),exploratory_family=sum(r['source_ncbi_taxid']=='2100421' and bool(r['exploratory_family']) for r in final)))
    (OUT/'reports/final_reconciliation_summary.json').write_text(json.dumps(brief,indent=2));print(json.dumps(brief,indent=2))
if __name__=='__main__':main()
