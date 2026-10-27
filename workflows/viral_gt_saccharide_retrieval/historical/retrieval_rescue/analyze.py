"""Read-only analysis of frozen CORE_ONLY vectors; no model imports or inference."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]
OUT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'analysis/bgcprophet_reference/.deps'))
import numpy as np
import csv,json,hashlib,collections,math,platform
BASE=ROOT/'analysis/viral_gt_mining_20261017/runs/baseline_594_v2'
PRIOR=ROOT/'analysis/viral_gt_interpretation_20261017'
REF=ROOT/'analysis/viral_glycan_benchmark_20261016'
PILOT=ROOT/'analysis/viral_glycan_benchmark_cluster7_ARF08726_ARF08739_20261016'
CUTS=[.4,.5,.6]
def read(path):
    with Path(path).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def write(name,rows,fields=None):
    path=OUT/name;path.parent.mkdir(parents=True,exist_ok=True)
    fields=fields or list(dict.fromkeys(k for r in rows for k in r))
    with path.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
def save(name,obj):(OUT/name).write_text(json.dumps(obj,indent=2)+'\n',encoding='utf-8')
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()
def yes(x):return str(x).lower() in {'true','1','yes'}
def verify_inputs():
    counts={}
    for directory,manifest in [(BASE,PRIOR/'freeze/baseline_files.tsv'),(PRIOR,PRIOR/'artifact_checksums.tsv')]:
        rr=read(manifest)
        for r in rr:assert sha(directory/r['path'])==r['sha256'],str(directory/r['path'])
        counts[directory.name]=len(rr)
    return counts
def retrieve(q,refs,ids,labels):
    sims=refs@(q.astype(np.float64)/np.linalg.norm(q.astype(np.float64)))
    order=np.lexsort((ids,-sims))[:20]
    return {k:int(labels[order[:k]].sum()) for k in [5,10,20]},order,sims
def aggregate_status(window_ids,candidate_ids,inferred_ids):
    if set(window_ids)&candidate_ids:return 'candidate'
    return 'unresolved' if set(window_ids)-inferred_ids else 'not_selected'

def run():
    preservation=verify_inputs()
    windows=read(BASE/'viral_gt_windows.tsv');wm={r['window_id']:r for r in windows}
    loci=read(BASE/'candidate_loci.tsv');genomes=read(BASE/'genome_units.tsv')
    lm={w:l for l in loci for w in l['window_ids'].split(';')};gm={w:g for g in genomes for w in g['window_ids'].split(';')}
    tax={r['genome_unit_id']:r for r in read(PRIOR/'tables/genome_relatedness_metadata.tsv')}
    prior={r['window_id']:r for r in read(PRIOR/'tables/window_interpretation.tsv')}
    preds={r['window_id']:r for r in read(BASE/'viral_bgc_predictions.tsv')};loc={r['window_id']:r for r in read(BASE/'saccharide_localization.tsv')}
    anchors=collections.defaultdict(set)
    for r in read(BASE/'anchor_mapping.tsv'):
        if r['window_id']:anchors[r['window_id']].add(r['input_accession'])
    manifest=read(BASE/'embedding_manifest.tsv');vectors=np.load(BASE/'viral_bgc_embeddings.npy',allow_pickle=False)
    refmeta=read(REF/'mibig_matched_reference_metadata.tsv');refvec=np.load(REF/'mibig_CORE_ONLY_embeddings.npy',allow_pickle=False)
    assert sha(REF/'mibig_CORE_ONLY_embeddings.npy')=='78036e1636eaeeab90f4d157794c868d9b57033195939d83960d8e70bb831e22'
    assert sha(REF/'mibig_matched_reference_metadata.tsv')=='5184f8f3199da2d1fac90e85ff11e2a30b443b70aac60992136718856c18ba4c'
    assert vectors.shape==(737,320) and refvec.shape==(2502,320)
    assert np.isfinite(vectors).all() and np.isfinite(refvec).all()
    assert len({r['BGC'] for r in refmeta})==2502 and all(int(r['embedding_row'])==i for i,r in enumerate(refmeta))
    assert all(r['strategy']=='CORE_ONLY' for r in manifest)
    refs=refvec.astype(np.float64);refs/=np.linalg.norm(refs,axis=1)[:,None]
    ids=np.array([r['BGC'] for r in refmeta]);labels=np.array(['Saccharide' in r['label'].split() for r in refmeta]);assert int(labels.sum())==195
    rows=[];neighbors=[];retrieval_checks=0
    for i,m in enumerate(manifest):
        assert int(m['embedding_row'])==i
        wid=m['window_id'];w=wm[wid];p=preds[wid];l=lm[wid];g=gm[wid];t=tax[g['genome_unit_id']]
        counts,order,sims=retrieve(vectors[i],refs,ids,labels)
        for k in [5,10,20]:
            assert abs(counts[k]/k-float(loc[wid][f'saccharide_fraction_top{k}']))<1e-12
            retrieval_checks+=1
        assert abs(float(prior[wid]['p10'])-counts[10]/10)<1e-12
        r=dict(window_id=wid,embedding_row=i,locus_id=l['locus_id'],genome_unit_id=g['genome_unit_id'],
            genome_accession=w['genome_accession'],source_genome_accessions=g['source_genome_accessions'],
            anchor_protein_id=w['anchor_protein_id'],cohort_anchor_accessions=';'.join(sorted(anchors[wid])),
            taxid=g['taxid'],source_organism_labels=g['source_organism_labels'],generic_source_label=yes(t['generic_source_label']),
            **{k:t[k] for k in ['family','family_taxid','genus','genus_taxid','species','species_taxid']},
            start_1based=w['start_1based'],end_1based=w['end_1based'],n_genes=int(w['n_genes']),complete=int(w['n_genes'])==15,
            predicted_labels=p['predicted_labels'],classifier_saccharide=yes(p['saccharide_predicted']),
            saccharide_score=float(p['saccharide_probability']),frozen_decision_threshold=float(p['threshold']),
            other_score=float(p['probability_Other']),context_signature=prior[wid]['context_signature'],
            **{f'saccharide_neighbors_top{k}':counts[k] for k in [5,10,20]},**{f'p{k}':counts[k]/k for k in [5,10,20]},
            reference_fraction=195/2502,descriptive_P10_enrichment=(counts[10]/10)/(195/2502),
            nearest_reference=str(ids[order[0]]),nearest_reference_cosine=float(sims[order[0]]),nearest_reference_labels=refmeta[order[0]]['label'])
        for c in CUTS:r[f'exploratory_B_P10_ge_{c}']=not r['classifier_saccharide'] and counts[10]>=round(c*10)
        rows.append(r)
        neighbors.extend(dict(window_id=wid,rank=j+1,reference_BGC=str(ids[o]),reference_labels=refmeta[o]['label'],
            reference_saccharide=bool(labels[o]),cosine=float(sims[o])) for j,o in enumerate(order))
    assert len(rows)==737 and sum(r['classifier_saccharide'] for r in rows)==246
    write('tables/all_inferred_windows.tsv',rows);write('tables/top20_reference_neighbors.tsv',neighbors)
    neg=[r for r in rows if not r['classifier_saccharide']];assert len(neg)==491
    write('tables/all_classifier_negative_windows.tsv',neg)
    discordant=sorted([r for r in neg if r['saccharide_neighbors_top10']>=3],key=lambda r:(-r['p10'],-r['p20'],-r['saccharide_score'],r['window_id']))
    write('tables/discordant_candidates_P10_ge_0.3.tsv',discordant)
    distribution=[]
    for k in [5,10,20]:
        for n in range(k+1):
            rr=[r for r in neg if r[f'saccharide_neighbors_top{k}']==n]
            distribution.append(dict(k=k,n_saccharide_neighbors=n,p_at_k=n/k,n_classifier_negative_windows=len(rr),
                fraction_of_491=len(rr)/491,n_loci=len({r['locus_id'] for r in rr}),n_genomes=len({r['genome_unit_id'] for r in rr})))
    write('tables/classifier_negative_Pk_distribution.tsv',distribution)
    thresholdcounts=[]
    for k in [5,10,20]:
        for c in [.3,.4,.5,.6,.8]:
            rr=[r for r in neg if r[f'saccharide_neighbors_top{k}']>=math.ceil(c*k-1e-10)]
            thresholdcounts.append(dict(k=k,cutoff=c,min_neighbors=math.ceil(c*k-1e-10),n_windows=len(rr),
                fraction_classifier_negative=len(rr)/491,n_loci=len({r['locus_id'] for r in rr}),n_genomes=len({r['genome_unit_id'] for r in rr})))
    write('tables/classifier_negative_cutoff_counts.tsv',thresholdcounts)
    inf=set(preds);A={r['window_id'] for r in rows if r['classifier_saccharide']}
    entities={'window':{w:{w} for w in wm},'locus':{l['locus_id']:set(l['window_ids'].split(';')) for l in loci},
        'genome':{g['genome_unit_id']:set(g['window_ids'].split(';')) for g in genomes}}
    status_rows=[];summary=[];memberships={};additional=[];cutsets={}
    for c in CUTS:
        B={r['window_id'] for r in neg if r['saccharide_neighbors_top10']>=round(c*10)};sets={'A':A,'B':B,'C':A|B};cutsets[c]=sets
        assert not A&B
        for level,ee in entities.items():
            hit={s:{e for e,ww in ee.items() if ww&ss} for s,ss in sets.items()};memberships[(c,level)]=hit
            for s,ss in sets.items():
                statuses={e:aggregate_status(ww,ss,inf) for e,ww in ee.items()};count=collections.Counter(statuses.values())
                summary.append(dict(P10_cutoff=c,set=s,level=level,total_entities=len(ee),any_inferred=sum(bool(ww&inf) for ww in ee.values()),
                    fully_inferred=sum(ww<=inf for ww in ee.values()),candidates=count['candidate'],not_selected=count['not_selected'],
                    unresolved=count['unresolved'],outcome_resolved=count['candidate']+count['not_selected'],
                    candidate_fraction_among_any_inferred=count['candidate']/sum(bool(ww&inf) for ww in ee.values()),
                    B_overlap_A=len(hit['B']&hit['A']),additional_B_absent_A=len(hit['B']-hit['A']),union_gain_over_A=len(hit['C'])-len(hit['A'])))
                for e,ww in ee.items():status_rows.append(dict(P10_cutoff=c,set=s,level=level,entity_id=e,status=statuses[e],n_windows=len(ww),
                    n_inferred_windows=len(ww&inf),n_candidate_windows=len(ww&ss),A_candidate=e in hit['A'],B_candidate=e in hit['B'],
                    C_candidate=e in hit['C'],newly_prioritized_by_retrieval=e in hit['B']-hit['A']))
            if level in {'locus','genome'}:
                for e in sorted(hit['B']-hit['A']):
                    ww=ee[e];gg=gm[next(iter(ww))];t=tax[gg['genome_unit_id']];rr=[r for r in rows if r['window_id'] in ww&B]
                    additional.append(dict(P10_cutoff=c,level=level,entity_id=e,genome_unit_id=gg['genome_unit_id'],
                        source_genome_accessions=gg['source_genome_accessions'],source_organism_labels=gg['source_organism_labels'],taxid=gg['taxid'],
                        baseline_A_status=aggregate_status(ww,A,inf),**{k:t[k] for k in ['family','family_taxid','genus','genus_taxid','generic_source_label']},
                        triggering_window_ids=';'.join(sorted(r['window_id'] for r in rr)),trigger_anchor_proteins=';'.join(sorted({r['anchor_protein_id'] for r in rr})),
                        trigger_cohort_accessions=';'.join(sorted({a for r in rr for a in anchors[r['window_id']]})),maximum_P10=max(r['p10'] for r in rr),
                        maximum_P20=max(r['p20'] for r in rr),maximum_saccharide_score=max(r['saccharide_score'] for r in rr)))
    write('tables/three_set_counts.tsv',summary);write('tables/entity_set_status.tsv',status_rows);write('tables/additional_loci_and_genomes.tsv',additional)
    write('tables/unresolved_windows.tsv',[dict(window_id=wid,locus_id=lm[wid]['locus_id'],genome_unit_id=gm[wid]['genome_unit_id'],
        genome_accession=wm[wid]['genome_accession'],anchor_protein_id=wm[wid]['anchor_protein_id'],reason=wm[wid]['inference_exclusion_reason'],
        status='unresolved_for_classifier_and_retrieval') for wid in sorted(set(wm)-inf)])
    taxrows=[];breadth=[];novel=[]
    for view in ['recorded','non_generic_source_sensitivity']:
        for rank in ['family','genus']:
            groups=collections.defaultdict(list)
            for g in genomes:
                t=tax[g['genome_unit_id']]
                key=('Generic source; rank withheld','') if view!='recorded' and yes(t['generic_source_label']) else (t[rank] or 'Unresolved',t[rank+'_taxid'])
                groups[key].append(g)
            for c in CUTS:
                ghit=memberships[(c,'genome')];lhit=memberships[(c,'locus')];lineages={s:set() for s in 'ABC'}
                for (name,tid),gg in sorted(groups.items()):
                    gset={g['genome_unit_id'] for g in gg};lset={l['locus_id'] for l in loci if l['genome_unit_id'] in gset}
                    tr=dict(view=view,rank=rank,taxon=name,rank_taxid=tid,P10_cutoff=c,total_genomes=len(gg),
                        any_inferred_genomes=sum(bool(set(g['window_ids'].split(';'))&inf) for g in gg),
                        fully_inferred_genomes=sum(set(g['window_ids'].split(';'))<=inf for g in gg),total_loci=len(lset),
                        any_inferred_loci=sum(bool(entities['locus'][l]&inf) for l in lset),total_windows=sum(int(g['n_windows']) for g in gg),
                        inferred_windows=sum(len(set(g['window_ids'].split(';'))&inf) for g in gg),
                        generic_source_genomes=sum(yes(tax[g['genome_unit_id']]['generic_source_label']) for g in gg))
                    for s in 'ABC':
                        tr[s+'_candidate_genomes']=len(gset&ghit[s]);tr[s+'_candidate_loci']=len(lset&lhit[s])
                        tr[s+'_unresolved_genomes']=sum(aggregate_status(g['window_ids'].split(';'),cutsets[c][s],inf)=='unresolved' for g in gg)
                        if tid and gset&ghit[s]:lineages[s].add(tid)
                    tr['additional_genomes']=len(gset&(ghit['B']-ghit['A']));tr['additional_loci']=len(lset&(lhit['B']-lhit['A']))
                    tr['B_rank_absent_A']=bool(tid and tr['B_candidate_genomes'] and not tr['A_candidate_genomes']);taxrows.append(tr)
                    if tr['B_rank_absent_A']:novel.append(tr)
                known={g['genome_unit_id'] for (name,tid),gg in groups.items() if tid for g in gg}
                for s in 'ABC':
                    candidates=ghit[s]
                    breadth.append(dict(view=view,rank=rank,P10_cutoff=c,set=s,distinct_candidate_taxa=len(lineages[s]),candidate_genomes=len(candidates),
                        candidate_genomes_rank_resolved=len(candidates&known),candidate_genomes_rank_unresolved_or_withheld=len(candidates-known),
                        newly_represented_taxa_vs_A=len(lineages[s]-lineages['A']),candidate_genomes_generic_source=sum(yes(tax[g]['generic_source_label']) for g in candidates)))
    write('tables/taxonomy_by_candidate_set.tsv',taxrows);write('tables/taxonomic_breadth.tsv',breadth)
    write('tables/taxa_retrieval_B_absent_classifier_A.tsv',novel,list(taxrows[0]))
    clusterrows=read(PRIOR/'tables/genome_similarity_clusters.tsv');redundancy=[];quality=[]
    for c in CUTS:
        for j in [.2,.5,.8]:
            assignment={r['genome_unit_id']:r['similarity_cluster'] for r in clusterrows if float(r['Jaccard_cutoff'])==j}
            ghit=memberships[(c,'genome')];components={s:{assignment[g] for g in ghit[s]} for s in 'ABC'}
            redundancy.append(dict(P10_cutoff=c,Jaccard_cutoff=j,total_components=len(set(assignment.values())),
                any_inferred_components=len({assignment[g['genome_unit_id']] for g in genomes if set(g['window_ids'].split(';'))&inf}),
                **{s+'_candidate_components':len(components[s]) for s in 'ABC'},additional_components_absent_A=len(components['B']-components['A']),
                components_containing_additional_genomes=len({assignment[g] for g in ghit['B']-ghit['A']})))
        for label,rr in [('all',rows),('complete_15gene',[r for r in rows if r['complete']]),('incomplete',[r for r in rows if not r['complete']]),
                         ('unique_contexts',list({r['context_signature']:r for r in rows}.values()))]:
            for s in 'ABC':
                selected=[r for r in rr if r['window_id'] in cutsets[c][s]]
                quality.append(dict(P10_cutoff=c,subset=label,set=s,total_contexts=len(rr),candidate_contexts=len(selected),
                    candidate_loci=len({r['locus_id'] for r in selected}),candidate_genomes=len({r['genome_unit_id'] for r in selected})))
    write('tables/genome_redundancy_sensitivity.tsv',redundancy);write('tables/window_quality_and_overlap_sensitivity.tsv',quality)
    concordance=[]
    for c in CUTS:
        for k in [5,20]:
            s10={r['window_id'] for r in neg if r['p10']>=c};sk={r['window_id'] for r in neg if r[f'p{k}']>=c}
            concordance.append(dict(cutoff=c,comparison_k=k,P10_B_windows=len(s10),Pk_B_windows=len(sk),both=len(s10&sk),P10_only=len(s10-sk),Pk_only=len(sk-s10)))
    write('tables/retrieval_k_sensitivity.tsv',concordance)
    pilot=check_pilot(refs,ids,labels,refmeta,manifest,wm)
    provenance=[]
    for directory,names in [(BASE,['viral_gt_windows.tsv','candidate_loci.tsv','genome_units.tsv','viral_bgc_predictions.tsv',
            'saccharide_localization.tsv','viral_bgc_embeddings.npy','embedding_manifest.tsv','anchor_mapping.tsv']),
        (PRIOR,['tables/window_interpretation.tsv','tables/genome_relatedness_metadata.tsv','tables/genome_similarity_clusters.tsv']),
        (REF,['mibig_CORE_ONLY_embeddings.npy','mibig_matched_reference_metadata.tsv','reference_validation.json'])]:
        provenance.extend(dict(path=str(p),sha256=sha(p),bytes=p.stat().st_size) for p in [directory/n for n in names])
    write('input_provenance.tsv',provenance)
    save('validation.json',dict(status='PASS',frozen_input_manifests_verified=preservation,baseline_windows=757,inferred_windows=737,
        classifier_negative_windows=491,recomputed_Pk_values_matching_frozen=retrieval_checks,model_inference_reruns=0,ESM_encoding_reruns=0,structure_predictions=0,
        unique_base_BGC_references=2502,saccharide_references=195,reference_augmentation_rows_used=0,no_binary_biological_boundary_claim=True,
        python=platform.python_version(),numpy=np.__version__,cluster7=pilot,baseline_and_prior_unchanged_after_analysis=verify_inputs()))
    save('summary.json',dict(three_set_counts=summary,negative_cutoff_counts=thresholdcounts,taxonomic_breadth=breadth,new_taxa=novel,pilot=pilot))
    print(json.dumps(dict(status='PASS',P10_counts=[r for r in thresholdcounts if r['k']==10],three_sets=summary,
        new_taxa=[dict(view=r['view'],rank=r['rank'],taxon=r['taxon'],cutoff=r['P10_cutoff'],additional_genomes=r['additional_genomes']) for r in novel],pilot=pilot),indent=2))

def check_pilot(refs,ids,labels,refmeta,baseline_manifest,wm):
    selected=['cluster7_embeddings.npy','embedding_metadata.tsv','classification_scores.tsv','input_strategy_manifest.tsv','inference_validation.json',
        'three_viral_loci_comparison.tsv','cluster7_core_gene_manifest.tsv','scripts/cluster7.py','reference_reuse_integrity.tsv']
    hashes={r['path']:r['sha256'] for r in read(PILOT/'artifact_checksums.tsv')}
    for name in selected:assert sha(PILOT/name)==hashes[name],name
    meta=read(PILOT/'embedding_metadata.tsv');i=next(i for i,r in enumerate(meta) if r['strategy']=='CORE_ONLY' and r['variant']=='primary')
    assert int(meta[i]['embedding_row'])==i and int(meta[i]['input_gene_count'])==14
    m=next(r for r in read(PILOT/'input_strategy_manifest.tsv') if r['strategy']=='CORE_ONLY')
    assert set(m['membership_mask'].split())=={'1'} and len(m['ordered_protein_ids'].split())==14
    proof=json.loads((PILOT/'inference_validation.json').read_text())
    assert proof['status']=='PASS' and proof['weights_unchanged'] and proof['hook_scores_unchanged']
    assert proof['checkpoint_sha256']==baseline_manifest[0]['classifier_sha256']
    assert proof['esm_checkpoint_sha256']==baseline_manifest[0]['esm_sha256'] and proof['ESM_layer']==6 and proof['truncation']==1022
    q=np.load(PILOT/'cluster7_embeddings.npy',allow_pickle=False)[i];counts,order,sims=retrieve(q,refs,ids,labels)
    score=next(r for r in read(PILOT/'classification_scores.tsv') if r['strategy']=='CORE_ONLY')
    historical=next(r for r in read(PILOT/'three_viral_loci_comparison.tsv') if r['locus_id'].startswith('Cluster7') and r['strategy']=='CORE_ONLY')
    for k in [5,10,20]:assert counts[k]/k==float(historical[f'saccharide_P_at_{k}'])
    assert not yes(score['saccharide_predicted'])
    write('tables/cluster7_CORE_ONLY_neighbors.tsv',[dict(rank=j+1,BGC=str(ids[o]),labels=refmeta[o]['label'],saccharide=bool(labels[o]),cosine=float(sims[o])) for j,o in enumerate(order)])
    proteins=set(m['ordered_protein_ids'].split());basalproteins={p for r in baseline_manifest for p in r['ordered_protein_ids'].split()}
    baseline_predictions={r['window_id']:r for r in read(BASE/'viral_bgc_predictions.tsv')}
    baseline_localization={r['window_id']:r for r in read(BASE/'saccharide_localization.tsv')}
    mappings=read(BASE/'anchor_mapping.tsv');overlap=[]
    for bm in baseline_manifest:
        shared=proteins&set(bm['ordered_protein_ids'].split())
        if bm['genome_accession']!=meta[i]['genome_id'] or not shared:continue
        wid=bm['window_id'];bp=baseline_predictions[wid];bl=baseline_localization[wid]
        overlap.append(dict(window_id=wid,anchor_protein_id=bm['anchor_protein_id'],
            cohort_anchor_accessions=';'.join(sorted({r['input_accession'] for r in mappings if r['window_id']==wid})),
            n_genes=int(bm['n_genes']),shared_pilot_core_proteins=len(shared),
            baseline_only_proteins=';'.join(sorted(set(bm['ordered_protein_ids'].split())-proteins)),
            pilot_only_proteins=';'.join(sorted(proteins-set(bm['ordered_protein_ids'].split()))),
            predicted_labels=bp['predicted_labels'],saccharide_score=float(bp['saccharide_probability']),
            **{f'P_at_{k}':float(bl[f'saccharide_fraction_top{k}']) for k in [5,10,20]}))
    write('tables/cluster7_overlapping_baseline_windows.tsv',overlap)
    result=dict(status='VERIFIED',locus_id=meta[i]['locus_id'],genome_accession=meta[i]['genome_id'],strategy='CORE_ONLY',n_input_genes=14,
        saccharide_score=float(score['saccharide_score']),predicted_labels=score['predicted_labels'],**{f'P_at_{k}':counts[k]/k for k in [5,10,20]},
        baseline_genome_accession_present=meta[i]['genome_id'] in {r['genome_accession'] for r in wm.values()},
        core_protein_accessions_overlapping_baseline_inputs=sorted(proteins&basalproteins),included_in_baseline_counts=False,cohort_membership_not_assumed=True,
        overlapping_baseline_windows=overlap,
        pooling='all-position mean of the 14 actual core genes; additive all-one mask',caveat='manually selected illustrative pilot; not verified glycan BGC ground truth')
    save('cluster7_benchmark.json',result);write('cluster7_input_provenance.tsv',[dict(path=str(PILOT/n),sha256=sha(PILOT/n)) for n in selected])
    return result
if __name__=='__main__':run()
