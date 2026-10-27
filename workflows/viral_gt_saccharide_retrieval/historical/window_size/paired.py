"""Paired target analysis, reference geometry, fixed-region recovery and taxonomy."""
from common import *
import time

def reference_metrics(vectors,metadata,condition):
    rmeta=read(REF/'mibig_matched_reference_metadata.tsv');ref=np.load(REF/'mibig_CORE_ONLY_embeddings.npy').astype(float)
    ref/=np.linalg.norm(ref,axis=1)[:,None];ids=np.array([r['BGC'] for r in rmeta]);sac=np.array(['Saccharide' in r['label'].split() for r in rmeta])
    assert len(ref)==2502 and int(sac.sum())==195 and len(set(ids))==2502
    centroid=ref[sac].mean(0);centroid/=np.linalg.norm(centroid)
    keys=list(vectors);q=np.stack([vectors[k] for k in keys]).astype(float);q/=np.linalg.norm(q,axis=1)[:,None]
    similarities=q@ref.T;result={};neighbors=[]
    for i,key in enumerate(keys):
        sim=similarities[i];order=np.lexsort((ids,-sim));first_sac=next(j for j in order if sac[j]);first_non=next(j for j in order if not sac[j])
        row=dict(nearest_saccharide_reference=str(ids[first_sac]),nearest_saccharide_cosine=float(sim[first_sac]),
            mean_saccharide_cosine=float(sim[sac].mean()),saccharide_centroid_cosine=float(q[i]@centroid),
            nearest_saccharide_minus_non_saccharide_cosine=float(sim[first_sac]-sim[first_non]))
        for k in [5,10,20]:
            row[f'p{k}']=float(sac[order[:k]].mean())
            row[f'top{k}_saccharide_reference_ids']=';'.join(str(ids[j]) for j in order[:k] if sac[j])
        row['top20_reference_ids']=';'.join(str(ids[j]) for j in order[:20]);result[key]=row
        neighbors.extend(dict(condition=condition,window_id=key,rank=rank+1,reference_BGC=str(ids[j]),reference_label=rmeta[j]['label'],
            saccharide=bool(sac[j]),cosine=float(sim[j])) for rank,j in enumerate(order[:20]))
    return result,neighbors

def signature(genes,gid):
    return extraction.digest(json.dumps([gid,[(g['start_1based'],g['end_1based'],g['strand'],g['sequence_sha256']) for g in genes]],separators=(',',':')))

def mean_stats(rows,subset,metric,clustering='genome_unit_id'):
    # Unique input contexts -> mean per fixed locus -> equal-weight genome means.
    dedup=list({r['context_signature15']:r for r in rows}.values());bylocus=collections.defaultdict(list);lg={}
    group_for={r['genome_unit_id']:r[clustering] for r in dedup}
    for r in dedup:bylocus[r['baseline_locus_id']].append(float(r['delta_'+metric]));lg[r['baseline_locus_id']]=r['genome_unit_id']
    bygenome=collections.defaultdict(list)
    for lid,vs in bylocus.items():bygenome[lg[lid]].append(float(np.mean(vs)))
    bygroup=collections.defaultdict(list)
    for gid,vs in bygenome.items():bygroup[group_for[gid]].append(float(np.mean(vs)))
    b=bootstrap([np.mean(vs) for vs in bygroup.values()])
    return dict(subset=subset,metric=metric,n_raw_pairs=len(rows),n_unique_context_pairs=len(dedup),n_fixed_loci=len(bylocus),n_resampling_groups=len(bygroup),
        resampling=clustering,window_mean_delta=float(np.mean([r['delta_'+metric] for r in rows])) if rows else '',
        genome_locus_balanced_delta=b['mean'],bootstrap95_low=b['low'],bootstrap95_high=b['high'],bootstrap_draws=5000)

def run():
    started=time.monotonic();preserve=verify_prior()
    oldp,oldloc,oldvec=load_results(BASE);newp,newloc,newvec=load_results(RUN)
    oldgenes=grouped_genes(BASE);newgenes=grouped_genes(RUN);mapping=read(RUN/'window15_window11_mapping.tsv')
    wm={r['window_id']:r for r in read(BASE/'viral_gt_windows.tsv')};w11={r['window_id']:r for r in read(RUN/'viral_gt_windows.tsv')}
    gmeta={r['genome_unit_id']:r for r in read(PRIOR/'tables/genome_relatedness_metadata.tsv')}
    gt_hashes={g['sequence_sha256'] for gs in oldgenes.values() for g in gs if int(g['relative_position'])==0}
    om,on=reference_metrics(oldvec,None,'15');nm,nn=reference_metrics(newvec,None,'11')
    write(OUT/'tables/reference_neighbors_top20.tsv',on+nn)
    pairs=[];allrows=[];propsrows=[]
    for m in mapping:
        old=m['window15_id'];new=m['window11_id'];t=gmeta[m['genome_unit_id']]
        removed=[g for g in oldgenes[old] if abs(int(g['relative_position']))>5]
        properties={tag:props(gs,gt_hashes) for tag,gs in [('15',oldgenes[old]),('11',newgenes[new]),('removed',removed)]}
        p15=oldp.get(old);p11=newp.get(new)
        state='paired' if p15 and p11 else 'newly_eligible_11' if p11 else 'lost_eligibility_11' if p15 else 'excluded_both'
        r=dict(m,pair_status=state,source_organism_labels=t['source_label'],taxid=t['taxid'],generic_source_label=t['generic_source_label'],
            **{k:t[k] for k in ['family','family_taxid','genus','genus_taxid']},
            context_signature15=signature(oldgenes[old],m['genome_unit_id']),context_signature11=signature(newgenes[new],m['genome_unit_id']),
            missing_positions15=wm[old]['missing_positions'],missing_positions11=w11[new]['missing_positions'],
            exclusion_reason15=wm[old]['inference_exclusion_reason'],exclusion_reason11=w11[new]['inference_exclusion_reason'],
            removed_protein_ids=';'.join(g['protein_id'] for g in removed),removed_products=';'.join(g['product'] for g in removed))
        for tag,pr in properties.items():r.update({tag+'_'+key:value for key,value in pr.items()})
        for suffix,p,loc,geometry in [('15',p15,oldloc.get(old),om.get(old)),('11',p11,newloc.get(new),nm.get(new))]:
            if p is None:continue
            r['saccharide_score'+suffix]=float(p['saccharide_probability']);r['classifier_positive'+suffix]=yes(p['saccharide_predicted'])
            r['predicted_labels'+suffix]=p['predicted_labels'];r['decision_threshold'+suffix]=float(p['threshold'])
            for label in LABELS:r['probability_'+label+'_'+suffix]=float(p['probability_'+label])
            for k in [5,10,20]:assert geometry[f'p{k}']==float(loc[f'saccharide_fraction_top{k}'])
            r.update({key+suffix:value for key,value in geometry.items()})
        if state=='paired':
            for metric in ['saccharide_score','classifier_positive','p5','p10','p20','nearest_saccharide_cosine','mean_saccharide_cosine',
                           'saccharide_centroid_cosine','nearest_saccharide_minus_non_saccharide_cosine']:
                r['delta_'+metric]=float(r[metric+'11'])-float(r[metric+'15'])
            for label in LABELS:r['delta_probability_'+label]=r['probability_'+label+'_11']-r['probability_'+label+'_15']
            r['embedding_cosine_15_vs_11']=float(np.dot(oldvec[old].astype(float),newvec[new].astype(float))/(np.linalg.norm(oldvec[old].astype(float))*np.linalg.norm(newvec[new].astype(float))))
            for k in [5,10,20]:
                a=set(filter(None,r[f'top{k}_saccharide_reference_ids15'].split(';')));b=set(filter(None,r[f'top{k}_saccharide_reference_ids11'].split(';')))
                r[f'top{k}_gained_saccharide_references']=';'.join(sorted(b-a));r[f'top{k}_lost_saccharide_references']=';'.join(sorted(a-b))
                r[f'top{k}_shared_saccharide_count']=len(a&b);r[f'top{k}_saccharide_Jaccard']=len(a&b)/len(a|b) if a|b else ''
            a=set(r['top20_reference_ids15'].split(';'));b=set(r['top20_reference_ids11'].split(';'));r['top20_neighbor_Jaccard']=len(a&b)/len(a|b)
            pairs.append(r)
        allrows.append(r);propsrows.append({k:v for k,v in r.items() if k in m or k.startswith(('15_','11_','removed_')) or k=='pair_status'})
    assert len(allrows)==757 and len(pairs)==len(oldp)==737 and not any(r['pair_status']=='lost_eligibility_11' for r in allrows)
    write(OUT/'tables/paired_window_comparison.tsv',pairs);write(OUT/'tables/all_occurrence_status.tsv',allrows)
    write(OUT/'tables/window_gene_content_properties.tsv',propsrows)
    write(OUT/'tables/newly_eligible_or_excluded.tsv',[r for r in allrows if r['pair_status']!='paired'])
    metrics=['saccharide_score','p5','p10','p20','nearest_saccharide_cosine','mean_saccharide_cosine','saccharide_centroid_cosine',
        'nearest_saccharide_minus_non_saccharide_cosine','classifier_positive']
    distributions=[];deltahist=[];transitions=[];direction=[]
    for metric in metrics:
        x=np.array([r['delta_'+metric] for r in pairs]);tol=1e-6 if metric=='saccharide_score' else 1e-12
        distributions.append(dict(metric=metric,n_pairs=len(x),n_increase=int((x>tol).sum()),n_decrease=int((x<-tol).sum()),n_unchanged_within_tolerance=int((np.abs(x)<=tol).sum()),
            mean=float(x.mean()),sd=float(x.std()),**{'quantile_'+str(q):float(np.quantile(x,q)) for q in [0,.05,.25,.5,.75,.95,1]},sign_tolerance=tol))
    for k in [5,10,20]:
        for n in range(-k,k+1):deltahist.append(dict(metric='p'+str(k),delta=n/k,n_pairs=sum(round(r['delta_p'+str(k)]*k)==n for r in pairs)))
    write(OUT/'tables/paired_change_distributions.tsv',distributions);write(OUT/'tables/paired_Pk_delta_frequencies.tsv',deltahist)
    scorebins=np.linspace(-1,1,41);hist,edges=np.histogram([r['delta_saccharide_score'] for r in pairs],bins=scorebins)
    write(OUT/'tables/paired_score_delta_histogram.tsv',[dict(lower=edges[i],upper=edges[i+1],n_pairs=int(h)) for i,h in enumerate(hist)])
    for metric in ['classifier_positive','p10_ge_0.5']:
        for a in [False,True]:
            for b in [False,True]:
                rr=[r for r in pairs if (r['classifier_positive15'] if metric=='classifier_positive' else r['p1015']>=.5)==a and
                    (r['classifier_positive11'] if metric=='classifier_positive' else r['p1011']>=.5)==b]
                transitions.append(dict(metric=metric,condition15=a,condition11=b,n_pairs=len(rr),n_loci=len({r['baseline_locus_id'] for r in rr}),n_genomes=len({r['genome_unit_id'] for r in rr})))
    for a in [-1,0,1]:
        for b in [-1,0,1]:
            def sign(x,tol):return 1 if x>tol else -1 if x<-tol else 0
            rr=[r for r in pairs if sign(r['delta_saccharide_score'],1e-6)==a and sign(r['delta_p10'],1e-12)==b]
            direction.append(dict(score_direction=a,P10_direction=b,n_pairs=len(rr)))
    write(OUT/'tables/paired_binary_transitions.tsv',transitions);write(OUT/'tables/score_retrieval_direction_concordance.tsv',direction)
    summaries=[]
    for label,rr in [('all_paired',pairs),('complete15_to_complete11',[r for r in pairs if int(r['n15'])==15 and int(r['n11'])==11]),
        ('incomplete15',[r for r in pairs if int(r['n15'])<15]),('zero_genes_removed',[r for r in pairs if int(r['n_removed'])==0]),
        ('exclude_generic_taxid2100421',[r for r in pairs if r['taxid']!='2100421'])]:
        for metric in metrics:summaries.append(mean_stats(rr,label,metric))
    assign={r['genome_unit_id']:r['similarity_cluster'] for r in read(PRIOR/'tables/genome_similarity_clusters.tsv') if float(r['Jaccard_cutoff'])==.5}
    clusterpairs=[dict(r,sequence_cluster=assign[r['genome_unit_id']]) for r in pairs]
    for metric in metrics:summaries.append(mean_stats(clusterpairs,'sequence_component_sensitivity',metric,'sequence_cluster'))
    write(OUT/'tables/genome_block_paired_statistics.tsv',summaries)
    candidates=analyze_candidates(allrows,oldp,oldloc,newp,newloc,mapping)
    pilot=analyze_pilot(allrows)
    save(OUT/'paired_summary.json',dict(status='PASS',n_pairs=len(pairs),occurrence_status=dict(collections.Counter(r['pair_status'] for r in allrows)),
        transitions=transitions,distributions=distributions,genome_statistics=summaries,candidate_counts=candidates,pilot=pilot,
        prior_files_unchanged=preserve,elapsed_seconds=time.monotonic()-started))
    print(json.dumps(dict(status='PASS',n_pairs=len(pairs),transitions=transitions,
        genome_statistics=[r for r in summaries if r['subset'] in {'all_paired','complete15_to_complete11'} and r['metric'] in {'saccharide_score','p10'}],
        candidate_counts=candidates,pilot=pilot),indent=2),flush=True)

def analyze_candidates(allrows,oldp,oldloc,newp,newloc,mapping):
    old2new={r['window15_id']:r['window11_id'] for r in mapping};new2old={v:k for k,v in old2new.items()}
    loci=read(BASE/'candidate_loci.tsv');genomes=read(BASE/'genome_units.tsv');wm={r['window_id']:r for r in read(BASE/'viral_gt_windows.tsv')}
    universes={'window':{w:[w] for w in wm},'locus':{r['locus_id']:r['window_ids'].split(';') for r in loci},'genome':{r['genome_unit_id']:r['window_ids'].split(';') for r in genomes}}
    sets={};statuses={};countrows=[];entityrows=[];comparison=[]
    for condition,pp,ll in [('15',oldp,oldloc),('11',newp,newloc)]:
        pred={w if condition=='15' else new2old[w]:r for w,r in pp.items()};local={w if condition=='15' else new2old[w]:r for w,r in ll.items()}
        A={w for w,r in pred.items() if yes(r['saccharide_predicted'])};B={w for w,r in pred.items() if w not in A and float(local[w]['saccharide_fraction_top10'])>=.5}
        for s,ids in [('A',A),('B',B),('C',A|B)]:
            flags={w:w in ids for w in pred};sets[(condition,s)]=ids
            for level,ee in universes.items():
                statuses[(condition,s,level)]={e:status(ww,flags) for e,ww in ee.items()};cc=collections.Counter(statuses[(condition,s,level)].values())
                countrows.append(dict(condition=condition,region_basis='frozen_15gene_loci',set=s,level=level,total=len(ee),any_inferred=sum(any(w in pred for w in ww) for ww in ee.values()),
                    fully_inferred=sum(all(w in pred for w in ww) for ww in ee.values()),candidates=cc['candidate'],not_selected=cc['not_selected'],unresolved=cc['unresolved']))
                entityrows.extend(dict(condition=condition,set=s,level=level,entity_id=e,status=v) for e,v in statuses[(condition,s,level)].items())
    byocc={r['window15_id']:r for r in allrows}
    for s in 'ABC':
        for level,ee in universes.items():
            for e,ww in ee.items():
                a=statuses[('15',s,level)][e];b=statuses[('11',s,level)][e]
                category='both' if a==b=='candidate' else '11_only_observed_gain' if b=='candidate' and a=='not_selected' else '11_candidate_baseline_unresolved' if b=='candidate' else '15_only_observed_loss' if a=='candidate' and b=='not_selected' else '15_candidate_11_unresolved' if a=='candidate' else 'neither_observed' if a==b=='not_selected' else 'unresolved_without_candidate'
                first=byocc[ww[0]]
                comparison.append(dict(set=s,level=level,entity_id=e,status15=a,status11=b,comparison=category,genome_unit_id=first['genome_unit_id'],
                    source_organism_labels=first['source_organism_labels'],genome_accessions=';'.join(sorted({wm[w]['genome_accession'] for w in ww})),
                    **{k:first[k] for k in ['taxid','family','family_taxid','genus','genus_taxid','generic_source_label']},
                    window15_ids=';'.join(ww),qualifying_window11_ids=';'.join(old2new[w] for w in ww if w in sets[('11',s)])))
    write(OUT/'tables/candidate_counts_by_condition.tsv',countrows);write(OUT/'tables/entity_condition_status.tsv',entityrows);write(OUT/'tables/candidate_recovery_comparison.tsv',comparison)
    taxonomy=[];breadth=[];novel=[]
    for view in ['recorded','non_generic_source_sensitivity']:
        for rank in ['family','genus']:
            groups=collections.defaultdict(list)
            for g in genomes:
                r=byocc[g['window_ids'].split(';')[0]]
                key=('Generic source; rank withheld','') if view!='recorded' and yes(r['generic_source_label']) else (r[rank] or 'Unresolved',r[rank+'_taxid'])
                groups[key].append(g)
            represented={condition:set() for condition in ['15','11']}
            for (name,tid),gg in sorted(groups.items()):
                gid={g['genome_unit_id'] for g in gg};ls={l['locus_id'] for l in loci if l['genome_unit_id'] in gid}
                tr=dict(view=view,rank=rank,taxon=name,rank_taxid=tid,total_genomes=len(gg),total_fixed_loci=len(ls))
                for c in ['15','11']:
                    tr['union_candidate_genomes'+c]=sum(statuses[(c,'C','genome')][g]=='candidate' for g in gid)
                    tr['union_candidate_loci'+c]=sum(statuses[(c,'C','locus')][l]=='candidate' for l in ls)
                    tr['unresolved_genomes'+c]=sum(statuses[(c,'C','genome')][g]=='unresolved' for g in gid)
                    tr['any_inferred_genomes'+c]=sum(any(w in (oldp if c=='15' else {new2old[w] for w in newp}) for w in g['window_ids'].split(';')) for g in gg)
                    if tid and tr['union_candidate_genomes'+c]:represented[c].add(tid)
                tr['newly_represented_at_11']=bool(tid and tr['union_candidate_genomes11'] and not tr['union_candidate_genomes15'])
                tr['represented_only_at_15']=bool(tid and tr['union_candidate_genomes15'] and not tr['union_candidate_genomes11'])
                taxonomy.append(tr)
                if tr['newly_represented_at_11'] or tr['represented_only_at_15']:novel.append(tr)
            for c in ['15','11']:
                selected={g['genome_unit_id'] for g in genomes if statuses[(c,'C','genome')][g['genome_unit_id']]=='candidate'}
                known={g['genome_unit_id'] for (name,tid),gg in groups.items() if tid for g in gg}
                breadth.append(dict(view=view,rank=rank,condition=c,distinct_union_taxa=len(represented[c]),candidate_genomes=len(selected),
                    rank_resolved_candidate_genomes=len(selected&known),rank_unresolved_or_withheld_candidate_genomes=len(selected-known)))
    write(OUT/'tables/union_taxonomy_comparison.tsv',taxonomy);write(OUT/'tables/union_taxonomic_breadth.tsv',breadth)
    write(OUT/'tables/new_or_lost_union_taxa.tsv',novel,list(taxonomy[0]))
    # Native 11-gene connected components are secondary: shortening can split a baseline region.
    windows11=read(RUN/'viral_gt_windows.tsv');genes11=read(RUN/'viral_gt_window_genes.tsv')
    native,_=baseline_analysis.merge_loci(windows11,genes11,newp);_,newgenomes=baseline_analysis.group_genomes(windows11,newp)
    write(RUN/'native_11gene_candidate_loci.tsv',native);write(RUN/'genome_units.tsv',newgenomes)
    crosswalk=[]
    for l in native:
        parent={byocc[new2old[w]]['baseline_locus_id'] for w in l['window_ids'].split(';')};assert len(parent)==1
        crosswalk.append(dict(native11_locus_id=l['locus_id'],baseline15_locus_id=next(iter(parent)),genome_unit_id=l['genome_unit_id'],n_windows=l['n_windows']))
    write(OUT/'tables/native11_to_fixed15_locus_crosswalk.tsv',crosswalk)
    nativecounts=[]
    for s in 'ABC':
        flags={w:(new2old[w] in sets[('11',s)]) for w in newp};cc=collections.Counter(status(l['window_ids'].split(';'),flags) for l in native)
        nativecounts.append(dict(condition='11',region_basis='native_11gene_connected_overlap',set=s,total_loci=len(native),candidates=cc['candidate'],not_selected=cc['not_selected'],unresolved=cc['unresolved']))
    write(OUT/'tables/native11_locus_counts.tsv',nativecounts)
    return countrows

def analyze_pilot(rows):
    r=next(r for r in rows if r['anchor_protein_id']=='ARF08732.1' and r['genome_accession']=='KY684083.1')
    historical=json.loads((RESCUE/'cluster7_benchmark.json').read_text())
    result=[]
    for c in ['15','11']:
        result.append(dict(condition=c,input_definition='systematic GT-centered +/-'+('7' if c=='15' else '5'),n_genes=int(r['n'+c]),
            window_id=r['window15_id' if c=='15' else 'window11_id'],saccharide_score=r['saccharide_score'+c],predicted_labels=r['predicted_labels'+c],
            p5=r['p5'+c],p10=r['p10'+c],p20=r['p20'+c],part_of_systematic_cohort=True))
    result.append(dict(condition='14',input_definition='manual ARF08726-ARF08739 core',n_genes=14,window_id=historical['locus_id'],
        saccharide_score=historical['saccharide_score'],predicted_labels=historical['predicted_labels'],p5=historical['P_at_5'],p10=historical['P_at_10'],p20=historical['P_at_20'],part_of_systematic_cohort=False))
    write(OUT/'tables/Catovirus_11_14_15_comparison.tsv',result)
    return result
if __name__=='__main__':run()
