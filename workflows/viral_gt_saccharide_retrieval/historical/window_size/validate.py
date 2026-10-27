"""Independent extraction, score-head, retrieval and paired-accounting audit."""
from common import *
import csv,time,shutil

def audit():
    started=time.monotonic();preserved=verify_prior();checks=0
    def check(ok):
        nonlocal checks
        assert ok
        checks+=1
    oldw={r['window_id']:r for r in read(BASE/'viral_gt_windows.tsv')};neww={r['window_id']:r for r in read(RUN/'viral_gt_windows.tsv')}
    oldg=grouped_genes(BASE);newg=grouped_genes(RUN);mapping=read(RUN/'window15_window11_mapping.tsv')
    check(len(mapping)==757);check(set(r['window15_id'] for r in mapping)==set(oldw));check(set(r['window11_id'] for r in mapping)==set(neww))
    removed=collections.defaultdict(list)
    for r in read(OUT/'tables/removed_outer_genes.tsv'):removed[r['baseline_window_id']].append(r)
    for m in mapping:
        a=oldw[m['window15_id']];b=neww[m['window11_id']];aa=oldg[a['window_id']];bb=newg[b['window_id']];idx=int(a['anchor_gene_index'])
        expected=list(range(max(0,idx-5),min(int(a['n_cds']),idx+6)))
        check([int(g['gene_index']) for g in bb]==expected)
        check(int(b['anchor_gene_index'])==idx and a['anchor_gene_id']==b['anchor_gene_id'] and a['anchor_protein_id']==b['anchor_protein_id'])
        check(sum(yes(g['is_anchor']) for g in bb)==1)
        byindex={int(g['gene_index']):g for g in aa}
        for g in bb:
            old=byindex[int(g['gene_index'])]
            for key in ['sequence','sequence_sha256','protein_id','gene_id','location','location_parts','source_feature_index','start_1based','end_1based','strand','product','translation_status']:
                check(g[key]==old[key])
            check(int(g['relative_position'])==int(g['gene_index'])-idx)
        kept={g['gene_id'] for g in bb};rm=removed[a['window_id']]
        check({g['gene_id'] for g in rm}=={g['gene_id'] for g in aa}-kept)
        check(len(rm)+len(bb)==len(aa));check(all(abs(int(g['relative_position']))>5 for g in rm))
        missing=[i for i in range(-5,6) if idx+i<0 or idx+i>=int(a['n_cds'])]
        check(b['missing_positions']==','.join(map(str,missing)))
        check(int(b['n_genes'])==len(bb) and int(b['n_genes'])<=11)
        expected_eligible=all(g['translation_status'] in {'verified','verified_partial'} for g in bb) and b['viral_status']!='unverified'
        check(yes(b['inference_eligible'])==expected_eligible)
    sourcefiles={r['source_path']:r['source_sha256'] for r in oldw.values()}
    for path,h in sourcefiles.items():check(sha(path)==h)
    write(RUN/'source_annotation_checksums.tsv',[dict(path=p,sha256=h) for p,h in sorted(sourcefiles.items())])
    # Independently evaluate the frozen linear sigmoid head from saved pooled vectors.
    import torch
    state=torch.load(REF/'sources/classifier.pt',map_location='cpu',weights_only=True)
    weight=state['classifier.1.weight'].numpy().astype(float);bias=state['classifier.1.bias'].numpy().astype(float)
    refmeta=read(REF/'mibig_matched_reference_metadata.tsv');ref=np.load(REF/'mibig_CORE_ONLY_embeddings.npy').astype(float)
    ref/=np.linalg.norm(ref,axis=1)[:,None];ids=np.array([r['BGC'] for r in refmeta]);sac=np.array(['Saccharide' in r['label'].split() for r in refmeta])
    validation_rows=[];head_errors=[]
    for directory in [RUN,CTRL]:
        predictions,local,vectors=load_results(directory);manifest=read(directory/'embedding_manifest.tsv');genes=grouped_genes(directory)
        check(len(vectors)==(740 if directory==RUN else 1552))
        matrix=np.stack(list(vectors.values())).astype(float);norm=matrix/np.linalg.norm(matrix,axis=1)[:,None];similarities=norm@ref.T
        scores=1/(1+np.exp(-(matrix@weight.T+bias)))
        for i,m in enumerate(manifest):
            wid=m['window_id'];p=predictions[wid];gg=genes[wid];loc=local[wid]
            expected=np.array([float(p['probability_'+label]) for label in LABELS]);error=float(np.max(np.abs(scores[i]-expected)));head_errors.append(error)
            check(error<2e-6);check(np.isfinite(vectors[wid]).all() and vectors[wid].shape==(320,))
            check(m['ordered_sequence_sha256'].split()==[g['sequence_sha256'] for g in gg]);check(m['ordered_protein_ids'].split()==[g['protein_id'] for g in gg])
            check(m['membership_mask'].split()==['1']*len(gg));check(m['classifier_sha256']==inference.CLASSIFIER_SHA and m['esm_sha256']==inference.ESM_SHA)
            threshold=max(float(p['probability_Other']),.5);check(threshold==float(p['threshold']))
            check(yes(p['saccharide_predicted'])==(float(p['saccharide_probability'])>=threshold))
            order=np.lexsort((ids,-similarities[i]))
            for k in [5,10,20]:check(float(sac[order[:k]].mean())==float(loc[f'saccharide_fraction_top{k}']))
            check(str(ids[order[0]])==loc['nearest_reference'])
            with np.load(directory/'inference_cache/windows'/(m['input_cache_key']+'.npz'),allow_pickle=False) as cache:
                check(np.array_equal(cache['embedding'],vectors[wid]));check(np.array_equal(cache['embedding'],cache['membership_embedding']))
                check(np.array_equal(cache['scores'],expected.astype(np.float32)))
            validation_rows.append(dict(run=directory.name,window_id=wid,head_max_abs_error=error,n_genes=len(gg),pooling_valid=True,retrieval_valid=True,status='PASS'))
        provenance=json.loads((directory/'inference_provenance.json').read_text())
        for name,h in provenance['input_checksums'].items():check(sha(directory/name)==h)
        for name,h in provenance['output_checksums'].items():check(sha(directory/name)==h)
    write(OUT/'tables/independent_inference_validation.tsv',validation_rows)
    a,al,_=load_results(BASE);b,bl,_=load_results(RUN);old2new={m['window15_id']:m['window11_id'] for m in mapping}
    pairs=read(OUT/'tables/paired_window_comparison.tsv');check(len(pairs)==737)
    for r in pairs:
        aa=a[r['window15_id']];bb=b[r['window11_id']]
        check(r['window11_id']==old2new[r['window15_id']])
        for label in LABELS:check(abs(float(r['delta_probability_'+label])-(float(bb['probability_'+label])-float(aa['probability_'+label])))<1e-12)
        check(abs(float(r['delta_saccharide_score'])-(float(bb['saccharide_probability'])-float(aa['saccharide_probability'])))<1e-12)
        for k in [5,10,20]:check(abs(float(r['delta_p'+str(k)])-(float(bl[r['window11_id']][f'saccharide_fraction_top{k}'])-float(al[r['window15_id']][f'saccharide_fraction_top{k}'])))<1e-12)
        if int(r['n_removed'])==0:check(float(r['delta_saccharide_score'])==0 and float(r['delta_p10'])==0)
    universes={'window':{w:[w] for w in oldw},'locus':{r['locus_id']:r['window_ids'].split(';') for r in read(BASE/'candidate_loci.tsv')},
        'genome':{r['genome_unit_id']:r['window_ids'].split(';') for r in read(BASE/'genome_units.tsv')}}
    expectedsets={};counts=read(OUT/'tables/candidate_counts_by_condition.tsv')
    for c in ['15','11']:
        pred=a if c=='15' else {old:b[new] for old,new in old2new.items() if new in b}
        local=al if c=='15' else {old:bl[new] for old,new in old2new.items() if new in bl}
        aa={w for w,p in pred.items() if yes(p['saccharide_predicted'])};bb={w for w,p in pred.items() if w not in aa and float(local[w]['saccharide_fraction_top10'])>=.5}
        for s,selected in [('A',aa),('B',bb),('C',aa|bb)]:
            for level,ee in universes.items():
                classified={e:'candidate' if set(ww)&selected else 'not_selected' if all(w in pred for w in ww) else 'unresolved' for e,ww in ee.items()}
                expectedsets[(c,s,level)]=classified;cc=collections.Counter(classified.values())
                saved=next(r for r in counts if r['condition']==c and r['set']==s and r['level']==level)
                for name,column in [('candidate','candidates'),('not_selected','not_selected'),('unresolved','unresolved')]:check(cc[name]==int(saved[column]))
                check(sum(cc.values())==len(ee))
    for r in read(OUT/'tables/entity_condition_status.tsv'):check(r['status']==expectedsets[(r['condition'],r['set'],r['level'])][r['entity_id']])
    for r in read(OUT/'tables/candidate_recovery_comparison.tsv'):
        check(r['status15']==expectedsets[('15',r['set'],r['level'])][r['entity_id']]);check(r['status11']==expectedsets[('11',r['set'],r['level'])][r['entity_id']])
    for view in ['recorded','non_generic_source_sensitivity']:
        for rank in ['family','genus']:
            tx=[r for r in read(OUT/'tables/union_taxonomy_comparison.tsv') if r['view']==view and r['rank']==rank]
            check(sum(int(r['total_genomes']) for r in tx)==466)
            for c in ['15','11']:
                check(sum(int(r['union_candidate_genomes'+c]) for r in tx)==sum(v=='candidate' for v in expectedsets[(c,'C','genome')].values()))
                check(sum(int(r['union_candidate_loci'+c]) for r in tx)==sum(v=='candidate' for v in expectedsets[(c,'C','locus')].values()))
    # Target-control contrast arithmetic and finite-group bootstrap inputs.
    for r in read(OUT/'tables/background_genome_pair_values.tsv'):
        d=(float(r['target11'])-float(r['target15']))-(float(r['control11'])-float(r['control15']))
        check(abs(d-float(r['difference_in_differences']))<1e-12)
    for folder in [RUN,CTRL,OUT/'tables']:
        for path in folder.glob('*.tsv'):
            with path.open(encoding='utf-8-sig') as f:header=next(csv.reader(f,delimiter='\t'))
            check(len(header)==len(set(x.casefold() for x in header)))
    # Store exact frozen implementation references rather than modifying/copying the primary run.
    source_paths=[BASE/'implementation_snapshot/src/viral_gt_mining/inference.py',BASE/'implementation_snapshot/src/viral_gt_mining/extract.py',
        ROOT/'sources/upstream/bgc_prophet/train/classifier.py',ROOT/'sources/upstream/bgc_prophet/train/model.py',REF/'scripts/benchmark.py',
        REF/'sources/classifier.pt',ROOT/'analysis/bgcprophet_sequence_recovery_20260912/model_cache/hub/checkpoints/esm2_t6_8M_UR50D.pt',
        REF/'mibig_CORE_ONLY_embeddings.npy',REF/'mibig_matched_reference_metadata.tsv']
    write(RUN/'frozen_implementation_and_reference_checksums.tsv',[dict(path=str(p),sha256=sha(p),bytes=p.stat().st_size) for p in source_paths])
    validation=dict(status='PASS',independent_checks=checks,score_head_checks=len(head_errors),max_score_head_abs_error=max(head_errors),
        source_annotation_files_reverified=len(sourcefiles),paired_occurrences=len(pairs),baseline_occurrences=757,window11_inferred=740,control11_inferred=1552,
        compatibility_regression=json.loads((RUN/'compatibility_gate.json').read_text()),
        frozen_inputs_verified_before=preserved,frozen_inputs_verified_after=verify_prior(),elapsed_seconds=time.monotonic()-started)
    save(OUT/'validation.json',validation);save(RUN/'validation.json',validation)
    print(json.dumps({k:validation[k] for k in ['status','independent_checks','score_head_checks','max_score_head_abs_error','elapsed_seconds']},indent=2),flush=True)

def manifests():
    for folder in [RUN,CTRL,OUT]:
        name='artifact_checksums.tsv';rows=[]
        for p in sorted(folder.rglob('*')):
            if not p.is_file() or p.name==name or any(x in {'__pycache__','.git','cache'} for x in p.relative_to(folder).parts):continue
            rows.append(dict(path=p.relative_to(folder).as_posix(),bytes=p.stat().st_size,sha256=sha(p)))
        write(folder/name,rows)
    # Include all newly persisted protein cache vectors/settings in a separate complete inventory.
    cache=OUT/'cache/esm'
    if cache.exists():write(OUT/'protein_cache_checksums.tsv',[dict(path=p.relative_to(OUT).as_posix(),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(cache.rglob('*')) if p.is_file()])
    # Add that inventory to the downstream manifest after it is written.
    rows=read(OUT/'artifact_checksums.tsv');p=OUT/'protein_cache_checksums.tsv'
    rows=[r for r in rows if r['path']!=p.name]+[dict(path=p.name,bytes=p.stat().st_size,sha256=sha(p))];write(OUT/'artifact_checksums.tsv',rows)
    for folder in [RUN,CTRL,OUT]:
        rows=read(folder/'artifact_checksums.tsv')
        assert all(sha(folder/r['path'])==r['sha256'] for r in rows)
        print(folder.name,'checksums PASS',len(rows),flush=True)
if __name__=='__main__':
    if '--manifests-only' in sys.argv:manifests()
    else:audit()
