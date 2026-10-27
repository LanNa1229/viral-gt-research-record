"""Extract exact +/-5 CDS contexts, verify compatibility, then frozen inference."""
from common import *
from Bio import SeqIO
import time,shutil,contextlib,subprocess,importlib.metadata

def parse_sources(windows):
    bypath=collections.defaultdict(list);result={}
    for w in windows:bypath[w['source_path']].append(w)
    for i,(path,ww) in enumerate(sorted(bypath.items())):
        checksum=sha(path);assert all(w['source_sha256']==checksum for w in ww)
        wanted={w['genome_accession'] for w in ww}
        with open(path) as f:
            for record in SeqIO.parse(f,'genbank'):
                if record.id not in wanted:continue
                info,genes=extraction.parse_record(record);result[(path,record.id)]=(info,genes)
        assert all((path,a) in result for a in wanted)
        if i%100==0:print('Verified source annotations',i+1,'/',len(bypath),flush=True)
    return result

def shortened(old,members,newid,anchor_index,existing_control=False):
    selected=[g for g in members if abs(int(g['gene_index'])-anchor_index)<=5]
    lo=max(int(members[0]['gene_index']) if existing_control else 0,anchor_index-5)
    hi=min(int(members[-1]['gene_index'])+1 if existing_control else int(old['n_cds']),anchor_index+6)
    assert selected and [int(g['gene_index']) for g in selected]==list(range(lo,hi))
    center=next(g for g in selected if int(g['gene_index'])==anchor_index)
    missing=[i for i in range(-5,6) if not lo<=anchor_index+i<hi]
    reasons=sorted({g['translation_status'] for g in selected if g['translation_status'] not in {'verified','verified_partial'}})
    if old['viral_status']=='unverified':reasons.append('viral_lineage_unverified')
    row=dict(old,window_id=newid,baseline_window_id=old['window_id'],anchor_gene_index=anchor_index,anchor_gene_id=center['gene_id'],
        anchor_protein_id=center['protein_id'],anchor_start_1based=center['start_1based'],anchor_end_1based=center['end_1based'],anchor_strand=center['strand'],
        n_genes=len(selected),missing_positions=','.join(map(str,missing)),missing_upstream=sum(i<0 for i in missing),missing_downstream=sum(i>0 for i in missing),
        inference_eligible=str(not reasons).lower(),inference_exclusion_reason=';'.join(reasons),
        start_1based=min(int(g['start_1based']) for g in selected),end_1based=max(int(g['end_1based']) for g in selected))
    genes=[dict(g,window_id=newid,relative_position=int(g['gene_index'])-anchor_index,is_anchor=int(g['gene_index'])==anchor_index) for g in selected]
    return row,genes

def prepare_main():
    started=time.monotonic();preserve=verify_prior();RUN.mkdir(parents=True,exist_ok=True)
    anchors=read(BASE/'input_snapshot/baseline_594_anchors.tsv');assert len(anchors)==len({r['accession'] for r in anchors})==594
    assert all(r['representative']=='A0A3Q8Q3U2' for r in anchors)
    windows=read(BASE/'viral_gt_windows.tsv');oldgenes=grouped_genes(BASE);source=parse_sources(windows)
    rows=[];genes=[];mapping=[];removed=[];checks=[]
    loci={w:l for l in read(BASE/'candidate_loci.tsv') for w in l['window_ids'].split(';')}
    genomes={w:g for g in read(BASE/'genome_units.tsv') for w in g['window_ids'].split(';')}
    for old in windows:
        wid=old['window_id'];info,allgenes=source[(old['source_path'],old['genome_accession'])];idx=int(old['anchor_gene_index'])
        assert info['genome_sequence_sha256']==old['genome_sequence_sha256'] and str(info['n_cds'])==old['n_cds']
        assert allgenes[idx]['gene_id']==old['anchor_gene_id'] and allgenes[idx]['protein_id']==old['anchor_protein_id']
        # Verify the full archived 15-gene context against the original sorted all-CDS annotation index.
        expected=allgenes[max(0,idx-7):idx+8];archived=oldgenes[wid]
        assert len(expected)==len(archived)
        for a,b in zip(expected,archived):
            for key in ['gene_id','gene_index','source_feature_index','protein_id','sequence','sequence_sha256','start_1based','end_1based','strand','location','product','translation_status']:
                assert str(a[key])==str(b[key]),(wid,key)
        newid='vgt11_'+extraction.digest(json.dumps([old['genome_accession'],old['genome_sequence_sha256'],old['anchor_gene_id'],5,'coordinate_no_wrap_v1'],separators=(',',':')))[:24]
        row,gg=shortened(old,archived,newid,idx)
        retained=allgenes[max(0,idx-5):idx+6]
        assert [g['gene_id'] for g in gg]==[g['gene_id'] for g in retained]
        assert all(str(g['sequence'])==r['sequence'] for g,r in zip(gg,retained))
        row.update(baseline_locus_id=loci[wid]['locus_id'],genome_unit_id=genomes[wid]['genome_unit_id'])
        rows.append(row);genes+=gg
        rm=[g for g in archived if abs(int(g['relative_position']))>5]
        removed.extend(dict(g,baseline_window_id=wid,window11_id=newid,locus_id=loci[wid]['locus_id'],genome_unit_id=genomes[wid]['genome_unit_id']) for g in rm)
        mapping.append(dict(window15_id=wid,window11_id=newid,genome_unit_id=genomes[wid]['genome_unit_id'],baseline_locus_id=loci[wid]['locus_id'],
            genome_accession=old['genome_accession'],anchor_protein_id=old['anchor_protein_id'],anchor_gene_id=old['anchor_gene_id'],
            anchor_gene_index=idx,n15=len(archived),n11=len(gg),n_removed=len(rm),eligible15=yes(old['inference_eligible']),eligible11=yes(row['inference_eligible']),
            source_sha256=old['source_sha256'],extraction_verified_from_source_CDS_index=True))
        checks.append(dict(window15_id=wid,window11_id=newid,n_source_CDS=len(allgenes),anchor_index=idx,
            original_15gene_members_match_source=True,new_11gene_members_match_source=True,sequence_coordinates_annotations_unchanged=True,
            relative_positions=' '.join(str(g['relative_position']) for g in gg),missing_positions=row['missing_positions'],status='PASS'))
    assert len(rows)==len({r['window_id'] for r in rows})==len(mapping)==757
    write(RUN/'viral_gt_windows.tsv',rows);write(RUN/'viral_gt_window_genes.tsv',genes)
    write(RUN/'window15_window11_mapping.tsv',mapping);write(RUN/'extraction_validation.tsv',checks)
    write(OUT/'tables/removed_outer_genes.tsv',removed)
    for name in ['anchor_mapping.tsv','taxonomy_cache.tsv']:
        shutil.copy2(BASE/name,RUN/('baseline_'+name))
    shutil.copy2(BASE/'input_snapshot/baseline_594_anchors.tsv',RUN/'baseline_594_anchors.tsv')
    with (RUN/'viral_gt_window_proteins.faa').open('w') as f:
        for g in genes:
            if g['sequence']:f.write(f">{g['window_id']}|{int(g['relative_position']):+d}|{g['gene_id']} protein_id={g['protein_id']}\n{g['sequence']}\n")
    save(RUN/'extraction_provenance.json',dict(status='PASS',parameters=dict(upstream=5,downstream=5,CDS_policy='all annotated CDS sorted by coordinates/source index',
        circular_policy='no wrap, same as baseline',no_padding=True,no_repetition=True,no_overlap_merge_before_inference=True),
        original_anchors=594,occurrences=757,complete_windows=sum(int(r['n_genes'])==11 for r in rows),eligible_windows=sum(yes(r['inference_eligible']) for r in rows),
        n_gene_rows=len(genes),removed_gene_rows=len(removed),baseline_files_verified=preserve,source_files=len({r['source_path'] for r in rows}),
        baseline_inputs={n:sha(BASE/n) for n in ['viral_gt_windows.tsv','viral_gt_window_genes.tsv','input_snapshot/baseline_594_anchors.tsv']},
        extraction_script_sha256=sha(Path(__file__)),elapsed_seconds=time.monotonic()-started))
    print('11-gene extraction:',len(rows),'windows;',sum(yes(r['inference_eligible']) for r in rows),'eligible;',sum(int(r['n_genes'])==11 for r in rows),'complete',flush=True)

def seed_new_cache(runs):
    # Copy only required cached proteins; never write to an old cache.
    needed={g['sequence_sha256'] for run in runs for _,gg in inference.prepare_windows(run)[0] for g in gg}
    sources=[]
    for parent in [PRIMARY/'cache/esm',PRIOR/'cache/esm']:
        sources.extend(p.parent for p in parent.glob('*/settings.json'))
    directories=[];copies=[];found={}
    for source in sources:
        settings=json.loads((source/'settings.json').read_text())
        assert settings['esm_checkpoint_sha256']==inference.ESM_SHA and settings['esm_layer']==6 and settings['esm_max_residues']==1022
        assert settings['device']=='cpu' and settings['dtype']=='float32' and settings['esm_batch_size']==1
        dest=OUT/'cache/esm'/source.name;dest.mkdir(parents=True,exist_ok=True)
        if (dest/'settings.json').exists():assert json.loads((dest/'settings.json').read_text())==settings
        else:shutil.copy2(source/'settings.json',dest/'settings.json')
        n=0
        for key in sorted(needed):
            p=source/(key+'.npy')
            if not p.exists():continue
            vector=np.load(p,allow_pickle=False);assert inference.valid_vector(vector,(320,))
            if key in found:assert np.array_equal(found[key],vector)
            found[key]=vector;n+=1;d=dest/p.name
            if d.exists():assert np.array_equal(np.load(d),vector)
            else:shutil.copy2(p,d)
            copies.append(dict(sequence_sha256=key,source=str(p),source_vector_sha256=sha(p),destination=str(d)))
        directories.append(dict(source=str(source),settings_sha256=sha(source/'settings.json'),matching_required_proteins=n))
    write(OUT/'cache_reuse_provenance.tsv',copies)
    save(OUT/'cache_reuse_summary.json',dict(required_unique_sequences=len(needed),unique_existing_cache_hits=len(found),source_directories=directories,
        policy='read old caches, validate vector/settings identity, copy required vectors only into new downstream cache; historical pilot caches additionally seeded by frozen inference'))
    print('Required proteins',len(needed),'; validated existing cache hits',len(found),flush=True)

def gate():
    proof=inference.verify_reference(ROOT);frozen=inference.load_frozen(ROOT);torch=frozen.torch_setup();model=frozen.load_model(torch)
    original={k:v.detach().clone() for k,v in model.state_dict().items()}
    oldp,oldloc,oldvec=load_results(BASE);gg=grouped_genes(BASE)
    cache,_=inference.seed_cache(ROOT)
    for source in [PRIMARY/'cache/esm',PRIOR/'cache/esm']:
        for settings in source.glob('*/settings.json'):
            for wid in [next(iter(oldp)),next(w for w,r in oldp.items() if r['anchor_protein_id']=='ARF08732.1'),
                next(w for w,r in oldp.items() if yes(r['saccharide_predicted']) and int(r['n_genes'])==15),
                min(oldp,key=lambda w:int(oldp[w]['n_genes']))]:
                for g in gg[wid]:
                    path=settings.parent/(g['sequence_sha256']+'.npy')
                    if path.exists():cache[g['sequence_sha256']]=np.load(path)
    rows=[]
    chosen=list(dict.fromkeys([next(iter(oldp)),next(w for w,r in oldp.items() if r['anchor_protein_id']=='ARF08732.1'),
        next(w for w,r in oldp.items() if yes(r['saccharide_predicted']) and int(r['n_genes'])==15),min(oldp,key=lambda w:int(oldp[w]['n_genes']))]))
    for wid in chosen:
        score,vec,core=frozen.infer(torch,model,[cache[g['sequence_sha256']] for g in gg[wid]],[1]*len(gg[wid]),True)
        expected=np.array([float(oldp[wid]['probability_'+label]) for label in LABELS])
        ds=float(np.max(np.abs(score-expected)));dv=float(np.max(np.abs(vec-oldvec[wid])))
        assert ds<=1e-6 and dv<=1e-6 and np.array_equal(vec,core)
        assert frozen.predicted(score)[0]==oldp[wid]['predicted_labels']
        rows.append(dict(window_id=wid,n_genes=len(gg[wid]),score_max_abs_error=ds,embedding_max_abs_error=dv,
            actual_pool_equals_original_pool=True,labels_match=True,status='PASS'))
    assert all(torch.equal(original[k],v) for k,v in model.state_dict().items())
    write(RUN/'compatibility_regression.tsv',rows);save(RUN/'compatibility_gate.json',dict(status='PASS',reference=proof,
        baseline_regression_windows=len(rows),model_weights_unchanged=True,regression_tolerance=1e-6,all_scores_and_vector_coordinates_checked=True))
    print('Compatibility gate PASS:',len(rows),'baseline inputs reproduce frozen scores and embeddings',flush=True)

def infer_one(run):
    assert json.loads((RUN/'compatibility_gate.json').read_text())['status']=='PASS'
    with (run/'inference_execution.log').open('w') as log,contextlib.redirect_stdout(log):
        result=inference.run_inference(ROOT,run,esm_batch_size=1,with_localization=True,sequence_cache_root=OUT/'cache/esm')
    print(run.name,result['status'],result['n_success'],'inferred;',result['n_excluded'],'excluded;',
        result['encoding']['newly_encoded'],'new protein encodings;',round(result['elapsed_seconds'],1),'seconds',flush=True)

def prepare_controls():
    CTRL.mkdir(parents=True,exist_ok=True);oldwindows=[];oldgenes={};origin={}
    for name in ['controls','controls_GT_union']:
        directory=PRIOR/name
        for w in read(directory/'viral_gt_windows.tsv'):oldwindows.append(w);origin[w['window_id']]=name
        oldgenes.update(grouped_genes(directory))
    source=parse_sources(oldwindows);rows=[];genes=[];mapping=[]
    for w in oldwindows:
        wid=w['window_id'];old=oldgenes[wid];idx=next(int(g['gene_index']) for g in old if int(g['relative_position'])==0)
        sourcegenes=source[(w['source_path'],w['genome_accession'])][1]
        for g in old:
            a=sourcegenes[int(g['gene_index'])]
            assert all(str(a[k])==str(g[k]) for k in ['gene_id','protein_id','sequence_sha256','sequence','location','product'])
        cid='ctrl11_'+extraction.digest(wid+'|minus5_plus5_v1')[:24];row,gg=shortened(w,old,cid,idx,existing_control=True)
        rows.append(row);genes+=gg
        mapping.append(dict(control15_id=wid,control11_id=cid,origin=origin[wid],genome_unit_id=w['genome_unit_id'],control_type=w['control_type'],
            target_locus_id=w.get('target_locus_id',''),n15=len(old),n11=len(gg),overlaps_baseline=w.get('overlaps_baseline','False')))
    write(CTRL/'viral_gt_windows.tsv',rows);write(CTRL/'viral_gt_window_genes.tsv',genes);write(CTRL/'control15_control11_mapping.tsv',mapping)
    gt_hashes={g['sequence_sha256'] for vv in grouped_genes(BASE).values() for g in vv if int(g['relative_position'])==0}
    newtarget=grouped_genes(RUN);targetmap={r['window15_id']:r['window11_id'] for r in read(RUN/'window15_window11_mapping.tsv')};controlmap={r['control15_id']:r['control11_id'] for r in mapping}
    newcontrol=grouped_genes(CTRL);checks=[]
    for p in read(PRIOR/'tables/GT_union_matched_pairs.tsv'):
        t=props(newtarget[targetmap[p['target_window_id']]],gt_hashes);c=props(newcontrol[controlmap[p['control_window_id']]],gt_hashes)
        reasons=[]
        if t['n_genes']!=c['n_genes']:reasons.append('gene_count')
        if t['gt_proxy_count']!=c['gt_proxy_count'] or t['gt_proxy_count']<1:reasons.append('GT_proxy_count')
        if not .75<=c['total_aa']/max(1,t['total_aa'])<=1.25:reasons.append('AA_length')
        if t['n_long']!=c['n_long']:reasons.append('long_protein_count')
        if abs(t['n_reverse']-c['n_reverse'])>2:reasons.append('reverse_count')
        checks.append(dict(p,target11_id=targetmap[p['target_window_id']],control11_id=controlmap[p['control_window_id']],
            passes_11gene_matching=not reasons,failing_constraints=';'.join(reasons),
            **{'target_'+k:v for k,v in t.items()},**{'control_'+k:v for k,v in c.items()}))
    write(OUT/'tables/control_matching_after_shortening.tsv',checks)
    save(CTRL/'extraction_provenance.json',dict(status='PASS',n_windows=len(rows),n_genes=len(genes),
        control_counts=dict(collections.Counter(w['control_type'] for w in rows)),revalidated_GT_matched_pairs=sum(yes(r['passes_11gene_matching']) for r in checks),
        method='same existing random/matched 15-gene control center; keep relative positions -5..+5; no resampling or score-based selection',
        inputs={n:{f:sha(PRIOR/n/f) for f in ['viral_gt_windows.tsv','viral_gt_window_genes.tsv','viral_bgc_predictions.tsv','viral_bgc_embeddings.npy']} for n in ['controls','controls_GT_union']}))
    print('Control extraction:',len(rows),'windows;',sum(yes(r['passes_11gene_matching']) for r in checks),'strict pairs still matched at 11 genes',flush=True)

def execute():
    OUT.mkdir(exist_ok=True);prepare_main();seed_new_cache([RUN])
    def child(stage):subprocess.run([sys.executable,'-B',str(Path(__file__).resolve()),stage],check=True,cwd=ROOT)
    child('gate');child('infer-main')
    # Primary extraction/inference precede the background experiment.
    prepare_controls();seed_new_cache([RUN,CTRL]);child('infer-controls')
    save(OUT/'experiment_execution.json',dict(status='PASS',primary_run=str(RUN),control_run=str(CTRL),preservation_after_execution=verify_prior()))
    print('Execution complete; primary baseline and prior downstream artifacts unchanged',flush=True)
if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['execute','prepare','gate','infer-main','prepare-controls','infer-controls']);a=p.parse_args()
    {'execute':execute,'prepare':prepare_main,'gate':gate,'infer-main':lambda:infer_one(RUN),'prepare-controls':prepare_controls,'infer-controls':lambda:infer_one(CTRL)}[a.stage]()
