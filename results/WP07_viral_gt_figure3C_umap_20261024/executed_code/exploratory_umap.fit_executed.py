"""Frozen CORE_ONLY input audit and one joint exploratory UMAP. No new inference/clustering."""
import argparse, collections, csv, hashlib, importlib.metadata, json, os, platform, sys, time
from pathlib import Path
OUT=Path(__file__).resolve().parents[1]; ROOT=OUT.parents[1]
if os.name=='nt': sys.path.insert(0,str(ROOT/'analysis/bgcprophet_reference/.deps'))
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'): os.environ[k]='1'
os.environ['MPLCONFIGDIR']=str(OUT/'cache/matplotlib'); os.environ['NUMBA_CACHE_DIR']=str(OUT/'cache/numba')
import numpy as np
BASE=ROOT/'analysis/viral_gt_mining_20261017/runs/baseline_594_v2'
REF=ROOT/'analysis/viral_glycan_benchmark_20261016'; SIZE=ROOT/'analysis/viral_gt_window_size_20261017'
CLASSES=('Alkaloid','Terpene','NRP','Polyketide','RiPP','Saccharide','Other')
COLORS=dict(zip(CLASSES,('#AA79A2','#6E9F8B','#718EB5','#C58B70','#99A56E','#C8A44C','#999999')))
PARAMS=dict(n_components=2,n_neighbors=30,min_dist=.1,metric='cosine',random_state=42)
HASHES={}; CHECKS=collections.Counter()
def check(v,n):
    if not v: raise ValueError('QC failed: '+n)
    CHECKS[n]+=1
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''): h.update(b)
    return h.hexdigest()
def track(p,expected=None):
    p=Path(p); h=sha(p)
    if expected: check(h==expected,'frozen_hash')
    HASHES[p]=h; return p
def read(p):
    with track(p).open(encoding='utf-8-sig',newline='') as f: return list(csv.DictReader(f,delimiter='\t'))
def jr(p): return json.loads(track(p).read_text(encoding='utf-8'))
def save(n,v): (OUT/n).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def write(n,rows):
    p=OUT/n; p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,list(rows[0]),lineterminator='\n'); w.writeheader(); w.writerows(rows)
def yes(s): return str(s).lower()=='true'
def unit(x):
    x=np.asarray(x,dtype=np.float64); return x/np.linalg.norm(x,axis=1,keepdims=True)
def audit():
    p=jr(BASE/'inference_provenance.json'); s=p['settings']; gate=p['reference_compatibility_gate']
    check(p['status']==gate['status']=='PASS' and s['strategy']=='CORE_ONLY','compatible_strategy')
    check(s['padding']==s['repetition']=='none','no_padding_repetition')
    check(s['membership']=='all float32 ones; original additive bias','mask_policy')
    check(s['esm_layer']==6 and s['esm_max_residues']==1022 and s['esm_pooling']=='residue mean excluding BOS/EOS','ESM_policy')
    for n,h in gate['fixed_artifact_sha256'].items(): track(REF/n,h)
    for n,h in {**p['output_checksums'],**p['input_checksums']}.items(): track(BASE/n,h)
    for path,h in [(REF/'sources/classifier.pt',s['classifier_sha256']),(ROOT/'sources/upstream/bgc_prophet/train/classifier.py',s['classifier_source_sha256']),
        (ROOT/'sources/upstream/bgc_prophet/train/model.py',s['encoder_source_sha256']),
        (ROOT/'analysis/bgcprophet_sequence_recovery_20260912/model_cache/hub/checkpoints/esm2_t6_8M_UR50D.pt',s['esm_checkpoint_sha256'])]: track(path,h)
    proof=jr(REF/'reference_validation.json'); pilot=jr(REF/'pilot_validation.json')
    check(proof['status']=='PASS' and proof['original_regression_passed'],'reference_regression')
    check(pilot['checkpoint_sha256']==s['classifier_sha256'] and pilot['weights_unchanged'],'same_checkpoint')
    ref=np.load(track(REF/'mibig_CORE_ONLY_embeddings.npy'),allow_pickle=False); rm=read(REF/'mibig_matched_reference_metadata.tsv')
    viral=np.load(track(BASE/'viral_bgc_embeddings.npy'),allow_pickle=False); vm=read(BASE/'embedding_manifest.tsv')
    for x,meta,name in [(ref,rm,'reference'),(viral,vm,'viral')]:
        check(x.shape==(len(meta),320),name+'_dimensions'); check(np.isfinite(x).all() and (np.linalg.norm(x,axis=1)>0).all(),name+'_valid')
        check([int(r['embedding_row']) for r in meta]==list(range(len(x))),name+'_row_alignment')
    check(len({r['BGC'] for r in rm})==len(rm)==proof['n_biological_bgcs'],'unique_reference_ids')
    check(len({r['window_id'] for r in vm})==len(vm),'unique_window_embeddings')
    with track(ROOT/'BGC_train_dataset_classify.csv').open(encoding='utf-8-sig',newline='') as f: train=list(csv.DictReader(f))
    groups=collections.defaultdict(list)
    for r in train: groups[r['ID'].rsplit('_',1)[0]].append(r)
    check(set(groups)=={r['BGC'] for r in rm},'reference_population')
    for m in rm:
        g=groups[m['BGC']]; core=g[0]['sentence'].split()
        check(len(g)==int(m['n_source_augmentations'])==5 and {r['ID'].rsplit('_',1)[1] for r in g}=={'-1','0','1','2','3'},'five_augmentations')
        check(len({r['sentence'] for r in g})==1 and {r['labels'] for r in g}=={m['label']},'identical_reference_core_labels')
        check(len(core)==len(set(core))==int(m['core_gene_count']),'reference_core_count'); check(set(m['label'].split())<=set(CLASSES),'seven_existing_classes')
    windows=read(BASE/'viral_gt_windows.tsv'); wins={r['window_id']:r for r in windows}; genes=collections.defaultdict(list)
    for g in read(BASE/'viral_gt_window_genes.tsv'): genes[g['window_id']].append(g)
    anchors=collections.defaultdict(set)
    for r in read(BASE/'anchor_mapping.tsv'):
        if r['window_id']: anchors[r['window_id']].add(r['input_accession'])
    ap=track(ROOT/'analysis/viral_gt_mining_20261017/inputs/baseline_594_anchors.txt'); aids=set(ap.read_text().split()); mapped=set().union(*anchors.values())
    check(mapped<=aids,'original_anchor_cohort'); loci=read(BASE/'candidate_loci.tsv')
    cand={r['locus_id']:r for r in read(ROOT/'analysis/viral_gt_figureB_taxonomy_20261018/tables/locus_audit.tsv')}
    reg={r['entity_id']:r for r in read(SIZE/'tables/candidate_recovery_comparison.tsv') if r['level']=='locus' and r['set']=='C'}
    check(set(cand)==set(reg)=={l['locus_id'] for l in loci},'same_locus_universe'); w2l={}
    for l in loci:
        c=cand[l['locus_id']]; r=reg[l['locus_id']]
        check(c['window15_ids']==l['window_ids'] and c['genome_unit_id']==l['genome_unit_id'],'frozen_locus_crosswalk')
        check((c['status15'],c['status11'])==(r['status15'],r['status11']),'frozen_candidate_labels')
        expected='candidate' if 'candidate' in (c['status15'],c['status11']) else 'not_selected' if c['status15']==c['status11']=='not_selected' else 'unresolved'
        check(c['combined_status']==expected,'candidate_union_consistency')
        for wid in l['window_ids'].split(';'):
            check(wid in wins and wid not in w2l,'unique_window_locus'); w2l[wid]=l['locus_id']
    check(set(w2l)==set(wins)==set(genes),'complete_window_mapping')
    qc={r['window_id']:r for r in read(BASE/'inference_qc.tsv')}; vmap={r['window_id']:r for r in vm}; cross=[]
    for wid,w in sorted(wins.items()):
        gs=sorted(genes[wid],key=lambda g:int(g['relative_position'])); n=int(w['n_genes']); ai=int(w['anchor_gene_index'])
        check(len(gs)==n and 1<=n<=15,'actual_gene_count')
        check([int(g['gene_index']) for g in gs]==list(range(max(0,ai-7),min(int(w['n_cds']),ai+8))),'contiguous_source_genes')
        check(len({g['gene_id'] for g in gs})==n,'no_repeated_gene_padding')
        check(all(int(g['relative_position'])==int(g['gene_index'])-ai for g in gs),'relative_position_identity')
        ag=[g for g in gs if int(g['relative_position'])==0]
        check(len(ag)==1 and ag[0]['gene_id']==w['anchor_gene_id'] and ag[0]['protein_id']==w['anchor_protein_id'] and bool(anchors[wid]),'anchor_identity')
        if wid in vmap:
            m=vmap[wid]; check(yes(w['inference_eligible']) and qc[wid]['status']=='PASS','embedding_eligibility')
            check(m['membership_mask'].split()==['1']*n and int(m['n_genes'])==n,'mask_length')
            check(m['strategy']=='CORE_ONLY' and m['classifier_sha256']==s['classifier_sha256'] and m['esm_sha256']==s['esm_checkpoint_sha256'],'per_embedding_model')
            check(m['ordered_protein_ids'].split()==[g['protein_id'] for g in gs] and m['ordered_sequence_sha256'].split()==[g['sequence_sha256'] for g in gs],'embedding_gene_order')
            check(yes(qc[wid]['pooled_matches_classifier_input']) and yes(qc[wid]['actual_and_published_pooling_equal']),'recorded_pooling')
            for g in gs: check(hashlib.sha256(g['sequence'].encode()).hexdigest()==g['sequence_sha256'],'gene_sequence_hash')
            with np.load(track(BASE/'inference_cache/windows'/(m['input_cache_key']+'.npz')),allow_pickle=False) as z:
                check(str(z['input_key'].item())==m['input_cache_key'],'cache_identity')
                check(np.array_equal(z['embedding'],viral[int(m['embedding_row'])]) and np.array_equal(z['embedding'],z['membership_embedding']),'cached_vector_pooling')
        else: check(not yes(w['inference_eligible']) and qc[wid]['status']=='excluded','missing_reason')
        cross.append(dict(window_id=wid,viral_locus_id=w2l[wid],original_anchor_accessions=';'.join(sorted(anchors[wid])),anchor_protein_id=w['anchor_protein_id'],
            genome_accession=w['genome_accession'],source_path=w['source_path'],source_sha256=w['source_sha256'],start_1based=w['start_1based'],end_1based=w['end_1based'],n_genes=n,
            has_embedding=wid in vmap,embedding_row=vmap.get(wid,{}).get('embedding_row',''),missing_reason=w['inference_exclusion_reason'],candidate_status=cand[w2l[wid]]['combined_status']))
    write('tables/window_locus_anchor_crosswalk.csv',cross)
    info=dict(reference_shape=list(ref.shape),viral_shape=list(viral.shape),anchors=len(aids),mapped_anchors=len(mapped),original_windows=len(windows),
        complete_windows=sum(int(w['n_genes'])==15 for w in windows),excluded_windows=len(windows)-len(vm),baseline_loci=len(loci),
        exclusion_reasons=dict(collections.Counter(w['inference_exclusion_reason'] for w in windows if w['window_id'] not in vmap)),
        baseline_status_counts=dict(collections.Counter(c['combined_status'] for c in cand.values())),reference_saccharide=sum('Saccharide' in m['label'].split() for m in rm),
        reference_multilabel=sum(len(m['label'].split())>1 for m in rm),classifier_sha256=s['classifier_sha256'],esm_sha256=s['esm_checkpoint_sha256'])
    return ref,rm,viral,vm,wins,loci,cand,info

def represent(viral,vm,wins,loci,cand):
    vmap={m['window_id']:m for m in vm}; norm=unit(viral); registry=[]; scores=[]; reps=[]
    for l in sorted(loci,key=lambda r:r['locus_id']):
        ids=sorted(l['window_ids'].split(';')); available=[w for w in ids if w in vmap]
        row=dict(viral_locus_id=l['locus_id'],genome_unit_id=l['genome_unit_id'],genome_accession=l['genome_accession'],source_genome_accessions=l['source_genome_accessions'],
            taxid=l['taxid'],organism=l['organism'],n_original_windows=len(ids),n_usable_windows=len(available),candidate_status=cand[l['locus_id']]['combined_status'],
            status15=cand[l['locus_id']]['status15'],status11=cand[l['locus_id']]['status11'],all_window_ids=';'.join(ids),included=bool(available))
        if available:
            ix=[int(vmap[w]['embedding_row']) for w in available]; local=norm[ix]@norm[ix].T; np.fill_diagonal(local,0)
            means=local.sum(axis=1)/(len(ix)-1) if len(ix)>1 else np.ones(1)
            # Numerical ties within 1e-12 use lexicographic stable IDs.
            chosen=int(np.flatnonzero(np.isclose(means,means.max(),atol=1e-12,rtol=0))[0]); wid=available[chosen]
            row.update(representative_window_id=wid,representative_embedding_row=ix[chosen],representative_n_genes=int(wins[wid]['n_genes']),mean_within_locus_cosine=float(means[chosen]),exclusion_reason='')
            reps.append(viral[ix[chosen]])
            for i,w in enumerate(available): scores.append(dict(viral_locus_id=l['locus_id'],window_id=w,embedding_row=ix[i],mean_cosine_to_other_windows=float(means[i]),selected=i==chosen))
        else: row.update(representative_window_id='',representative_embedding_row='',representative_n_genes='',mean_within_locus_cosine='',exclusion_reason='no_usable_15gene_workflow_embedding')
        registry.append(row)
    write('tables/locus_registry.csv',registry); write('tables/medoid_window_scores.csv',scores)
    return np.stack(reps),[r for r in registry if r['included']],registry

def neighbors(ref,rm,reps,selected):
    ids=np.array([r['BGC'] for r in rm]); pos=np.array(['Saccharide' in r['label'].split() for r in rm]); cos=unit(reps)@unit(ref).T
    old={r['window_id']:r for r in read(BASE/'saccharide_localization.tsv')}; rows=[]
    for i,r in enumerate(selected):
        order=np.lexsort((ids,-cos[i])); p10=float(pos[order[:10]].mean()); o=old[r['representative_window_id']]
        check(p10==float(o['saccharide_fraction_top10']),'P10_reproduction'); check(ids[order[0]]==o['nearest_reference'],'nearest_reproduction')
        r.update(P_at_10=p10,nearest_reference=ids[order[0]],nearest_reference_class=rm[order[0]]['label'],nearest_reference_saccharide=bool(pos[order[0]]),nearest_reference_cosine=float(cos[i,order[0]]))
        for rank,j in enumerate(order[:20],1): rows.append(dict(viral_locus_id=r['viral_locus_id'],representative_window_id=r['representative_window_id'],rank=rank,reference_id=ids[j],reference_class=rm[j]['label'],reference_saccharide=bool(pos[j]),cosine_320D=float(cos[i,j])))
    write('tables/representative_neighbors_top20.csv',rows); summaries=[]; frequencies=[]
    for group,keep in [('selected',[r for r in selected if r['candidate_status']=='candidate']),('not_selected',[r for r in selected if r['candidate_status']=='not_selected']),
        ('unresolved',[r for r in selected if r['candidate_status']=='unresolved']),('other_including_unresolved',[r for r in selected if r['candidate_status']!='candidate'])]:
        values=np.array([r['P_at_10'] for r in keep]); qs=np.quantile(values,[0,.25,.5,.75,1]) if len(values) else [None]*5
        summaries.append(dict(group=group,n=len(values),mean=float(values.mean()) if len(values) else None,minimum=qs[0],q25=qs[1],median=qs[2],q75=qs[3],maximum=qs[4],
            fraction_P10_ge_0_5=float(np.mean(values>=.5)) if len(values) else None,nearest_saccharide=sum(r['nearest_reference_saccharide'] for r in keep)))
        for k in range(11): frequencies.append(dict(group=group,P_at_10=k/10,n=int(np.count_nonzero(values==k/10))))
    write('tables/P10_group_summary.csv',summaries); write('tables/P10_group_frequencies.csv',frequencies)
    lengths=[]
    for group,keep in [('complete_15',[r for r in selected if r['representative_n_genes']==15]),('boundary_short',[r for r in selected if r['representative_n_genes']<15])]:
        lengths.append(dict(group=group,n=len(keep),min_genes=min(r['representative_n_genes'] for r in keep),max_genes=max(r['representative_n_genes'] for r in keep),mean_P10=float(np.mean([r['P_at_10'] for r in keep])),median_P10=float(np.median([r['P_at_10'] for r in keep]))))
    write('tables/gene_count_diagnostic.csv',lengths)
    count=sum(r['nearest_reference_saccharide'] for r in selected)
    return dict(nearest_reference_saccharide_count=count,nearest_reference_saccharide_fraction=count/len(selected),P10_groups=summaries,gene_count_groups=lengths)

def plot(xy,rm,vr):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.path import Path as MPath
    from PIL import Image,ImageOps
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.labelsize':11,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.7,'svg.fonttype':'none'})
    nr=len(rm); rxy=xy[:nr]; vxy=xy[nr:]; pos=np.array(['Saccharide' in r['label'].split() for r in rm]); chosen=np.array([r['candidate_status']=='candidate' for r in vr])
    unresolved=sum(r['candidate_status']=='unresolved' for r in vr); pad=np.ptp(xy,axis=0)*.055; low=xy.min(axis=0)-pad; high=xy.max(axis=0)+pad
    names=[]
    for panel in 'ABC':
        fig=plt.figure(figsize=(8.5,5.7)); ax=fig.add_axes([.10,.17,.61,.70]); ax.set(xlabel='UMAP 1',ylabel='UMAP 2',xlim=(low[0],high[0]),ylim=(low[1],high[1])); ax.set_aspect('equal',adjustable='box'); ax.tick_params(labelsize=9,width=.7)
        if panel=='A':
            groups=collections.defaultdict(list)
            for i,r in enumerate(rm): groups[tuple(c for c in CLASSES if c in r['label'].split())].append(i)
            for labels,ix in groups.items():
                for k,label in enumerate(labels):
                    marker='o'
                    if len(labels)>1:
                        angles=np.linspace(2*np.pi*k/len(labels),2*np.pi*(k+1)/len(labels),20)
                        marker=MPath(np.vstack(([0,0],np.column_stack((np.cos(angles),np.sin(angles))),[0,0])))
                    ax.scatter(*rxy[ix].T,s=12,c=COLORS[label],marker=marker,linewidths=0,alpha=.78,zorder=1)
            ax.scatter(*vxy.T,s=19,c='#151515',marker='^',linewidths=.15,edgecolors='white',alpha=.9,zorder=4)
            counts={c:sum(c in r['label'].split() for r in rm) for c in CLASSES}
            handles=[Line2D([],[],marker='o',linestyle='',color=COLORS[c],markersize=6,label=f'{c} ({counts[c]:,})') for c in CLASSES]
            handles.append(Line2D([],[],marker='^',linestyle='',color='#151515',markersize=6,label=f'Viral loci ({len(vr):,})'))
            ax.legend(handles=handles,loc='upper left',bbox_to_anchor=(1.02,1),frameon=False,title='Existing class memberships',fontsize=9,title_fontsize=10,labelspacing=.8)
            title='Reference BGC classes and viral genomic loci'; note='Split-color circles retain multilabel references; class membership counts overlap.'; name='Figure_A_reference_classes'
        elif panel=='B':
            ax.scatter(*rxy[~pos].T,s=12,c='#D5D5D5',linewidths=0,alpha=.7,zorder=1,label=f'Other references ({(~pos).sum():,})')
            ax.scatter(*rxy[pos].T,s=19,c='#C9A23B',linewidths=0,alpha=.9,zorder=2,label=f'Saccharide references ({pos.sum():,})')
            ax.scatter(*vxy[~chosen].T,s=23,c='#397BB1',marker='^',linewidths=.2,edgecolors='white',zorder=4,label=f'Other viral loci ({(~chosen).sum():,})')
            ax.scatter(*vxy[chosen].T,s=25,c='#C74640',marker='^',linewidths=.2,edgecolors='white',zorder=5,label=f'Selected viral loci ({chosen.sum():,})')
            ax.legend(loc='upper left',bbox_to_anchor=(1.02,1),frameon=False,fontsize=9,labelspacing=1)
            title='Saccharide references and prior candidate status'; note=f'Prior 11/15-gene union labels retained. Other viral loci include {unresolved} unresolved loci.'; name='Figure_B_candidate_status'
        else:
            ax.scatter(*rxy.T,s=12,c='#D5D5D5',linewidths=0,alpha=.7,zorder=1)
            vals=np.array([r['P_at_10'] for r in vr]); order=np.argsort(vals,kind='stable')
            sc=ax.scatter(*vxy[order].T,s=26,c=vals[order],cmap='viridis',vmin=0,vmax=1,marker='^',linewidths=.2,edgecolors='white',zorder=4)
            handles=[Line2D([],[],marker='o',linestyle='',color='#CFCFCF',markersize=6,label=f'Reference BGCs ({nr:,})'),Line2D([],[],marker='^',linestyle='',color='#555555',markersize=6,label=f'Viral loci ({len(vr):,})')]
            ax.legend(handles=handles,loc='upper left',bbox_to_anchor=(1.02,1),frameon=False,fontsize=9)
            cax=fig.add_axes([.76,.28,.025,.36]); bar=fig.colorbar(sc,cax=cax,ticks=np.linspace(0,1,6)); bar.set_label('Saccharide P@10',fontsize=10,labelpad=8); bar.ax.tick_params(labelsize=9)
            title='Saccharide neighbors of representative viral windows'; note='P@10 uses the plotted representative window and cosine neighbors in the original 320D space.'; name='Figure_C_saccharide_P10'
        fig.text(.025,.945,panel,weight='bold',fontsize=16); fig.text(.10,.945,title,fontsize=12)
        fig.text(.10,.905,f'{nr:,} reference BGCs + {len(vr):,} viral loci | joint cosine UMAP | seed 42',color='#555555',fontsize=9)
        fig.text(.10,.075,note,fontsize=8.5,color='#444444'); fig.text(.10,.040,'Exploratory visualization; UMAP separation does not define biological clusters.',fontsize=8.5,color='#555555')
        for ext in ('png','svg'): fig.savefig(OUT/'figures'/(name+'.'+ext),dpi=600,facecolor='white')
        names.append(name); plt.close(fig)
    canvas=Image.new('RGB',(1100,2214),'white')
    for i,name in enumerate(names):
        with Image.open(OUT/'figures'/(name+'.png')) as im: canvas.paste(ImageOps.contain(im.convert('RGB'),(1100,738)),(0,i*738))
    canvas.save(OUT/'figures/preview_ABC.png')

def main():
    ap=argparse.ArgumentParser(description=__doc__); ap.add_argument('--audit-only',action='store_true'); args=ap.parse_args(); started=time.monotonic()
    for d in ('figures','tables','cache/numba','cache/matplotlib'): (OUT/d).mkdir(parents=True,exist_ok=True)
    ref,rm,viral,vm,wins,loci,cand,info=audit(); reps,selected,registry=represent(viral,vm,wins,loci,cand)
    info.update(included_loci=len(selected),excluded_loci=len(registry)-len(selected),included_status_counts=dict(collections.Counter(r['candidate_status'] for r in selected)),
        excluded_status_counts=dict(collections.Counter(r['candidate_status'] for r in registry if not r['included'])),multiple_original_windows=sum(r['n_original_windows']>1 for r in registry),
        multiple_usable_windows=sum(r['n_usable_windows']>1 for r in registry),representative_complete15=sum(r['representative_n_genes']==15 for r in selected),
        representative_min_genes=min(r['representative_n_genes'] for r in selected),representative_max_genes=max(r['representative_n_genes'] for r in selected),
        reference_core_min_genes=min(int(r['core_gene_count']) for r in rm),reference_core_max_genes=max(int(r['core_gene_count']) for r in rm),unique_representative_vectors=len(np.unique(reps,axis=0)))
    info.update(neighbors(ref,rm,reps,selected)); info['checks_passed']=dict(CHECKS)
    write('tables/included_locus_representatives.csv',selected); np.save(OUT/'tables/representative_15gene_embeddings.npy',reps); save('input_and_representation_audit.json',info)
    print(json.dumps(info,indent=2),flush=True)
    if args.audit_only: return
    print('Input gates passed. Fitting one joint UMAP.',flush=True)
    import umap
    # Raw original 320D vectors: no PCA, standardization, or normalization before UMAP.
    joint=np.concatenate((ref,reps),axis=0).astype(np.float64); model=umap.UMAP(**PARAMS,n_jobs=1); xy=model.fit_transform(joint)
    check(xy.shape==(len(ref)+len(reps),2) and np.isfinite(xy).all(),'joint_coordinates'); np.save(OUT/'tables/joint_umap_coordinates.npy',xy)
    rows=[]
    for m,coord in zip(rm,xy[:len(ref)]): rows.append(dict(sample_id=m['BGC'],sample_type='reference_BGC',reference_class=m['label'],viral_locus_id='',representative_window_id='',candidate_status='',P_at_10='',UMAP1=float(coord[0]),UMAP2=float(coord[1])))
    for r,coord in zip(selected,xy[len(ref):]): rows.append(dict(sample_id=r['viral_locus_id'],sample_type='viral_locus',reference_class='',viral_locus_id=r['viral_locus_id'],representative_window_id=r['representative_window_id'],candidate_status=r['candidate_status'],P_at_10=r['P_at_10'],UMAP1=float(coord[0]),UMAP2=float(coord[1])))
    write('joint_umap_samples.csv',rows); plot(xy,rm,selected)
    versions={p:importlib.metadata.version(p) for p in ('numpy','scipy','scikit-learn','umap-learn','numba','pynndescent','matplotlib','pillow')}
    integrity=[]
    for p,h in sorted(HASHES.items(),key=lambda x:str(x[0])):
        after=sha(p); check(after==h,'input_unchanged'); integrity.append(dict(path=str(p.relative_to(ROOT)),sha256_before=h,sha256_after=after,unchanged=True))
    write('tables/input_integrity.csv',integrity)
    provenance=dict(status='PASS',script_sha256=sha(__file__),umap_parameters=model.get_params(),software_versions=versions,python=sys.version,executable=sys.executable,platform=platform.platform(),
        normalization='Raw 320D vectors promoted to float64; no UMAP input normalization/scaling/PCA. Retrieval uses float64 L2 normalization.',
        input_order='Reference embedding_row order, then lexicographic locus ID',fit_count=1,labels_used_for_fit=False,clustering_performed=False,input_files_unchanged=len(integrity),checks_passed=dict(CHECKS),elapsed_seconds=time.monotonic()-started)
    save('execution_provenance.json',provenance); save('software_versions.json',versions)
    import pickle
    with (OUT/'joint_umap_model.pkl').open('wb') as f: pickle.dump(model,f,protocol=5)
    write('output_checksums.csv',[dict(path=str(p.relative_to(OUT)),sha256=sha(p)) for p in sorted(OUT.rglob('*')) if p.is_file() and 'cache' not in p.relative_to(OUT).parts and p.name!='output_checksums.csv'])
    print('Completed one joint UMAP, A/B/C PNG+SVG, CSV, medoid mappings and 320D diagnostics.',flush=True)
if __name__=='__main__': main()
