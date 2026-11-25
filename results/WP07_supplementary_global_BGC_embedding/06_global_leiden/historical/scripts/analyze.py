"""Prespecified clustering in original 320D space; all annotation interpretation is post-fit."""
import collections,json
import common as c
import numpy as np
import engine as e
import contextualize
OLD=c.ROOT/c.CONFIG['previous_analysis'];BASE=c.ROOT/'analysis/viral_gt_mining_20261017/runs/baseline_594_v2';REF=c.ROOT/'analysis/viral_glycan_benchmark_20261016'
meta=c.read(OLD/'tables/included_locus_representatives.csv');x=np.load(OLD/'tables/representative_15gene_embeddings.npy',allow_pickle=False)
assert x.shape==(556,320) and len(meta)==556 and np.isfinite(x).all()
vm={r['window_id']:r for r in c.read(BASE/'embedding_manifest.tsv','\t')};raw=np.load(BASE/'viral_bgc_embeddings.npy',allow_pickle=False)
assert np.array_equal(x,np.stack([raw[int(vm[r['representative_window_id']]['embedding_row'])] for r in meta]))
ids=[r['viral_locus_id'] for r in meta];selected_idx=np.array([i for i,r in enumerate(meta) if r['candidate_status']=='candidate']);assert len(selected_idx)==287,'STOP: primary selected count discrepancy'
selected=[meta[i] for i in selected_idx];sx=x[selected_idx];sids=[ids[i] for i in selected_idx]
ref=np.load(REF/'mibig_CORE_ONLY_embeddings.npy',allow_pickle=False);rm=c.read(REF/'mibig_matched_reference_metadata.tsv','\t');rids=[r['BGC'] for r in rm]
assert ref.shape==(2502,320) and [int(r['embedding_row']) for r in rm]==list(range(2502))
sac=np.array(['Saccharide' in r['label'].split() for r in rm]);assert sac.sum()==195
# Exact identities and context equivalence are determined independently of labels/scores.
all_vector_keys=[c.hashlib.sha256(v.tobytes()).hexdigest() for v in x]
all_context_keys=[c.hashlib.sha256(vm[r['representative_window_id']]['ordered_sequence_sha256'].encode()).hexdigest() for r in meta]
vector_members=collections.defaultdict(list);context_members=collections.defaultdict(list)
for i,r in enumerate(meta):vector_members[all_vector_keys[i]].append(i);context_members[all_context_keys[i]].append(i)
vector_groups={key:f'DV{i+1:03d}' for i,key in enumerate(sorted(vector_members))};context_groups={key:f'DC{i+1:03d}' for i,key in enumerate(sorted(context_members))}
for i,r in enumerate(meta):r.update(duplicate_vector_group=vector_groups[all_vector_keys[i]],duplicate_context_group=context_groups[all_context_keys[i]],duplicate_vector_multiplicity=len(vector_members[all_vector_keys[i]]),duplicate_context_multiplicity=len(context_members[all_context_keys[i]]))
seen=set();unique_idx=[]
for i in selected_idx:
 if all_vector_keys[i] not in seen:seen.add(all_vector_keys[i]);unique_idx.append(i)
# Exact-vector and ordered-context collapse must be the same equivalence relation here.
for i in selected_idx:
 for j in selected_idx:assert (all_vector_keys[i]==all_vector_keys[j])==(all_context_keys[i]==all_context_keys[j]),'Vector/context duplicate equivalences differ; separate sensitivity required.'
full_idx=np.array([i for i in selected_idx if int(meta[i]['representative_n_genes'])==15]);assert len(full_idx)==243
c.write('tables/duplicate_context_registry.csv',[dict(viral_locus_id=r['viral_locus_id'],representative_window_id=r['representative_window_id'],vector_sha256=all_vector_keys[i],ordered_context_sha256=all_context_keys[i],duplicate_vector_group=r['duplicate_vector_group'],duplicate_context_group=r['duplicate_context_group'],vector_multiplicity=r['duplicate_vector_multiplicity'],context_multiplicity=r['duplicate_context_multiplicity']) for i,r in enumerate(meta)])
c.write('tables/short_representatives_excluded_from_full15.csv',[r for r in selected if int(r['representative_n_genes'])<15])
c.write('tables/selected_287_input_order.csv',[dict(selected_row=i,original_556_row=int(index),viral_locus_id=ids[index],representative_window_id=meta[index]['representative_window_id']) for i,index in enumerate(selected_idx)])
cohort_specs=[('selected287',sx,sids,'VC',True),('all556',x,ids,'AV',True),('selected_full15',x[full_idx],[ids[i] for i in full_idx],'FV',True),('selected_deduplicated',x[unique_idx],[ids[i] for i in unique_idx],'DV',True),('global_joint',np.concatenate((ref,sx)),rids+sids,'GC',False),('saccharide_joint',np.concatenate((ref[sac],sx)),[rids[i] for i in np.flatnonzero(sac)]+sids,'SC',False)]
assert len(cohort_specs[-2][2])==2789 and len(cohort_specs[-1][2])==482
results={name:e.run_cohort(name,vectors,nodeids,prefix,hd) for name,vectors,nodeids,prefix,hd in cohort_specs}
# Explicit genome/confounder sensitivity: retain all loci and alter only nuisance block weights.
related=c.read(c.ROOT/'analysis/viral_gt_interpretation_20261017/tables/genome_similarity_clusters.tsv','\t');cut=c.CONFIG['genome_sensitivity']['existing_sketch_Jaccard_component_cutoff']
components={r['genome_unit_id']:r['similarity_cluster'] for r in related if float(r['Jaccard_cutoff'])==cut}
assert all(r['genome_unit_id'] in components for r in selected)
for r in meta:r['genome_sketch_component']=components[r['genome_unit_id']]
for mode,key,prefix in [('genome_balanced','genome_unit_id','GB'),('genome_component_balanced','genome_sketch_component','CB')]:
 g=results['selected287']['graph'].copy();counts=collections.Counter(r[key] for r in selected)
 g.es['weight']=[float(edge['weight'])/np.sqrt(counts[selected[edge.source][key]]*counts[selected[edge.target][key]]) for edge in g.es]
 runs={};details=[]
 for seed in [0,1,2,3,4,42]:
  label,info=e.leiden(g,sids,1.0,seed,prefix);runs[f'k20_r1_s{seed}']=label;details.append(dict(seed=seed,**info))
 results[mode]=dict(name=mode,ids=sids,baseline=runs['k20_r1_s42'],runs=runs,graph=g)
 c.save(f'provenance/{mode}_iterations.json',details);c.write(f'assignments/{mode}_leiden_seeds.csv',[dict(sample_id=sid,**{run:labels[i] for run,labels in runs.items()}) for i,sid in enumerate(sids)])
contextualize.summarize(results,meta,x,selected_idx,rm,ref,vm)
summary={name:dict(n_samples=len(r['ids']),n_communities=len(set(r['baseline'])),community_sizes=dict(collections.Counter(r['baseline'])),silhouette=r.get('silhouette'),HDBSCAN=r.get('hdb_summary',[])) for name,r in results.items()}
c.save('RESULTS.json',summary);c.save('provenance/scientific_execution.json',dict(status='PASS',input_alignment_verified=True,original_representatives_unchanged=True,source_features='320D saved embeddings only; metadata used only for cohort restriction, documented confounder sensitivities and post-fit summaries',umap_refitted=False,embedding_inference_run=False,cluster_parameters_selected_using_labels=False,configuration=c.CONFIG,cohorts={name:len(r['ids']) for name,r in results.items()}))
print('SCIENTIFIC ANALYSIS COMPLETE',json.dumps({name:dict(n=r['n_samples'],communities=r['n_communities'],sizes=r['community_sizes']) for name,r in summary.items()}),flush=True)
