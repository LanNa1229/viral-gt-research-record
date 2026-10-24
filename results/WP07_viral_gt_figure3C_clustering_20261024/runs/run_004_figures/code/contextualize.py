"""Post-fit reference context, cross-cohort comparisons and confounder summaries."""
import collections,itertools
import common as c
import numpy as np
from sklearn.metrics import adjusted_rand_score
CLASSES=['Alkaloid','Terpene','NRP','Polyketide','RiPP','Saccharide','Other']
def compare(baseline_ids,baseline_labels,target_ids,target_labels,name):
 a=dict(zip(baseline_ids,baseline_labels));b=dict(zip(target_ids,target_labels));shared=sorted(set(a)&set(b));aa=np.array([a[s] for s in shared]);bb=np.array([b[s] for s in shared]);rows=[]
 for cluster in sorted(set(aa)):
  members={s for s in shared if a[s]==cluster};matches=[]
  for target in sorted(set(bb)):
   other={s for s in shared if b[s]==target};overlap=len(members&other);matches.append((overlap/len(members|other),target,overlap,len(other)))
  best=max(matches,key=lambda z:(z[0],z[1]));rows.append(dict(comparison=name,Leiden_cluster=cluster,n_retained_members=len(members),best_target_cluster=best[1],best_Jaccard=best[0],overlap=best[2],target_size_on_shared_loci=best[3],retained_fraction_in_best_target=best[2]/len(members)))
 return dict(comparison=name,n_shared=len(shared),ARI=float(adjusted_rand_score(aa,bb)),n_primary_communities_on_shared=len(set(aa)),n_target_communities_on_shared=len(set(bb))),rows

def summarize(results,meta,x,selected_idx,rm,ref,vm):
 selected=[meta[i] for i in selected_idx];ids=[r['viral_locus_id'] for r in selected];primary=results['selected287'];labels=primary['baseline'];hlabels,hprob=primary['hdb']['mcs15_ms5'];n=len(ids)
 sample_lookup={r['sample_id']:r for r in c.read(c.ROOT/c.CONFIG['previous_analysis']/'joint_umap_samples.csv')}
 lookup={r['viral_locus_id']:r for r in meta};index={r['viral_locus_id']:i for i,r in enumerate(meta)}
 # 320D reference comparison is post-fit; old metadata remains unchanged.
 def unit(z):
  z=z.astype(float);return z/np.linalg.norm(z,axis=1,keepdims=True)
 cosine=np.clip(unit(x[selected_idx])@unit(ref).T,-1,1);rids=np.array([r['BGC'] for r in rm]);sac=np.array(['Saccharide' in r['label'].split() for r in rm]);refrows=[];refstats=[]
 for i,r in enumerate(selected):
  order=np.lexsort((rids,-cosine[i]));si=order[sac[order]];ni=order[~sac[order]];p10=float(sac[order[:10]].mean())
  assert p10==float(r['P_at_10']) and rids[order[0]]==r['nearest_reference']
  assert abs(float(cosine[i,order[0]])-float(r['nearest_reference_cosine']))<1e-12
  refstats.append(dict(viral_locus_id=ids[i],Leiden_cluster=labels[i],nearest_saccharide_reference=rids[si[0]],nearest_saccharide_cosine=float(cosine[i,si[0]]),nearest_non_saccharide_reference=rids[ni[0]],nearest_non_saccharide_cosine=float(cosine[i,ni[0]]),saccharide_minus_non_saccharide_nearest_cosine=float(cosine[i,si[0]]-cosine[i,ni[0]]),mean_top10_saccharide_reference_cosine=float(np.mean(cosine[i,si[:10]]))))
  for rank,j in enumerate(order[:20],1):refrows.append(dict(viral_locus_id=ids[i],Leiden_cluster=labels[i],rank=rank,reference_id=rids[j],reference_class=rm[j]['label'],reference_saccharide=bool(sac[j]),cosine_320D=float(cosine[i,j])))
 c.write('tables/selected_reference_neighbors_top20.csv',refrows);c.write('tables/selected_saccharide_context.csv',refstats)
 # Baseline communities are never relabeled from reference composition.
 mappings=[];comparisonrows=[];cross=[]
 for name,target in results.items():
  if name=='selected287':continue
  target_ids=target['ids'];target_labels=target['baseline']
  if name=='selected_deduplicated':
   vector_to_label={lookup[sid]['duplicate_vector_group']:lab for sid,lab in zip(target_ids,target_labels)}
   target_ids=ids;target_labels=np.array([vector_to_label[r['duplicate_vector_group']] for r in selected])
  stat,jaccards=compare(ids,labels,target_ids,target_labels,name);comparisonrows.append(stat);mappings.extend(jaccards)
  b=dict(zip(target_ids,target_labels))
  for i,sid in enumerate(ids):
   if sid in b:cross.append(dict(comparison=name,viral_locus_id=sid,Leiden_cluster=labels[i],comparison_cluster=b[sid]))
 c.write('tables/cross_cohort_partition_agreement.csv',comparisonrows);c.write('tables/community_sensitivity_Jaccard.csv',mappings);c.write('tables/cross_cohort_locus_memberships.csv',cross)
 # Per-community stability against every prespecified primary grid partition.
 stability=[];community_flags={};sensitivity={(r['comparison'],r['Leiden_cluster']):r for r in mappings}
 for cluster in sorted(set(labels)):
  members=np.flatnonzero(labels==cluster);member_ids={ids[i] for i in members};values=[];seedvalues=[]
  for run,partition in primary['runs'].items():
   if run=='k20_r1_s42':continue
   best=max(len(member_ids&{ids[i] for i in range(n) if partition[i]==other})/len(member_ids|{ids[i] for i in range(n) if partition[i]==other}) for other in set(partition))
   values.append(best)
   if run.startswith('k20_r1_s'):seedvalues.append(best)
   stability.append(dict(Leiden_cluster=cluster,run_id=run,best_Jaccard=best,baseline_members=len(members),baseline_k_resolution=run.startswith('k20_r1_s')))
  flags=c.CONFIG['descriptive_stability_flags'];seedmed=float(np.median(seedvalues));gridmed=float(np.median(values));sens=[sensitivity.get((name,cluster)) for name in ['all556','selected_full15','selected_deduplicated']]
  assessed=all(row is not None for row in sens) and sens[1]['n_retained_members']>=flags['minimum_retained_full15_members']
  stable_seed=seedmed>=flags['seed_median_best_Jaccard_min'];stable_grid=gridmed>=flags['grid_median_best_Jaccard_min'];stable_sens=assessed and all(row['best_Jaccard']>=flags['all556_full15_dedup_best_Jaccard_min'] for row in sens)
  community_flags[cluster]=dict(seed_median_best_Jaccard=seedmed,grid_median_best_Jaccard=gridmed,grid_q10_best_Jaccard=float(np.quantile(values,.1)),seed_stability_flag=stable_seed,parameter_stability_flag=stable_grid,sensitivity_assessable=assessed,core_sensitivity_flag=stable_sens,operational_robust_flag=bool(stable_seed and stable_grid and stable_sens))
 c.write('stability/primary_community_all_runs_Jaccard.csv',stability)
 # Reference compositions: per-query NN edges, and joint-community membership separately.
 compositions=[]
 for cluster in sorted(set(labels)):
  neighbors=[r for r in refrows if r['Leiden_cluster']==cluster and r['rank']<=10];denom=len(neighbors)
  for label in CLASSES:
   count=sum(label in r['reference_class'].split() for r in neighbors)
   compositions.append(dict(Leiden_cluster=cluster,existing_reference_class=label,positive_neighbor_edges=count,total_top10_neighbor_edges=denom,fraction=count/denom,definition='marginal label membership; hybrid labels count in multiple classes; repeated query-neighbor edges retained'))
 c.write('tables/community_reference_neighbor_class_composition.csv',compositions)
 jointsummary=[];jointclasses=[]
 for name in ['global_joint','saccharide_joint']:
  target=results[name];mapping=dict(zip(target['ids'],target['baseline']));refmeta={r['BGC']:r for r in rm}
  for community in sorted(set(target['baseline'])):
   members=[sid for sid,lab in mapping.items() if lab==community];viral=[sid for sid in members if sid in lookup];refs=[sid for sid in members if sid in refmeta]
   jointsummary.append(dict(analysis=name,joint_community=community,n_viral=len(viral),n_reference=len(refs),reference_to_viral_ratio=len(refs)/len(viral) if viral else None,viral_fraction=len(viral)/len(members),saccharide_references=sum('Saccharide' in refmeta[sid]['label'].split() for sid in refs),viral_only_communities=';'.join(sorted({labels[ids.index(sid)] for sid in viral}))))
   for label in CLASSES:
    count=sum(label in refmeta[sid]['label'].split() for sid in refs);jointclasses.append(dict(analysis=name,joint_community=community,reference_class=label,n_reference_with_label=count,n_reference_total=len(refs),fraction=count/len(refs) if refs else None,hybrid_rule='marginal memberships; fractions can sum above 1'))
 c.write('tables/joint_community_membership_summary.csv',jointsummary);c.write('tables/joint_reference_class_composition.csv',jointclasses)
 # Full locus records preserve original metadata columns exactly.
 all_labels=dict(zip(results['all556']['ids'],results['all556']['baseline']));selected_labels=dict(zip(ids,labels));selected_hdb=dict(zip(ids,hlabels));selected_probs=dict(zip(ids,hprob));allrows=[];primaryrows=[]
 selected_vector_counts=collections.Counter(r['duplicate_vector_group'] for r in selected);selected_context_counts=collections.Counter(r['duplicate_context_group'] for r in selected)
 for r in meta:
  sid=r['viral_locus_id'];coord=sample_lookup[sid];row=dict(r,Leiden_cluster=selected_labels.get(sid,''),HDBSCAN_cluster=selected_hdb.get(sid,''),HDBSCAN_probability=float(selected_probs[sid]) if sid in selected_probs else '',all556_Leiden_cluster=all_labels[sid],UMAP1=coord['UMAP1'],UMAP2=coord['UMAP2'],selected_context_multiplicity=selected_context_counts.get(r['duplicate_context_group'],0),selected_vector_multiplicity=selected_vector_counts.get(r['duplicate_vector_group'],0))
  allrows.append(row)
  if sid in selected_labels:primaryrows.append(row)
 c.write('tables/per_locus_selected287_clustering.csv',primaryrows);c.write('tables/per_locus_all556_clustering.csv',allrows)
 c.write('tables/Leiden_HDBSCAN_contingency.csv',[dict(Leiden_cluster=a,HDBSCAN_cluster=b,n_loci=int(np.sum((labels==a)&(hlabels==b)))) for a in sorted(set(labels)) for b in sorted(set(hlabels))])
 summaries=[];genome_detail=[];distance_groups=[]
 for cluster in sorted(set(labels)):
  ix=np.flatnonzero(labels==cluster);rs=[selected[i] for i in ix];p=np.array([float(r['P_at_10']) for r in rs]);gc=collections.Counter(r['genome_unit_id'] for r in rs);cc=collections.Counter(r['genome_sketch_component'] for r in rs);hc=collections.Counter(hlabels[ix]);sac_stats=[refstats[i] for i in ix];j=community_flags[cluster]
  for key,count in sorted(gc.items()):genome_detail.append(dict(Leiden_cluster=cluster,genome_unit_id=key,n_loci=count,fraction=count/len(ix)))
  summary=dict(Leiden_cluster=cluster,n_loci=len(ix),n_genome_units=len(gc),largest_genome_loci=max(gc.values()),largest_genome_fraction=max(gc.values())/len(ix),n_existing_genome_sketch_components=len(cc),largest_relatedness_component_fraction=max(cc.values())/len(ix),full15_n=sum(int(r['representative_n_genes'])==15 for r in rs),full15_fraction=np.mean([int(r['representative_n_genes'])==15 for r in rs]),minimum_genes=min(int(r['representative_n_genes']) for r in rs),unique_ordered_contexts=len({r['duplicate_context_group'] for r in rs}),duplicate_context_fraction=np.mean([selected_context_counts[r['duplicate_context_group']]>1 for r in rs]),P10_mean=float(p.mean()),P10_median=float(np.median(p)),P10_q25=float(np.quantile(p,.25)),P10_q75=float(np.quantile(p,.75)),P10_min=float(p.min()),P10_max=float(p.max()),nearest_reference_saccharide_fraction=np.mean(['Saccharide' in r['nearest_reference_class'].split() for r in rs]),top10_reference_saccharide_fraction=np.mean([r['reference_saccharide'] for r in refrows if r['Leiden_cluster']==cluster and r['rank']<=10]),nearest_saccharide_cosine_median=float(np.median([r['nearest_saccharide_cosine'] for r in sac_stats])),nearest_saccharide_cosine_q25=float(np.quantile([r['nearest_saccharide_cosine'] for r in sac_stats],.25)),nearest_saccharide_cosine_q75=float(np.quantile([r['nearest_saccharide_cosine'] for r in sac_stats],.75)),HDBSCAN_noise_fraction=float(np.mean(hlabels[ix]=='Noise')),HDBSCAN_largest_group_fraction=max(hc.values())/len(ix),**j)
  for name in ['all556','selected_full15','selected_deduplicated','global_joint','saccharide_joint','genome_balanced','genome_component_balanced']:
   match=sensitivity.get((name,cluster));summary[name+'_best_Jaccard']=match['best_Jaccard'] if match else None
  summaries.append(summary)
  distances=primary['d'][np.ix_(ix,ix)];vals=distances[np.triu_indices(len(ix),1)]
  if len(vals):distance_groups.append(dict(Leiden_cluster=cluster,within_pair_count=len(vals),within_cosine_distance_median=float(np.median(vals)),within_cosine_distance_q25=float(np.quantile(vals,.25)),within_cosine_distance_q75=float(np.quantile(vals,.75))))
 c.write('tables/community_summary.csv',summaries);c.write('tables/community_genome_composition.csv',genome_detail);c.write('tables/community_distance_summary.csv',distance_groups)
 # Joint graph reference-degree context directly measures reference-poor neighborhoods at baseline.
 contact=[]
 for name in ['global_joint','saccharide_joint']:
  target=results[name];g=target['graph'];nodeidx={sid:i for i,sid in enumerate(target['ids'])}
  for sid in ids:
   neighbors=g.neighbors(nodeidx[sid]);refneighbors=[j for j in neighbors if target['ids'][j].startswith('BGC')]
   contact.append(dict(analysis=name,viral_locus_id=sid,Leiden_cluster=selected_labels[sid],joint_degree=len(neighbors),reference_neighbor_count=len(refneighbors),reference_neighbor_fraction=len(refneighbors)/len(neighbors) if neighbors else 0,zero_reference_neighbor=len(refneighbors)==0))
 c.write('tables/joint_graph_viral_reference_contacts.csv',contact)
 c.save('provenance/context_summary.json',dict(reference_P10_and_nearest_identity_reproduced=287,annotations_used_only_after_fits=True,restricted_saccharide_reference_composition_is_not_enrichment=True,community_flags_are_prespecified_operational_descriptions=True,duplicate_vector_and_context_equivalence_identical=True,selected_unique_vectors=len(results['selected_deduplicated']['ids'])))
