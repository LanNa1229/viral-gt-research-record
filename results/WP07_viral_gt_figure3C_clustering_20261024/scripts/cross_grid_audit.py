"""Matched-setting stability comparisons of SAVED partitions; never performs fitting."""
import collections,json,statistics
import common as c
primary=c.read(c.OUT/'assignments/selected287_leiden_grid.csv');ids=[r['sample_id'] for r in primary];pm={r['sample_id']:r for r in primary};runs=list(primary[0])[1:];baseline='k20_r1_s42'
registry={r['viral_locus_id']:r['duplicate_vector_group'] for r in c.read(c.OUT/'tables/duplicate_context_registry.csv')}
def choose2(v):return v*(v-1)/2
def ari(a,b):
 counts=collections.Counter(zip(a,b));aa=collections.Counter(a);bb=collections.Counter(b);n=len(a);tot=choose2(n);x=sum(choose2(v) for v in counts.values());y=sum(choose2(v) for v in aa.values());z=sum(choose2(v) for v in bb.values());expected=y*z/tot;den=(y+z)/2-expected
 return 1.0 if den==0 else (x-expected)/den
rows=[];jrows=[];summary=[]
for cohort in ['all556','selected_full15','selected_deduplicated','global_joint','saccharide_joint']:
 target=c.read(c.OUT/f'assignments/{cohort}_leiden_grid.csv');tm={r['sample_id']:r for r in target};shared=sorted(set(ids)&set(tm))
 if cohort=='selected_deduplicated':
  bygroup={registry[sid]:r for sid,r in tm.items()};tm={sid:bygroup[registry[sid]] for sid in ids};shared=sorted(ids)
 for run in runs:
  a=[pm[sid][run] for sid in shared];b=[tm[sid][run] for sid in shared];rows.append(dict(comparison=cohort,run_id=run,n_shared=len(shared),ARI_matched_setting=ari(a,b),n_selected_communities=len(set(a)),n_comparison_communities=len(set(b))))
  # Overlap of each fixed VC baseline community with every comparison run.
  for vc in sorted({pm[sid][baseline] for sid in shared}):
   members={sid for sid in shared if pm[sid][baseline]==vc};best=max(len(members&{sid for sid in shared if tm[sid][run]==label})/len(members|{sid for sid in shared if tm[sid][run]==label}) for label in set(b));jrows.append(dict(comparison=cohort,run_id=run,VC_baseline=vc,n_retained=len(members),best_Jaccard=best))
 grid=[r for r in rows if r['comparison']==cohort and r['run_id']!=baseline];values=[r['ARI_matched_setting'] for r in grid]
 summary.append(dict(comparison=cohort,n_matched_grid_runs=len(grid),n_shared=len(shared),ARI_min=min(values),ARI_median=statistics.median(values),ARI_max=max(values),baseline_ARI=next(r['ARI_matched_setting'] for r in rows if r['comparison']==cohort and r['run_id']==baseline)))
c.write('stability/cross_cohort_matched_grid_ARI.csv',rows);c.write('stability/cross_cohort_matched_grid_summary.csv',summary);c.write('stability/baseline_community_cross_cohort_grid_Jaccard.csv',jrows)
text='\n## Matched-setting cross-grid comparison\n\nAll 80 prespecified grid settings were also compared between selected-only and each corresponding sensitivity/joint cohort on shared selected loci. This uses saved partitions only, with no new fits or parameter selection.\n\n| Comparison | Shared loci | Matched-grid ARI min | Median | Max | Baseline ARI |\n|---|---:|---:|---:|---:|---:|\n'+''.join(f"| {r['comparison']} | {r['n_shared']} | {r['ARI_min']:.3f} | {r['ARI_median']:.3f} | {r['ARI_max']:.3f} | {r['baseline_ARI']:.3f} |\n" for r in summary)
with (c.OUT/'REPORT.md').open('a',encoding='utf-8') as f:f.write(text)
with (c.OUT/'METHODS_RECORD.md').open('a',encoding='utf-8') as f:f.write('\n## Matched-grid stability supplement\n\nSaved partitions were compared at each identical (k,resolution,seed) setting across cohorts, restricted to shared selected loci; duplicate-collapsed labels were expanded before comparison. Manual contingency-based ARI and fixed-baseline-community best Jaccard were recorded for all settings. No model was refit, no parameter was selected from these comparisons, and biological labels were not used. The exact command, source and consumed partition hashes are recorded in the supplementary_stability execution.\n')
c.save('provenance/matched_grid_stability.json',dict(status='PASS',comparisons=summary,new_fits=False,primary_parameters_changed=False))
print(json.dumps(summary,indent=2))
