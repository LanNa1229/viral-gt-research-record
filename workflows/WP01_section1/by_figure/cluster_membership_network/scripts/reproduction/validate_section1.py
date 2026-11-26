#!/usr/bin/env python3
"""AUDIT-ONLY portable validation of frozen Section 1 numerical evidence.

This new utility is not an original analysis script. It executes no Foldseek
searches. Default checks require Python standard library only; --bootstrap
additionally requires NumPy. Original input TSV identities are supplied in
figure_inputs.json metadata and the release manifests.
"""
import argparse,bisect,collections,json,math,statistics
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--source-data',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--bootstrap',action='store_true');a=p.parse_args()
D=json.loads(a.source_data.read_text(encoding='utf-8'));checks=[]
def test(name,value,expected,tol=0):
    ok=math.isclose(value,expected,rel_tol=tol,abs_tol=tol) if isinstance(value,(int,float)) and isinstance(expected,(int,float)) else value==expected
    checks.append(dict(check=name,observed=value,expected=expected,status='VERIFIED' if ok else 'BLOCKED'))
sizes=collections.Counter(r['representative_id'] for r in D['membership']);ids={r['protein_id'] for r in D['membership']}
test('unique baseline entries',len(ids),1727);test('clusters',len(sizes),761);test('multi-member clusters',sum(n>1 for n in sizes.values()),115);test('singleton clusters',sum(n==1 for n in sizes.values()),646);test('dominant cluster',sizes['A0A3Q8Q3U2'],594)
base={r['protein_id'] for r in D['membership'] if r['representative_id']=='A0A3Q8Q3U2'};new={r['protein_id'] for r in D['spike_membership']}-ids
spike={r['protein_id'] for r in D['spike_membership'] if r['representative_id']=='A0A3Q8Q3U2'};gt=new-{'7SP6_B','7SP6_C'}
test('GT entries',len(gt),7);test('distinct control proteins',len({'L230' if x in ['Q5UQC3_model','Q5UQC3_apo_model'] else x for x in gt}),6);test('GT recovery',len(gt&spike),7);test('spike cluster',len(spike),483);test('original giant retained',len(base&spike),439);test('other baseline recruited',len((spike&ids)-base),37)
cc=collections.Counter(r['panel_C_class'] for r in D['classes']);test('representatives',len({r['query_id'] for r in D['classes']}),114)
for label,n in [('GT-like',4),('non-GT',23),('ambiguous/unresolved',87)]:test(label,cc[label],n)
test('relevant original batch alignments in final hit audit',sum(int(r['n_total_hits_merged']) for r in D['classes']),29359)
test('representatives with original selected batch matches',sum(int(r['n_batches_with_hits'])>0 for r in D['classes']),112)
coh=D['coherence'];pairs=D['pairs'];test('small clusters',len(coh),63);test('small-cluster members',sum(int(r['cluster_size']) for r in coh),185);test('unordered pairs',len({tuple(sorted([r['protein_1'],r['protein_2']])) for r in pairs}),383)
support=collections.Counter()
for r in pairs:support[r['cluster_representative']]+=int(min(float(r['qtmscore']),float(r['ttmscore']))>=.5 and min(float(r['qcov']),float(r['tcov']))>=.5)
test('supported pairs',sum(support.values()),84);test('clusters with no supported pair',sum(support[r['cluster_representative']]==0 for r in coh),54)
master=D['master'];single=[r for r in master if r['cluster_group']=='singleton'];test('QC singletons',len(single),646)
for short,conf,n in [(True,False,328),(True,True,120),(False,False,105),(False,True,93)]:test(f'stratum short={short};confident={conf}',sum((float(r['analysis_sequence_length'])<122)==short and (float(r['median_plddt'])>=70)==conf for r in single),n)
def arr(metric,group):return [float(r[metric]) for r in master if (r['cluster_group']!='singleton' if group=='all_multimember' else r['cluster_group']==group) and r[metric] not in ('','NA','nan')]
def mw(x,y):
    sy=sorted(y);u=sum(bisect.bisect_left(sy,v)+(bisect.bisect_right(sy,v)-bisect.bisect_left(sy,v))/2 for v in x)
    n,m=len(x),len(y);N=n+m;t=collections.Counter(x+y);var=n*m/12*(N+1-sum(k**3-k for k in t.values())/(N*(N-1)))
    return u,min(1,math.erfc((max(u,n*m-u)-n*m/2-.5)/math.sqrt(2*var))),2*u/(n*m)-1
stats=[]
for r in D['statistics']:
    if r['test']!='Mann-Whitney U':continue
    g,h=r['groups_compared'].split(' vs ');u,pval,delta=mw(arr(r['metric'],g),arr(r['metric'],h));test(r['metric']+' '+g+'/'+h+' U',u,float(r['statistic']));test(r['metric']+' '+g+'/'+h+' delta',delta,float(r['effect_size']),1e-12)
    stats.append(dict(original=r,p=pval,delta=delta))
order=sorted(range(len(stats)),key=lambda i:stats[i]['p']);running=1.
for rank in range(len(order),0,-1):
    i=order[rank-1];running=min(running,stats[i]['p']*len(order)/rank);stats[i]['adjusted_p']=running
test('BH family size',len(stats),24)
for r in stats:test(r['original']['metric']+' '+r['original']['groups_compared']+' BH relative error',abs(r['adjusted_p']/float(r['original']['adjusted_p_value'])-1),0,1e-10)
bootstrap_status='NOT APPLICABLE: optional flag not requested'
if a.bootstrap:
    import numpy as np
    for m in ['fasta_length','median_plddt','frac_plddt_lt70']:
        x=np.array(arr(m,'singleton'));y=np.array(arr(m,'other_multimember'));rng=np.random.default_rng(20260712);ds=[]
        for _ in range(500):
            u=x[rng.integers(len(x),size=len(x))];v=np.sort(y[rng.integers(len(y),size=len(y))]);ds.append(float(np.sum(np.searchsorted(v,u,'left')+np.searchsorted(v,u,'right')-len(v))/(len(u)*len(v))))
        low,high=np.quantile(ds,[.025,.975]);r=next(r['original'] for r in stats if r['original']['metric']==m and r['original']['groups_compared']=='singleton vs other_multimember')
        test(m+' bootstrap lower',float(low),float(r['effect_ci_lower']),1e-12);test(m+' bootstrap upper',float(high),float(r['effect_ci_upper']),1e-12)
    bootstrap_status='VERIFIED:500resamples;seed20260712;three displayed intervals'
result={'status':'VERIFIED' if all(r['status']=='VERIFIED' for r in checks) else 'BLOCKED','checks':checks,'bootstrap':bootstrap_status,'scope':'frozen processed numerical evidence; not an end-to-end Foldseek or initial-search rerun','historical_gaps':['original MGAT2 web-search export (non-blocking for verified downstream baseline)','final PanelC merge/classifier source','exact historical execution hash/environment for some scripts'],'non_blocking_scope_limitation':'Separate truncated-core files unavailable on another computer; no search/rerun or independent equivalence claim; author decision 2026-10-24'}
a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,indent=2));print(json.dumps({'status':result['status'],'checks':len(checks),'failed':[r['check'] for r in checks if r['status']!='VERIFIED']}))
raise SystemExit(0 if result['status']=='VERIFIED' else 1)
