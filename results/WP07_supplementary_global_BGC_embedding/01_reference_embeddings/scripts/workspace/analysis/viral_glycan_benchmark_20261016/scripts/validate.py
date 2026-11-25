"""Independent checks of sources, input edge cases, tensors and retrieval tables."""
import argparse
import json
import numpy as np
from benchmark import OUT,ROOT,RECOVERY,read,tsv,sha,require,save_json,strategies,LABELS


def main(final=False):
    checks=[]
    def check(name,condition,detail=''):
        checks.append(dict(check=name,passed=bool(condition),detail=detail))
        require(condition,name+': '+detail)
    expected={
        ROOT/'BGC_train_dataset_classify.csv':'78dab471e8ef7dcdec16f6b4eb82ef15ef01b8bf900743c14238d26186a1d7f9',
        ROOT/'analysis/bgcprophet_reference/data/mibig_12510_embeddings.npy':'9a747f4a73e4d76074c95e6d2e8fbeb329e688e7ed27bda11d8a5ed8efb7d038',
        ROOT/'analysis/bgcprophet_reference/data/mibig_12510_probabilities.npy':'2d7bd10c71cd396b9da43313b11b1f55a9425a696dac969d0b69d14c85f025e0',
        RECOVERY/'model_cache/hub/checkpoints/esm2_t6_8M_UR50D.pt':'46f002a9870c9bdecd0ea887acb1f9a38a6b561e8f8bf8a6990b679b9d31b928',
        RECOVERY/'lmdb_snapshot/data.mdb':'ff1cb70d4df7132a6463c0b34fac56bced97ea81b383a7818475c7c15fd7062a',
        ROOT/'recovered_historical_scripts/extract_bgc_embeddings.py':'33e6a9eb5d07e7634b3033d96667c450be1b3b7da9dcf8417600dbdd4ff3f616',
    }
    integrity=[]
    for path,digest in expected.items():
        actual=sha(path)
        check('historical_hash_'+path.name,actual==digest,str(path))
        integrity.append(dict(path=str(path),expected_historical_sha256=digest,actual_sha256=actual,match=actual==digest))
    tsv(OUT/'source_integrity.tsv',integrity)
    testids=[f'g{i}' for i in range(200)]
    check('short_genome_A_unavailable',strategies(testids[:127],testids[10:30])[0][2].startswith('unavailable'))
    check('128_genome_exact_window',strategies(testids[:128],testids[10:30])[0][1]==testids[:128])
    for core in [testids[:20],testids[90:110],testids[-20:]]:
        a,b,c=strategies(testids,core)
        check('A_authentic_contiguous_'+core[0],len(a[1])==128 and set(core)<=set(a[1]) and
              a[1]==testids[testids.index(a[1][0]):testids.index(a[1][0])+128])
        check('B_C_membership_'+core[0],b[1]==(core*7)[:128] and c[1]==core)
    check('wide_core_A_unavailable',strategies(testids,[testids[0],testids[150]])[0][2].startswith('unavailable'))
    meta=read(OUT/'mibig_matched_reference_metadata.tsv')
    b=np.load(OUT/'mibig_REPEAT_128_embeddings.npy',allow_pickle=False)
    c=np.load(OUT/'mibig_CORE_ONLY_embeddings.npy',allow_pickle=False)
    check('matched_reference_shapes',b.shape==c.shape==(2502,320))
    check('reference_finite',np.isfinite(b).all() and np.isfinite(c).all())
    selected=[i for i,m in enumerate(meta) if 128%int(m['core_gene_count'])==0]
    d=np.abs(b[selected]-c[selected])
    check('actual_MIBiG_uniform_repeat_equivalence',np.allclose(b[selected],c[selected],atol=1e-4,rtol=1e-4),
          f'n={len(selected)}; max_abs={d.max()}; atol=rtol=1e-4')
    tsv(OUT/'mibig_uniform_repetition_controls.tsv',[
        dict(BGC=meta[i]['BGC'],core_gene_count=meta[i]['core_gene_count'],max_abs=float(np.abs(b[i]-c[i]).max()),
             passed=bool(np.allclose(b[i],c[i],atol=1e-4,rtol=1e-4))) for i in selected])
    raw=np.load(ROOT/'analysis/bgcprophet_reference/data/mibig_12510_embeddings.npy',allow_pickle=False)
    rawmeta=read(ROOT/'analysis/bgcprophet_reference/data/mibig_12510_metadata.tsv')
    groups={m['BGC']:[] for m in meta}
    for i,m in enumerate(rawmeta):groups[m['ID'].rsplit('_',1)[0]].append(i)
    consensus=np.stack([raw[groups[m['BGC']]].astype(np.float64).mean(0) for m in meta])
    a=np.load(ROOT/'analysis/bgcprophet_reference/results/consensus_embeddings.npy',allow_pickle=False)
    check('original_consensus_arithmetic',np.array_equal(a,consensus))
    q=np.load(OUT/'viral_bgc_embeddings.npy',allow_pickle=False)
    qm=read(OUT/'embedding_metadata.tsv')
    check('viral_finite_and_shape',q.shape==(2,320) and np.isfinite(q).all())
    scores=read(OUT/'classification_scores.tsv')
    for r in scores:
        v=np.array([float(r['score_'+name]) for name in LABELS])
        threshold=max(.5,v[-1])
        predictions=[name for name,value in zip(LABELS,v) if value>=threshold]
        if len(predictions)>1 and 'Other' in predictions:predictions.remove('Other')
        check('threshold_'+r['strategy'],r['predicted_labels']==(' '.join(predictions) or 'Unknown'))
        check('rank_'+r['strategy'],int(r['saccharide_rank'])==1+sum(v>v[5]))
    if final:
        neighbors=read(OUT/'mibig_nearest_neighbors.tsv')
        summary=read(OUT/'strategy_comparison.tsv')
        ids=[m['BGC'] for m in meta]
        for i,m in enumerate(qm):
            ref=b if m['strategy']=='REPEAT_128' else c
            sims=[float(np.dot(q[i].astype(float),v.astype(float))/(np.linalg.norm(q[i].astype(float))*np.linalg.norm(v.astype(float)))) for v in ref]
            order=sorted(range(2502),key=lambda k:(-sims[k],ids[k]))[:20]
            rows=[r for r in neighbors if r['strategy']==m['strategy']]
            check('nearest_neighbor_order_'+m['strategy'],[r['mibig_bgc'] for r in rows]==[ids[k] for k in order])
            check('cosine_values_'+m['strategy'],np.allclose([float(r['cosine_similarity']) for r in rows],[sims[k] for k in order],atol=1e-12))
            s=next(r for r in summary if r['strategy']==m['strategy'])
            for k in [5,10,20]:
                fraction=sum('Saccharide' in meta[j]['label'].split() for j in order[:k])/k
                check(f'P{k}_'+m['strategy'],float(s[f'saccharide_P_at_{k}'])==fraction)
        for name in ['pca_comparison.png','umap_comparison.png']:
            check('figure_exists_'+name,(OUT/name).stat().st_size>10000)
    tsv(OUT/('final_validation.tsv' if final else 'validation_before_comparison.tsv'),checks)
    save_json(OUT/('final_validation.json' if final else 'validation_before_comparison.json'),
              dict(status='PASS',checks=len(checks),mibig_divisor_length_controls=len(selected),
                   mibig_divisor_max_abs=float(d.max()),stage='final' if final else 'before_comparison'))
    print('PASS',len(checks),'checks',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--final',action='store_true');a=p.parse_args();main(a.final)
