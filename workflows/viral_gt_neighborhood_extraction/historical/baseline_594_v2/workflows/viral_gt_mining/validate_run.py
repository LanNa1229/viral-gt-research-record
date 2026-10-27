"""Independent source-coordinate, alignment, finite-output and 320D retrieval checks."""
import argparse
import collections
import json
from pathlib import Path
from Bio import SeqIO
import numpy as np
from viral_gt_mining.extract import read_tsv, write_tsv, json_write, sha_file


def validate(root,run,pilot=False):
    root,run=Path(root),Path(run);checks=[]
    def check(name,value):
        checks.append(dict(check=name,passed=bool(value)))
        if not value:raise ValueError('Validation failed: '+name)
    windows=read_tsv(run/'viral_gt_windows.tsv');genes=read_tsv(run/'viral_gt_window_genes.tsv')
    grouped=collections.defaultdict(list)
    for g in genes:grouped[g['window_id']].append(g)
    originals={}
    comparison=[]
    for w in windows:
        path=w['source_path']
        if path not in originals:
            with open(path) as f:originals[path]={r.id:r for r in SeqIO.parse(f,'genbank')}
        record=originals[path][w['genome_accession']]
        cds=sorted([(i,f) for i,f in enumerate(record.features) if f.type=='CDS'],
                   key=lambda x:(int(x[1].location.start),int(x[1].location.end),x[0]))
        anchor=int(w['anchor_gene_index']);expected=list(range(max(0,anchor-7),min(len(cds),anchor+8)))
        actual=sorted(grouped[w['window_id']],key=lambda g:int(g['relative_position']))
        check(w['window_id']+':all_expected_CDSs',[int(g['gene_index']) for g in actual]==expected)
        check(w['window_id']+':offsets',[int(g['relative_position']) for g in actual]==[i-anchor for i in expected])
        check(w['window_id']+':missing',w['missing_positions']==','.join(str(i) for i in range(-7,8) if not 0<=anchor+i<len(cds)))
        for g,i in zip(actual,expected):
            feature=cds[i][1]
            check(w['window_id']+':source_gene:'+str(i),int(feature.location.start)+1==int(g['start_1based'])
                and int(feature.location.end)==int(g['end_1based']) and (feature.location.strand or 0)==int(g['strand'])
                and feature.qualifiers.get('protein_id',[''])[0]==g['protein_id']
                and ''.join(feature.qualifiers.get('translation',[''])[0].split()).upper()==g['sequence'])
        if pilot:
            ids=[g['protein_id'] for g in actual]
            old=([f'QBP{i:05d}.1' for i in range(5878,5898)] if w['genome_accession']=='MK016663.1'
                 else [f'ARF{i:05d}.1' for i in range(8726,8740)] if w['genome_accession']=='KY684083.1' else [])
            comparison.append(dict(window_id=w['window_id'],anchor=w['anchor_protein_id'],genome=w['genome_accession'],
                start_1based=w['start_1based'],end_1based=w['end_1based'],actual_genes=' '.join(ids),
                manual_genes=' '.join(old),same_inputs=ids==old,n_shared=len(set(ids)&set(old)),
                only_new=' '.join(x for x in ids if x not in old),only_manual=' '.join(x for x in old if x not in ids)))
    manifest=read_tsv(run/'embedding_manifest.tsv');pred=read_tsv(run/'viral_bgc_predictions.tsv')
    array=np.load(run/'viral_bgc_embeddings.npy',allow_pickle=False)
    check('embedding_dimensions',array.shape==(len(manifest),320))
    check('embedding_finite',np.isfinite(array).all())
    check('array_manifest_prediction_alignment',[r['window_id'] for r in manifest]==[r['window_id'] for r in pred])
    check('all_extraction_windows_accounted',set(r['window_id'] for r in read_tsv(run/'inference_qc.tsv'))==set(w['window_id'] for w in windows))
    local=read_tsv(run/'saccharide_localization.tsv')
    refs=np.load(root/'analysis/viral_glycan_benchmark_20261016/mibig_CORE_ONLY_embeddings.npy',allow_pickle=False).astype(float)
    meta=read_tsv(root/'analysis/viral_glycan_benchmark_20261016/mibig_matched_reference_metadata.tsv')
    check('reference_population',refs.shape==(2502,320) and len(meta)==2502)
    positive=np.array(['Saccharide' in r['label'].split() for r in meta]);ids=np.array([r['BGC'] for r in meta])
    for i,(p,m,l) in enumerate(zip(pred,manifest,local)):
        scores=np.array([float(p['probability_'+name]) for name in ['Alkaloid','Terpene','NRP','Polyketide','RiPP','Saccharide','Other']])
        check(m['window_id']+':finite_scores',np.isfinite(scores).all() and ((scores>=0)&(scores<=1)).all())
        check(m['window_id']+':row_id',int(m['embedding_row'])==i and l['window_id']==m['window_id'])
        q=array[i].astype(float);similarity=(refs@q)/(np.linalg.norm(refs,axis=1)*np.linalg.norm(q))
        order=np.lexsort((ids,-similarity))
        check(m['window_id']+':independent_nearest',ids[order[0]]==l['nearest_reference'])
        check(m['window_id']+':independent_P10',abs(positive[order[:10]].mean()-float(l['saccharide_fraction_top10']))<1e-12)
        check(m['window_id']+':independent_cosine',abs(similarity[order[0]]-float(l['nearest_reference_cosine']))<1e-10)
    integrity=[]
    for name in ['viral_glycan_benchmark_20261016','viral_glycan_benchmark_A0A0K0N728_20261016','viral_glycan_benchmark_cluster7_ARF08726_ARF08739_20261016']:
        base=root/'analysis'/name
        for r in read_tsv(base/'artifact_checksums.tsv'):
            path=base/r['path'];actual=sha_file(path)
            check('historical:'+str(path),actual==r['sha256'])
            integrity.append(dict(path=str(path),sha256=actual))
    if pilot:
        check('pilot_7_success',len(pred)==len(windows)==7)
        check('six_complete_one_boundary',sum(int(w['n_genes'])==15 for w in windows)==6 and sum(int(w['n_genes'])==11 for w in windows)==1)
        write_tsv(run/'manual_region_comparison.tsv',comparison)
    write_tsv(run/'independent_validation.tsv',checks)
    write_tsv(run/'historical_integrity_verified.tsv',integrity)
    result=dict(status='PASS',checks=len(checks),windows=len(windows),inferred=len(pred),historical_artifacts_unchanged=len(integrity))
    json_write(run/'validation.json',result);print(json.dumps(result,indent=2));return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',required=True);p.add_argument('--run',required=True);p.add_argument('--pilot',action='store_true')
    a=p.parse_args();validate(a.root,a.run,a.pilot)
