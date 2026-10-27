#!/usr/bin/env python3
"""Pinned historical Step2 functions with portable resource and cold-cache plumbing."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import argparse,csv,json,hashlib,importlib.util,importlib.metadata,urllib.request,tarfile,io,time,platform
PACKAGE=Path(__file__).resolve().parents[1]
def digest(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def rows(p):
    with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def save(p,x):Path(p).write_text(json.dumps(x,indent=2)+'\n',encoding='utf8')
def write(p,items):
    with Path(p).open('w',encoding='utf8',newline='') as f:
        w=csv.DictWriter(f,items[0].keys(),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(items)
def registry():return json.loads((PACKAGE/'resources.json').read_text())
def checked(p,h):
    p=Path(p)
    if not p.is_file() or digest(p)!=h:raise ValueError('Missing or checksum-mismatched resource: '+str(p))
    return p
def download(url):
    with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'Viral-GT-Step2/1.0'}),timeout=120) as r:return r.read()
def put(p,data,h):
    p=Path(p)
    if hashlib.sha256(data).hexdigest()!=h:raise ValueError('Downloaded checksum mismatch: '+str(p))
    p.parent.mkdir(parents=True,exist_ok=True)
    if p.exists():checked(p,h);return
    with p.open('xb') as f:f.write(data)
def prepare(cache,step1_archive=None,step2_data=None,offline=False):
    import zipfile
    cache=Path(cache).resolve();reg=registry();cache.mkdir(parents=True,exist_ok=True)
    for item in reg['files']:
        p=cache/item['cache_path']
        if p.exists():checked(p,item['sha256']);continue
        if offline:raise ValueError('Offline resource missing: '+str(p))
        data=download(item['url'])
        if item['kind']=='tar_member':
            if hashlib.sha256(data).hexdigest()!=item['archive_sha256']:raise ValueError('Model archive checksum mismatch')
            with tarfile.open(fileobj=io.BytesIO(data),mode='r:gz') as t:
                members=[m for m in t.getmembers() if m.isfile() and Path(m.name).name==item['member_basename']]
                if len(members)!=1:raise ValueError('Ambiguous model archive member')
                data=t.extractfile(members[0]).read()
        put(p,data,item['sha256'])
    missing=[x for x in reg['step1']['inputs'] if not (cache/x['cache_path']).is_file()]
    if missing:
        archive=Path(step1_archive) if step1_archive else cache/'downloads'/reg['step1']['archive_name']
        if not archive.exists():
            if offline:raise ValueError('Supply the existing Step1 archive')
            put(archive,download(reg['step1']['archive_url']),reg['step1']['archive_sha256'])
        checked(archive,reg['step1']['archive_sha256'])
        with zipfile.ZipFile(archive) as z:
            for x in missing:put(cache/x['cache_path'],z.read(x['member']),x['sha256'])
    for x in reg['step1']['inputs']:checked(cache/x['cache_path'],x['sha256'])
    for x in reg['step2_outputs']:
        p=cache/x['cache_path']
        if p.exists():checked(p,x['sha256']);continue
        if not step2_data:raise ValueError('Supply --step2-data extracted dataset directory')
        source=checked(Path(step2_data)/x['dataset_path'],x['sha256']);put(p,source.read_bytes(),x['sha256'])
    return dict(status='PASS',resources_verified=len(reg['files'])+len(reg['step1']['inputs'])+len(reg['step2_outputs']),historical_pilot_caches_required=False)
def verify(cache):
    import numpy as np
    cache=Path(cache).resolve();reg=registry()
    for x in reg['files']+reg['step1']['inputs']+reg['step2_outputs']:checked(cache/x['cache_path'],x['sha256'])
    result=[]
    for size,n,positive in [('15',737,246),('11',740,265)]:
        d=cache/'runs'/size;windows=rows(d/'viral_gt_windows.tsv');genes=rows(d/'viral_gt_window_genes.tsv');pred=rows(d/'viral_bgc_predictions.tsv');man=rows(d/'embedding_manifest.tsv');qc=rows(d/'inference_qc.tsv');v=np.load(d/'viral_bgc_embeddings.npy',mmap_mode='r',allow_pickle=False)
        if len(windows)!=757 or len(pred)!=n or len(man)!=n or v.shape!=(n,320) or v.dtype!=np.float32 or not np.isfinite(v).all():raise ValueError('Frozen accounting failed: '+size)
        if sum(x['saccharide_predicted'].lower()=='true' for x in pred)!=positive:raise ValueError('Saccharide count changed')
        if [x['window_id'] for x in pred]!=[x['window_id'] for x in man]:raise ValueError('Prediction/embedding row mismatch')
        grouped={}
        for g in genes:
            if hashlib.sha256(g['sequence'].encode()).hexdigest()!=g['sequence_sha256']:raise ValueError('Input sequence hash mismatch')
            grouped.setdefault(g['window_id'],[]).append(g)
        for i,m in enumerate(man):
            ordered=sorted(grouped[m['window_id']],key=lambda g:int(g['relative_position']))
            if int(m['embedding_row'])!=i or m['ordered_sequence_sha256'].split()!=[g['sequence_sha256'] for g in ordered] or m['ordered_protein_ids'].split()!=[g['protein_id'] for g in ordered]:raise ValueError('Ordered input correspondence failed')
            if m['membership_mask'].split()!=['1']*len(ordered):raise ValueError('Unexpected membership/padding')
        result.append(dict(size=size,windows=757,CDS=len(genes),success=n,saccharide_positive=positive,excluded=sum(x['status']!='PASS' for x in qc),embedding_shape=list(v.shape),row_correspondence=True))
    return dict(status='PASS',runs=result)
def load(path,name):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m
def regress(cache,output,selection):
    import numpy as np
    cache=Path(cache).resolve();output=Path(output).resolve();verify(cache)
    if output.exists() or output.is_relative_to(cache) or output.is_relative_to(PACKAGE):raise ValueError('Use a new output directory outside inputs/software')
    wanted={'torch':'2.5.1+cpu','numpy':'1.26.4','fair-esm':'2.0.0'};versions={k:importlib.metadata.version(k) for k in wanted}
    if versions!=wanted:raise ValueError('Use documented runtime: '+str(versions))
    pairs=json.loads(Path(selection).read_text())['pairs']
    if not 1<=len(pairs)<=8:raise ValueError('Representative regression is limited to eight paired occurrences')
    adapter=load(cache/'historical/inference.py','step2_pinned_inference');wrapper=load(cache/'analysis/viral_glycan_benchmark_20261016/scripts/benchmark.py','step2_pinned_benchmark')
    prepared={};sequences={};expected={};arrays={};indexes={}
    for size in ['15','11']:
        d=cache/'runs'/size;eligible,_,_=adapter.prepare_windows(d);lookup={w['window_id']:(w,g) for w,g in eligible};prepared[size]=[]
        for pair in pairs:
            wid=pair['window'+size+'_id'];prepared[size].append((pair['reason'],*lookup[wid]))
            for g in lookup[wid][1]:sequences[g['sequence_sha256']]=g['sequence']
        expected[size]={x['window_id']:x for x in rows(d/'viral_bgc_predictions.tsv')};arrays[size]=np.load(d/'viral_bgc_embeddings.npy',allow_pickle=False);indexes[size]={x['window_id']:int(x['embedding_row']) for x in rows(d/'embedding_manifest.tsv')}
    output.mkdir(parents=True);(output/'protein_cache').mkdir();started=time.monotonic();torch=wrapper.torch_setup();vectors={}
    encoding=adapter.encode_missing(torch,cache/'models/esm2_t6_8M_UR50D.pt',sequences,vectors,output/'protein_cache',1)
    if encoding['historical_cache_hits'] or encoding['persistent_cache_hits']:raise ValueError('Unexpected historical protein-cache reuse')
    model=wrapper.load_model(torch);before={k:v.detach().clone() for k,v in model.state_dict().items()};comparisons=[]
    for size in ['15','11']:
        regenerated=[]
        for reason,w,genes in prepared[size]:
            scores,pooled,core=wrapper.infer(torch,model,[vectors[g['sequence_sha256']] for g in genes],[1]*len(genes),True);labels,threshold=wrapper.predicted(scores);old=expected[size][w['window_id']]
            oldscores=np.array([float(old['probability_'+k]) for k in adapter.LABELS]);oldvec=arrays[size][indexes[size][w['window_id']]];ds=float(np.max(np.abs(scores-oldscores)));dv=float(np.max(np.abs(pooled-oldvec)));equal=bool(np.array_equal(pooled,core));passed=ds<=1e-6 and dv<=1e-6 and equal and labels==old['predicted_labels']
            comparisons.append(dict(size=size,window_id=w['window_id'],selection_reason=reason,n_genes=len(genes),model_input_shape=f'(1,{len(genes)},320)',contextual_pool_shape='(1,320)',only_actual_positions=True,position_count_below_128=len(genes)<128,pooled_equals_classifier_input=True,membership_pool_exactly_equal=equal,score_max_abs_error=ds,embedding_max_abs_error=dv,predicted_labels_match=labels==old['predicted_labels'],threshold_error=abs(threshold-float(old['threshold'])),passed=passed));regenerated.append(pooled)
        np.save(output/('regenerated_'+size+'_embeddings.npy'),np.stack(regenerated),allow_pickle=False)
    unchanged=all(torch.equal(before[k],v) for k,v in model.state_dict().items());write(output/'regression_comparisons.tsv',comparisons)
    result=dict(status='PASS' if all(x['passed'] for x in comparisons) and unchanged else 'FAIL',representative_pairs=len(pairs),windows_inferred=len(comparisons),full_cohort_inference=False,protein_encoding=encoding,all_proteins_encoded_from_sequences=True,historical_pilot_cache_or_LMDB_used=False,model_weights_unchanged=unchanged,max_score_abs_error=max(x['score_max_abs_error'] for x in comparisons),max_embedding_abs_error=max(x['embedding_max_abs_error'] for x in comparisons),tolerance=1e-6,all_actual_position_pooling=True,positional_encoding_quirk_preserved=True,python=platform.python_version(),packages=versions,elapsed_seconds=time.monotonic()-started,model_checkpoint_sha256=digest(cache/'analysis/viral_glycan_benchmark_20261016/sources/classifier.pt'),esm_checkpoint_sha256=digest(cache/'models/esm2_t6_8M_UR50D.pt'),source_adapter_sha256=digest(cache/'historical/inference.py'),source_wrapper_sha256=digest(cache/'analysis/viral_glycan_benchmark_20261016/scripts/benchmark.py'))
    save(output/'regression_receipt.json',result)
    if result['status']!='PASS':raise ValueError('Representative regression differs; inspect new comparison evidence')
    return result
def infer_full(cache,output,size,full_cohort=False):
    if not full_cohort:raise ValueError('Full-cohort inference requires the explicit --full-cohort flag and separate execution approval')
    import numpy as np
    cache=Path(cache).resolve();output=Path(output).resolve();verify(cache)
    if output.exists() or output.is_relative_to(cache) or output.is_relative_to(PACKAGE):raise ValueError('Use a new output directory outside software/scientific inputs')
    versions={k:importlib.metadata.version(k) for k in ['torch','numpy','fair-esm']}
    if versions!={'torch':'2.5.1+cpu','numpy':'1.26.4','fair-esm':'2.0.0'}:raise ValueError('Use the documented frozen runtime')
    adapter=load(cache/'historical/inference.py','step2_full_pinned_inference');wrapper=load(cache/'analysis/viral_glycan_benchmark_20261016/scripts/benchmark.py','step2_full_pinned_benchmark')
    eligible,qc,sequences=adapter.prepare_windows(cache/'runs'/size);output.mkdir(parents=True);(output/'protein_cache').mkdir();started=time.monotonic();torch=wrapper.torch_setup();vectors={}
    encoding=adapter.encode_missing(torch,cache/'models/esm2_t6_8M_UR50D.pt',sequences,vectors,output/'protein_cache',1)
    model=wrapper.load_model(torch);before={k:v.detach().clone() for k,v in model.state_dict().items()};pred=[];manifest=[];emb=[]
    settings=json.loads((cache/'runs'/size/'inference_provenance.json').read_text())['settings'];protein_settings={k:settings[k] for k in ['esm_checkpoint_sha256','esm_layer','esm_max_residues','esm_pooling','device','dtype','versions']};protein_settings['esm_batch_size']=1
    for window,genes in eligible:
        scores,pooled,core=wrapper.infer(torch,model,[vectors[g['sequence_sha256']] for g in genes],[1]*len(genes),True)
        if not np.array_equal(pooled,core):raise ValueError('Actual-position pooling check failed')
        labels,threshold=wrapper.predicted(scores);common=dict(window_id=window['window_id'],genome_accession=window['genome_accession'],anchor_protein_id=window['anchor_protein_id'],n_genes=len(genes),strategy='CORE_ONLY')
        pred.append(dict(**common,**{'probability_'+k:float(v) for k,v in zip(adapter.LABELS,scores)},saccharide_probability=float(scores[5]),saccharide_rank=1+int((scores>scores[5]).sum()),predicted_labels=labels,threshold=threshold,saccharide_predicted='Saccharide' in labels.split(),probability_semantics='raw frozen-model sigmoid outputs; seven scores need not sum to one'))
        key=adapter.json_digest(dict(settings=settings,protein_settings=protein_settings,ordered_sequence_sha256=[g['sequence_sha256'] for g in genes]))
        manifest.append(dict(**common,embedding_row=len(emb),pooling=settings['pooling'],membership_mask=' '.join(['1']*len(genes)),ordered_protein_ids=' '.join(g['protein_id'] for g in genes),ordered_sequence_sha256=' '.join(g['sequence_sha256'] for g in genes),input_cache_key=key,cache_reused=False,classifier_sha256=adapter.CLASSIFIER_SHA,esm_sha256=adapter.ESM_SHA));emb.append(pooled)
        qc.append(dict(window_id=window['window_id'],status='PASS',reason='',n_genes=len(genes),input_shape=f'(1,{len(genes)},320)',finite=True,pooled_matches_classifier_input=True,actual_and_published_pooling_equal=True,cache_reused=False,n_truncated_proteins=sum(len(g['sequence'])>1022 for g in genes),n_sequence_cache_hits=0))
    if not all(torch.equal(before[k],v) for k,v in model.state_dict().items()):raise ValueError('Model state changed')
    np.save(output/'viral_bgc_embeddings.npy',np.stack(emb).astype(np.float32),allow_pickle=False)
    adapter.write_tsv(output/'viral_bgc_predictions.tsv',pred,adapter.PREDICTION_FIELDS);adapter.write_tsv(output/'embedding_manifest.tsv',manifest,adapter.MANIFEST_FIELDS);adapter.write_tsv(output/'inference_qc.tsv',qc,adapter.QC_FIELDS)
    result=dict(status='GENERATED_NOT_YET_COMPARED',size=size,successful_windows=len(pred),excluded=len(qc)-len(pred),encoding=encoding,cold_cache=True,model_weights_unchanged=True,historical_sources_unchanged=True,elapsed_seconds=time.monotonic()-started,limitation='Compare this new directory against frozen outputs before claiming full reproduction; portable full-cohort execution was not run for release preparation')
    save(output/'portable_execution.json',result);return result

def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    for name in ['prepare','verify','regress','infer']:
        s=sub.add_parser(name);s.add_argument('--cache',type=Path,required=True)
        if name=='prepare':s.add_argument('--step1-archive',type=Path);s.add_argument('--step2-data',type=Path);s.add_argument('--offline',action='store_true')
        if name=='infer':s.add_argument('--size',choices=['15','11'],required=True);s.add_argument('--output',type=Path,required=True);s.add_argument('--full-cohort',action='store_true')
        if name=='regress':s.add_argument('--output',type=Path,required=True);s.add_argument('--selection',type=Path,default=PACKAGE/'tests/representative_windows.json')
    a=p.parse_args()
    if a.command=='prepare':result=prepare(a.cache,a.step1_archive,a.step2_data,a.offline)
    elif a.command=='verify':result=verify(a.cache)
    elif a.command=='infer':result=infer_full(a.cache,a.output,a.size,a.full_cohort)
    else:result=regress(a.cache,a.output,a.selection)
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
