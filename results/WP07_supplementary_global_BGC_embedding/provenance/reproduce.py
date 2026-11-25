"""Verify and replay the documented baseline from publicly archived inputs."""
import argparse
from collections import Counter
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import tempfile

PACKAGE = Path(__file__).resolve().parents[1]
REPO = PACKAGE.parents[1]
FIG = PACKAGE/'07_supplementary_figure'
FIGURE_NAME = 'SuppFig_global_leiden_biosynthetic_classes_minimalist_v2'
PREFIX = 'analysis/viral_gt_clustering_20261024/'
REF = 'analysis/viral_glycan_benchmark_20261016/'
LOC = 'analysis/viral_gt_umap_20261024/'
BASE = 'analysis/viral_gt_mining_20261017/runs/baseline_594_v2/'

def read(path):
    with path.open(encoding='utf-8-sig',newline='') as f:
        return list(csv.DictReader(f,delimiter='\t' if path.suffix=='.tsv' else ','))

def sha(p):
    with p.open('rb') as stream:
        return hashlib.file_digest(stream,'sha256').hexdigest()

def git_show(commit, path):
    return subprocess.check_output(['git','-C',str(REPO),'show',f'{commit}:{path}'])

def verify_checksum_manifest():
    for r in read(PACKAGE/'provenance/ARCHIVE_SHA256.tsv'):
        f=PACKAGE/r['relative_path']
        if f.stat().st_size!=int(r['size_bytes']) or sha(f)!=r['sha256']:
            raise ValueError(f'Checksum mismatch for {f}')
    print('Public package file checksums verified')

def original_files(temp):
    """Materialize a verified public-only workspace for scientific inputs."""
    files={}
    for r in read(PACKAGE/'provenance/INPUT_MANIFEST.tsv'):
        rel=r['original_workspace_path']
        path=Path(r['repository_path'])
        assert not path.is_absolute() and '..' not in path.parts
        if r['public_commit']=='RELEASE_COMMIT':
            source=REPO/path
        else:
            source=temp/'pinned'/path
            source.parent.mkdir(parents=True,exist_ok=True)
            source.write_bytes(git_show(r['public_commit'],path.as_posix()))
        assert sha(source)==r['sha256'] and source.stat().st_size==int(r['bytes'])
        files[rel]=source
    return files

def graph_baseline(replay):
    # Use the original distance and union-kNN arithmetic. BLAS threads are bounded before NumPy import.
    os.environ['OPENBLAS_NUM_THREADS']='1'
    os.environ['OMP_NUM_THREADS']='1'
    os.environ['MKL_NUM_THREADS']='1'
    import numpy as np
    with tempfile.TemporaryDirectory(prefix='bgc-public-inputs-') as scratch:
        temp=Path(scratch)
        files=original_files(temp)
        ref=np.load(files[REF+'mibig_CORE_ONLY_embeddings.npy'],allow_pickle=False)
        rep=np.load(files[LOC+'tables/representative_15gene_embeddings.npy'],allow_pickle=False)
        raw=np.load(files[BASE+'viral_bgc_embeddings.npy'],allow_pickle=False)
        rmeta=read(files[REF+'mibig_matched_reference_metadata.tsv'])
        vmeta=read(files[LOC+'tables/included_locus_representatives.csv'])
        selected=read(files[PREFIX+'tables/selected_287_input_order.csv'])
        row_map={r['window_id']:int(r['embedding_row']) for r in read(files[BASE+'embedding_manifest.tsv'])}
        assert ref.shape==(2502,320) and rep.shape==(556,320) and raw.shape==(737,320)
        assert np.array_equal(rep,np.stack([raw[row_map[r['representative_window_id']]] for r in vmeta]))
        idx=[int(r['original_556_row']) for r in selected]
        assert len(set(idx))==287
        assert [int(r['selected_row']) for r in selected]==list(range(287))
        assert [r['viral_locus_id'] for r in selected]==[vmeta[i]['viral_locus_id'] for i in idx]
        assert [r['representative_window_id'] for r in selected]==[vmeta[i]['representative_window_id'] for i in idx]
        assert [int(r['embedding_row']) for r in rmeta]==list(range(2502))
        ids=np.array([r['BGC'] for r in rmeta]+[r['viral_locus_id'] for r in selected])
        matrix=np.concatenate((ref,rep[idx]))
        assert matrix.shape==(2789,320) and len(set(ids))==2789
        unit=matrix.astype(np.float64)
        unit/=np.linalg.norm(unit,axis=1,keepdims=True)
        cosine=np.clip(unit@unit.T,-1,1)
        cosine=(cosine+cosine.T)/2
        np.fill_diagonal(cosine,1)
        _,inv=np.unique(matrix,axis=0,return_inverse=True)
        inv=inv.reshape(-1)
        cosine[inv[:,None]==inv[None,:]]=1
        distance=np.clip(1-cosine,0,2)
        np.fill_diagonal(distance,0)
        edges=set()
        for i in range(len(ids)):
            ranked=np.lexsort((ids,distance[i]))
            edges.update((min(i,j),max(i,j)) for j in ranked[ranked!=i][:20] if cosine[i,j]>0)
        observed={(str(ids[i]),str(ids[j])):cosine[i,j] for i,j in edges}
        archived={(r['source'],r['target']):float(r['weight']) for r in read(files[PREFIX+'graphs/global_joint_k20_edges.csv'])}
        if set(archived)!=set(observed):
            raise ValueError(f'Graph edge identities differ: {len(set(archived)^set(observed))} edges')
        max_error=float(max(abs(v-observed[k]) for k,v in archived.items()))
        assert max_error<1e-12
        grid=read(files[PREFIX+'assignments/global_joint_leiden_grid.csv'])
        saved={r['sample_id']:r['k20_r1_s42'] for r in grid}
        assert set(saved)==set(ids) and len(set(saved.values()))==10
        degrees=Counter()
        internal=0.0
        total=sum(archived.values())
        for (left,right),weight in archived.items():
            degrees[saved[left]]+=weight
            degrees[saved[right]]+=weight
            if saved[left]==saved[right]:
                internal+=weight
        rb_quality=2*internal-sum(value**2 for value in degrees.values())/(2*total)
        baseline=next(r for r in json.loads(files[PREFIX+'provenance/global_joint_optimizer_iterations.json'].read_text()) if r['run_id']=='k20_r1_s42')
        assert abs(rb_quality-baseline['rb_objective'])<1e-7
        assert abs(rb_quality/(2*total)-baseline['weighted_modularity_gamma1'])<1e-12
        counts=Counter()
        for r in rmeta:
            counts.update((saved[r['BGC']],c) for c in r['label'].split())
        table=read(files[PREFIX+'tables/joint_reference_class_composition.csv'])
        global_rows=[r for r in table if r['analysis']=='global_joint']
        assert len(global_rows)==70
        assert all(counts[r['joint_community'],r['reference_class']]==int(r['n_reference_with_label']) for r in global_rows)
        coordinates=read(files[LOC+'joint_umap_samples.csv'])
        coord={r['sample_id']:r for r in coordinates}
        assert len(coord)==3058
        display=read(FIG/'source_data/FIGURE_SOURCE_DATA.tsv')
        assert len(display)==2789 and {r['sample_id'] for r in display}==set(ids)
        for r in display:
            sid=r['sample_id']
            assert r['global_leiden_community']==saved[sid]
            assert (r['UMAP1'],r['UMAP2'],r['reference_class_original'])==(coord[sid]['UMAP1'],coord[sid]['UMAP2'],coord[sid]['reference_class'])
        result=dict(reference_shape=list(ref.shape),representative_shape=list(rep.shape),global_input_shape=list(matrix.shape),
                    weighted_edges=len(edges),max_cosine_weight_error=max_error,archived_communities=10,class_counts_verified=70,
                    coordinates_verified=len(display),weighted_modularity=baseline['weighted_modularity_gamma1'],
                    rb_quality=baseline['rb_objective'],Leiden_optimizer_replayed=False)
        if replay:
            # Optional reader action: run the saved baseline only, using its historical engine and config.
            # This codepath is not part of the public release validation.
            engine=files[PREFIX+'scripts/engine.py']
            cfg=json.loads(files[PREFIX+'config.json'].read_text(encoding='utf-8-sig'))
            runtime=temp/'engine'
            (runtime/'scripts').mkdir(parents=True)
            for name in ['common.py','engine.py']:
                shutil.copyfile(files[PREFIX+'scripts/'+name],runtime/'scripts'/name)
            cfg.update(workspace_root=str(temp),graph_runtime='',analysis_dependencies='')
            (runtime/'config.json').write_text(json.dumps(cfg),encoding='utf-8')
            sys.path.insert(0,str(runtime/'scripts'))
            import engine
            # Replay with the exact archived floating-point weights after confirming the graph from embeddings.
            index={sid:i for i,sid in enumerate(ids)}
            pairs=sorted((index[left],index[right]) for left,right in archived)
            g=engine.ig.Graph(n=len(ids),edges=pairs,directed=False)
            g.es['weight']=[archived[ids[i],ids[j]] for i,j in pairs]
            labels,_=engine.leiden(g,ids,cfg['baseline']['resolution'],cfg['baseline']['seed'],'GC')
            assert all(saved[sid]==label for sid,label in zip(ids,labels))
            result['Leiden_optimizer_replayed']=True
        print(json.dumps(result,indent=2))

def plot(font,output):
    if not font.is_file():
        raise FileNotFoundError(font)
    output=output.resolve()
    output.mkdir(parents=True,exist_ok=False)
    runtime=output/'replay'/'figure'
    runtime.mkdir(parents=True)
    for name in ['plot_integrated_umap.py','config.json']:
        src=FIG/'scripts'/name
        shutil.copyfile(src,runtime/name)
    for name in ['FIGURE_SOURCE_DATA.tsv','COMMUNITY_ANNOTATION_KEY.tsv','REFERENCE_CLASS_MARGINAL_COUNTS.tsv','COLOR_PALETTE.tsv']:
        shutil.copyfile(FIG/'source_data'/name,runtime/name)
    config=json.loads((runtime/'config.json').read_text())
    config['font_path']=str(font.resolve())
    (runtime/'config.json').write_text(json.dumps(config),encoding='utf-8')
    subprocess.run([sys.executable,'-B',str(runtime/'plot_integrated_umap.py')],cwd=runtime,check=True)
    from PIL import Image, ImageChops
    original=Image.open(FIG/'figures'/f'{FIGURE_NAME}.png').convert('RGB')
    exported=Image.open(runtime/'figures'/f'{FIGURE_NAME}.png').convert('RGB')
    assert original.size==exported.size
    delta=ImageChops.difference(original,exported)
    print(json.dumps(dict(plot_resolution=list(exported.size),pixel_identical=delta.getbbox() is None,
                          font_sha256=sha(font),output_directory=str(runtime/'figures')),indent=2))

def restore(output):
    """Materialize historical workspace paths for readers rebuilding upstream stages."""
    output=output.resolve()
    output.mkdir(parents=True,exist_ok=False)
    with tempfile.TemporaryDirectory(prefix='bgc-public-restore-') as scratch:
        for relative,source in original_files(Path(scratch)).items():
            rel=Path(relative)
            assert not rel.is_absolute() and '..' not in rel.parts,relative
            target=output/rel
            target.parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(source,target)
    bundle_rows=read(PACKAGE/'provenance/BUNDLED_SCIENTIFIC_INPUTS.tsv')
    for bundle in {r['archive_path'] for r in bundle_rows}:
        expected={r['member']:r for r in bundle_rows if r['archive_path']==bundle}
        with tarfile.open(PACKAGE/bundle) as ar:
            for member in ar:
                original=expected.pop(member.name)
                rel=Path(member.name)
                assert member.isfile() and not rel.is_absolute() and '..' not in rel.parts
                target=output/rel
                assert not target.exists()
                target.parent.mkdir(parents=True,exist_ok=True)
                with ar.extractfile(member) as stream, target.open('xb') as f:
                    shutil.copyfileobj(stream,f)
                assert sha(target)==original['sha256']
        assert not expected
    print('Restored the workspace and bundled inputs. See EXTERNAL_DEPENDENCIES.tsv for publisher-hosted inputs.')

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parts=parser.add_subparsers(dest='step',required=True)
    parts.add_parser('checksums')
    graph=parts.add_parser('graph')
    graph.add_argument('--replay-leiden',action='store_true')
    restore_cmd=parts.add_parser('restore')
    restore_cmd.add_argument('--output',type=Path,required=True)
    draw=parts.add_parser('figure')
    draw.add_argument('--font',type=Path,required=True)
    draw.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.step=='checksums': verify_checksum_manifest()
    elif args.step=='graph': graph_baseline(args.replay_leiden)
    elif args.step=='restore': restore(args.output)
    else: plot(args.font,args.output)

if __name__=='__main__':
    main()
