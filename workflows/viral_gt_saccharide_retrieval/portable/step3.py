#!/usr/bin/env python3
"""Step 3 retrieval-only adapter. Original historical mains are never executed."""
from pathlib import Path
import argparse, ast, collections, csv, hashlib, importlib.metadata
import json, sys, urllib.request, zipfile

PACKAGE = Path(__file__).resolve().parents[1]
ATOL = RTOL = 1e-12
PRIMARY = 'public/analysis/viral_gt_mining_20261017/runs/'
REFERENCE = 'public/analysis/viral_glycan_benchmark_20261016/'
RUNS = {'15': 'baseline_594_v2', '11': 'window11_594_v1'}

def sha(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()

def read(path):
    with Path(path).open(encoding='utf-8-sig', newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))

def write(path, rows, fields=None):
    rows = list(rows)
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    fields = fields or list(rows[0])
    with path.open('w', encoding='utf8', newline='') as handle:
        writer = csv.DictWriter(handle, fields, delimiter='\t', lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)

def save(path, value):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n', encoding='utf8')

def require(condition, message):
    if not condition:
        raise ValueError(message)

def spec():
    return json.loads((PACKAGE/'resources.json').read_text())

def request(url):
    with urllib.request.urlopen(urllib.request.Request(url, headers={'User-Agent':'Viral-GT-Step3/1.0'}), timeout=90) as response:
        return response.read()

def materialize(cache, archive_arguments):
    """Fetch only the pinned public Step 1/2 resources; never download model weights."""
    cache = Path(cache).resolve(); cache.mkdir(parents=True, exist_ok=True)
    supplied = dict(item.split('=', 1) for item in archive_arguments)
    manifest = spec(); opened = {}
    try:
        for row in manifest['resources']:
            dest = cache/row['key']; dest.parent.mkdir(parents=True, exist_ok=True)
            if dest.exists():
                require(sha(dest) == row['sha256'], 'Existing cache hash mismatch: '+row['key'])
                continue
            if row['provider'] == 'GitHub':
                data = request(row['url'])
            else:
                aid = row['archive_id']; archive = manifest['archives'][aid]
                if aid not in opened:
                    archive_path = Path(supplied[aid]).resolve() if aid in supplied else cache/'archives'/archive['name']
                    if not archive_path.exists():
                        archive_path.parent.mkdir(parents=True, exist_ok=True)
                        # A failed/partial download is retained for inspection, not silently replaced.
                        with urllib.request.urlopen(archive['url'], timeout=90) as response, archive_path.open('xb') as handle:
                            while True:
                                chunk = response.read(1024*1024)
                                if not chunk: break
                                handle.write(chunk)
                    require(sha(archive_path) == archive['sha256'], 'Public archive hash mismatch')
                    opened[aid] = zipfile.ZipFile(archive_path)
                data = opened[aid].read(row['public_path'])
            require(hashlib.sha256(data).hexdigest() == row['sha256'], 'Downloaded/member hash mismatch')
            with dest.open('xb') as handle: handle.write(data)
    finally:
        for archive in opened.values(): archive.close()
    return verify(cache, None)

def verify(cache, dataset):
    manifest = spec(); checked = []
    for row in manifest['resources']:
        path = Path(cache)/row['key']
        require(path.is_file() and path.stat().st_size == row['bytes'] and sha(path) == row['sha256'],
                'Missing or changed public input: '+row['key'])
        checked.append(row['key'])
    for row in manifest['historical_sources']:
        require(sha(PACKAGE/row['file']) == row['sha256'], 'Historical source changed: '+row['file'])
    if dataset is not None:
        for row in manifest['frozen_dataset_files']:
            require(sha(Path(dataset)/row['file']) == row['sha256'], 'Frozen output changed: '+row['file'])
        # Final packages also carry full metadata/evidence integrity manifests.
        sums = Path(dataset)/'SHA256SUMS.tsv'
        if sums.exists():
            for row in read(sums):
                require(sha(Path(dataset)/row['path']) == row['sha256'], 'Dataset integrity failure: '+row['path'])
    return {'passed': True, 'public_resources_verified': len(checked), 'historical_scripts_verified': 3,
            'frozen_outputs_verified': 4 if dataset is not None else 0, 'model_resources_required': 0}

def paired_function(reference_dir):
    """Load only the pinned historical reference_metrics function through AST extraction."""
    import numpy as np
    path = PACKAGE/'historical/window_size/paired.py'
    expected = next(row['sha256'] for row in spec()['historical_sources'] if row['file'].endswith('/paired.py'))
    require(sha(path) == expected, 'Historical paired.py identity changed')
    tree = ast.parse(path.read_text(), filename=str(path))
    nodes = [node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == 'reference_metrics']
    require(len(nodes) == 1, 'Expected exactly one historical reference_metrics function')
    namespace = {'np': np, 'read': read, 'REF': reference_dir}
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(path), 'exec'), namespace)
    return namespace['reference_metrics']

def rank_one(reference_unit, reference_ids, query):
    """Independent per-query expression matching published inference.py:localize."""
    import numpy as np
    query = np.asarray(query, dtype=np.float64)
    norm = np.linalg.norm(query)
    require(np.isfinite(query).all() and norm > 0, 'Invalid query vector')
    similarity = reference_unit @ (query/norm)
    order = np.lexsort((reference_ids, -similarity))
    return order, similarity

def status(ids, inferred, selected):
    if set(ids) & selected: return 'candidate'
    return 'not_selected' if set(ids) <= inferred else 'unresolved'

def support(a, b):
    return 'Both' if a == b == 'candidate' else '11-only' if b == 'candidate' else '15-only' if a == 'candidate' else ''

def compare_counts(cache, dataset, out, evidence):
    """Reproduce the original fixed-15-locus accounting without changing candidate records."""
    base = Path(cache)/(PRIMARY+RUNS['15'])
    mapping = read(Path(cache)/(PRIMARY+RUNS['11'])/'window15_window11_mapping.tsv')
    old_to_new = {row['window15_id']: row['window11_id'] for row in mapping}
    require(len(old_to_new) == len(set(old_to_new.values())) == 757, 'Occurrence crosswalk changed')
    windows = read(base/'viral_gt_windows.tsv')
    loci = read(base/'candidate_loci.tsv'); genomes = read(base/'genome_units.tsv')
    entities = {'window': {row['window_id']:[row['window_id']] for row in windows},
                'locus': {row['locus_id']:row['window_ids'].split(';') for row in loci},
                'genome': {row['genome_unit_id']:row['window_ids'].split(';') for row in genomes}}
    require([len(entities[key]) for key in ['window','locus','genome']] == [757,565,466], 'Grouping population mismatch')
    states = {}; rows = []
    for condition in ['15','11']:
        data = evidence[condition]; inferred = set(data)
        a = {w for w,r in data.items() if r['classifier_positive']}
        b = {w for w,r in data.items() if not r['classifier_positive'] and r['p10'] >= .5}
        for name, selected in [('A',a),('B',b),('C',a|b)]:
            for level, units in entities.items():
                mapped = {key: ids if condition == '15' else [old_to_new[w] for w in ids] for key,ids in units.items()}
                classified = {key:status(ids,inferred,selected) for key,ids in mapped.items()}
                states[condition,name,level] = classified
                counts = collections.Counter(classified.values())
                rows.append(dict(condition=condition,region_basis='frozen_15gene_loci',set=name,level=level,total=len(units),
                    any_inferred=sum(bool(set(ids)&inferred) for ids in mapped.values()),
                    fully_inferred=sum(set(ids)<=inferred for ids in mapped.values()),
                    candidates=counts['candidate'],not_selected=counts['not_selected'],unresolved=counts['unresolved']))
    frozen = read(Path(dataset)/'outputs/window_size/candidate_counts_by_condition.tsv')
    lookup = {(r['condition'],r['set'],r['level']):r for r in rows}
    require(len(frozen) == len(rows) == 18, 'Unexpected frozen count table')
    for row in frozen:
        reproduced = lookup[row['condition'],row['set'],row['level']]
        require(all(str(reproduced[key]) == value for key,value in row.items()), 'Historical count table mismatch')
    write(out/'candidate_counts_by_condition_reproduced.tsv',rows)

    snapshot = read(Path(dataset)/'retrospective_sensitivity/historical_candidate_membership.tsv')
    frozen_candidates = {r['genome_unit_id']:r for r in snapshot}
    exclusions = read(Path(dataset)/'retrospective_sensitivity/historical_figure_exclusions.tsv')
    excluded = {r['genome_unit_id'] for r in exclusions}
    require(len(excluded) == 3 and set(';'.join(r['source_genome_accessions'] for r in exclusions).split(';')) ==
            {'KC511810.1','KX960781.1','KY189087.1'}, 'Historical figure exclusions changed')
    sensitivity = []; boundaries = []
    window_meta = {}
    for condition in ['15','11']:
        window_meta[condition] = {r['window_id']:r for r in read(Path(cache)/(PRIMARY+RUNS[condition])/'viral_gt_windows.tsv')}
    for genome in genomes:
        gid = genome['genome_unit_id']; ids15 = genome['window_ids'].split(';')
        ids = {'15':ids15,'11':[old_to_new[w] for w in ids15]}
        observed = {}; strict = {}
        for condition in ['15','11']:
            data = evidence[condition]; inferred = set(data)
            observed[condition] = states[condition,'C','genome'][gid]
            selected = {w for w,r in data.items() if r['classifier_positive'] or r['p10'] > .5}
            strict[condition] = status(ids[condition],inferred,selected)
            for wid in ids[condition]:
                if wid in data and data[wid]['p10'] == .5:
                    wm = window_meta[condition][wid]; item = data[wid]
                    boundaries.append(dict(condition=condition,window_id=wid,genome_unit_id=gid,
                        source_genome_accessions=genome['source_genome_accessions'],anchor_protein_id=wm['anchor_protein_id'],
                        p10=.5,classifier_positive=item['classifier_positive'],
                        retained_as_window_under_strict=item['classifier_positive']))
        historical_selected = 'candidate' in observed.values(); strict_selected = 'candidate' in strict.values()
        historical_support = support(observed['15'],observed['11']); strict_support = support(strict['15'],strict['11'])
        if historical_selected:
            require(gid in frozen_candidates and frozen_candidates[gid]['historical_window_support'] == historical_support,
                    'Frozen 245-candidate snapshot mismatch')
            require((gid not in excluded) == (frozen_candidates[gid]['included_in_frozen_242_sankey'].lower() == 'true'),
                    'Frozen 242 membership mismatch')
        else: require(gid not in frozen_candidates, 'Unexpected historical candidate')
        sensitivity.append(dict(genome_unit_id=gid,source_genome_accessions=genome['source_genome_accessions'],
            historical_selected=historical_selected,strict_selected=strict_selected,historical_window_support=historical_support,
            strict_window_support=strict_support,included_in_frozen_sankey=historical_selected and gid not in excluded,
            included_under_strict_sensitivity=strict_selected and gid not in excluded,
            any_classifier_positive=any(evidence[c].get(w,{}).get('classifier_positive',False) for c in ['15','11'] for w in ids[c]),
            interpretation='Retrospective comparator sensitivity only; historical cohort unchanged'))
    lost = [r for r in sensitivity if r['historical_selected'] and not r['strict_selected']]
    changes = [r for r in sensitivity if r['strict_selected'] and r['historical_window_support'] != r['strict_window_support']]
    result = dict(historical_comparator='>= 0.5',sensitivity_comparator='> 0.5',historical_master=245,
        historical_sankey=242,strict_master=sum(r['strict_selected'] for r in sensitivity),
        strict_sankey=sum(r['included_under_strict_sensitivity'] for r in sensitivity),
        lost_genomes=len(lost),retained_genomes_with_support_change=len(changes),
        boundary_windows={c:sum(r['condition']==c for r in boundaries) for c in ['15','11']},
        classifier_positive_boundary_windows={c:sum(r['condition']==c and r['classifier_positive'] for r in boundaries) for c in ['15','11']},
        historical_support=dict(collections.Counter(r['historical_window_support'] for r in sensitivity if r['historical_selected'])),
        strict_support=dict(collections.Counter(r['strict_window_support'] for r in sensitivity if r['strict_selected'])),
        historical_cohort_modified=False,new_candidate_selection_performed=False)
    require((len(snapshot),result['strict_master'],result['strict_sankey'],len(lost),len(changes)) == (245,230,227,15,14),
            'Comparator sensitivity mismatch')
    require(result['boundary_windows'] == {'15':34,'11':26} and result['classifier_positive_boundary_windows'] == {'15':2,'11':4},
            'Boundary accounting mismatch')
    folder = out/'retrospective_sensitivity'
    write(folder/'genome_comparator_sensitivity.tsv',sensitivity)
    write(folder/'boundary_windows.tsv',boundaries)
    write(folder/'lost_genomes_under_strict_comparator.tsv',lost)
    write(folder/'retained_genomes_with_support_changes.tsv',changes)
    save(folder/'summary.json',result)
    return result

def reproduce(cache, dataset, output):
    import numpy as np
    cache = Path(cache).resolve(); dataset = Path(dataset).resolve(); output = Path(output).resolve()
    require(not output.exists(), 'Use a completely new output directory')
    for protected in [PACKAGE,cache,dataset]:
        require(output != protected and protected not in output.parents, 'Output must be separate from packages and inputs')
    verification = verify(cache,dataset)
    output.mkdir(parents=True)
    reference_dir = cache/REFERENCE
    reference = np.load(reference_dir/'mibig_CORE_ONLY_embeddings.npy',allow_pickle=False)
    metadata = read(reference_dir/'mibig_matched_reference_metadata.tsv')
    require(reference.shape == (2502,320) and np.isfinite(reference).all(), 'Reference matrix mismatch')
    require([int(r['embedding_row']) for r in metadata] == list(range(2502)), 'Reference row order mismatch')
    ids = np.asarray([r['BGC'] for r in metadata])
    positive = np.asarray(['Saccharide' in r['label'].split() for r in metadata],dtype=bool)
    require(len(set(ids)) == 2502 and int(positive.sum()) == 195 and all(r['label'].strip() for r in metadata),
            'Reference membership/labels mismatch')
    reference64 = reference.astype(np.float64); norms = np.linalg.norm(reference64,axis=1)
    require((norms>0).all(), 'Zero reference vector')
    unit = reference64/norms[:,None]
    duplicates = collections.defaultdict(list)
    for row,vector in zip(metadata,reference):
        duplicates[hashlib.sha256(vector.tobytes()).hexdigest()].append(row['BGC'])
    duplicate_groups = [v for v in duplicates.values() if len(v)>1]
    reference_report = dict(rows=2502,dimensions=320,saccharide_accessions=195,
        ordered_accession_sha256=hashlib.sha256(('\n'.join(ids.tolist())+'\n').encode()).hexdigest(),
        unique_exact_vectors=len(duplicates),duplicate_vector_groups=duplicate_groups,
        duplicates_retained_as_separate_accessions=True,missing_labels=0,
        reference_matrix_sha256=sha(reference_dir/'mibig_CORE_ONLY_embeddings.npy'),
        metadata_sha256=sha(reference_dir/'mibig_matched_reference_metadata.tsv'))
    calculate = paired_function(reference_dir)
    expected_neighbors = read(dataset/'outputs/window_size/reference_neighbors_top20.tsv')
    neighbor_groups = collections.defaultdict(list)
    for row in expected_neighbors: neighbor_groups[row['condition'],row['window_id']].append(row)
    require(len(expected_neighbors) == 29540, 'Frozen neighbor row count mismatch')
    comparisons = []; metrics = []; evidence = {}; summaries = {}
    for condition in ['15','11']:
        folder = cache/(PRIMARY+RUNS[condition])
        vectors = np.load(folder/'viral_bgc_embeddings.npy',allow_pickle=False)
        manifest = read(folder/'embedding_manifest.tsv'); predictions = read(folder/'viral_bgc_predictions.tsv')
        frozen_local = read(folder/'saccharide_localization.tsv') if condition == '15' else read(dataset/'outputs/window11_594_v1/saccharide_localization.tsv')
        n = 737 if condition == '15' else 740
        require(vectors.shape == (n,320) and len(manifest) == len(predictions) == len(frozen_local) == n, 'Query population mismatch')
        query_ids = [r['window_id'] for r in manifest]
        require(len(set(query_ids)) == n and not (set(query_ids)&set(ids.tolist())), 'Unexpected query identity collision')
        require(query_ids == [r['window_id'] for r in predictions] == [r['window_id'] for r in frozen_local] and
                [int(r['embedding_row']) for r in manifest] == list(range(n)), 'Window/embedding correspondence mismatch')
        mapped = {wid:vectors[i] for i,wid in enumerate(query_ids)}
        original_metrics, original_neighbors = calculate(mapped,None,condition)
        by_original = collections.defaultdict(list)
        for row in original_neighbors: by_original[row['window_id']].append(row)
        evidence[condition] = {}; reproduced_rows = []; max_neighbor_error = max_local_error = 0.0
        for index,(wid,pred,prior) in enumerate(zip(query_ids,predictions,frozen_local)):
            order, similarity = rank_one(unit,ids,vectors[index])
            saved = sorted(neighbor_groups[condition,wid],key=lambda r:int(r['rank']))
            computed = by_original[wid]
            require(len(saved) == len(computed) == 20, 'Top-20 row mismatch')
            per_error = 0.0
            for rank,(row,historical) in enumerate(zip(saved,computed),1):
                j = order[rank-1]
                require(row['window_id'] == historical['window_id'] == wid and int(row['rank']) == historical['rank'] == rank,
                        'Neighbor rank/window mismatch')
                require(row['reference_BGC'] == historical['reference_BGC'] == str(ids[j]) and
                        row['reference_label'] == historical['reference_label'] == metadata[j]['label'] and
                        (row['saccharide'].lower()=='true') == historical['saccharide'] == bool(positive[j]), 'Neighbor identity/label mismatch')
                for number in [historical['cosine'],float(similarity[j])]:
                    error=abs(number-float(row['cosine']));per_error=max(per_error,error)
                    require(np.isclose(number,float(row['cosine']),atol=ATOL,rtol=RTOL), 'Cosine tolerance failure')
            values = {k:float(positive[order[:k]].mean()) for k in [5,10,20]}
            for k,value in values.items():
                require(value == original_metrics[wid]['p'+str(k)] == float(prior['saccharide_fraction_top'+str(k)]), 'P@k exact mismatch')
            require(prior['nearest_reference'] == str(ids[order[0]]) and prior['nearest_reference_label'] == metadata[order[0]]['label'],
                    'Nearest reference mismatch')
            nearest_positive = next(int(j) for j in order if positive[j])
            require(prior['nearest_saccharide_reference'] == str(ids[nearest_positive]), 'Nearest Saccharide identity mismatch')
            for name,value in [('nearest_reference_cosine',float(similarity[order[0]])),
                               ('nearest_saccharide_cosine',float(similarity[nearest_positive])),
                               ('descriptive_saccharide_enrichment_top10',values[10]/positive.mean()),
                               ('reference_saccharide_fraction',float(positive.mean()))]:
                err=abs(value-float(prior[name]));max_local_error=max(max_local_error,err)
                require(np.isclose(value,float(prior[name]),atol=ATOL,rtol=RTOL), 'Localization scalar tolerance failure')
            classifier = pred['saccharide_predicted'].lower() == 'true'
            require((prior['predicted_saccharide'].lower()=='true') == classifier and
                    float(prior['saccharide_probability']) == float(pred['saccharide_probability']), 'Classifier comparison data changed')
            evidence[condition][wid] = dict(classifier_positive=classifier,p5=values[5],p10=values[10],p20=values[20])
            max_neighbor_error=max(max_neighbor_error,per_error)
            comparisons.append(dict(condition=condition,window_id=wid,embedding_row=index,
                exact_neighbor_id_rank_label_match=True,exact_P5_P10_P20_match=True,max_neighbor_cosine_abs_error=per_error,status='PASS'))
            reproduced_rows.append(dict(window_id=wid,embedding_row=index,nearest_reference=str(ids[order[0]]),
                nearest_reference_cosine=float(similarity[order[0]]),p5=values[5],p10=values[10],p20=values[20],
                classifier_positive=classifier,retrieval_ge_0p5=values[10]>=.5,union_positive=classifier or values[10]>=.5))
        require(len(neighbor_groups) >= n, 'Missing frozen query neighbor groups')
        counts = collections.Counter((r['classifier_positive'],r['p10']>=.5) for r in evidence[condition].values())
        actual = dict(evaluated=n,retrieval_positive=sum(r['p10']>=.5 for r in evidence[condition].values()),
            classifier_positive=sum(r['classifier_positive'] for r in evidence[condition].values()),
            classifier_only=counts[True,False],retrieval_only=counts[False,True],both=counts[True,True],neither=counts[False,False],
            union_positive=n-counts[False,False])
        expected = {'15':(325,246,6,85,240,406,331),'11':(337,265,7,79,258,396,344)}[condition]
        require(tuple(actual[k] for k in ['retrieval_positive','classifier_positive','classifier_only','retrieval_only','both','neither','union_positive']) == expected,
                'Expected historical screening counts not reproduced')
        summaries[condition] = dict(**actual,exact_Pk_comparisons=3*n,exact_top20_neighbor_rows=20*n,
            max_neighbor_cosine_abs_error=max_neighbor_error,max_localization_scalar_abs_error=max_local_error)
        write(output/('retrieval_'+condition+'.tsv'),reproduced_rows)
        write(output/('reference_neighbors_top20_'+condition+'.tsv'),original_neighbors)
    require(len(neighbor_groups) == 1477, 'Unexpected frozen neighbor query groups')
    write(output/'window_comparisons.tsv',comparisons)
    save(output/'reference_membership_validation.json',reference_report)
    sensitivity = compare_counts(cache,dataset,output,evidence)
    forbidden_modules = [name for name in sys.modules if name.split('.')[0] in {'torch','esm','lmdb','transformers'}]
    require(not forbidden_modules, 'Inference library unexpectedly imported')
    result = dict(passed=True,scope='Full frozen-query retrieval reproduction only; no model inference or embedding generation',
        evaluated_windows=1477,reference=reference_report,conditions=summaries,neighbor_rows_checked=29540,
        exact_Pk_values_checked=4431,neighbor_id_rank_label_comparisons='exact',
        cosine_comparison=dict(atol=ATOL,rtol=RTOL,byte_identity_claim=False),
        historical_count_table_rows_exact=18,threshold_sensitivity=sensitivity,verification=verification,
        original_historical_mains_executed=False,only_historical_function_executed='paired.py:reference_metrics',
        independent_numeric_check='inference.py:localize-style per-query normalized dot product',
        model_inference_calls=0,embedding_generation_calls=0,inference_modules_imported=forbidden_modules,
        environment=dict(python=sys.version,numpy=np.__version__,platform=sys.platform))
    save(output/'reproduction_validation.json',result)
    return result

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    subs=parser.add_subparsers(dest='command',required=True)
    fetch=subs.add_parser('fetch');fetch.add_argument('--cache',type=Path,required=True);fetch.add_argument('--archive',action='append',default=[],help='step1=/path/to/public.zip or step2=/path/to/public.zip')
    check=subs.add_parser('verify');check.add_argument('--cache',type=Path,required=True);check.add_argument('--dataset',type=Path,required=True)
    run=subs.add_parser('reproduce');run.add_argument('--cache',type=Path,required=True);run.add_argument('--dataset',type=Path,required=True);run.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.command=='fetch':result=materialize(args.cache,args.archive)
    elif args.command=='verify':result=verify(args.cache,args.dataset)
    else:result=reproduce(args.cache,args.dataset,args.output)
    print(json.dumps(result,indent=2,allow_nan=False))
if __name__=='__main__':
    main()
