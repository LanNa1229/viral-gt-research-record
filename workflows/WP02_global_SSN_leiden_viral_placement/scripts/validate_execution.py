#!/usr/bin/env python3
"""Validate completed preserved/fresh workflow runs without changing research data."""
import argparse
import ast
from collections import Counter, defaultdict
import hashlib
import importlib.metadata as metadata
import importlib.util
import json
import math
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd
from packaging.requirements import Requirement


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read(path):
    return pd.read_csv(path, sep='\t' if path.suffix == '.tsv' else ',', dtype=str, keep_default_na=False)


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def semantic_differences(a, b, key):
    if set(a.columns) != set(b.columns) or set(a[key]) != set(b[key]) or a[key].duplicated().any() or b[key].duplicated().any():
        return {'schema_or_ID_difference': True}
    a = a.set_index(key).sort_index()
    b = b.set_index(key).sort_index()
    differences = {}
    for column in a.columns:
        av, bv = a[column], b[column]
        same = av.eq(bv).to_numpy()
        an, bn = pd.to_numeric(av, errors='coerce'), pd.to_numeric(bv, errors='coerce')
        same |= np.isclose(an, bn, rtol=1e-10, atol=1e-12, equal_nan=False)
        if not same.all():
            differences[column] = int((~same).sum())
    return differences


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root', type=Path, required=True)
    parser.add_argument('--preserved-run', type=Path, required=True)
    parser.add_argument('--fresh-run', type=Path, required=True)
    parser.add_argument('--report', type=Path, required=True)
    args = parser.parse_args()
    workflow = Path(__file__).resolve().parents[1]
    checks = []

    def check(name, condition, observed=None):
        checks.append(dict(test=name, status='PASS' if bool(condition) else 'FAIL', observed=observed))

    analyzer = load(workflow/'original/analyze_ssn_leiden.py', 'ssn_analyzer')
    plotter = load(workflow/'original/make_figure2b_alluvial.py', 'ssn_figure2b')
    for path in sorted((workflow/'original').glob('*.py')) + sorted((workflow/'scripts').glob('*.py')):
        ast.parse(path.read_text())
    check('Python_syntax', True)
    original_manifest=read(workflow/'manifest.tsv').set_index('path')
    for name in ['analyze_ssn_leiden.py','classify_viral_structural_placement_v3.py']:
        check(name+'_unchanged',sha(workflow/'original'/name)==original_manifest.loc['original/'+name,'sha256'])
    for name in ['analyze_ssn_leiden.py', 'classify_viral_structural_placement_v3.py']:
        result = subprocess.run([sys.executable, '-B', str(workflow/'original'/name), '--help'], capture_output=True, text=True)
        check(name+'_CLI', result.returncode == 0 and '--nodes' in result.stdout)
    dependencies = {}
    for line in (workflow/'environment/requirements.lock.txt').read_text().splitlines():
        if not line.strip() or line.startswith('#'):
            continue
        req = Requirement(line)
        version = metadata.version(req.name)
        dependencies[req.name] = version
        check('dependency_'+req.name, version in req.specifier, version)

    # Independent expected values exercise the important graph-construction boundaries.
    fixture = pd.DataFrame([
        ['A','B',.9,.8,.85,.8,.8,1e-8,100],
        ['B','A',.4,.9,.7,.7,.7,1e-9,90],
        ['A','C',.7,.8,.8,.4,.9,1e-7,80],
        ['A','A',1,1,1,1,1,0,200],
        ['X','A',.9,.9,.9,.9,.9,1e-8,100],
    ], columns=['query','target','qtmscore','ttmscore','alntmscore','qcov','tcov','evalue','bits'])
    ids = {'A','B','C','D'}
    edges, stats = analyzer.clean_undirected_edges(fixture, ids, False, None, None)
    weights = {(r.source,r.target):r.weight_global_min_tm for r in edges.itertuples()}
    check('fixture_directional_mean_of_row_minima', math.isclose(weights['A','B'], .6) and math.isclose(weights['A','C'], .7), {'AB':weights['A','B'],'AC':weights['A','C']})
    check('fixture_self_and_unknown_removal', stats['self_hits_removed']==1 and stats['alignment_rows_missing_node']==1 and len(edges)==2)
    reciprocal, _ = analyzer.clean_undirected_edges(fixture, ids, True, None, None)
    coverage, _ = analyzer.clean_undirected_edges(fixture, ids, False, .5, None)
    score, _ = analyzer.clean_undirected_edges(fixture, ids, False, None, .75)
    check('fixture_reciprocal_filter', len(reciprocal)==1 and reciprocal.iloc[0]['target']=='B')
    check('fixture_coverage_filter', len(coverage)==1 and coverage.iloc[0]['target']=='B')
    check('fixture_optional_alignment_TM_filter', math.isclose(score.iloc[0]['weight_global_min_tm'], .8))
    graph, _, _ = analyzer.build_graph(pd.DataFrame({'name':['A','B','C','D']}), edges, 'name')
    partition = analyzer.la.find_partition(graph, analyzer.la.RBConfigurationVertexPartition, weights=edges.weight_global_min_tm.tolist(), resolution_parameter=1.0, n_iterations=-1, seed=1)
    check('fixture_weighted_Leiden_and_isolate_retention', graph.vcount()==4 and graph.degree(3)==0 and len(partition.membership)==4 and math.isfinite(partition.quality()) and list(partition.membership).count(partition.membership[3])==1)

    data = args.data_root
    primary = data/'outputs/primary'
    ref_nodes = read(primary/'leiden/cytoscape_nodes_analyzed.csv')
    ref_edges = read(primary/'leiden/edges_clean_undirected.csv')
    scientific = {}
    for label, run, nodes_path, edges_path in [
        ('preserved', args.preserved_run, primary/'leiden/cytoscape_nodes_analyzed.csv', primary/'leiden/edges_clean_undirected.csv'),
        ('fresh', args.fresh_run, args.fresh_run/'leiden/cytoscape_nodes_analyzed.csv', args.fresh_run/'leiden/edges_clean_undirected.csv'),
    ]:
        nodes, edges = read(nodes_path), read(edges_path)
        viral = read(run/'placement/viral_structural_placement.tsv')
        all_nodes = read(run/'placement/cytoscape_nodes_with_figure2b_classes.tsv')
        check(label+'_all_nodes_and_unique_IDs', set(nodes.name)==set(ref_nodes.name) and not nodes.name.duplicated().any())
        check(label+'_classifier_all_node_IDs', set(all_nodes.name)==set(nodes.name) and not all_nodes.name.duplicated().any())
        check(label+'_viral_IDs', set(viral.name)==set(nodes.loc[nodes.is_virus.str.lower().eq('true'),'name']) and not viral.name.duplicated().any())
        check(label+'_community_assignments_complete', nodes.community.ne('').all() and all_nodes.community.ne('').all())
        pairs = [tuple(sorted((a,b))) for a,b in zip(edges.source, edges.target)]
        weights = pd.to_numeric(edges.weight_global_min_tm)
        check(label+'_valid_undirected_edges', len(set(pairs))==len(pairs) and all(a!=b for a,b in pairs) and set(edges.source)|set(edges.target)<=set(nodes.name) and weights.between(0,1,inclusive='right').all())
        degree = Counter()
        strength = defaultdict(float)
        adjacent = defaultdict(list)
        for a,b,w in zip(edges.source,edges.target,weights):
            degree.update([a,b]);strength[a]+=w;strength[b]+=w
            adjacent[a].append((b,w));adjacent[b].append((a,w))
        check(label+'_classifier_degrees', all(int(r.degree)==degree[r.name] for r in all_nodes.itertuples()))
        check(label+'_classifier_strengths', all(math.isclose(float(r.weighted_degree),strength[r.name],rel_tol=1e-10,abs_tol=1e-12) for r in all_nodes.itertuples()))
        families = dict(zip(all_nodes.name,all_nodes.cazy_family_normalized))
        communities = dict(zip(nodes.name,nodes.community))
        thresholds = json.loads((workflow/'configs/historical/primary_classification_run.json').read_text())['thresholds']
        support_errors=[]
        for row in viral.itertuples():
            fs=defaultdict(float);fc=Counter()
            for other,weight in adjacent[row.name]:
                if communities[other]==row.community and families[other]:
                    fs[families[other]]+=weight;fc[families[other]]+=1
            ranked=sorted(fs,key=lambda f:(-fs[f],f));total=sum(fs.values());best=ranked[0] if ranked else ''
            second=fs[ranked[1]] if len(ranked)>1 else 0
            passed=bool(best) and fc[best]>=thresholds['family_support_min_neighbors'] and fs[best]/total>=thresholds['family_support_min_strength_fraction'] and (fs[best]-second)/total>=thresholds['family_support_min_margin']
            expected=best if passed and row.cazy_status=='CAZy_unassigned' else ''
            if expected!=row.directly_supported_cazy_family:support_errors.append(row.name)
        check(label+'_direct_support_recomputed_from_graph', not support_errors, support_errors)
        plot_data=plotter.prepare_data(run/'placement/viral_structural_placement.tsv')
        check(label+'_Figure2B_mapping_complete', len(plot_data)==len(viral) and plot_data[['annotation','placement','outcome']].notna().all().all())
        counts=Counter(nodes.community)
        scientific[label]=dict(nodes=len(nodes),edges=len(edges),isolates=sum(degree[n]==0 for n in nodes.name),communities=len(counts),multi_protein_communities=sum(n>1 for n in counts.values()),viral_proteins=len(viral),cazy_unassigned=int(viral.cazy_status.eq('CAZy_unassigned').sum()),direct_family_supported=int(viral.directly_supported_cazy_family.ne('').sum()),dominant_family_concordant=int(viral.unassigned_family_support_call.eq('direct_support_concordant_with_community').sum()),alternative_family_supported=int(viral.unassigned_family_support_call.eq('direct_support_discordant_with_community').sum()))
        if label=='preserved':
            diff=semantic_differences(read(primary/'placement/viral_structural_placement.tsv'),viral,'name')
            check('preserved_all_protein_columns_reproduced',not diff,dict(columns=len(viral.columns),differences=diff))
            for filename,key in [('figure2b_alluvial_flow_counts.tsv',['annotation','placement','outcome']),('figure2b_node_counts.tsv',['stage','category'])]:
                a=read(data/'figures/figure2b'/filename);b=read(run/'figure2b'/filename)
                a['_key']=a[key].agg('|'.join,axis=1);b['_key']=b[key].agg('|'.join,axis=1)
                differences=semantic_differences(a,b,'_key')
                check('preserved_'+filename,not differences,differences)
    fresh=args.fresh_run/'leiden'
    check('fresh_full_edge_table_byte_identical',sha(fresh/'edges_clean_undirected.csv')==sha(primary/'leiden/edges_clean_undirected.csv'))
    cfg=json.loads((fresh/'run_config.json').read_text())
    check('fresh_primary_parameters',cfg['seeds']==[1] and cfg['resolutions']==[1.0] and cfg['arguments']['n_iterations']==-1 and not cfg['arguments']['reciprocal_only'] and cfg['arguments']['analysis_min_coverage'] is None and cfg['arguments']['analysis_min_alntmscore'] is None and cfg['preferred_method_for_cytoscape_node_metrics']=='global_min_tm')
    check('fresh_isolates_retained',scientific['fresh']['isolates']==scientific['preserved']['isolates'])
    fresh_result=json.loads((args.fresh_run/'reproduction_result.json').read_text())
    comparison=fresh_result['comparison']
    result=dict(passed=all(c['status']=='PASS' for c in checks),python=sys.version,dependencies=dependencies,scientific_counts=scientific,fresh_partition_ARI_vs_preserved=comparison['ARI_vs_original_primary'],fresh_is_exact_historical_reproduction=fresh_result['exact_scientific_reproduction'],checks=checks)
    args.report.parent.mkdir(parents=True,exist_ok=True)
    args.report.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['checks','dependencies']},indent=2))
    for c in checks:
        if c['status']=='FAIL':print(json.dumps(c))
    raise SystemExit(0 if result['passed'] else 1)


if __name__=='__main__':
    main()
