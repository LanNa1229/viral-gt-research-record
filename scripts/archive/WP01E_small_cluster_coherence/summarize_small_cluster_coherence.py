#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import statistics
from collections import defaultdict, deque
from pathlib import Path

COLS = [
    'query','target','qlen','tlen','alnlen','qcov','tcov',
    'qstart','qend','tstart','tend','fident','evalue','bits',
    'qtmscore','ttmscore','alntmscore','lddt'
]


def num(x):
    try:
        return float(x)
    except Exception:
        return math.nan


def avg(values):
    vals = [x for x in values if not math.isnan(x)]
    return statistics.fmean(vals) if vals else math.nan


def med(values):
    vals = [x for x in values if not math.isnan(x)]
    return statistics.median(vals) if vals else math.nan


def fmt(x):
    if isinstance(x, float):
        return '' if math.isnan(x) else f'{x:.4f}'
    return x


def fasta_lengths(path):
    out, current, seq = {}, None, []
    def save():
        if current is not None:
            out[current] = len(''.join(seq).replace('-', '').replace('.', ''))
    with open(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith('>'):
                save()
                current = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
        save()
    return out


def components(nodes, adj):
    unseen, comps = set(nodes), []
    while unseen:
        start = unseen.pop()
        comp, q = {start}, deque([start])
        while q:
            x = q.popleft()
            for y in adj[x]:
                if y in unseen:
                    unseen.remove(y)
                    comp.add(y)
                    q.append(y)
        comps.append(comp)
    return sorted(comps, key=len, reverse=True)


def diameter(nodes, adj):
    if len(components(nodes, adj)) != 1:
        return None
    dmax = 0
    for start in nodes:
        dist, q = {start: 0}, deque([start])
        while q:
            x = q.popleft()
            for y in adj[x]:
                if y not in dist:
                    dist[y] = dist[x] + 1
                    q.append(y)
        dmax = max(dmax, max(dist.values()))
    return dmax


def open_wedges(nodes, adj):
    count = 0
    for center in nodes:
        nbrs = sorted(adj[center])
        for i in range(len(nbrs)):
            for j in range(i + 1, len(nbrs)):
                if nbrs[j] not in adj[nbrs[i]]:
                    count += 1
    return count


def label(size, observed_fraction, density, ncomp, diam, wedges):
    if observed_fraction < 0.90:
        return 'incomplete_pairwise_search'
    if size == 2:
        return 'supported_pair' if density == 1 else 'weak_pair'
    if ncomp > 1:
        return 'fragmented'
    if density == 1:
        return 'compact_clique'
    if density >= 0.80 and diam is not None and diam <= 2:
        return 'compact_dense'
    if wedges > 0:
        return 'chained_connected'
    return 'connected_sparse'


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--membership', required=True, type=Path)
    p.add_argument('--pairwise', required=True, type=Path)
    p.add_argument('--fasta', required=True, type=Path)
    p.add_argument('--outdir', required=True, type=Path)
    p.add_argument('--tm-threshold', type=float, default=0.50)
    p.add_argument('--coverage-threshold', type=float, default=0.50)
    p.add_argument('--lddt-threshold', type=float, default=0.0)
    a = p.parse_args()
    a.outdir.mkdir(parents=True, exist_ok=True)

    clusters, member_to_rep = defaultdict(list), {}
    with open(a.membership) as fh:
        for line in fh:
            rep, member = line.rstrip('\n').split('\t')[:2]
            clusters[rep].append(member)
            member_to_rep[member] = rep

    lengths = fasta_lengths(a.fasta)
    missing_lengths = sorted(set(member_to_rep) - set(lengths))
    if missing_lengths:
        raise SystemExit(f'Missing FASTA lengths: {missing_lengths[:10]}')

    directed = defaultdict(list)
    with open(a.pairwise) as fh:
        for line_no, row in enumerate(csv.reader(fh, delimiter='\t'), 1):
            if len(row) != len(COLS):
                raise SystemExit(f'{a.pairwise}:{line_no}: expected {len(COLS)} columns, got {len(row)}')
            r = dict(zip(COLS, row))
            q, t = r['query'], r['target']
            if q == t or q not in member_to_rep or t not in member_to_rep:
                continue
            if member_to_rep[q] != member_to_rep[t]:
                continue
            directed[tuple(sorted((q, t)))].append(r)

    pairs = {}
    pair_rows = []
    for (x, y), records in sorted(directed.items()):
        def score(r):
            return (min(num(r['qtmscore']), num(r['ttmscore'])),
                    min(num(r['qcov']), num(r['tcov'])), num(r['lddt']))
        best = max(records, key=score)
        qtm, ttm = num(best['qtmscore']), num(best['ttmscore'])
        qcov, tcov = num(best['qcov']), num(best['tcov'])
        tm_min, cov_min, lddt = min(qtm, ttm), min(qcov, tcov), num(best['lddt'])
        supported = (tm_min >= a.tm_threshold and cov_min >= a.coverage_threshold and
                     (a.lddt_threshold <= 0 or lddt >= a.lddt_threshold))
        out = {
            'cluster_representative': member_to_rep[x], 'protein_1': x, 'protein_2': y,
            'protein_1_length': lengths[x], 'protein_2_length': lengths[y],
            'direction_count': len(records), 'aligned_length': int(float(best['alnlen'])),
            'qcov': qcov, 'tcov': tcov, 'coverage_min': cov_min,
            'qtmscore': qtm, 'ttmscore': ttm, 'tm_min': tm_min,
            'tm_mean': avg([qtm, ttm]), 'alignment_tmscore': num(best['alntmscore']),
            'lddt': lddt, 'fident': num(best['fident']), 'evalue': best['evalue'],
            'bits': num(best['bits']), 'supported_edge': int(supported)
        }
        pairs[(x, y)] = out
        pair_rows.append(out)

    pair_fields = list(pair_rows[0]) if pair_rows else []
    with open(a.outdir / 'pairwise_within_cluster.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, pair_fields, delimiter='\t', lineterminator='\n')
        if pair_fields:
            w.writeheader()
            for r in pair_rows:
                w.writerow({k: fmt(v) for k, v in r.items()})

    summaries, nodes, missing = [], [], []
    for rep in sorted(clusters, key=lambda r: (-len(clusters[r]), r)):
        members = sorted(clusters[rep])
        n = len(members)
        expected = n * (n - 1) // 2
        adj = {m: set() for m in members}
        scores = {m: [] for m in members}
        rows = []

        for i in range(n):
            for j in range(i + 1, n):
                x, y = members[i], members[j]
                r = pairs.get(tuple(sorted((x, y))))
                if r is None:
                    missing.append((rep, x, y))
                    scores[x].append(0.0); scores[y].append(0.0)
                    continue
                rows.append(r)
                scores[x].append(r['tm_min']); scores[y].append(r['tm_min'])
                if r['supported_edge']:
                    adj[x].add(y); adj[y].add(x)

        observed = len(rows)
        observed_fraction = observed / expected if expected else 1.0
        supported = sum(r['supported_edge'] for r in rows)
        density = supported / expected if expected else 1.0
        comps = components(members, adj)
        diam = diameter(members, adj)
        wedges = open_wedges(members, adj)
        medoid = max(members, key=lambda m: (avg(scores[m]), m))
        cls = label(n, observed_fraction, density, len(comps), diam, wedges)

        tmvals = [r['tm_min'] for r in rows]
        lvals = [r['lddt'] for r in rows]
        cvals = [r['coverage_min'] for r in rows]
        avals = [r['aligned_length'] for r in rows]
        lens = [lengths[m] for m in members]

        summaries.append({
            'cluster_representative': rep, 'cluster_size': n,
            'mean_protein_length': avg(lens), 'median_protein_length': med(lens),
            'min_protein_length': min(lens), 'max_protein_length': max(lens),
            'expected_pairs': expected, 'observed_pairs': observed,
            'observed_pair_fraction': observed_fraction,
            'supported_pairs': supported, 'supported_pair_density': density,
            'median_tm_min': med(tmvals), 'min_tm_min': min(tmvals) if tmvals else math.nan,
            'mean_tm_min': avg(tmvals), 'median_lddt': med(lvals),
            'median_coverage_min': med(cvals), 'median_aligned_length': med(avals),
            'connected_components': len(comps), 'largest_component_size': len(comps[0]),
            'graph_diameter': '' if diam is None else diam, 'open_wedges': wedges,
            'original_rep_mean_tm': avg(scores[rep]),
            'structural_medoid': medoid, 'medoid_mean_tm': avg(scores[medoid]),
            'coherence_class': cls
        })
        for m in members:
            nodes.append({
                'node_id': m, 'cluster_representative': rep, 'cluster_size': n,
                'protein_length': lengths[m], 'is_original_representative': int(m == rep),
                'is_structural_medoid': int(m == medoid),
                'mean_tm_min_to_cluster': avg(scores[m]), 'cluster_coherence_class': cls
            })

    def write_dicts(path, rows):
        fields = list(rows[0]) if rows else []
        with open(path, 'w', newline='') as fh:
            w = csv.DictWriter(fh, fields, delimiter='\t', lineterminator='\n')
            if fields:
                w.writeheader()
                for r in rows:
                    w.writerow({k: fmt(v) for k, v in r.items()})

    write_dicts(a.outdir / 'cluster_coherence_summary.tsv', summaries)
    write_dicts(a.outdir / 'small_cluster_node_metadata.tsv', nodes)

    with open(a.outdir / 'missing_within_cluster_pairs.tsv', 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n')
        w.writerow(['cluster_representative','protein_1','protein_2'])
        w.writerows(missing)

    with open(a.outdir / 'supported_edges_for_cytoscape.tsv', 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n')
        w.writerow(['source','target','interaction','cluster_representative','tm_min','coverage_min','lddt','aligned_length'])
        for r in pair_rows:
            if r['supported_edge']:
                w.writerow([r['protein_1'], r['protein_2'], 'pairwise_structural_support',
                            r['cluster_representative'], fmt(r['tm_min']), fmt(r['coverage_min']),
                            fmt(r['lddt']), r['aligned_length']])

    print(f'Clusters: {len(clusters)}')
    print(f'Proteins: {len(member_to_rep)}')
    print(f'Observed within-cluster pairs: {len(pair_rows)}')
    print(f'Missing within-cluster pairs: {len(missing)}')
    print(f'Output: {a.outdir}')

if __name__ == '__main__':
    main()
