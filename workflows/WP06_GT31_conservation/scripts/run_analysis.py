#!/usr/bin/env python3
"""Reproduce GT31 reference-core definition and F2Y1U4 mapping offline."""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))


def check(condition, message):
    if not condition:
        raise ValueError(message)


def verify_files():
    data_manifest = ROOT / 'INPUT_SHA256SUMS'
    if not data_manifest.is_file():
        raise SystemExit(
            'Required input data are distributed separately. Extract '
            'GT31_conservation_inputs_20261022.zip into this workflow directory '
            '(beside README.md), then run this command again. See README.md.')
    for manifest in [ROOT / 'SHA256SUMS', data_manifest]:
        for line in manifest.read_text().splitlines():
            digest, relpath = line.split('  ', 1)
            path = ROOT / relpath
            check(path.is_file(), f'Missing package file: {relpath}')
            check(hashlib.sha256(path.read_bytes()).hexdigest() == digest,
                  f'File differs from the validated package: {relpath}')


def class_counts(path):
    return dict(Counter(str(int(float(line[60:66])))
                        for line in path.read_text().splitlines()
                        if line.startswith('ATOM  ') and line[12:16].strip() == 'CA'))


def check_preserved_atoms(original, annotated):
    def atoms(path):
        return [line for line in path.read_text().splitlines()
                if line.startswith(('ATOM  ', 'HETATM'))]
    before, after = atoms(original), atoms(annotated)
    check(len(before) == len(after), 'Class annotation changed atom count')
    check(all(a[:60] == b[:60] and a[66:] == b[66:] for a, b in zip(before, after)),
          'Class annotation changed fields outside B factor')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', type=Path, default=ROOT / 'output')
    args = parser.parse_args()
    started = time.perf_counter()
    verify_files()
    out = args.outdir.resolve()
    check(not out.exists() or (out.is_dir() and not any(out.iterdir())),
          'Output directory must be new or empty; use --outdir with another path.')
    out.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault('MPLBACKEND', 'Agg')
    os.environ.setdefault('MPLCONFIGDIR', str(out / '.matplotlib_cache'))

    import numpy as np
    import pandas as pd
    from vgt_constraints.common import load_config, read_fasta, read_tsv, structure_qc
    from vgt_constraints.medoids import choose_medoids, standardize_pairwise
    from vgt_constraints.pan_core import compute_pan_core_workflow
    from make_query_core_classes import build as annotate_query
    from make_pymol_pml import build as make_view

    config = load_config(ROOT / 'config.yaml')
    members = read_tsv(ROOT / 'inputs/reference_members.tsv')
    check(len(members) == 17 and members.protein_id.is_unique, 'Expected 17 unique references')
    check('F2Y1U4' not in set(members.protein_id), 'Query must be excluded from reference selection')
    members['structure_path'] = members.structure_path.map(lambda p: str((ROOT / p).resolve()))
    members['foldmason_id'] = members.protein_id
    members['exclude'] = members.exclude.astype(str).str.lower().isin(['true', '1', 'yes'])
    qc = pd.DataFrame([structure_qc(row.structure_path, row.chain_id or None)
                       for row in members.itertuples(index=False)])
    members['mean_pLDDT'] = qc.mean_structure_bfactor
    members['n_residues'] = qc.n_residues
    members['eligible_for_medoid'] = (
        (members.mean_pLDDT >= config['quality']['min_model_plddt'])
        & (members.n_residues >= config['quality']['min_model_residues'])
        & ~members.exclude)
    # The archived alignments use biological IDs. Do not supply a duplicate
    # FoldMason alias to the original pan-core header resolver.
    members.drop(columns='foldmason_id').to_csv(out / 'prepared_members.tsv', sep='\t', index=False)

    print('1/3 Selecting representative structures from the complete reference matrix', flush=True)
    raw = read_tsv(ROOT / 'inputs/foldseek/reference17_pairwise.tsv')
    ids = set(members.protein_id)
    check(len(raw) == 17 * 17 and not raw.duplicated(['query', 'target']).any()
          and set(raw['query']) == ids and set(raw['target']) == ids,
          'Pairwise input must contain all 289 directional reference comparisons')
    pairwise = standardize_pairwise(raw)
    check(len(pairwise) == 136 and (pairwise.n_directional_rows == 2).all(),
          'Both directions are required for every non-self pair')
    medoids, centrality = choose_medoids(members, pairwise, config)
    expected_medoids = read_tsv(ROOT / 'expected/medoids.tsv')
    pd.testing.assert_frame_equal(
        medoids[['cluster_id', 'medoid_id']].reset_index(drop=True),
        expected_medoids.reset_index(drop=True))
    medoids.to_csv(out / 'medoids.tsv', sep='\t', index=False)
    centrality.to_csv(out / 'medoid_centrality.tsv', sep='\t', index=False)

    # Select rows from the archived reference19 MSA, retaining every column.
    aligns = ROOT / 'inputs/alignments'
    reference_paths = []
    for kind in ['aa', '3di']:
        records = read_fasta(aligns / f'C012_GT31_reference19_{kind}.fa')
        check(len(records) == 19, 'Expected original reference19 alignment')
        check(ids <= set(records), 'Reference IDs are missing from the alignment')
        path = out / f'reference17_common_{kind}.fa'
        path.write_text(''.join(f'>{pid}\n{records[pid]}\n' for pid in members.protein_id))
        check({len(records[pid]) for pid in ids} == {966}, 'Reference alignment columns changed')
        reference_paths.append(path)

    print('2/3 Defining reference-only core and mapping it onto F2Y1U4', flush=True)
    pan = out / 'pan_core'
    summary = compute_pan_core_workflow(
        members_path=out / 'prepared_members.tsv', medoids_path=out / 'medoids.tsv',
        reference_aa_path=reference_paths[0], reference_di_path=reference_paths[1],
        output_dir=pan, config=config, anchor_medoid_id='Q9N294', query_id='F2Y1U4',
        query_aa_path=aligns / 'C012_GT31_reference19_plus_F2Y1U4_aa.fa',
        query_di_path=aligns / 'C012_GT31_reference19_plus_F2Y1U4_3di.fa',
        query_structure=ROOT / 'inputs/structures/F2Y1U4.pdb')

    print('3/3 Writing class PDBs, PyMOL view, and validating against archived assignments', flush=True)
    query_class = pan / 'query_mapping/F2Y1U4_pan_core_classes.pdb'
    annotation = annotate_query(ROOT / 'inputs/structures/F2Y1U4.pdb',
                                pan / 'query_mapping/F2Y1U4_pan_core_mapping.tsv', query_class)
    check_preserved_atoms(ROOT / 'inputs/structures/F2Y1U4.pdb', query_class)
    counts = {'Q18515': class_counts(pan / 'Q18515_pan_core_classes.pdb'),
              'F2Y1U4': class_counts(query_class)}
    check(counts == {'Q18515': {'100': 98, '60': 49, '0': 242},
                     'F2Y1U4': {'100': 93, '60': 41, '0': 97}}, 'Class counts changed')
    comparisons = [
        ('pan_shared_core_columns.tsv', 'pan_shared_core_columns.tsv'),
        ('pan_core_medoid_anchor_metrics.tsv', 'pan_core_medoid_anchor_metrics.tsv'),
        ('query_mapping/F2Y1U4_pan_core_mapping.tsv', 'F2Y1U4_pan_core_mapping.tsv')]
    for current, archived in comparisons:
        pd.testing.assert_frame_equal(
            pd.read_csv(pan / current, sep='\t'),
            pd.read_csv(ROOT / 'expected' / archived, sep='\t'),
            check_dtype=False, check_exact=False, rtol=1e-10, atol=1e-10)
    expected_summary = json.loads((ROOT / 'expected/pan_core_run_summary.json').read_text())
    for key, value in expected_summary.items():
        if isinstance(value, float):
            check(np.isclose(summary[key], value, rtol=1e-10, atol=1e-10), f'Summary changed: {key}')
        else:
            check(summary[key] == value, f'Summary changed: {key}')
    make_view(pan / 'Q18515_pan_core_classes.pdb', query_class, out / 'view_pan_core_comparison.pml')
    versions = {name: importlib.metadata.version(name)
                for name in ['numpy', 'pandas', 'biopython', 'PyYAML', 'matplotlib', 'seaborn', 'scipy']}
    validation = {
        'status': 'PASS', 'core_calls_match_archive': True,
        'anchor_metrics_match_archive': True, 'query_mapping_matches_archive': True,
        'automatic_medoids_match_archive': True, 'input_sha256_verified': True,
        'query_only_b_factor_changed': annotation['only_B_factor_field_changed'],
        'class_counts': counts, 'summary': summary, 'python': sys.version,
        'packages': versions, 'elapsed_seconds': round(time.perf_counter() - started, 2),
        'upstream_search_and_alignment_rerun': False,
        'pymol_rendering_tested': False,
    }
    (out / 'validation.json').write_text(json.dumps(validation, indent=2) + '\n')
    print(json.dumps(validation, indent=2))


if __name__ == '__main__':
    main()
