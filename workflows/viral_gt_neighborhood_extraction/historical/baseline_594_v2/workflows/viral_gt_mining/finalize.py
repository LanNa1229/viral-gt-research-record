"""Write a verified milestone receipt, output checksums and scientific caveats."""
import argparse
import collections
import importlib.metadata
import json
import subprocess
from pathlib import Path
from viral_gt_mining.extract import read_tsv,write_tsv,json_write,sha_file


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--task',type=Path,required=True)
    p.add_argument('--run-name',default='baseline_594_v2');a=p.parse_args();task=a.task;run=task/'runs'/a.run_name
    def readj(name):return json.loads((run/name).read_text())
    ext,inf,analysis,validation=[readj(n) for n in ['extraction_provenance.json','inference_provenance.json','analysis_provenance.json','validation.json']]
    assert validation['status']=='PASS' and inf['status']=='PASS'
    summaries=read_tsv(run/'three_level_summary.tsv')
    mapping=read_tsv(run/'anchor_mapping.tsv');unresolved=[r for r in mapping if r['status']=='unresolved']
    qc=read_tsv(run/'extraction_qc.tsv');reasons=collections.Counter(r['reason'] for r in qc if r['status']!='PASS')
    code_root=task/'checkout'
    code=[p for directory in ['src/viral_gt_mining','tests/viral_gt_mining','workflows/viral_gt_mining']
          for p in (code_root/directory).rglob('*') if p.is_file() and '__pycache__' not in p.parts]
    write_tsv(task/'implementation_checksums.tsv',[dict(path=str(p.relative_to(code_root)),sha256=sha_file(p),bytes=p.stat().st_size) for p in sorted(code)])
    for directory in [run,task/'runs/pilot_v1']:
        files=[p for p in directory.iterdir() if p.is_file() and p.name!='artifact_checksums.tsv']
        write_tsv(directory/'artifact_checksums.tsv',[dict(path=p.name,sha256=sha_file(p),bytes=p.stat().st_size) for p in sorted(files)])
    versions={x:importlib.metadata.version(x) for x in ['biopython','numpy','torch','fair-esm']}
    json_write(task/'software_versions.json',versions)
    git=dict(development_branch=subprocess.check_output(['git','-C',str(code_root),'branch','--show-current'],text=True).strip(),
        base_commit=subprocess.check_output(['git','-C',str(code_root),'rev-parse','HEAD'],text=True).strip(),
        status=subprocess.check_output(['git','-C',str(code_root),'status','--short'],text=True),
        original_status=subprocess.check_output(['git','-C','/home/rahnn/projects/viral-gt-research-record','status','--short'],text=True))
    json_write(task/'git_receipt.json',git)
    lines=['# Baseline 594-member viral GT neighborhood run','',
        'The original WP01B cluster A0A3Q8Q3U2 contains exactly 594 unique anchors, including the representative.',
        'Source: `/home/rahnn/projects/viral-gt-research-record/results/WP01B_coarse_clustering/canonical/cluster_membership.tsv`.',
        'Headerless column 1 is representative; column 2 is accession. Baseline SHA256: `cd731e42c1bffd7946427540d4cb61d72f706f494a93df2a8be11831f554fd2f`.',
        'This is the original baseline, not a spike-in rerun. The seven-anchor pilot is a separate cohort.','',
        f'Latest full run: [`runs/{a.run_name}`](runs/{a.run_name}/).',
        '[Commands, dependencies, methods, output schemas and interpretation](checkout/workflows/viral_gt_mining/README.md).',
        '[Frozen inference audit](audit_inference.md). [Exact cohort receipt](inputs/cohort_validation.json).','',
        '## Executed and verified','',
        f'- Input anchors: {ext["input_rows"]}; mapped anchors: {ext["successful_input_rows"]}; unresolved anchors: {ext["unresolved_input_rows"]}.',
        f'- Occurrence-specific windows: {ext["windows"]}; complete 15-gene windows: {ext["complete_windows"]}; gene rows: {ext["gene_rows"]}.',
        f'- Frozen inferences: {inf["n_success"]}; excluded windows: {inf["n_excluded"]}. Overlapping windows remain separate.',
        f'- Independent validation: {validation["checks"]} checks passed; {validation["historical_artifacts_unchanged"]} historical artifacts retain their recorded content hashes.',
        f'- Extraction: {ext["wall_seconds"]:.2f} s; inference/encoding: {inf["elapsed_seconds"]:.2f} s; peak inference RSS: {inf["peak_process_rss_kib"]/1024:.1f} MiB.',
        f'- ESM unique proteins: {inf["encoding"]["unique_sequences"]}; new encodings: {inf["encoding"]["newly_encoded"]}; cached: {inf["encoding"]["historical_cache_hits"]+inf["encoding"]["persistent_cache_hits"]}.',
        f'- Development branch: `{git["development_branch"]}`. Original repository status is {"clean" if not git["original_status"].strip() else "recorded in git_receipt.json"}. No push/publication performed.','',
        '## Classification counts','',
        '| Unit | Total | Analyzed at least once | Positive | Fully analyzed negative | Unknown | Fraction among resolved |',
        '|---|---:|---:|---:|---:|---:|---:|']
    for s in summaries:
        f=float(s['positive_fraction_resolved']) if s['positive_fraction_resolved'] else None
        lines.append('| '+ ' | '.join([s['level'],s['n_units_total'],s['n_any_analyzed'],s['n_positive'],s['n_negative_fully_analyzed'],s['n_unknown'],f'{f:.1%}' if f is not None else 'NA'])+' |')
    lines+=['', 'A locus/genome is positive if any analyzed constituent window meets the original classifier Saccharide threshold. A negative requires all constituent windows to be analyzed; otherwise its outcome remains unknown.',
        'Identical nucleotide-sequence hash plus source taxid collapses GenBank/RefSeq aliases for genome counts. Loci merge only connected overlapping nucleotide intervals after inference.',
        'Unmapped anchors have no known genome and remain outside genome denominators. Fractions describe this selected structural cohort, not all-virus prevalence. Wilson intervals are descriptive because observations overlap and were selected.',
        '', '## Mapping and exclusion scope','',
        'All remote mappings are bounded by explicit UniProt/UniSave genomic cross-references and exact source-CDS protein sequence matches; this is not an exhaustive public-genome/IPG search.',
        'Every annotated CDS is retained in genomic order. Missing/uncertain translations exclude inference without shifting gene positions. Linear/contig ends remain incomplete; no wrapping or compensation occurs.','']
    for reason,n in sorted(reasons.items()):lines.append(f'- {n} excluded windows: `{reason}`')
    lines+=['','Unresolved anchor details are in `anchor_mapping.tsv`, supplemented by `inputs/retrieval_anchor_status.tsv`.','',
        '## Pilot and model interpretation','',
        'The pilot had 7/7 mapped and inferred windows, six complete and one 11-gene boundary case; 544 independent pilot checks passed. S-H68 QBP05883.1 gives Saccharide 0.999795. The five Catovirus windows span Saccharide scores 0.024441–0.996528, with only two receiving the Saccharide label. The Catovirus ARF08732.1-centered window is not the historical ARF08726–ARF08739 input; see `runs/pilot_v1/manual_region_comparison.tsv`.',
        'Primary representation uses actual window genes, no artificial padding/repetition, the exact frozen checkpoint and the compatible 2,502-BGC CORE_ONLY reference. Seven sigmoid outputs need not sum to one. Original all-actual-position mean pooling equals actual-gene pooling in this input. Reference startup identity/proof gates protect comparability.',
        'Classification and 320D cosine/P@k retrieval remain separate. P@10 enrichment is descriptive over 195/2502 Saccharide-labeled references; no binary localization cutoff or taxonomic enrichment test was chosen. Neither result confirms glycan biosynthesis.',
        'Proteins above 1,022 aa retain complete source sequences but use the original first-1,022-residue encoder convention. Viral domain shift and supervised training-reference geometry remain limitations.','',
        '## Files','',
        f'- [Anchor mappings](runs/{a.run_name}/anchor_mapping.tsv)',
        f'- [Window manifest](runs/{a.run_name}/viral_gt_windows.tsv)',
        f'- [Ordered window genes](runs/{a.run_name}/viral_gt_window_genes.tsv)',
        f'- [Full source protein FASTA](runs/{a.run_name}/viral_gt_window_proteins.faa)',
        f'- [Predictions](runs/{a.run_name}/viral_bgc_predictions.tsv)',
        f'- [Embedding manifest](runs/{a.run_name}/embedding_manifest.tsv)',
        f'- [320D localization](runs/{a.run_name}/saccharide_localization.tsv)',
        f'- [Candidate loci](runs/{a.run_name}/candidate_loci.tsv)',
        f'- [Taxonomy summaries](runs/{a.run_name}/viral_taxonomy_summary.tsv)',
        f'- [Independent validation](runs/{a.run_name}/validation.json)',
        f'- [Output checksums](runs/{a.run_name}/artifact_checksums.tsv)',
        '', '## Figures','',
        f'![Reference/viral overlay](runs/{a.run_name}/reference_viral_embedding_overlay.png)',
        f'![Saccharide-associated fractions](runs/{a.run_name}/saccharide_associated_fractions.png)',
        f'![Taxonomic distribution](runs/{a.run_name}/viral_taxonomic_distribution.png)',
        '', 'PNG exports are 300 dpi and paired PDFs are vector. Visual review status is recorded alongside the figures. Quantitative retrieval uses the original 320D space, never PCA coordinates.','']
    (task/'README.md').write_text('\n'.join(lines))
    print(json.dumps(dict(report=str(task/'README.md'),run=str(run),mapped=ext['successful_input_rows'],inferred=inf['n_success'],validation=validation),indent=2))


if __name__=='__main__':main()
