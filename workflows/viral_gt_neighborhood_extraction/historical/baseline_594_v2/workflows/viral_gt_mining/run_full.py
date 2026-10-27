"""Resume the prespecified baseline cohort after a verified independent pilot."""
import argparse
import json
import os
import resource
import shutil
import sys
import time
from pathlib import Path
from viral_gt_mining.extract import extract,read_tsv,write_tsv,json_write,sha_file
from viral_gt_mining.inference import run_inference,prepare_windows,seed_cache
from viral_gt_mining.taxonomy import resolve
from viral_gt_mining.analysis import analyze
from validate_run import validate


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',required=True);p.add_argument('--task',required=True)
    p.add_argument('--run-name',default='baseline_594_v2');a=p.parse_args();root=Path(a.root);task=Path(a.task)
    started=time.perf_counter();inputs=task/'inputs';run=task/'runs'/a.run_name
    assert json.loads((task/'runs/pilot_v1/validation.json').read_text())['status']=='PASS','Independent pilot gate'
    assert json.loads((inputs/'cohort_validation.json').read_text())['selected_count']==594,'Baseline population gate'
    assert json.loads((task/'retrieval_progress.json').read_text())['stage']=='retrieval_complete','Retrieval not complete'
    run.mkdir(parents=True,exist_ok=True)
    snapshot_manifest=[]
    for directory in ['src/viral_gt_mining','tests/viral_gt_mining','workflows/viral_gt_mining']:
        for source in (task/'checkout'/directory).rglob('*'):
            if not source.is_file() or '__pycache__' in source.parts:continue
            rel=source.relative_to(task/'checkout');dest=run/'implementation_snapshot'/rel
            dest.parent.mkdir(parents=True,exist_ok=True)
            if dest.exists():assert sha_file(dest)==sha_file(source),'Code changed: choose a new versioned run name'
            else:shutil.copyfile(source,dest)
            snapshot_manifest.append(dict(path=str(dest.relative_to(run)),sha256=sha_file(dest)))
    for source in inputs.iterdir():
        if not source.is_file():continue
        dest=run/'input_snapshot'/source.name;dest.parent.mkdir(exist_ok=True)
        if dest.exists():assert sha_file(dest)==sha_file(source),'Inputs changed: choose a new versioned run name'
        else:shutil.copyfile(source,dest)
        snapshot_manifest.append(dict(path=str(dest.relative_to(run)),sha256=sha_file(dest)))
    write_tsv(run/'run_snapshot_checksums.tsv',snapshot_manifest)
    paths=[r['path'] for r in read_tsv(inputs/'full_genbanks.tsv')]
    aliases={}
    alias_rows=read_tsv(inputs/'retrieved_aliases.tsv')
    if (inputs/'supplement_aliases.tsv').exists():
        alias_rows+=read_tsv(inputs/'supplement_aliases.tsv')
    if (inputs/'supplement_genbanks.tsv').exists():
        paths+= [r['path'] for r in read_tsv(inputs/'supplement_genbanks.tsv')]
    for r in alias_rows:
        assert sha_file(r['evidence_path'])==r['evidence_sha256'],'Alias evidence changed'
        aliases.setdefault(r['input_accession'],[]).append(r)
    json_write(task/'full_run_progress.json',dict(stage='extraction',run=str(run),pid=os.getpid()))
    extraction=extract(inputs/'baseline_594_anchors.txt',None,paths,run,task/'cache',aliases)
    # Attach granular retrieval causes without treating unresolved anchors as windows.
    statuses={r['input_accession']:r for r in read_tsv(inputs/'retrieval_anchor_status.tsv')}
    if (inputs/'supplement_status.tsv').exists():
        statuses.update({r['input_accession']:r for r in read_tsv(inputs/'supplement_status.tsv')})
    mapping=read_tsv(run/'anchor_mapping.tsv')
    for r in mapping:
        m=statuses[r['input_accession']]
        r['retrieval_status']=m['status'];r['retrieval_reason']=m.get('reason','')
    write_tsv(run/'anchor_mapping.tsv',mapping)
    eligible,qc,sequences=prepare_windows(run)
    historical,_=seed_cache(root)
    new_sequences={k:v for k,v in sequences.items() if k not in historical}
    estimate=dict(windows=len(eligible),unique_sequences=len(sequences),historical_sequence_hits=len(sequences)-len(new_sequences),
        potentially_new_sequences=len(new_sequences),long_proteins_truncated=sum(len(v)>1022 for v in new_sequences.values()),
        output_array_bytes=len(eligible)*320*4,unique_sequence_vector_bytes=len(sequences)*320*4,
        cpu_threads=4,GPU_available=False,esm_planning_seconds_low=len(new_sequences)*.2,
        esm_planning_seconds_high=len(new_sequences)*2.,
        estimate_basis='pilot two new ESM calls averaged 0.218 seconds; broad planning range allows longer proteins; shared caches may reduce work',
        source_annotation_bytes=sum(Path(p).stat().st_size for p in paths),
        caveat='not a calibrated length-stratified benchmark; actual elapsed time recorded after completion')
    json_write(run/'resource_estimate.json',estimate)
    print(json.dumps(dict(extraction=extraction,resource_estimate=estimate),indent=2),flush=True)
    json_write(task/'full_run_progress.json',dict(stage='inference',run=str(run),extraction=extraction,pid=os.getpid()))
    inference=run_inference(root,run,esm_batch_size=1,with_localization=True,sequence_cache_root=task/'cache/esm')
    tax=resolve(task,run)
    json_write(task/'full_run_progress.json',dict(stage='analysis',run=str(run),inference_success=inference['n_success'],pid=os.getpid()))
    # The existing WSL environment has inference dependencies but no matplotlib.
    # Scientific plotting uses the existing Windows packages via render_windows.py.
    analysis=analyze(root,run,run/'taxonomy_cache.tsv',make_plots=False)
    validation=validate(root,run)
    result=dict(stage='complete',status='PASS',run=str(run),extraction=extraction,
        inference_success=inference['n_success'],inference_excluded=inference['n_excluded'],taxonomy=tax,
        analysis=analysis,validation=validation,elapsed_seconds=time.perf_counter()-started,
        peak_process_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    json_write(task/'full_run_progress.json',result);print(json.dumps(result,indent=2),flush=True)
    json_write(run/'pipeline_provenance.json',result)


if __name__=='__main__':main()
