#!/usr/bin/env python3
"""Run preserved code into a NEW directory; never run Foldseek automatically."""
import argparse,csv,hashlib,importlib.metadata,json,os,subprocess,sys
from pathlib import Path
def table(p,sep=','):
    with p.open(newline='') as f:return list(csv.DictReader(f,delimiter=sep))
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--data-root',type=Path,required=True);p.add_argument('--outdir',type=Path,required=True);p.add_argument('--mode',choices=['frozen','primary'],default='frozen');p.add_argument('--render-figure2b',action='store_true');a=p.parse_args()
    if a.outdir.exists():p.error('Use a new output directory; refusing to overwrite a previous run')
    wf=Path(__file__).resolve().parents[1];data=a.data_root.resolve();out=a.outdir.resolve()
    if out.is_relative_to(data) or out==wf or out.is_relative_to(wf/'original'):p.error('Output must be outside immutable data and original-code directories')
    out.mkdir(parents=True);commands=[]
    env=os.environ.copy();env['MPLBACKEND']='Agg';env['MPLCONFIGDIR']=str(out/'matplotlib_cache');env['PYTHONDONTWRITEBYTECODE']='1'
    def run(args,name):
        command=[sys.executable,'-B']+[str(x) for x in args];commands.append(dict(step=name,argv=command))
        with (out/(name+'.log')).open('w') as f:r=subprocess.run(command,stdout=f,stderr=subprocess.STDOUT,env=env)
        if r.returncode:raise RuntimeError(f'{name} failed; see {out/(name+".log")}')
    leiden=data/'outputs/primary/leiden';comparison={}
    if a.mode=='primary':
        original=json.loads((wf/'configs/historical/primary_leiden_run_config.json').read_text())['arguments'];args=[wf/'original/analyze_ssn_leiden.py']
        original.update(nodes=str(data/'inputs/metadata/005_GTA_node_metadata_completed_20260727.tsv'),alignments=str(data/'inputs/foldseek/005_GTA_all_vs_all_TM065_cov050_rich.id_corrected.m8'),outdir=str(out/'leiden'))
        for key,value in original.items():
            if value is None or value is False:continue
            args.append('--'+key.replace('_','-'))
            if value is not True:args.append(str(value))
        run(args,'primary_leiden');leiden=out/'leiden'
        ref={r['name']:r['community'] for r in table(data/'outputs/primary/leiden/cytoscape_nodes_analyzed.csv')};new=table(leiden/'cytoscape_nodes_analyzed.csv');from sklearn.metrics import adjusted_rand_score
        comparison=dict(ARI_vs_original_primary=adjusted_rand_score([ref[r['name']] for r in new],[r['community'] for r in new]),raw_community_label_differences=sum(ref[r['name']]!=r['community'] for r in new),method='new environment; original algorithm and primary settings; no Foldseek run')
    nodes=leiden/'cytoscape_nodes_analyzed.csv';edges=leiden/'edges_clean_undirected.csv';cfg=json.loads((wf/'configs/historical/primary_classification_run.json').read_text());args=[wf/'original/classify_viral_structural_placement_v3.py','--nodes',nodes,'--edges',edges,'--outdir',out/'placement','--node-id-column','name','--community-column','community','--virus-column','is_virus','--cazy-family-column','CAZy_family','--gtfam-column','GTFam','--weight-column','weight_global_min_tm']
    for key,value in cfg['thresholds'].items():args+=['--'+key.replace('_','-'),str(value)]
    run(args,'viral_placement')
    ref={r['name']:r for r in table(data/'outputs/primary/placement/viral_structural_placement.tsv','\t')};new=table(out/'placement/viral_structural_placement.tsv','\t');cols=['community','cazy_status','placement_major','unassigned_family_support_call','directly_supported_cazy_family','cazy_relation_to_community','figure2b_group'];delta=[]
    for r in new:
        for col in cols:
            if ref[r['name']][col]!=r[col]:delta.append(dict(name=r['name'],column=col,original=ref[r['name']][col],reproduced=r[col]))
    comparison['placement_scientific_field_differences']=len(delta);comparison['placement_rows']=len(new);comparison['placement_differences']=delta
    # Always derive the manuscript display from the ORIGINAL primary full graph.
    run([wf/'scripts/derive_figure2a_display.py','--data-root',data,'--outdir',out/'figure2a'],'figure2a_display')
    original_edges={(r['source'],r['target']) for r in table(data/'figures/figure2a/historical_source/figure2a_edges_display_backbone.csv')};newedges={(r['source'],r['target']) for r in table(out/'figure2a/reproduced_figure2a_display_edges.csv')};comparison['figure2a_edge_set_equal']=newedges==original_edges
    if a.render_figure2b:run([wf/'original/make_figure2b_alluvial.py',out/'placement/viral_structural_placement.tsv','--output-dir',out/'figure2b'],'figure2b_render')
    versions={name:importlib.metadata.version(name) for name in ['pandas','numpy','scipy','scikit-learn','igraph','leidenalg','matplotlib','seaborn']};result=dict(mode=a.mode,python=sys.version,versions=versions,commands=commands,comparison=comparison,exact_scientific_reproduction=(not delta and comparison['figure2a_edge_set_equal'] and comparison.get('ARI_vs_original_primary',1)==1),historical_environment=False)
    (out/'reproduction_result.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2));raise SystemExit(0 if result['exact_scientific_reproduction'] else 2)
if __name__=='__main__':main()
