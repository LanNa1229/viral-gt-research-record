#!/usr/bin/env python3
"""Portable edge derivation only; no reconstruction of historical colors or labels."""
import argparse,importlib.util,json
from pathlib import Path
import pandas as pd
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--data-root',required=True,type=Path);p.add_argument('--outdir',required=True,type=Path);a=p.parse_args()
    if a.outdir.exists():p.error('Use a new output directory')
    root=Path(__file__).resolve().parents[1];spec=importlib.util.spec_from_file_location('original_figure2a',root/'original/prepare_figure2a_cytoscape.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    e=m.read_table(a.data_root/'outputs/primary/leiden/edges_clean_undirected.csv');n=m.read_table(a.data_root/'outputs/primary/leiden/cytoscape_nodes_analyzed.csv')
    output,stats=m.build_backbone(e,n['name'].astype(str).tolist(),'source','target','weight_global_min_tm',3)
    a.outdir.mkdir(parents=True);output.to_csv(a.outdir/'reproduced_figure2a_display_edges.csv',index=False)
    (a.outdir/'reproduction.json').write_text(json.dumps(dict(stats=stats,role='later computational reproduction; not historical Cytoscape export',code_sha256=m.sha256(root/'original/prepare_figure2a_cytoscape.py')),indent=2)+'\n');print(json.dumps(stats))
if __name__=='__main__':main()
