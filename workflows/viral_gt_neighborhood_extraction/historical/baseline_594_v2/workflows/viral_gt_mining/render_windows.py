"""Render with the existing Windows scientific packages; no dependency install."""
import argparse
from pathlib import Path
import sys
p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--run',type=Path,required=True)
a=p.parse_args()
sys.path.insert(0,str(a.root/'analysis/bgcprophet_reference/.deps'))
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'src'))
from viral_gt_mining.analysis import analyze
analyze(a.root,a.run,a.run/'taxonomy_cache.tsv',make_plots=True)
