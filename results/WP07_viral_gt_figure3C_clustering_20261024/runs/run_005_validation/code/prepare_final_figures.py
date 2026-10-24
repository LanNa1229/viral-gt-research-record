import shutil
from pathlib import Path
OUT=Path(__file__).resolve().parents[1];dest=OUT/'runs/run_002_figures/figure_exports';dest.mkdir(parents=True,exist_ok=True)
for p in (OUT/'figures').iterdir():
 if p.is_file():shutil.copyfile(p,dest/p.name)
s=(OUT/'scripts/figures.py').read_text(encoding='utf-8-sig');s=s.replace('bottom=.13,top=.87,hspace=.52','bottom=.18,top=.87,hspace=.60')
s=s.replace("titles=['All 556','Full 15 genes','Exact-context collapse','Genome-balanced','Sketch-component-balanced']","titles=['All 556','Full 15','Deduplicated','Genome-balanced','Sketch-balanced']")
lines=s.splitlines()
for i,line in enumerate(lines):
 if line.startswith('for i,r in enumerate(hd):ax.text'):
  lines[i]='for i,r in enumerate(hd):ax.text(i,max(5,float(r["noise_fraction"])*100-7),str(r["n_clusters"])+" groups",ha="center",fontsize=8,color="white")'
(OUT/'scripts/figures_final.py').write_text('\n'.join(lines)+'\n',encoding='utf-8')
print('Preserved initial figure exports and prepared final layout correction; no rendering or fits executed.')
