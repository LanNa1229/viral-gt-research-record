#!/usr/bin/env python3
"""Replot frozen HMM results plus validated AF3 update; no searches.
Run update_structural_analysis.py first. Requires matplotlib/numpy, Arial."""
from pathlib import Path
import csv,json,math,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
R=Path(__file__).resolve().parents[2];O=Path(__file__).resolve().parent
def read(p):return list(csv.DictReader(p.open(),delimiter='\t'))
def table(name,rows):
 fields=list(dict.fromkeys(k for r in rows for k in r))
 with (O/name).open('w') as f:
  w=csv.DictWriter(f,fields,delimiter='\t');w.writeheader();w.writerows(rows)
u=json.loads((O/'results/update_summary.json').read_text())
assert all(hashlib.sha256((R/n).read_bytes()).hexdigest()==v for n,v in u['source_hashes'].items())
for n in ['arial.ttf','arialbd.ttf','ariali.ttf','arialbi.ttf']:font_manager.fontManager.addfont('/mnt/c/Windows/Fonts/'+n)
plt.rcParams.update({'font.family':'Arial','font.size':7.5,'text.color':'black','axes.labelcolor':'black','xtick.color':'black','ytick.color':'black','axes.edgecolor':'black','axes.linewidth':.6,'xtick.major.width':.6,'ytick.major.width':.6,'xtick.major.size':3,'ytick.major.size':3,'pdf.fonttype':42,'svg.fonttype':'none','mathtext.fontset':'custom','mathtext.rm':'Arial','mathtext.it':'Arial:italic','mathtext.bf':'Arial:bold','mathtext.cal':'Arial','mathtext.sf':'Arial','mathtext.tt':'Arial','mathtext.fallback':None})
h=read(R/'results/hmm_chma68_domains.tsv');seen=set();domains=[]
raw=[l.split(maxsplit=22) for l in (R/'results/chma68.domtbl').read_text().splitlines() if l and not l.startswith('#')]
assert len(raw)==len(h)==207
for r,x in zip(h,raw):
 assert r['protein_id']==x[0] and float(r['full_evalue'])==float(x[6])
 assert math.isclose(float(r['profile_coverage']),(int(x[16])-int(x[15])+1)/int(x[5]))
 k=tuple(r[c] for c in ['protein_id','hmm_start','hmm_end','target_start','target_end'])
 if k not in seen:domains.append(r);seen.add(k)
lead=next(r for r in domains if r['protein_id']=='YP_001950147.1')
assert float(lead['full_evalue'])==.0063 and len({r['protein_id'] for r in domains})==191
tm=read(O/'results/foldseek_tm_updated.tsv');best={}
for r in tm:
 k=r['protein_id']
 if k not in best or float(r['qtmscore'])>float(best[k]['qtmscore']):best[k]=r
assert len(tm)==471 and len(best)==157
assert all(len([r for r in tm if r['protein_id']==pid])==3 for pid in best)
ctrl=[r for r in read(O/'results/controls_global.tsv') if r['query']!=r['target']]
assert len(ctrl)==6
fig=plt.figure(figsize=(180/25.4,95/25.4),facecolor='white')
ax=fig.add_axes([.087,.18,.365,.70]);bx=fig.add_axes([.60,.18,.365,.70])
for axis,label,title in [(ax,'A','Sequence homology'),(bx,'B','Structural similarity')]:
 axis.spines[['top','right']].set_visible(False)
 axis.text(-.18,1.10,label,transform=axis.transAxes,fontsize=10,weight='bold')
 axis.set_title(title,loc='left',fontsize=8.5,pad=10)
ax.scatter([100*float(r['profile_coverage']) for r in domains],[-math.log10(float(r['full_evalue'])) for r in domains],s=13,c='#74899A',alpha=.72,linewidths=0)
x=100*float(lead['profile_coverage']);y=-math.log10(float(lead['full_evalue']))
ax.scatter([x],[y],s=26,c='#B57966',edgecolors='black',linewidths=.35,zorder=4)
ax.annotate('ORF270\nE = 0.0063; coverage 12.5%',xy=(x,y),xytext=(33,3.1),fontsize=7,arrowprops={'arrowstyle':'-','lw':.5,'color':'black'})
ax.axhline(5,c='#999999',ls=(0,(4,3)),lw=.7)
ax.text(98,5.13,'E = $10^{-5}$',ha='right',fontsize=7)
ax.set(xlim=(0,100),ylim=(-2.3,5.8),xticks=[0,25,50,75,100],yticks=[-2,0,2,4,5])
ax.set_xlabel('ChmA HMM profile coverage (%)',fontsize=8)
ax.set_ylabel('−log$_{10}$(full-sequence E-value)',fontsize=8)
models=sorted(best.values(),key=lambda r:r['protein_id'])
jit=np.random.default_rng(20261117).uniform(-.14,.14,len(models))
bx.scatter(jit,[float(r['qtmscore']) for r in models],s=12,c='#74899A',alpha=.72,linewidths=0)
i=next(i for i,r in enumerate(models) if r['protein_id']=='YP_001950147.1')
score=float(models[i]['qtmscore']);assert score==u['orf270_max']
bx.scatter([jit[i]],[score],s=30,c='#B57966',edgecolors='black',linewidths=.35,zorder=5)
bx.annotate(f'ORF270 (AF3)\nTM = {score:.3f}',xy=(jit[i],score),xytext=(.27,.35),fontsize=7,arrowprops={'arrowstyle':'-','lw':.5,'color':'black'})
cxs=np.linspace(.9,1.1,6)
bx.scatter(cxs,[float(r['qtmscore']) for r in ctrl],s=24,c='#789586',edgecolors='black',linewidths=.3)
bx.set(xlim=(-.4,1.4),ylim=(0,1),yticks=np.arange(0,1.01,.2),xticks=[0,1],xticklabels=['phiRSL1 models','ChmA controls'])
bx.set_ylabel('TM-score normalized to resolved ChmA',fontsize=8)
for t in fig.findobj(matplotlib.text.Text):t.set_color('#000000')
for ext in ['pdf','svg','png']:fig.savefig(O/f'Supplementary_Figure_ChmA_updated.{ext}',dpi=600,facecolor='white')
fig.savefig(O/'preview.png',dpi=170)
table('plotting_data_panelA.tsv',[dict(r,x_profile_coverage_pct=100*float(r['profile_coverage']),y_neglog10_fullE=-math.log10(float(r['full_evalue'])),source='results/hmm_chma68_domains.tsv') for r in domains])
table('plotting_data_panelB.tsv',[dict(r,group='phiRSL1',display_x=float(jit[i]),plotted_tm=float(r['qtmscore']),normalization='resolved query ChmA',method='global TM-align',source='results/foldseek_tm_updated.tsv') for i,r in enumerate(models)]+[dict(r,group='ChmA control',display_x=float(x),plotted_tm=float(r['qtmscore']),normalization='resolved query ChmA',method='global TM-align',source='results/controls_global.tsv') for x,r in zip(cxs,ctrl)])
(O/'FIGURE_LEGEND.md').write_text("""# Supplementary Figure | Sequence and structural assessment of canonical ChmA homology in phiRSL1
**A,** ChmA profile-HMM search of the complete annotated phiRSL1 proteome (NC_010811.2; 343 proteins), using the frozen 68-homologue profile. Points represent 207 distinct reported domain alignments from 191 proteins; alternative domains retain their own coverage and share the protein's full-sequence E-value. Proteins without reported alignments are not assigned values. Coverage is the aligned HMM span divided by profile length. ORF270 (YP_001950147.1), the strongest weak sequence candidate, is highlighted in terracotta (E = 0.0063; coverage 12.5%). The dashed line marks E = 10⁻⁵; no protein meets this threshold.
**B,** Maximum global TM-align score for each of 157 phiRSL1 protein models across three resolved experimental ChmA chains (7SQR, 7SQV and 8FNE; 471 comparisons), including the newly supplied AF3 model of ORF270 (highlighted; maximum 0.131). Scores are normalized to the resolved ChmA query chain. Controls show six directed non-self global comparisons among the same three ChmA structures; reciprocal comparisons are not independent. All groups use global Foldseek TM-align. Horizontal jitter separates observations.
Models cover 157 of 343 proteins; 70 of the original 156 models have mean pLDDT below 70, whereas the added ORF270 model has mean pLDDT 90.2. No convincing canonical ChmA homolog was detected by sequence- or structure-based analysis. These results do not establish biological absence of a phage nucleus or demonstrate replacement of nuclear-shell protection by DNA glycosylation.
""",encoding='utf-8')
lines=['# Updated figure QC','',f"- Selected model: {u['model_path']}",f"- SHA-256: {u['model_sha256']}",f"- Selection: supplied deterministic preference and highest supplied ranking; sample 1. Supplied full-precision ranking = {u['ranking_score_user_full_precision']}; original summary_confidences.json rounds it to {u['ranking_score_metadata']}. Both retained.",f"- Sequence: exact match to all 265 aa of YP_001950147.1; mean CA pLDDT {u['mean_plddt']:.6f}.",'- Copied original CIF and summary metadata preserved under inputs/. PDB conversion preserves coordinates and confidence B-factors.','- Method: identical Foldseek executable and global --alignment-type 1 --exhaustive-search 1 -a 1 -e 100 --max-seqs 1000 --threads 4 for added model and recomputed controls. Exact argv, version, stdout and stderr retained in logs/. Original 468 global comparisons reused unchanged; no sequence searches or original-model searches rerun.','- References: unchanged resolved monomers, 7SQR/7SQV/8FNE, lengths 522/523/525; same MBP exclusion in 8FNE.',f"- Updated count: {u['models']} unique proteins, {u['comparisons']} comparisons. One representative model per protein. Only three added ORF270 comparisons. Six non-self global controls (nine including self, all retained).",'- ORF270 scores (ChmA-query normalization / target normalization):']
for r in u['orf270_scores']:lines.append(f"  - {r['query']}: {r['qtmscore']} / {r['ttmscore']}; RMSD {r['rmsd']}; lDDT {r['lddt']}; query/target coverage {r['qcov']}/{r['tcov']}.")
lines += [f"- Maximum ORF270 score: {u['orf270_max']} against {u['orf270_best_ref']}.",f"- New global query-normalized control range: {u['global_control_range']}. Previous figure used local controls and both-chain normalization; new figure uses six directed global query-normalized scores throughout. Old control values were not reused.",'- ORF333 remains the strongest original and overall model match (0.277, DNA topoisomerase II); it is not highlighted to prioritize ORF270.','- Interpretation unchanged: no convincing canonical ChmA homology; not proof of absence of a nucleus or of other compartmentalization functions.','- Panel A unchanged observed data: 207 distinct domains, 191 proteins; raw HMM coordinates/E-values cross-checked. No artificial observations. Six leave-one-out controls remain in original results (E range 1.5e-198 to 1.8e-24), not in the rendered plot.','- Original source hashes verified before plotting; originals and previous figure untouched. Source hashes are in results/update_summary.json.','- Presentation: 180 x 95 mm; Arial exclusively; all written text black; vector PDF/ editable-text SVG; PNG 600 dpi. Requested clutter removed.','- Reproduce computations: python update_structural_analysis.py; reproduce figure: python plot_chma_supplementary_updated.py (WSL environment from original METHODS; matplotlib/numpy/Biopython). No download or model generation.','- Visual inspection: pending.']
(O/'FIGURE_QC.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
print(json.dumps(dict(orf270_max=score,models=len(models),controls=u['global_control_range']),indent=2))
