#!/usr/bin/env python3
"""RECONSTRUCTED plotting/assembly, 2026-10-24. Not historical executed code.

Uses preserved table values in figure_inputs.json. No structural searches,
model predictions, functional reassignment or image-derived measurements.
"""
import argparse,collections,json,math
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

p=argparse.ArgumentParser();p.add_argument('--source-data',type=Path,required=True);p.add_argument('--outdir',type=Path,required=True);p.add_argument('--only-figure1',action='store_true',help='Render only the Figure 1 schematic');a=p.parse_args()
D=json.loads(a.source_data.read_text(encoding='utf-8'));a.outdir.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.titlesize':10,'axes.labelsize':8,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white'})
BLUE='#347DA3';GOLD='#D39A25';GREEN='#429584';GREY='#9CA3AA';COLORS=[BLUE,GOLD,GREEN]
def save(fig,name):
    for ext in ['pdf','svg','png']:fig.savefig(a.outdir/(name+'.'+ext),dpi=300,bbox_inches='tight',pad_inches=.12)
    plt.close(fig)
def tag(ax,label):ax.text(-.16,1.10,label,transform=ax.transAxes,fontweight='bold',fontsize=13,va='top')
def clean(ax):ax.tick_params(labelsize=7);ax.set_axisbelow(True)
def box(ax,values,labels,ylabel,threshold=None,log=False):
    b=ax.boxplot(values,tick_labels=labels,patch_artist=True,widths=.45,showfliers=False,medianprops={'color':'#111111','linewidth':1.3})
    rng=np.random.default_rng(20261024)
    for i,(v,c) in enumerate(zip(values,COLORS),1):
        b['boxes'][i-1].set(facecolor=c,alpha=.25)
        ax.scatter(i+rng.uniform(-.17,.17,len(v)),v,s=3,c=c,alpha=.22,linewidths=0,rasterized=False)
        med=float(np.median(v));ax.annotate(f'{med:g}',(i,med),xytext=(0,8),textcoords='offset points',ha='center',fontsize=8,fontweight='bold',bbox={'facecolor':'white','alpha':.8,'edgecolor':'none','pad':1})
    if threshold is not None:ax.axhline(threshold,color=BLUE,linestyle='--',linewidth=.9)
    if log:ax.set_yscale('log');ax.set_yticks([50,100,200,500,1000],labels=['50','100','200','500','1000'])
    ax.set_ylabel(ylabel);clean(ax)

# Figure 1: explicit stratification instead of a functional KEEP/DROP decision.
fig,ax=plt.subplots(figsize=(10,7));ax.set(xlim=(0,10),ylim=(0,7));ax.axis('off')
def node(x,y,w,h,text,fill='#EAF1F5'):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=.08',facecolor=fill,edgecolor='#71818E',linewidth=.8));ax.text(x+w/2,y+h/2,text,ha='center',va='center',fontsize=9)
def arrow(x,y,u,v):ax.annotate('',xy=(u,v),xytext=(x,y),arrowprops={'arrowstyle':'->','color':'#425767','lw':1})
node(.35,5.45,2.35,1.0,'447-residue MGAT2\nFoldseek search\nBaseline source')
node(3.3,5.45,2.35,1.0,'Foldseek / BFVD retrieval\n1,727 preserved structures\nDownstream baseline')
node(6.25,5.45,3.3,1.0,'Permissive coarse clustering\n761 clusters\n115 multi-member + 646 singletons')
arrow(2.8,5.95,3.15,5.95);arrow(5.75,5.95,6.1,5.95)
node(.35,3.05,2.7,1.25,'Dominant cluster\n594 structural entries\nPrioritized for downstream\nstructure-similarity analysis')
node(3.65,3.05,2.7,1.25,'Other multi-member clusters\n114 clusters / 487 entries\n4 GT-like, 23 non-GT,\n87 unresolved representatives')
node(6.95,3.05,2.7,1.25,'Singletons\n646 structural entries\nLength / confidence QC\n93 reference-length + confident')
for x in [1.7,5,8.3]:arrow(7.9,5.35,x,4.42)
node(.35,.95,4.25,1.1,'Positive-control benchmark\n6 proteins / 7 GT structural entries\n7/7 recovered in a 483-entry spike-in cluster','#FFF4DA')
node(5.2,.95,4.45,1.1,'Independent evaluation retained\n63 small clusters / 383 structural pairs\nAll singleton and non-giant evidence preserved','#EBF4EF')
ax.text(.35,.2,'Structural clustering and QC strata do not establish biochemical function. No population is discarded solely by cluster size.',fontsize=9)
ax.set_title('Structure-guided retrieval and structural stratification',fontsize=15,pad=12)
save(fig,'Figure_1_reconstructed')
if a.only_figure1: raise SystemExit(0)

# Supplementary 2: data-based benchmark and direct alignment panels; qualitative
# panel B is the archived plate with its own provenance, not newly simulated.
fig=plt.figure(figsize=(10.5,12));gs=fig.add_gridspec(3,1,height_ratios=[1,1.25,2.1],hspace=.35)
ax=fig.add_subplot(gs[0]);ax.axis('off');tag(ax,'A')
ax.text(.02,.83,'Baseline: 1,727 entries â†’ 761 clusters; dominant cluster = 594',transform=ax.transAxes,fontsize=11)
ax.text(.02,.55,'Add seven GT structural entries representing six distinct proteins',transform=ax.transAxes,fontsize=11)
ax.text(.02,.24,'Spike-in cluster: 439 + 37 + 7 = 483 entries; GT recovery = 7/7',transform=ax.transAxes,fontsize=12,fontweight='bold',bbox={'facecolor':'#FFF0CB','edgecolor':'none','pad':7})
ax.text(.02,0,'7SP6_B and 7SP6_C are excluded from the GT denominator; two L230 models remain in that denominator.',transform=ax.transAxes,fontsize=8)
ax=fig.add_subplot(gs[1]);ax.axis('off');tag(ax,'B')
asset=a.source_data.parent/'assets/archived_supplementary_2B_plate.png'
if asset.exists():
    ax.imshow(plt.imread(asset));ax.set_title('Archived six-control structure plate',loc='left',fontsize=10)
else:
    ax.text(.05,.5,'Panel B plate withheld pending model/image redistribution clearance.\nNumerical panels A/C and their Source Data are preserved.',transform=ax.transAxes,fontsize=10)
ax.text(0,-.08,'Graphic key: A197 green; CvHAS cyan; gp110 magenta; A064R yellow; L230 salmon; A111/114R purple/blue.',transform=ax.transAxes,fontsize=8)
ax=fig.add_subplot(gs[2]);tag(ax,'C')
rows=D['controls'];queries=list(dict.fromkeys(r['query'] for r in rows));names={'2C0N':'A197','7SP6_A':'CvHAS','Che8_gp110_model':'gp110','Q84432_model':'A111/114R','pbcv1_a064r_1_model':'A064R','Q5UQC3_apo_model':'L230 (apo)'}
ys=[];labels=[]
for k,q in enumerate(queries):
    subset=[r for r in rows if r['query']==q]
    qlen=int(subset[0]['qlen']);base=(len(queries)-k-1)*3
    ax.hlines(base+.5,0,qlen,color='#E7EAEC',linewidth=9)
    for r,off,c in [(next(r for r in subset if r['comparison_type']=='representative'),1,BLUE),(next(r for r in subset if r['comparison_type']=='best_nonrepresentative_giant_member'),0,GOLD)]:
        y=base+off;start=int(r['qstart']);end=int(r['qend']);ax.hlines(y,start,end,color=c,linewidth=5)
        ax.text(945,y,f"{end-start+1}/{qlen}   {float(r['qcov']):.3f}   {float(r['tcov']):.3f}   {float(r['min_qt_tm']):.3f}",fontsize=7,va='center')
    ys.append(base+.5);labels.append(names[q])
ax.set_yticks(ys,labels=labels);ax.set(xlim=(0,1450),ylim=(-1,18),xlabel='Control query residue coordinate (colored span = qend âˆ’ qstart + 1)')
ax.set_xticks([0,200,400,600,800]);ax.text(945,17.5,'Span/length   qcov    tcov    minTM',fontsize=8,fontweight='bold')
ax.plot([],[],color=BLUE,lw=4,label='Versus representative');ax.plot([],[],color=GOLD,lw=4,label='Versus best non-representative giant member');ax.legend(loc='upper left',bbox_to_anchor=(0,1.13),ncol=2,frameon=False,fontsize=8);clean(ax)
save(fig,'Supplementary_Figure_2_reconstructed')

# Supplementary 3: all six panels, using model-level records or cluster-level
# summaries according to the original metric definition.
fig,axs=plt.subplots(3,2,figsize=(11,11));fig.subplots_adjust(hspace=.52,wspace=.38)
sizes=[int(r['cluster_size']) for r in D['cluster_sizes'] if int(r['cluster_size'])>1];ct=collections.Counter(sizes)
ax=axs[0,0];xs=sorted(ct);ax.bar(range(len(xs)),[ct[x] for x in xs],color=BLUE);ax.set_xticks(range(len(xs)),xs);ax.set(xlabel='Exact cluster size (structural entries)',ylabel='Number of clusters',title='Multi-member cluster size distribution')
for i,x in enumerate(xs):ax.text(i,ct[x]+.5,str(ct[x]),ha='center',fontsize=7)
tag(ax,'a');clean(ax)
master=D['master'];byrep=collections.defaultdict(list)
for r in master:
    if r['cluster_group']!='singleton':byrep[r['representative_id']].append(float(r['analysis_sequence_length']))
groups=[('2',2,2),('3',3,3),('4â€“5',4,5),('6â€“10',6,10),('11â€“50',11,50),('594',594,594)]
vals=[]
for label,lo,hi in groups:vals.append([float(np.median(v)) for rep,v in byrep.items() if lo<=len(v)<=hi])
ax=axs[0,1];bp=ax.boxplot(vals,tick_labels=[g[0] for g in groups],showfliers=False,patch_artist=True,medianprops={'color':'#222222'})
rng=np.random.default_rng(20261024)
for i,v in enumerate(vals,1):
    ax.scatter(i+rng.uniform(-.14,.14,len(v)),v,s=10,color=BLUE,alpha=.45)
    ax.scatter(i,float(np.mean(v)),s=19,marker='D',facecolor='white',edgecolor='#222222',linewidth=.6,zorder=5)
for b in bp['boxes']:b.set(facecolor='#DCE8EF')
ax.set(xlabel='Cluster-size class',ylabel='Cluster median model length (aa)',title='Cluster-level length distributions');tag(ax,'b');clean(ax)
ax=axs[1,0];out=D['outlier_counts'];x=np.arange(len(out));sh=[r['short_le102']/r['n'] for r in out];lg=[r['long_ge763']/r['n'] for r in out]
ax.bar(x-.18,sh,.36,color=BLUE,label='Short â‰¤102 aa');ax.bar(x+.18,lg,.36,color=GOLD,label='Long â‰¥763 aa')
for i,r in enumerate(out):
    ax.text(i-.18,sh[i]+.003,f"{r['short_le102']}/{r['n']}",ha='center',fontsize=6);ax.text(i+.18,lg[i]+.003,f"{r['long_ge763']}/{r['n']}",ha='center',fontsize=6)
ax.set_xticks(x,[r['size_class'] for r in out]);ax.set(xlabel='Cluster-size class',ylabel='Fraction of structural entries',ylim=(0,.20),title='Protein-weighted length outliers');ax.legend(frameon=False,fontsize=7);tag(ax,'c');clean(ax)
ax=axs[1,1];pairs=D['pairs'];special='A0A6J5N5Z7'
for is_special,marker in [(False,'^'),(True,'o')]:
    rs=[r for r in pairs if (r['cluster_representative']==special)==is_special]
    sc=ax.scatter([float(r['coverage_min']) for r in rs],[float(r['tm_min']) for r in rs],c=[float(r['lddt']) for r in rs],cmap='viridis',vmin=.2,vmax=1,s=13,marker=marker,edgecolor='#333333' if is_special else 'none',linewidth=.3)
ax.axhline(.5,color=GREY,ls='--',lw=.8);ax.axvline(.5,color=GREY,ls='--',lw=.8);ax.set(xlabel='min(qcov, tcov)',ylabel='min(qTM, tTM)',xlim=(0,1),ylim=(0,1),title='383 unique within-cluster pairs');fig.colorbar(sc,ax=ax,label='LDDT',fraction=.045);tag(ax,'d');clean(ax)
ax=axs[2,0];coh=D['coherence'];seen=collections.Counter()
for r in coh:
    size=int(r['cluster_size']);den=float(r['supported_pair_density']);key=(size,den);off=seen[key];seen[key]+=1
    if r['cluster_representative']==special:ax.scatter(size,den,s=60,marker='*',color=GOLD,edgecolor='black',zorder=4)
    else:ax.scatter(size+(off%7-3)*.045,den,s=13,color=GREY if den==0 else BLUE,alpha=.65)
ax.set(xlabel='Cluster size',ylabel='Supported-pair density',ylim=(-.04,1.1),title='63 clusters; 54 with no supported pair');tag(ax,'e');clean(ax)
ax=axs[2,1];counts=collections.Counter(r['panel_C_class'] for r in D['classes']);left=0
for name,c in [('GT-like',GOLD),('non-GT',BLUE),('ambiguous/unresolved',GREY)]:
    n=counts[name];ax.barh([0],[n],left=left,color=c,height=.45,label=f'{name}: {n}')
    ax.text(left+n/2,.28,str(n),ha='center',fontsize=9);left+=n
ax.set(xlim=(0,114),ylim=(-.6,.8),xlabel='Representatives (n = 114)',title='Final representative-level assignments');ax.set_yticks([]);ax.legend(loc='lower center',bbox_to_anchor=(.5,-.45),frameon=False,fontsize=7,ncol=1);tag(ax,'f');clean(ax)
save(fig,'Supplementary_Figure_3_reconstructed')

# Supplementary 4: corrected non-overlapping layout; actual medians are shown.
fig,axs=plt.subplots(2,2,figsize=(11,8.5));fig.subplots_adjust(hspace=.45,wspace=.5,bottom=.12)
gn=['giant_cluster','other_multimember','singleton'];labels=['Dominant\n(n=594)','Other multi-member\n(n=487)','Singleton\n(n=646)']
ax=axs[0,0];box(ax,[[float(r['analysis_sequence_length']) for r in master if r['cluster_group']==g] for g in gn],labels,'Structure-represented length (aa)',122,True);ax.set_title('Structure-represented length');tag(ax,'A')
ax=axs[0,1];box(ax,[[float(r['median_plddt']) for r in master if r['cluster_group']==g] for g in gn],labels,'Median residue pLDDT',70);ax.set(ylim=(0,105),title='Prediction confidence');tag(ax,'B')
ax=axs[1,0];metrics=['fasta_length','median_plddt','frac_plddt_lt70'];stats=[next(r for r in D['statistics'] if r['metric']==m and r['groups_compared']=='singleton vs other_multimember') for m in metrics]
for i,r in enumerate(stats):
    value=float(r['effect_size']);low=float(r['effect_ci_lower']);high=float(r['effect_ci_upper']);ax.errorbar(value,2-i,xerr=[[value-low],[high-value]],fmt='o',color=BLUE,capsize=3)
    ax.text(value,2-i+.18,f'{value:.3f}',ha='center',fontsize=8)
ax.axvline(0,color=GREY,ls='--',lw=.8);ax.set(yticks=[2,1,0],yticklabels=['Model length','Median pLDDT','Fraction pLDDT <70'],xlabel="Cliffâ€™s Î´ (95% bootstrap CI)",xlim=(-.85,.45),ylim=(-.5,2.6),title='Singleton versus other multi-member');tag(ax,'C');clean(ax)
ax=axs[1,1];strata=D['singleton_strata'];names=['Short + pLDDT <70','Short + pLDDT â‰¥70','Length â‰¥122 + pLDDT <70','Length â‰¥122 + pLDDT â‰¥70'];ns=[r['count'] for r in strata];ax.barh(range(4),ns,color=BLUE);ax.set_yticks(range(4),names);ax.invert_yaxis();ax.set(xlim=(0,390),xlabel='Singleton structural entries',title='Length and confidence strata (n = 646)')
for i,n in enumerate(ns):ax.text(n+5,i,f'{n} ({n/646:.1%})',va='center',fontsize=8)
tag(ax,'D');clean(ax)
save(fig,'Supplementary_Figure_4_reconstructed')
(a.outdir/'reconstructed_plot_environment.json').write_text(json.dumps({'classification':'RECONSTRUCTED','matplotlib':matplotlib.__version__,'numpy':np.__version__,'source_data':a.source_data.name,'new_structural_computation':False},indent=2))
print('Created four reconstructed figure sets (PDF/SVG/PNG).')
