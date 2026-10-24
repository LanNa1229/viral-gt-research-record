"""Publication figures from saved clustering outputs and UNCHANGED original UMAP coordinates."""
import collections
import common as c
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from PIL import Image,ImageOps
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.labelsize':11,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','axes.linewidth':.7})
PALETTE=['#3C6C9D','#C17042','#648C54','#9D688D','#C4A43E','#5A9EA2','#8C7565','#7D82B4','#AD6571','#839A88','#B99876','#786A96','#667D89','#9C963C','#B67D96','#61927C']
old=c.read(c.ROOT/c.CONFIG['previous_analysis']/'joint_umap_samples.csv');xy=np.array([[float(r['UMAP1']),float(r['UMAP2'])] for r in old]);byid={r['sample_id']:i for i,r in enumerate(old)}
rows=c.read(c.OUT/'tables/per_locus_selected287_clustering.csv');communities=c.read(c.OUT/'tables/community_summary.csv');ids=[r['viral_locus_id'] for r in rows];ix=[byid[sid] for sid in ids];ref=np.array([r['sample_type']=='reference_BGC' for r in old]);selected=set(ids);other=np.array([r['sample_type']=='viral_locus' and r['sample_id'] not in selected for r in old]);pad=np.ptp(xy,axis=0)*.055;low=xy.min(0)-pad;high=xy.max(0)+pad
exports=[]
def save(fig,name,size):
 for ext in ['png','svg']:fig.savefig(c.OUT/'figures'/(name+'.'+ext),dpi=600,facecolor='white')
 exports.append(dict(name=name,figure_inches=list(size),dpi=600));plt.close(fig)
for panel,field,title in [('C1','Leiden_cluster','Selected viral loci: 320D Leiden communities'),('C2','HDBSCAN_cluster','Selected viral loci: 320D HDBSCAN comparison')]:
 size=(9,6);fig=plt.figure(figsize=size);ax=fig.add_axes([.09,.16,.62,.70]);ax.set(xlabel='UMAP 1',ylabel='UMAP 2',xlim=(low[0],high[0]),ylim=(low[1],high[1]));ax.set_aspect('equal',adjustable='box')
 ax.scatter(*xy[ref].T,c='#DCDCDC',s=11,linewidths=0,alpha=.65,zorder=1);ax.scatter(*xy[other].T,c='#AEBAC3',s=13,marker='^',linewidths=0,alpha=.6,zorder=2)
 labels=sorted({r[field] for r in rows if r[field]!='Noise'});handles=[]
 for j,label in enumerate(labels):
  keep=np.array([i for i,r in enumerate(rows) if r[field]==label]);color=PALETTE[j%len(PALETTE)];ax.scatter(*xy[np.asarray(ix)[keep]].T,s=29,c=color,marker='^',edgecolors='white',linewidths=.25,zorder=4);handles.append(Line2D([],[],marker='^',linestyle='',color=color,markersize=6,label=f'{label} ({len(keep)})'))
 noise=np.array([i for i,r in enumerate(rows) if r[field]=='Noise'])
 if len(noise):ax.scatter(*xy[np.asarray(ix)[noise]].T,s=25,c='#262626',marker='x',linewidths=.8,zorder=5);handles.append(Line2D([],[],marker='x',linestyle='',color='#262626',markersize=6,label=f'Noise ({len(noise)})'))
 handles.extend([Line2D([],[],marker='o',linestyle='',color='#CCCCCC',markersize=5,label='Reference BGCs (2,502)'),Line2D([],[],marker='^',linestyle='',color='#AEBAC3',markersize=5,label='Other viral loci (269)')]);ax.legend(handles=handles,loc='upper left',bbox_to_anchor=(1.02,1),frameon=False,fontsize=9,title='Clustering assignments')
 fig.text(.025,.945,panel,fontsize=15,weight='bold');fig.text(.09,.945,title,fontsize=12);fig.text(.09,.902,'287 prior-selected loci | clusters inferred in original 320D | saved UMAP coordinates',fontsize=9,color='#555555')
 note='Baseline: union cosine kNN k=20; RB Leiden resolution=1.0; seed=42.' if panel=='C1' else 'HDBSCAN: cosine distance; min_cluster_size=15; min_samples=5 (includes self).'
 fig.text(.09,.075,note,fontsize=8.5,color='#444444');fig.text(.09,.038,'UMAP is visualization only; community colors imply no functional or pathway identity.',fontsize=8.5,color='#555555');save(fig,'Figure_'+panel+'_saved_UMAP',size)
# Cross-cohort membership relations on exactly the shared selected loci.
relations=c.read(c.OUT/'tables/cross_cohort_locus_memberships.csv');vcs=sorted(r['Leiden_cluster'] for r in communities)
size=(12,5.8);fig,axes=plt.subplots(1,2,figsize=size,gridspec_kw={'width_ratios':[1,1]});fig.subplots_adjust(left=.07,right=.92,bottom=.22,top=.77,wspace=.35)
for ax,name,title in zip(axes,['global_joint','saccharide_joint'],['Global joint: 287 viral + 2,502 references','Saccharide joint: 287 viral + 195 references']):
 data=[r for r in relations if r['comparison']==name];cols=sorted({r['comparison_cluster'] for r in data});counts=np.array([[sum(r['Leiden_cluster']==a and r['comparison_cluster']==b for r in data) for b in cols] for a in vcs]);fractions=counts/counts.sum(axis=1,keepdims=True)
 mesh=ax.pcolormesh(np.arange(len(cols)+1),np.arange(len(vcs)+1),fractions,cmap='Blues',vmin=0,vmax=1,edgecolors='white',linewidth=.5);ax.invert_yaxis();ax.set_xticks(np.arange(len(cols))+.5,cols,rotation=45,ha='right');ax.set_yticks(np.arange(len(vcs))+.5,vcs);ax.set(xlabel='Separate joint-community IDs',ylabel='Viral-only baseline community',title=title)
 for i in range(len(vcs)):
  for j in range(len(cols)):
   if counts[i,j]:ax.text(j+.5,i+.5,str(counts[i,j]),ha='center',va='center',fontsize=9,color='white' if fractions[i,j]>.6 else '#222222')
bar=fig.colorbar(mesh,ax=list(axes),fraction=.025,pad=.03);bar.set_label('Fraction of viral-only community');bar.solids.set_rasterized(False)
fig.text(.02,.945,'C3',fontsize=15,weight='bold');fig.text(.09,.945,'Viral-only and joint community relationships',fontsize=13);fig.text(.09,.89,'Cell labels: selected-locus counts; colors: row fractions. Joint communities without viral loci are omitted here.',fontsize=9,color='#555555');fig.text(.07,.065,'All clustering uses original 320D vectors. Reference annotations are post-fit descriptors.',fontsize=9,color='#555555');save(fig,'Figure_C3_community_relationships',size)
# Four complementary diagnostics with no parameter tuning from these plots.
size=(12,9);fig,axes=plt.subplots(2,2,figsize=size);fig.subplots_adjust(left=.08,right=.94,bottom=.13,top=.87,hspace=.52,wspace=.4)
ax=axes[0,0];ax.bar(vcs,[int(r['n_loci']) for r in communities],color=[PALETTE[i%len(PALETTE)] for i in range(len(vcs))]);ax.set(ylabel='Selected viral loci',title='Prespecified baseline community sizes');ax.tick_params(axis='x',rotation=45)
for i,r in enumerate(communities):ax.text(i,int(r['n_loci'])+1,r['n_loci'],ha='center',fontsize=9)
grid=c.read(c.OUT/'stability/selected287_grid_summary.csv');ks=c.CONFIG['grid']['k'];res=c.CONFIG['grid']['resolution'];mat=np.array([[np.median([float(r['ARI_to_cohort_baseline']) for r in grid if int(r['k'])==k and float(r['resolution'])==v and int(r['seed'])!=42]) for v in res] for k in ks]);ax=axes[0,1];mesh=ax.pcolormesh(np.arange(len(res)+1),np.arange(len(ks)+1),mat,cmap='viridis',vmin=0,vmax=1,edgecolors='white');ax.invert_yaxis();ax.set_xticks(np.arange(len(res))+.5,res);ax.set_yticks(np.arange(len(ks))+.5,ks);ax.set(xlabel='Leiden resolution',ylabel='kNN neighbors',title='Median ARI to baseline across five seeds')
for i in range(len(ks)):
 for j in range(len(res)):ax.text(j+.5,i+.5,f'{mat[i,j]:.2f}',ha='center',va='center',fontsize=9,color='white' if mat[i,j]<.55 else '#111111')
bar=fig.colorbar(mesh,ax=ax,fraction=.045);bar.solids.set_rasterized(False)
ax=axes[1,0];settings=['all556','selected_full15','selected_deduplicated','genome_balanced','genome_component_balanced'];titles=['All 556','Full 15 genes','Exact-context collapse','Genome-balanced','Sketch-component-balanced'];values=np.array([[float(r[name+'_best_Jaccard']) for name in settings] for r in communities]);mesh=ax.pcolormesh(np.arange(6),np.arange(len(vcs)+1),values,cmap='Blues',vmin=0,vmax=1,edgecolors='white');ax.invert_yaxis();ax.set_xticks(np.arange(5)+.5,titles,rotation=35,ha='right',fontsize=8.5);ax.set_yticks(np.arange(len(vcs))+.5,vcs);ax.set(title='Best community Jaccard on shared selected loci')
for i in range(len(vcs)):
 for j in range(5):ax.text(j+.5,i+.5,f'{values[i,j]:.2f}',ha='center',va='center',fontsize=8,color='white' if values[i,j]>.6 else '#222222')
ax=axes[1,1];hd=c.read(c.OUT/'stability/selected287_HDBSCAN_summary.csv');names=[f"{r['min_cluster_size']}/{r['min_samples']}" for r in hd];ax.bar(names,[float(r['noise_fraction'])*100 for r in hd],color='#808890');ax.set(ylim=(0,100),ylabel='Noise (%)',xlabel='min_cluster_size / min_samples',title='HDBSCAN density-based comparison');ax.tick_params(axis='x',rotation=35)
for i,r in enumerate(hd):ax.text(i,float(r['noise_fraction'])*100+2,f"{r['n_clusters']} groups",ha='center',fontsize=8)
fig.text(.02,.955,'C4',fontsize=15,weight='bold');fig.text(.09,.955,'Community size, parameter stability and confounder sensitivity',fontsize=13);fig.text(.08,.035,'Diagnostics are descriptive; stability thresholds were fixed before fitting. Agreement is not biochemical validation.',fontsize=9,color='#555555');save(fig,'Figure_C4_size_stability_QC',size)
# Preview uses saved exports, not additional scientific analysis.
thumbs=[]
for name in ['Figure_C1_saved_UMAP','Figure_C2_saved_UMAP','Figure_C3_community_relationships','Figure_C4_size_stability_QC']:
 with Image.open(c.OUT/'figures'/(name+'.png')) as im:thumbs.append(ImageOps.contain(im.convert('RGB'),(1200,900)))
canvas=Image.new('RGB',(2400,1800),'white')
for i,im in enumerate(thumbs):canvas.paste(im,((i%2)*1200,(i//2)*900))
canvas.save(c.OUT/'figures/preview_C1_C4.png');c.save('provenance/figure_exports.json',dict(plots=exports,coordinates_source=str(c.ROOT/c.CONFIG['previous_analysis']/'joint_umap_samples.csv'),coordinates_sha256=c.sha(c.ROOT/c.CONFIG['previous_analysis']/'joint_umap_samples.csv'),umap_refitted=False))
print('Saved C1-C4 PNG/SVG; reused original coordinates unchanged.')
