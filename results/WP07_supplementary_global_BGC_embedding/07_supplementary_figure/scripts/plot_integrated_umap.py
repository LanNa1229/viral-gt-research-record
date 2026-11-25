"""Minimalist v2: saved UMAP coordinates with full axes and restrained hollow markers."""
import csv,hashlib,json,os,sys
from pathlib import Path
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[1]
os.environ['MPLCONFIGDIR']=str(OUT/'runtime/mplconfig');os.environ['MPLBACKEND']='Agg'
sys.path.append(str(ROOT/'analysis/bgcprophet_reference/.deps'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
import matplotlib.font_manager as fm
CONFIG=json.loads((OUT/'config.json').read_text())
fm.fontManager.addfont(CONFIG['font_path'])
assert fm.FontProperties(fname=CONFIG['font_path']).get_name()=='Arial'
plt.rcParams.update({'font.family':'Arial','font.size':8.5,'axes.labelsize':9,'axes.linewidth':.5,'axes.edgecolor':'black','axes.labelcolor':'black','text.color':'black','xtick.color':'black','ytick.color':'black','svg.fonttype':'none','pdf.fonttype':42,'ps.fonttype':42,'figure.facecolor':'white','savefig.facecolor':'white','axes.grid':False,'axes.titleweight':'normal'})
MATRICES={
 'protanopia':[[.152286,1.052583,-.204868],[.114503,.786281,.099216],[-.003882,-.048116,1.051998]],
 'deuteranopia':[[.367322,.860646,-.227968],[.280085,.672501,.047413],[-.011820,.042940,.968881]],
 'tritanopia':[[1.255528,-.076749,-.178779],[-.078411,.930809,.147602],[.004733,.691367,.303900]]
}
def read(n):
 with (OUT/n).open(encoding='utf-8',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def save(n,x):(OUT/n).write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
def write(n,rows):
 with (OUT/n).open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def rgb(h):return np.array([int(h[i:i+2],16)/255 for i in [1,3,5]])
def lin(a):return np.where(a<=.04045,a/12.92,((a+.055)/1.055)**2.4)
def enc(a):
 a=np.clip(a,0,1);return np.where(a<=.0031308,12.92*a,1.055*a**(1/2.4)-.055)
def color(h,mode):
 a=rgb(h)
 if mode in MATRICES:return enc(lin(a)@np.asarray(MATRICES[mode]).T)
 if mode=='grayscale':
  y=float(lin(a)@np.array([.2126729,.7151522,.0721750]));return enc(np.array([y,y,y]))
 return a
def build(data,mode='normal'):
 W=CONFIG['width_mm'];H=CONFIG['height_mm']
 fig=plt.figure(figsize=(W/25.4,H/25.4))
 ax=fig.add_axes([CONFIG['axes_left_fraction'],CONFIG['axes_bottom_fraction'],CONFIG['axes_width_fraction'],CONFIG['axes_max_height_fraction']])
 ax.set_xlim(CONFIG['axis_limits']['UMAP1']);ax.set_ylim(CONFIG['axis_limits']['UMAP2'])
 ax.set_aspect('equal',adjustable='box')
 ax.set_xlabel('UMAP 1',fontsize=CONFIG['axis_titles_pt'],labelpad=4)
 ax.set_ylabel('UMAP 2',fontsize=CONFIG['axis_titles_pt'],labelpad=4)
 ax.xaxis.set_major_locator(MaxNLocator(nbins=6,steps=[1,2,2.5,5,10]))
 ax.yaxis.set_major_locator(MaxNLocator(nbins=6,steps=[1,2,2.5,5,10]))
 ax.tick_params(axis='both',which='major',bottom=True,left=True,top=False,right=False,labelbottom=True,labelleft=True,length=2.5,width=.5,labelsize=CONFIG['tick_labels_pt'],pad=3,direction='out')
 for s in ax.spines.values():s.set_visible(True);s.set_color('black');s.set_linewidth(CONFIG['frame_linewidth_pt']);s.set_bounds(None,None)
 ax.grid(False)
 groups=[];glyphs=[]
 def draw(rows,name,marker,face,edge,area,lw,alpha,z):
  xy=np.array([[float(r['UMAP1']),float(r['UMAP2'])] for r in rows])
  coll=ax.scatter(xy[:,0],xy[:,1],s=area,marker=marker,facecolors=face,edgecolors=edge,linewidths=lw,alpha=alpha,zorder=z,rasterized=False)
  coll.set_gid('samples_'+name)
  groups.append(dict(group_id='samples_'+name,n_samples=len(rows),marker=marker,facecolor='none' if isinstance(face,str) else 'explicit palette RGB',edgecolor='none' if isinstance(edge,str) and edge=='none' else 'explicit neutral RGB',area_pt2=area,alpha=alpha,linewidth_pt=lw))
  for r in rows:glyphs.append(dict(sample_id=r['sample_id'],group_id='samples_'+name,rendering_group=r['rendering_group'],marker=r['plot_marker'],UMAP1=r['UMAP1'],UMAP2=r['UMAP2'],global_leiden_community=r['global_leiden_community'],original_class_labels_json=r['reference_class_labels_json']))
 for c in CONFIG['classes']:
  rr=[r for r in data if r['rendering_group']==c]
  draw(rr,c,'o',color(CONFIG['palette'][c],mode),'none',CONFIG['reference_area_pt2'],0,CONFIG['reference_alpha'],2)
 rr=[r for r in data if r['plot_marker']=='hollow_circle']
 draw(rr,'Hybrid','o','none',color(CONFIG['hybrid_color'],mode),CONFIG['hybrid_area_pt2'],CONFIG['hybrid_linewidth_pt'],CONFIG['hybrid_alpha'],3)
 rr=[r for r in data if r['sample_type']=='viral_locus']
 draw(rr,'SelectedViral','^','none',color(CONFIG['viral_color'],mode),CONFIG['viral_area_pt2'],CONFIG['viral_linewidth_pt'],CONFIG['viral_alpha'],4)
 handles=[Line2D([],[],marker='o',markersize=4,markerfacecolor=color(CONFIG['palette'][c],mode),markeredgewidth=0,linestyle='none',label=c) for c in CONFIG['classes']]
 handles.append(Line2D([],[],marker='o',markersize=4,markerfacecolor='none',markeredgecolor=color(CONFIG['hybrid_color'],mode),markeredgewidth=.6,linestyle='none',label='Hybrid\n(multiple labels)'))
 handles.append(Line2D([],[],marker='^',markersize=4.5,markerfacecolor='none',markeredgecolor=color(CONFIG['viral_color'],mode),markeredgewidth=.6,linestyle='none',label='Selected viral\nGT-associated\nloci (287)'))
 legend=fig.legend(handles=handles,loc='center left',bbox_to_anchor=(.835,.53),bbox_transform=fig.transFigure,fontsize=CONFIG['legend_pt'],title='Reference BGCs\n(2,502)',title_fontsize=CONFIG['legend_title_pt'],frameon=False,ncol=1,handlelength=.8,handletextpad=.45,labelspacing=.65,borderaxespad=0,borderpad=0,alignment='left')
 legend.get_title().set_fontweight('normal')
 fig.canvas.draw();renderer=fig.canvas.get_renderer()
 bb=ax.get_position();fw,fh=fig.canvas.get_width_height()
 extents=[]
 for kind,tt in [('axis_title',[ax.xaxis.label,ax.yaxis.label]),('x_tick',ax.get_xticklabels()),('y_tick',ax.get_yticklabels()),('legend',[legend.get_title(),*legend.get_texts()])]:
  for t in tt:
   if not t.get_visible():continue
   b=t.get_window_extent(renderer)
   # Matplotlib may instantiate locator ticks beyond the displayed limits.
   if kind=='x_tick' and not ax.get_xlim()[0]<=t.get_position()[0]<=ax.get_xlim()[1]:continue
   if kind=='y_tick' and not ax.get_ylim()[0]<=t.get_position()[1]<=ax.get_ylim()[1]:continue
   extents.append(dict(kind=kind,text=t.get_text(),font_pt=t.get_fontsize(),x0=float(b.x0),x1=float(b.x1),y0=float(b.y0),y1=float(b.y1)))
 receipt=dict(scientific_axes=len(fig.axes),fig_width_mm=W,fig_height_mm=H,canvas_pixels=[fw,fh],canvas_dpi=fig.dpi,axes_bounds_figure_fraction=list(bb.bounds),axes_width_fraction=bb.width,axes_width_mm=bb.width*W,axes_height_mm=bb.height*H,axis_limits={k:[float(x) for x in v] for k,v in [('UMAP1',ax.get_xlim()),('UMAP2',ax.get_ylim())]},visible_x_ticks=[float(x) for x in ax.get_xticks() if ax.get_xlim()[0]<=x<=ax.get_xlim()[1]],visible_y_ticks=[float(y) for y in ax.get_yticks() if ax.get_ylim()[0]<=y<=ax.get_ylim()[1]],all_four_spines_visible=all(s.get_visible() for s in ax.spines.values()),grid_visible=any(a.get_visible() for a in [*ax.get_xgridlines(),*ax.get_ygridlines()]),equal_axis_units=ax.get_aspect()==1.0,community_text_labels=0,community_leader_lines=0,community_annotation_table=False,text=extents)
 return fig,glyphs,groups,receipt
def main():
 data=read('FIGURE_SOURCE_DATA.tsv')
 for n in ['figures','qc']:(OUT/n).mkdir(exist_ok=True)
 save('SOFTWARE_VERSIONS.json',dict(python=sys.version,numpy=np.__version__,matplotlib=matplotlib.__version__,matplotlib_file=matplotlib.__file__,font_resolved=fm.findfont(fm.FontProperties(family='Arial')),backend=matplotlib.get_backend(),PDF_fonttype=42,SVG_fonttype='none',no_model_execution=True))
 fig,glyphs,groups,layout=build(data)
 stem=OUT/'figures/SuppFig_global_leiden_biosynthetic_classes_minimalist_v2'
 fig.savefig(stem.with_suffix('.pdf'),metadata={'Title':'Reference biosynthetic classes and selected viral GT-associated loci','Subject':'Presentation revision of saved joint UMAP; no new analysis'})
 fig.savefig(stem.with_suffix('.svg'),metadata={'Title':'Reference biosynthetic classes and selected viral GT-associated loci'})
 fig.savefig(stem.with_suffix('.png'),dpi=CONFIG['dpi'])
 fig.savefig(OUT/'qc/preview_183mm_150dpi.png',dpi=150)
 fig.savefig(OUT/'qc/preview_183mm_96dpi.png',dpi=96)
 plt.close(fig)
 write('PLOTTED_SAMPLE_RECEIPTS.tsv',glyphs);write('RENDER_GROUP_RECEIPTS.tsv',groups);save('LAYOUT_RECEIPTS.json',layout)
 for mode in [*MATRICES,'grayscale']:
  f,_,_,_=build(data,mode);f.savefig(OUT/('qc/'+mode+'_preview.png'),dpi=150);plt.close(f)
 save('PLOT_PROVENANCE.json',dict(version=CONFIG['version'],original_projection_points=3058,displayed_samples=2789,reference_BGCs=2502,single_label_colored_references=2035,neutral_hollow_hybrid_references=467,hollow_viral_triangles=287,original_UMAP_analysis=CONFIG['original_UMAP_analysis'],Leiden_archival_commit=CONFIG['Leiden_archival_commit'],Leiden_baseline=CONFIG['Leiden_baseline'],community_ids_unchanged_in_source_data=True,coordinate_strings_unchanged=True,original_class_memberships_unchanged=True,reference_original_classes_not_reassigned=True,hybrid_rendering_is_not_an_additional_biological_class=True,no_UMAP_refit=True,no_Leiden_rerun=True,no_model_inference=True,PDF_SVG_fully_vector=True,one_scientific_panel=True))
 print(json.dumps(dict(status='GENERATED',samples=len(glyphs),groups=len(groups),axes_width_fraction=layout['axes_width_fraction'],axes_size_mm=[layout['axes_width_mm'],layout['axes_height_mm']],dimensions_mm=[CONFIG['width_mm'],CONFIG['height_mm']]),indent=2),flush=True)
if __name__=='__main__':main()
