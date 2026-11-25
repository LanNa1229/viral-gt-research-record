"""Two local figure designs; only saved coordinates and frozen table labels are used."""
import csv,hashlib,itertools,json,os,platform,sys
from pathlib import Path
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[1]
os.environ['MPLCONFIGDIR']=str(OUT/'runtime/mplconfig');os.environ['MPLBACKEND']='Agg'
sys.path.append(str(ROOT/'analysis/bgcprophet_reference/.deps'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator
from PIL import Image,ImageDraw,ImageFont
C=json.loads((OUT/'config.json').read_text())
for p in [C['font_path'],C['bold_font_path']]:fm.fontManager.addfont(p)
plt.rcParams.update({'font.family':'Arial','font.size':7.5,'axes.linewidth':.5,'axes.labelcolor':'#202020','text.color':'#202020','xtick.color':'#202020','ytick.color':'#202020','svg.fonttype':'none','pdf.fonttype':42,'ps.fonttype':42,'figure.facecolor':'white','savefig.facecolor':'white','axes.grid':False,'axes.titleweight':'normal','svg.hashsalt':'viral-gt-figure3C-comparison-20261028'})
MATRICES={
 'protanopia':[[.152286,1.052583,-.204868],[.114503,.786281,.099216],[-.003882,-.048116,1.051998]],
 'deuteranopia':[[.367322,.860646,-.227968],[.280085,.672501,.047413],[-.011820,.042940,.968881]],
 'tritanopia':[[1.255528,-.076749,-.178779],[-.078411,.930809,.147602],[.004733,.691367,.303900]]}
def read(n):
 with (OUT/n).open(encoding='utf-8',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def save(n,x):(OUT/n).write_text(json.dumps(x,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
def write(n,rows):
 with (OUT/n).open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rgb(h):return np.array([int(h[i:i+2],16)/255 for i in [1,3,5]])
def lin(a):return np.where(a<=.04045,a/12.92,((a+.055)/1.055)**2.4)
def enc(a):
 a=np.clip(a,0,1);return np.where(a<=.0031308,12.92*a,1.055*a**(1/2.4)-.055)
def shown(h,mode):
 a=rgb(h)
 if mode in MATRICES:return enc(lin(a)@np.asarray(MATRICES[mode]).T)
 if mode=='grayscale':return enc(np.repeat(float(lin(a)@np.array([.2126729,.7151522,.0721750])),3))
 return a
def lab(a):
 xyz=lin(a)@np.array([[.4124564,.3575761,.1804375],[.2126729,.7151522,.0721750],[.0193339,.1191920,.9503041]]).T/np.array([.95047,1,1.08883]);d=6/29
 f=np.where(xyz>d**3,np.cbrt(xyz),xyz/(3*d*d)+4/29)
 return np.array([116*f[1]-16,500*(f[0]-f[1]),200*(f[1]-f[2])])
def build(data,version,mode='normal'):
 W=C['width_mm'];H=C[version+'_height_mm'];fig=plt.figure(figsize=(W/25.4,H/25.4));fig.dpi=100
 axes=[];colls=[];glyphs=[];groups=[];texts=[];symbols=[]
 def text(x,y,s,size,weight='normal'):
  t=fig.text(x/W,y/H,s,fontsize=size,ha='left',va='center',fontweight=weight);texts.append(t);return t
 def symbol(x,y,marker,hc,fill='solid',size=4.2,alpha=1):
  col=shown(hc,mode);face='none' if fill=='none' else col
  edge=col if fill=='none' else ('#6F6F6F' if marker=='o' else 'white')
  lw=.65 if fill=='none' else (.2 if marker=='o' else .25)
  l=Line2D([x/W],[y/H],marker=marker,markersize=size,markerfacecolor=face,markeredgecolor=edge,markeredgewidth=lw,alpha=alpha,linestyle='none',transform=fig.transFigure)
  fig.add_artist(l);symbols.append(l)
 def axis(bounds,panel,title=None):
  x,y,w,h=bounds;ax=fig.add_axes([x/W,y/H,w/W,h/H]);axes.append((panel,ax))
  ax.set_xlim(C['axis_limits']['UMAP1']);ax.set_ylim(C['axis_limits']['UMAP2']);ax.set_aspect('equal',adjustable='box')
  ax.xaxis.set_major_locator(FixedLocator(C['x_ticks']));ax.yaxis.set_major_locator(FixedLocator(C['y_ticks']))
  ax.tick_params(labelsize=C['tick_pt'],length=2.3,width=.5,pad=2.5,direction='out')
  ax.set_xlabel('UMAP 1',fontsize=C['axis_titles_pt'],labelpad=3);ax.set_ylabel('UMAP 2',fontsize=C['axis_titles_pt'],labelpad=3)
  for s in ax.spines.values():s.set_linewidth(.5);s.set_color('#303030');s.set_visible(True)
  ax.grid(False)
  if title:text(x,y+h+5.3,title,C['panel_title_pt'])
  return ax
 def draw(ax,panel,rows,gid,marker,hc,fill,area,alpha,lw,z):
  assert rows
  xy=np.array([[float(r['UMAP1']),float(r['UMAP2'])] for r in rows]);col=shown(hc,mode)
  face='none' if fill=='none' else col
  edge=col if fill=='none' else (shown('#666666',mode) if marker=='o' and lw else ('white' if marker=='^' else 'none'))
  coll=ax.scatter(xy[:,0],xy[:,1],s=area,marker=marker,facecolors=face,edgecolors=edge,linewidths=lw,alpha=alpha,zorder=z,rasterized=False)
  gid=panel+'_'+gid;coll.set_gid(gid);colls.append((panel,coll,rows,hc,fill,marker))
  groups.append(dict(version=version,panel=panel,gid=gid,n_samples=len(rows),marker=marker,HEX=hc,fill=fill,area_pt2=area,alpha=alpha,linewidth_pt=lw,zorder=z))
 def reference_layers(ax,panel,small=False,neutral_sac=False,viral_context=False):
  bg=[r for r in data if r['sample_type']=='reference_BGC' and r['is_saccharide_reference']=='False']
  draw(ax,panel,bg,'nonSac','o',C['non_saccharide_HEX'],'solid',2.3 if small else 3.2,.65,0,1)
  if viral_context:
   vv=[r for r in data if r['sample_type']=='viral_locus'];draw(ax,panel,vv,'viral_context','^',C['context_viral_HEX'],'none',5,.24,.38,2)
  if neutral_sac:
   ss=[r for r in data if r['is_saccharide_reference']=='True'];draw(ax,panel,ss,'Sac_landmarks','o',C['context_saccharide_HEX'],'solid',5,.85,.12,3)
  else:
   for c in C['reference_categories']:
    rr=[r for r in data if r['reference_chemotype_visualization_group']==c['name']]
    draw(ax,panel,rr,'chem'+str(c['order']),'o',c['HEX'],'solid',6.5 if small else 10,.98,.18,3)
 def viral_layers(ax,panel,small=False):
  for c in C['viral_categories']:
   vv=[r for r in data if r['viral_only_leiden_community']==c['name']]
   draw(ax,panel,vv,c['name'],'^',c['HEX'],c['fill'],7.5 if small else 10.5,1,.65 if c['fill']=='none' else .27,5)
 def viral_legend(x,y,spacing=5):
  text(x,y,'Viral Leiden community',C['legend_title_pt'])
  for i,c in enumerate(C['viral_categories']):
   yy=y-5-spacing*i;symbol(x+1.2,yy,'^',c['HEX'],c['fill'],4.5);text(x+4.4,yy,c['name']+' (n = '+str(c['n'])+')',C['legend_pt'])
 if version=='A':
  ax=axis(C['A_axes_mm'],'A');reference_layers(ax,'A');viral_layers(ax,'A')
  text(12,30,'Reference BGC chemotype',C['legend_title_pt'])
  for i,c in enumerate(C['reference_categories']):
   x=12 if i<4 else 75;y=24.5-5.2*(i if i<4 else i-4)
   symbol(x+1.2,y,'o',c['HEX']);text(x+4.3,y,c['legend_label'],C['legend_pt'])
  symbol(76.2,8.9,'o',C['non_saccharide_HEX']);text(79.3,8.9,'Other reference BGCs (2,307)',C['legend_pt'])
  viral_legend(140,30)
 else:
  ax=axis(C['B_axes_mm'][0],'B_left','Reference glycan chemotypes');reference_layers(ax,'B_left',small=True,viral_context=True)
  ax=axis(C['B_axes_mm'][1],'B_right','Viral GT-associated communities');reference_layers(ax,'B_right',small=True,neutral_sac=True);viral_layers(ax,'B_right',small=True)
  text(12,40,'Reference BGC chemotype',C['legend_title_pt'])
  for i,c in enumerate(C['reference_categories']):
   y=34.8-4.7*i;symbol(13.2,y,'o',c['HEX']);text(16.3,y,c['legend_label'],C['legend_pt'])
  symbol(13.2,1.9,'o',C['non_saccharide_HEX']);text(16.3,1.9,'Other reference BGCs (2,307)',C['legend_pt'])
  viral_legend(104,40,4.7)
  symbol(105.2,9,'o',C['context_saccharide_HEX']);text(108.3,9,'Saccharide reference landmarks',C['legend_pt'])
  symbol(105.2,4.3,'^',C['context_viral_HEX'],'none',alpha=.4);text(108.3,4.3,'Viral locus context, left map',C['legend_pt'])
 fig.canvas.draw();renderer=fig.canvas.get_renderer();fw,fh=fig.canvas.get_width_height()
 for panel,coll,rows,hc,fill,marker in colls:
  actual=np.asarray(coll.get_offsets());expected=np.array([[float(r['UMAP1']),float(r['UMAP2'])] for r in rows]);assert np.array_equal(actual,expected)
  display=coll.get_offset_transform().transform(actual)
  for r,xy,pix in zip(rows,actual,display):
   glyphs.append(dict(version=version,panel=panel,sample_id=r['sample_id'],gid=coll.get_gid(),UMAP1=r['UMAP1'],UMAP2=r['UMAP2'],actual_offset_x=float(xy[0]),actual_offset_y=float(xy[1]),svg_x_pt=float(pix[0]*72/fig.dpi),svg_y_pt=float((fig.get_figheight()*fig.dpi-pix[1])*72/fig.dpi),HEX=hc,fill=fill,marker=marker,viral_only_leiden_community=r['viral_only_leiden_community'],reference_chemotype_visualization_group=r['reference_chemotype_visualization_group']))
 extents=[]
 def extent(t,kind):
  b=t.get_window_extent(renderer);extents.append(dict(kind=kind,text=t.get_text(),font_pt=t.get_fontsize(),bbox=[float(b.x0),float(b.y0),float(b.x1),float(b.y1)]))
 for t in texts:extent(t,'legend_or_panel_title')
 for panel,ax in axes:
  for t in [ax.xaxis.label,ax.yaxis.label]:extent(t,'axis_title')
  for t in [*ax.get_xticklabels(),*ax.get_yticklabels()]:extent(t,'tick')
 layout=dict(version=version,width_mm=W,height_mm=H,canvas_pixels=[fw,fh],canvas_dpi=fig.dpi,axes=[dict(panel=p,bounds_mm=[float(x) for x in (ax.get_position().x0*W,ax.get_position().y0*H,ax.get_position().width*W,ax.get_position().height*H)],UMAP1=[float(x) for x in ax.get_xlim()],UMAP2=[float(x) for x in ax.get_ylim()],aspect=float(ax.get_aspect()),x_ticks=[float(x) for x in ax.get_xticks()],y_ticks=[float(x) for x in ax.get_yticks()],visible_spines=all(s.get_visible() for s in ax.spines.values()),grid=False) for p,ax in axes],text=extents)
 return fig,glyphs,groups,layout
def contact_sheet():
 # Each preview remains exactly 183 mm wide at 150 dpi. No shrinking one design relative to the other.
 a=Image.open(OUT/'qc/A_preview_183mm_150dpi.png').convert('RGB');b=Image.open(OUT/'qc/B_preview_183mm_150dpi.png').convert('RGB')
 pad=30;head=60;canvas=Image.new('RGB',(a.width+b.width+pad*3,max(a.height,b.height)+head+pad*2),'white');d=ImageDraw.Draw(canvas);font=ImageFont.truetype(C['font_path'],23)
 d.text((pad,pad),'A  Single integrated UMAP',fill='#202020',font=font);d.text((a.width+pad*2,pad),'B  Paired UMAPs',fill='#202020',font=font)
 canvas.paste(a,(pad,head));canvas.paste(b,(a.width+pad*2,head));canvas.save(OUT/'Figure3C_design_comparison.png',dpi=(150,150))
 save('CONTACT_SHEET_RECEIPT.json',dict(source_A_sha256=sha(OUT/'qc/A_preview_183mm_150dpi.png'),source_B_sha256=sha(OUT/'qc/B_preview_183mm_150dpi.png'),same_physical_width_mm=183,dpi=150,resampling=False,positions_px={'A':[pad,head],'B':[a.width+pad*2,head]},dimensions_px=list(canvas.size)))
def main():
 for n in ['qc']:(OUT/n).mkdir(exist_ok=True)
 data=read('FIGURE_SOURCE_DATA.tsv');glyphs=[];groups=[];layouts=[]
 save('SOFTWARE_VERSIONS.json',dict(python=sys.version,platform=platform.platform(),numpy=np.__version__,matplotlib=matplotlib.__version__,Pillow=Image.__version__,matplotlib_file=matplotlib.__file__,font_resolved=fm.findfont(fm.FontProperties(family='Arial')),font_SHA256=sha(Path(C['font_path'])),backend=matplotlib.get_backend(),PDF_fonttype=42,SVG_fonttype='none',no_model_execution=True))
 for version,stem in [('A','Figure3C_A_single_integrated_umap'),('B','Figure3C_B_paired_umaps')]:
  fig,g,gp,layout=build(data,version);glyphs.extend(g);groups.extend(gp);layouts.append(layout)
  fig.savefig(OUT/(stem+'.pdf'),metadata={'Title':'Figure 3C alternative '+version,'Subject':'Frozen joint UMAP; reference chemotypes and viral-only Leiden communities; exploratory design'})
  fig.savefig(OUT/(stem+'.svg'),metadata={'Title':'Figure 3C alternative '+version})
  fig.savefig(OUT/(stem+'.png'),dpi=C['dpi'])
  for dpi in [150,96]:fig.savefig(OUT/('qc/'+version+'_preview_183mm_'+str(dpi)+'dpi.png'),dpi=dpi)
  plt.close(fig)
  for mode in [*MATRICES,'grayscale']:
   fig,_,_,_=build(data,version,mode);fig.savefig(OUT/('qc/'+version+'_'+mode+'_preview.png'),dpi=150);plt.close(fig)
 write('PLOTTED_SAMPLE_RECEIPTS.tsv',glyphs);write('RENDER_GROUP_RECEIPTS.tsv',groups);save('LAYOUT_RECEIPTS.json',layouts)
 diagnostics=[];pal=read('FINAL_PALETTE.tsv')
 for mode in ['normal',*MATRICES,'grayscale']:
  for x,y in itertools.combinations([r for r in pal if r['section']!='Context'],2):
   diagnostics.append(dict(mode=mode,category_1=x['category'],category_2=y['category'],same_section=str(x['section']==y['section']),CIE76_deltaE=float(np.linalg.norm(lab(shown(x['HEX'],mode))-lab(shown(y['HEX'],mode)))),redundant_cue='different shape' if x['marker']!=y['marker'] else ('filled versus open' if x['fill']!=y['fill'] else 'color and position only')))
 write('PALETTE_CVD_DIAGNOSTICS.tsv',diagnostics);save('CVD_SIMULATION_PARAMETERS.json',dict(method='Machado severity-100 matrices applied in linear sRGB, gamut-clipped; grayscale uses D65 relative luminance',matrices=MATRICES,certifies_accessibility=False))
 contact_sheet()
 save('PLOT_PROVENANCE.json',dict(source_data_SHA256=sha(OUT/'FIGURE_SOURCE_DATA.tsv'),original_projection_samples=3058,original_references=2502,original_viral_loci=556,displayed_unique_samples=2789,reference_BGCs=2502,Saccharide_references=195,selected_viral_loci=287,omitted_unselected_viral_loci=269,A_maps=1,B_maps=2,each_map_samples=2789,viral_Leiden_space_dimensions=320,viral_Leiden_baseline='k20_r1_s42',no_new_UMAP=True,no_new_clustering=True,no_new_inference=True,no_new_chemical_annotation=True,no_new_nearest_neighbors=True,no_jitter=True,axes_equal_in_all_maps=True,identical_underlying_data=True,committed=False,pushed=False))
 print(json.dumps(dict(status='GENERATED',unique_samples=len(data),rendered_glyphs=len(glyphs),A_dimensions_mm=[183,150],B_dimensions_mm=[183,111],render_groups=len(groups)),indent=2))
if __name__=='__main__':main()
