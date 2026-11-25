"""One 100% stacked horizontal bar chart from verified frozen source tables."""
import csv,hashlib,itertools,json,os,sys
from pathlib import Path
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[1]
os.environ['MPLCONFIGDIR']=str(OUT/'runtime/mplconfig');os.environ['MPLBACKEND']='Agg'
sys.path.append(str(ROOT/'analysis/bgcprophet_reference/.deps'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.patches import Patch
from matplotlib.ticker import FixedLocator,PercentFormatter
CONFIG=json.loads((OUT/'config.json').read_text())
fm.fontManager.addfont(CONFIG['font_path'])
plt.rcParams.update({'font.family':'Arial','font.size':8,'text.color':'black','axes.labelcolor':'black','axes.edgecolor':'black','axes.linewidth':.5,'axes.grid':False,'figure.facecolor':'white','savefig.facecolor':'white','pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none'})
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
def shown(h,mode):
 a=rgb(h)
 if mode in MATRICES:return enc(lin(a)@np.asarray(MATRICES[mode]).T)
 if mode=='grayscale':
  y=float(lin(a)@np.array([.2126729,.7151522,.0721750]));return enc(np.repeat(y,3))
 return a
def lab(a):
 xyz=lin(a)@np.array([[.4124564,.3575761,.1804375],[.2126729,.7151522,.0721750],[.0193339,.1191920,.9503041]]).T
 xyz=xyz/np.array([.95047,1,1.08883]);d=6/29
 f=np.where(xyz>d**3,np.cbrt(xyz),xyz/(3*d*d)+4/29)
 return np.array([116*f[1]-16,500*(f[0]-f[1]),200*(f[1]-f[2])])
def build(rows,mode='normal'):
 W=CONFIG['width_mm'];H=CONFIG['height_mm'];fig=plt.figure(figsize=(W/25.4,H/25.4))
 left,bottom,width,height=CONFIG['axes_bounds_mm'];ax=fig.add_axes([left/W,bottom/H,width/W,height/H])
 yy={g:2-i for i,g in enumerate(CONFIG['populations'])}
 patches=[];receipts=[]
 for r in rows:
  c=CONFIG['categories'][int(r['category_order'])-1];p=ax.barh(yy[r['population']],float(r['percentage']),left=float(r['left_percentage']),height=CONFIG['bar_height'],color=shown(c['HEX'],mode),edgecolor='white',linewidth=CONFIG['separator_linewidth_pt'],zorder=3)[0]
  p.set_gid('segment_P'+r['population_order']+'_C'+r['category_order'])
  patches.append(p);receipts.append(dict(population=r['population'],category_order=r['category_order'],chemotype_group=r['chemotype_group'],n_BGCs=r['n_BGCs'],denominator=r['denominator_saccharide_BGCs'],left_percentage=p.get_x(),width_percentage=p.get_width(),right_percentage=p.get_x()+p.get_width(),y=p.get_y(),height=p.get_height(),gid=p.get_gid(),HEX=c['HEX']))
 sizes={r['population']:int(r['denominator_saccharide_BGCs']) for r in rows}
 ax.set_xlim(CONFIG['x_limits']);ax.set_ylim(-.6,2.6)
 ax.xaxis.set_major_locator(FixedLocator(CONFIG['x_ticks']));ax.xaxis.set_major_formatter(PercentFormatter(xmax=100,decimals=0))
 ax.set_yticks([2,1,0],[g+' (n = '+str(sizes[g])+')' for g in CONFIG['populations']],fontsize=CONFIG['bar_label_pt'])
 ax.set_xlabel(CONFIG['x_title'],fontsize=CONFIG['axis_title_pt'],labelpad=5)
 ax.tick_params(axis='x',labelsize=CONFIG['tick_pt'],length=3,width=.5,pad=3)
 ax.tick_params(axis='y',length=0,pad=7)
 ax.spines[['top','right','left']].set_visible(False);ax.spines['bottom'].set_linewidth(.5);ax.grid(False)
 handles=[Patch(facecolor=shown(c['HEX'],mode),edgecolor='none',label=c['short_label']) for c in CONFIG['categories']]
 leg=fig.legend(handles=handles,ncol=2,loc='lower center',bbox_to_anchor=(.5,.035),frameon=False,fontsize=CONFIG['legend_pt'],handlelength=1.15,handleheight=.85,handletextpad=.65,columnspacing=2.1,labelspacing=.7,borderaxespad=0,borderpad=0)
 fig.canvas.draw();renderer=fig.canvas.get_renderer()
 texts=[]
 for kind,items in [('axis_title',[ax.xaxis.label]),('x_tick',ax.get_xticklabels()),('bar_label',ax.get_yticklabels()),('legend',leg.get_texts())]:
  for t in items:
   b=t.get_window_extent(renderer);texts.append(dict(kind=kind,text=t.get_text(),font_pt=t.get_fontsize(),x0=float(b.x0),x1=float(b.x1),y0=float(b.y0),y1=float(b.y1)))
 actual=ax.get_position();fw,fh=fig.canvas.get_width_height()
 layout=dict(dimensions_mm=[W,H],canvas_pixels=[fw,fh],canvas_dpi=fig.dpi,axes_bounds=list(actual.bounds),axes_bounds_mm=[actual.x0*W,actual.y0*H,actual.width*W,actual.height*H],x_limits=list(ax.get_xlim()),x_ticks=[float(x) for x in ax.get_xticks()],y_order_top_to_bottom=CONFIG['populations'],text=texts,scientific_axes=len(fig.axes),in_bar_text_count=len(ax.texts),legend_columns=2,zero_width_segments=sum(float(r['width_percentage'])==0 for r in receipts),minimum_positive_segment_percent=min(float(r['width_percentage']) for r in receipts if float(r['width_percentage'])>0))
 return fig,receipts,layout
def main():
 status=json.loads((OUT/'INPUT_VALIDATION.json').read_text());assert status['status']=='PASS' and status['counts_match_brief']
 for n in ['figures','qc']:(OUT/n).mkdir(exist_ok=True)
 rows=read('FIGURE_SOURCE_DATA.tsv');assert len(rows)==21
 save('SOFTWARE_VERSIONS.json',dict(python=sys.version,numpy=np.__version__,matplotlib=matplotlib.__version__,matplotlib_file=matplotlib.__file__,backend=matplotlib.get_backend(),resolved_font=fm.findfont(fm.FontProperties(family='Arial')),PDF_fonttype=42,SVG_fonttype='none',no_models_or_annotation=True))
 fig,receipts,layout=build(rows)
 stem=OUT/'figures/SuppFig2_saccharide_reference_chemotype_distribution'
 fig.savefig(stem.with_suffix('.pdf'),metadata={'Title':'Chemotype distribution of Saccharide reference BGCs','Subject':'Descriptive display of frozen MIBiG 3.1 proposed annotations; no new analysis'})
 fig.savefig(stem.with_suffix('.svg'),metadata={'Title':'Chemotype distribution of Saccharide reference BGCs'})
 fig.savefig(stem.with_suffix('.png'),dpi=CONFIG['dpi'])
 fig.savefig(OUT/'qc/preview_183mm_150dpi.png',dpi=150);fig.savefig(OUT/'qc/preview_183mm_96dpi.png',dpi=96);plt.close(fig)
 write('PLOTTED_SEGMENT_RECEIPTS.tsv',receipts);save('LAYOUT_RECEIPTS.json',layout)
 diagnostics=[]
 for mode in ['normal',*MATRICES,'grayscale']:
  if mode!='normal':
   f,_,_=build(rows,mode);f.savefig(OUT/('qc/'+mode+'_preview.png'),dpi=150);plt.close(f)
  for a,b in itertools.combinations(CONFIG['categories'],2):
   difference=float(np.linalg.norm(lab(shown(a['HEX'],mode))-lab(shown(b['HEX'],mode))))
   diagnostics.append(dict(view=mode,category_a=a['name'],category_b=b['name'],DeltaE76=difference,scope='Color appearance diagnostic only; not universal accessibility threshold'))
 write('PALETTE_DISTINGUISHABILITY.tsv',diagnostics)
 save('PLOT_PROVENANCE.json',dict(input_source_data_sha256=hashlib.sha256((OUT/'FIGURE_SOURCE_DATA.tsv').read_bytes()).hexdigest(),n_populations=3,n_categories=7,n_source_BGCs=195,bar_denominators=status['population_sizes'],bar_segments_including_structural_zeros=21,all_counts_reconstructed=True,expected_values_used_only_for_validation=True,not_biological_reannotation=True,no_new_clustering=True,no_enrichment_tests=True,reference_version='MIBiG 3.1',color_checks='Machado severity-100 matrices in linear sRGB plus grayscale; exact matrices in script',surface_definition='Categories 1-4 only; not the broader phosphonoglycan-inclusive flag'))
 print(json.dumps(dict(status='GENERATED',segments=len(receipts),minimum_positive_segment_percent=layout['minimum_positive_segment_percent'],dimensions_mm=layout['dimensions_mm']),indent=2),flush=True)
if __name__=='__main__':main()
