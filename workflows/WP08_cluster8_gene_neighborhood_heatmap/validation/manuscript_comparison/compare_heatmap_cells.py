"""Decode all cells using solid interior pixels and independently located image grids.
No Foldseek, original builder, or plotting script is executed.
"""
import ast,csv,hashlib,json,math,shutil
from collections import Counter
from pathlib import Path
import numpy as np
from PIL import Image,ImageDraw,ImageFont
BASE=Path(__file__).resolve().parent
OUT=BASE/'WP08_manuscript_comparison_20261102'
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 matrix=OUT/'historical/cluster8_neighborhood_matrix.tsv'
 with matrix.open(newline='') as f: rows=list(csv.DictReader(f,delimiter='\t'))
 # Read class column order as a literal, without importing or executing historical code.
 code=BASE/'WP08_recovery_staging/originals/260630_Temp/plot_cluster8_panelB_heatmap_clean.py'
 tree=ast.parse(code.read_text())
 modules=next(ast.literal_eval(n.value) for n in tree.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='MODULES' for t in n.targets))
 classes=[c for cols in modules.values() for c in cols]
 assert len(rows)==11 and len(classes)==33
 groups={'C5-MT + 5mYOX-associated':1,'dCMP hydroxymethylase-associated':2,'ThyX/ThyS / nucleotide-metabolism-associated':3,'additional-GT/glycan-associated':4,'ambiguous/no DNA-context':5}
 rows.sort(key=lambda r:(groups[r['primary_group']],-int(r['DNA_context_score']),-int(r['glycan_score']),-int(r['virion_score']),r['candidate_id']))
 files={'historical':OUT/'historical/PanelB_cluster8_global_neighborhood_heatmap_clean.png','manuscript':OUT/'manuscript_images/Viral_GT_draftv1.4_261022_LanNa_references_added/word/media/image19.png'}
 # Grid bounds from native-resolution image runs; right/bottom are exclusive.
 bounds={'historical':(286,364,3806,1452),'manuscript':(114,145,1512,578)}
 arrays={n:np.asarray(Image.open(p).convert('RGB')) for n,p in files.items()}
 palette=np.array([[242,242,242],[198,219,239],[107,174,214],[8,81,156]])
 palette_checks={}; measured={}; evidence=[]
 for name,a in arrays.items():
  x0,y0,x1,y1=bounds[name]; w=(x1-x0)/33;h=(y1-y0)/11
  # Validate the palette independently against the colorbar from bottom=0 to top=3.
  color_x=3910 if name=='historical' else 1553
  bar=[np.median(a[int(y1-(s+.5)*(y1-y0)/4)-3:int(y1-(s+.5)*(y1-y0)/4)+4,color_x-3:color_x+4].reshape(-1,3),axis=0).astype(int).tolist() for s in range(4)]
  palette_checks[name]=dict(sampled_colorbar_RGB=bar,expected_RGB=palette.tolist(),exact_match=bool(np.array_equal(np.array(bar),palette)))
  assert palette_checks[name]['exact_match']
  measured[name]=[]
  for ri,row in enumerate(rows):
   for ci,cls in enumerate(classes):
    cx=x0+(ci+.5)*w;cy=y0+(ri+.5)*h
    # A 40%-width/height interior patch avoids gridlines; test center offsets too.
    patches=[]
    for dx,dy in [(0,0),(-.12,0),(.12,0),(0,-.12),(0,.12)]:
     xa=round(cx+(dx-.18)*w);xb=round(cx+(dx+.18)*w)
     ya=round(cy+(dy-.18)*h);yb=round(cy+(dy+.18)*h)
     pixels=a[ya:yb,xa:xb].reshape(-1,3).astype(float)
     distances=np.sqrt(((pixels[:,None,:]-palette[None,:,:])**2).sum(axis=2))
     nearest=distances.argmin(axis=1)
     counts=np.bincount(nearest,minlength=4); score=int(counts.argmax())
     patches.append(dict(score=score,purity=float(counts[score]/len(pixels)),median_RGB=np.median(pixels,axis=0).tolist(),median_distance=float(np.median(distances[:,score]))))
    certain=all(p['score']==patches[0]['score'] and p['purity']>=.99 and p['median_distance']<=3 for p in patches)
    rec=dict(locus=row['candidate_id'],gene_class=cls,TSV_score=int(row[cls]),image_score=patches[0]['score'] if certain else None,
       status='confirmed' if certain else 'uncertain',center_x=round(cx,3),center_y=round(cy,3),min_patch_purity=min(p['purity'] for p in patches),median_RGB=patches[0]['median_RGB'])
    measured[name].append(rec)
 # Separate grouping strip from score grid. Read central strip color per row.
 group_palette={'dCMP hydroxymethylase-associated':(31,120,180),'ThyX/ThyS / nucleotide-metabolism-associated':(51,160,44),'additional-GT/glycan-associated':(255,127,0)}
 grouping=[]
 for ri,row in enumerate(rows):
  result=dict(locus=row['candidate_id'],TSV_group=row['primary_group'])
  for name,a in arrays.items():
   x0,y0,x1,y1=bounds[name];ya=round(y0+(ri+.12)*(y1-y0)/11);yb=round(y0+(ri+.88)*(y1-y0)/11);x=195 if name=='historical' else 78
   patch=a[ya:yb,x-10:x+11].reshape(-1,3)
   # Text overlaps the historical group strip. Its modal fill, away from the
   # horizontal separators, identifies the category without mistaking black ink.
   rgb,n=Counter(map(tuple,patch)).most_common(1)[0]
   result[name+'_modal_fill_fraction']=n/len(patch)
   result[name+'_RGB']=[int(n) for n in rgb]
   result[name+'_group']=next((g for g,c in group_palette.items() if rgb==c),'uncertain')
  grouping.append(result)
 combined=[]
 for h,m in zip(measured['historical'],measured['manuscript']):
  combined.append(dict(locus=m['locus'],gene_class=m['gene_class'],original_TSV_score=m['TSV_score'],historical_PNG_score=h['image_score'],manuscript_figure_score=m['image_score'],comparison='uncertain' if m['image_score'] is None or h['image_score'] is None else 'agree' if m['TSV_score']==m['image_score']==h['image_score'] else 'discrepancy',manuscript_patch_purity=m['min_patch_purity'],manuscript_RGB=m['median_RGB'],manuscript_center_x=m['center_x'],manuscript_center_y=m['center_y']))
 def table(p,recs):
  with p.open('w',newline='') as f:
   w=csv.DictWriter(f,fieldnames=list(recs[0]));w.writeheader();w.writerows(recs)
 table(OUT/'all_363_score_cell_comparisons.csv',combined)
 table(OUT/'row_group_comparisons.csv',grouping)
 (OUT/'image_sampling_details.json').write_text(json.dumps(dict(bounds=bounds,palette_checks=palette_checks,measurements=measured),indent=2))
 # Preserve a figure-derived matrix independently from the historical data.
 with (OUT/'manuscript_figure_decoded_matrix.tsv').open('w',newline='') as f:
  writer=csv.DictWriter(f,fieldnames=['candidate_id']+classes+['primary_group'],delimiter='\t');writer.writeheader()
  for ri,row in enumerate(rows):
   decoded={c:combined[ri*33+ci]['manuscript_figure_score'] for ci,c in enumerate(classes)}
   writer.writerow(dict(candidate_id=row['candidate_id'],**decoded,primary_group=grouping[ri]['manuscript_group']))
 # Compare common downsampling algorithms; mismatch alone cannot establish editing.
 original=Image.open(files['historical']).convert('RGB');target=arrays['manuscript']
 resize=[]
 for name in ['NEAREST','BOX','BILINEAR','HAMMING','BICUBIC','LANCZOS']:
  candidate=np.asarray(original.resize((target.shape[1],target.shape[0]),resample=getattr(Image.Resampling,name)))
  diff=np.abs(candidate.astype(float)-target.astype(float))
  resize.append(dict(method=name,pixel_equal=bool(np.array_equal(candidate,target)),mean_absolute_channel_error=float(diff.mean()),fraction_pixels_identical=float(np.all(diff==0,axis=2).mean())))
 best=min(resize,key=lambda r:r['mean_absolute_channel_error'])
 resized=original.resize((target.shape[1],target.shape[0]),resample=getattr(Image.Resampling,best['method']))
 resized.save(OUT/'comparison_resized_historical.png')
 # Enlarged border-free score grid: top original TSV, middle historical pixels, bottom manuscript pixels.
 sheet=Image.new('RGB',(33*28+150,11*28*3+90),'white');draw=ImageDraw.Draw(sheet)
 for panel,name in enumerate(['TSV','historical PNG','manuscript PNG']):
  ybase=panel*(11*28+30);draw.text((2,ybase),name,fill='black')
  for ri,row in enumerate(rows):
   draw.text((2,ybase+25+ri*28),row['candidate_id'],fill='black')
   for ci,cls in enumerate(classes):
    score=int(row[cls]) if panel==0 else measured['historical' if panel==1 else 'manuscript'][ri*33+ci]['image_score']
    color=tuple(palette[score]) if score is not None else (255,0,255)
    xx=150+ci*28;yy=ybase+22+ri*28;draw.rectangle((xx,yy,xx+26,yy+26),fill=color)
 sheet.save(OUT/'score_grid_comparison.png')
 summary=dict(cells=363,agreements=sum(r['comparison']=='agree' for r in combined),
  confirmed_discrepancies=[r for r in combined if r['comparison']=='discrepancy'],uncertain=[r for r in combined if r['comparison']=='uncertain'],
  grouping_agrees=all(r['TSV_group']==r['historical_group']==r['manuscript_group'] for r in grouping),
  resize_comparisons=resize,best_resize=best,
  preservation=dict(historical_matrix_sha256=sha(matrix),historical_png_sha256=sha(files['historical']),manuscript_png_sha256=sha(files['manuscript'])),
  method='Each cell decoded from five border-excluding interior patches and calibrated to the image colorbar; certainty requires >=99% nearest-color agreement in every patch, median RGB distance <=3 and the same decoded score across offsets. No scoring code or Foldseek executed.')
 (OUT/'comparison_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
 shutil.copy2(__file__,OUT/Path(__file__).name)
 print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
