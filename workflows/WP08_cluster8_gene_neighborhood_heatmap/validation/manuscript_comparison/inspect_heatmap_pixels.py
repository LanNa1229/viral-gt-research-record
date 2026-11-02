import io,json,zipfile
from pathlib import Path
from collections import Counter
from PIL import Image,ImageDraw
import numpy as np
BASE=Path(__file__).resolve().parent
OUT=BASE/'WP08_manuscript_comparison_20261102'
records=json.loads((OUT/'document_inventory.json').read_text())
for ri in [0,2]:
 r=records[ri]; ims=[]
 with zipfile.ZipFile(r['path']) as z:
  for entry in r['images']:
   if entry.get('format')!='PNG': continue
   im=Image.open(io.BytesIO(z.read(entry['member']))).convert('RGB'); im.thumbnail((280,160))
   ims.append((entry['member'],im))
 sheet=Image.new('RGB',(5*300,((len(ims)+4)//5)*190),'white'); draw=ImageDraw.Draw(sheet)
 for i,(label,im) in enumerate(ims):
  x=(i%5)*300;y=(i//5)*190;draw.text((x+5,y+2),label,fill='black');sheet.paste(im,(x+5,y+22))
 sheet.save(OUT/f'contact_sheet_{ri}.png')
paths=[OUT/'historical/PanelB_cluster8_global_neighborhood_heatmap_clean.png',OUT/'manuscript_images/Viral_GT_draftv1.4_261022_LanNa_references_added/word/media/image19.png']
for p in paths:
 im=Image.open(p).convert('RGB'); a=np.asarray(im)
 print(p.name,im.size,Counter(map(tuple,a.reshape(-1,3))).most_common(12))
 # Blue pixels: locate horizontal heatmap body independently of labels or colorbar.
 mask=(a[:,:,2]>a[:,:,0]+20)&(a[:,:,1]<180)&(a[:,:,0]<160)
 counts=mask.sum(axis=1); ys=np.where(counts>im.width*.10)[0]
 print('blue rows',ys.min(),ys.max())
 for y in [int(ys.min()+12),int((ys.min()+ys.max())/2)]:
  row=mask[y]; edges=np.flatnonzero(np.diff(np.pad(row.astype(int),(1,1))))
  print('y',y,'segments',list(zip(edges[::2].tolist(),edges[1::2].tolist())))
