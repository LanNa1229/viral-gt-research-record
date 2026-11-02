"""Test image-coordinate conventions without editing either source image."""
import json
from pathlib import Path
import numpy as np
from PIL import Image
B=Path(__file__).resolve().parent/'WP08_manuscript_comparison_20261102'
a=np.asarray(Image.open(B/'historical/PanelB_cluster8_global_neighborhood_heatmap_clean.png').convert('RGB'))
b=np.asarray(Image.open(B/'manuscript_images/Viral_GT_draftv1.4_261022_LanNa_references_added/word/media/image19.png').convert('RGB'))
sh,sw=a.shape[:2];th,tw=b.shape[:2]
results=[]
for convention in ['edge_floor','center_floor','endpoints_round']:
 for ox in [-.5,-.25,0,.25,.5]:
  for oy in [-.5,-.25,0,.25,.5]:
   if convention=='edge_floor':
    xx=np.floor(np.arange(tw)*sw/tw+ox);yy=np.floor(np.arange(th)*sh/th+oy)
   elif convention=='center_floor':
    xx=np.floor((np.arange(tw)+.5)*sw/tw+ox);yy=np.floor((np.arange(th)+.5)*sh/th+oy)
   else:
    xx=np.round(np.arange(tw)*(sw-1)/(tw-1)+ox);yy=np.round(np.arange(th)*(sh-1)/(th-1)+oy)
   c=a[np.clip(yy.astype(int),0,sh-1)[:,None],np.clip(xx.astype(int),0,sw-1)[None,:]]
   exact=float(np.all(c==b,axis=2).mean());mae=float(np.abs(c.astype(float)-b.astype(float)).mean())
   results.append(dict(convention=convention,x_offset=ox,y_offset=oy,fraction_pixels_identical=exact,mean_absolute_channel_error=mae,pixel_equal=bool(np.array_equal(c,b))))
results.sort(key=lambda r:r['mean_absolute_channel_error'])
(B/'resampling_coordinate_tests.json').write_text(json.dumps(results,indent=2))
print(json.dumps(results[:5],indent=2))
