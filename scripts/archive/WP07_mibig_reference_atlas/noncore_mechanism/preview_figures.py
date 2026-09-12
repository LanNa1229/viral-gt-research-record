from common import *
from PIL import Image
import io,base64
ims=[]
for name in sys.argv[1:]:
 im=Image.open(ROOT/f"figures/{name}.png").convert("RGB");im.thumbnail((1500,650));ims.append(im)
out=Image.new("RGB",(1500,sum(im.height+15 for im in ims)),"white");y=0
for im in ims:out.paste(im,((1500-im.width)//2,y));y+=im.height+15
b=io.BytesIO();out.save(b,format="PNG");print(base64.b64encode(b.getvalue()).decode())
