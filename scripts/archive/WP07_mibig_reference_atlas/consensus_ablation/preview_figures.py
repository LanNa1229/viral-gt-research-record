from common import *
from PIL import Image,ImageOps
import io,base64
images=[]
for name in sys.argv[1:]:
 im=Image.open(ROOT/"figures"/(name+".png")).convert("RGB")
 im.thumbnail((1350,650));images.append(im)
canvas=Image.new("RGB",(1350,sum(im.height+12 for im in images)),"white")
y=0
for im in images:canvas.paste(im,((1350-im.width)//2,y));y+=im.height+12
b=io.BytesIO();canvas.save(b,format="PNG")
print(base64.b64encode(b.getvalue()).decode())
