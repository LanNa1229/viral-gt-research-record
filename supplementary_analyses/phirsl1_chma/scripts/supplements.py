from pipeline import *
import re
s=(ROOT/'sources/PMC10299810.xml').read_text()
for path in sorted(set(re.findall(r'<\?cloudpmc-path ([^?]+)\?>',s))):
 if path.endswith(('.pdf','.xlsx','.csv')):
  b=fetch('https://cdn.ncbi.nlm.nih.gov/pmc/blobs/'+path,'sources/PMC10299810/'+path.split('/')[-1],True)
  print(path.split('/')[-1],len(b) if b else None,flush=True)
