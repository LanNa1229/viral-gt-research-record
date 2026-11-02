"""Read-only DOCX media extraction and targeted Cluster 8 figure audit."""
import csv, hashlib, io, json, shutil, zipfile
from pathlib import Path
from xml.etree import ElementTree as ET
from PIL import Image

BASE=Path(__file__).resolve().parent
OUT=BASE/'WP08_manuscript_comparison_20261102'
KNOWN='8a938331cff92bd986ec142b63c18cf316ad6557cf6f2172701ee9178764e5bb'
DOCS=[
 Path('D:/Projects/Viral_GT_Section1_Structural_Clustering_WorkFlow/Viral_GT_draftv1.3_261019_LanNa_references_added_2.docx'),
 Path('D:/Projects/Viral_GT_Section1_Structural_Clustering_WorkFlow/Viral_GT_draftv1.4_261022_LanNa_references_added.docx'),
 Path('C:/Users/rahnn/Downloads/Viral_GT_draftv1.5_261027_SuppleFig_embedding_analysis.docx'),
 Path('C:/Users/rahnn/Downloads/Viral_GT_draftv1.5_261027_LanNa_embedding_supple_final.docx')]
def sha_bytes(b): return hashlib.sha256(b).hexdigest()
def sha(p): return sha_bytes(p.read_bytes())
def main():
 OUT.mkdir(exist_ok=False)
 src=BASE/'WP08_recovery_staging/originals/260630_Temp'
 originals=OUT/'historical'; originals.mkdir()
 for p in [src/'cluster8_global_neighborhood_panelB/cluster8_neighborhood_matrix.tsv',src/'cluster8_global_neighborhood_panelB/figures_clean/PanelB_cluster8_global_neighborhood_heatmap_clean.png']:
  shutil.copy2(p,originals/p.name)
  assert sha(p)==sha(originals/p.name)
 records=[]
 ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main','a':'http://schemas.openxmlformats.org/drawingml/2006/main','r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships'}
 for doc in DOCS:
  st=doc.stat(); before=sha(doc); record=dict(path=str(doc),sha256=before,bytes=st.st_size,mtime_ns=st.st_mtime_ns,images=[],drawings=[])
  with zipfile.ZipFile(doc) as z:
   root=ET.fromstring(z.read('word/document.xml'))
   relroot=ET.fromstring(z.read('word/_rels/document.xml.rels'))
   rels={r.attrib['Id']:r.attrib.get('Target','') for r in relroot}
   paragraphs=root.findall('.//w:p',ns)
   for i,p in enumerate(paragraphs):
    for blip in p.findall('.//a:blip',ns):
     rid=blip.get('{'+ns['r']+'}embed'); target=rels.get(rid,'')
     record['drawings'].append(dict(paragraph=i,target=target,xml=ET.tostring(p,encoding='unicode'),
       surrounding_text=[''.join(t.text or '' for t in q.findall('.//w:t',ns)) for q in paragraphs[max(0,i-2):i+5]]))
   record['core_properties']=z.read('docProps/core.xml').decode() if 'docProps/core.xml' in z.namelist() else None
   for n in z.namelist():
    if not n.startswith('word/media/'): continue
    b=z.read(n); entry=dict(member=n,sha256=sha_bytes(b),bytes=len(b))
    try:
     im=Image.open(io.BytesIO(b)); entry.update(size=list(im.size),format=im.format)
    except Exception: pass
    if entry['sha256']==KNOWN or entry.get('size')==[1646,779]:
     dest=OUT/'manuscript_images'/doc.stem/n
     dest.parent.mkdir(parents=True,exist_ok=True); dest.write_bytes(b)
     entry['preserved_path']=dest.relative_to(OUT).as_posix()
    record['images'].append(entry)
   text=[''.join(t.text or '' for t in p.findall('.//w:t',ns)) for p in paragraphs]
   record['cluster8_relevant_text']=[dict(paragraph=i,text=p) for i,p in enumerate(text) if any(s in p.lower() for s in ['cluster 8','cluster8','tet-like','neighborhood heatmap','5myox'])]
  record['source_unchanged']=sha(doc)==before and doc.stat().st_mtime_ns==st.st_mtime_ns
  assert record['source_unchanged']; records.append(record)
 (OUT/'document_inventory.json').write_text(json.dumps(records,indent=2)+'\n')
 (OUT/'historical_provenance.json').write_text(json.dumps([dict(original_source='/home/rahnn/260630_Temp/'+p.relative_to(src).as_posix(),preserved_source=str(p),sha256=sha(p)) for p in [src/'cluster8_global_neighborhood_panelB/cluster8_neighborhood_matrix.tsv',src/'cluster8_global_neighborhood_panelB/figures_clean/PanelB_cluster8_global_neighborhood_heatmap_clean.png']],indent=2))
 shutil.copy2(__file__,OUT/Path(__file__).name)
 print(json.dumps([dict(document=r['path'],matches=[i for i in r['images'] if i['sha256']==KNOWN or i.get('size')==[1646,779]]) for r in records],indent=2))
if __name__=='__main__': main()
