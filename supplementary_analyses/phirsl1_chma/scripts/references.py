#!/usr/bin/env python3
"""Acquire native ChmA references and published conserved-gene supplement."""
from pipeline import *
import xml.etree.ElementTree as ET,zipfile
def main():
 rows=[];seqs=[]
 for acc in ['B3FIW8','A0A482GDX1','F8SJT5']:
  j=json.loads(fetch(f'https://rest.uniprot.org/uniprotkb/{acc}.json',f'sources/uniprot_{acc}.json'))
  seq=j['sequence']['value'];seqs.append(f'>{acc}\n{seq}\n')
  rows.append(dict(accession=acc,organism=j['organism']['scientificName'],length=len(seq),genes=json.dumps(j.get('genes')),refseq=','.join(x['id'] for x in j['uniProtKBCrossReferences'] if x['database']=='RefSeq')))
 # Wider sequence set, annotated as chimallin; retain evidence distinction from characterized references.
 query='protein_name:chimallin'
 fetch('https://rest.uniprot.org/uniprotkb/search?'+urlencode(dict(query=query,format='fasta',size=500)),'inputs/chma_uniprot_expanded.faa')
 for name,acc in [('phiKZ','NC_004629.1'),('RAY','KU886224.1')]:
  efetch('nuccore',acc,'gbwithparts',f'inputs/{acc}.gb'); rec=parse_proteome(f'inputs/{acc}.gb',name)
  wanted=54 if name=='phiKZ' else 222
  found=[]
  for f in rec.features:
   if f.type!='CDS':continue
   q=f.qualifiers;labels=' '.join(q.get('gene',[])+q.get('locus_tag',[])+q.get('product',[]))
   if ((name=='phiKZ' and q.get('locus_tag')==['PHIKZ054']) or (name=='RAY' and q.get('locus_tag')==['RAY_222'])):
    found.append(q)
  print(name,'ChmA label candidates',[(x.get('protein_id'),x.get('gene'),x.get('locus_tag'),x.get('product')) for x in found],flush=True)
  for q in found:
   pid=q['protein_id'][0];seq=q['translation'][0];seqs.append(f'>{pid}\n{seq}\n');rows.append(dict(accession=pid,organism=name,length=len(seq),genes=json.dumps(q.get('gene')),refseq=pid))
 (ROOT/'inputs/chma_characterized.faa').write_text(''.join(seqs));table('results/characterized_references.tsv',rows)
 for pmc in ['PMC10299810','PMC9365700']:
  root=ET.parse(ROOT/f'sources/{pmc}.xml').getroot()
  links=[]
  for el in root.iter():
   h=el.attrib.get('{http://www.w3.org/1999/xlink}href','')
   if re.search(r'\.(?:csv|xlsx|zip|txt|fasta|fa|pdf)$',h) and 'supp' in h:links.append(h)
  print(pmc,'supplement links',sorted(set(links)),flush=True)
 # Europe PMC supplementary endpoint supplies original files in a ZIP.
 for pmc in ['PMC10299810','PMC9365700']:
  b=fetch(f'https://www.ebi.ac.uk/europepmc/webservices/rest/{pmc}/supplementaryFiles',f'sources/{pmc}_supplements.zip',True)
  if b and zipfile.is_zipfile(io.BytesIO(b)):
   with zipfile.ZipFile(io.BytesIO(b)) as z:
    for n in z.namelist():
     if n.lower().endswith(('.csv','.xlsx','.pdf','.txt','.fasta','.fa','.docx')):
      dest=ROOT/'sources'/pmc/Path(n).name;dest.parent.mkdir(exist_ok=True);dest.write_bytes(z.read(n));print('Extracted',dest.name,flush=True)
if __name__=='__main__':main()
