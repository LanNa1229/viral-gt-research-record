#!/usr/bin/env python3
"""phiRSL1 reproducible ChmA analysis. No uploads."""
from pathlib import Path
import csv,hashlib,json,subprocess,time,sys,io,re
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime,timezone
from urllib.request import Request,urlopen
from urllib.parse import urlencode
from Bio import SeqIO
ROOT=Path(__file__).resolve().parents[1]
for d in ['inputs','inputs/pdb','inputs/afdb','results','logs','work','sources']:
 (ROOT/d).mkdir(parents=True,exist_ok=True)
def stamp(): return datetime.now(timezone.utc).isoformat()
def digest(b): return hashlib.sha256(b).hexdigest()
def table(path,rows,fields=None):
 rows=list(rows); fields=fields or (list(rows[0]) if rows else [])
 with (ROOT/path).open('w') as f:
  w=csv.DictWriter(f,fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)
def fetch(url,dest,optional=False):
 dest=ROOT/dest
 if dest.exists(): return dest.read_bytes()
 err=None
 for attempt in range(3):
  try:
   with urlopen(Request(url,headers={'User-Agent':'phiRSL1-ChmA-research/1.0'}),timeout=45) as r:b=r.read()
   dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(b)
   with (ROOT/'logs/downloads.jsonl').open('a') as f:f.write(json.dumps(dict(url=url,path=str(dest.relative_to(ROOT)),utc=stamp(),bytes=len(b),sha256=digest(b)))+'\n')
   return b
  except Exception as e:
   err=str(e)
   if '404' in err or '400' in err:break
   time.sleep(1+attempt)
 with (ROOT/'logs/download_errors.jsonl').open('a') as f:f.write(json.dumps(dict(url=url,path=str(dest.relative_to(ROOT)),utc=stamp(),error=err))+'\n')
 if not optional:raise RuntimeError(f'{url}: {err}')
def run(args,log,stdout=None):
 args=list(map(str,args));start=stamp();p=subprocess.run(args,cwd=ROOT,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
 (ROOT/f'logs/{log}.stdout').write_text(p.stdout);(ROOT/f'logs/{log}.stderr').write_text(p.stderr)
 if stdout:(ROOT/stdout).write_text(p.stdout)
 with (ROOT/'logs/commands.jsonl').open('a') as f:f.write(json.dumps(dict(argv=args,utc=start,end=stamp(),exit_code=p.returncode))+'\n')
 if p.returncode:raise RuntimeError(f'{args[0]}: {p.stderr[-1000:]}')
 return p.stdout
def efetch(db,acc,typ,dest):return fetch('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?'+urlencode(dict(db=db,id=acc,rettype=typ,retmode='text',tool='phirsl1_chma_research')),dest)
def parse_proteome(gb,name):
 rec=SeqIO.read(ROOT/gb,'genbank');rows=[];seqs=[];checks=[]
 for f in rec.features:
  if f.type!='CDS':continue
  q=f.qualifiers;aa=q.get('translation',[''])[0];pid=q.get('protein_id',[''])[0]
  row=dict(protein_id=pid,locus_tag=q.get('locus_tag',[''])[0],gene=q.get('gene',[''])[0],product=q.get('product',[''])[0],length_aa=len(aa),location=str(f.location),pseudo=('pseudo' in q or 'pseudogene' in q),sequence_sha256=digest(aa.encode()));rows.append(row)
  if aa:
   seqs.append(f'>{pid} {row["locus_tag"]} {row["product"]}\n{aa}\n')
   calc=str(f.extract(rec.seq)[int(q.get('codon_start',['1'])[0])-1:].translate(table=int(q.get('transl_table',['11'])[0]))).rstrip('*')
   good=calc==aa or (len(calc)==len(aa) and calc[1:]==aa[1:] and aa[0]=='M')
   checks.append(dict(protein_id=pid,translation_matches=good))
 (ROOT/f'inputs/{name}.faa').write_text(''.join(seqs));table(f'results/{name}_annotations.tsv',rows);table(f'results/{name}_translation_validation.tsv',checks)
 info=dict(accession=rec.id,length_bp=len(rec.seq),taxonomy=rec.annotations.get('taxonomy'),cds_count=len(rows),translated_cds=len(seqs),translation_mismatches=sum(not x['translation_matches'] for x in checks))
 (ROOT/f'results/{name}_genome_summary.json').write_text(json.dumps(info,indent=2));print(info,flush=True);return rec
def acquire():
 efetch('nuccore','NC_010811.2','gbwithparts','inputs/NC_010811.2.gb');parse_proteome('inputs/NC_010811.2.gb','phiRSL1')
 efetch('nuccore','NC_010811.2','fasta_cds_aa','inputs/NC_010811.2_ncbi_cds.faa')
 rows=[]
 for pdb in ['7SQR','7SQV','8FNE']:
  j=json.loads(fetch(f'https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb}/1',f'sources/{pdb}_entity1.json'))
  fetch(f'https://files.rcsb.org/download/{pdb}.cif',f'inputs/pdb/{pdb}.cif')
  ids=j['rcsb_polymer_entity_container_identifiers'];print(pdb,ids.get('uniprot_ids'),flush=True)
  seq=j['entity_poly']['pdbx_seq_one_letter_code_can'].replace('\n','');(ROOT/f'inputs/{pdb}_entity1.faa').write_text(f'>{pdb}_entity1\n{seq}\n')
  rows.append(dict(pdb=pdb,entity=1,length=len(seq),identifiers=json.dumps(ids)))
 table('results/pdb_references.tsv',rows)
 fetch('https://rest.uniprot.org/uniprotkb/search?size=500&'+urlencode(dict(query='organism_name:"Ralstonia phage phiRSL1"',format='tsv',fields='accession,id,protein_name,length,sequence,xref_refseq,organism_id')),'sources/phiRSL1_uniprot.tsv')
 for pmc in ['PMC10299810','PMC9365700']:
  fetch(f'https://www.ebi.ac.uk/europepmc/webservices/rest/{pmc}/fullTextXML',f'sources/{pmc}.xml',optional=True)
def inventory():
 target={str(r.seq):r.id for r in SeqIO.parse(ROOT/'inputs/phiRSL1.faa','fasta')}
 uni=list(csv.DictReader((ROOT/'sources/phiRSL1_uniprot.tsv').open(),delimiter='\t'));print('UniProt:',len(uni),flush=True)
 def one(row):
  acc=row['Entry'];seq=row['Sequence'];pid=target.get(seq,'')
  if not pid:return dict(uniprot=acc,protein_id='',sequence_exact=False,models=0)
  data=fetch(f'https://alphafold.ebi.ac.uk/api/prediction/{acc}',f'sources/afdb_{acc}.json',True);models=json.loads(data) if data else []
  if isinstance(models,dict):models=[]
  for m in models:
   if m.get('uniprotSequence')==seq:fetch(m['pdbUrl'],f'inputs/afdb/{acc}.pdb',True)
  return dict(uniprot=acc,protein_id=pid,sequence_exact=True,models=len(models),downloaded=(ROOT/f'inputs/afdb/{acc}.pdb').exists(),product=row['Protein names'],length=len(seq))
 with ThreadPoolExecutor(max_workers=6) as ex:rows=list(ex.map(one,uni))
 table('results/structure_inventory.tsv',rows,['uniprot','protein_id','sequence_exact','models','downloaded','product','length'])
 print('Structures:',sum(r.get('downloaded',False) for r in rows),'exact sequences',sum(r['sequence_exact'] for r in rows),flush=True)
def checksums():
 table('checksums.tsv',[dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,sha256=digest(p.read_bytes())) for p in sorted(ROOT.rglob('*')) if p.is_file() and 'work' not in p.relative_to(ROOT).parts and p.name!='checksums.tsv'])
if __name__=='__main__':globals()[sys.argv[1]]()
