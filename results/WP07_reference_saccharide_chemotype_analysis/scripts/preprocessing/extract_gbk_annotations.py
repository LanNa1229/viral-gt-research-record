"""Lossless-source GenBank annotation extraction. No sequence analyses."""
import collections,csv,hashlib,json,re,sys,tarfile
from pathlib import Path
ROOT=Path(__file__).resolve().parent
def parse_gbk(text):
 headers={};references=[];features=[];current_ref=None;feature=None;last_h=None;last_q=None;in_features=False
 for line in text.splitlines():
  if line.startswith('ORIGIN') or line.startswith('//'):break
  if line.startswith('FEATURES'):in_features=True;continue
  if not in_features:
   key=line[:12].strip(); value=line[12:].strip()
   if key=='REFERENCE': current_ref={};references.append(current_ref);last_h=None
   elif key in ('AUTHORS','TITLE','JOURNAL','PUBMED','REMARK') and current_ref is not None:
    current_ref[key]=value;last_h=('reference',key)
   elif key:
    headers[key]=value if key not in headers else headers[key]+' '+value;last_h=('header',key)
   elif last_h:
    owner,k=last_h;target=current_ref if owner=='reference' else headers
    target[k]=target.get(k,'')+' '+value
   continue
  match=re.match(r'^ {5}(\S+)\s+(.+)',line)
  if match:
   feature=dict(type=match[1],location=match[2],qualifiers={});features.append(feature);last_q=None;continue
  if feature is None:continue
  tail=line[21:].strip() if len(line)>21 else ''
  if tail.startswith('/'):
   k,sep,v=tail[1:].partition('=');feature['qualifiers'].setdefault(k,[]).append(v if sep else 'true');last_q=k
  elif tail and last_q:
   q=feature['qualifiers'][last_q]
   q[-1]+= (' ' if last_q!='translation' else '')+tail
  elif tail:feature['location']+=tail
 for feature in features:
  for k,vs in feature['qualifiers'].items():
   feature['qualifiers'][k]=[(v[1:-1].replace('""','"') if len(v)>1 and v.startswith('"') and v.endswith('"') else v) for v in vs]
 organism=headers.get('ORGANISM','')
 # ORGANISM continuation is the documented lineage, not newly retrieved taxonomy.
 return dict(headers=headers,references=references,features=features)
def main():
 records={};manifest=[];priority=json.loads((ROOT/'PRIORITY_SOURCE_RECORDS.json').read_text(encoding='utf-8'))
 wanted={p.stem for p in (ROOT/'raw_sources/mibig_json_3.1').glob('*.json')}
 with tarfile.open(ROOT/'raw_sources/mibig_gbk_3.1.tar.gz','r:gz') as tar:
  for member in tar:
   accession=Path(member.name).stem
   if member.isfile() and accession in wanted and member.name.endswith('.gbk'):
    b=tar.extractfile(member).read();assert accession not in records
    records[accession]=parse_gbk(b.decode('utf-8'))
    manifest.append(dict(mibig_accession=accession,member=member.name,size_bytes=len(b),sha256=hashlib.sha256(b).hexdigest()))
 assert set(records)==wanted and len(records)==2502
 (ROOT/'GBK_PARSED_ANNOTATIONS.json').write_text(json.dumps(records,ensure_ascii=False,separators=(',',':'))+'\n',encoding='utf-8')
 with (ROOT/'GBK_RECORD_MANIFEST.tsv').open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(manifest[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(manifest)
 cards=[]
 for p in priority:
  acc=p['accession'];record=records[acc]
  cds=[f for f in record['features'] if f['type']=='CDS']
  cards.append(dict(accession=acc,GC=p['GC'],original_label=p['original_label'],
   products=[c['compound'] for c in p['raw_record']['cluster']['compounds']],
   source_definition=record['headers'].get('DEFINITION',''),source_comment=record['headers'].get('COMMENT',''),
   publications=p['raw_record']['cluster']['publications'],references=record['references'],
   JSON_genes=p['raw_record']['cluster'].get('genes',{}),saccharide=p['raw_record']['cluster'].get('saccharide',{}),
   GBK_gene_products=[dict(id=f['qualifiers'].get('protein_id',f['qualifiers'].get('locus_tag',[])),gene=f['qualifiers'].get('gene',[]),product=f['qualifiers'].get('product',[])) for f in cds]))
 (ROOT/'PRIORITY_EVIDENCE_CARDS.json').write_text(json.dumps(cards,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
 print('Parsed',len(records),'records;',sum(sum(f['type']=='CDS' for f in r['features']) for r in records.values()),'CDS source annotations.')
 for c in cards:
  print(c['accession'],c['GC'],'; '.join(c['products']),'TITLE:', '; '.join(r.get('TITLE','') for r in c['references'])[:450])
if __name__=='__main__':
 sys.stdout.reconfigure(encoding='utf-8');main()
