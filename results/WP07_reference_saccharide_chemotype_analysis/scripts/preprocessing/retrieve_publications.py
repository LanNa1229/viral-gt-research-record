"""Retrieve primary-publication metadata/abstracts from NCBI PubMed; no annotations inferred."""
import datetime,hashlib,json,time,urllib.parse,urllib.request,xml.etree.ElementTree as ET
from pathlib import Path
ROOT=Path(__file__).resolve().parent; out=ROOT/'raw_sources/pubmed';out.mkdir(exist_ok=True)
cards=json.loads((ROOT/'PRIORITY_SOURCE_RECORDS.json').read_text(encoding='utf-8'))
ids=sorted({p.split(':',1)[1] for r in cards for p in r['raw_record']['cluster']['publications'] if p.startswith('pubmed:')},key=int)
manifest=[];records={}
for batch_index,start in enumerate(range(0,len(ids),100),1):
 selected=ids[start:start+100]
 url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?'+urllib.parse.urlencode(dict(db='pubmed',id=','.join(selected),retmode='xml',tool='ViralGT_source_curation'))
 path=out/f'priority_publications_batch_{batch_index:02d}.xml'
 if not path.exists():
  req=urllib.request.Request(url,headers={'User-Agent':'ViralGT-curation/1.0 (research record retrieval)'})
  with urllib.request.urlopen(req,timeout=40) as response:
   b=response.read();status=response.status
  root=ET.fromstring(b);assert root.tag=='PubmedArticleSet'
  path.write_bytes(b)
 else:b=path.read_bytes();status='cached'
 fetched=datetime.datetime.now(datetime.timezone.utc).isoformat()
 manifest.append(dict(url=url,requested_pmids=selected,path=str(path.relative_to(ROOT)),retrieval_utc=fetched,status=status,sha256=hashlib.sha256(b).hexdigest(),bytes=len(b)))
 root=ET.fromstring(b)
 for article in root.findall('PubmedArticle'):
  pmid=article.findtext('MedlineCitation/PMID')
  def text(e):return ''.join(e.itertext()) if e is not None else ''
  abstract='\n'.join((n.get('Label','')+': ' if n.get('Label') else '')+text(n) for n in article.findall('MedlineCitation/Article/Abstract/AbstractText'))
  doi=[n.text for n in article.findall('PubmedData/ArticleIdList/ArticleId') if n.get('IdType')=='doi']
  records[pmid]=dict(pmid=pmid,title=text(article.find('MedlineCitation/Article/ArticleTitle')),abstract=abstract,doi=doi,
   publication_types=[text(n) for n in article.findall('MedlineCitation/Article/PublicationTypeList/PublicationType')],
   url='https://pubmed.ncbi.nlm.nih.gov/'+pmid+'/',retrieval_batch=str(path.relative_to(ROOT)))
 print('Retrieved batch',batch_index,'records',len(records),flush=True);time.sleep(.4)
(ROOT/'PUBMED_RETRIEVAL_MANIFEST.json').write_text(json.dumps(dict(requested=len(ids),returned=len(records),unreturned=sorted(set(ids)-set(records)),batches=manifest),indent=2)+'\n',encoding='utf-8')
(ROOT/'PUBMED_PRIMARY_RECORDS.json').write_text(json.dumps(records,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print('DONE',len(ids),'requested,',len(records),'returned.')
