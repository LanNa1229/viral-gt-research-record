"""Cache public authoritative metadata; never modify screening inputs."""
import pathlib,csv,json,hashlib,urllib.request,urllib.parse,time,datetime,xml.etree.ElementTree as ET,zipfile
OUT=pathlib.Path(__file__).resolve().parents[1];ROOT=OUT.parents[1]
for sub in ['sources/http','sources/ictv','tables','reports','inputs','runs','env','db','logs']:(OUT/sub).mkdir(parents=True,exist_ok=True)
def read(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,rows):
    with p.open('w',encoding='utf8',newline='') as f:
        w=csv.DictWriter(f,list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def fetch(url,p):
    if p.exists() and p.stat().st_size:return p
    for attempt in range(4):
        try:
            with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'viral-taxonomy-reconciliation/1.0'}),timeout=90) as r:data=r.read();headers=dict(r.headers);final=r.url
            p.write_bytes(data)
            pathlib.Path(str(p)+'.json').write_text(json.dumps(dict(url=url,final_url=final,headers=headers,bytes=len(data),sha256=hashlib.sha256(data).hexdigest(),retrieved_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()),indent=2))
            print('Fetched',p.name,len(data),flush=True);time.sleep(.4);return p
        except Exception as e:print('Retry',attempt,url[:100],str(e),flush=True);time.sleep(2**attempt)
    raise RuntimeError(url)
def cached(url):return fetch(url,OUT/'sources/http'/hashlib.sha256(url.encode()).hexdigest())
def xlsx_sheets(p):
    ns={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
    with zipfile.ZipFile(p) as z:
        strings=[]
        if 'xl/sharedStrings.xml' in z.namelist():strings=[''.join(n.itertext()) for n in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('s:si',ns)]
        rels={r.attrib['Id']:r.attrib['Target'] for r in ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))}
        for sheet in ET.fromstring(z.read('xl/workbook.xml')).findall('s:sheets/s:sheet',ns):
            target=rels[sheet.attrib['{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id']];target=target.lstrip('/') if target.startswith('/') else 'xl/'+target;rows=[]
            for row in ET.fromstring(z.read(target)).findall('s:sheetData/s:row',ns):
                values={}
                for c in row.findall('s:c',ns):
                    idx=0
                    for letter in ''.join(x for x in c.attrib['r'] if x.isalpha()):idx=idx*26+ord(letter)-64
                    v=c.findtext('s:v',default='',namespaces=ns)
                    if c.attrib.get('t')=='s':v=strings[int(v)]
                    elif c.attrib.get('t')=='inlineStr':v=''.join(c.find('s:is',ns).itertext())
                    values[idx-1]=v
                rows.append([values.get(i,'') for i in range(max(values,default=-1)+1)])
            yield sheet.attrib['name'],rows
def main():
    master=read(ROOT/'analysis/viral_gt_figureB_taxonomy_20261018/tables/master_genome_dataset.tsv')
    assert len(master)==466 and len({r['genome_unit_id'] for r in master})==466
    accessions=sorted({a for r in master for a in r['source_genome_accessions'].split(';')});groups={}
    for r in master:
        key=(r['source_ncbi_taxid'],r['source_organism_labels']);g=groups.setdefault(key,dict(source_taxid=key[0],organism_labels=key[1],generic_source=r['generic_source_label'],genome_units=0,candidates=0,accessions=set()))
        g['genome_units']+=1;g['candidates']+=r['combined_status']=='candidate';g['accessions'].update(r['source_genome_accessions'].split(';'))
    rows=[]
    for g in groups.values():g['accessions']=';'.join(sorted(g['accessions']));rows.append(g)
    write(OUT/'tables/source_taxid_organism_groups.tsv',rows)
    for name,url in [('VMR_MSL41.v1.20260729.xlsx','https://ictv.global/vmr/current'),('ICTV_Master_Species_List_2025_MSL41.v1.xlsx','https://ictv.global/msl/current')]:
        p=fetch(url,OUT/'sources/ictv'/name)
        for sheet,rs in xlsx_sheets(p):
            safe=''.join(c if c.isalnum() else '_' for c in sheet)
            with (OUT/'sources/ictv'/(p.stem+'__'+safe+'.tsv')).open('w',newline='',encoding='utf8') as f:csv.writer(f,delimiter='\t').writerows(rs)
            print(name,sheet,len(rs),'HEAD',rs[:3],flush=True)
    summaries=[]
    for start in range(0,len(accessions),80):
        url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?'+urllib.parse.urlencode(dict(db='nuccore',id=','.join(accessions[start:start+80]),retmode='json',version='2.0',tool='viral_taxonomy_reconciliation'))
        p=cached(url);doc=json.loads(p.read_text())
        for uid in doc.get('result',{}).get('uids',[]):
            r=doc['result'][uid];r['_source_file']=str(p.relative_to(OUT));summaries.append(r)
    (OUT/'sources/ncbi_accession_summaries.json').write_text(json.dumps(summaries,indent=2))
    tids=sorted({str(r['taxid']) for r in summaries}|{r['source_ncbi_taxid'] for r in master})
    for start in range(0,len(tids),50):cached('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?'+urllib.parse.urlencode(dict(db='taxonomy',id=','.join(tids[start:start+50]),retmode='xml',tool='viral_taxonomy_reconciliation')))
    print('Metadata downloaded:',len(master),'units',len(accessions),'accessions',len(summaries),'summaries',len(tids),'taxids',flush=True)
if __name__=='__main__':main()
