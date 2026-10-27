"""Resolve explicit taxonomic ranks from cached NCBI taxonomy XML."""
import argparse
import json
import urllib.parse
import xml.etree.ElementTree as ET
from pathlib import Path
from .extract import read_tsv,write_tsv,sha_file
from .retrieve import Client


def resolve(task,run):
    task,run=Path(task),Path(run)
    taxids=sorted({x for r in read_tsv(run/'viral_gt_windows.tsv') for x in r.get('taxid','').split(';') if x.isdigit()})
    client=Client(task/'cache/http');lookup={};errors=[]
    for start in range(0,len(taxids),50):
        batch=taxids[start:start+50]
        url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?'+urllib.parse.urlencode(
            dict(db='taxonomy',id=','.join(batch),retmode='xml',tool='viral_gt_mining'))
        try:
            path=client.get(url,lambda data: ET.fromstring(data))
            doc=ET.parse(path)
            for item in doc.getroot().findall('Taxon'):
                canonical=item.findtext('TaxId')
                identities={canonical}|{x.text for x in item.findall('./AkaTaxIds/TaxId')}
                row=dict(resolved_taxid=canonical,source_path=str(path.resolve()),source_sha256=sha_file(path),status='resolved')
                lineage=item.findall('./LineageEx/Taxon')+[item]
                for rank in ['family','genus','species']:
                    entries=[x for x in lineage if x.findtext('Rank')==rank]
                    row[rank]=entries[-1].findtext('ScientificName') if entries else ''
                    row[rank+'_taxid']=entries[-1].findtext('TaxId') if entries else ''
                for ident in identities:
                    if ident in batch:lookup[ident]=dict(row,taxid=ident)
        except Exception as exc:
            errors.append(dict(batch=','.join(batch),reason=str(exc)))
    rows=[lookup.get(t,dict(taxid=t,status='unresolved',resolved_taxid='')) for t in taxids]
    fields=['taxid','resolved_taxid','status','family','family_taxid','genus','genus_taxid','species','species_taxid','source_path','source_sha256']
    write_tsv(run/'taxonomy_cache.tsv',rows,fields)
    write_tsv(run/'taxonomy_retrieval_errors.tsv',errors,['batch','reason'])
    return dict(requested_taxids=len(taxids),resolved=sum(r['status']=='resolved' for r in rows),errors=len(errors))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--task',required=True);p.add_argument('--run',required=True)
    a=p.parse_args();print(json.dumps(resolve(a.task,a.run),indent=2))
