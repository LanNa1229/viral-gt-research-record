"""Cached accession crosswalk and batch NCBI retrieval; no structural searches.

UniProt/UniSave EMBL cross-references bound mapping scope. Protein sequences must
match source CDSs exactly in extraction. Missing cross-references remain unresolved;
this is not an exhaustive sequence search or an IPG occurrence census.
"""
import argparse
import concurrent.futures
import io
import json
import re
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from Bio import SeqIO, SwissProt
from .extract import digest, sha_file, read_tsv, write_tsv, json_write, load_anchors


class Client:
    def __init__(self, cache):
        self.cache=Path(cache);self.cache.mkdir(parents=True,exist_ok=True)
        self.lock=threading.Lock();self.next_request={}

    def get(self,url,validator):
        key=digest(url);path=self.cache/key;meta=path.with_suffix('.json')
        if path.exists() and meta.exists():
            m=json.loads(meta.read_text())
            if m['sha256']==sha_file(path):
                validator(path.read_bytes());return path
        domain=urllib.parse.urlsplit(url).netloc
        for attempt in range(4):
            with self.lock:
                delay=max(0,self.next_request.get(domain,0)-time.monotonic())
                self.next_request[domain]=time.monotonic()+delay+(.37 if 'ncbi' in domain else .15)
            if delay:time.sleep(delay)
            try:
                req=urllib.request.Request(url,headers={'User-Agent':'viral-gt-mining/0.1 scientific-research'})
                with urllib.request.urlopen(req,timeout=120) as response:data=response.read()
                validator(data)
                tmp=path.with_suffix('.part');tmp.write_bytes(data);tmp.replace(path)
                json_write(meta,dict(url=url,sha256=sha_file(path),bytes=len(data),
                    retrieved_utc=datetime.now(timezone.utc).isoformat()))
                return path
            except urllib.error.HTTPError as e:
                if e.code in (400,404,410):raise
                if attempt==3:raise
                time.sleep(min(30,2**attempt*2))
            except Exception:
                if attempt==3:raise
                time.sleep(min(30,2**attempt*2))


def validate_json(data):
    obj=json.loads(data)
    if 'messages' in obj or 'error' in obj:raise ValueError('Remote JSON error')


def validate_swiss(data):
    if not list(SwissProt.parse(io.StringIO(data.decode()))):raise ValueError('No UniSave records')


def validate_genbank(data):
    if not list(SeqIO.parse(io.StringIO(data.decode()),'genbank')):raise ValueError('No genomic records')


def crosslinks(entry):
    rows=[]
    for x in entry.get('uniProtKBCrossReferences',[]):
        if x['database']=='EMBL':
            p={v['key']:v['value'] for v in x.get('properties',[])}
            pid=p.get('ProteinId','')
            if pid and pid!='-':rows.append((x['id'],pid))
        elif x['database']=='RefSeq':
            for v in x.get('properties',[]):
                if v['key']=='NucleotideSequenceId' and v['value']!='-':rows.append((v['value'],x['id']))
    return sorted(set(rows))


def resolve_anchor(accession,client,local_records):
    """Use current record, or the most recent saved UniSave entry on deletion."""
    evidence=None;method='current_uniprot';error=''
    local=local_records.get(accession,{})
    try:
        if local.get('json'):
            evidence=local['json'];entry=json.loads(evidence.read_text())
        else:
            evidence=client.get(f'https://rest.uniprot.org/uniprotkb/{accession}.json',validate_json)
            entry=json.loads(evidence.read_text())
        if 'sequence' not in entry:raise ValueError('Inactive accession without current sequence')
        sequence=entry['sequence']['value'];links=crosslinks(entry)
        taxid=entry.get('organism',{}).get('taxonId','')
    except Exception as exc:
        error=str(exc);method='archived_unisave_latest_available'
        try:
            evidence=local.get('unisave') or client.get(f'https://rest.uniprot.org/unisave/{accession}?format=txt',validate_swiss)
            with evidence.open() as f:records=list(SwissProt.parse(f))
            records=[r for r in records if accession in r.accessions]
            if not records:raise ValueError('No matching accession in archive')
            entry=max(records,key=lambda r:int(r.annotation_update[1]))
            sequence=entry.sequence;taxid=';'.join(entry.taxonomy_id)
            links=sorted({(x[1],x[2]) for x in entry.cross_references
                          if x[0]=='EMBL' and len(x)>2 and x[2]!='-'})
        except Exception as archive_exc:
            return dict(input_accession=accession,status='unresolved',reason=f'current: {error}; archive: {archive_exc}'),[]
    links=[(g,p) for g,p in links if re.fullmatch(r'[A-Za-z_]+\d+(?:\.\d+)?',g) and re.fullmatch(r'[A-Za-z_]+\d+(?:\.\d+)?',p)]
    common=dict(input_accession=accession,evidence_path=str(evidence.resolve()),evidence_sha256=sha_file(evidence),
        sequence_sha256=digest(sequence),taxid=taxid,mapping_method=method)
    rows=[dict(common,genome_accession=g,protein_id=p) for g,p in links]
    return dict(common,status='crosslinks_available' if links else 'unresolved',reason='' if links else 'no_usable_genomic_crossreferences',
                n_genome_links=len(links)),rows


def run(root,task,workers=4):
    root,task=Path(root),Path(task);inputs=task/'inputs';inputs.mkdir(exist_ok=True)
    client=Client(task/'cache/http');started=time.perf_counter()
    accessions=load_anchors(inputs/'baseline_594_anchors.txt')
    assert len(accessions)==len(set(accessions))==594
    local={}
    dirs=[root/'analysis'/d/'sources' for d in ['viral_glycan_benchmark_20261016',
        'viral_glycan_benchmark_A0A0K0N728_20261016','viral_glycan_benchmark_cluster7_ARF08726_ARF08739_20261016']]
    for directory in dirs:
        for p in directory.glob('*.json'):
            if p.stem in accessions:local.setdefault(p.stem,{})['json']=p
        for p in directory.glob('*_unisave.txt'):
            local.setdefault(p.name.split('_')[0],{})['unisave']=p
    metadata_dir=task/'cache/anchor_metadata';metadata_dir.mkdir(parents=True,exist_ok=True)
    metadata=[];aliases=[]
    def worker(acc):
        saved=metadata_dir/(acc+'.json')
        if saved.exists():
            obj=json.loads(saved.read_text())
            evidence=obj['metadata'].get('evidence_path')
            if evidence and Path(evidence).exists() and sha_file(evidence)==obj['metadata']['evidence_sha256']:
                return obj['metadata'],obj['aliases']
        m,a=resolve_anchor(acc,client,local)
        json_write(saved,dict(metadata=m,aliases=a));return m,a
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        futures={pool.submit(worker,acc):acc for acc in accessions}
        for future in concurrent.futures.as_completed(futures):
            m,a=future.result();metadata.append(m);aliases.extend(a)
            if len(metadata)%20==0:
                print(f'Crosswalk {len(metadata)}/594; mapped links {sum(x["status"]=="crosslinks_available" for x in metadata)}',flush=True)
                json_write(task/'retrieval_progress.json',dict(stage='uniprot_crosswalk',completed=len(metadata),total=594))
    metadata.sort(key=lambda r:r['input_accession']);aliases.sort(key=lambda r:(r['input_accession'],r['genome_accession'],r['protein_id']))
    write_tsv(inputs/'retrieval_anchor_status.tsv',metadata)
    write_tsv(inputs/'retrieved_aliases.tsv',aliases,['input_accession','genome_accession','protein_id','sequence_sha256',
        'evidence_path','evidence_sha256','taxid','mapping_method'])
    # Reuse validated local genomes; never refetch the same accession repeatedly.
    genomes={}
    for directory in dirs:
        for path in directory.glob('*.gb'):
            with path.open() as f:
                for rec in SeqIO.parse(f,'genbank'):genomes[rec.id.split('.')[0]]=str(path.resolve())
    needed=sorted({r['genome_accession'] for r in aliases if r['genome_accession'].split('.')[0] not in genomes})
    genome_dir=task/'cache/genomes';genome_dir.mkdir(parents=True,exist_ok=True)
    for path in genome_dir.glob('*.gb'):
        with path.open() as f:
            for rec in SeqIO.parse(f,'genbank'):genomes[rec.id.split('.')[0]]=str(path.resolve())
    needed=[g for g in needed if g.split('.')[0] not in genomes]
    failures=[];batches=[needed[i:i+5] for i in range(0,len(needed),5)]
    def download(batch):
        url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?'+urllib.parse.urlencode(
            dict(db='nuccore',id=','.join(batch),rettype='gbwithparts',retmode='text',tool='viral_gt_mining'))
        try:
            raw=client.get(url,validate_genbank)
            with raw.open() as f:records=list(SeqIO.parse(f,'genbank'))
            found={}
            for rec in records:
                if rec.id.split('.')[0] not in {g.split('.')[0] for g in batch}:continue
                path=genome_dir/(rec.id+'.gb')
                temp=path.with_suffix('.tmp')
                with temp.open('w') as f:SeqIO.write(rec,f,'genbank')
                temp.replace(path);found[rec.id.split('.')[0]]=str(path.resolve())
            return found,[dict(genome_accession=g,reason='not_returned_by_efetch') for g in batch if g.split('.')[0] not in found]
        except Exception as exc:
            # One invalid accession must not discard every other member of its batch.
            if len(batch)>1:
                found={};errors=[]
                for g in batch:
                    f,e=download([g]);found.update(f);errors.extend(e)
                return found,errors
            return {},[dict(genome_accession=batch[0],reason=str(exc))]
    with concurrent.futures.ThreadPoolExecutor(max_workers=min(workers,3)) as pool:
        for i,(found,errors) in enumerate(pool.map(download,batches),1):
            genomes.update(found);failures.extend(errors)
            print(f'Genome batches {i}/{len(batches)}; available {len(genomes)}; failed {len(failures)}',flush=True)
            json_write(task/'retrieval_progress.json',dict(stage='ncbi_genomes',completed_batches=i,total_batches=len(batches),available_genomes=len(genomes),failed_genomes=len(failures)))
    wanted={r['genome_accession'].split('.')[0] for r in aliases}
    paths=sorted({p for g,p in genomes.items() if g in wanted})
    write_tsv(inputs/'full_genbanks.tsv',[dict(path=p) for p in paths],['path'])
    write_tsv(inputs/'genome_retrieval_failures.tsv',failures,['genome_accession','reason'])
    result=dict(stage='retrieval_complete',anchors=594,with_crosslinks=sum(r['status']=='crosslinks_available' for r in metadata),
        unresolved_crosslinks=sum(r['status']=='unresolved' for r in metadata),alias_links=len(aliases),genome_files=len(paths),
        failed_genome_accessions=len(failures),elapsed_seconds=time.perf_counter()-started,
        mapping_scope='all_explicit_EMBL_and_current_RefSeq_links; sequence_verified_at_extraction; no_exhaustive_IPG_search',
        http_cache=str(client.cache.resolve()))
    json_write(task/'retrieval_progress.json',result);print(json.dumps(result,indent=2),flush=True)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',required=True);p.add_argument('--task',required=True)
    p.add_argument('--workers',type=int,default=4);a=p.parse_args();run(a.root,a.task,a.workers)


if __name__=='__main__':main()
