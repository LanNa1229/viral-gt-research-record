import sys
sys.dont_write_bytecode=True
import hashlib,json,collections,re
from Bio import SeqIO
import pandas as pd
from fetch_metadata import OUT,ROOT,read,write
def main():
    master=read(OUT/'tables/phase1_genome_taxonomy_master.tsv')
    windows=read(ROOT/'analysis/viral_gt_mining_20261017/runs/baseline_594_v2/viral_gt_windows.tsv')
    ws={w['window_id']:w for w in windows};summaries={r['accessionversion']:r for r in json.loads((OUT/'sources/ncbi_accession_summaries.json').read_text())}
    taxonomy={}
    import xml.etree.ElementTree as ET
    for p in (OUT/'sources/http').iterdir():
        if p.suffix!='.json' and p.read_bytes().lstrip().startswith(b'<?xml'):
            for t in ET.parse(p).getroot().findall('Taxon'):
                taxonomy[t.findtext('TaxId')]=[(x.findtext('Rank'),x.findtext('ScientificName')) for x in t.findall('LineageEx/Taxon')+[t]]
    streams={};proteins=collections.defaultdict(list);lengths=collections.defaultdict(list);quality=[];grows=[]
    for name in ['all_genomes','unresolved_genomes']:
        streams[name]=(OUT/'inputs'/(name+'.fna')).open('w')
    for domain in ['prokaryotes','eukaryotes']:
        for mode in ['full','no_gt']:streams[domain+'_'+mode]=(OUT/'inputs'/(domain+'_'+mode+'.faa')).open('w')
    for g in master:
        gid=g['genome_unit_id'];w0=ws[g['window15_ids'].split(';')[0]];a=w0['genome_accession'];p=__import__('pathlib').Path(w0['source_path'])
        record=next(r for r in SeqIO.parse(p,'genbank') if r.id==a);seq=str(record.seq).upper()
        assert hashlib.sha256(seq.encode()).hexdigest()==g['genome_sequence_sha256'],gid
        title=summaries[a]['title'];anchors={(int(ws[w]['anchor_start_1based']),int(ws[w]['anchor_end_1based']),int(ws[w]['anchor_strand'])) for w in g['window15_ids'].split(';')}
        lin=taxonomy.get(summaries[a]['taxid'] if isinstance(summaries[a]['taxid'],str) else str(summaries[a]['taxid']),[]);names={n for r,n in lin}
        # No host-domain assertion from a generic family label: use supported higher ranks only.
        if names & {'Caudoviricetes','Tokiviricetes','Leviviricetes','Inoviridae','Microviricetes','Naldaviricetes','Faserviricetes'}:domains=['prokaryotes']
        elif names & {'Nucleocytoviricota','Megaviricetes','Peploviricota','Papillomaviridae','Anelloviridae','Naldaviricetes'}:domains=['eukaryotes']
        else:domains=['prokaryotes','eukaryotes']
        # Naldaviricetes are arthropod viruses, not prokaryotic viruses.
        if 'Naldaviricetes' in names:domains=['eukaryotes']
        pending=g['needs_genome_classification']=='True'
        streams['all_genomes'].write('>'+gid+'\n'+seq+'\n')
        if pending:streams['unresolved_genomes'].write('>'+gid+'\n'+seq+'\n')
        ncds=0;nt=0;gtmatches=0;partial=0
        for i,f in enumerate(record.features):
            if f.type!='CDS':continue
            ncds+=1;aa=f.qualifiers.get('translation',[''])[0].replace(' ','').rstrip('*');partial+=int('<' in str(f.location) or '>' in str(f.location))
            if not aa:continue
            nt+=1;pid=gid+'__cds'+str(i);coords=(int(f.location.start)+1,int(f.location.end),f.location.strand or 0);isgt=coords in anchors;gtmatches+=isgt
            grows.append(dict(protein_id=pid,genome_unit_id=gid,source_accession=a,source_protein_id=';'.join(f.qualifiers.get('protein_id',[])),start=coords[0],end=coords[1],strand=coords[2],location=str(f.location),selected_gt_anchor=isgt,translation_length=len(aa),sha256=hashlib.sha256(aa.encode()).hexdigest()))
            if pending:
                for domain in domains:
                    for mode in ['full','no_gt']:
                        if mode=='no_gt' and isgt:continue
                        key=domain+'_'+mode;streams[key].write('>'+pid+'\n'+aa+'\n');proteins[key].append(dict(protein_id=pid,genome_id=gid,keywords='None'))
        if pending:
            for domain in domains:
                for mode in ['full','no_gt']:lengths[domain+'_'+mode].append(dict(genome_id=gid,**{'Size (Kb)':len(seq)/1000}))
        completeness='record_claimed_complete_genome' if re.search(r'complete genome',title,re.I) else 'record_claimed_coding_complete' if re.search(r'complete (?:sequence|coding)',title,re.I) else 'fragment_or_scaffold' if re.search(r'scaffold|contig|partial',title,re.I) else 'unknown_assembly_extent'
        quality.append(dict(genome_unit_id=gid,source_accession=a,cached_genbank=str(p.relative_to(ROOT)),cached_genbank_sha256=hashlib.file_digest(p.open('rb'),'sha256').hexdigest(),nucleotide_sha256=hashlib.sha256(seq.encode()).hexdigest(),genome_length=len(seq),n_ambiguous=sum(c not in 'ACGT' for c in seq),topology=record.annotations.get('topology',''),source_completeness=completeness,n_cds=ncds,n_cached_translations=nt,n_partial_cds=partial,n_selected_gt_anchors=len(anchors),n_gt_proteins_identified=gtmatches,gt_ablation_complete=gtmatches==len(anchors),classification_pending=pending,vcontact_reference_domains=';'.join(domains),genbank_definition=record.description,latest_ncbi_title=title))
    for f in streams.values():f.close()
    for key,rows in proteins.items():pd.DataFrame(rows).to_parquet(OUT/'inputs'/(key+'_gene2genome.parquet'));pd.DataFrame(lengths[key]).to_parquet(OUT/'inputs'/(key+'_lengths.parquet'))
    write(OUT/'tables/full_genome_quality.tsv',quality);write(OUT/'tables/cached_full_genome_proteins.tsv',grows)
    print(json.dumps(dict(all_genomes=len(quality),pending=sum(r['classification_pending'] for r in quality),domains={d:sum(r['classification_pending'] and d in r['vcontact_reference_domains'] for r in quality) for d in ['prokaryotes','eukaryotes']},nt_total=sum(r['genome_length'] for r in quality),proteins=len(grows),gt_ablation_complete=sum(r['gt_ablation_complete'] for r in quality)),indent=2))
if __name__=='__main__':main()
