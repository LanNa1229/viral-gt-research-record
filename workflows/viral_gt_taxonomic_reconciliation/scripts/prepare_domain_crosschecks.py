import sys
sys.dont_write_bytecode=True
import json
import pandas as pd
from Bio import SeqIO
from fetch_metadata import OUT,read,write
markers=read(OUT/'tables/genomad_genome_marker_evidence.tsv');quality={r['genome_unit_id']:r for r in read(OUT/'tables/full_genome_quality.tsv')};checks=[]
for r in markers:
    c=r['no_gt_class'];strong=float(r['no_gt_class_weight_fraction'])>=.5 and int(r['no_gt_class_marker_count'])>=3
    desired='prokaryotes' if c in ['Caudoviricetes','Microviricetes','Faserviricetes','Tokiviricetes','Leviviricetes'] else 'eukaryotes' if c in ['Megaviricetes','Naldaviricetes','Papovaviricetes','Revtraviricetes','Herviviricetes','Poxviricetes'] else ''
    if strong and desired and desired not in quality[r['genome_unit_id']]['vcontact_reference_domains'].split(';'):checks.append(dict(genome_unit_id=r['genome_unit_id'],original_domains=quality[r['genome_unit_id']]['vcontact_reference_domains'],additional_domain=desired,marker_class=c,marker_weight_fraction=r['no_gt_class_weight_fraction'],marker_gene_count=r['no_gt_class_marker_count']))
write(OUT/'tables/reference_domain_crosschecks.tsv',checks)
config=[]
for domain in ['prokaryotes','eukaryotes']:
    ids={r['genome_unit_id'] for r in checks if r['additional_domain']==domain}
    if not ids:continue
    other='eukaryotes' if domain=='prokaryotes' else 'prokaryotes'
    for mode in ['full','no_gt']:
        key=domain+'_crosscheck_'+mode;src=other+'_'+mode
        genes=pd.read_parquet(OUT/'inputs'/(src+'_gene2genome.parquet'));genes=genes[genes['genome_id'].isin(ids)];genes.to_parquet(OUT/'inputs'/(key+'_gene2genome.parquet'))
        lengths=pd.read_parquet(OUT/'inputs'/(src+'_lengths.parquet'));lengths[lengths['genome_id'].isin(ids)].to_parquet(OUT/'inputs'/(key+'_lengths.parquet'))
        proteins=set(genes['protein_id'])
        with (OUT/'inputs'/(key+'.faa')).open('w') as f:SeqIO.write((r for r in SeqIO.parse(OUT/'inputs'/(src+'.faa'),'fasta') if r.id in proteins),f,'fasta')
        config.append(dict(run='vcontact3_'+key,domain=domain,mode=mode,genome_unit_ids=sorted(ids)))
(OUT/'reports/additional_reference_domain_checks.json').write_text(json.dumps(config,indent=2));print('Prepared',len(checks),'cross-domain genomes in',len(config),'full/ablation comparisons')
