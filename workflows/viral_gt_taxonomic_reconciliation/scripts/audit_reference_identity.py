"""Identify exact reference-sequence matches so they cannot count as independent support."""
import sys
sys.dont_write_bytecode=True
import hashlib,json,collections
from fetch_metadata import OUT,read,write
master=read(OUT/'tables/phase1_genome_taxonomy_master.tsv');byhash=collections.defaultdict(list)
for r in master:byhash[r['genome_sequence_sha256']].append(r['genome_unit_id'])
config=json.loads((OUT/'db/vcontact3/230.json').read_text())['230']['Domains'];matches=[];counts={}
for domain in ['prokaryotes','eukaryotes']:
    p=OUT/'db/vcontact3'/config[domain]['nucleotides'];rid='';hasher=None;length=0;n=0
    def finish():
        if hasher:
            for gid in byhash.get(hasher.hexdigest(),[]):matches.append(dict(genome_unit_id=gid,reference_genome_id=rid,domain=domain,length=length,nucleotide_sha256=hasher.hexdigest(),source_file=str(p.relative_to(OUT))))
    with p.open() as f:
        for line in f:
            if line.startswith('>'):
                finish();rid=line[1:].split()[0];hasher=hashlib.sha256();length=0;n+=1
            else:
                s=line.strip().upper();hasher.update(s.encode());length+=len(s)
        finish()
    counts[domain]=n
fields=['genome_unit_id','reference_genome_id','domain','length','nucleotide_sha256','source_file']
with (OUT/'tables/exact_reference_sequence_matches.tsv').open('w',encoding='utf8',newline='') as f:
    import csv
    w=csv.DictWriter(f,fields,delimiter='\t');w.writeheader();w.writerows(matches)
(OUT/'reports/reference_sequence_identity_audit.json').write_text(json.dumps(dict(reference_genomes_scanned=counts,exact_match_genome_units=len({r['genome_unit_id'] for r in matches}),matches=len(matches)),indent=2))
print('Reference identity audit:',counts,'exact matches',len(matches))
