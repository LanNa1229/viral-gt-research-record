#!/usr/bin/env python3
"""Published 68 core gene groups / 21 specific groups (Prichard 2023 Table S3)."""
from pipeline import *
from searches import MUSCLE,BFIELDS,buildsearch
def markers():
 efetch('nuccore','NC_048170.1','gbwithparts','inputs/NC_048170.1.gb');parse_proteome('inputs/NC_048170.1.gb','Goslar')
 text=(ROOT/'sources/PMC10299810/publisher_1.txt').read_text().split('Table S3. Core genome blocks by function,')[1].split('Table S4.')[0]
 text=text[text.index('cg1 gp192'):]
 matches=list(re.finditer(r'\bcg(\d+)\b',text));mapping=[];seqs=[];groups={}
 genomes={name:SeqIO.read(ROOT/'inputs'/gb,'genbank') for name,gb in [('RAY','KU886224.1.gb'),('phiKZ','NC_004629.1.gb'),('Goslar','NC_048170.1.gb')]}
 def getfeature(name,n):
  rec=genomes[name];wanted={'RAY':f'RAY_{n:03d}','phiKZ':f'PHIKZ{n:03d}','Goslar':f'Goslar_{n:05d}'}[name]
  variants={wanted}
  if name=='RAY':variants.add(f'RAY_{n}')
  hits=[f for f in rec.features if f.type=='CDS' and any(v in variants for key in ['locus_tag','old_locus_tag'] for v in f.qualifiers.get(key,[]))]
  if not hits:
   hits=[f for f in rec.features if f.type=='CDS' and any(re.fullmatch(re.escape(wanted),v,re.I) for k in ['gene','product'] for v in f.qualifiers.get(k,[]))]
  return hits
 for i,m in enumerate(matches):
  cg=int(m.group(1));block=text[m.end():matches[i+1].start() if i+1<len(matches) else len(text)]
  labels=re.findall(r'gp(\d+)\*?',block);counts={5:(1,1,2),14:(1,0,1),15:(1,2,1),27:(1,3,2),42:(3,4,3)}.get(cg,(1,1,1))
  if sum(counts)!=len(labels):raise ValueError((cg,labels,block))
  specific='*' in block;pos=0;groups[cg]=[]
  for name,count in zip(['Goslar','phiKZ','RAY'],counts):
   for num in labels[pos:pos+count]:
    hits=getfeature(name,int(num))
    if len(hits)!=1:
     print('Mapping issue',cg,name,num,len(hits),flush=True)
     mapping.append(dict(core_group=f'cg{cg}',specific_2023=specific,phage=name,gp=num,protein_id='',status='unmapped'));continue
    f=hits[0];q=f.qualifiers;pid=q['protein_id'][0];aa=q['translation'][0];ident=f'cg{cg}__{name}__{pid}'
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    groups[cg].append(SeqRecord(Seq(aa),id=ident,description=''))
    mapping.append(dict(core_group=f'cg{cg}',specific_2023=specific,phage=name,gp=num,protein_id=pid,status='mapped',length=len(aa),product=q['product'][0]))
   pos+=count
 table('results/core_gene_reference_mapping.tsv',mapping)
 specificgroups={int(r['core_group'][2:]) for r in mapping if r['specific_2023']}
 print('Core groups',len(groups),'specific groups',len(specificgroups),'mapped',sum(r['status']=='mapped' for r in mapping),'of',len(mapping),flush=True)
 allseq=[s for g in groups.values() for s in g];SeqIO.write(allseq,ROOT/'inputs/core_gene_references.faa','fasta')
 run(['blastp','-query','inputs/core_gene_references.faa','-db','work/phiRSL1','-evalue','100','-outfmt','6 '+BFIELDS,'-max_target_seqs','500','-num_threads','4','-out','results/core_blast.raw.tsv'],'core_blast')
 table('results/core_blast.tsv',csv.DictReader((ROOT/'results/core_blast.raw.tsv').open(),fieldnames=BFIELDS.split(),delimiter='\t'),BFIELDS.split())
 for cg in sorted(specificgroups):
  name=f'core_cg{cg}';SeqIO.write(groups[cg],ROOT/f'inputs/{name}.faa','fasta')
  run([MUSCLE,'-align',f'inputs/{name}.faa','-output',f'inputs/{name}.aln.faa','-threads','4'],'muscle_'+name)
  buildsearch(f'inputs/{name}.aln.faa',name,'inputs/phiRSL1.faa')
 print('Specific core profile searches completed',flush=True)
if __name__=='__main__':markers()
