#!/usr/bin/env python3
from pipeline import *
def main():
 s=(ROOT/'sources/PMC9365700/41586_2022_5013_MOESM1_ESM.txt').read_text()
 section=s.split('SI Table 7 |')[1].split('SI Table 8 |')[0]
 accs=re.findall(r'\b(?:[A-Z]{2}_\d+|[A-Z]{3}\d+)\.\d+',section)
 table('results/published_chma_accessions.tsv',[dict(accession=a,source='Laughlin et al. 2022 SI Table 7') for a in accs])
 efetch('protein',','.join(accs),'fasta','inputs/chma_published_table7.faa')
 efetch('protein','QFP93061.1','fasta','inputs/PCH45_chma.faa')
 refs=list(SeqIO.parse(ROOT/'inputs/chma_characterized.faa','fasta'))
 if 'QFP93061.1' not in [r.id for r in refs]:
  refs+=list(SeqIO.parse(ROOT/'inputs/PCH45_chma.faa','fasta'));SeqIO.write(refs,ROOT/'inputs/chma_characterized.faa','fasta')
 print('Published ChmA',len(accs),'received',len(list(SeqIO.parse(ROOT/'inputs/chma_published_table7.faa','fasta'))),flush=True)
 # Publisher originals are an independent distribution route for the same study supplements.
 for ext,num in [('pdf',1),('xlsx',3),('csv',2)]:
  for base in ['https://ars.els-cdn.com/content/image/1-s2.0-S2211124723004436-mmc','https://europepmc.org/articles/PMC10299810/bin/NIHMS1905355-supplement-']:
   url=base+str(num)+'.'+ext
   b=fetch(url,f'sources/PMC10299810/publisher_{num}.{ext}',True)
   if b:print(url,len(b),b[:15],flush=True);break
if __name__=='__main__':main()
