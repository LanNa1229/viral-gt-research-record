#!/usr/bin/env python3
"""Additional checks: expanded BLAST, weak lead context, and model quality."""
from pipeline import *
from searches import BFIELDS
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
def main():
 run(['blastp','-query','inputs/chma_published_table7.faa','-db','work/phiRSL1','-evalue','100','-outfmt','6 '+BFIELDS,'-max_target_seqs','500','-num_threads','4','-out','results/blastp68.raw.tsv'],'blastp68')
 table('results/blastp68.tsv',csv.DictReader((ROOT/'results/blastp68.raw.tsv').open(),fieldnames=BFIELDS.split(),delimiter='\t'),BFIELDS.split())
 refs={r.id:str(r.seq) for r in SeqIO.parse(ROOT/'inputs/phiRSL1.faa','fasta')}
 inv=list(csv.DictReader((ROOT/'results/structure_inventory.tsv').open(),delimiter='\t'));rows=[]
 parser=PDBParser(QUIET=True)
 for r in inv:
  if r['downloaded']!='True':continue
  p=ROOT/f'inputs/afdb/{r["uniprot"]}.pdb';structure=parser.get_structure(r['uniprot'],p)
  residues=[res for res in next(structure.get_models()).get_residues() if 'CA' in res and res.id[0]==' ']
  seq=''.join(seq1(res.resname) for res in residues);plddt=[res['CA'].bfactor for res in residues]
  rows.append(dict(protein_id=r['protein_id'],uniprot=r['uniprot'],model_residues=len(seq),protein_length=len(refs[r['protein_id']]),model_coverage=len(seq)/len(refs[r['protein_id']]),model_sequence_exact=seq==refs[r['protein_id']],mean_plddt=sum(plddt)/len(plddt),fraction_plddt70=sum(x>=70 for x in plddt)/len(plddt),product=r['product']))
 table('results/structure_model_qc.tsv',rows)
 ann=list(csv.DictReader((ROOT/'results/phiRSL1_annotations.tsv').open(),delimiter='\t'))
 j=next(i for i,r in enumerate(ann) if r['protein_id']=='YP_001950147.1')
 table('results/ORF270_neighborhood.tsv',ann[j-5:j+6])
 # Examine an annotation-independent six-frame translation as a sensitivity extension.
 # Stops delimit ORFs; min 30 aa. No length filter is applied to the annotated proteome.
 from Bio.SeqRecord import SeqRecord
 rec=SeqIO.read(ROOT/'inputs/NC_010811.2.gb','genbank');orfs=[]
 for strand,nt in [('+',rec.seq),('-',rec.seq.reverse_complement())]:
  for frame in range(3):
   pep=str(nt[frame:len(nt)-((len(nt)-frame)%3)].translate(table=11));pos=0
   for aa in pep.split('*'):
    if len(aa)>=30:
     start=frame+3*pos;end=start+3*len(aa)
     if strand=='-':start,end=len(rec)-end,len(rec)-start
     from Bio.Seq import Seq
     orfs.append(SeqRecord(Seq(aa),id=f'strand{strand}_frame{frame}_nt{start+1}_{end}',description='stop-delimited; no start-codon requirement'))
    pos+=len(aa)+1
 SeqIO.write(orfs,ROOT/'inputs/phiRSL1_sixframe_min30.faa','fasta')
 run(['hmmsearch','--max','--cpu','4','-E','100','--domE','100','--tblout','results/chma68_sixframe.tbl','--domtblout','results/chma68_sixframe.domtbl','inputs/chma68.hmm','inputs/phiRSL1_sixframe_min30.faa'],'hmmsearch_sixframe')
 print('Model QC',len(rows),'exact',sum(r['model_sequence_exact'] for r in rows),'sixframe ORFs',len(orfs),flush=True)
if __name__=='__main__':main()
