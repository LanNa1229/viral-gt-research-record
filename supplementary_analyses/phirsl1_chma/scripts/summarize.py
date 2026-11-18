#!/usr/bin/env python3
"""Summarize observed results and independently validate principal claims."""
from pipeline import *
from collections import Counter
import statistics,platform,shutil
def readtsv(path):return list(csv.DictReader((ROOT/path).open(),delimiter='\t'))
def domrows(path):
 rows=[]
 for line in (ROOT/path).read_text().splitlines():
  if not line or line.startswith('#'):continue
  x=line.split(maxsplit=22)
  rows.append(dict(protein_id=x[0],target_length=int(x[2]),profile=x[3],profile_length=int(x[5]),full_evalue=float(x[6]),full_bits=float(x[7]),domain_ievalue=float(x[12]),domain_bits=float(x[13]),hmm_start=int(x[15]),hmm_end=int(x[16]),target_start=int(x[17]),target_end=int(x[18]),profile_coverage=(int(x[16])-int(x[15])+1)/int(x[5]),target_coverage=(int(x[18])-int(x[17])+1)/int(x[2])))
 return rows
def main():
 ann=readtsv('results/phiRSL1_annotations.tsv');annotations={r['protein_id']:r for r in ann}
 refs=list(SeqIO.parse(ROOT/'inputs/chma_characterized.faa','fasta'));pub=list(SeqIO.parse(ROOT/'inputs/chma_published_table7.faa','fasta'))
 ncbi=list(SeqIO.parse(ROOT/'inputs/NC_010811.2_ncbi_cds.faa','fasta'));target=list(SeqIO.parse(ROOT/'inputs/phiRSL1.faa','fasta'))
 blast=readtsv('results/blastp_standard.tsv');blast68=readtsv('results/blastp68.tsv');unmask=readtsv('results/blastp_unmasked.tsv')
 for rows in [blast,blast68,unmask]:
  for r in rows:
   r['query_span_coverage']=(int(r['qend'])-int(r['qstart'])+1)/int(r['qlen']);r['subject_span_coverage']=(int(r['send'])-int(r['sstart'])+1)/int(r['slen'])
 h6=domrows('results/chma_characterized.domtbl');h68=domrows('results/chma68.domtbl');sixframe=domrows('results/chma68_sixframe.domtbl')
 table('results/hmm_chma6_domains.tsv',h6);table('results/hmm_chma68_domains.tsv',h68);table('results/hmm_sixframe_domains.tsv',sixframe)
 fs=readtsv('results/foldseek.tsv');tm=readtsv('results/foldseek_tm.tsv');qc=readtsv('results/structure_model_qc.tsv');mapping=readtsv('results/core_gene_reference_mapping.tsv')
 inv=readtsv('results/structure_inventory.tsv')
 bestblast=min(blast,key=lambda r:float(r['evalue']));best68=min(blast68,key=lambda r:float(r['evalue']));bestunmask=min(unmask,key=lambda r:float(r['evalue']))
 besth68=min(h68,key=lambda r:(r['full_evalue'],r['domain_ievalue']));besth6=min(h6,key=lambda r:(r['full_evalue'],r['domain_ievalue']))
 bestfs=min(fs,key=lambda r:float(r['evalue']));besttm=max(tm,key=lambda r:float(r['qtmscore']))
 besttmhyp=max([r for r in tm if 'Uncharacterized' in r['product']],key=lambda r:float(r['qtmscore']))
 groups=[]
 for cg in sorted(set(r['core_group'] for r in mapping),key=lambda x:int(x[2:])):
  sub=[r for r in mapping if r['core_group']==cg];specific=sub[0]['specific_2023']=='True';br=[r for r in readtsv('results/core_blast.tsv') if r['qseqid'].startswith(cg+'__')];b=min(br,key=lambda r:float(r['evalue'])) if br else None
  hr=domrows('results/core_'+cg+'.domtbl') if specific else [];h=min(hr,key=lambda r:(r['full_evalue'],r['domain_ievalue'])) if hr else None
  groups.append(dict(core_group=cg,specific_2023=specific,mapped_references=sum(r['status']=='mapped' for r in sub),expected_references=len(sub),reference_proteins=';'.join(r['protein_id'] for r in sub if r['protein_id']),blast_best_target=b['sseqid'] if b else '',blast_evalue=b['evalue'] if b else '',blast_identity_pct=b['pident'] if b else '',blast_alignment_length=b['length'] if b else '',blast_query_coverage=(int(b['qend'])-int(b['qstart'])+1)/int(b['qlen']) if b else '',blast_target_coverage=(int(b['send'])-int(b['sstart'])+1)/int(b['slen']) if b else '',hmm_best_target=h['protein_id'] if h else '',hmm_evalue=h['full_evalue'] if h else '',hmm_profile_coverage=h['profile_coverage'] if h else '',hmm_target_coverage=h['target_coverage'] if h else '',blast_E1e5=(float(b['evalue'])<=1e-5) if b else False,hmm_E1e5=(h['full_evalue']<=1e-5) if h else False))
 table('results/core_gene_summary.tsv',groups)
 specific=[g for g in groups if g['specific_2023']]
 table('results/chimalliviridae_specific_markers.tsv',specific)
 # Reference metadata: native sequences; published homolog accessions matched by exact sequence.
 names=['201phi2-1','Goslar','phiPA3','phiKZ','RAY','PCH45'];gps=['gp105','gp189','gp53','gp54','gp222','gp202'];pdbs=['7SQR','7SQV','8FNE','','','']
 sources=['10.1038/s41586-022-05013-4','10.1038/s41586-022-05013-4','10.1038/s41467-023-36526-9','10.1016/j.celrep.2017.07.064','10.1016/j.celrep.2023.112432','10.1038/s41586-019-1786-y']
 refrows=[]
 for r,name,gp,pdb,doi in zip(refs,names,gps,pdbs,sources):
  refrows.append(dict(phage=name,protein=gp,query_accession=r.id,length_aa=len(r.seq),published_table7_accessions=';'.join(p.id for p in pub if p.seq==r.seq),experimental_atomic_pdb=pdb,evidence_doi=doi,native_sequence_sha256=digest(str(r.seq).encode())))
 table('results/reference_sequences.tsv',refrows)
 table('results/blastp_best_by_reference.tsv',[min([r for r in blast if r['qseqid']==ref.id],key=lambda x:float(x['evalue'])) for ref in refs])
 # Preserve weak leads across methods without promoting chance alignments to homologs.
 chosen={bestblast['sseqid'],best68['sseqid'],bestunmask['sseqid'],besth68['protein_id'],besth6['protein_id'],bestfs['protein_id'],besttm['protein_id'],besttmhyp['protein_id']}
 candidates=[]
 for pid in sorted(chosen):
  b=min([r for r in blast68 if r['sseqid']==pid],key=lambda x:float(x['evalue']),default=None);h=min([r for r in h68 if r['protein_id']==pid],key=lambda x:(x['full_evalue'],x['domain_ievalue']),default=None);f=min([r for r in fs if r['protein_id']==pid],key=lambda x:float(x['evalue']),default=None);t=max([r for r in tm if r['protein_id']==pid],key=lambda x:float(x['qtmscore']),default=None)
  r=annotations[pid]
  candidates.append(dict(protein_id=pid,locus_tag=r['locus_tag'],product=r['product'],length_aa=r['length_aa'],blast68_reference=b['qseqid'] if b else '',blast68_identity_pct=b['pident'] if b else '',blast68_alignment_aa=b['length'] if b else '',blast68_query_coverage=b['query_span_coverage'] if b else '',blast68_target_coverage=b['subject_span_coverage'] if b else '',blast68_evalue=b['evalue'] if b else '',hmm68_evalue=h['full_evalue'] if h else '',hmm68_profile_coverage=h['profile_coverage'] if h else '',hmm68_target_coverage=h['target_coverage'] if h else '',foldseek_evalue=f['evalue'] if f else '',foldseek_reference=f['query'] if f else '',foldseek_local_query_tm=f['qtmscore'] if f else '',foldseek_local_target_tm=f['ttmscore'] if f else '',global_tm_reference=t['query'] if t else '',global_query_tm=t['qtmscore'] if t else '',global_target_tm=t['ttmscore'] if t else '',assessment='weak partial profile lead; no model; insufficient for ChmA assignment' if pid=='YP_001950147.1' else 'non-significant match; no evidence of ChmA homology'))
 table('results/candidate_review.tsv',candidates)
 table('results/foldseek_top20_by_evalue.tsv',sorted(fs,key=lambda r:float(r['evalue']))[:20])
 table('results/foldseek_tm_top20.tsv',sorted(tm,key=lambda r:float(r['qtmscore']),reverse=True)[:20])
 summary=dict(proteins=len(target),genome_bp=231255,characterized_references=len(refs),published_chma_references=len(pub),models=len(qc),models_uncharacterized=sum('Uncharacterized' in r['product'] for r in qc),models_coverage=len(qc)/len(target),models_full_exact=sum(r['model_sequence_exact']=='True' for r in qc),mean_model_plddt_median=statistics.median(float(r['mean_plddt']) for r in qc),models_mean_plddt_below70=sum(float(r['mean_plddt'])<70 for r in qc),blast6_best=bestblast,blast68_best=best68,blast_unmasked_best=bestunmask,hmm6_best=besth6,hmm68_best=besth68,sixframe_orfs=sum(1 for r in SeqIO.parse(ROOT/'inputs/phiRSL1_sixframe_min30.faa','fasta')),sixframe_best=min(sixframe,key=lambda r:(r['full_evalue'],r['domain_ievalue'])),foldseek_best=bestfs,foldseek_global_best=besttm,foldseek_global_best_uncharacterized=besttmhyp,specific_groups=len(specific),specific_groups_significant=sum(g['hmm_E1e5'] or g['blast_E1e5'] for g in specific),marker_reference_mapping=[sum(r['status']=='mapped' for r in mapping),len(mapping)],conclusion='Unlikely to use the canonical ChmA-based phage nucleus; moderate confidence. Biological absence is not established.')
 (ROOT/'results/summary.json').write_text(json.dumps(summary,indent=2))
 checks=[]
 def check(name,value,detail=''):
  checks.append(dict(check=name,pass_=bool(value),detail=detail))
 check('343 annotated and extracted translated CDS',len(ann)==len(target)==343)
 check('Fresh NCBI CDS FASTA agrees with GenBank translations',Counter(str(r.seq) for r in target)==Counter(str(r.seq) for r in ncbi))
 check('Every annotated CDS independently translates correctly',all(r['translation_matches']=='True' for r in readtsv('results/phiRSL1_translation_validation.tsv')))
 check('All six native references occur in published homolog set',all(any(r.seq==p.seq for p in pub) for r in refs))
 check('All downloaded structure sequences equal the requested full protein',len(qc)==156 and all(r['model_sequence_exact']=='True' for r in qc))
 check('468 exhaustive pairs in each structural method',len(fs)==len(tm)==3*156)
 check('68 core groups and 21 starred groups recovered',len(groups)==68 and len(specific)==21)
 check('Specific reference panel fully mapped',all(g['mapped_references']==g['expected_references'] for g in specific))
 check('Six withheld ChmA references detected by independent leave-one-out profiles',len(readtsv('results/hmm_leave_one_out.tsv'))==6 and all(float(r['evalue'])<=1e-5 for r in readtsv('results/hmm_leave_one_out.tsv')))
 controls=readtsv('results/foldseek_controls.tsv')
 check('Non-self known ChmA folds match strongly',all(float(r['evalue'])<1e-20 and float(r['qtmscore'])>.8 for r in controls if r['query']!=r['target']))
 check('No ChmA sequence hits at E<=1e-5',not any(float(r['evalue'])<=1e-5 for r in blast+blast68+unmask) and not any(r['full_evalue']<=1e-5 for r in h6+h68))
 check('No specific core hits at E<=1e-5',summary['specific_groups_significant']==0)
 table('results/validation.tsv',checks)
 (ROOT/'results/environment.json').write_text(json.dumps(dict(utc=stamp(),python=sys.version,platform=platform.platform(),biopython=__import__('Bio').__version__,executables={x:shutil.which(x) for x in ['blastp','hmmsearch','hmmbuild']},muscle='<USER_HOME>/miniforge3/envs/foldtree/bin/muscle',foldseek='<USER_HOME>/apps/foldseek/foldseek/bin/foldseek'),indent=2))
 print(json.dumps(summary,indent=2),flush=True)
 print('Validation failures:',[r for r in checks if not r['pass_']],flush=True)
 if not all(r['pass_'] for r in checks):raise RuntimeError('Validation failed')
if __name__=='__main__':main()
