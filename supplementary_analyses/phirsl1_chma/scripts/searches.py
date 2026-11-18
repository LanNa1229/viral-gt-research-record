#!/usr/bin/env python3
"""BLASTP, HMMER and Foldseek against all annotated phiRSL1 proteins."""
from pipeline import *
from Bio.PDB import MMCIF2Dict
from Bio.SeqUtils import seq1
MUSCLE='<USER_HOME>/miniforge3/envs/foldtree/bin/muscle'
FOLDSEEK='<USER_HOME>/apps/foldseek/foldseek/bin/foldseek'
BFIELDS='qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovhsp'
def sequence():
 for tool,args in [('blastp',['blastp','-version']),('hmmer',['hmmsearch','-h']),('muscle',[MUSCLE,'-version']),('foldseek',[FOLDSEEK,'version'])]:run(args,'version_'+tool)
 run(['makeblastdb','-in','inputs/phiRSL1.faa','-dbtype','prot','-out','work/phiRSL1'],'makeblastdb')
 for mode,extra in [('standard',[]),('unmasked',['-seg','no','-comp_based_stats','0'])]:
  run(['blastp','-query','inputs/chma_characterized.faa','-db','work/phiRSL1','-evalue','100','-outfmt','6 '+BFIELDS,'-max_target_seqs','500','-num_threads','4','-out',f'results/blastp_{mode}.raw.tsv']+extra,'blastp_'+mode)
  rows=list(csv.DictReader((ROOT/f'results/blastp_{mode}.raw.tsv').open(),fieldnames=BFIELDS.split(),delimiter='\t'))
  for r in rows:r['target_span_coverage']=round((abs(int(r['send'])-int(r['sstart']))+1)/int(r['slen']),4)
  table(f'results/blastp_{mode}.tsv',rows,BFIELDS.split()+['target_span_coverage'])
 run([MUSCLE,'-align','inputs/chma_characterized.faa','-output','inputs/chma_characterized.aln.faa','-threads','4'],'muscle_chma')
 buildsearch('inputs/chma_characterized.aln.faa','chma_characterized','inputs/phiRSL1.faa')
 buildsearch('inputs/chma_characterized.aln.faa','chma_characterized_control','inputs/chma_characterized.faa')
 refs=list(SeqIO.parse(ROOT/'inputs/chma_characterized.faa','fasta'));loo=[]
 for i,ref in enumerate(refs):
  SeqIO.write([r for r in refs if r.id!=ref.id],ROOT/f'work/loo_{i}.faa','fasta')
  SeqIO.write([ref],ROOT/f'work/holdout_{i}.faa','fasta')
  run([MUSCLE,'-align',f'work/loo_{i}.faa','-output',f'work/loo_{i}.aln','-threads','4'],f'muscle_loo_{i}')
  buildsearch(f'work/loo_{i}.aln',f'loo_{i}',f'work/holdout_{i}.faa')
  raw=(ROOT/f'results/loo_{i}.tbl').read_text()
  for line in raw.splitlines():
   if line and not line.startswith('#'):
    x=line.split();loo.append(dict(held_out=ref.id,evalue=x[4],bitscore=x[5]))
 table('results/hmm_leave_one_out.tsv',loo)
 print('Sequence searches completed',flush=True)
def buildsearch(aln,name,target):
 run(['hmmbuild','--amino','-n',name,f'inputs/{name}.hmm',aln],'hmmbuild_'+name)
 run(['hmmsearch','--max','--cpu','4','-E','100','--domE','100','--tblout',f'results/{name}.tbl','--domtblout',f'results/{name}.domtbl',f'inputs/{name}.hmm',target],'hmmsearch_'+name)
def expanded():
 run([MUSCLE,'-align','inputs/chma_published_table7.faa','-output','inputs/chma68.aln.faa','-threads','4'],'muscle_chma68')
 buildsearch('inputs/chma68.aln.faa','chma68','inputs/phiRSL1.faa')
 buildsearch('inputs/chma68.aln.faa','chma68_control','inputs/chma_characterized.faa')
 print('Expanded HMM completed',flush=True)
def monomers():
 out=ROOT/'inputs/chma_monomers';out.mkdir(exist_ok=True)
 rows=[]
 for pdb in ['7SQR','7SQV','8FNE']:
  d=MMCIF2Dict.MMCIF2Dict(str(ROOT/f'inputs/pdb/{pdb}.cif'));lines=[];res=set();serial=0
  for i,chain in enumerate(d['_atom_site.label_asym_id']):
   if chain!='A' or d['_atom_site.group_PDB'][i]!='ATOM':continue
   pos=int(d['_atom_site.label_seq_id'][i])
   if pdb=='8FNE' and pos<395:continue
   atom=d['_atom_site.label_atom_id'][i]
   alt=d['_atom_site.label_alt_id'][i]
   if alt not in ('.','A'):continue
   serial+=1;res.add(pos);rn=d['_atom_site.label_comp_id'][i]
   x,y,z=[float(d[k][i]) for k in ['_atom_site.Cartn_x','_atom_site.Cartn_y','_atom_site.Cartn_z']]
   num=pos-394 if pdb=='8FNE' else pos
   element=d['_atom_site.type_symbol'][i]
   lines.append(f'ATOM  {serial:5d} {atom:^4s} {rn:>3s} A{num:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 50.00          {element:>2s}\n')
  (out/f'{pdb}.pdb').write_text(''.join(lines)+'TER\nEND\n')
  rows.append(dict(pdb=pdb,chain='A',resolved_residues=len(res),min_label_seq_id=min(res),max_label_seq_id=max(res),tag_excluded=(pdb=='8FNE')))
 table('results/structural_reference_preparation.tsv',rows)
def structure_tm():
 fields='query,target,fident,alnlen,qstart,qend,qlen,tstart,tend,tlen,evalue,bits,qcov,tcov,alntmscore,qtmscore,ttmscore,rmsd,lddt'
 run([FOLDSEEK,'search','work/chma_struct','work/afdb','work/fs_tm','work/fs_tm_tmp','--alignment-type','1','--exhaustive-search','1','-a','1','-e','100','--max-seqs','1000','--threads','4'],'foldseek_tm_search')
 run([FOLDSEEK,'convertalis','work/chma_struct','work/afdb','work/fs_tm','results/foldseek_tm.raw.tsv','--format-output',fields],'foldseek_tm_convert')
 rows=list(csv.DictReader((ROOT/'results/foldseek_tm.raw.tsv').open(),fieldnames=fields.split(','),delimiter='\t'))
 inventory={r['uniprot']:r for r in csv.DictReader((ROOT/'results/structure_inventory.tsv').open(),delimiter='\t')}
 for r in rows:r.update({k:inventory.get(r['target'],{}).get(k,'') for k in ['protein_id','uniprot','product']})
 table('results/foldseek_tm.tsv',rows,fields.split(',')+['protein_id','uniprot','product'])
 print('TM-align comparisons',len(rows),flush=True)
def structure():
 import gzip,shutil
 raw=ROOT/'inputs/afdb_compressed';raw.mkdir(exist_ok=True)
 for p in (ROOT/'inputs/afdb').glob('*.pdb'):
  b=p.read_bytes()
  if b[:2]==b'\x1f\x8b':
   (raw/(p.name+'.gz')).write_bytes(b);p.write_bytes(gzip.decompress(b))
 monomers()
 run([FOLDSEEK,'createdb',*[str(p.relative_to(ROOT)) for p in sorted((ROOT/'inputs/afdb').glob('*.pdb'))],'work/afdb','--threads','4'],'foldseek_createdb_af')
 run([FOLDSEEK,'createdb',*[str(p.relative_to(ROOT)) for p in sorted((ROOT/'inputs/chma_monomers').glob('*.pdb'))],'work/chma_struct','--threads','4'],'foldseek_createdb_chma')
 # Exhaustive comparison bypasses prefilter for this small target set.
 run([FOLDSEEK,'search','work/chma_struct','work/afdb','work/fs_result_trace','work/fs_tmp_trace','-a','1','--exhaustive-search','1','-e','100','--max-seqs','1000','--threads','4'],'foldseek_search')
 fields='query,target,fident,alnlen,qstart,qend,qlen,tstart,tend,tlen,evalue,bits,qcov,tcov,alntmscore,qtmscore,ttmscore,rmsd,lddt'
 run([FOLDSEEK,'convertalis','work/chma_struct','work/afdb','work/fs_result_trace','results/foldseek.raw.tsv','--format-output',fields],'foldseek_convert')
 rows=list(csv.DictReader((ROOT/'results/foldseek.raw.tsv').open(),fieldnames=fields.split(','),delimiter='\t'))
 inventory={r['uniprot']:r for r in csv.DictReader((ROOT/'results/structure_inventory.tsv').open(),delimiter='\t')}
 for r in rows:
  acc=r['target'].replace('.pdb','').split('_')[0];r.update({k:inventory.get(acc,{}).get(k,'') for k in ['protein_id','uniprot','product']})
 table('results/foldseek.tsv',rows,fields.split(',')+['protein_id','uniprot','product'])
 # Positive controls: three known ChmA folds against each other.
 run([FOLDSEEK,'search','work/chma_struct','work/chma_struct','work/fs_control_trace','work/fs_control_trace_tmp','-a','1','--exhaustive-search','1','-e','100','--threads','4'],'foldseek_control_search')
 run([FOLDSEEK,'convertalis','work/chma_struct','work/chma_struct','work/fs_control_trace','results/foldseek_controls.raw.tsv','--format-output',fields],'foldseek_control_convert')
 table('results/foldseek_controls.tsv',csv.DictReader((ROOT/'results/foldseek_controls.raw.tsv').open(),fieldnames=fields.split(','),delimiter='\t'),fields.split(','))
 print('Structural comparisons completed',len(rows),'reported alignments',flush=True)
if __name__=='__main__':globals()[sys.argv[1]]()
