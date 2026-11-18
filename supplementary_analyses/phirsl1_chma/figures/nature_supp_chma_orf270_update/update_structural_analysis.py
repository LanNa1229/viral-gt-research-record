from pathlib import Path
import csv,json,hashlib,shutil,subprocess,sys
from Bio import SeqIO
from Bio.PDB import MMCIFParser,PDBIO
from Bio.SeqUtils import seq1
R=Path(__file__).resolve().parents[2];O=Path(__file__).resolve().parent
for d in ['inputs','results','logs','work']: (O/d).mkdir(exist_ok=True)
S=Path('<USER_HOME>/Downloads/phirsl1_orf270/seed-1_sample-1/model.cif')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return list(csv.DictReader(p.open(),delimiter='\t'))
def table(p,rows):
 fields=list(dict.fromkeys(k for r in rows for k in r))
 with p.open('w') as f:
  w=csv.DictWriter(f,fields,delimiter='\t');w.writeheader();w.writerows(rows)
for n in ['README.md','METHODS.md','results/candidate_review.tsv','results/hmm_chma68_domains.tsv','results/reference_sequences.tsv','results/foldseek_tm.tsv','results/foldseek_controls.tsv','results/structure_model_qc.tsv','results/validation.tsv']:
 (R/n).read_text()
frozen={r['path']:r['sha256'] for r in read(R/'checksums.tsv')}
assert all(sha(R/n)==v for n,v in frozen.items() if n in ['README.md','METHODS.md','results/foldseek_tm.tsv','results/hmm_chma68_domains.tsv'])
shutil.copy2(S,O/'inputs/ORF270_AF3_sample1.cif')
summary=json.loads((S.parent/'summary_confidences.json').read_text())
shutil.copy2(S.parent/'summary_confidences.json',O/'inputs/summary_confidences.json')
ranking=summary['ranking_score'];assert abs(ranking-.8052457233370919)<.005
st=MMCIFParser(QUIET=True).get_structure('ORF270',str(S))
res=[r for r in st.get_residues() if r.id[0]==' ' and 'CA' in r]
seq=''.join(seq1(r.resname) for r in res)
native=next(r for r in SeqIO.parse(R/'inputs/phiRSL1.faa','fasta') if r.id=='YP_001950147.1')
assert seq==str(native.seq) and len(seq)==265
io=PDBIO();io.set_structure(st);io.save(str(O/'inputs/ORF270_AF3_sample1.pdb'))
for ref in ['7SQR','7SQV','8FNE']:shutil.copy2(R/f'inputs/chma_monomers/{ref}.pdb',O/f'inputs/{ref}.pdb')
FS='<USER_HOME>/apps/foldseek/foldseek/bin/foldseek'
logs=[]
def run(args,name):
 p=subprocess.run([FS]+args,cwd=O,capture_output=True,text=True)
 (O/f'logs/{name}.stdout').write_text(p.stdout);(O/f'logs/{name}.stderr').write_text(p.stderr)
 logs.append(dict(argv=[FS]+args,exit_code=p.returncode))
 (O/'logs/commands.json').write_text(json.dumps(logs,indent=2))
 if p.returncode:raise RuntimeError(p.stderr)
run(['version'],'version')
run(['createdb','inputs/7SQR.pdb','inputs/7SQV.pdb','inputs/8FNE.pdb','work/refs','--threads','4'],'refs')
run(['createdb','inputs/ORF270_AF3_sample1.pdb','work/orf270','--threads','4'],'orf270')
fields='query,target,fident,alnlen,qstart,qend,qlen,tstart,tend,tlen,evalue,bits,qcov,tcov,alntmscore,qtmscore,ttmscore,rmsd,lddt'
for name,db in [('orf270','work/orf270'),('controls','work/refs')]:
 run(['search','work/refs',db,f'work/{name}_aln',f'work/{name}_tmp','--alignment-type','1','--exhaustive-search','1','-a','1','-e','100','--max-seqs','1000','--threads','4'],name+'_search')
 run(['convertalis','work/refs',db,f'work/{name}_aln',f'results/{name}.raw.tsv','--format-output',fields],name+'_convert')
 rows=list(csv.DictReader((O/f'results/{name}.raw.tsv').open(),fieldnames=fields.split(','),delimiter='\t'))
 table(O/f'results/{name}_global.tsv',rows)
old=read(R/'results/foldseek_tm.tsv');assert len(old)==468 and len({r['protein_id'] for r in old})==156
new=read(O/'results/orf270_global.tsv');assert len(new)==3 and {r['query'] for r in new}=={'7SQR','7SQV','8FNE'}
for r in new:r.update(protein_id='YP_001950147.1',uniprot='B2ZY82',product='hypothetical protein; ORF270 AF3',model_source=str(S))
combined=old+new;table(O/'results/foldseek_tm_updated.tsv',combined)
controls=[r for r in read(O/'results/controls_global.tsv') if r['query']!=r['target']]
assert len(controls)==6
best=max(new,key=lambda r:float(r['qtmscore']))
data=dict(model_path=str(S),model_sha256=sha(S),ranking_score_metadata=ranking,ranking_score_user_full_precision=.8052457233370919,sequence_exact=True,length=265,mean_plddt=sum(r['CA'].bfactor for r in res)/265,models=157,comparisons=471,orf270_scores=new,orf270_max=float(best['qtmscore']),orf270_best_ref=best['query'],global_control_range=[min(float(r['qtmscore']) for r in controls),max(float(r['qtmscore']) for r in controls)],source_hashes={n:sha(R/n) for n in ['README.md','METHODS.md','results/foldseek_tm.tsv','results/hmm_chma68_domains.tsv']})
(O/'results/update_summary.json').write_text(json.dumps(data,indent=2))
print(json.dumps(data,indent=2))
