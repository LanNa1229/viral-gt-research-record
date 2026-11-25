"""Finalize report/asset checks and immutable output inventories after visual QA."""
from common import *
import report,validate,csv,struct

def run():
    report.run();checks=0
    def check(ok):
        nonlocal checks
        assert ok
        checks+=1
    for folder in [OUT/'tables',RUN,CTRL]:
        for p in folder.glob('*.tsv'):
            with p.open(encoding='utf-8-sig') as f:head=next(csv.reader(f,delimiter='\t'))
            check(len(head)==len(set(x.casefold() for x in head)))
    summary=read(OUT/'tables/candidate_recovery_summary.tsv')
    check(int(next(r['candidate_in_either_condition'] for r in summary if r['set']=='C' and r['level']=='locus'))==288)
    check(int(next(r['candidate_in_either_condition'] for r in summary if r['set']=='C' and r['level']=='genome'))==245)
    controls=read(OUT/'tables/background_genome_pair_values.tsv');contrast=read(OUT/'tables/window_length_background_contrasts.tsv')
    for r in contrast:
        rr=[c for c in controls if c['method']==r['method'] and c['metric']==r['metric']]
        check(len(rr)==int(r['n_genomes']));check(len({c['genome_unit_id'] for c in rr})==len(rr))
        check(abs(np.mean([float(c['difference_in_differences']) for c in rr])-float(r['difference_in_differences']))<1e-12)
    for k in [5,10,20]:
        rows=[r for r in read(OUT/'tables/neighbor_identity_turnover_frequencies.tsv') if int(r['k'])==k]
        check(sum(int(r['n_pairs_gaining_exactly_n']) for r in rows)==737)
        check(sum(int(r['n_pairs_losing_exactly_n']) for r in rows)==737)
    figures=[]
    for p in sorted((OUT/'figures').glob('*.png')):
        data=p.read_bytes();check(data.startswith(b'\x89PNG'));w,h=struct.unpack('>II',data[16:24])
        check(w>1000 and h>1000)
        for ext in ['pdf','svg']:
            other=p.with_suffix('.'+ext);check(other.exists())
            check(other.read_bytes().startswith(b'%PDF') if ext=='pdf' else '<svg' in other.read_text(encoding='utf-8'))
        figures.append(dict(name=p.stem,width_px=w,height_px=h,PNG_sha256=sha(p),PDF_sha256=sha(p.with_suffix('.pdf')),SVG_sha256=sha(p.with_suffix('.svg'))))
    check(len(figures)==6)
    for name in ['REPORT.md','README.md','METHODS.md']:
        content=(OUT/name).read_text(encoding='utf-8')
        for target in re.findall(r'\[[^\]]+\]\(([^)]+)\)',content):
            if not target.startswith(('http:','https:','#')):check((OUT/target.split('#')[0]).exists())
    save(OUT/'figure_exports.json',dict(status='VISUALLY_REVIEWED_PNG_EXPORTS',figures=figures,review='figures/visual_review.md',
        PDF_SVG='same Matplotlib objects; headers checked; not separately rasterized'))
    save(OUT/'final_delivery_validation.json',dict(status='PASS',delivery_checks=checks,frozen_input_verification=verify_prior(),
        independent_inference_and_accounting_checks=json.loads((OUT/'validation.json').read_text())['independent_checks'],
        figures=6,model_inference_not_repeated_for_finalization=True))
    validate.manifests()
    # Reverify the separate protein inventory as well as its own protected manifest hash.
    proteins=read(OUT/'protein_cache_checksums.tsv');check(all(sha(OUT/r['path'])==r['sha256'] for r in proteins))
    print(json.dumps(dict(status='PASS',delivery_checks=checks,protein_cache_files_verified=len(proteins),
        target_run=str(RUN),report=str(OUT/'REPORT.md')),indent=2),flush=True)
if __name__=='__main__':run()
