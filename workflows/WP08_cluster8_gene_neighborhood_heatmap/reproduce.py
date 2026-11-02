"""Offline validation using recovered code without editing archived data."""
import argparse, csv, hashlib, importlib.util, json, os, platform, re, subprocess, sys
from collections import Counter
from pathlib import Path

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):
    spec=importlib.util.spec_from_file_location(p.stem,p); m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--root',type=Path,default=Path(__file__).resolve().parent)
    ap.add_argument('--out',type=Path,required=True); ap.add_argument('--require-exact-render',action='store_true'); a=ap.parse_args()
    root=a.root.resolve(); out=a.out.resolve()
    if out==root or root/'originals'==out or root/'originals' in out.parents: raise ValueError('Output cannot overwrite originals')
    out.mkdir(parents=True,exist_ok=False)
    import pandas as pd, numpy as np, matplotlib
    from PIL import Image
    src=root/'originals/260630_Temp'; orig=src/'cluster8_global_neighborhood_panelB'
    historical=pd.read_csv(orig/'cluster8_neighborhood_matrix.tsv',sep='\t')
    manifest=pd.read_csv(src/'cluster8_analysis_manifest.fixed.tsv',sep='\t')
    assert len(manifest)==11 and manifest.candidate_id.is_unique
    mod=load(src/'build_cluster8_neighborhood_matrix.py')
    input_audit=[]; all_evidence=[]; rows=[]
    for rec in manifest.to_dict('records'):
        p=src/rec['candidate_id']/'final/analysis_v2/10_database_top_hit_names.tsv'
        df=pd.read_csv(p,sep='\t')
        scores, evidence=mod.infer_scores_for_candidate(df,rec['candidate_id'],rec['cluster8_query_id'])
        dna,glycan,virion=mod.compute_module_scores(scores)
        row=dict(candidate_id=rec['candidate_id'],cluster8_query_id=rec['cluster8_query_id'],analysis_dir=rec['analysis_dir'],**scores,
                 DNA_context_score=dna,glycan_score=glycan,virion_score=virion)
        row['primary_group']=mod.infer_primary_group(row); row['inferred_model']=mod.infer_model(row)
        rows.append(row); all_evidence+=evidence
        input_audit.append(dict(candidate_id=rec['candidate_id'],cluster8_query_id=rec['cluster8_query_id'],
            query_exact_match=bool((df.query_id==rec['cluster8_query_id']).any()),
            query_ids=sorted(df.query_id.unique()),rows=len(df),sha256=sha(p)))
    rebuilt=pd.DataFrame(rows); order=load(src/'plot_cluster8_panelB_heatmap_clean.py').GROUP_ORDER
    rebuilt['group_rank']=rebuilt.primary_group.map(order).fillna(99)
    rebuilt=rebuilt.sort_values(['group_rank','DNA_context_score','glycan_score','virion_score','candidate_id'],ascending=[True,False,False,False,True]).reset_index(drop=True)
    rebuilt.to_csv(out/'rebuilt_matrix.tsv',sep='\t',index=False)
    differences=[]
    for col in historical:
        for i,(x,y) in enumerate(zip(historical[col],rebuilt[col])):
            if x!=y: differences.append(dict(candidate_id=historical.iloc[i].candidate_id,column=col,historical=x,rebuilt=y))
    old_e=pd.read_csv(orig/'cluster8_neighborhood_evidence_long.tsv',sep='\t').fillna('')
    new_e=pd.DataFrame(all_evidence).fillna('')[old_e.columns]
    new_e.to_csv(out/'rebuilt_evidence.tsv',sep='\t',index=False)
    def counter(df):
        return Counter(tuple(str(v) for v in row) for row in df.itertuples(index=False,name=None))
    # Read both serialized TSVs identically to normalize pandas int/float formatting.
    new_e=pd.read_csv(out/'rebuilt_evidence.tsv',sep='\t').fillna('')
    ec_equal=counter(old_e)==counter(new_e)
    field_order=['candidate_id','query_id','protein_id','gene_class','score','database','rank','target_id','hit_name','probability','evalue','q_coverage','t_coverage','aln_length']
    # Order of multi-class evidence is nondeterministic because original classify_hit returns a set.
    hit_issues=[]
    for cand in manifest.candidate_id:
        ev=old_e[old_e.candidate_id==cand]
        for query,g in ev.groupby('query_id'):
            strong=sorted(g.loc[g.score>=2,'gene_class'].unique())
            if len(strong)>1: hit_issues.append(dict(candidate_id=cand,query_id=query,strong_classes=strong))
    pd.DataFrame(input_audit).to_json(out/'input_audit.json',orient='records',indent=2)
    (out/'ambiguous_multiclass_ORFs.json').write_text(json.dumps(hit_issues,indent=2))
    special=old_e[(old_e.candidate_id=='A0A6J5M498') & (old_e.gene_class=='TET_5mYOX')]
    special.to_csv(out/'A0A6J5M498_TET_evidence.tsv',sep='\t',index=False)
    plots=[]
    env=dict(os.environ,MPLBACKEND='Agg',PYTHONDONTWRITEBYTECODE='1',MPLCONFIGDIR=str(out/'mplconfig'))
    original_png=orig/'figures_clean/PanelB_cluster8_global_neighborhood_heatmap_clean.png'
    old=np.array(Image.open(original_png))
    for title in [True,False]:
        dest=out/('plot_with_title' if title else 'plot_no_title')
        cmd=[sys.executable,'-B',str(src/'plot_cluster8_panelB_heatmap_clean.py'),'--matrix',str(orig/'cluster8_neighborhood_matrix.tsv'),'--outdir',str(dest)]
        if not title: cmd+=['--no-title']
        r=subprocess.run(cmd,env=env,capture_output=True,text=True)
        (out/('plot_with_title.log' if title else 'plot_no_title.log')).write_text(r.stdout+r.stderr)
        if r.returncode: raise RuntimeError(r.stderr)
        png=dest/original_png.name; new=np.array(Image.open(png)); same=old.shape==new.shape
        plots.append(dict(title=title,shape=list(new.shape),original_shape=list(old.shape),sha256=sha(png),
                          byte_equal=sha(png)==sha(original_png),pixel_equal=same and bool(np.array_equal(old,new)),
                          different_pixels=int(np.any(old!=new,axis=2).sum()) if same else None))
    provenance=json.loads((root/'audit/provenance.json').read_text())
    preservation=[dict(path=x['archive_path'],matches=sha(root/x['archive_path'])==x['sha256']) for x in provenance['files']]
    pdf_name=original_png.with_suffix('.pdf').name
    original_pdf=(orig/'figures_clean'/pdf_name).read_bytes()
    reproduced_pdf=(out/'plot_with_title'/pdf_name).read_bytes()
    normalize_pdf=lambda b: re.sub(rb'/CreationDate [(][^)]*[)]',b'/CreationDate (DATE)',b)
    pdf_compare=dict(byte_equal=original_pdf==reproduced_pdf,
                     equal_except_creation_date=normalize_pdf(original_pdf)==normalize_pdf(reproduced_pdf),
                     original_creation_dates=[v.decode() for v in re.findall(rb'/CreationDate [(][^)]*[)]',original_pdf)],
                     reproduced_creation_dates=[v.decode() for v in re.findall(rb'/CreationDate [(][^)]*[)]',reproduced_pdf)])
    # Independent aggregation of the saved evidence, without the original scoring helper.
    independent=old_e.groupby(['candidate_id','gene_class']).score.max().unstack(fill_value=0).reindex(index=historical.candidate_id,columns=mod.GENE_CLASSES,fill_value=0)
    independent_equal=bool(np.array_equal(independent.to_numpy(),historical[mod.GENE_CLASSES].to_numpy()))
    import PIL, matplotlib.ft2font
    results=dict(matrix_rows=len(historical),gene_classes=len(mod.GENE_CLASSES),matrix_equal=not differences,matrix_differences=differences,
        independent_evidence_max_equals_matrix=independent_equal,pdf=pdf_compare,
        evidence_rows=len(old_e),rebuilt_evidence_rows=len(new_e),evidence_multiset_equal=ec_equal,
        evidence_missing=sum((counter(old_e)-counter(new_e)).values()),evidence_extra=sum((counter(new_e)-counter(old_e)).values()),
        group_counts=historical.primary_group.value_counts().to_dict(),input_count=len(input_audit),
        query_mismatches=[x for x in input_audit if not x['query_exact_match']],plots=plots,
        originals_preserved=all(x['matches'] for x in preservation),preservation_count=len(preservation),
        environment=dict(python=sys.version,executable=sys.executable,platform=platform.platform(),pandas=pd.__version__,numpy=np.__version__,matplotlib=matplotlib.__version__,pillow=PIL.__version__,freetype=matplotlib.ft2font.__freetype_version__))
    (out/'validation.json').write_text(json.dumps(results,indent=2)); print(json.dumps(results,indent=2))
    assert results['matrix_equal'] and ec_equal and independent_equal and results['originals_preserved']
    if a.require_exact_render:
        assert plots[0]['byte_equal'] and plots[0]['pixel_equal'] and pdf_compare['equal_except_creation_date'], 'Data reproduced; exact render requires the documented matching environment'

if __name__=='__main__': main()
