"""Exhaustive frozen-vs-fresh comparisons; every change and schema omission is reported."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import csv,json,hashlib,collections,datetime
D=Path(__file__).resolve().parents[1];TMP=Path((D/'temporary_workspace.txt').read_text().strip());DATA=TMP/'dataset/step1';OUT=TMP/'regenerated';diffs=[];schemas=[];summaries=[];window_details=[];cds_details=[];file_results=[]
def read(p):
    with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def write(p,rows,fields):
    with Path(p).open('w',encoding='utf8',newline='') as f:
        w=csv.DictWriter(f,fields,delimiter='\t');w.writeheader();w.writerows(rows)
relocation={r['original_path']:r['relocated_path'] for r in read(DATA/'portable_inputs/relocation_manifest.tsv')}
def classify(field,e,a):
    if field in ['source_path','path'] and e in relocation and Path(a).resolve()==(DATA/relocation[e]).resolve():return 'expected_path_relocation'
    if field=='cache_reused':return 'fresh_cache_provenance_difference'
    return 'unexpected_value_difference'
def compare(ep,ap,label,keys,kind='metadata'):
    er,ar=read(ep),read(ap);ek=set(er[0]) if er else set(next(csv.reader(Path(ep).open(),delimiter='\t')));ak=set(ar[0]) if ar else set(next(csv.reader(Path(ap).open(),delimiter='\t')))
    key=lambda r:tuple(Path(r[x]).name if x=='path' else r[x] for x in keys)
    ee={key(r):r for r in er};aa={key(r):r for r in ar};assert len(ee)==len(er) and len(aa)==len(ar),(label,'duplicate keys')
    missing=sorted(ee.keys()-aa.keys());extra=sorted(aa.keys()-ee.keys())
    for k in missing:diffs.append(dict(table=label,key=json.dumps(k),field='ROW',expected='present',actual='absent',category='missing_row'))
    for k in extra:diffs.append(dict(table=label,key=json.dumps(k),field='ROW',expected='absent',actual='present',category='extra_row'))
    for c in sorted(ek-ak):schemas.append(dict(table=label,column=c,status='not_emitted_by_portable_adapter',affected_rows=len(ar),meaning=('model-status field intentionally omitted; inference not rerun' if 'native_loci' in label or 'genome_units' in label else 'retrieval enrichment metadata not appended by extraction adapter' if label=='15/anchor_mapping' else 'pair-map identity/eligibility/verification metadata omitted; core fields checked in window and CDS tables' if label=='11/size_pair_mapping' else 'historical auxiliary metadata') if kind not in ['window','CDS'] else 'unexpected scientific-field omission'))
    for c in sorted(ak-ek):schemas.append(dict(table=label,column=c,status='new_column',affected_rows=len(ar),meaning='new adapter metadata'))
    nchanged=0;unexpected=0;rows_exact=0;rows_exceptpaths=0;field_total=0;field_same=0;sequence_n=0;sequence_matched=0
    for k in sorted(ee.keys()&aa.keys()):
        e,a=ee[k],aa[k];changes=[]
        for f in sorted(ek&ak):
            field_total+=1
            if e[f]==a[f]:field_same+=1;continue
            cat=classify(f,e[f],a[f]);changes.append((f,cat));diffs.append(dict(table=label,key=json.dumps(k),field=f,expected=e[f],actual=a[f],category=cat));nchanged+=1;unexpected+=cat=='unexpected_value_difference'
        exact=not changes;exceptpaths=all(cat=='expected_path_relocation' for _,cat in changes);rows_exact+=exact;rows_exceptpaths+=exceptpaths
        if kind=='window':
            gkeys=['genome_accession','genome_sequence_sha256','anchor_protein_id','anchor_gene_id','anchor_gene_index','anchor_start_1based','anchor_end_1based','anchor_strand','n_genes','n_cds','missing_positions','missing_upstream','missing_downstream','inference_eligible','inference_exclusion_reason','start_1based','end_1based']
            window_details.append(dict(size=label.split('/')[0],window_id=e['window_id'],source_accession=e['genome_accession'],anchor=e['anchor_protein_id'],frozen_gene_count=e['n_genes'],regenerated_gene_count=a['n_genes'],all_scientific_fields_equal=all(e[x]==a[x] for x in gkeys),all_nonpath_fields_equal=exceptpaths,source_path_relocated=bool(changes),shared_fields=len(ek&ak),equal_fields=len(ek&ak)-len(changes),changed_columns=';'.join(x[0] for x in changes)))
        if kind=='CDS':
            applicable=bool(e['sequence']);sequence_n+=applicable;sequence_matched+=bool(applicable and e['sequence']==a['sequence'])
            cds_details.append(dict(size=label.split('/')[0],window_id=e['window_id'],gene_id=e['gene_id'],relative_position=e['relative_position'],gene_index=e['gene_index'],source_feature_index=e['source_feature_index'],protein_accession=e['protein_id'],all_fields_equal=exact,sequence_present=applicable,sequence_equal=e['sequence']==a['sequence'],sequence_sha256_equal=e['sequence_sha256']==a['sequence_sha256'],sequence_identity_fraction=1.0 if applicable and e['sequence']==a['sequence'] else '' if not applicable else 'different',coordinates_equal=all(e[x]==a[x] for x in ['start_1based','end_1based','strand','location','location_parts']),annotation_and_translation_status_equal=all(e[x]==a[x] for x in ['product','translation_status','translation_detail','partial']),changed_columns=';'.join(x[0] for x in changes)))
    order_e=collections.defaultdict(list);order_a=collections.defaultdict(list)
    if kind=='CDS':
        for r in er:order_e[r['window_id']].append(r['gene_id'])
        for r in ar:order_a[r['window_id']].append(r['gene_id'])
    result=dict(table=label,frozen_rows=len(er),fresh_rows=len(ar),missing_rows=len(missing),extra_rows=len(extra),row_keys_equal=not missing and not extra,rows_exact_on_shared_fields=rows_exact,rows_match_after_path_relocation=rows_exceptpaths,shared_columns=len(ek&ak),missing_columns=';'.join(sorted(ek-ak)),extra_columns=';'.join(sorted(ak-ek)),compared_cells=field_total,equal_cells=field_same,value_differences=nchanged,unexpected_value_differences=unexpected,source_protein_sequences_present=sequence_n,source_protein_sequences_identical=sequence_matched,within_window_gene_order_identical=(order_e==order_a) if kind=='CDS' else '',file_byte_identical=sha(ep)==sha(ap),frozen_sha256=sha(ep),fresh_sha256=sha(ap))
    summaries.append(result);file_results.append(dict(table=label,frozen_path=str(ep),fresh_path=str(ap),frozen_bytes=Path(ep).stat().st_size,fresh_bytes=Path(ap).stat().st_size,frozen_sha256=sha(ep),fresh_sha256=sha(ap),byte_identical=sha(ep)==sha(ap)))
    return result

def fasta(p):
    result={};header=None;seq=[]
    for line in Path(p).read_text().splitlines():
        if line.startswith('>'):
            if header is not None:assert header not in result;result[header]=''.join(seq)
            header=line[1:];seq=[]
        else:seq.append(line.strip())
    if header is not None:assert header not in result;result[header]=''.join(seq)
    return result

def main():
    for size,name in [('15','baseline_594_v2'),('11','window11_594_v1')]:
        old=DATA/'outputs'/name;new=OUT/('window'+size)
        compare(old/'viral_gt_windows.tsv',new/'viral_gt_windows.tsv',size+'/windows',['window_id'],'window')
        compare(old/'viral_gt_window_genes.tsv',new/'viral_gt_window_genes.tsv',size+'/CDS',['window_id','gene_id'],'CDS')
        ef,af=fasta(old/'viral_gt_window_proteins.faa'),fasta(new/'viral_gt_window_proteins.faa');assert ef==af,(size,'FASTA record mismatch')
        file_results.append(dict(table=size+'/protein_FASTA',frozen_path=str(old/'viral_gt_window_proteins.faa'),fresh_path=str(new/'viral_gt_window_proteins.faa'),frozen_bytes=(old/'viral_gt_window_proteins.faa').stat().st_size,fresh_bytes=(new/'viral_gt_window_proteins.faa').stat().st_size,frozen_sha256=sha(old/'viral_gt_window_proteins.faa'),fresh_sha256=sha(new/'viral_gt_window_proteins.faa'),byte_identical=sha(old/'viral_gt_window_proteins.faa')==sha(new/'viral_gt_window_proteins.faa')))
        loc=old/('candidate_loci.tsv' if size=='15' else 'native_11gene_candidate_loci.tsv')
        compare(loc,new/'structural_loci.tsv',size+'/native_loci',['locus_id'])
        compare(old/'genome_units.tsv',new/'structural_genome_units.tsv',size+'/genome_units',['genome_unit_id'])
    compare(DATA/'outputs/baseline_594_v2/anchor_mapping.tsv',OUT/'window15/anchor_mapping.tsv','15/anchor_mapping',['input_accession','status','window_id'])
    compare(DATA/'outputs/baseline_594_v2/extraction_qc.tsv',OUT/'window15/extraction_qc.tsv','15/extraction_QC',['window_id'])
    compare(DATA/'outputs/baseline_594_v2/source_manifest.tsv',OUT/'window15/source_manifest.tsv','15/source_manifest',['path'])
    compare(DATA/'outputs/window11_594_v1/window15_window11_mapping.tsv',OUT/'window11/window15_window11_mapping.tsv','11/size_pair_mapping',['window15_id'])
    # Execution-provenance changes are explicit; scientific parameters must still match.
    e=json.loads((DATA/'outputs/baseline_594_v2/extraction_provenance.json').read_text());a=json.loads((OUT/'window15/extraction_provenance.json').read_text());pd=[]
    for f in sorted(set(e)|set(a)):
        if e.get(f)==a.get(f):continue
        cat='expected_path_relocation' if f=='input_path' else 'expected_runtime_measurement' if f=='wall_seconds' else 'unexpected_value_difference'
        pd.append(dict(table='15/extraction_provenance',key='run',field=f,expected=json.dumps(e.get(f)),actual=json.dumps(a.get(f)),category=cat))
    diffs.extend(pd)
    write(D/'tables/window_comparison.tsv',window_details,list(window_details[0]));write(D/'tables/CDS_comparison.tsv',cds_details,list(cds_details[0]));write(D/'tables/table_comparison_summary.tsv',summaries,list(summaries[0]));write(D/'tables/value_differences.tsv',diffs,['table','key','field','expected','actual','category']);write(D/'tables/schema_differences.tsv',schemas,['table','column','status','affected_rows','meaning']);write(D/'tables/file_byte_comparison.tsv',file_results,list(file_results[0]))
    # Record every original output file, including intentionally unregenerated later-stage receipts.
    inventory=[]
    mapped_names={'candidate_loci.tsv':'structural_loci.tsv','native_11gene_candidate_loci.tsv':'structural_loci.tsv','genome_units.tsv':'structural_genome_units.tsv'}
    for size,name in [('15','baseline_594_v2'),('11','window11_594_v1')]:
        for p in sorted((DATA/'outputs'/name).iterdir()):
            if not p.is_file():continue
            target=OUT/('window'+size)/mapped_names.get(p.name,p.name)
            inventory.append(dict(size=size,frozen_file=p.name,regenerated_file=target.name if target.exists() else '',regenerated=target.exists(),note='Compared in detailed tables; structural summaries omit model classification fields' if target.exists() else 'Historical auxiliary/later-stage receipt or input copy not emitted by this extraction-only adapter; remains in packaged archive'))
    write(D/'tables/frozen_output_file_inventory.tsv',inventory,list(inventory[0]))
    unexpected=[r for r in diffs if r['category'] in ['missing_row','extra_row','unexpected_value_difference']]
    scientific_missing=[r for r in schemas if r['meaning']=='unexpected scientific-field omission']
    result=dict(compared_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),scientific_extraction_reproduced=not unexpected and not scientific_missing,window_rows_compared=len(window_details),window_rows_with_all_nonpath_fields_equal=sum(r['all_nonpath_fields_equal'] for r in window_details),CDS_rows_compared=len(cds_details),CDS_rows_with_all_fields_equal=sum(r['all_fields_equal'] for r in cds_details),nonempty_protein_sequences_compared=sum(r['sequence_present'] for r in cds_details),nonempty_protein_sequences_identical=sum(r['sequence_present'] and r['sequence_equal'] for r in cds_details),within_window_order_identical=all(r['within_window_gene_order_identical'] for r in summaries if r['within_window_gene_order_identical']!=''),unexpected_row_or_value_differences=len(unexpected),unexpected_scientific_schema_omissions=len(scientific_missing),reported_value_differences_by_category=dict(collections.Counter(r['category'] for r in diffs)),reported_schema_omissions=len(schemas),tables=summaries,exact_files=[r['table'] for r in file_results if r['byte_identical']],non_byte_identical_files=[r['table'] for r in file_results if not r['byte_identical']])
    (D/'qc/comparison_result.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='tables'},indent=2))
if __name__=='__main__':main()

