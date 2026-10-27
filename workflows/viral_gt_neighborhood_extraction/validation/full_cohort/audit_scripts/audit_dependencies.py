"""Classify actual file dependencies and retain a clear audit boundary statement."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import json,csv,collections,hashlib
D=Path(__file__).resolve().parents[1]
def write(p,rows,fields):
    with p.open('w',encoding='utf8',newline='') as f:
        w=csv.DictWriter(f,fields,delimiter='\t');w.writeheader();w.writerows(rows)
def main():
    g=json.loads((D/'qc/filesystem_guard_summary.json').read_text());root=Path(g['root']);code=Path(g['code_root']);data=Path(g['data_root']);rows=[]
    for r in g['file_accesses']:
        p=Path(r['path']);category=r['category'];logical=str(p.relative_to(data)) if category=='packaged_data' else str(p.relative_to(code)) if category=='packaged_code' else str(p)
        role='source annotation' if category=='packaged_data' and logical.startswith('sources/genbank/') else 'raw mapping evidence' if category=='packaged_data' and logical.startswith('sources/alias_evidence/') else 'packaged fixed-baseline identity crosswalk' if category=='packaged_data' and logical.startswith('relationships/') else 'packaged anchor or portable input manifest' if category=='packaged_data' else category
        rows.append(dict(path=r['path'],logical_path=logical,category=category,scientific_role=role,write=r['write'],access_count=r['count']))
    write(D/'tables/observed_file_dependencies.tsv',rows,list(rows[0]));pd=[r for r in rows if r['category']=='packaged_data'];pc=[r for r in rows if r['category']=='packaged_code'];denied=[r for r in g['denied_events'] if r['phase']=='extraction'];frozenreads=[r for r in pd if r['logical_path'].startswith('outputs/')];expectedreads=[r for r in pc if r['logical_path'].startswith('tests/fixture/')]
    assert g['succeeded'] and g['canary_blocked'] and not denied and not frozenreads and not expectedreads
    assert not any(x.split('.')[0] in ['torch','esm'] for x in g['imports'])
    assert all(not r['write'] or r['category']=='regenerated_output' for r in rows)
    counts=collections.Counter(r['category'] for r in rows);rolecounts=collections.Counter(r['scientific_role'] for r in pd)
    result=dict(passed=True,isolated_interpreter=g['isolated_interpreter'],user_site_disabled=g['user_site_disabled'],canary_read_outside_packages_blocked=True,algorithm_denied_external_read_attempts=len(denied),algorithm_external_scientific_reads=0,frozen_window_or_CDS_output_files_read_during_extraction=len(frozenreads),test_expected_files_read_during_extraction=len(expectedreads),writes_outside_fresh_output=0,model_imports=[],network_or_subprocess_events_during_extraction=[r for r in denied if r['event']!='open'],unique_access_entries_by_category=dict(counts),packaged_scientific_file_roles=dict(rolecounts),packaged_data_files=[r['logical_path'] for r in pd if r['scientific_role'] not in ['source annotation','raw mapping evidence']],documented_external_dependencies='Python3.12.11 / Biopython1.85 and transitive runtime/OS libraries from existing declared environment',identity_crosswalk_caveat='The archived packaged relationship table is explicitly read to preserve fixed baseline locus IDs; native loci are computed fresh. No frozen window/gene output table is read during extraction.',enforcement_limit='Python audit hooks enforce file opens, directory reads, SQLite database path, subprocess/network and model-import restrictions. This is not a separate kernel container or syscall tracer; no cross-machine execution is claimed.')
    (D/'qc/dependency_validation.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
