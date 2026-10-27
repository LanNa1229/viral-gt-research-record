"""Verify the original headerless baseline membership, then write immutable cohorts."""
import argparse
import csv
import json
import subprocess
from pathlib import Path
from viral_gt_mining.extract import sha_file, write_tsv, json_write


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',required=True);p.add_argument('--task',required=True)
    a=p.parse_args();root=Path(a.root);task=Path(a.task);task.mkdir(exist_ok=True,parents=True)
    checkout=task/'checkout'
    source=root/'../projects/viral-gt-research-record/results/WP01B_coarse_clustering/canonical/cluster_membership.tsv'
    source=source.resolve();copy=checkout/'results/WP01B_coarse_clustering/canonical/cluster_membership.tsv'
    assert sha_file(source)==sha_file(copy),'Checkout baseline differs from original'
    with source.open() as f:
        rows=list(csv.reader(f,delimiter='\t'))
    assert all(len(r)==2 for r in rows)
    cluster=[r[1] for r in rows if r[0]=='A0A3Q8Q3U2']
    assert len(cluster)==len(set(cluster))==594 and 'A0A3Q8Q3U2' in cluster
    counts={r[0]:sum(x[0]==r[0] for x in rows) for r in rows}
    assert max(counts,key=counts.get)=='A0A3Q8Q3U2'
    cohort=task/'inputs';cohort.mkdir(exist_ok=True)
    text='\n'.join(sorted(cluster))+'\n';dest=cohort/'baseline_594_anchors.txt'
    if dest.exists():assert dest.read_text()==text
    else:dest.write_text(text)
    write_tsv(cohort/'baseline_594_anchors.tsv',[dict(accession=x,representative='A0A3Q8Q3U2',source_column='2 (member; headerless)',source_sha256=sha_file(source)) for x in sorted(cluster)])
    # Independent prespecified benchmark/control cases; never expand the full cohort.
    pilot=['QBP05883.1','ARF08732.1','ARF08726.1','ARF08731.1','ARF08736.1','ARF08737.1','AKJ72542.1']
    (cohort/'pilot_anchors.txt').write_text('\n'.join(pilot)+'\n')
    paths=[root/'analysis/viral_glycan_benchmark_20261016/sources/MK016663.1.gb',
        root/'analysis/viral_glycan_benchmark_cluster7_ARF08726_ARF08739_20261016/sources/KY684083.1.gb',
        root/'analysis/viral_glycan_benchmark_A0A0K0N728_20261016/sources/KR063281.1.gb']
    write_tsv(cohort/'pilot_genbanks.tsv',[dict(path=str(x)) for x in paths])
    report=dict(status='PASS',source=str(source),source_sha256=sha_file(source),
        representative_column='1 (headerless)',accession_column='2 (headerless)',representative='A0A3Q8Q3U2',
        n_rows=len(rows),n_clusters=len(counts),selected_count=len(cluster),unique_selected_count=len(set(cluster)),
        representative_included=True,largest_cluster=True,baseline_provenance='WP01B README excludes positive-control spike-in reruns',
        accession_list_sha256=sha_file(dest),pilot_anchors=pilot,
        baseline_commit=subprocess.check_output(['git','-C',str(checkout),'rev-parse','HEAD'],text=True).strip())
    json_write(cohort/'cohort_validation.json',report);print(json.dumps(report,indent=2))


if __name__=='__main__':main()
