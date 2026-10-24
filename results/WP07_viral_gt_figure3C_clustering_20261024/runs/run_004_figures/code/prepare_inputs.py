"""Prepare and hash the fixed scientific inputs, without clustering."""
import common as c
from pathlib import Path
OLD=c.ROOT/c.CONFIG['previous_analysis'];rows=[]
relative=[
'analysis/viral_gt_umap_20261024/tables/representative_15gene_embeddings.npy',
'analysis/viral_gt_umap_20261024/tables/included_locus_representatives.csv',
'analysis/viral_gt_umap_20261024/joint_umap_samples.csv',
'analysis/viral_gt_umap_20261024/tables/locus_registry.csv',
'analysis/viral_gt_umap_20261024/tables/window_locus_anchor_crosswalk.csv',
'analysis/viral_gt_umap_20261024/tables/representative_neighbors_top20.csv',
'analysis/viral_gt_umap_20261024/output_checksums.csv',
'analysis/viral_glycan_benchmark_20261016/mibig_CORE_ONLY_embeddings.npy',
'analysis/viral_glycan_benchmark_20261016/mibig_matched_reference_metadata.tsv',
'analysis/viral_gt_mining_20261017/runs/baseline_594_v2/viral_bgc_embeddings.npy',
'analysis/viral_gt_mining_20261017/runs/baseline_594_v2/embedding_manifest.tsv',
'analysis/viral_gt_interpretation_20261017/tables/genome_relatedness_metadata.tsv',
'analysis/viral_gt_interpretation_20261017/tables/genome_similarity_clusters.tsv',
'analysis/viral_gt_interpretation_20261017/artifact_checksums.tsv']
prior={r['path'].replace('\\','/'):r['sha256'] for r in c.read(c.ROOT/'analysis/viral_gt_interpretation_20261017/artifact_checksums.tsv','\t')}
for name in relative:
 p=c.ROOT/name;h=c.sha(p)
 if name.startswith('analysis/viral_gt_interpretation_20261017/tables/'):
  assert h==prior[p.relative_to(c.ROOT/'analysis/viral_gt_interpretation_20261017').as_posix()]
 rows.append(dict(path=str(p),workspace_relative_path=name,role='scientific input or provenance gate',size_bytes=p.stat().st_size,sha256=h))
for p in sorted(OLD.rglob('*')):
 if p.is_file() and 'cache' not in p.relative_to(OLD).parts and str(p) not in {r['path'] for r in rows}:
  rows.append(dict(path=str(p),workspace_relative_path=p.relative_to(c.ROOT).as_posix(),role='protected previous analytical artifact; not a clustering feature',size_bytes=p.stat().st_size,sha256=c.sha(p)))
c.write('INPUT_MANIFEST.csv',rows)
print('Prepared',len(rows),'input/protected-file identities; no fits.')
