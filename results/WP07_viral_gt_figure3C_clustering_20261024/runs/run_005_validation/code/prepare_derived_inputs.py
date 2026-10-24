import common as c
names=['tables/per_locus_selected287_clustering.csv','tables/community_summary.csv','tables/cross_cohort_locus_memberships.csv','stability/selected287_grid_summary.csv','stability/selected287_HDBSCAN_summary.csv']
c.write('DERIVED_INPUT_MANIFEST.csv',[dict(path=str(c.OUT/name),workspace_relative_path=c.OUT.relative_to(c.ROOT).as_posix()+'/'+name,role='post-fit derived input; never a clustering feature',size_bytes=(c.OUT/name).stat().st_size,sha256=c.sha(c.OUT/name)) for name in names])
import ast
for p in (c.OUT/'scripts').glob('*.py'):ast.parse(p.read_text(encoding='utf-8-sig'))
print('All current source syntax passed; exact derived input hashes prepared.')
