"""Render A/B/C from the saved joint coordinates without fitting UMAP."""
import exploratory_umap as a
xy=a.np.load(a.OUT/'tables/joint_umap_coordinates.npy')
rm=a.read(a.REF/'mibig_matched_reference_metadata.tsv')
with (a.OUT/'tables/included_locus_representatives.csv').open(encoding='utf-8-sig',newline='') as f: rows=list(a.csv.DictReader(f))
for r in rows: r['P_at_10']=float(r['P_at_10'])
a.plot(xy,rm,rows)
print('Re-rendered saved coordinates; no UMAP fit or data changes.')

