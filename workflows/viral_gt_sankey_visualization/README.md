# Viral GT Figure B: definitive 242-genome Sankeys

This workflow archives only the final Compact (18 terminal groups) and Full (31 groups) views of the same 242 computationally prioritized viral GT-associated glycan-biosynthetic neighborhoods. Counts describe a GT-selected cohort, not all-virus prevalence or experimentally verified glycan pathways. The only figure exclusions from the 245 candidates are KC511810.1, KX960781.1 and KY189087.1, each a single-CDS source. Support is 187 Both, 29 11-only and 26 15-only.

## Existing code and required data

The exact original final generator is already public: [build_sankeys.py at a72e62a](https://github.com/LanNa1229/viral-gt-research-record/blob/a72e62ad7a5efc23b474564a28e8c3c85ed016f7/workflows/viral_gt_candidate_selection/historical/build_sankeys.py), SHA256 `770e6c20050e0f84c62ae54b1e2c6bf8f82ab63de3cc7cae2e240237a37bd335`. It is referenced, not copied into this directory. `historical/legacy_geometry.py` is its missing original dependency, copied unchanged (SHA256 `aee3fdcaaf562f60221e5ba96b766ab4139bf6c0585618426e4e8b30c1a43237`). Its historical `main()`/input-loading functions are not invoked by the replay below.

Upstream taxonomy/screening is the unchanged [Step 6 integrated Master](https://doi.org/10.5281/zenodo.23701289); the exclusion audit is already archived in [Step 4](https://doi.org/10.5281/zenodo.23685561). Neither is repackaged. The companion Step 7 dataset is proposed locally and has no DOI yet. Its 12 files are six approved SVG/PDF/PNG exports, `FigureB_242_genome_to_display_groups.tsv`, Compact/Full node and link TSVs, and `FigureB_style_specification.json`.

The mapping is a 242-row, 14-column projection of the frozen figure membership table: genome/accession IDs, support, lineage/evidence, accepted family/genus and statuses, validated path, and both terminal IDs/labels. Every retained cell is unchanged. It avoids duplicating the upstream Master and unrelated taxonomy/classifier fields. Node/link tables retain exact memberships, counts, coordinates and colors. The style JSON preserves exact heights, typography and 600-dpi settings.

## Replay the frozen final geometry

Use Python 3.11.17 and Matplotlib 3.11.2, preferably the recorded existing environment. An independent environment can use `python -m pip install matplotlib==3.11.2`; this installation command was not executed during preparation. Arial must be legally available to Matplotlib for the original SVG/PNG appearance. The renderer checks the original WSL Windows-font paths. PDF uses built-in Helvetica/Helvetica-Bold and embeds no proprietary font programs.

Download the Step 7 files into a read-only `DATA` directory. From this workflow directory, stage the two original modules in a **new** output directory; never run against the approved export folder:

```sh
DATA=/absolute/path/to/downloaded-step7-files
WORK=/absolute/path/to/new-empty-replay
mkdir -p "$WORK/scripts" "$WORK/qc"
cp historical/legacy_geometry.py "$WORK/scripts/"
curl -fL https://raw.githubusercontent.com/LanNa1229/viral-gt-research-record/a72e62ad7a5efc23b474564a28e8c3c85ed016f7/workflows/viral_gt_candidate_selection/historical/build_sankeys.py -o "$WORK/scripts/build_sankeys.py"
cp "$DATA"/FigureB_compact_nodes.tsv "$DATA"/FigureB_compact_links.tsv \
   "$DATA"/FigureB_full_nodes.tsv "$DATA"/FigureB_full_links.tsv \
   "$DATA"/FigureB_style_specification.json "$WORK/"
cd "$WORK"
python -B - <<'PY'
from pathlib import Path
import importlib.util,json,hashlib
w=Path.cwd()
assert hashlib.sha256((w/'scripts/build_sankeys.py').read_bytes()).hexdigest()=='770e6c20050e0f84c62ae54b1e2c6bf8f82ab63de3cc7cae2e240237a37bd335'
assert hashlib.sha256((w/'scripts/legacy_geometry.py').read_bytes()).hexdigest()=='aee3fdcaaf562f60221e5ba96b766ab4139bf6c0585618426e4e8b30c1a43237'
s=importlib.util.spec_from_file_location('final_sankey',w/'scripts/build_sankeys.py')
m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
style=json.loads((w/'FigureB_style_specification.json').read_text())
assert m.v.P==style['palette'] and m.v.FONTS==style['font_sizes_pt']
assert m.v.K==style['points_per_genome'] and m.v.NODE_W==style['node_width_pt']
for view,nterm in [('compact',18),('full',31)]:
    nodes={}
    for r in m.read(w/f'FigureB_{view}_nodes.tsv'):
        nodes[r['node_id']]=dict(id=r['node_id'],label=r['label'].replace(' / ','\n'),rank=r['taxonomic_rank'],kind=r['node_kind'],stage=float(r['stage']),major=r['major_lineage'],n=int(r['candidate_count']),x=float(r['x_pt']),y=float(r['y_pt']),height=float(r['height_pt']),center=float(r['y_pt'])+float(r['height_pt'])/2,node_color=r['node_hex'],ribbon_color=r['ribbon_hex'])
    links=[dict(source=r['source_node'],target=r['target_node'],n=int(r['candidate_count']),thickness=float(r['thickness_pt']),source_y=float(r['source_y_pt']),target_y=float(r['target_y_pt'])) for r in m.read(w/f'FigureB_{view}_links.tsv')]
    assert nodes['root_viruses']['n']==242
    terminal=[n for k,n in nodes.items() if k.startswith('assignment_')]
    assert len(terminal)==nterm and sum(n['n'] for n in terminal)==242
    m.render(view,nodes,links,style['versions'][view]['height_pt'])
PY
```

This tested entry point calls the unchanged public `render()` on frozen 242-genome geometry; it does not execute the historical `main()`, recompute grouping, or require earlier 245-genome/iTOL intermediates. The historical full main additionally required five unpublished historical table inputs; it is not claimed to run from this minimal dataset. Original source/inputs remain locally preserved.

Expected exports are 180 × 172.303722 mm (Compact) and 180 × 210.058 mm (Full), with 600-dpi PNGs of 4251 × 4070 and 4251 × 4962 pixels. On the recorded environment the replayed PNGs matched the frozen exports pixel-for-pixel. Vector creation timestamps/internal SVG IDs may differ, so regenerated copies do not replace the approved files.

Minimal numerical checks: mapping IDs are unique and identical between views; terminal counts total 242; node/link membership-set sizes equal their counts; incoming/outgoing membership partitions conserve each genome exactly once; heights/thicknesses are 0.65 pt per genome. Prepared validation also checked all stage cuts, source CDS exclusions, 85/111 text boxes, editable SVG text, vector PDFs and actual raster exports. No clipped/overlapping labels or visual changes were observed.

Original independent code uses the [existing MIT license text](https://github.com/LanNa1229/viral-gt-research-record/blob/1f7ff50e5dd982bbb8214672ba06f20ae35a8a7e/workflows/viral_gt_neighborhood_extraction/LICENSES/MIT.txt). Original figures/derived source data are proposed under CC BY 4.0 with upstream-source attribution/qualifications retained. No earlier artwork, font files, database, model weights or scientific intermediate datasets are bundled.
