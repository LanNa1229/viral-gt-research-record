# Supplementary computational methods

## Structural inputs and metadata

The recovered curated directory contains 4,854 unique structure stems: 2,656 PDB and 2,198 CIF files. Its current hashes and identifiers are joined to primary node metadata in `inputs/structures/structure_source_manifest.tsv`. These identify the recovered snapshot; no contemporaneous July coordinate-file hash inventory was recovered. Metadata preserve canonical name, file, protein/accession fields, recorded source URLs, taxonomy, family assignment, structure scope and curation notes. Do not infer unrecorded structure-selection rules or assert a complete upstream curation script exists.

The dominant coarse-cluster snapshot contains 594 proteins, all retained in the network: 568 viral-labeled and 26 taxonomy-unresolved. A further 119 viral-labeled nodes lie outside it, comprising 105 baseline-reference-source entries and 14 other source entries. Of the latter 14, six are CAZy-assigned and eight unassigned; source membership is distinct from annotation status. The analysis therefore has 687 viral labels and 4,167 other labels, not 4,167 independently verified nonviral references.

## Foldseek execution and alignment schema

Original successful run identifier: `local_20260727_021051`; Foldseek commit `19800e921de4dcadd646feac40d83b4114137515`; 60 threads. `createdb` generated 4,854 entries. The script constructed all 23,561,316 ordered pairs and imported them with `tsv2db --output-dbtype 7`. Alignment used `structurealign --alignment-type 1 -a 1 -e inf --tmscore-threshold 0.65 --tmscore-threshold-mode 0 -c 0.5 --cov-mode 0`. The `-a 1` backtrace flag was required for rich export. The archived batch script is a submission template; the successful execution's exact local launch wrapper was not recovered.

Rich export used `convertalis --format-mode 4` with the explicit order:

```
query,target,fident,alnlen,evalue,bits,alntmscore,qtmscore,ttmscore,qcov,tcov,qlen,tlen,qstart,qend,tstart,tend
```

The raw file contains a header; the corrected file is headerless. Four ID repairs change 282 data rows and 286 endpoint occurrences, with all other fields unchanged. `fident` is a fraction, not percent. Coverage is fractional; alignment lengths and coordinates are amino-acid positions. The Foldseek TM threshold is not an extra filter requiring every global-minimum weight to be at least 0.65. No E-value cutoff was imposed on the primary graph.

## Undirected graph and community detection

The corrected file has 454,207 directional data rows, including 4,854 self hits. Removal leaves 449,353 rows collapsed into 224,825 unordered pairs. There are 224,528 reciprocal and 297 one-way pairs, no duplicate rows within an orientation and no metadata-unmatched IDs.

For each directional row, the global score is `min(qtmscore, ttmscore)`. The primary undirected weight is the mean of these row scores; a one-way pair uses its only observation. The comparison aligned-core weight is mean `alntmscore`; the third comparison is unweighted. E-value uses the minimum; bit score has mean/max; mutual coverage is mean/min of row-wise `min(qcov, tcov)`. No reciprocal-only or additional analysis TM/coverage filter was used. If duplicated orientations occur in another input, the original code averages rows rather than first balancing orientations.

`python-igraph` constructs the undirected graph including all metadata nodes. `leidenalg.find_partition` uses `RBConfigurationVertexPartition`, resolution 1.0, seed 1 and `n_iterations=-1`. Labels are assigned by community size and are not family identifiers. CAZy/taxonomy are not optimization labels. Original primary totals are 660 communities: 80 multi-protein and 580 singleton communities; the 580 singletons are graph isolates. There are 4,274 connected nodes and a 4,012-node giant component.

Weighted degree is incident-weight sum. Within-community strength sums weights to members of the same community. Within-strength fraction divides by total strength; participation is `1 - sum_c (s_ic/s_i)^2`. Zero-strength fractions and participation are zero in the analyzer. Family bases are parsed from GTFam; subgroup labels are retained. Viral labels are derived from the original Virus group or a viral species-name keyword; no manual override file was supplied.

## Placement v3.0.0

The classifier input is the exact primary full node table plus the full undirected edge table, using community column `community`. Its historical input SHA256 values are retained in `classification_run.json`. No external curated viral-reference metadata was supplied. Thus `unclassified_viral` is an evidence-metadata status, not a statement that no known viral enzymes occur in the dataset.

All thresholds are in `configs/historical/primary_classification_run.json`; the runner passes them explicitly. Reliable family composition requires at least five annotated members and purity at least 0.80, with additional mixed-family rules in the original classifier. Small/large community boundaries are 4/20. Peripheral classification uses within-community strength and the community's 0.25 quantile. Connector screening uses degree, internal/external neighbors, external-strength fraction and participation, followed by bridge subtypes. Viral enrichment is an upper-tail hypergeometric test against all 4,854 nodes/687 viral labels; BH correction is performed across the community summary. Enrichment threshold is q<0.05.

For an unassigned node, direct support requires at least two same-family annotated within-community neighbors, their best-family share of annotated-neighbor strength at least 0.67, and best-minus-second-family share at least 0.20. There are 130 passing nodes: 110 support the community modal family and 20 another family. The 110 are not necessarily all in high-purity communities. Only 57 of the 130 have reliable community-family context (53 concordant, 4 discordant).

Mutually exclusive major placement follows precedence: singleton/isolate; intercommunity connector; peripheral reliable-family context; other reliable-family context; genuine mixed-family; sparse; virus-enriched dark; other annotation-depleted; composition ambiguous; unresolved. The code and exact configuration govern all boundary cases; biological interpretation does not override computational calls.

## Figure source data

Figure 2A's 9,942-edge display is the union of a maximum-weight spanning forest (4,216 edges) and each node's three strongest incident edges (union 9,817); 125 forest edges supplement top-three. Selection sorts by descending global-minimum weight then source/target identifiers. It is a rendering-only subset of the 224,825-edge graph. `derive_figure2a_display.py` calls the unmodified original selection function, avoiding unrecovered historical color/community-label inputs. It makes no claim to reconstruct the original Cytoscape session. The two source maps display 4,274 connected nodes; isolates are separately tabulated.

Figure 2B has 687 protein rows: 304 CAZy-assigned and 383 unassigned. Its plotted unassigned outcomes are 110 concordant direct supports, 73 limited/ambiguous/alternative cases (20 discordant direct +11 ambiguous +42 insufficient-neighbor) and 200 without annotated neighbors in their own community. Assigned outcomes 213/91 compare the assigned family with the modal community family; they are distinct from direct-neighbor support. All 26 nonzero three-stage groups were independently reconciled. See `DATA_DICTIONARY.tsv` and the complete protein-level evidence table in the data archive.

## Sensitivity and historical environment

The stored sweep covers three methods, six resolutions (0.5,0.75,1,1.25,1.5,2) and seeds 1–5. At each method/resolution, the code selects maximum Leiden quality; first-listed seed wins an exact tie. Global-minimum resolution 1 selects seed 5. All-node ARI against original seed 1 is 0.9694049152; 102 proteins, including 36 viral-labeled, change after optimal one-to-one community matching. The stored five-seed ARI summaries cannot be independently recomputed because per-seed partitions/qualities were not retained. No full sweep was rerun for this package.

Historical Python/igraph/leidenalg versions are unknown. The separately locked later environment replays placement from frozen assignments exactly, but fresh clustering has ARI 0.94517455697, despite a byte-identical graph; it yields 126/105/21 total/concordant/alternative direct supports. Those later outputs are validation evidence, not replacement manuscript results. Exact downstream reproduction therefore explicitly uses archived primary assignments. No new Foldseek run was performed.


## Structure-model sources and terms

The source/license map identifies 420 AlphaFold DB models (Monomer v2.0), 2,194 team-generated AlphaFold 3 models (1,797 v3.0.0, 160 v3.0.1, 237 v3.0.2), three experimental PDB references and 2,237 unresolved files. Sources use existing metadata and coordinate headers; no source is inferred from file extensions alone. Exact terms and attributions are in `THIRD_PARTY_NOTICES.md` and `licenses/`. Unresolved-permission files remain prepared locally and are not cleared for public redistribution.
