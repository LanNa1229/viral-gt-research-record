"""Build Methods, evidence matrix and readout from recorded executions and saved results."""
import collections,datetime,json
import common as c
runs=sorted((c.OUT/'runs').glob('run_*'));before=json.loads((runs[0]/'before.json').read_text());after=json.loads((runs[0]/'after.json').read_text())
assert after['exit_status']==0 and after['inputs_unchanged'] and after['executable_sources_unchanged']
results=json.loads((c.OUT/'RESULTS.json').read_text());communities=c.read(c.OUT/'tables/community_summary.csv');comparison=c.read(c.OUT/'tables/cross_cohort_partition_agreement.csv');hd=c.read(c.OUT/'stability/selected287_HDBSCAN_summary.csv');graph=next(r for r in c.read(c.OUT/'graphs/selected287_diagnostics.csv') if r['k']=='20');seed=next(r for r in c.read(c.OUT/'stability/selected287_seed_ARI.csv') if r['k']=='20' and float(r['resolution'])==1);grid=c.read(c.OUT/'stability/selected287_grid_summary.csv');joint=c.read(c.OUT/'tables/joint_community_membership_summary.csv');contacts=c.read(c.OUT/'tables/joint_graph_viral_reference_contacts.csv');validation=json.loads((c.OUT/'VALIDATION.json').read_text());assert validation['status']=='PASS'
bycomp={r['comparison']:r for r in comparison};sizes=', '.join(r['Leiden_cluster']+'='+r['n_loci'] for r in communities);mcs15=next(r for r in hd if r['run_id']=='mcs15_ms5');small=next(r for r in hd if r['run_id']=='mcs10_ms5');robust=[r['Leiden_cluster'] for r in communities if r['operational_robust_flag']=='True'];minc=min(int(r['n_communities']) for r in grid);maxc=max(int(r['n_communities']) for r in grid)
summarytable='| Community | Loci | Genome units | Full 15 genes | Median P@10 | Seed Jaccard | All-556 Jaccard |\n|---|---:|---:|---:|---:|---:|---:|\n'+''.join(f"| {r['Leiden_cluster']} | {r['n_loci']} | {r['n_genome_units']} | {100*float(r['full15_fraction']):.1f}% | {float(r['P10_median']):.1f} | {float(r['seed_median_best_Jaccard']):.3f} | {float(r['all556_best_Jaccard']):.3f} |\n" for r in communities)
contactrows=[]
for cluster in sorted({r['Leiden_cluster'] for r in communities}):
 for name in ['global_joint','saccharide_joint']:
  rs=[r for r in contacts if r['analysis']==name and r['Leiden_cluster']==cluster];values=sorted(float(r['reference_neighbor_fraction']) for r in rs);n=len(values);med=(values[(n-1)//2]+values[n//2])/2
  contactrows.append(dict(Leiden_cluster=cluster,comparison=name,n_loci=n,zero_reference_neighbor_loci=sum(r['zero_reference_neighbor']=='True' for r in rs),median_reference_neighbor_fraction=med))
c.write('tables/community_reference_graph_contacts_summary.csv',contactrows)
report=f'''# Initial 320D viral-locus clustering report

The prespecified Leiden baseline detects **five exploratory communities** ({sizes}), but the combined evidence does **not support treating this five-way partition as a reproducible set of distinct biological groups**. No community passes every prespecified operational stability flag. The primary graph is connected, HDBSCAN mostly rejects discrete density groups, and inclusion of unselected loci changes community membership.

{summarytable}

## Stability and algorithm comparison

At k=20 and resolution=1.0, pairwise seed ARI across seeds 0-4 ranges {float(seed['pairwise_seed_ARI_min']):.3f}-{float(seed['pairwise_seed_ARI_max']):.3f}, median {float(seed['pairwise_seed_ARI_median']):.3f}. Across the prespecified grid, the primary cohort yields {minc}-{maxc} communities. VC03 and VC04 have the strongest same-setting seed median best-Jaccard values (0.923 and 0.971), but both change substantially in the all-556 sensitivity. VC01, VC02 and VC05 have lower same-setting seed stability. Higher stability at a different setting was not used to replace the prespecified baseline.

Baseline HDBSCAN (min_cluster_size=15, min_samples=5, sklearn convention including self) finds **zero clusters and 100% noise**. Five of six HDBSCAN settings assign all points to noise. The mcs10/ms5 setting finds three small groups of sizes 16, 14 and 12, with **245/287 (85.4%) noise**. Its ARI to Leiden is {float(small['ARI_to_Leiden_including_noise']):.4f} when noise is a category. The nonnoise-only ARI is {float(small['ARI_to_Leiden_nonnoise']):.3f}, but covers only 42 loci (14.6%): it is not broad algorithmic agreement. HDBSCAN membership strengths are not probabilities of biochemical function.

## Sensitivities

- **All 556:** ARI on the 287 shared selected loci is {float(bycomp['all556']['ARI']):.3f}. VC01 remains intact but acquires other selected members (best Jaccard 0.840); VC02 is partly retained (0.662). VC03, VC04 and VC05 change more (0.526, 0.306 and 0.553). The selection-conditioned five-way partition does not persist as a whole.
- **Full 15 genes:** retaining 243 selected loci, excluding all 44 shorter representatives, gives ARI {float(bycomp['selected_full15']['ARI']):.3f}; per-community best Jaccard is 0.861-1.000. Short windows do not by themselves explain the principal baseline pattern. Extreme short-input/hubness examples are explicitly listed in tables/extreme_short_window_QC.csv.
- **Duplicate collapse:** 287 loci contain 263 distinct vectors/ordered contexts. The two exact-duplicate equivalence relations coincide. Thirty-eight selected loci belong to 14 repeated-context sets; collapsing and mapping assignments back gives ARI {float(bycomp['selected_deduplicated']['ARI']):.3f}, with Jaccard 0.763-0.947. Duplicates influence details but do not account for all structure; biological loci were retained separately in the primary analysis.
- **Genome/relatedness:** communities contain 35-63 distinct genome units and 30-58 pre-existing sketch components. No single genome contributes more than 7.3% of a community, and the largest relatedness component contributes at most 14.7%. Genome-balanced and existing k21-sketch-component-balanced graph weights give ARIs {float(bycomp['genome_balanced']['ARI']):.3f} and {float(bycomp['genome_component_balanced']['ARI']):.3f}. These controls change some boundaries; they do not establish independence of related contexts. Sketch components are not ANI/species definitions.

## Graph and distance QC

The baseline union kNN graph has {graph['n_edges']} edges, density {float(graph['density']):.4f}, **one connected component** and no isolated nodes. Weighted modularity is {float(next(r for r in grid if r['run_id']=='k20_r1_s42')['weighted_modularity_gamma1']):.4f}; silhouette in original cosine space is {results['selected287']['silhouette']:.3f}. A positive modularity and modest silhouette do not establish a unique natural partition. Maximum directed kNN in-degree is {graph['max_knn_indegree']} (mean k=20); in-degree Gini is {float(graph['indegree_gini']):.3f}. Exact-duplicate edges are {graph['duplicate_pair_edges']}/{graph['n_edges']} ({100*float(graph['duplicate_edge_fraction']):.2f}%). Twenty-nine nodes have exact ties at the k=20 boundary, handled by stable IDs. Full similarity, nearest-distance and k-th-distance distributions are saved; absolute cosine magnitude is not interpreted as biological separation.

## Global and Saccharide reference context

The global joint dataset (2,789 samples) has 10 baseline communities; six contain viral loci. Its labels on selected viral loci have ARI {float(bycomp['global_joint']['ARI']):.3f} with viral-only labels. GC08 combines most VC01, all VC05 and part of VC04, with 147 viral loci and 105 references (100 carrying Saccharide). GC05 contains 107 viral loci and 53 references, mixing VC02/VC03 with a VC01 locus. Viral-only and global partitions therefore differ substantially.

The deliberately Saccharide-focused joint dataset (482 samples) has seven baseline communities; six mix viral loci with Saccharide references. ARI on shared viral loci is {float(bycomp['saccharide_joint']['ARI']):.3f}. VC03 falls entirely into SC03, VC04 mostly into SC07, and VC05 divides across several Saccharide-focused communities. VC01 and VC05 are comparatively close to Saccharide references (median nearest-Saccharide cosine 0.982 and 0.991); VC02-VC04 have lower medians (0.953-0.959). These are relative descriptive comparisons, not calibrated distance thresholds.

The baseline joint graph's per-community reference contacts and zero-reference-neighbor counts are in tables/community_reference_graph_contacts_summary.csv. Reference-poor local neighborhoods can occur even when the containing joint community includes references; counts are explicitly tied to the chosen kNN graph, not inferred from UMAP islands. Hybrid reference classes are retained as overlapping memberships, so marginal class fractions may sum above one. The Saccharide-only comparison is not an enrichment test.

## Scope and interpretation

Candidate ascertainment already used the 11/15-gene classifier and Saccharide P@10, and the reference BGCs were used in supervised training. Functional labels, P@10, taxonomy, annotations and UMAP positions were not clustering features or parameter-selection criteria. Agreement with Saccharide references and algorithmic agreement are not independent functional validation. The robust result is the **sensitivity pattern**, not a validated five-family classification. Inspect these communities together before selecting any biological follow-up.

No new embeddings, BGC-Prophet inference, UMAP fits, functional annotation, pathway assignment or representative-neighborhood selection was performed. The original UMAP analysis and archive remain unchanged.
'''
(c.OUT/'REPORT.md').write_text(report,encoding='utf-8')
# Parameter evidence always points to pre-run configuration AND actual execution outputs.
evidence=[]
def add(parameter,value,source,kind='EXECUTION-VERIFIED',qualification='No additional qualification'):
 evidence.append(dict(evidence_id=f'CE{len(evidence)+1:03d}',parameter_or_decision=parameter,observed_value=value if isinstance(value,str) else json.dumps(value,sort_keys=True),evidence_type=kind,exact_execution_evidence=source,qualification=qualification))
prefix='runs/run_001_science/before.json';post='runs/run_001_science/after.json'
for key in ['command','cwd','python','executable','platform','configuration_sha256']:
 add(key,before[key],prefix+' $.'+key)
for key,value in before['software_versions'].items():add('distribution version '+key,value,prefix+' $.software_versions.'+key+'; runs/run_001_science/dependency_file_hashes.csv')
for key in ['baseline','grid','graph','leiden','hdbscan_baseline','hdbscan_grid','hdbscan','genome_sensitivity','descriptive_stability_flags','parameter_selection']:
 add(key,before['configuration'][key],prefix+' $.configuration.'+key+'; assignments/ and provenance/* actual executed records')
for name,result in results.items():
 add(name+' actual sample/community counts',dict(n=result['n_samples'],communities=result['n_communities'],sizes=result['community_sizes']),'RESULTS.json $.'+name+'; assignments/'+name+'_leiden_grid.csv' if name not in ['genome_balanced','genome_component_balanced'] else 'RESULTS.json $.'+name+'; assignments/'+name+'_leiden_seeds.csv')
add('Primary representative/UMAP alignment','556 saved medoids; 287 original candidates; IDs/window IDs unchanged','PREFLIGHT.json; provenance/scientific_execution.json; VALIDATION.json')
add('Raw input normalization','float64 L2 normalization for cosine only; no PCA/standardization','runs/run_001_science/code/engine.py distances(); config.json graph.normalization')
add('Duplicate handling','All primary loci retained; exact-vector pairs use distance zero; 263 unique selected contexts in separate mapped-back sensitivity','tables/duplicate_context_registry.csv; provenance/context_summary.json; run_001 science code/analyze.py')
add('Nearest-neighbor search','Exhaustive float64 pairwise cosine; lexsort distance then stable sample ID; no self-edges','runs/run_001_science/code/engine.py graph(); graphs/selected287_directed_k20_neighbors.csv')
add('Graph symmetrization and weights','Union, not mutual kNN; max(0,cosine); omit zero weights; no self-loops','config.json graph; graphs/*_k20_edges.csv; VALIDATION.json independent graph/objective checks')
add('Leiden objective','RBConfigurationVertexPartition; weighted configuration null model; gamma=resolution','runs/run_001_science/code/engine.py leiden(); provenance/*_optimizer_iterations.json; VALIDATION.json independent_RB_configuration_objective')
add('Leiden actual iterations','Recorded per fit; one optimizer iteration per call, stop improvement <=1e-12, cap 100','provenance/*_optimizer_iterations.json; provenance/genome*_iterations.json; stability/*_grid_summary.csv')
add('HDBSCAN min_samples convention','sklearn includes query point; differs from contrib implementation by one','runs/run_001_science/before.json $.api.HDBSCAN_signature; provenance/*HDBSCAN_parameters.json; https://scikit-learn.org/stable/modules/generated/sklearn.cluster.HDBSCAN.html')
add('HDBSCAN actual baseline','0 clusters; 287/287 noise; all membership strengths zero','stability/selected287_HDBSCAN_summary.csv mcs15_ms5; assignments/selected287_HDBSCAN_grid.csv')
add('Reference annotations','Used post-fit only; Saccharide includes any hybrid membership; no supplied y','runs/run_001_science/code/analyze.py and contextualize.py; provenance/context_summary.json')
add('Reference P@10 context','287 representative-specific original-space values reproduce prior metadata; never max across windows','tables/selected_reference_neighbors_top20.csv; provenance/context_summary.json reference_P10_and_nearest_identity_reproduced')
add('ARI and Jaccard','Label-invariant shared-locus comparison; deduplicated assignments expanded to original loci; best Jaccard allows many-to-one matching','tables/cross_cohort_partition_agreement.csv; tables/community_sensitivity_Jaccard.csv; VALIDATION.json manual_ARI_*')
add('All-556 selection sensitivity',bycomp['all556']['ARI'],'tables/cross_cohort_partition_agreement.csv comparison=all556')
add('Full15 selected sensitivity',bycomp['selected_full15']['ARI'],'tables/short_representatives_excluded_from_full15.csv; tables/cross_cohort_partition_agreement.csv comparison=selected_full15')
add('Duplicate collapse sensitivity',bycomp['selected_deduplicated']['ARI'],'tables/duplicate_context_registry.csv; cross_cohort_partition_agreement.csv selected_deduplicated')
add('UMAP coordinates','Original joint_umap_samples.csv, unchanged; no fitting here','INPUT_MANIFEST.csv; tables/per_locus_selected287_clustering.csv UMAP1,UMAP2; provenance/figure_exports.json')
add('Initial figure-wrapper integrity warning','Child exit 0; wrapper failed after unused validator file changed during render','runs/run_002_figures/after.json executable_sources_unchanged=false; preserved code/ and figure_exports/; logs/codex_command_records.jsonl','EXECUTION-VERIFIED','Final clean re-render is recorded separately; original run is not rewritten.')
add('First validation failure','Child exit 1; validator syntax error before checks','runs/run_003_validation/after.json; stderr.txt; code/validate_extended.py')
add('Successful scientific exit and integrity',dict(exit_status=after['exit_status'],inputs_unchanged=after['inputs_unchanged'],sources_unchanged=after['executable_sources_unchanged']),post)
add('Original UMAP effective optimization epochs','NOT RECORDED','Previous UMAP METHODS_RECORD audit; no original UMAP refit or historical default inference in this task','NOT RECORDED')
add('Original prior-inference raw execution transcript','NOT RECORDED in this clustering run','Archived earlier UMAP record distinguishes stored prior-inference proofs from raw execution evidence','NOT RECORDED')
add('Community functional identities','NOT RECORDED; not inferred or assigned','Scope and stopping point; neutral VC/GC/SC labels only','NOT RECORDED')
add('Robustness interpretation','No five-way partition passes every prespecified core flag; HDBSCAN mostly noise; inference remains exploratory','tables/community_summary.csv; REPORT.md; Figures C1-C4','RETROSPECTIVE-ONLY','Interpretation of executed results, not independent biochemical validation.')
c.write('EXECUTION_EVIDENCE_MATRIX.csv',evidence)
# Also emit the requested evidence matrix as TSV.
with (c.OUT/'EXECUTION_EVIDENCE_MATRIX.tsv').open('w',encoding='utf-8',newline='') as f:
 import csv
 w=csv.DictWriter(f,list(evidence[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(evidence)
methods=f'''# Methods record: executed original-320D clustering

This record is based on **pre-execution records and saved outputs**, not inferred current defaults. Source snapshots, configuration, input hashes, interpreter and versions were recorded before each scientific subprocess. See EXECUTION_EVIDENCE_MATRIX.tsv and EXECUTION_LOG.md for exact evidence locations.

## Datasets and ascertainment

We reused the archived 556x320 observed 15-gene-workflow medoids, preserving locus/window identifiers and row alignment. Primary clustering used only the 287 previously selected evaluable loci. Candidate status was not recalculated; selection was the earlier classifier/P@10 union across 11/15-gene windows. This cohort restriction creates ascertainment bias. Existing compatible CORE_ONLY reference embeddings contain 2,502 unique BGCs, including 195 with any Saccharide annotation. Original hashes were matched to archival commit 7b35093eefa99a4704e7b64039f51b7564788fd4. Prior model inference and the original UMAP were not rerun.

## Cosine graph

The executed engine uses float64 L2-normalized vectors and exhaustive pairwise dot products. Cosine similarity is clipped to [-1,1]; distances are 1-similarity, with self and exact-vector duplicate distances set to zero. No feature scaling, PCA or whitening is applied. For each node, k other nodes are selected by ascending distance, breaking exact ties by stable sample ID. The undirected graph is the union of directed neighbor relations, with weight max(0,cosine); zero weights and self-loops are omitted. Every biological locus is retained in the primary graph, including exact duplicates. Stable-ID tie handling can influence duplicate/hub edges; it is explicitly audited. The baseline graph has {graph['n_edges']} edges, density {float(graph['density']):.6f}, one component and no isolates.

## Leiden and stability

The executed objective is leidenalg.RBConfigurationVertexPartition with singleton initialization, no maximum community-size cap, configured resolution and explicit RNG seed. One optimizer iteration is called repeatedly until improvement <=1e-12, capped at 100 calls; the actual counts/improvements and convergence flag are saved for every fit. The baseline is k=20, resolution=1.0, seed=42. The grid crosses k in 10/15/20/30, resolution in 0.5/1.0/1.5/2.0 and seeds 0-4, giving 80 grid fits plus the baseline for each of six cohorts. Parameter choice did not use annotations, P@10, desired counts or UMAP appearance. Neutral cluster identifiers are ordered by lexicographically smallest member ID, not by function.

We compare partitions with ARI and best community Jaccard on shared loci; matching is label-invariant and may be many-to-one. Seed pairwise ARIs, all-grid primary ARIs and coassignment frequencies are retained. Operational flags were fixed before fitting: seed median Jaccard >=0.8, grid median >=0.6, and all-556/full15/deduplication best-Jaccard >=0.6, with at least five retained full15 members. They are descriptive, uncalibrated flags and do not include a biochemical validation claim or select parameters. None of the baseline communities passes all flags.

## HDBSCAN

The independent comparison uses sklearn.cluster.HDBSCAN on a dense precomputed cosine-distance matrix with algorithm=brute, alpha=1, EOM cluster selection, epsilon=0, allow_single_cluster=false, n_jobs=1 and copy=true. The original distance matrix is checked unchanged. We cross min_cluster_size=10/15/25 with min_samples=5/10; the fixed reporting baseline is 15/5. This sklearn min_samples counts the query point itself, unlike the contrib package's convention. All actual constructor values and memberships/probabilities are saved. Noise is reported explicitly. ARI with noise as one category is distinguished from nonnoise-only ARI and coverage; membership strengths are not functional probabilities.

## Global and Saccharide reference comparisons

We preserve viral-only labels when summarizing original-space nearest-reference identities and marginal class memberships, retaining hybrids. Separate label-free joint Leiden grids use 2,789 samples (287 viral plus all references) and 482 samples (287 viral plus 195 Saccharide references). Global GC and focused SC IDs remain distinct from viral VC IDs. Tables report shared-viral partition correspondence, joint reference/viral ratios and graph reference contacts. P@10 remains that of each saved representative; all 287 prior values and nearest IDs are reproduced. The Saccharide-restricted comparison is deliberately enriched and is not an enrichment test.

## Sensitivity and confounders

The identical primary grid is repeated on all 556 evaluable loci, on the 243 selected full15 representatives, and on 263 unique selected vectors/ordered contexts. All 44 excluded short representatives are listed. Exact-vector and ordered-context equivalence coincide here; the collapse retains one deterministic observed vector per equivalence class and maps labels back to all original loci, without averaging. Additional baseline/seed controls retain all 287 loci and divide edge weights by sqrt(block_multiplicity_i*block_multiplicity_j), using genome-unit and pre-existing k21-sketch Jaccard-0.2 component blocks separately. These are confounder controls, not functional features. The sketch components are coarse proxies, not ANI/species or independent biological replicates. We report genome diversity, repeated contexts, hubness, exact ties, components, distance distributions and extreme short-window statistics.

## Visualization and validation

C1 and C2 reuse the original joint UMAP coordinates. C3 shows shared-viral memberships across separate partitions; C4 shows sizes/stability/sensitivity. UMAP positions do not enter clustering. Final exports are PNG at 600 dpi and SVG with vector marks/color bars. Independent validation uses SciPy cosine distances, graph/objective arithmetic, manual contingency-based ARI, metadata/coordinate equality and export structure; it does not rerun the clustering. Input integrity and code equality are checked after each recorded subprocess. Initial unsuccessful setup, figure-wrapper and validation attempts remain preserved rather than replaced.

## Software and execution

The successful scientific command, working directory, exact executable source SHA256 and code snapshots are in runs/run_001_science/before.json. Python is {before['python'].splitlines()[0]}. Distribution versions were recorded at execution: {json.dumps(before['software_versions'],sort_keys=True)}. Full relevant Python/native dependency-file fingerprints are retained for the scientific run. Postprocessing records its actual loaded entrypoints, native files and versions and links the primary fingerprint. Explicit numerical-library thread variables are one. Complete command stdout/stderr, child exit status, inputs-before/after and outputs-after hashes are retained. No current default is substituted for an unrecorded historical parameter.

## Nature Microbiology Methods draft

Previously generated 320-dimensional BGC-Prophet CORE_ONLY embeddings were reused to examine representation-space structure among viral glycosyltransferase-associated loci. Primary analyses included 287 previously selected evaluable loci, each represented by its existing observed 15-gene-workflow cosine medoid; candidate selection and locus representatives were not recalculated. L2-normalized embeddings were compared by cosine similarity without PCA or feature standardization. An undirected graph was constructed from the union of directed k-nearest-neighbor relations, excluding self-edges and assigning positive cosine similarities as edge weights. Exact ties were resolved by stable sample identifier and duplicate loci were retained.

Leiden optimization used the weighted RB configuration-model objective, with a prespecified baseline of k=20, resolution=1.0 and seed=42. Stability was assessed across k=10,15,20,30, resolution=0.5,1.0,1.5,2.0 and seeds 0-4. Optimizer passes continued until improvement <=1e-12 or 100 passes, with actual iteration counts recorded. HDBSCAN was applied independently to precomputed cosine distances, using EOM selection and min_cluster_size=10,15,25 crossed with min_samples=5,10 under the scikit-learn self-inclusive convention. Noise, membership strengths and agreement with Leiden were reported without forcing concordance.

Reference context was characterized by original-space nearest neighbors and separate label-free joint partitions with all 2,502 reference BGCs or the 195 Saccharide-annotated references, retaining hybrid annotations. Sensitivities included all 556 evaluable viral loci, full15 representatives, exact-context collapse and graph-weight controls for genome-unit and pre-existing sequence-sketch-component multiplicity. Partitions were compared on shared loci using adjusted Rand indices and community Jaccard overlap. The original UMAP was reused only for display. Input/code/configuration identities, software versions, execution logs and output checksums were retained. Community structure and reference agreement were treated as exploratory properties of a selected, supervised learned representation, not as independent functional validation.

## NOT RECORDED and limits

The effective optimization epoch count of the **old UMAP** remains NOT RECORDED; it is irrelevant to new 320D clustering and was not reconstructed from installed defaults. A raw original prior-inference execution transcript is NOT RECORDED in this clustering execution; saved source/embedding identities and prior archival proofs are used. No functional identity, pathway identity or new annotation is recorded or inferred. The original archive's historical limitations remain in force. Current clustering iteration counts, parameters and versions are recorded contemporaneously.
'''
(c.OUT/'METHODS_RECORD.md').write_text(methods,encoding='utf-8')
# Index all retained executions; distinguish child success from wrapper guard failure.
log=['# Execution log index\n\nEvery recorded scientific/postprocessing attempt retains before.json, code snapshots, stdout.txt, stderr.txt, after.json and output manifests. Dependency setup attempts are retained under logs/.\n\n| Run | Stage | Child exit | Input unchanged | Recorded-source unchanged | Interpretation |\n|---|---|---:|---|---|---|\n']
for run in runs:
 p=run/'after.json'
 if not p.exists():continue
 r=json.loads(p.read_text());interpretation='Successful recorded subprocess' if r['exit_status']==0 and r['executable_sources_unchanged'] else 'FAILED child' if r['exit_status'] else 'Child completed; wrapper failed integrity guard after unused validator edit'
 log.append(f"| {run.name} | {r['stage']} | {r['exit_status']} | {r['inputs_unchanged']} | {r['executable_sources_unchanged']} | {interpretation} |\n")
log.append('\nRun 001 science is the only clustering-fitting execution. The initial figure child completed, but the wrapper detected an edit to an unused validator; its source snapshots and initial figure assets remain in run_002_figures. Run 003 failed on validator syntax before checks. A corrected validator and clean final figure render are separate subsequent executions. Parent-wrapper traceback evidence is preserved in logs/codex_command_records.jsonl. This index does not relabel the unsuccessful attempts as successful.\n\nEvery record and log has a SHA256 in OUTPUT_MANIFEST.tsv / archive manifest. The pre-execution records, not this later index or current source, establish executed configuration and versions. Raw local session metadata/private system material is not redistributed.\n')
(c.OUT/'EXECUTION_LOG.md').write_text(''.join(log),encoding='utf-8')
print('Methods,  evidence matrix, execution index and report created from retained records; no fitting.')
