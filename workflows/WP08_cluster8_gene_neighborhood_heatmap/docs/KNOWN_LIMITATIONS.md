# Historical limitations and separate interpretation correction

## TET-positive loci hidden by historical grouping

A0A6J5M498 has TET_5mYOX=3 and C5_MT=0, but is labeled additional-GT/glycan-associated because the first row-class branch requires both C5-MT and TET, and there is no standalone TET branch. Its `inferred_model` already says hybrid DNA-modification + glycan module. This is a mismatch in the grouping scheme, not evidence that the TET-like score should be removed.

The saved evidence specifically identifies UFOVP405_24, afdb-swissprot rank 14, AF-D2W6T1-F1-model_v6, “Tet-like dioxygenase 1”, probability 1.0, E-value 5.287e−9, query coverage 0.9183, alignment length 257. The per-locus table also contains stronger-ranked TET/JBP structural matches with other naming forms, illustrating incomplete name matching. Structural association does not establish a particular DNA substrate or biochemical activity.

A0A3G1L3B7 likewise has TET_5mYOX=2, C5_MT=0; it is historically grouped under broad nucleotide metabolism. The separate correction applies the same rule consistently to both loci. It adds `TET/5mYOX-associated (C5-MT unsupported)` for TET≥2, C5-MT<2 and dCMP-HMase<2. It changes primary_group and display group_rank only (to 2.5); all scores, query IDs, evidence, inferred models and source paths remain unchanged. The corrected plot uses the recovered plotting code with this added group color/order/legend entry. This is an explicitly proposed interpretation, not a new biological validation or replacement of the historical analysis.

## Broad nucleotide and glycan labels

The historical “ThyX/ThyS / nucleotide-metabolism-associated” group accepts a dCMP/dTMP/pyrimidine-related match or a nucleotidyltransferase/pyrophosphorylase match, even without ThyX or ThyS. All 11 ThyX values are zero; all five rows assigned this broad group also have ThyS=0. The shortened figure legend “ThyX/ThyS” is therefore misleading if interpreted literally. Nucleotide/sugar metabolism is not evidence by itself for DNA glycosylation or hydroxymethyl-base production.

Generic reductase, ligase, acyltransferase, epimerase-like and APS/PAPS annotations are broad. A single hit may populate several classes and affect several totals; database redundancy does not constitute independent replication. TET_associated_GT_like uses a broad `tet.*gt` expression as well as explicit phrases, so this column requires hit-level review before mechanistic interpretation.

## Ambiguous ORFs and evidence limitations

Maximum-per-class aggregation can assign mutually alternative functions to the same ORF from different hits, or count broad domains as whole-protein function. Target coverage is unused in the score, and functional specificity, catalytic residues, substrate identity, structural superposition and genomic orientation/distance are not independently tested by the matrix builder. `validation/ambiguous_multiclass_ORFs.json` inventories ORFs with two or more score≥2 classes; multi-class status is a review flag, not proof of an erroneous assignment. Many valid multi-domain/overlapping annotations exist.

Unknown/hypothetical and function-assigned hits can coexist. Zero is a lack of a qualifying named match, not confirmed absence. Neighborhoods contain differing query counts, and this is a selected 11-locus panel rather than an unbiased prevalence estimate. The heatmap summarizes context, not gene order or a demonstrated pathway.

## Accession and query-ID mismatches

| Candidate | Fixed-manifest focal accession | Query namespace in saved table |
|---|---|---|
| B2ZYD6 | BAG41682.1 | BAG_### |
| A0A5P8D8M7 | QFP95377.1 | SEA_YUUY_## |
| A0A3G1L3B7 | ATW62693.1 | SCBWM1_gp# |
| A0A8F2DCD3 | QWS69713.1 | SEA_WILLIAMSTRONG_## |
| A0A159B7X0 | AKJ73160.1 | AKJ_### |

The builder requires exact equality to recognize the focal Cluster 8 ORF. These five focal IDs do not occur in their tables, and all five have historical cluster8_GT_like=0 despite GT-like neighborhood hits. Their focal GT may instead contribute to additional_GT. Related manifests recover local model/query aliases but do not, by themselves, establish all accession-to-query mappings by sequence. No accession was silently remapped. BAG41682/BAG41686 or any other cross-session identity must not be equated without sequence/source evidence.

## Historical provenance limits

The completed raw combined-hit TSVs, annotation outputs and code are recovered; the exact command line and generating script version for every historical Foldseek submission are not conclusively established. Multiple batch script variants exist; all recovered annotation script copies share one hash, while the June 17 batch differs from later copies. Recovery preserves variants and reports the distinction rather than choosing one as the proven generator. Absolute paths in preserved manifests are historical, not portable paths. API database versions and server state cannot be recreated by rerunning an online search today; no new searches were made. Related manifests sometimes name other loci or earlier templates, so directory/name proximity is not treated as confirmed lineage.

## Rendering and publication

Exact PNG reproduction was verified in the existing Linux environment; Windows Matplotlib 3.9.1 produced a one-pixel width difference. The PDF differs only in CreationDate in the matching Linux environment. The recovered layout has some header/row-label crowding and was preserved. This package validates computational reproducibility with caveats, not biochemical function or publication rights. Creators, licensing and final manuscript figure associations still require author confirmation in release materials.
