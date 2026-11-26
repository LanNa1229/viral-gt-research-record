# WP01B — Foldseek coarse clustering and Supplementary Figure X-1-1

## Purpose

WP01B documents the permissive Foldseek coarse-clustering step used to
stratify the 1,727 GT-A-query-derived viral structure entries assembled in
WP01A. It also preserves the complete membership network underlying:

> Supplementary Figure X-1-1. Complete membership network of Foldseek
> multi-member clusters recovered from the GT-A structural search.

The record connects the archived clustering script, raw Foldseek assignments,
biological identifier mapping, canonical cluster membership, representative-
centered alignment output, reconstructed Cytoscape source tables, saved
Cytoscape session, and final SVG.

## Scope

WP01B covers:

```text
1,727 WP01A structure entries
  → Foldseek structure database
  → permissive Foldseek clustering
  → final representative–member assignments
  → 761 final clusters
  → 115 multi-member clusters and 646 singletons
  → complete multi-member membership network
  → Cytoscape session and Supplementary Figure X-1-1 SVG
```

The following analyses are outside WP01B:

- functional classification of non-giant representatives;
- pairwise structural-coherence testing of small clusters;
- positive-control spike-in benchmarking;
- singleton length, confidence, rescue, and trimming analyses;
- high-resolution all-versus-all structure-similarity networking.

These belong to the linked downstream work packages.

## Authorship and analytical responsibility

- Analysis conception and design: Lan Na
- Structural dataset curation: Lan Na
- Foldseek clustering execution: Lan Na
- Cluster-membership reconstruction and validation: Lan Na
- Complete membership-network analysis and visualization: Lan Na
- Scientific interpretation and manuscript integration: Lan Na

This repository record documents the analytical provenance and contribution
history of WP01B. Authorship decisions remain subject to applicable journal
and institutional policies.

## Scientific rationale

The initial MGAT2–BFVD search was intentionally permissive and was expected to
recover both GT-A-like structures and false-positive or weakly resolved
structures. Coarse clustering was introduced as a structural stratification
step: recurrent structural signal should concentrate in one or more
multi-member populations, whereas weak, fragmented, unrelated, or highly
divergent entries should be enriched in the long tail and singleton fraction.

The contemporaneous hypothesis recorded by Lan Na was:

> Key hypothesis:  
> Main argument is why I choose large cluster? Why not others?  
> Why I choose coarse clustering? -> Remove false positive.  
> Big, giant cluster = GT-A like protein.  
> Others = false positive.

This language is retained as the original hypothesis. The technically revised
interpretation is narrower: coarse-cluster membership is not a functional
annotation and does not by itself prove that a protein is or is not a GT.
Instead, it provides a reproducible structural partition for downstream
validation. The dominant cluster is treated as a candidate GT-A-enriched core;
non-dominant clusters and singletons remain a mixture of non-GT, divergent
GT-like, fragmented, low-confidence, domain-level, and unresolved entries.

## Input linkage to WP01A

The clustering input was the 1,727-PDB viral structure dataset assembled by
the initial MGAT2 query against BFVD and stored at:

```text
/home/ln72030/virus_pdb/virus_initial_hits
```

The corresponding WP01A retrospective notebook documents the query, retained
Foldseek web export, downloaded PDB count, and unresolved difference between
1,755 BFVD alignment targets and 1,727 downloaded structures.

## Original clustering implementation

Archived source:

```text
scripts/archive/WP01B_coarse_clustering/run_cluster_clean_virus.py
```

Original HPC path:

```text
/home/ln72030/foldseek_script/run_cluster_clean_virus.py
```

Important implementation note: despite its `.py` suffix, this archived file is
a Bash script. The original filename is preserved for provenance.

Archived script SHA-256:

```text
96570c629f220806b11c4908ae9e1aa642971a97c359ee07125d6d5eb436a32a
```

Recorded clustering command:

```bash
foldseek cluster \
  "$DB_NAME" \
  "$CLUSTER_DIR/cluster" \
  "$CLUSTER_TMP" \
  --min-seq-id 0.0 \
  -c 0.05 \
  --cov-mode 0 \
  --threads 64
```

Parameter interpretation:

- `--min-seq-id 0.0` imposed no minimum sequence-identity threshold.
- `-c 0.05 --cov-mode 0` applied a permissive 5% bidirectional coverage
  threshold.
- `--threads 64` controlled parallel execution and not the biological
  criterion.

The Foldseek build observed during retrospective verification was:

```text
19800e921de4dcadd646feac40d83b4114137515
```

This is the build observed in 2026 and is not represented as definitive proof
of the exact executable used in May 2025.

## Membership reconstruction

The archived script prints a path for `cluster_cluster.tsv`, but the surviving
copy does not contain the `foldseek createtsv` command that would have produced
that table. The raw table therefore remains a direct surviving output, while
its exact original export command is unresolved.

The repository preserves:

```text
results/WP01B_coarse_clustering/raw/cluster_cluster.tsv
results/WP01B_coarse_clustering/raw/protein_db.source
results/WP01B_coarse_clustering/canonical/cluster_membership.tsv
```

The reconstruction maps Foldseek numeric identifiers through
`protein_db.source`, converts the raw assignment into biological identifiers,
and verifies that every input entry is assigned once to one final
representative.

Portable validation:

```bash
bash scripts/reproducible/WP01B_coarse_clustering/validate_membership_reconstruction.sh
```

## Verified clustering result

| Quantity | Verified value |
|---|---:|
| Input structure entries | 1,727 |
| Final clusters | 761 |
| Multi-member clusters | 115 |
| Proteins in multi-member clusters | 1,081 |
| Singleton clusters/proteins | 646 |
| Dominant-cluster representative | `A0A3Q8Q3U2` |
| Dominant-cluster size | 594 |
| Second-largest cluster size | 43 |

These values describe the baseline 1,727-entry clustering. They must not be
interchanged with counts from later positive-control spike-in reruns.

## Supplementary Figure X-1-1 network model

The final figure displays all 1,081 proteins in the 115 multi-member clusters.
Singletons are retained in the complete node table for accounting but are not
shown in the exported multi-member network.

### Node roles

| Node role | Count | Definition |
|---|---:|---|
| Representative | 115 | Final Foldseek cluster representative |
| Direct member | 620 | Non-self representative–member alignment is present in `aln_result.m8` |
| Indirect member | 346 | Assigned to the final cluster but no explicit alignment to the final representative is present in the retained m8 |
| Singleton | 646 | One-member final cluster; excluded from the displayed multi-member network |

### Edge roles

| Edge type | Count | Meaning |
|---|---:|---|
| `direct_alignment` | 620 | Explicit representative–member alignment in the retained m8 |
| `cluster_membership` | 346 | Membership-only edge added to reconstruct the complete final cluster |

The complete multi-member network therefore contains 1,081 nodes and 966
edges. Each non-representative multi-member protein contributes exactly one
representative-centered edge:

```text
1,081 nodes - 115 representatives = 966 edges
620 direct-alignment edges + 346 membership-only edges = 966 edges
```

Membership-only edges do not assert an additional pairwise structural
alignment. Layout positions and edge lengths do not encode structural
similarity, evolutionary distance, or alignment strength.

## Figure source files

Representative-centered alignment evidence:

```text
results/WP01B_coarse_clustering/raw/alignment/aln_result.m8
```

Complete network source:

```text
results/WP01B_coarse_clustering/figure_X-1-1/network_source/
├── cytoscape_complete_nodes.tsv
├── cytoscape_complete_edges.tsv
└── cytoscape_edge_style_legend.tsv
```

Network accounting:

```text
results/WP01B_coarse_clustering/figure_X-1-1/accounting/
├── cluster_cytoscape_accounting.tsv
├── cytoscape_absent_cluster_sizes.tsv
├── cytoscape_absent_reps.txt
└── cytoscape_visible_reps.txt
```

Saved visualization:

```text
figures/WP01B_coarse_clustering/figure_X-1-1/
├── 250710_Coarse_Clustering_GTA_CC01_01.cys
├── CC01_01.svg
└── SHA256SUMS
```

## Cytoscape visual encoding

The source tables and saved session encode:

- blue nodes: final representatives;
- orange nodes: direct members;
- purple nodes: indirect members;
- solid green edges: explicit representative–member alignments;
- dashed green edges: membership-only reconstruction edges.

The saved `.cys` session is the authoritative editable visualization record.
The SVG is the archived figure export.

## Integrity anchors

| Artifact | SHA-256 |
|---|---|
| `aln_result.m8` | `46c09ece565970a9e9660e5158c9f4c1446c11844aabdf3c536333c8f1f8bd26` |
| `cytoscape_complete_nodes.tsv` | `42b94616734b1cd62d553fb4c375c7c5a19c49ee17de81dfe25d1f926e787c14` |
| `cytoscape_complete_edges.tsv` | `55ddce7e149b781b8a48f1e8267b463697cb94da3cd998184a7217e8bcd8e5e9` |
| `cytoscape_edge_style_legend.tsv` | `b04937e8cb5dfc327894912b627d2b9fe8ccead55fd7a8915e89077fe4678377` |
| Cytoscape session | `83f81ea299d6ff58c60d2abfab6b5352e05af315b4feeb646211a61ee01d5d87` |
| SVG export | `178fb82bb96b1c3ef12ce8b4e40dd4f6f9c4fdac90e293220bb8c6339c572eb1` |

The complete evidence inventory is maintained in:

```text
results/WP01B_coarse_clustering/evidence_manifest.tsv
```

## Reproducibility status

### Complete

- Raw cluster assignments are preserved.
- Numeric-to-biological identifier mapping is preserved.
- Canonical membership reconstruction is validated.
- The retained representative-centered m8 is preserved.
- Complete node and edge tables are preserved.
- Network accounting files are preserved.
- The editable Cytoscape session and final SVG are preserved.
- Repository and figure artifacts are protected by SHA-256 checksums.

### Incomplete or unresolved

- The exact `foldseek createtsv` invocation that produced the surviving
  `cluster_cluster.tsv` is not present in the archived clustering script.
- The separate code that generated the complete Cytoscape node, edge, and
  accounting tables has not yet been recovered.
- Consequently, artifact-level reproduction of Supplementary Figure X-1-1 is
  complete, but fully automated code-level regeneration of the complete
  network tables remains incomplete.
- The exact May 2025 Foldseek executable build is not independently fixed by
  the surviving evidence.

These limitations do not alter the validated membership counts or the
integrity of the archived figure source files.

## Repository history

| Commit | Record |
|---|---|
| `e3f4890` | Validated WP01B coarse-clustering membership record |
| `622c052` | Figure X-1-1 m8, complete network source, and accounting files |
| `a5a2bef` | Cytoscape session, final SVG, and figure checksums |

## Scientific conclusion and downstream decision

The baseline clustering produced a highly asymmetric structural partition:
one 594-member population was separated from a long tail whose second-largest
cluster contained 43 proteins. This topology was consistent with, but did not
by itself prove, enrichment of recurrent GT-A-like signal in the dominant
cluster.

The dominant cluster was therefore advanced to positive-control recovery and
high-resolution structure-similarity analysis. Non-giant clusters and
singletons were not declared negative; they were retained for independent
functional, structural-coherence, prediction-confidence, and rescue analyses.

