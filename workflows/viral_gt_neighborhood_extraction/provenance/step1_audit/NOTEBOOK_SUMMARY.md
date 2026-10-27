# Step 1 — verified notebook summary

**Audited 26 October 2026; no extraction or inference rerun.** The input was the original WP01B largest coarse structural cluster, represented by **A0A3Q8Q3U2**, drawn from the WP01A MGAT2/GT-A–BFVD viral-structure discovery workflow. Headerless membership column 2 supplied **594 distinct accessions** (592 UniProt and 2 UniParc), including the representative; this was not the later spike-in cohort. Structural selection is not proof that every protein is a GT.

The frozen run retained all 594 requested anchors: **593 mapped**, while **UPI001FBC14AE** remained unresolved because only a tagged PDB construct lacked an authoritative genomic crosslink. The 593 mapped accessions produced **757 occurrence-specific windows**, including all qualifying source occurrences and aliases. Every annotated CDS counted in ascending genomic-coordinate order, regardless of strand. The anchor occupied relative position 0; nominal 15-gene and 11-gene contexts used **±7 and ±5 CDS**, respectively. Record ends were neither wrapped nor compensated by shifting, padding or repetition. Invalid translations remained in the extracted tables and excluded model inference rather than changing CDS order.

The 15-gene run produced **616 complete + 141 shortened windows** (10,688 CDS rows); the later 11-gene experiment retained the same 757 occurrences, producing **658 complete + 99 shortened windows** (7,967 CDS rows). Downstream grouping gave **565 baseline overlap loci** and **466 nucleotide-hash-plus-taxid genome units** from 515 source accessions. Native 11-gene regrouping gives **575 loci**; paired comparisons retained the 565 baseline identities. Inference is a later stage: **737/757** 15-gene inputs and **740/757** 11-gene inputs were inferred, with 737 shared pairs and 3 eligibility rescues.

Exact references (all paths below share `/home/rahnn/260908_MiBG/` unless stated):

- Original cohort: `/home/rahnn/projects/viral-gt-research-record/results/WP01B_coarse_clustering/canonical/cluster_membership.tsv`.
- Actual input: `analysis/viral_gt_mining_20261017/inputs/baseline_594_anchors.txt`; matching frozen copy in `runs/baseline_594_v2/input_snapshot/`.
- Exact 15-gene source: `analysis/viral_gt_mining_20261017/runs/baseline_594_v2/implementation_snapshot/src/viral_gt_mining/extract.py`.
- Exact 11-gene source: `analysis/viral_gt_window_size_20261017/scripts/experiment.py` (`prepare_main`, `shortened`), with `common.py`.
- Outputs: `analysis/viral_gt_mining_20261017/runs/{baseline_594_v2,window11_594_v1}/viral_gt_windows.tsv` and `viral_gt_window_genes.tsv`; original mapping is `baseline_594_v2/anchor_mapping.tsv`.

Source/input/code SHA256 identities and file sizes are in the accompanying `reproducibility_manifest.tsv`. The historical base commit is verified, but the extraction implementation was uncommitted; its frozen source hash, not that commit alone, identifies the executed code.
