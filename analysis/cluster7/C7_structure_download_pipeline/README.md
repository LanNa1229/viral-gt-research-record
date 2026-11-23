# Cluster 7 structural-homologue discovery and FoldTree input provenance

This directory documents the `C7_structure_download_pipeline.zip` research package created for Cluster 7 structural-phylogeny preparation.

**Archive status:** The source archive is not stored in this repository as of this README commit. Preserve the original archive and upload it alongside this README to archive the runnable script and source selection tables.

## Contents of source archive

- `C7_structure_download_pipeline/c7_structure_downloader.py`
- `C7_structure_download_pipeline/README.md` (original pipeline instructions)
- `C7_structure_download_pipeline/requirements.txt`
- `C7_structure_download_pipeline/C7_novel_priority_candidates.tsv`
- `C7_structure_download_pipeline/C7_LPS_Oantigen_Foldseek_hits.tsv`
- `C7_structure_download_pipeline/stage1_selection_preview.tsv`
- `C7_structure_download_pipeline/stage2_neighborhood_selection_preview.tsv`
- `C7_structure_download_pipeline/Tree1A_edge_informed_35_ids.txt`
- `C7_structure_download_pipeline/Tree1A_edge_informed_35_manifest.tsv`

## Analysis provenance

Stage 1 used three Cluster 7 structural queries (A0A1V0SAH7, A0A346FKF4, A0A6J5KTT4) to select 15 top MGnify structural matches and eight annotated LPS/O-antigen GT reference proteins. The eight LPS anchors were selected with support from at least two of three queries, minimum Foldseek probability of 0.85, query coverage of 0.85, and target coverage of 0.60.

Stage 2 defined a separate 35-member edge-informed neighbourhood panel. **This Stage-2 panel was not the complete final FoldTree input**: the final 143-protein analysis used 121 structures from the existing Cluster 7 network, plus 15 MGnify hits and eight LPS/O-antigen anchors, deduplicating Q5JBG6.

The later 143-protein merge and PDB-format normalization scripts are separate workflows and are **not** included in this source archive.

## Source archive integrity

- Filename: `C7_structure_download_pipeline.zip`
- SHA-256: `01114ede891b016db381e3059576a0a4a72efd161a38ebdc6e8a289188447ebd`
- Size: 26,679 bytes (check archive on upload)

## Recommended archive path

`analysis/cluster7/C7_structure_download_pipeline/C7_structure_download_pipeline.zip`

Upload the original ZIP without recompressing it if the above checksum is to remain valid.
