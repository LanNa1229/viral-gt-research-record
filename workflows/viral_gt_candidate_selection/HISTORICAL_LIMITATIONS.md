# Historical versus portable reproducibility

Both historical scripts are preserved byte-for-byte and are not modified or executed. build_figure_b.py includes taxonomy/iTOL generation and validates broad earlier artifact manifests. Literal execution also requires the original taxonomy_cache.tsv, five NCBI taxonomy responses and provenance JSON sidecars, generic-source metadata and earlier integrity caches. build_sankeys.py imports legacy_geometry.py and depends on later taxonomy/figure inputs and plotting/font resources. Those dependencies are unnecessary for the portable selection/exclusion computation and are deliberately not newly distributed.

The portable code independently reproduces selection and grouping from public frozen inputs. It reproduces every one of 33 selection/grouping fields for 466 genomes, all 565 locus rows and exact 245/242 memberships. It does not reconstruct the 14 historical taxonomy/display metadata fields or redo later taxonomy reconciliation; those original fields are preserved byte-for-byte in the data master. Thus full historical Figure B execution is not claimed to be self-contained.

The 12 representative model-window reproductions in Step 2 remain distinct from its full frozen inference cohort; this Step 4 release does not perform or imply another model rerun. The 2,502 reference embeddings, 1,477 scored viral inputs, cosine space and >=0.5 retrieval criterion are inherited from Steps 2-3, not recomputed.

GT-selected frequencies are not all-virus prevalence. Overlap grouping and exact accession-alias collapse do not eliminate all related-genome redundancy. Both is not independent replication. Incomplete inputs and six unresolved outcomes are retained as historically defined. Original generic-source taxonomy must not be called accession-specific validation.

The manual 14-gene Cluster 7 pilot (score 0.329664) is separate from the systematic Catovirus 15/11 pair. Hypothetical strict P@10 >0.5 memberships (230/227) are retrospective Step 3 sensitivity results and are not used here.
