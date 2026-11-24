# Archive notes

## Files and provenance

The final trees are the three `*_struct_tree.PP.nwk.rooted.final` files. Each metric has one labeled text distance matrix and its matching full-precision NumPy array. The NumPy arrays follow the row and column order of the corresponding text file. `finalset.csv` is FASTA. `allvall_1.csv` is the original headerless, tab-delimited Foldseek export.

The accession manifest records the input cohorts, paths, and SHA-256 hashes of the original selected structures. The 121 local structures match the coordinate members in the restricted [WP02 Zenodo dataset](https://doi.org/10.5281/zenodo.23720756). Its archive and per-file verification instructions are in [WP02 Rights and Retrieval](../../../workflows/WP02_global_SSN_leiden_viral_placement/RIGHTS_AND_RETRIEVAL.md). The other 22 structures match the corresponding historical Stage 1 download records. The discovery selection workflow is documented in the [adjacent discovery package](../C7_structure_download_pipeline/). The overlapping accession Q5JBG6 used the local network structure.

The public package contains the numerical results, sequences, and selected execution records. Coordinate files, normalized or intermediate copies of coordinates, and transient databases are excluded. The accession manifest uses `/home/REDACTED` in place of the personal home-directory component. The execution log uses `REDACTED_HOST` in place of the workstation name. `SOURCE_SHA256.tsv` records the original and published hashes for these and the other retained files. `SHA256SUMS` covers the public files, with the checksum list itself excluded.

## Methods and recorded environment

The run used 16 Snakemake cores. The logged Foldseek job used one thread. The Snakemake workflow requests version 7.8.0 or later. Its archived environment specification requests Python 3.10 and Toytree 2.0.1 and lists QuickTree and Foldseek without package versions. The Foldseek output reports the MMseqs engine identifier `755d3fa84b8259d4d30aa92782f0bfc2bfc3363f`. Exact installed Snakemake and QuickTree versions were not verified for the historical run.

The wrapper records the historical invocation template. It uses the spelling `customstructs=True`. The executed workflow reads `custom_structs`, and the completed log independently confirms `custom_structs: True`. See [Commands](COMMANDS.md) for the observed settings.

The workflow averages the two directional similarities for each pair, takes one minus that average, and applies its truncated Tajima transformation. Its `bfactor` is 0.93 for Fident and 1 for AlnTM and lDDT. The workflow uses QuickTree on the corresponding labeled distance matrix. The included FoldTree workflow and scripts are snapshots from the retained checkout at `d22e4ca08a2903f995ab2bb8eefa4820e226160f`. It contains a locally modified `grabplddt.py`. These snapshots document the retained workflow, without proving the checkout contents on the run date.

The final rule is named `mad_root_post`. Its script reads and writes a tree through Toytree. All three final files are byte-identical to their corresponding postprocessed `.PP.nwk` files in the source run. The retained records do not independently establish that an MAD rerooting occurred.

## Limitations

- Foldseek was run exhaustively with `-e inf` and no separate coverage cutoff. Short and weak alignments are retained. Interpret individual matches alongside alignment length, query and target spans, and structural scores.
- Full input structures were used. Domain extensions and alignment-register differences can affect inferred relationships.
- The three metrics yield different topologies. The archived trees have no internal bootstrap support. Protein structural relationships alone cannot establish organismal phylogeny or the direction of transfer.
- The saved alignment table includes local lDDT scores and alignment endpoints but does not include gapped sequences or a CIGAR string.

## Attribution and terms

- [Moi et al., FoldTree, Nature Structural & Molecular Biology (2025)](https://doi.org/10.1038/s41594-025-01649-8). The included upstream scripts carry the accompanying MIT license.
- Original source coordinates include [AlphaFold DB](https://alphafold.ebi.ac.uk/download) from EMBL-EBI and Google DeepMind, the [ESM Metagenomic Atlas](https://esmatlas.com/about), and [wwPDB entry 2C0N](https://www.rcsb.org/structure/2C0N). The original input selection and source tags are recorded in the manifest.
- The computation also used 23 locally generated AlphaFold 3 models. Cite [Abramson et al., Nature (2024)](https://doi.org/10.1038/s41586-024-07487-w). Numerical and tree outputs were derived by comparing those models with other inputs. To the extent applicable to those derivatives, the [AlphaFold 3 Output Terms](notices/ALPHAFOLD3_OUTPUT_TERMS_OF_USE.md) and [binding notice](notices/Legally%20Binding%20Terms%20of%20Use.txt) accompany the public package.
