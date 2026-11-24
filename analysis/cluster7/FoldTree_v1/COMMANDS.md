# Historical commands and parameters

The completed run log is the primary execution record. The commands below describe the archived computation. The wrapper and normalization script are historical copies.

```bash
foldseek easy-search C7_gated143_foldtree_run/structs/ C7_gated143_foldtree_run/structs/ \
  C7_gated143_foldtree_run/allvall_1.csv C7_gated143_foldtree_run/tmp \
  --format-output 'query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,lddt,lddtfull,alntmscore' \
  --exhaustive-search --alignment-type 2 -e inf --threads 1

quicktree -i m C7_gated143_foldtree_run/foldtree_fastmemat.txt > C7_gated143_foldtree_run/foldtree_struct_tree.nwk
quicktree -i m C7_gated143_foldtree_run/alntmscore_fastmemat.txt > C7_gated143_foldtree_run/alntmscore_struct_tree.nwk
quicktree -i m C7_gated143_foldtree_run/lddt_fastmemat.txt > C7_gated143_foldtree_run/lddt_struct_tree.nwk
```

The effective custom structure setting, filter setting, and configuration are documented in [Archive notes](ARCHIVE_NOTES.md). The [historical wrapper](scripts/run_C7_gated143_FoldTree.sh) records the Snakemake flags, including `--use-conda`, `--rerun-incomplete` and `--latency-wait 60`. A portable release of the coordinate inputs is outside this archive.
