# Main Figure 7C | Viral GT communities and reference glycan chemotypes

Main Figure 7C maps viral glycosyltransferase-associated loci relative to chemically annotated biosynthetic gene clusters (BGCs) in the learned BGC-Prophet representation. The original joint UMAP included 2,502 MIBiG reference BGCs and 556 evaluable viral loci. Each panel displays all 2,502 references and the same 287 previously selected viral loci using the original saved coordinates. Among the references, 195 carry the original Saccharide label, including 65 hybrid BGCs. The other 2,307 references provide context.

Viral communities were obtained in the original 320-dimensional embedding space. The workflow L2-normalizes the selected representative vectors and constructs a union-symmetrized 20-nearest-neighbor graph with positive cosine edge weights. Weighted Leiden optimization uses the RB configuration-model objective, resolution 1.0 and seed 42. Community identifiers follow the lexicographically smallest member in each community. VC01 through VC05 contain 68, 55, 74, 35 and 55 loci, respectively.

The left panel highlights seven groups derived from frozen MIBiG 3.1 chemotypes, with viral loci in the background. The right highlights the five viral communities against Saccharide reference landmarks. Stable identifiers link the frozen chemotypes, community assignments and UMAP coordinates.

The dependency manifest and selected-row crosswalk are under `inputs/`. The original clustering engine and historical plotting script are under `scripts/`, and the graph, assignments and verification records are under `outputs/`. The unified table, original PDF, editable SVG, 600-dpi PNG and legend are under `figures/`. Follow [REPRODUCTION.md](REPRODUCTION.md) to check the public inputs, replay the baseline and render the paired figure.

The projection is descriptive. Its geometry and community assignments share the learned representation, and the historical sensitivity analysis found that none of the communities met every operational stability criterion. Spatial correspondence with a reference chemotype does not establish which products the viral loci synthesize.
