# Gene-class assignments

Human-readable audit of all 33 classes. The exact regular expressions and GT exclusion/overlap conditions remain in the recovered `classify_hit` function; this table is explanatory, not replacement code. All rules search the concatenated fields described in Methods.

| Class | Recovered matching logic |
|---|---|
| C5_MT | C5/cytosine/5-methylcytosine methyltransferase text patterns. |
| TET_5mYOX | 5mYOX, 5-methylpyrimidine dioxygenase, tet-like/tet family, 2OG/Fe or 2-oxoglutarate dioxygenase; excludes TET-associated glycosyltransferase text. Does not comprehensively recognize all TET/JBP names. |
| dCMP_hydroxymethylase | dCMP hydroxymethylase, thymidylate synthase/dCMP hydroxymethylase or dUTPase/hydroxymethyl patterns. |
| hmC_hmdU_related | 5hmC/5hmdU, hydroxymethylcytosine/uracil/DNA, hmdUTP/hmdCTP. |
| ThyX | thyx word or flavin-dependent thymidylate synthase. |
| ThyS | thys word or thymidylate synthase; can overlap ThyX and dCMP-HMase. |
| dCMP_dTMP_enzyme | dCMP/dTMP/dUMP/dUTP/dCTP/dTTP, pyrimidine, thymidine, deoxycytidine or deoxyuridine. |
| nucleotidyltransferase | Nucleotidyltransferase, pyrophosphorylase, nucleoside-diphosphate, NTP transferase or sugar-1-phosphate nucleotidyltransferase. |
| DNA_polymerase | DNA polymerase, polymerase III/I, family A/B polymerase. |
| helicase | Helicase, DnaB, UvrD or RecA-like helicase. |
| ligase | Any ligase pattern; not limited to DNA ligase. |
| nuclease | Nuclease, endo/exonuclease, RNase H, restriction endonuclease, HNH. |
| primase | Primase. |
| recombinase | Recombinase, integrase, resolvase or broad rec[a-z] pattern. |
| cluster8_GT_like | Exact query_id == manifest cluster8_query_id, plus GT-specific or TET-associated GT text. |
| TET_associated_GT_like | TET-associated glycosyltransferase phrases or broad tet.*gt regex; independent of exact cluster8 query match. |
| additional_GT | GT-specific text on a non-exact cluster8 query, excluding generic transferase/kinase/ligase/polymerase/NTase terms. |
| NDP_sugar_enzyme | NDP/nucleotide/UDP/GDP/CDP/TDP sugar; UDP-glucose/GlcNAc/galactose, GDP-mannose, dTDP or rml/rfb A–D. |
| epimerase | Epimerase/dehydratase/dehydrogenase/isomerase/mutase AND sugar or nucleotide-name pattern; thus broader than epimerases. |
| dehydratase | Any dehydratase. |
| reductase | Any reductase. |
| aminotransferase | Aminotransferase or transaminase. |
| methyltransferase | Any methyltransferase; overlaps C5_MT. |
| acetyltransferase | Acetyltransferase or acyltransferase. |
| sulfotransferase | Sulfotransferase. |
| APS_PAPS_kinase | APS kinase/adenylyl-sulfate kinase/PAPS/sulfate adenylyltransferase/atps/adenosine phosphosulfate; broader than one kinase. |
| capsid | Capsid or HK97. |
| portal | Portal. |
| tail | Tail, baseplate or fiber. |
| tape_measure | Tape measure/tapemeasure. |
| terminase | Terminase. |
| conserved_hypothetical | Conserved hypothetical, hypothetical protein or uncharacterized protein. |
| DUF | DUF followed by digits or domain of unknown function. |
