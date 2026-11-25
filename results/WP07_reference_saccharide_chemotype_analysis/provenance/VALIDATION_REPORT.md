# Validation report

## Accession and composition checks

The independent utility retrieved the three immutable public files and both official MIBiG 3.1 archives without repository credentials. All matched the historical sizes and SHA-256 values. The 2,502 accession identities, JSON member digests, original metadata labels and saved Leiden assignments agreed with the reference master. The cohort includes 195 unique Saccharide-associated references, including 65 hybrids. The 229 detailed decisions matched the expected priority population and primary chemotype assignments.

| Display category | GC08 | GC05 | Other GCs |
| --- | ---: | ---: | ---: |
| Capsular polysaccharide | 31 | 4 | 1 |
| Extracellular/matrix polysaccharide | 19 | 2 | 1 |
| LPS/LOS/O-antigen | 18 | 3 | 1 |
| S-layer/other surface glycan | 4 | 0 | 0 |
| Aminoglycoside/aminocyclitol | 20 | 7 | 7 |
| Glycosylated specialized metabolite | 4 | 5 | 51 |
| Other/unresolved | 4 | 3 | 10 |
| **Total** | **100** | **24** | **71** |

All 21 figure-source rows matched their independent counts, fractions, percentages, segment limits and exact accession sets. The four-category surface/extracellular fractions were 72/100 (72.0%), 9/24 (37.5%) and 3/71 (4.2%). Historical GC08 totals and the broader 74/100 flag remain explicitly documented in the archive.

## Replay checks

- **Figure.** The original plotting code and font produced PDF, editable SVG and PNG in scratch. The PNG was byte-identical and the SVG page dimensions, 21 segment paths and styles, and 16 text labels matched. The three published figure exports retain their original hashes.
- **Reference preparation and decisions.** The pinned GenBank parser rebuilt the 2,502-record manifest and priority evidence cards. The historical curation rules and recorded evidence regenerated the 229 frozen decisions. JSON comparisons account for line endings and the publication redaction of contact strings.
- **Reference tables.** The historical reference construction regenerated the master, compound, gene, gene-observation, chemical-feature and curation-audit TSVs. All matched the staged publication files after the same contact-string redaction. The historical viral-link branch was outside the replay.
- **Integrity.** The release rehashed the archived originals, screened release copies for private paths, contact details and credentials, and sealed the public copies with SHA-256.

The tests replay established outputs. They do not reproduce the original model-assisted interpretation as an independent judgment exercise or experimentally validate reference or viral products. [Reference replay receipt](REFERENCE_REPLAY.json) · [Reproduction commands](../REPRODUCTION.md).
