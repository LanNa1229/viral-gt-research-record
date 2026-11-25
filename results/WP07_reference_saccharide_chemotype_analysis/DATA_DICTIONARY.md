# Data dictionary

| File or directory | Grain and key | Interpretation |
| --- | --- | --- |
| [reference_2502_chemotype_master.tsv](outputs/reference_chemotype_master/reference_2502_chemotype_master.tsv) | One biological reference per mibig_accession | Source class string and token list, fixed embedding row, community, product, curated fields, evidence and raw JSON digest. |
| [reference_compounds.tsv](outputs/reference_chemotype_master/reference_compounds.tsv) | Compound relationship within mibig_accession | Source product names, formula, mass and structure where supplied. A parent BGC chemotype is distinct from a compound-specific decision. |
| [reference_biosynthetic_genes.tsv](outputs/reference_chemotype_master/reference_biosynthetic_genes.tsv) | Gene record within mibig_accession | Original GenBank and MIBiG product annotations, mapping status, function evidence and source record digests. |
| [reference_biosynthetic_gene_observations.tsv](outputs/reference_chemotype_master/reference_biosynthetic_gene_observations.tsv) | Individual source gene or annotation observation | Retains record-level evidence and source JSON locations. |
| [reference_chemical_features.tsv](outputs/reference_chemotype_master/reference_chemical_features.tsv) | Source-supported feature within mibig_accession | Reported sugar chemistry with evidence fields. |
| [REFERENCE_COMMUNITY_CROSSWALK.tsv](inputs/community_assignments/REFERENCE_COMMUNITY_CROSSWALK.tsv) | One row per reference accession | Exact projection of original class, row index and community from the frozen master. |
| [CURATION_DECISIONS.json](annotation_evidence/curation_records/CURATION_DECISIONS.json) | One decision per priority reference accession | Primary category, tags, evidence source, confidence, review state and rationale. |
| [ACCESSION_LEVEL_CHEMOTYPE_MAPPING.tsv](figures/source_data/ACCESSION_LEVEL_CHEMOTYPE_MAPPING.tsv) | One row per plotted Saccharide accession | Master fields plus disjoint figure population, display group, order and four-category membership. |
| [CHEMOTYPE_GROUP_MAPPING.tsv](figures/source_data/CHEMOTYPE_GROUP_MAPPING.tsv) | One row per historical primary category | Defines the seven display groups and operational surface membership. |
| [FIGURE_SOURCE_DATA.tsv](figures/source_data/FIGURE_SOURCE_DATA.tsv) | One row per population and display category, including zeros | Count, Saccharide denominator, fraction (0–1), percentage (0–100), cumulative percentage limits and exact accession list. |
| [SURFACE_EXTRACELLULAR_SUMMARY.tsv](outputs/community_composition/SURFACE_EXTRACELLULAR_SUMMARY.tsv) | One row per disjoint figure population | Numerator from the first four figure categories and actual Saccharide denominator. |
| [DATA_LINEAGE.tsv](provenance/DATA_LINEAGE.tsv) | One row per archived historical artifact | Source and archive checksums, source path relative to the historical workspace, and release treatment. |

The [full reference data dictionary](outputs/reference_chemotype_master/CHEMOTYPE_DATA_DICTIONARY.md) and [field definitions](outputs/reference_chemotype_master/FIELD_DEFINITIONS.json) document historical fields and missingness. Other/unresolved is a display category spanning several minor assigned groups and unresolved records. The dedicated confidence and review fields retain those distinctions.
