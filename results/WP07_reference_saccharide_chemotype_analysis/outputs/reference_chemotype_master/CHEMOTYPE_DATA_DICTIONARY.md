# Chemotype data dictionary

The master has 91 fields and exactly 2,502 original biological BGC rows. UTF-8 TSV; source lists/structured annotations use JSON cells. Original labels and original embedding-row strings are retained. Missingness and review states are explicit.

## Missingness and review codes

- NOT_REPORTED: the inspected source/curated evidence does not supply this fact; not an absence claim.
- NOT_REVIEWED: detailed chemical curation has not been performed; never a biological negative.
- NOT_DETERMINED / Unresolved: evidence does not support a precise conclusion or sources conflict.
- True/False surface/extracellular flags describe the documented product scope only, not the absence of genes or all possible functions.
- Human manual verification was not performed. Methods distinguish rule-based source annotation from evidence-guided record/literature curation.

## Master fields

| Field | Definition | Source |
| --- | --- | --- |
| mibig_accession | Original biological reference accession; unique primary key. | frozen CORE_ONLY metadata:BGC |
| original_embedding_row | Original zero-based CORE_ONLY row, retained exactly; augmentation rows are not observations. | frozen metadata:embedding_row |
| original_bgcprophet_class | Original label string used in clustering; never replaced by source chemotypes. | frozen metadata:label |
| original_class_multilabel | Original class tokens in original order as a JSON list; overlapping hybrid classes retained. | frozen metadata:label |
| original_is_hybrid | More than one original class token. | original label tokens |
| global_leiden_community | Existing global joint baseline, k20_r1_s42; no recomputation. | global_joint_leiden_grid.csv:k20_r1_s42 |
| is_saccharide_reference | Any Saccharide token in original labels, including hybrids. | original labels and saved global assignment |
| is_GC08 | GC08 membership of the reference. | original labels and saved global assignment |
| is_GC05 | GC05 membership of the reference. | original labels and saved global assignment |
| is_saccharide_outside_GC08 | Original Saccharide reference outside GC08. | original labels and saved global assignment |
| reference_core_gene_count | Frozen CORE_ONLY reference gene count; distinct from all annotated CDS count. | frozen metadata:core_gene_count |
| is_priority_curation_reference | Union of all Saccharide references and all GC08/GC05 references; 229 unique BGCs. | original labels and saved assignments |
| organism_name | Original MIBiG organism string; no current taxonomy substitution. | JSON:/cluster/organism_name |
| ncbi_taxid | Recorded taxon ID kept as text. | JSON:/cluster/ncbi_tax_id |
| genus | Lexical genus component of recorded organism name where unambiguous; not a new taxonomic determination. | recorded organism string; conservative lexical parsing |
| species | Recorded binomial component where unambiguous; strain text omitted without taxonomy lookup. | recorded organism string; conservative lexical parsing |
| broader_taxonomy | Documented GenBank lineage, retained separately from JSON organism. | frozen GBK:ORGANISM continuation |
| taxonomy_provenance | Distinguishes recorded taxonomy and derived name components from new retrievals. | curation method |
| genbank_organism_name | Frozen GenBank organism qualifier, including any disagreement with JSON. | GBK:source/organism |
| genbank_ncbi_taxid | Taxon cross-references in GenBank source qualifier. | GBK:source/db_xref |
| source_record_accession | Original nucleotide record accession/version, not the BGC accession. | JSON:/cluster/loci/accession |
| bgc_start | Original JSON BGC start coordinate; no coordinate-system conversion or inferred default. | JSON:/cluster/loci/start_coord |
| bgc_end | Original JSON BGC end coordinate; no inferred default. | JSON:/cluster/loci/end_coord |
| bgc_completeness | Original locus completeness, retaining Unknown separately from missing. | JSON:/cluster/loci/completeness |
| genbank_definition | Frozen sequence description, retained for source conflict audit. | GBK:DEFINITION |
| genbank_comment | Original sequence warnings/comments; sequence replacement warnings are retained. | GBK:COMMENT |
| compound_names_original | All original compound names, ordered JSON list; no synonym merging. | JSON:/cluster/compounds/*/compound |
| normalized_product_names | Unicode NFC plus whitespace normalization only; no chemical name reassignment. | original compound names |
| number_of_annotated_compounds | Number of source BGC-compound relationship records, including all congeners. | JSON:/cluster/compounds |
| chemical_formula | Original formula per source compound. | JSON:/cluster/compounds/*/molecular_formula |
| molecular_mass | Original source mol_mass values; no mass computation or unit reinterpretation. | JSON:/cluster/compounds/*/mol_mass |
| SMILES | Original source chemical structure strings; no structure generation. | JSON:/cluster/compounds/*/chem_struct |
| chemical_structure_availability | Availability of recorded structures, not new predictions. | source chem_struct presence |
| source_mibig_biosynthetic_classes | MIBiG 3.1 classes retained independently of frozen training labels. | JSON:/cluster/biosyn_class |
| source_saccharide_subclass | Exact structured source Saccharide subclass; not replaced by curated primary category. | JSON:/cluster/saccharide/subclass |
| product_architecture | Source-supported gross architecture; NOT_REVIEWED is not a negative call. | curation decisions |
| primary_chemotype | Mutually exclusive curated category; independent from original model classes. | CURATION_DECISIONS.json |
| secondary_chemotype_tags | Overlapping descriptors and ambiguity flags; do not add category counts as independent observations. | curation decisions |
| product_glycan_or_glycosylated_metabolite | Distinguishes assembled glycan, glycometabolite, sugar precursor and non-glycosylated final product. | source-supported curation |
| polymeric_or_nonpolymeric | Gross source-supported architecture, not a computed molecular descriptor. | source-supported curation |
| glycan_association_type | Scope-specific relation; final-product chemistry and biosynthetic glycosylation are separate. | curation decisions and cited sources |
| biological_localization_or_target | Only documented glycan location scope; not inferred from a generic biosynthetic label. | curation decisions |
| bacterial_surface_glycan_related | Explicit source-supported curation descriptor. True/False apply to documented product scope; NOT_DETERMINED and NOT_REVIEWED are not biological negatives. | CURATION_DECISIONS.json |
| bacterial_extracellular_glycan_related | Explicit source-supported curation descriptor. True/False apply to documented product scope; NOT_DETERMINED and NOT_REVIEWED are not biological negatives. | CURATION_DECISIONS.json |
| surface_or_extracellular_glycan_combined | Explicit source-supported curation descriptor. True/False apply to documented product scope; NOT_DETERMINED and NOT_REVIEWED are not biological negatives. | CURATION_DECISIONS.json |
| LPS_LOS_O_antigen_combined | Explicit source-supported curation descriptor. True/False apply to documented product scope; NOT_DETERMINED and NOT_REVIEWED are not biological negatives. | CURATION_DECISIONS.json |
| chemotype_annotation_confidence | Explicit source-supported curation descriptor. True/False apply to documented product scope; NOT_DETERMINED and NOT_REVIEWED are not biological negatives. | CURATION_DECISIONS.json |
| chemotype_evidence_source | Explicit source-supported curation descriptor. True/False apply to documented product scope; NOT_DETERMINED and NOT_REVIEWED are not biological negatives. | CURATION_DECISIONS.json |
| known_monosaccharide_identities | Only explicit curated source claims; missing means not reported/reviewed in inspected sources, not chemical absence. See long evidence table for scope, level and source. | reference_chemical_features.tsv |
| unusual_or_rare_monosaccharides | Only explicit curated source claims; missing means not reported/reviewed in inspected sources, not chemical absence. See long evidence table for scope, level and source. | reference_chemical_features.tsv |
| amino_sugars | Only explicit curated source claims; missing means not reported/reviewed in inspected sources, not chemical absence. See long evidence table for scope, level and source. | reference_chemical_features.tsv |
| deoxy_sugars | Only explicit curated source claims; missing means not reported/reviewed in inspected sources, not chemical absence. See long evidence table for scope, level and source. | reference_chemical_features.tsv |
| sugar_nucleotide_precursors | Only explicit curated source claims; missing means not reported/reviewed in inspected sources, not chemical absence. See long evidence table for scope, level and source. | reference_chemical_features.tsv |
| reported_glycosidic_linkages | Only explicit curated source claims; missing means not reported/reviewed in inspected sources, not chemical absence. See long evidence table for scope, level and source. | reference_chemical_features.tsv |
| branching_or_polymeric_architecture | Only explicit curated source claims; missing means not reported/reviewed in inspected sources, not chemical absence. See long evidence table for scope, level and source. | reference_chemical_features.tsv |
| methylation | Only explicit curated source claims; missing means not reported/reviewed in inspected sources, not chemical absence. See long evidence table for scope, level and source. | reference_chemical_features.tsv |
| acetylation | Only explicit curated source claims; missing means not reported/reviewed in inspected sources, not chemical absence. See long evidence table for scope, level and source. | reference_chemical_features.tsv |
| sulfation | Only explicit curated source claims; missing means not reported/reviewed in inspected sources, not chemical absence. See long evidence table for scope, level and source. | reference_chemical_features.tsv |
| phosphorylation | Only explicit curated source claims; missing means not reported/reviewed in inspected sources, not chemical absence. See long evidence table for scope, level and source. | reference_chemical_features.tsv |
| other_glycan_tailoring_modifications | Only explicit curated source claims; missing means not reported/reviewed in inspected sources, not chemical absence. See long evidence table for scope, level and source. | reference_chemical_features.tsv |
| chemical_composition_evidence_level | Chemical evidence levels distinguish primary experimental chemistry, literature annotation and review claims. | curated chemical features |
| sugar_specificity_annotations | Recorded enzyme substrate specificities and evidence; not silently treated as demonstrated product composition. | JSON:/cluster/saccharide/glycosyltransferases |
| number_of_annotated_genes | Number of CDS source annotations in the frozen GenBank BGC record; not the count of experimentally characterized enzymes. | frozen GBK:CDS features |
| number_of_explicit_MIBiG_gene_annotations | Number of explicitly curated JSON gene annotations; distinguish from all CDS. | JSON:/cluster/genes/annotations |
| glycosyltransferase_genes | Index of explicit source function descriptions. Roles may overlap; no sequence-based annotation and no absence claim. | reference_biosynthetic_genes.tsv |
| nucleotide_sugar_biosynthesis_genes | Index of explicit source function descriptions. Roles may overlap; no sequence-based annotation and no absence claim. | reference_biosynthetic_genes.tsv |
| sugar_modification_genes | Index of explicit source function descriptions. Roles may overlap; no sequence-based annotation and no absence claim. | reference_biosynthetic_genes.tsv |
| glycan_polymerization_or_assembly_genes | Index of explicit source function descriptions. Roles may overlap; no sequence-based annotation and no absence claim. | reference_biosynthetic_genes.tsv |
| glycan_export_or_transport_genes | Index of explicit source function descriptions. Roles may overlap; no sequence-based annotation and no absence claim. | reference_biosynthetic_genes.tsv |
| other_biosynthetic_genes | Index of explicit source function descriptions. Roles may overlap; no sequence-based annotation and no absence claim. | reference_biosynthetic_genes.tsv |
| experimentally_characterized_enzymes | Only genes with explicit source Activity assay evidence; knock-outs alone are not direct catalytic characterization. | MIBiG gene-level evidence |
| pathway_description | Source descriptions only; no new pathway reconstruction. | frozen JSON and GBK |
| pathway_evidence_level | Exact locus/BGC evidence list; not transferred to every enzyme or product claim. | JSON:/cluster/loci/evidence |
| pathway_source_annotations_json | Lossless source pathway-class annotations including predictions; not newly calculated chemistry. | MIBiG 3.1 class-specific sections |
| mibig_record_status | Original record status. | JSON:/cluster/status |
| mibig_minimal_record | Original minimal-record flag. | JSON:/cluster/minimal |
| experimental_evidence_types | Explicit source experimental/genetic/chemical evidence labels, excluding prediction, homology and inference. Aggregate labels do not imply every enzyme is characterized. | explicit JSON evidence arrays, experimental terms only |
| source_evidence_types_all | Complete original source evidence vocabulary, retaining predictions/inferences separately from experimental_evidence_types. | all explicit JSON evidence arrays |
| publication_PMIDs | Original source PubMed IDs; no replacement by retrieved identifiers. | JSON:/cluster/publications |
| publication_DOIs | Original source DOI strings. | JSON:/cluster/publications |
| publication_identifiers_original | All original publication identifier strings, unchanged. | JSON:/cluster/publications |
| retrieved_publication_DOIs | Supplemental PubMed metadata, kept separate from original DOI annotations. | 473-record PubMed retrieval |
| source_record_URL | Pinned official snapshot URL with archive member; current live records are not substituted. | frozen official MIBiG 3.1 |
| source_record_SHA256 | SHA-256 of exact raw JSON record bytes. | raw snapshot member |
| source_version | Single release used across all reference annotations. | frozen source manifest |
| curated_chemotype_status | Rule-based or evidence-guided status is explicit; no human manual verification is claimed. | curation decision method |
| manual_review_required | Records needing further human/scientific review, including unreviewed background entries. | curation decisions |
| curation_notes | Accession-specific reasoning and unresolved evidence limits. | curation decisions |
| original_vs_mibig_class_disagreement | Exact original/source class token-set disagreement; both remain unchanged. | frozen metadata versus MIBiG 3.1 |
| JSON_GBK_organism_disagreement | Exact source organism string disagreement; may reflect nomenclature or strain differences. | JSON organism versus GBK source organism |
| annotation_conflict_flag | Recorded source or chemical-scope conflict; inspect notes and retained source values. | source cross-checks |

## Supporting tables

- reference_compounds.tsv: one BGC–source-compound relationship per row, keyed by compound_relationship_id. Original compound JSON and all names/structures are retained. The parent BGC chemotype does not assign a chemical category to every congener.
- reference_biosynthetic_genes.tsv: one source-documented gene entity per row, keyed by gene_record_id. GBK CDS rows use source ordinal; explicit unmatched MIBiG gene IDs remain separate. Mapping ambiguity is explicit. Protein, locus and gene aliases are preserved. JSON function arrays and GBK qualifiers retain exact documented functions.
- Gene indexed_function_categories are transparent lexical indexes of existing annotations, not new sequence/domain predictions. Empty categories mean unclassified/not reported. Genes may appear in several master function lists.
- Only explicit gene-level activity-assay evidence sets experimental_enzyme_activity_recorded=True. Genetic evidence, sequence predictions and unspecified source-function annotations retain separate levels.
- reference_biosynthetic_gene_observations.tsv: normalized gene references in class-specific source sections, with exact JSON pointers and lossless annotation fragments. Module/domain or multi-gene evidence is not silently treated as whole-gene catalytic validation.
- chemotype_curation_audit.tsv: all 2,502 records, including 229 priority records with explicit source/evidence/rationale/status and 2,273 NOT_REVIEWED background records.
- reference_chemical_features.tsv: one accession–field–source claim. Values are JSON lists with evidence_level, source, claim_scope and notes. No chemistry is inferred from broad BGC classes.
- viral_287_reference_chemotype_links.tsv: 5,740 original Top-20 associations. original_320D_cosine_similarity is copied verbatim; community labels and P@10 come unchanged from the existing viral master. These associations are not viral product predictions.
- chemotype_group_comparison.tsv: descriptive category counts and percentages; GC05 Saccharide is a subset of non-GC08 Saccharide, not an independent comparison sample.
- GC08_accession_level_audit.tsv: source-based inclusion and notes for every GC08 Saccharide reference.
