# Historical reference curation task specification

Sections 2–6 are preserved below. These sections describe the intended reference curation. The final decisions and validation records document what the workflow actually produced. SHA-256 of the original task specification is b8a4974d489714fd450c1c9fb640b1fd5e10a43135a79214708df346c0a0ff08.

## 2. REFERENCE POPULATION

Create a master table containing **all 2,502 unique reference BGCs**.

Every row must represent one original biological reference BGC, not an augmentation row.

Include all seven existing BGC-Prophet biosynthetic classes, retaining hybrid annotations.

Identify the following subsets:

**A. All Saccharide references**

Expected: 195 BGCs, including hybrid annotations.

**B. GC08 references**

Expected: 105 BGCs, including 100 Saccharide-annotated references.

**C. GC05 references**

Expected: 53 BGCs, including 24 Saccharide-annotated references.

**D. Saccharide references outside GC08**

Expected: 95 BGCs.

These subsets may overlap. Do not count overlapping records as independent biological observations.

Detailed chemical curation must prioritize:

- All 195 Saccharide references.
- All non-Saccharide references within GC08 and GC05.

The expected combined priority set is 229 unique BGCs.

Automatically extract available source annotations for the full 2,502-reference population.

Use an explicit `NOT_REVIEWED` designation for records not yet subjected to detailed manual or evidence-guided chemotype curation. Do not classify them as biologically negative.

---

## 3. MIBiG DATA ACQUISITION AND PROVENANCE

Use authoritative MIBiG reference records.

First determine whether a frozen MIBiG JSON dataset already exists in the workspace.

If a new MIBiG dataset must be retrieved, obtain it from the official MIBiG database or its official repository.

Record:

- Source URL
- MIBiG version or Git commit SHA
- Retrieval date
- File hashes
- MIBiG accession
- Record status
- Annotation completeness
- Supporting publication identifiers

Use a consistent, pinned source snapshot for the analysis.

Do not silently mix different MIBiG releases.

If the MIBiG annotation differs from the original BGC-Prophet training/reference annotation, preserve both values and explicitly flag the disagreement.

Do not substitute current MIBiG labels for the original frozen reference labels used during model training or clustering.

If a reference accession cannot be matched to a MIBiG record, retain the reference row and mark its annotation as unavailable.

---

## 4. CREATE THE MASTER REFERENCE TABLE

Primary output:

`reference_2502_chemotype_master.tsv`

Exactly one row per reference BGC.

Include the following information.

### A. Reference identity

- mibig_accession
- original_embedding_row
- original_bgcprophet_class
- original_class_multilabel
- global_leiden_community
- is_saccharide_reference
- is_GC08
- is_GC05
- is_saccharide_outside_GC08
- reference_core_gene_count

### B. Source organism

- organism_name
- ncbi_taxid
- genus
- species
- broader_taxonomy, where available
- source_record_accession
- bgc_start
- bgc_end
- bgc_completeness

Distinguish explicitly documented source taxonomy from newly retrieved taxonomy.

### C. Chemical product identity

- compound_names_original
- normalized_product_names
- number_of_annotated_compounds
- chemical_formula, if available
- molecular_mass, if available
- SMILES, if available
- chemical_structure_availability

For multiple compounds, preserve all original names and create a separate compound-level table.

Do not invent chemical structures from product names.

### D. Curated chemical subclass

Create an independent, evidence-supported chemotype classification.

Include:

- product_architecture
- primary_chemotype
- secondary_chemotype_tags
- product_glycan_or_glycosylated_metabolite
- polymeric_or_nonpolymeric
- glycan_association_type
- biological_localization_or_target
- bacterial_surface_glycan_related
- bacterial_extracellular_glycan_related
- surface_or_extracellular_glycan_combined
- chemotype_annotation_confidence
- chemotype_evidence_source

Suggested controlled categories:

1. Capsular polysaccharide
2. Extracellular/matrix polysaccharide
3. Lipopolysaccharide (LPS)
4. Lipooligosaccharide (LOS)
5. O-antigen
6. S-layer glycan
7. Other cell-surface glycan
8. Glycopeptidolipid
9. Teichuronic acid
10. Phosphonoglycan
11. Aminoglycoside / aminocyclitol-related metabolite
12. Glycosylated specialized metabolite
13. Other carbohydrate-related product
14. Non-carbohydrate-associated product
15. Unresolved

Separate LPS, LOS and O-antigen when source evidence permits, while retaining a broader combined category for comparisons.

Do not force ambiguous products into a specific category.

Maintain a mutually exclusive primary category and optional overlapping secondary tags.

### E. Chemical composition and modifications

Where source evidence permits, extract:

- Known monosaccharide identities
- Unusual or rare monosaccharides
- Amino sugars
- Deoxy sugars
- Sugar-nucleotide precursors
- Reported glycosidic linkages
- Branching or polymeric architecture
- Methylation
- Acetylation
- Sulfation
- Phosphorylation
- Other glycan-tailoring modifications

Distinguish:

- Experimentally established chemistry
- Literature-supported chemical annotation
- Genomically predicted chemistry
- Not reported / unresolved

Do not infer linkage stereochemistry, sugar identity or modification chemistry from broad BGC classes alone.

### F. Biosynthetic functions and machinery

Extract available gene-level information from MIBiG and supporting records.

Include:

- number_of_annotated_genes
- glycosyltransferase_genes
- nucleotide_sugar_biosynthesis_genes
- sugar_modification_genes
- glycan_polymerization_or_assembly_genes
- glycan_export_or_transport_genes
- other_biosynthetic_genes
- experimentally_characterized_enzymes
- pathway_description
- pathway_evidence_level

Preserve gene identifiers and documented enzyme functions.

Do not interpret absence of annotation as demonstrated absence of the corresponding enzyme.

Detailed gene-level information should also be saved in a separate long-format table rather than compressed exclusively into one master-cell string.

### G. Evidence and curation

Include:

- mibig_record_status
- mibig_minimal_record
- experimental_evidence_types
- publication_PMIDs
- publication_DOIs
- source_record_URL
- source_record_SHA256
- curated_chemotype_status
- manual_review_required
- curation_notes
- annotation_conflict_flag

Use an explicit evidence hierarchy.

A recorded product identity must not automatically imply that every enzyme in the corresponding BGC has been experimentally characterized.

Preserve uncertainty.

---

## 5. CREATE SUPPORTING NORMALIZED TABLES

In addition to the reference master, generate:

**A. Compound-level table**

`reference_compounds.tsv`

One row per reference BGC–compound relationship, retaining compound names, structural information, original annotations and source evidence.

**B. Gene-level table**

`reference_biosynthetic_genes.tsv`

One row per documented biosynthetic gene or annotated enzyme, retaining gene IDs, functions, enzyme categories, evidence and parent BGC accessions.

**C. Chemotype curation audit**

`chemotype_curation_audit.tsv`

Include:

- MIBiG accession
- Original product designation
- Proposed curated category
- Exact supporting record or publication
- Curation confidence
- Reasoning for ambiguous cases
- Review status

Ensure all 195 Saccharide references and GC08/GC05 exceptions receive an explicit curation status.

Do not claim manual verification where only automated extraction or rule-based classification was performed.

---

## 6. SPECIAL QUALITY CONTROL — GC08 AND GC05

Independently verify the preliminary observations.

### GC08

Expected reference composition:

105 total references.

100 Saccharide-annotated references.

Previous exploratory chemotype counts:

- Capsular polysaccharides: 29
- Extracellular/matrix polysaccharides: 21
- LPS/LOS/O-antigen: 19
- S-layer glycans: 3
- Other categories: 28

Previously reported surface/extracellular glycan-associated total:

72/100 Saccharide references.

These historical numbers are **audit targets, not mandatory outputs**.

Reconstruct the classifications from source records without forcing agreement.

Report any discrepancy, including its accession-level cause.

### GC05

Expected:

53 total references.

24 Saccharide-annotated references.

Curate the actual product categories of all 53 references, including non-Saccharide classes.

### Non-GC08 Saccharide references

Expected:

95 Saccharide-annotated references.

Apply exactly the same chemical classification rules used for GC08.

This dataset will serve as the comparison population in a subsequent chemical-specificity experiment.

Do not perform that enrichment experiment yet.

### Important edge cases

Review entries such as:

- BGC0000727, whose recorded Saccharide label and annotated product require careful interpretation.
- Hybrid Saccharide/Other pathways, which may encode glycosylated specialized metabolites rather than polysaccharide assembly.
- Other-labeled records associated with carbohydrate metabolism.

Preserve original MIBiG labels even when product chemistry suggests a different interpretation.

Do not silently relabel or exclude these records.

---

