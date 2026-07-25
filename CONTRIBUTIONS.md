# Contribution Documentation Standard

## 1. Purpose

This document defines how scientific and technical contributions are recorded across the Viral GT project. The objective is transparent, specific, contemporaneous attribution tied to reproducible evidence.

Contribution documentation must describe what a person actually contributed to a defined research output. It must not be used to infer motives, assign blame, or unilaterally determine authorship.

## 2. CRediT roles used

The project uses the 14-role CRediT taxonomy:

1. Conceptualization
2. Methodology
3. Software
4. Formal analysis
5. Investigation
6. Data curation
7. Validation
8. Visualization
9. Writing – original draft
10. Writing – review and editing
11. Project administration
12. Supervision
13. Resources
14. Funding acquisition

A contributor may hold multiple roles. A role may be shared. When useful, degree may be recorded as `lead`, `equal`, or `supporting`.

CRediT describes contributions and does not independently determine authorship. Authorship decisions must also consider applicable journal, disciplinary, institutional, and collaboration standards.

## 3. Required writing pattern

Each contribution statement should contain five elements:

```text
Contributor + action + scientific object + intellectual or technical purpose + evidence or output
```

Preferred example:

```text
Lan designed the weighted Leiden comparison strategy to test whether Foldseek edge strength improves separation of CAZy-heterogeneous structural communities; the implementation and outputs are recorded in NB-VGT-WP05-20260723-001 and commit <hash>.
```

Insufficient example:

```text
Lan performed analysis.
```

## 4. Technical execution versus intellectual contribution

### 4.1 Technical execution

Technical execution describes implementation of a defined task.

Examples:

- ran Foldseek using an approved configuration;
- converted alignment output into a Cytoscape edge table;
- downloaded AlphaFold models from a defined list;
- generated a phylogenetic tree using a supplied alignment and model;
- formatted a figure using established panel specifications.

Technical execution should be recorded under roles such as Software, Formal analysis, Data curation, Validation, Visualization, or Investigation, depending on the work.

### 4.2 Intellectual contribution

Intellectual contribution changes the scientific framing, analytical design, interpretation, or experimental direction.

Examples:

- formulated the structure-first viral GT discovery strategy;
- selected and justified the structural filtering criteria;
- recognized that the Cytoscape layout did not constitute an objective clustering method;
- designed the weighted Leiden comparison and resolution analysis;
- interpreted Cluster 8 as a DNA-modification-associated lineage;
- connected a multi-GT neighborhood to a sulfated DNA-glycan hypothesis;
- selected the experimental system and designed staged validation;
- formulated manuscript claims and figure architecture.

Intellectual contributions are commonly documented under Conceptualization, Methodology, Formal analysis, Validation, Writing, or Visualization. A single activity may legitimately support more than one role, but the specific contribution must be stated separately for each role.

## 5. Role-specific drafting guidance

### Conceptualization

Record the idea, question, hypothesis, or research objective that was formulated or materially evolved.

Good:

- `Lan conceived the structure-first strategy for recovering viral GTs missed by sequence annotation.`
- `Lan formulated the hypothesis that Cluster 8 proteins function in hydroxymethyl-pyrimidine-associated DNA glycan extension.`

Avoid:

- `Lan had the idea.`

### Methodology

Record design or development of methods, models, criteria, or analytical workflows.

Good:

- `Lan selected min(qTM,tTM) as a conservative global structural-similarity metric and designed the resolution series used for weighted Leiden sensitivity analysis.`

Distinguish from Software: methodology is the scientific design; software is the implementation.

### Software

Record creation, implementation, testing, maintenance, or substantive modification of code.

Good:

- `Lan implemented the reproducible weighted Leiden workflow, including directional-hit collapse, reciprocal-edge audit, viral metadata overrides, and Cytoscape-ready exports.`

Do not imply conceptual ownership merely from coding unless the contributor also designed the method.

### Formal analysis

Record application of computational, statistical, mathematical, or other formal analytical techniques.

Good:

- `Lan quantified community composition, CAZy purity, viral enrichment, and partition stability across edge definitions and Leiden resolutions.`

### Investigation

Record performance of experiments or collection of evidence, including computational investigation when it constitutes research execution.

Good:

- `Lan investigated Cluster 8 gene neighborhoods across viral genomes and classified recurring hydroxymethyl-pyrimidine-associated architectures.`

### Data curation

Record annotation, cleaning, normalization, metadata construction, deduplication, and maintenance enabling reuse.

Good:

- `Lan curated viral origin, CAZy family, host, taxonomy, structure source, and model-quality metadata for the cross-domain structure network.`

### Validation

Record verification, replication, benchmarking, control analysis, or confirmation of reproducibility.

Good:

- `Lan designed and performed positive-control recovery analysis using characterized viral GT structures and audited non-giant clusters for off-target behavior.`

### Visualization

Record conceptual design and generation of figures, panels, network views, or visual encodings.

Good:

- `Lan designed Figure 2B to classify viral GT placement as known-family embedded, CAZy-unassigned embedded, peripheral, bridging, virus-specific dark cluster, or unresolved.`

### Writing – original draft

Record drafting of substantive manuscript or report text.

Good:

- `Lan drafted the initial Results section describing structure-first mining, Cluster 8 neighborhood architectures, and the jumbo-phage sulfated-DNA hypothesis.`

### Writing – review and editing

Record substantive revision, critical interpretation, or rewriting.

Good:

- `Contributor X revised the interpretation of structural phylogeny limitations and clarified the distinction between network placement and evolutionary ancestry.`

### Project administration

Record coordination, milestone management, documentation systems, and integration of work packages.

Good:

- `Lan established the repository, notebook, claim-evidence, and transfer-manifest system and maintained cross-work-package traceability.`

### Supervision

Record scientific oversight, mentoring, or responsibility for a contributor's work.

Good:

- `Contributor X supervised trainee Y's neighborhood annotation and reviewed the resulting classifications.`

Supervision must not be used as a substitute for describing the supervised contributor's actual work.

### Resources

Record provision of materials, data access, computing allocations, instruments, samples, or other resources.

Good:

- `Contributor X provided access to the viral proteome structure collection and associated metadata under the project data-use terms.`

### Funding acquisition

Record acquisition of financial support specifically enabling the work.

Good:

- `Contributor X obtained NSF funding that supported the computational infrastructure and personnel used in this project.`

## 6. Evidence requirements

Each ledger entry should point to one or more of:

- notebook ID;
- Git commit;
- code or configuration file;
- output file;
- dated email;
- meeting summary;
- presentation;
- manuscript version;
- HPC job log or SLURM output;
- file timestamp and checksum;
- data-transfer manifest;
- validation report.

Evidence establishes that a contribution was documented; it does not automatically establish exclusivity or priority. Where priority matters, use the earliest reliable contemporaneous evidence and state any uncertainty.

## 7. Verification status

Use one of:

- `draft-self-recorded`
- `corroborated-by-contemporaneous-record`
- `reviewed-by-contributor`
- `reviewed-by-project-lead`
- `confirmed-by-collaborators`
- `superseded`
- `contested`
- `uncertain`

Do not use `verified` without identifying who verified it and on what evidence.

## 8. Corrections and disagreements

- Never delete a contribution entry solely because an interpretation changed.
- Add a new entry or update the status with a dated note.
- Preserve the original Git history.
- Record factual disagreement neutrally and link both supporting records.
- Separate contribution evidence from authorship conclusions.
- When another person's role is uncertain, record `TBD` rather than guessing.

## 9. Notebook integration

Each notebook must contain a CRediT table with:

| Contributor | CRediT role | Degree | Specific contribution | Evidence | Status |
|---|---|---|---|---|---|

The master ledger is the project-wide index. Notebook entries are the analysis-level record. Both should point to the same Notebook ID, Entry ID, and Git commit.

## 10. Periodic review

Review the master ledger:

- at completion of each work package;
- before each internal presentation;
- before manuscript circulation;
- before submission;
- after major revision;
- when a contributor joins or leaves;
- when an analysis is reassigned or reused in a different manuscript scope.
