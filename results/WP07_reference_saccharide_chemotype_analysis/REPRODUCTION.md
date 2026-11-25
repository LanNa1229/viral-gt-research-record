# Reproduction

Run the commands from this archive directory. Python 3.11 or later provides checksum and data verification with its standard library. Use a scratch directory outside the archive. The replay commands create new files there. The archive manifest covers all archive files except itself.

~~~sh
python -B provenance/reproduce.py archive
python -B provenance/reproduce.py composition ../chemotype_replay
~~~

The composition command downloads the three files in [PUBLIC_DEPENDENCIES.tsv](provenance/PUBLIC_DEPENDENCIES.tsv) anonymously by commit and checks size and SHA-256. It also verifies the official MIBiG 3.1 JSON archive from [EXTERNAL_DEPENDENCIES.tsv](provenance/EXTERNAL_DEPENDENCIES.tsv). It checks every reference accession, original class and community. It derives membership from the original labels, joins the frozen primary decisions and applies the archived seven-group map. It independently reconstructs the 21 counts, fractions, segment limits and accession sets and checks against the historical figure source data. The numerical targets serve as validation assertions.

## 1. Regenerate the figure

Use the recorded NumPy 2.3.5 and Matplotlib 3.11.1 with a licensed copy of Arial Regular. An exact font digest is recorded for the original rendering environment. Replace the example font path with your local file.

~~~sh
python -m pip install numpy==2.3.5 matplotlib==3.11.1
python -B provenance/reproduce.py figure ../chemotype_replay --font /path/to/arial.ttf
~~~

The wrapper runs the archived plotting script from a scratch directory with the frozen source table and configuration. The SVG comparison checks the 21 segment shapes and styles, the labels and page dimensions. The PNG check reports byte equality. PDF and SVG metadata and identifiers can vary by rendering environment or time. The files in [figures/final](figures/final/) are the original exports.

## 2. Recalculate composition

Run the composition command above. Its JSON output gives independently recalculated rows and summaries, including the surface/extracellular counts and fractions. The [archived reference crosswalk](inputs/community_assignments/REFERENCE_COMMUNITY_CROSSWALK.tsv) also lets readers examine the fixed accessions, labels and communities offline.

## 3. Replay the historical reference construction

~~~sh
python -B provenance/reproduce.py reference ../chemotype_replay
~~~

This command adds the official GenBank archive to the dependency cache. To use a local copy, pass its path as `--gbk /path/to/mibig_gbk_3.1.tar.gz`. Use a new scratch directory to replay again. The scratch run extracts annotations from the pinned inputs and executes the historical curation script with its existing source rules, accession decisions and frozen evidence. The wrapper adapts workstation paths and executes the historical reference table construction without the separate viral-neighbor branch. It compares the evidence cards, record manifest, frozen decisions, master and supporting reference tables against the archive with the same contact-string redaction applied.

The replay reconstructs the recorded workflow from the encoded decisions. The original literature reading, model-assisted judgment and the preliminary manual accession assignments cannot be regenerated independently from the saved records. The archive keeps the scientific task specification and decision rationale. It supplies PMID/DOI references and the historical retrieval manifest in place of raw abstracts.
