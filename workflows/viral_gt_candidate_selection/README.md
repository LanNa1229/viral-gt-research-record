# Viral GT Step 4 — candidate selection (v1.0.0)

Creator: Lan Na (given name Lan; family name Na). Published v1.0.0 on 2026-10-27 following explicit creator approval.

This selection-only package combines the already frozen Step 1–3 records. It reconstructs exact nucleotide-record-SHA256/source-taxid genome units and the original 15-gene overlap components, applies the frozen classifier-positive OR Saccharide P@10 >= 0.5 rule separately to 11/15-gene inputs, and consolidates evidence into 245 candidates among 466 selected genome units. Three explicitly named single-CDS records are excluded only from a separate 242-genome Sankey membership.

Expected: 331/344 positive 15/11 inputs; 388 positive paired occurrences; 288 positive fixed loci; 245/215/6 candidate/negative/unresolved genome units; support 190 Both, 29 11-only, 26 15-only. Sankey support: 187/29/26.

## Original and portable code

- historical/build_figure_b.py and historical/build_sankeys.py are byte-for-byte original archives. They are evidence, not standalone entry points.
- portable/step4.py, portable/fetch_resources.py and validation/run_isolated.py are new standard-library reproduction/verification code.
- Existing Step 1 grouping code and Step 3 paired.py are referenced with immutable paths/checksums, not copied.
- No source sequences, model weights, predictions, embeddings, retrieval databases or old membership snapshots are bundled.

See EXECUTION_GUIDE.md, METHODS.md, HISTORICAL_LIMITATIONS.md and THIRD_PARTY_NOTICES.md. Runtime and documentation resources are pinned in references/public_resources.tsv; transport ZIPs are pinned separately in references/public_archives.tsv.

The companion Step4_candidate_selection_v1.0.0_DATA.zip preserves the original 466-row screening master (including inherited historical source taxonomy) and detailed three-record exclusion evidence. Its master is used as an expected-output comparison after selection is independently computed; it is not used to determine candidate membership.

Existing publications: https://doi.org/10.5281/zenodo.23680510 ; https://doi.org/10.5281/zenodo.23682612 ; https://doi.org/10.5281/zenodo.23683720 .

Predicted association is exploratory prioritization, not verified glycan biosynthesis. Both denotes two related context sizes, not independent experimental replication. No all-virus prevalence inference is justified.

## Published dataset and immutable software snapshot

Dataset: https://doi.org/10.5281/zenodo.23685561 ; public record: https://zenodo.org/records/23685561 . Exact approved 22-file software snapshot: https://github.com/LanNa1229/viral-gt-research-record/commit/a72e62ad7a5efc23b474564a28e8c3c85ed016f7 .

The one public dataset file is Step4_candidate_selection_v1.0.0_DATA.zip (eight internal files; SHA256 66f4a2c97bd4f21237cfb1bb6afcc8eb7f6ad6d0f10623a5e50d544ea7f56b42). The reviewed software ZIP SHA256 is 1ab29faaa419c5b7c96109cd75d7256c46acdb758df75e5f42cb1e046e00fc4d; its individual files are in the exact software commit, without a redundant repository/Zenodo ZIP.

Explicit approval activates MIT only for original software and CC BY 4.0 only for original derived data/documentation; external source rights/terms remain separate. Pending/proposed wording inside immutable preparation snapshots records their original status. This documentation-only follow-up links the DOI and code; it changes no historical script, executable selection code, validation result, original dataset or cohort.

Selection-only reproducibility was verified; complete historical Figure B execution was not independently reproduced. Original 466 screened / 245 candidate / 242 Sankey accounting, three figure-only exclusions, inherited taxonomy/display limitations and missing historical cache dependencies remain as documented. All 26 pre-existing public references and Steps 1–3 remain unchanged.
