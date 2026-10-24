# Section 1 structural mining and coarse clustering

Version section1-v1.0.0. Sole record creator: Lan Na, University of Georgia, https://orcid.org/0000-0003-2788-6123 . Original code: MIT. Original data/documentation: CC BY 4.0. Third-party terms remain in THIRD_PARTY_NOTICES.md.

This is the selected Section 1 snapshot supporting Figure 1 and Supplementary Figures 1–4, not a license or DOI assignment for the wider repository. The 1,727-structure baseline came from the original 447-residue MGAT2 Foldseek search. Missing truncated-core files are non-blocking historical provenance; no independently verified output equivalence is claimed.

Run from the repository root (or the root of the selected software ZIP):

```bash
python -B scripts/reproducible/WP01_section1/validate_section1.py --source-data results/section1_source_data/figure_inputs.json --out /path/to/new_validation/result.json --bootstrap
MPLCONFIGDIR=/path/to/new_plot_cache python -B scripts/reproducible/WP01_section1/reconstruct_figures.py --source-data results/section1_source_data/figure_inputs.json --outdir /path/to/new_figures
```

The default numerical validator needs only Python; --bootstrap requires NumPy. Reconstructed figures require NumPy and Matplotlib. The frozen-input run passed 106 checks, including 24 Mann–Whitney/BH comparisons and three displayed 500-resample intervals (seed 20260712). The prior Source Data verification covered 18 sheets and 79,674 data cells without mismatches.

Historical scientific scripts remain byte-identical in scripts/archive/. ORIGINAL execution-supported sources, ARCHIVED surviving sources, RECONSTRUCTED plotting code and AUDIT-ONLY verification code are distinguished in reports/section1_archive/. Historical absolute defaults in original code are retained as necessary provenance; use isolated working directories and explicit inputs/environment variables. Account-specific prefixes are redacted from public metadata, with changed-file identities recorded in REDACTED_METADATA.tsv.

The public dataset withholds institutional AF3 model/source files and structural databases containing them because redistribution authority was not confirmed. The archived model plate is also withheld conservatively; new plotting code shows a disclosure placeholder for panel B while reproducing numerical A/C. Exact exclusions and hashes are in PUBLIC_EXCLUSIONS.tsv. The complete private rc2 archive remains unchanged. Final classifier source, L230 copy lineage and historical-runtime limitations remain disclosed; no new investigation is required for the nonessential truncated-core issue.

Software DOI: https://doi.org/10.5281/zenodo.23626496 . Dataset DOI: https://doi.org/10.5281/zenodo.23626546 . GitHub release: https://github.com/LanNa1229/viral-gt-research-record/releases/tag/section1-v1.0.0 . Record pages establish actual publication status.
