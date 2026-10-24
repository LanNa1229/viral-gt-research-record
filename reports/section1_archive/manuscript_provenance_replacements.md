# Proposed manuscript provenance replacements — author review

Original DOCX and merged PDF are unchanged. These replacements implement the author's scope decision without modifying reported downstream values.

**Results, first structural-retrieval paragraph:** “A Foldseek search of viral protein structures using a 447-residue MGAT2 query yielded the 1,727-structure dataset used for all downstream coarse clustering and analyses described here.”

**Methods, Structure-guided retrieval and coarse clustering:** “The fixed input dataset for downstream analyses comprised 1,727 structures retrieved through the original Foldseek search using a 447-residue MGAT2 structure as the query. Structure-file identities were verified against the preserved historical checksum manifest. A separate truncated GT-A core search was performed on another computer; its original files were not available for this archive. No independently verified equivalence between its output and the MGAT2-derived dataset is asserted.”

Use the revised Figure 1 legend in figures/CANDIDATE_CAPTIONS.md. Remove any clause attributing this fixed downstream baseline to the truncated-core search or asserting identical search outputs unless independently supported. The separate search may be omitted from the main workflow description or retained as this qualified historical note. This limitation is non-blocking and does not require additional retrieval or computation.
