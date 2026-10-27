# Step 4 Methods and exact semantics

The 757 paired mapped anchor occurrences are kept as separate 15- and 11-gene inputs, including overlapping neighborhoods and incomplete inputs with valid frozen inference. No inference or nearest-neighbor search is performed here.

For each successfully inferred window, positivity is the stored Saccharide classifier flag OR saccharide_fraction_top10 >= 0.5. The flag is preserved; raw scores are not rethresholded. Historical classifier labels used score >= max(Other, 0.5). P@10 >= 0.5 is the historical exploratory boundary, not a validated biological classifier.

Genome keys are ("sequence_and_taxid", available nucleotide-record SHA256, nonempty source taxid). Missing hash/taxid keeps an accession-specific key. IDs are vgg_ plus the first 24 hexadecimal SHA256 characters of compact JSON. Identical sequence with different taxids does not merge; names/ANI/species are not grouping criteria. These are source-record units, not guaranteed complete genomes or independent isolates.

The new adapter independently rebuilds the original grouping algorithm from published windows and CDS locations, then checks every identity/member against published frozen groupings. Loci are transitive connected components of overlapping inclusive nucleotide segments within a genome key. Adjacency alone does not merge. Source origin-spanning CDS parts do not fill the intervening genomic gap. All extracted windows participate, so an unscored/nonpositive bridging window can connect a component. Locus IDs hash compact JSON [genome_key, merged_segments].

The same 565 baseline 15-gene loci are used for both conditions. Native 11-gene components (575 in the historical secondary analysis) are not used for final paired consolidation. Multiple anchors and overlapping inputs remain individual evidence rows before aggregation; repeat/alias inputs are not treated as independent genomes.

A unit is candidate if any observed constituent window is positive. It is not_selected only when all its windows have predictions and all are nonpositive; otherwise it is unresolved. Across sizes, either candidate wins; both nonselected are negative; all remaining combinations are unresolved. Both requires candidate evidence in each size at genome level, not necessarily the same occurrence. All 190 observed Both units also have a jointly supported fixed locus.

Source-accession aliases, 678 coordinate-defined occurrences after alias collapse, 757 accession-specific paired occurrence IDs, 565 loci and 466 genome units are distinct accounting levels. Positive unions are 388 occurrences and 288 loci, with 245 genome units (268 source accessions). Four candidate units have partial inference coverage but observed positive support. Missing predictions are never negatives.

The canonical master is the original 466-row master_genome_dataset.tsv; later taxonomy-enriched tables and the 242 figure cohort do not replace it. Its original inherited taxonomy/display fields are retained unchanged and are not revalidated here. The reproduction compares 33 independently generated selection/grouping fields on all 466 rows, all 565 locus rows, the original paired registry, and both published cohort memberships.

Figure-only rule: remove units whose source accession aliases intersect {KC511810.1, KX960781.1, KY189087.1} from the 245 candidate list. Published Step 1 GenBank records independently verify exactly one annotated CDS and one-gene 11/15 inputs. Preserve the three units as candidate/Both in the canonical master. This yields 242 (187 Both, 29 11-only, 26 15-only).

The Catovirus ARF08732.1/KY684083.1 locus is negative at 15 genes and positive at 11, with score 0.024441218003630638 -> 0.8780882358551025 and P@10 0.2 -> 0.8; the genome remains Both through other loci. This is context-sensitive model behavior, not proof of improved biological specificity.
