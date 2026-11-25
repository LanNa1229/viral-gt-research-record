"""Assemble inspectable result summaries and a concise scientific report."""
from common import *

def table(headers,rows):
    return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+['| '+' | '.join(str(v) for v in r)+' |' for r in rows])
def fmt(x,n=3):return f'{float(x):.{n}f}'
def run():
    pair=read(OUT/'tables/paired_window_comparison.tsv');allrows=read(OUT/'tables/all_occurrence_status.tsv')
    counts=read(OUT/'tables/candidate_counts_by_condition.tsv');recovery=read(OUT/'tables/candidate_recovery_comparison.tsv')
    stats=read(OUT/'tables/genome_block_paired_statistics.tsv');background=read(OUT/'tables/window_length_background_contrasts.tsv')
    distribution=read(OUT/'tables/paired_change_distributions.tsv');tax=read(OUT/'tables/union_taxonomic_breadth.tsv')
    pilot=read(OUT/'tables/Catovirus_11_14_15_comparison.tsv');validation=json.loads((OUT/'validation.json').read_text())
    inference11=json.loads((RUN/'inference_provenance.json').read_text());control11=json.loads((CTRL/'inference_provenance.json').read_text())
    # Additional complete distributions of neighbor-identity turnover and all-label transitions.
    identity=[];frequencies=[]
    for k in [5,10,20]:
        gained=[len(list(filter(None,r[f'top{k}_gained_saccharide_references'].split(';')))) for r in pair]
        lost=[len(list(filter(None,r[f'top{k}_lost_saccharide_references'].split(';')))) for r in pair]
        j=[float(r[f'top{k}_saccharide_Jaccard']) for r in pair if r[f'top{k}_saccharide_Jaccard']!='']
        identity.append(dict(k=k,n_pairs=len(pair),n_same_saccharide_neighbor_set=sum(a==b==0 for a,b in zip(gained,lost)),
            n_changed_nearest_saccharide_reference=sum(r['nearest_saccharide_reference15']!=r['nearest_saccharide_reference11'] for r in pair),
            saccharide_set_Jaccard_n_defined=len(j),**{f'Jaccard_q{q}':float(np.quantile(j,q)) for q in [0,.25,.5,.75,1]}))
        for n in range(k+1):frequencies.append(dict(k=k,n_references=n,n_pairs_gaining_exactly_n=gained.count(n),n_pairs_losing_exactly_n=lost.count(n)))
    write(OUT/'tables/neighbor_identity_turnover_summary.tsv',identity);write(OUT/'tables/neighbor_identity_turnover_frequencies.tsv',frequencies)
    labelcounts=collections.Counter((r['predicted_labels15'],r['predicted_labels11']) for r in pair)
    write(OUT/'tables/all_label_transition_counts.tsv',[dict(labels15=a,labels11=b,n_pairs=n) for (a,b),n in sorted(labelcounts.items())])
    # Attribute yields to paired inputs versus the three eligibility rescues, without imputing absent predictions.
    matched15={r['window15_id']:r for r in pair};oldw={r['window_id']:r for r in read(BASE/'viral_gt_windows.tsv')}
    entities={'window':{w:[w] for w in oldw},'locus':{r['locus_id']:r['window_ids'].split(';') for r in read(BASE/'candidate_loci.tsv')},'genome':{r['genome_unit_id']:r['window_ids'].split(';') for r in read(BASE/'genome_units.tsv')}}
    matched_counts=[]
    for c in ['15','11']:
        aa={w for w,r in matched15.items() if yes(r['classifier_positive'+c])};bb={w for w,r in matched15.items() if w not in aa and float(r['p10'+c])>=.5}
        for s,selected in [('A',aa),('B',bb),('C',aa|bb)]:
            for level,ee in entities.items():
                flags={w:w in selected for w in matched15};cc=collections.Counter(status(ww,flags) for ww in ee.values())
                matched_counts.append(dict(condition=c,set=s,level=level,n_scored_matched_windows=737,total_entities=len(ee),candidates=cc['candidate'],not_selected=cc['not_selected'],unresolved=cc['unresolved']))
    write(OUT/'tables/paired_only_candidate_counts.tsv',matched_counts)
    # Compact recovery summary keeps uncertain states distinct from observed losses.
    recsummary=[]
    for s in 'ABC':
        for level in ['window','locus','genome']:
            cc=collections.Counter(r['comparison'] for r in recovery if r['set']==s and r['level']==level)
            recsummary.append(dict(set=s,level=level,**dict(cc),candidate_in_either_condition=sum(n for name,n in cc.items() if name in {'both','11_only_observed_gain','11_candidate_baseline_unresolved','15_only_observed_loss','15_candidate_11_unresolved'})))
    write(OUT/'tables/candidate_recovery_summary.tsv',recsummary)
    def cn(c,s,level,key='candidates'):return next(r[key] for r in counts if r['condition']==c and r['set']==s and r['level']==level)
    counttable=[[s,cn('15',s,'window'),cn('11',s,'window'),cn('15',s,'locus'),cn('11',s,'locus'),cn('15',s,'genome'),cn('11',s,'genome')] for s in 'ABC']
    changetable=[[r['metric'],r['n_increase'],r['n_decrease'],r['n_unchanged_within_tolerance'],fmt(r['mean']),fmt(r['quantile_0.5'])] for r in distribution if r['metric'] in {'saccharide_score','p5','p10','p20'}]
    stattable=[]
    for subset in ['all_paired','complete15_to_complete11']:
        for metric in ['saccharide_score','p10']:
            r=next(r for r in stats if r['subset']==subset and r['metric']==metric)
            stattable.append([subset,metric,r['n_resampling_groups'],fmt(r['genome_locus_balanced_delta']),f"[{fmt(r['bootstrap95_low'])}, {fmt(r['bootstrap95_high'])}]"])
    bgtable=[]
    for method in ['complete_target_vs_same_genome_random','complete_target_vs_nonoverlap_random','strict_matched_retained_at11']:
        for metric in ['saccharide_score','p10']:
            r=next(r for r in background if r['method']==method and r['metric']==metric)
            bgtable.append([method,metric,r['n_genomes'],fmt(r['difference_in_differences']),f"[{fmt(r['bootstrap95_low'])}, {fmt(r['bootstrap95_high'])}]"])
    text='''# 11-gene versus frozen 15-gene viral GT windows

**Shortening changes candidate recovery substantially in both directions; it does not establish uniformly better glycan-module detection.** Catovirus shows a strong favorable example, but the overall genome-balanced change is small and uncertain, some 15-gene candidates disappear, and strict GT/property-matched controls do not demonstrate improved specificity. Retain the frozen 15-gene screen as Data Point 1 and the separate 11-gene run as Data Point 2.

## Execution and exact denominators

The same **594 original GT-A anchors** and **757 mapped genomic occurrences** were used. Every ±5-ORF extraction was checked against the original sorted all-CDS index and the archived ±7 context, including protein accessions, source sequences, coordinates, product annotations and the identical center gene. Overlapping windows remained separate for inference. Ends were neither wrapped nor compensated, and genes were not repeated or padded.

Data Point 2 contains **658 complete 11-gene windows and 99 incomplete windows**, 7,967 gene rows, and 2,721 removed outer-gene rows (occurrence rows, not distinct proteins). Frozen inference succeeded for **740 windows**, with **17 excluded**; all 737 baseline-inferred windows remain inferred, giving **737 properly matched pairs**. Three baseline exclusions became eligible and are reported separately: MK072508.1/AYV86512.1 (Polyketide; P@10=0), MG779342.1/AUV58485.1 (Saccharide; score 0.9813, P@10=1), and UPSH01000001.1/VBB17815.1 (Other; P@10=0). No baseline-inferred window lost eligibility.

Four representative original inputs reproduced their frozen scores and all 320 embedding coordinates within 10⁻⁶ before full inference. The exact checkpoint, ESM2 layer-6 residue pooling/truncation, original additive all-one mask, actual-gene contextual pooling and classifier rule max(Other,0.5) were retained. References remain the same **2,502 biological CORE_ONLY BGCs**, including 195 Saccharide-labeled BGCs. Only **33 previously uncached proteins** were encoded for newly supported inputs; all 1,552 shortened controls used existing protein vectors. No retraining, recalibration, structure prediction or Foldseek search was performed.

## Paired changes: all 737 matched occurrences

'''+table(['Metric','Increase','Decrease','Unchanged','Window mean Δ','Median Δ'],changetable)+'''

Changes are 11 minus 15. Classifier score direction uses tolerance 10⁻⁶; retrieval fractions use their exact discrete neighbor counts. **53 classifier-negative windows become positive, while 35 positives become negative**; 211 stay positive and 438 stay negative. **56 P@10-low windows reach ≥0.5, while 45 previously high windows fall below 0.5**. Among the 358 windows with changed P@10, score and retrieval move in the same direction for 297 (83.0%) and in opposite directions for 61. The median P@10 change is zero. The nearest Saccharide reference changes identity in **442/737 pairs**; the top-10 Saccharide neighbor set is unchanged in 253/737. Full score and retrieval distributions, all seven class scores, reference cosines, nearest Saccharide identities, gained/lost top-5/10/20 Saccharide neighbors and neighbor-set turnover are in the [paired table](tables/paired_window_comparison.tsv), [distribution tables](tables/paired_change_distributions.tsv) and [identity-turnover summary](tables/neighbor_identity_turnover_summary.tsv). Empty/empty neighbor sets have undefined Jaccard, not an assigned zero.

The 22 pairs with no removed genes reproduce their original scores and retrieval exactly. Removing genes changes reference geometry as well as labels; a higher score is not evidence that the removed genes were biologically unrelated.

## Genome-grouped uncertainty and background comparison

Duplicate contexts were collapsed, context changes averaged within fixed baseline loci, loci averaged within genomes, and whole genomes resampled 5,000 times. The complete 15→11 subset contains **600 raw pairs, 531 unique contexts and 356 genomes**.

'''+table(['Subset','Metric','Genomes','Balanced mean Δ','95% bootstrap interval'],stattable)+'''

The primary genome-balanced intervals include zero for both endpoints. Dropping the overrepresented generic taxid 2100421 leaves this conclusion unchanged. Grouping related genomes by the previous k21 sketch components gives a small positive score interval but an uncertain P@10 interval; these are operational sensitivity groups, not species or a phylogeny.

The **same 1,252 random 15-gene control centers** were shortened to 11 genes, along with the 239 source-annotation-matched and 61 stricter GT-proxy-matched controls. All 1,552 were inferred with the same frozen pipeline. Strict matching was checked before control inference: **26 pairs across 20 genomes** still satisfy equal gene/GT-proxy counts, AA totals within 25%, equal long-protein counts and reverse-strand counts within two after shortening. The remaining 35 strict pairs are retained as a separate baseline-matched sensitivity, not silently treated as still matched.

The contrast below is **(target11−target15)−(control11−control15)**, averaged within genome:

'''+table(['Comparison','Metric','Genomes','Difference in changes','95% interval'],bgtable)+'''

Random controls decline slightly on shortening, so the target–random contrast improves modestly. For P@10 this evidence is weak: its interval crosses zero after restricting random controls to nonoverlapping windows and after grouping genomes into sequence components (component contrast +0.023, interval −0.008 to +0.053). Strict GT/property-matched controls show **no demonstrated advantage**, with wide intervals and negative point estimates. Thus there is evidence of context-sensitive model behavior and a modest selected-cohort contrast against random neighborhoods, not established improvement in biological specificity beyond GT content.

## Candidate recovery and taxonomy

A=frozen-rule Saccharide-positive windows; B=classifier-negative windows with exploratory P@10≥0.5; C=A∪B. The threshold is not a validated biosynthetic boundary. The main locus comparison uses the **same 565 baseline regions**, preventing shortening-induced fragmentation from appearing as new recovery. Genome units retain exact nucleotide/taxid alias collapsing.

'''+table(['Set','15 windows','11 windows','15 fixed loci','11 fixed loci','15 genomes','11 genomes'],counttable)+'''

Data Point 2 has any inference for **559/565 fixed loci and 461/466 genomes**; 556 loci and 458 genomes are fully inferred. C retains **17 unresolved windows, seven unresolved fixed loci and six unresolved genomes**. These are selected-cohort frequencies only, not all-virus prevalence.

For union C, **217 fixed loci are recovered by both conditions**, **42 observed noncandidate loci become candidates at 11 genes**, one additional locus was previously unresolved, and **28 loci are candidates only at 15 genes**. At genome level, 190 are recovered by both; 28 observed noncandidate genomes and one previously unresolved genome enter at 11; 25 observed candidate genomes lose the criterion and one 15-gene candidate genome becomes unresolved. An unresolved 11-gene status is not an observed biological loss. All identities and statuses are in the [recovery registry](tables/candidate_recovery_comparison.tsv).

One of the three eligibility rescues contributes an additional classifier/union-positive locus and genome. Restricting both conditions to the same 737 scored occurrences gives 264 classifier-positive and 343 union windows at 11 genes; [paired-only entity counts](tables/paired_only_candidate_counts.tsv) separate this from full-run eligibility recovery. Independently merging the shorter intervals yields **575 native 11-gene components**, with 204 classifier-positive and 262 union-positive components; these are secondary counts because splitting changes the region denominator.

Recorded union breadth changes from **19 to 17 families and 44 to 41 genera**. **No new recorded family enters**; Kleczkowskaviridae and Orpheoviridae cease meeting the 11-gene union criterion. New genera are **Aspduovirus, Hokovirus and Megavirus**; Megavirus enters through the newly eligible MG779342.1 window rather than a scored baseline pair. Alphaorpheovirus, Betabaculovirus, Cuauhnahuacvirus, Elerivirus and Oceanusvirus become observed nonselected, while Yasminevirus has no observed 11-gene candidate but remains unresolved. With generic-source inherited ranks withheld, breadth changes from **15 to 13 families and 40 to 38 genera**. At 11 genes, only 61/219 union genomes have a non-generic resolved family and 57/219 a non-generic resolved genus; uncertainty is substantial. The [complete taxon table](tables/union_taxonomy_comparison.tsv) retains exact per-condition analyzed denominators and unresolved groups.

Combining candidate evidence from either window size would prioritize **288 fixed loci and 245 genomes**, an exploratory superset rather than a validated ensemble. The nonnested gains and losses favor retaining both data points for biological review instead of replacing Data Point 1 solely because one total count increases.

## What was removed, and the Catovirus example

For the 600 complete paired windows, 55 lose at least one GT-proxy gene, 66 lose a carbohydrate-annotated gene, and 357 lose at least three unknown/hypothetical-annotated genes. These overlapping categories use original product text and cohort-sequence proxies, not new functional predictions. Genome-block comparisons of change between annotation strata all have intervals spanning zero. They do not establish that removing unknown genes improves association or that carbohydrate-gene removal explains losses. Every removed gene is retained with its original functional annotation in the [outer-gene table](tables/removed_outer_genes.tsv); proxy definitions and stratum intervals are supplied.

'''+table(['Input','Genes','Saccharide score','Label','P@5','P@10','P@20'],[[r['input_definition'],r['n_genes'],fmt(r['saccharide_score'],6),r['predicted_labels'],r['p5'],r['p10'],r['p20']] for r in pilot])+'''

At the same ARF08732.1 occurrence, the 11-gene window **does recover stronger association**. This is a locus-level rescue in a genome already represented by other baseline candidates, not a newly discovered genus. However, its four removals include **ARF08726.1, annotated as a group-1 glycosyltransferase**, alongside ARF08725.1 (hypothetical), ARF08738.1 (DUF3626) and ARF08739.1 (nucleotidyl transferase). This is not a clean removal of four known irrelevant genes. The original 14-gene manually selected core remains separate from systematic-cohort counts. The preserved forward computes positional encoding but does not feed it into the encoder; changed gene composition affects attention and mean pooling, and this experiment does not isolate those mechanisms.

## Manuscript use and reproducibility

The supported claim is **window-boundary sensitivity and complementary candidate recovery**, with a strong illustrative Catovirus rescue. Greater glycan biosynthetic accuracy, an optimal window size, autonomous pathway completeness and actual glycan production remain unverified. Scores and P@k are correlated outputs of the same supervised representation, and all intervals are exploratory/pointwise.

- [Figure 1: paired scores and retrieval](figures/Fig1_paired_score_and_retrieval.pdf).
- [Figure 2: full paired change distributions](figures/Fig2_paired_change_distributions.pdf).
- [Figure 3: classifier/retrieval/union yields](figures/Fig3_candidate_yield_comparison.pdf).
- [Figure 4: grouped uncertainty and background contrasts](figures/Fig4_dilution_and_background_contrasts.pdf).
- [Figure 5: recorded family/genus distribution](figures/Fig5_union_taxonomy_comparison.pdf).
- [Figure S1: Catovirus boundary and removed GT](figures/FigS1_Catovirus_boundary_comparison.pdf).

PNG, vector PDF and editable SVG exports accompany the TSVs. [Methods](METHODS.md), [validation](validation.json), run provenance, source/model/reference SHA256 records and output manifests make the experiment reproducible. **The 676 frozen baseline files, 97 prior interpretation artifacts and 47 retrieval-rescue artifacts remain unchanged.** No results were pushed or published.
'''
    (OUT/'REPORT.md').write_text(text,encoding='utf-8')
    (OUT/'README.md').write_text('# Paired viral GT window-size experiment\n\nStart with [REPORT.md](REPORT.md) and [METHODS.md](METHODS.md).\n\nData Point 1: `analysis/viral_gt_mining_20261017/runs/baseline_594_v2` (unchanged).\n\nData Point 2: `analysis/viral_gt_mining_20261017/runs/window11_594_v1`.\n\nShortened existing controls: `analysis/viral_gt_mining_20261017/runs/controls11_existing_v1`.\n\nThis downstream directory contains paired comparisons, statistical contrasts, taxonomy, figures and independent validation.\n',encoding='utf-8')
    save(OUT/'analysis_manifest.json',dict(status='COMPLETE',data_point1=str(BASE),data_point2=str(RUN),control_run=str(CTRL),
        all_mapped_occurrences=757,paired_inference_occurrences=737,newly_eligible_11=3,excluded_both=17,
        independent_checks=validation['independent_checks'],new_target_protein_encodings=inference11['encoding']['newly_encoded'],
        new_control_protein_encodings=control11['encoding']['newly_encoded'],
        scripts=[dict(path=str(p.relative_to(OUT)),sha256=sha(p)) for p in sorted((OUT/'scripts').glob('*.py'))],
        hypothesis_status='context-sensitive candidate recovery; no demonstrated optimal length or biological specificity gain',published=False))
    print('Report, identity-turnover tables and paired-only candidate accounting saved',flush=True)
if __name__=='__main__':run()
