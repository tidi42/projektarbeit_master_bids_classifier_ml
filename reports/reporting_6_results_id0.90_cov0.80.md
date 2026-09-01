# Reporting 6 — Results (ground-truth profile: identity ≥ 0.90, coverage ≥ 0.80)

First full pipeline run (`--gt fixed`) completed 2026-08-17 in 35.9 h, 0 errors.
All numbers below come from
[hypothesis_tests.tsv](../results/gt_fixed_id0.90_cov0.80/hypothesis_tests.tsv)
and the companion tables in `results/gt_fixed_id0.90_cov0.80/`.

## Abstract

**Background.** Metagenomic long-read sequencing can detect bacterial pathogens
directly in clinical samples, but interpretation is confounded when the human
host background is high and bacterial abundance is low — the regime in which
taxonomic classifiers emit false-positive hits that are hard to separate from
genuine detections. We asked whether supervised machine learning can distinguish
true bacterial reads from background taxonomic hits better than a conventional
fixed score threshold, and where along the abundance range any advantage lies.

**Methods.** Eight donors were sequenced on Oxford Nanopore with live
adaptive-sampling host depletion after a defined community (ZymoBIOMICS Gut
Microbiome Standard) was spiked into human clinical background across a five-step
10× titration (c1 highest → c5 lowest) plus a no-template control. Per-read
ground truth was assigned by a minimap2 reference competition against the human
genome. Taxon-agnostic features derived from BLASTn alignment margins and Kraken2
assignments were used to train five models — a fixed-threshold baseline, logistic
regression (GLM), a generalised linear mixed model with donor/species random
effects (GLMM), random forest and XGBoost — under nested leave-one-donor-out
cross-validation (8 folds), with leave-one-taxon-out used to probe generalisation
to unseen organisms. Read-level performance was summarised by AUPRC, precision at
fixed recall, false-discovery rate and Brier calibration; twelve pre-selected
hypotheses were tested with paired Wilcoxon signed-rank tests and Holm–Šídák
correction across the six-member primary family.

**Results.** All trained models substantially outperformed the fixed threshold
(read-level AUPRC ≈ 0.999 for the tree ensembles vs 0.866; median ΔAUPRC ≈ 0.14,
corrected p = 0.038), and the advantage widened as bacterial input fell: model
precision remained ≈ 100 % across the whole titration while the threshold
collapsed from ≈ 90 % to a few percent. The discriminative signal resided in the
BLASTn alignment-margin features — removing them roughly doubled the
false-discovery rate — whereas adding Kraken2 on top of BLAST yielded no
measurable gain. At a 95 %-recall operating point, XGBoost classified 4.7 M
held-out reads at 99.98 % precision (412 false positives) and 97.1 % accuracy,
and the learned probabilities were far better calibrated than the threshold
(Brier ≈ 0.003 vs ≈ 0.14). Generalisation to taxa unseen in training cost only
≈ 0.037 AUPRC, and a per-**species** — not per-donor — random effect supplied the
transferable structure.

**Conclusions.** In a realistic ~99 %-human, low-biomass nanopore matrix,
supervised machine learning driven by alignment-margin features reliably
separated true bacterial reads from background taxonomic noise, with the largest
and most clinically relevant gains at low abundance and with well-calibrated
probabilities suitable for decision cut-offs. The findings are bounded by a
near-perfect read-level ceiling and an eight-donor design and await confirmation
under data-driven ground-truth cut-offs and external validation.

## All hypotheses at a glance

| H | Family | What was tested (in plain words) | Why we test it |
|---|---|---|---|
| H1 | p | Does adding Kraken2 features on top of BLAST improve the read calls? | Is the extra classifier arm worth including? |
| H2 | p | Does a trained model beat a single fixed score cut-off? | The core reason to use ML at all. |
| H3 | p | Is the model's edge over the cut-off bigger at low concentrations? | Low-abundance detection is the hardest and most clinically important case. |
| H4 | p | Do tree models (random forest / XGBoost) beat a linear model (GLM)? | Justify the extra complexity of trees. |
| H5 | p | Do the BLAST-margin + human-competitor features drive the accuracy? (drop them and re-score) | Find out which features actually carry the signal. |
| H5b | s | Do genome breadth + evenness help once reads are pooled per taxon? | Check the value of coverage-shape features at the sample level. |
| H6 | p | Does one combined call per taxon beat judging reads one at a time? | The intended clinical read-out is per taxon. |
| H7 | s | Are reads carrying sequencing adapters over-represented among false positives? | See whether adapter contamination drives errors. *(deferred — not run)* |
| H8 | s | Would a feed-forward neural network beat XGBoost? | Confirm nothing better was left out. *(dropped — not run)* |
| H9 | s | Does a per-donor effect help the model transfer to a new, unseen donor? | Test whether it generalises across people. |
| H9b | s | How much does a per-species baseline add on top of the donor effect (species from truth vs from the classifier)? | Separate donor from species structure, and check the result is robust. |
| H10 | s | How much accuracy is lost on a species the model never saw in training? | Test whether it generalises to new organisms. |
| H11 | s | Do short "unblocked" / truncated adaptive-sampling reads behave differently and cause false positives? | Check that a known ONT artefact isn't fooling the classifier. |
| H12 | s | Are the model's probabilities better calibrated than the threshold's? | Reliable probabilities matter for choosing clinical cut-offs. |

**Conventions.** LOEO = leave-one-donor-out (8 donor folds); the read-level
metric is **AUPRC**; primary tests use a **paired Wilcoxon signed-rank** across
the 8 donors with **Holm–Šídák** correction over the 6-member primary family
(α = 0.05); secondary tests report an effect + donor-bootstrap 95 % CI. Each
section adds a **Why this test** note (unit of replication, data distribution,
sample size) and a compact results table.

**Orientation — the ceiling.** Every ML model is near-perfect at read level, so
gaps between strong models are small even when significant:

| combined arm (overall) | AUPRC |
|---|---|
| xgboost | 0.9999 |
| ranger (RF) | 0.9998 |
| glmmTMB (1&#124;donor)+(1&#124;species) | 0.9985 |
| glm | 0.985 |
| fixed_threshold (baseline) | 0.866 |

BLAST-only and combined arms both reach ~0.9999 (xgboost); the kraken2-only arm
is weaker (xgboost 0.9954, glm 0.930).

---

## Results

The results below follow the five figures in order, one section each, followed by
a formal summary of the pre-selected hypotheses. All numbers are read-level
unless stated, use the **LOEO** scheme (8 held-out-donor folds), and refer to the
combined feature arm by default.

### 1. Dataset parameters (Figure 1)

The analysed dataset comprises eight donors (D01 and D03–D09), each sequenced
across a five-step ZymoBIOMICS titration (c1 = highest input mass down to
c5 = lowest) plus a no-template control (NC). After host depletion only the
non-human (“classified”) reads were retained, and their number is highly uneven
across the design (Figure 1A): sequencing depth is concentrated at the high-input
level c1 — up to 772k reads for donor D08 — and falls sharply toward c4, c5 and
NC, so the low-abundance conditions are also the read-poor, statistically hardest
cases. The ground-truth labels, assigned by a minimap2 competition between the
Zymo reference and the human genome, split each level into positive (true Zymo),
negative (false-positive), and a small ambiguous/indeterminate remainder
(Figure 1B). At c1 the positive class dominates (~2.9M reads), but the positive
block shrinks level by level until c4, c5 and NC are almost entirely negative.
Quantifying this, the mean positive-read fraction across the eight donors
collapses monotonically from ≈80 % at c1 to ≈42 % (c2), ≈8 % (c3), ≈2 % (c4),
≈1 % (c5) and ≈0 % at NC (Figure 1C). This built-in prevalence gradient is the
central challenge of the task: a naïve score threshold tuned at c1 must still
control false positives when true reads become vanishingly rare. The gradient
also motivates the LOEO evaluation and the low-abundance hypothesis (H3), because
both accuracy and its cost change dramatically along the titration. In total the
run produced 116.3M read-level predictions, from which the per-read call table
used below (4.7M combined-arm reads) was derived.

### 2. Feature-set comparison (Figure 2)

The first modelling question is whether the Kraken2 feature block adds anything
on top of the BLAST features. Three learners — logistic regression (GLM), random
forest, and XGBoost — were each trained on three feature arms: BLAST-only,
Kraken2-only, and the combined set (Figure 2). At read level the two tree
ensembles are essentially saturated: they reach an AUPRC of ≈1.000 on both the
BLAST-only and combined arms, and only slightly lower (≈0.993) on Kraken2-only
(Figure 2A). The GLM is more sensitive to the feature set, scoring ≈0.982 on
BLAST-only/combined but dropping to ≈0.922 (with a wide fold-to-fold spread) on
Kraken2-only. Expressed as a within-fold difference, adding Kraken2 to BLAST
(combined − BLAST-only) moves AUPRC by ≈0 for every learner, whereas replacing
BLAST with Kraken2 (Kraken2 − BLAST-only) is uniformly negative and worst for the
GLM (≈−0.059) (Figure 2B). The same ordering holds at the clinical operating
point of 95 % recall (Figure 2C): trees achieve ≈100 % precision on
BLAST/combined, the GLM ≈95 %, and Kraken2-only is clearly worst (GLM ≈82 %,
wide spread). Together this shows the BLAST features already carry the
discriminative signal; Kraken2 is a useful fallback when BLAST is unavailable but
contributes nothing measurable once BLAST is present. This is exactly the pattern
formalised in hypothesis H1, and it explains why the combined arm — although kept
as the default — ties rather than beats BLAST-only.

### 3. ML model performance (Figure 3)

Figure 3 reports the main performance comparison on the combined feature set,
covering all five models including the fixed-threshold baseline. Cross-validated
read-level AUPRC increases steeply from the baseline to the learned models
(Figure 3A): the fixed threshold averages ≈0.852 with a wide spread across folds,
the GLM ≈0.982, and the GLMM, random forest and XGBoost all reach ≈1.000. The
advantage is stable across operating points (Figure 3B): XGBoost (and the GLMM)
hold ≈100 % precision at the 90 %, 95 % and 99 % recall targets, whereas the GLM
erodes from 96 % to 88 % as the target tightens and the fixed threshold stays far
below (77 % to 68 %). The per-fold AUPRC heatmap confirms this is not driven by a
single easy donor (Figure 3C): XGBoost and the random forest score 1.000 in every
one of the eight LOEO folds, the GLMM 0.996–0.999, the GLM 0.948–0.994, and the
fixed threshold is both lower and noisier (0.802–0.898). Two conclusions follow.
First, any trained model decisively beats the single-threshold baseline
(hypothesis H2), and it does so consistently across held-out donors rather than
on average only. Second, the tree ensembles edge out the linear GLM
(hypothesis H4), but the margin is small (≈0.016 AUPRC) precisely because the
problem is near its performance ceiling. At the chosen 95 %-recall operating
point the combined XGBoost model classifies the 4.7M-read evaluation set with
99.98 % precision and 97.1 % accuracy (412 false positives against 2.57M true
positives), i.e. it recovers almost all true Zymo reads while emitting very few
false calls.

### 4. Titration dependence and feature-block ablation (Figure 4)

Figure 4 addresses where the machine-learning advantage comes from and which
features produce it. Splitting precision at 95 % recall by titration level shows
that XGBoost holds ≈100 % precision all the way from c1 to c5, while the fixed
threshold degrades from ≈90 % at c1 to only a few percent at c4–c5 (Figure 4A).
The gap between model and baseline therefore widens as input mass falls, meaning
the learned model helps most exactly where detection is hardest and clinically
most important — the low-abundance regime (hypothesis H3). The same picture
appears in AUPRC (Figure 4B): the tree ensembles and the GLMM stay flat near 1.0
across all five levels, whereas the GLM decays to ≈0.38 by c5 and the fixed
threshold to ≈0.3, so it is the linear and threshold-based methods that break
down under dilution. A feature-block ablation on XGBoost identifies the
responsible features (Figure 4C): using the false-discovery rate at 99 % recall
(a stricter operating point that resolves the otherwise sub-0.001 differences),
the full combined set gives an FDR of ≈0.11 %. Removing the BLAST-margin block,
or the BLAST-margin together with the human-competitor features, roughly doubles
the FDR to ≈0.21–0.22 %, whereas removing only the human-competitor feature
leaves it essentially unchanged (≈0.11 %). The BLAST-margin features — the score
gap between the best Zymo hit and the best competing hit — therefore carry most
of the read-level signal (hypothesis H5). This is mechanistically sensible: in a
~99 %-human background, the margin by which a read prefers a Zymo reference over
the human genome is the most direct evidence that it is genuinely microbial.

### 5. Statistical results (Figure 5)

Figure 5 summarises the formal statistics, probability calibration, and the
random-effect decomposition. The effect-size forest plot (Figure 5A) orders every
hypothesis by its median difference on a pseudo-log axis so that sub-0.001
effects stay legible; the largest effects are the low-abundance interaction H3
(+0.520), ML-over-baseline H2 (+0.141), calibration H12 (+0.133) and per-taxon
aggregation H6 (+0.110), while H5, H1 and the two non-significant results H9
(+5.8e-04) and H11 (−7.4e-06) sit near zero. Significance (green vs grey) and
family (primary ● vs secondary ▲) are encoded by colour and shape, giving a
one-glance map of which comparisons matter and by how much. Calibration, measured
as the Brier score where lower is better (Figure 5B), separates the models
cleanly: the fixed threshold is worst and noisiest (≈0.14), the GLM intermediate
(≈0.045), and the GLMM, random forest and XGBoost cluster near ≈0.003. The
trained models therefore return far more trustworthy probabilities than the
threshold, which matters when a downstream clinical cut-off must be chosen
(hypothesis H12). The random-effect decomposition (Figure 5C) explains where the
GLMM's generalisation comes from: median AUPRC is 0.9837 for the GLM, essentially
unchanged at 0.9843 when only a per-donor random effect is added, but rises to
0.9988 (species from the truth) or 0.9975 (species from the classifier) once a
per-species baseline is included. In other words, a donor effect alone buys
nothing for an unseen donor, whereas a species baseline provides the real lift —
and that lift is robust to whether “species” is taken from the ground truth or
inferred by the classifier itself (hypotheses H9/H9b).

### 6. Pre-selected hypotheses (Supplementary Table S2)

The pre-selected hypotheses were evaluated with a fixed analysis plan and are
reported in full in **Supplementary Table S2**
([hypothesis_results_table.tsv](hypothesis_results_table.tsv)). Each read-level
primary hypothesis (H1–H6) was tested with a paired Wilcoxon signed-rank test
across the eight LOEO donor folds and Holm–Šídák correction over the six-member
primary family (α = 0.05); secondary hypotheses report an effect with a
donor-bootstrap 95 % CI. All six primary tests are significant at the corrected
level (Holm–Šídák p = 0.0385): richer features tie BLAST (H1), ML beats the
threshold (H2), the ML advantage grows at low abundance (H3), trees beat the GLM
(H4), the BLAST-margin/human features drive the read-level signal (H5), and
per-taxon aggregation beats per-read calling (H6). Because n = 8, every primary
test lands exactly on the exact-test p-floor of 0.0078, so significance reflects
all eight donor folds agreeing in sign rather than a large per-fold effect; the
effect sizes and the mixed-model supplement are the more informative inference.
H1 is the clearest example of this caveat — it is technically significant but its
effect is ≈7e-06 and the donor×level mixed model returns p = 0.73, so the
combined arm should be read as tying, not beating, BLAST-only. Among the
secondary hypotheses, the species random-effect contribution is confirmed under
both truth and classifier sources (H9b, +0.014 and +0.012), the model generalises
to an unseen taxon with only ≈0.037 AUPRC degradation (H10), and the models are
better calibrated than the threshold (H12, +0.133). Two secondary hypotheses are
not significant: a per-donor random effect alone does not improve transfer to a
new donor (H9, p = 0.15), and truncated/unblocked adaptive-sampling reads neither
differ meaningfully in score nor carry any measurable share of the false
positives (H11, FP share = 0.000). Two further pre-listed items, H7 (adapter
enrichment) and H8 (neural-network comparison), were not run and are marked as
such in Supplementary Table S2. Overall the practically important claims — ML
over baseline, the low-abundance advantage, calibration, per-taxon aggregation
and species structure — are all supported, while the near-ceiling differences
(H1, H5) are statistically but not practically meaningful. Full per-hypothesis
statistics (effect, 95 % CI, raw and corrected p, significance, and method notes)
are listed in Supplementary Table S2 and reproduced in the table below.

## Hypothesis test results (all statistics)

Exact statistics for every hypothesis (**Supplementary Table S2**); the visual
companion is **Figure 5A** (the effect-size forest plot). Primary hypotheses
(H1–H6) are corrected for multiple testing with **Holm–Šídák**; secondary
hypotheses report an effect with a donor-bootstrap 95 % CI. Effects are
read-level AUPRC differences unless noted. This table is also written
machine-readably to
[hypothesis_results_table.tsv](hypothesis_results_table.tsv).

| H | Claim | Test | Effect (median Δ) | 95 % CI | p | Holm–Šídák p | Significant |
|---|---|---|---|---|---|---|---|
| H1 | Combined features beat BLAST-only | Paired Wilcoxon (Holm–Šídák) | 7.27e-06 | [6.81e-06, 1.41e-05] | 0.0078 | 0.0385 | yes |
| H2 | ML beats fixed-threshold baseline | Paired Wilcoxon (Holm–Šídák) | 0.141 | [0.108, 0.185] | 0.0078 | 0.0385 | yes |
| H3 † | ML gain larger at low titration | LMM method×level | 0.520 | — | 3.2e-21 | <1e-16 | yes |
| H4 | Tree ensembles beat linear GLM | Paired Wilcoxon (Holm–Šídák) | 0.016 | [7.85e-03, 0.021] | 0.0078 | 0.0385 | yes |
| H5 | Margin + human features drive read-level gain | Paired Wilcoxon (Holm–Šídák) | 7.97e-05 | [4.79e-05, 1.30e-04] | 0.0078 | 0.0385 | yes |
| H5b | Genome breadth helps (sample×taxon) | Paired Wilcoxon + bootstrap | 0.012 | [9.91e-03, 0.016] | 0.0078 | — | yes |
| H6 | Sample×taxon aggregation beats read threshold | Paired Wilcoxon (Holm–Šídák) | 0.110 | [0.074, 0.155] | 0.0078 | 0.0385 | yes |
| H7 | Adapter+ reads enriched in false positives | Fisher exact | — | — | — | — | not run (deferred) |
| H8 | Neural net not better than XGBoost | — | — | — | — | — | not run (dropped) |
| H9 | Donor random effect improves transfer | Paired Wilcoxon + bootstrap | 5.79e-04 | [-1.73e-04, 8.76e-04] | 0.1484 | — | no |
| H9b (truth) | Species-RE contribution (truth source) | Paired Wilcoxon + bootstrap | 0.014 | [6.88e-03, 0.016] | 0.0078 | — | yes |
| H9b (classifier) | Species-RE contribution (classifier source) | Paired Wilcoxon + bootstrap | 0.012 | [6.42e-03, 0.016] | 0.0078 | — | yes |
| H10 ‡ | Generalises to an unseen taxon (LOTO vs LOEO) | Median degradation + bootstrap | 0.037 | — | — | — | yes |
| H11 § | Truncated/unblocked reads drive false positives | Score-distribution | -7.37e-06 | — | <1e-16 | — | no |
| H12 | Models better calibrated than fixed threshold | Brier + paired Wilcoxon | 0.133 | [0.126, 0.149] | 0.0078 | — | yes |

† H3 was fit with `lm` (the `(1|donor)` term was dropped), so its p ignores
within-donor correlation and is over-optimistic — trust the direction/size, not
the tiny p. ‡ H10 = median AUPRC degradation (LOEO 1.000 → LOTO 0.963); the
stored CI describes the absolute-level distribution, not the paired difference,
so it is omitted here. § H11's p is inflated by the very large read count; the
practical result is a **false-positive share of 0.000** carried by truncated reads.

All six primary tests land on the exact n = 8 paired-Wilcoxon p-floor
(0.0078 → 0.0385 after correction), so significance rests on all 8 donor folds
agreeing in sign; the effect sizes and the mixed-model supplement are the more
informative inference (see the per-hypothesis sections and cross-cutting caveats).

---

## H1 — Do richer feature sets beat BLAST alone? (primary)
**Tested / method.** AUPRC of `combined(xgboost)` − `blast_only(xgboost)` per
donor, paired Wilcoxon (Holm–Šídák over the primary family).
**Why this test.** Unit of replication = **donor** (n = 8 paired folds; both arms
scored on the *same* held-out donor). Small n + bounded, ceiling-skewed AUPRC
differences → a paired **Wilcoxon signed-rank** (exact, distribution-free — no
normality assumed, uses sign + rank of the 8 differences); the 6 primary tests
share a **Holm–Šídák** family-wise correction.

| quantity | value |
|---|---|
| AUPRC combined (xgboost) | 0.9999 |
| AUPRC blast_only (xgboost) | 0.9999 |
| median Δ (per donor, n = 8) | **+7.3e‑06** |
| 95 % CI | [6.8e‑06, 1.4e‑05] |
| Wilcoxon p / Holm–Šídák | 0.0078 / 0.038 |
| significant | **yes** (see caveat) |

**Interpretation.** Adding the Kraken2 features on top of BLAST makes no real
difference — both arms already score ~0.9999, so there is nothing left to gain.
It only counts as "significant" because the tiny gap leans the same way in all
8 donors.
**Caveat.** The donor×level mixed-model supplement gives p = 0.73 (no effect) —
so H1's "significance" is a ceiling/unanimous-sign artefact, not a practical gain.

## H2 — Does ML beat the fixed-threshold baseline? (primary)
**Tested / method.** AUPRC of `combined(xgboost)` − `combined(fixed_threshold)`,
paired Wilcoxon.
**Why this test.** Same paired design as H1 — 8 donors, both scored per fold. With
n = 8 and non-normal, bounded AUPRC gaps, the exact paired Wilcoxon is the
robust choice; part of the Holm–Šídák primary family.

| quantity | value |
|---|---|
| AUPRC combined (xgboost) | 0.9999 |
| AUPRC fixed_threshold (baseline) | 0.866 |
| median Δ (per donor, n = 8) | **+0.141** |
| 95 % CI | [0.108, 0.185] |
| Wilcoxon p / Holm–Šídák | 0.0078 / 0.038 |
| mixed-model supplement (160 units) | p = 5e‑28 |
| significant | **yes** |

**Interpretation.** The main takeaway: a trained model beats a simple score
cut-off by a wide margin (0.9999 vs 0.866). A second, independent test agrees
strongly.
**Caveat.** The baseline is a single Platt-scaled score threshold; the large gap
partly reflects how weak that baseline is at low abundance (see H3).

## H3 — Is the ML advantage larger at low titration? (primary)
**Tested / method.** method×titration-level interaction on AUPRC
(`AUPRC ~ method*level + (1|donor)`).
**Why this test.** H3 is an **interaction** (does the ML−baseline gain change
across levels), not a single paired contrast — so it needs a *model* with a
method×level term fit across all donor×level cells, with `(1|donor)` for the
repeated measures within a donor. A (mixed) linear model yields that interaction
estimate and a continuous p.

| quantity | value |
|---|---|
| interaction estimate (ML gain: low − high level) | **+0.520** |
| model p | 3e‑21 |
| Holm–Šídák | ≈ 0 |
| model actually fit | `lm` (fallback) |
| significant | **yes** (see caveat) |

**Interpretation.** The model's edge over the cut-off gets bigger as there is
less material to detect. At the lowest concentrations (c4/c5) the simple
threshold falls apart while the model keeps working — so ML helps most exactly
where detection is hardest.
**Caveat.** The model fell back to plain `lm` (`note: fit=lm`) — the `(1|donor)`
term was dropped, so the p-value ignores within-donor correlation and is far too
small; trust the **direction/size**, not the tiny p.

## H4 — Do tree ensembles beat a linear GLM? (primary)
**Tested / method.** AUPRC of `xgboost` − `glm` on the combined arm (tree family
fixed by inner-CV), paired Wilcoxon.
**Why this test.** Paired over the 8 donors (both models, same folds); small-n,
non-normal differences → exact paired Wilcoxon; part of the primary family.

| quantity | value |
|---|---|
| AUPRC xgboost (combined) | 0.9999 |
| AUPRC glm (combined) | 0.985 |
| median Δ (per donor, n = 8) | **+0.016** |
| 95 % CI | [0.008, 0.021] |
| Wilcoxon p / Holm–Šídák | 0.0078 / 0.038 |
| mixed-model supplement (160 units) | p = 2e‑12 |
| significant | **yes** |

**Interpretation.** Trees do a little better than the linear model (0.9999 vs
0.985) — a small edge, but consistent across every donor and confirmed by the
second test. Being able to use feature *combinations* is what wins.
**Caveat.** Magnitude is small because the GLM is already strong; the gap widens
on the weaker kraken2-only arm (0.995 vs 0.930).

## H5 — Do margin + human-competitor features drive the read-level gain? (primary)
**Tested / method.** AUPRC of combined − `combined_minus_H5_key` (drop
BLAST-margin + human-competitor blocks) for xgboost, paired Wilcoxon. Breadth is
tested separately as **H5b**.
**Why this test.** Ablation is paired by construction (same model, same donor,
with vs without the block) → paired Wilcoxon over the 8 donors; exact/nonparametric
for small n; part of the primary family.

| quantity | value |
|---|---|
| AUPRC combined (all features) | 0.9999 |
| AUPRC combined − key blocks | 0.9998 |
| median Δ (per donor, n = 8) | **+8.0e‑05** |
| 95 % CI | [4.8e‑05, 1.3e‑04] |
| Wilcoxon p / Holm–Šídák | 0.0078 / 0.038 |
| significant | **yes** (negligible size) |

**Interpretation.** Removing the margin and human-competitor features barely
changes the score (0.9999 → 0.9998). The model is already near-perfect at read
level, so these features help but aren't essential — other features carry much
of the same signal.
**Caveat.** A near-zero ablation effect at the ceiling does **not** mean the
features are worthless; their value shows up in aggregation (H5b) and calibration.

## H5b — Does genome breadth help at the sample×taxon level? (secondary)
**Tested / method.** AUPRC of sample×taxon aggregate **with** vs **without**
`genome_breadth` + `coverage_evenness`, leave-one-donor-out.
**Why this test.** Same paired ablation logic, but at the sample×taxon level and
held out by donor → paired Wilcoxon over the 8 donors (nonparametric, small n).

| quantity | value |
|---|---|
| median Δ (full − minus-breadth) | **+0.012** |
| 95 % CI | [0.010, 0.016] |
| Wilcoxon p | 0.0078 |
| significant | **yes** |

**Interpretation.** Genome breadth and evenness only exist after pooling reads
per taxon, and they give a real, measurable boost — signal the read-level
features can't provide on their own.
**Caveat.** Sample×taxon has far fewer units than reads, so the CI rests on the
8 donors.

## H6 — Does sample×taxon aggregation beat read-level thresholding? (primary)
**Tested / method.** AUPRC of `sample_taxon(xgboost)` − read-level
`fixed_threshold`, paired Wilcoxon over donors.
**Why this test.** The two read-outs are compared on the **same 8 donors** →
paired; AUPRC differences are bounded/non-normal at small n → exact paired
Wilcoxon; part of the primary family.

| quantity | value |
|---|---|
| read-level fixed_threshold AUPRC | 0.866 |
| median Δ (sxt xgboost − read baseline) | **+0.110** |
| 95 % CI | [0.074, 0.155] |
| Wilcoxon p / Holm–Šídák | 0.0078 / 0.038 |
| significant | **yes** |

**Interpretation.** Combining all of a taxon's reads into one call works much
better than judging reads one at a time. The intended clinical read-out (a call
per taxon) is the stronger one.
**Caveat.** It compares an aggregated model score against a read-level baseline,
so part of the gain is the aggregation step itself, not only the model.

## H7 — Adapter-positive reads enriched among false positives? (secondary)
**Tested / method.** Fisher exact / odds-ratio of adapter content in FP vs TP.
**Why this test.** The planned data are a 2×2 count table (adapter± × FP/TP) with
possibly small cells → **Fisher's exact test** (valid without large-count
asymptotics), with the odds ratio as effect size.

| quantity | value |
|---|---|
| status | **not run** (`DEFERRED`, OI T1) |
| significant | — |

**Interpretation.** Not tested in this run, so there is nothing to conclude.
**Caveat.** Inactive **by design** in both GT profiles; requires the adapter
table before it can be enabled.

## H8 — Feed-forward NN vs XGBoost? (secondary)
**Tested / method.** Equivalence of a feed-forward NN to XGBoost (paired AUPRC).
**Why this test.** Would have reused the paired-over-donors Wilcoxon design (same
8 folds) had the NN been in scope.

| quantity | value |
|---|---|
| status | **not run** (`DROPPED`, note K) |
| significant | — |

**Interpretation.** Not tested — the neural network was left out on purpose.
XGBoost is the best model here.
**Caveat.** Permanent design decision, not a pending item.

## H9 — Does a donor random effect improve transfer to a new donor? (secondary)
**Tested / method.** AUPRC of `GLMM(1|donor)` − `GLM` on combined, held-out
donors; the species term is decomposed and source-swept in **H9b**.
**Why this test.** Paired over the 8 held-out donors (GLMM vs GLM on the same
fold) → paired Wilcoxon + donor-bootstrap CI; the species source is swept 3 ways
as a robustness axis (table below).

| species source | GLMM AUPRC | GLM AUPRC | Δ (GLMM−GLM) | 95 % CI | species contrib. | robust |
|---|---|---|---|---|---|---|
| none (1&#124;donor) | 0.9843 | 0.9837 | **+0.00058** | [−0.00017, 0.00088] | 0 | ✓ |
| truth | 0.9988 | 0.9837 | +0.01466 | [0.0065, 0.0170] | +0.0141 | ✓ |
| classifier | 0.9975 | 0.9837 | +0.01272 | [0.0064, 0.0164] | +0.0121 | ✓ |

H9 itself = the **none** row: Wilcoxon p = 0.15 → **not significant**.

**Interpretation.** Adding a per-donor effect alone does not help predict a new
donor. That makes sense: a brand-new donor was never seen in training, so the
model just falls back to the average — there is nothing donor-specific to apply.
**Caveat.** This is the honest donor-transfer term only; the useful structure is
the species baseline (H9b), not the donor RE.

## H9b — Species random-effect contribution (truth vs classifier source) (secondary)
**Tested / method.** AUPRC of `GLMM(1|donor)+(1|species)` − `GLMM(1|donor)`,
with species from minimap2 **truth** and from the **classifier** coarse rank.
**Why this test.** Again paired over the 8 donors (species-RE model vs donor-only
model on the same fold) → paired Wilcoxon + donor-bootstrap CI, run under two
species sources to check the lift isn't an artefact of the truth labels.

| species source | species contribution (Δ) | 95 % CI | Wilcoxon p | significant |
|---|---|---|---|---|
| truth | **+0.0141** | [0.0069, 0.0161] | 0.0078 | yes |
| classifier | **+0.0121** | [0.0064, 0.0155] | 0.0078 | yes |

**Interpretation.** What actually helps is adding a per-species effect (0.984 →
~0.998). And it barely matters whether "species" comes from the ground truth or
from the classifier's own guess — both give the same lift, so the finding is
robust.
**Caveat.** The truth source uses the same minimap2 signal that defines the
labels, so its lift is a mild upper bound; the classifier source (independent)
still gives essentially the same effect.

## H10 — Generalisation to an unseen taxon (LOTO vs LOEO) (secondary)
**Tested / method.** Median AUPRC degradation, leave-one-taxon-out vs
leave-one-donor-out (16 LOTO folds).
**Why this test.** LOEO (8 donor folds) and LOTO (16 taxon folds) have **different
fold structures**, so they are *not* one-to-one paired; the honest summary is a
descriptive **median degradation** with a resampling (bootstrap) CI over folds,
not a paired test.

| quantity | value |
|---|---|
| LOEO median AUPRC | 1.000 |
| LOTO median AUPRC | 0.963 |
| median degradation | **≈ 0.037** |
| LOTO folds | 16 |
| significant | yes |

**Interpretation.** Tested on a species it never saw in training, the model
loses only ~4 points (1.000 → 0.963). So it generalises to new species — exactly
what the taxon-agnostic feature design was meant to achieve.
**Caveat.** The `ci_lo/ci_hi` fields for this row (0.36, 0.98) report the
absolute-level distribution, **not** the paired-difference CI — read the ~0.037
degradation, not that interval.

## H11 — Do truncated/unblocked reads differ and drive false positives? (secondary)
**Tested / method.** Score-distribution comparison of unblocked vs
signal-positive negative reads + the FP share carried by truncated reads.
**Why this test.** This contrasts whole **score distributions** across two large
read subsets (many reads per group, not 8 paired summaries) → a rank-based
distribution test suits it; the FP share is a descriptive proportion.

| quantity | value |
|---|---|
| median score Δ (unblocked − signal-pos) | −7e‑06 |
| FP share carried by truncated reads | **0.000** |
| significant | **no** |

**Interpretation.** The short "unblocked" reads from adaptive sampling are not
causing false positives — they carry essentially none of them.
**Caveat.** With an FP share of exactly 0, the test is near-degenerate (little to
detect), so this is best read as "no signal" rather than a powered null.

## H12 — Are models better calibrated than the fixed threshold? (secondary)
**Tested / method.** Brier score of Platt-scaled `fixed_threshold` − Brier of
`xgboost` (positive = model better calibrated), paired over donors.
**Why this test.** Brier is a **proper scoring rule** for probability calibration;
the per-donor Brier differences are then compared with the same paired Wilcoxon
over the 8 donors (small-n, non-normal → exact/nonparametric).

| quantity | value |
|---|---|
| Brier fixed_threshold (Platt) | ≈ 0.13 |
| Brier xgboost | ≈ 0.003 |
| median Δ (per donor, n = 8) | **+0.133** |
| 95 % CI | [0.126, 0.149] |
| Wilcoxon p | 0.0078 |
| significant | **yes** |

**Interpretation.** The model's probabilities are much more trustworthy than the
threshold's (Brier ~0.003 vs ~0.13) — that matters when you need to pick a
reliable cut-off for clinical use.
**Caveat.** Both scores were Platt-scaled on train folds; calibration is reported
as Brier only (reliability curves are a further check).

---

## Cross-cutting caveats (apply to all sections)
1. **Read-level ceiling / partial circularity.** ML AUPRC is ~0.98–0.9999
   because the ground truth is a minimap2 competition and the BLAST/Kraken2
   features correlate with it; "significant" differences between strong models
   are often <0.001 and not practically meaningful (clearest in H1, H5).
2. **n = 8 power floor.** The exact paired Wilcoxon p-floor at 8 donors is
   0.0078; every primary test lands **exactly** there, i.e. significance rests on
   all 8 folds agreeing in sign. Effect sizes + the mixed-model supplement
   (~160 donor×level units) are the more informative inference.
3. **run = donor confound (F8).** One flowcell per donor, so donor biology and
   batch cannot be separated.
4. **Single GT profile.** These are the **fixed** cutoffs (id ≥ 0.90, cov ≥ 0.80);
   the `--gt calculated` robustness run (data-driven cutoffs) is not yet done, so
   cutoff-sensitivity of these conclusions is still to be confirmed.

---

## Discussion

We show that a lightweight, taxon-agnostic supervised classifier reliably
separates true ZymoBIOMICS reads from false positives in an ~99 %-human Oxford
Nanopore adaptive-sampling matrix, and that its advantage over a single fixed
score cut-off is concentrated exactly where it matters clinically — at the lowest
microbial input. Three results carry the interpretation. First, the
discriminative signal already resides in the BLAST alignment-margin features:
adding Kraken2 on top of BLAST produced no measurable read-level gain (H1), while
ablating the BLAST-margin/human-competitor block roughly doubled the FDR at 99 %
recall (H5). This is mechanistically coherent — in a host-dominated background
the decisive evidence that a read is microbial is the margin by which it prefers
a Zymo reference over the human genome, which is precisely what these features
encode. Second, every trained model reaches a near-perfect read-level ceiling
(AUPRC ≈ 0.999 for the tree ensembles), so differences among strong models
(trees vs GLM, H4: ΔAUPRC ≈ 0.016) are statistically robust but practically
negligible; the meaningful gap is between any learned model and the fixed
threshold (H2: ΔAUPRC ≈ 0.14), and it widens from ~90 % to a few percent
precision as input mass falls across the titration (H3). Third, the models
generalise: leave-one-taxon-out degraded median AUPRC by only ~0.037 (H10), and a
per-**species** random effect — not a per-donor one — supplied the transferable
structure (H9/H9b), indicating the classifier learned organism-level alignment
behaviour rather than donor idiosyncrasies. At the 95 %-recall operating point
this translated into 99.98 % precision on 4.7 M held-out reads with only 412
false positives, i.e. a usable low-biomass detector rather than a merely
rank-ordering one — a property reinforced by the models' far lower Brier scores
than the threshold (H12).

Several methodological features temper these conclusions. The read-level ceiling
reflects a **partial circularity**: the ground truth is itself a minimap2
competition against the same reference set the BLAST features align to, so the
labels and the strongest features share a common signal, inflating absolute
performance and compressing between-model differences (clearest in H1/H5).
Inference rests on **eight donors**, so every primary paired test sits exactly on
the n = 8 exact-Wilcoxon p-floor (0.0078 → 0.0385 after Holm–Šídák); significance
therefore certifies a unanimous sign across donors, and the effect sizes plus the
~160-unit donor×level mixed-model supplement are the more informative readouts.
Because one flow cell was run per donor, donor biology and batch are inseparable
(**run = donor**), and all results derive from a single ground-truth profile
(fixed, id ≥ 0.90 / cov ≥ 0.80); the data-driven ("calculated") profile is still
running, so cut-off sensitivity is not yet established. Finally, the **hac**
basecaller caps single-read identity at ≈0.95–0.98, the GLMM was fit on a
250k-read label-stratified subsample for tractability, and two pre-listed
secondary items (H7 adapter enrichment, H8 neural network) were not executed —
none of which alters the primary claims but each of which bounds their scope.

The immediate priority is the robustness replicate under data-driven cut-offs,
followed by external validation on an independent donor cohort and sequencing run
to break the run = donor confound and to test whether the near-perfect ceiling
survives outside the training distribution. Testing the adapter-content
hypothesis (H7) and reporting reliability curves alongside the Brier scores would
further harden the calibration claim (H12) — the property most relevant to any
clinical deployment where probabilities, not just rankings, drive decisions.

## Conclusion

In a realistic ~99 %-human nanopore matrix, a taxon-agnostic classifier built on
BLAST alignment-margin features recovers essentially all true ZymoBIOMICS reads
while emitting very few false positives (99.98 % precision at 95 % recall on
4.7 M held-out reads), and it retains this precision into the low-abundance
regime where a fixed threshold collapses. The gains are genuine but bounded by a
near-perfect read-level ceiling and an eight-donor design; even so, the
practically important conclusions — machine learning over baseline (H2), a
widening advantage at low input (H3), well-calibrated probabilities (H12), and
species-level rather than donor-level generalisation (H9b/H10) — are all
supported and internally consistent across figures, folds and tests. Pending the
data-driven cut-off replicate and external validation, the method is a strong,
reproducible foundation for low-biomass pathogen-read discrimination in
host-dominated clinical metagenomics.

---

## Figures

Five publication figures are in
[figures/gt_fixed_id0.90_cov0.80/](../figures/gt_fixed_id0.90_cov0.80/) as PDF
(vector) and PNG (320 dpi), generated by
[create_project_figures.R](../create_project_figures.R) directly from the
pipeline summary tables. Style is GraphPad/Prism (bold black axes, outward ticks,
no grid) with the Okabe–Ito colourblind-safe palette; model and feature-set
colours are consistent across all figures. Unless noted, read-level panels use
the **LOEO** scheme (8 held-out-donor folds, all reads, truncated reads included)
and error bars are **mean ± SD across the 8 folds**.

### Figure 1 — Dataset design and ground-truth label structure

![Figure 1](../figures/gt_fixed_id0.90_cov0.80/Figure_1_dataset_and_labels.png)

Shows what the input actually looks like before any modelling: how many reads
survive host-depletion per donor and level, and how the positive/negative ground
truth is distributed.

- **(A) Dataset design.** Heatmap of classified (non-human) reads for each of the
  8 donors (D01, D03–D09) × titration level (c1–c5 = high→low Zymo input, NC =
  no-template control). Fill = read count (log₁₀, viridis); cell labels give the
  count. Depth concentrates at c1 (max D08 = 772k) and falls sharply as input
  mass drops, so the low-abundance levels are the read-poor, hard cases.
- **(B) Ground-truth label composition.** Stacked absolute read counts per class
  (positive / negative / ambiguous / indeterminate) at each level. c1 is the
  largest bar (~2.9 M reads) and mostly positive; the positive block shrinks
  level by level until c4–c5–NC are almost entirely negative.
- **(C) Positive-read fraction vs titration.** Mean ± SD across the 8 donors of
  the per-sample positive share. It collapses monotonically — ≈ 80 % at c1,
  ≈ 42 % at c2, ≈ 8 % at c3, ≈ 2 % at c4, ≈ 1 % at c5, ≈ 0 % at NC — i.e. the
  titration is a built-in prevalence gradient (motivates the H3 low-abundance
  test).

### Figure 2 — Feature-set comparison: BLAST-only vs Kraken2-only vs combined (H1)

![Figure 2](../figures/gt_fixed_id0.90_cov0.80/Figure_2_feature_set_arms.png)

Asks whether the Kraken2 feature block adds anything beyond BLAST. Three learners
(GLM, random forest, XGBoost) are each trained on the three feature arms.

- **(A) Read-level AUPRC by feature set.** Points = mean ± SD AUPRC. The tree
  models reach ≈ 1.000 for both BLAST-only and combined; kraken2-only is slightly
  lower (≈ 0.993). The GLM sits at ≈ 0.982 for BLAST-only/combined but drops to
  ≈ 0.922 (wide SD) on kraken2-only.
- **(B) Change in AUPRC vs BLAST-only.** Per-fold ΔAUPRC relative to the same
  learner's BLAST-only score. *Combined − BLAST* (green) is ≈ 0 for every model —
  adding Kraken2 on top of BLAST changes nothing; *Kraken2 − BLAST* (orange) is
  negative, worst for the GLM (≈ −0.059). Kraken2 helps only as a fallback, never
  as an addition.
- **(C) Precision at 95 % recall.** The clinical operating point. Trees ≈ 100 %
  on BLAST/combined; the GLM ≈ 95 %; kraken2-only is worst (GLM ≈ 82 %, wide SD).
  Consistent with H1: the combined arm ties BLAST-only, it does not beat it.

### Figure 3 — Machine-learning model performance, combined feature set (H2, H4)

![Figure 3](../figures/gt_fixed_id0.90_cov0.80/Figure_3_ml_model_performance.png)

The main performance figure: all five models on the combined arm, including the
fixed-threshold baseline.

- **(A) Cross-validated read-level AUPRC.** Dot ± horizontal SD per model. The
  baseline fixed threshold is ≈ 0.852 with a wide spread; the GLM ≈ 0.982; GLMM,
  random forest and XGBoost are all ≈ 1.000 (H2: ML clearly beats the threshold).
- **(B) Precision across recall targets (90 / 95 / 99 %).** XGBoost (and GLMM)
  stay ≈ 100 % at every operating point; the GLM erodes 96 % → 88 % as the recall
  target tightens; the fixed threshold is far below (77 % → 68 %).
- **(C) Per-fold AUPRC heatmap.** One column per held-out donor fold (1–8),
  viridis fill + printed value. XGBoost and RF are 1.000 in every fold; GLMM
  0.996–0.999; GLM ranges 0.948–0.994; the fixed threshold is both lower and
  noisier (0.802–0.898), showing the trees' edge over the GLM (H4) is small but
  perfectly consistent across donors.

### Figure 4 — Titration dependence and feature-block ablation (H3, H5)

![Figure 4](../figures/gt_fixed_id0.90_cov0.80/Figure_4_titration_and_ablation.png)

Where the ML advantage comes from: performance as input mass drops, and which
feature blocks carry the signal.

- **(A) Precision vs titration, model vs baseline.** XGBoost (orange) holds
  ≈ 100 % precision @ 95 % recall across c1→c5; the fixed threshold (grey) falls
  from ≈ 90 % at c1 to a few percent at c4–c5. The gap widens as abundance drops —
  the model's advantage is largest exactly where detection is hardest (H3).
- **(B) AUPRC vs titration level, all models.** Trees and GLMM stay flat at ≈ 1.0
  across all levels; the GLM decays to ≈ 0.38 by c5 and the fixed threshold to
  ≈ 0.3, so the linear/baseline methods are the ones that break down at low input.
- **(C) Feature ablation (XGBoost), FDR @ 99 % recall (lower = better).** Full
  combined set ≈ 0.11 %; removing the BLAST-margin block, or the margin + human
  block, roughly **doubles** the FDR (≈ 0.21–0.22 %); removing only the
  human-competitor feature leaves it ≈ 0.11 %. The BLAST-margin features are the
  ones that carry the read-level signal (H5).

### Figure 5 — Statistics, calibration and random-effect decomposition

![Figure 5](../figures/gt_fixed_id0.90_cov0.80/Figure_5_statistics_and_calibration.png)

The statistical summary companion to the results table above.

- **(A) Hypothesis effect sizes (forest plot).** Median difference per hypothesis
  on a pseudo-log x-axis (so sub-0.001 effects stay legible), with 95 % CI where
  a paired difference exists. Shape marks the family (● primary, ▲ secondary) and
  colour the outcome (green = significant, grey = n.s.). Effects rank from
  H3 (+0.520) and H2 (+0.141), H12 (+0.133), H6 (+0.110) down to the near-zero
  H5, H1 and the two n.s. results H9 (+5.8e-04) and H11 (−7.4e-06). This is the
  visual form of the *all statistics* table.
- **(B) Calibration (Brier score, lower = better).** The fixed threshold is worst
  and noisiest (≈ 0.14); the GLM ≈ 0.045; GLMM, RF and XGBoost cluster near
  ≈ 0.003 — the trained models give far more trustworthy probabilities (H12).
- **(C) H9 donor vs species random effects.** Median AUPRC (combined arm) for GLM
  (0.9837), GLMM with a donor effect only (0.9843), and GLMM with donor + species
  effects using species from the truth (0.9988) or from the classifier (0.9975).
  The donor effect alone adds nothing (H9 n.s.); the lift comes from the species
  baseline, and it is robust to how species is sourced (H9b).
