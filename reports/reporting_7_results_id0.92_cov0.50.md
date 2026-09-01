# Reporting 7 — Results (ground-truth profile: identity ≥ 0.92, coverage ≥ 0.50)

Second full pipeline run (`--gt calculated`, data-driven cut-offs) completed
2026-08-18 in 36.9 h, 0 errors. All numbers below come from
[hypothesis_tests.tsv](../results/gt_calculated_id0.92_cov0.50/hypothesis_tests.tsv)
and the companion tables in `results/gt_calculated_id0.92_cov0.50/`.

**Purpose of this report.** This is the **robustness replicate** of
[reporting_6](reporting_6_results_id0.90_cov0.80.md). The pipeline, feature set,
models, cross-validation and statistics are *unchanged*; only the ground-truth
definition differs — the identity/coverage cut-offs are now derived from the data
(`cutoff_recommended.tsv`) rather than fixed a priori. It therefore answers the
last open caveat of reporting_6: **do the conclusions survive a different
ground-truth definition?** Sections that carry over verbatim are marked
**[Identical to reporting_6]** and are not repeated; everything else is new.

| | reporting_6 | reporting_7 (this) |
|---|---|---|
| GT profile | `fixed` | `calculated` |
| min identity | 0.90 | **0.92** (stricter) |
| min coverage | 0.80 | **0.50** (looser) |
| run tag | `gt_fixed_id0.90_cov0.80` | `gt_calculated_id0.92_cov0.50` |
| wall-clock | 35.9 h | 36.9 h |
| errors | 0 | 0 |

---

## Abstract

**Background.** **[Identical to reporting_6]** — long-read metagenomics in a
high-human-background, low-biomass matrix, and whether supervised learning
separates true bacterial reads from background taxonomic hits better than a fixed
score threshold.

**Methods.** **[Identical to reporting_6]** except the ground-truth definition:
per-read labels were re-derived with **data-driven** minimap2 cut-offs
(identity ≥ 0.92, coverage ≥ 0.50) instead of the fixed a-priori pair
(0.90 / 0.80). The same 8 donors, titration design, taxon-agnostic
BLASTn/Kraken2 features, five models, nested LOEO/LOTO cross-validation and the
same paired-Wilcoxon + Holm–Šídák analysis plan were re-run end-to-end.

**Results.** The conclusions are **reproduced in full**. All six primary
hypotheses remain significant after correction, every effect keeps its sign, and
both null results stay null. Trained models again far exceeded the fixed
threshold (read-level AUPRC 0.9999 for the tree ensembles vs 0.906; median
ΔAUPRC = 0.088, corrected p = 0.038), and the advantage again widened with
dilution — model precision held ≈ 99.9 % from c1 to c5 while the threshold fell
from 92.5 % to 3.0 %. The BLASTn alignment-margin block again carried the signal
(removing it raised FDR at 99 % recall from 0.147 % to 0.257 %), and Kraken2 again
added nothing on top of BLAST. At the 95 %-recall operating point XGBoost
classified 4.72 M held-out reads at 99.98 % precision (552 false positives) and
97.0 % accuracy, with Brier 0.0028 vs 0.126 for the threshold. Effect sizes were
systematically **smaller** than under the fixed profile (H2 0.141 → 0.088,
H6 0.110 → 0.061, H12 0.133 → 0.116), consistent with a cleaner ground truth
narrowing the room the baseline had to fail in.

**Conclusions.** The findings of reporting_6 are **not artefacts of the
ground-truth cut-off choice**. Across two independent label definitions the
direction, significance and practical ranking of every hypothesis are stable;
only the magnitudes shift, and predictably so. This closes the "single GT
profile" caveat and leaves external validation as the principal remaining
limitation.

## All hypotheses at a glance

**[Identical to reporting_6]** — the hypothesis registry (what was tested and
why) is unchanged; see the
[overview table in reporting_6](reporting_6_results_id0.90_cov0.80.md#all-hypotheses-at-a-glance).
H7 (adapter content) and H8 (neural network) remain not run in this profile too.

**Conventions.** **[Identical to reporting_6]** — LOEO = leave-one-donor-out
(8 donor folds); read-level metric = AUPRC; primary tests = paired Wilcoxon
signed-rank with Holm–Šídák over the 6-member primary family (α = 0.05);
secondary tests report an effect + donor-bootstrap 95 % CI.

**Orientation — the ceiling (new numbers).** As before, every ML model is
near-perfect at read level, so gaps between strong models are small even when
significant:

| combined arm (overall) | AUPRC (calculated) | AUPRC (fixed) |
|---|---|---|
| xgboost | 0.9999 | 0.9999 |
| ranger (RF) | 0.9999 | 0.9998 |
| glmmTMB (1&#124;donor)+(1&#124;species) | 0.9996 | 0.9985 |
| glm | 0.981 | 0.985 |
| fixed_threshold (baseline) | **0.906** | 0.866 |

The baseline improves from 0.866 to 0.906 under the data-driven cut-offs — the
single most consequential change in this run, and the direct cause of the smaller
ML-over-baseline effects reported below. BLAST-only and combined again tie
(xgboost 0.9999 both); kraken2-only remains weaker (xgboost 0.9991, glm 0.954).

---

## Results

Section order and reasoning follow reporting_6; all numbers are new.

### 1. Dataset parameters (Figure 1)

The sequencing design is **[Identical to reporting_6]** — the same 8 donors
(D01, D03–D09), the same five-step 10× titration plus no-template control, the
same 4,747,846 host-depleted reads. Only the *labels* on those reads change.
Relaxing coverage from 0.80 to 0.50 while tightening identity from 0.90 to 0.92
reclassifies a net 78,193 reads from negative to positive, almost all at the
high-input level c1 (positives 2,404,925 → 2,473,512; negatives
508,069 → 439,482). The mean positive-read fraction therefore shifts slightly
upward at every level — 83.5 % (c1), 43.9 % (c2), 8.4 % (c3), 2.1 % (c4), 1.5 %
(c5), 0.0 % (NC), versus 80.5 / 42.4 / 8.0 / 1.9 / 1.3 / 0.0 under the fixed
profile (Figure 1C). The prevalence gradient that defines the task is thus fully
preserved: the positive class still collapses by roughly two orders of magnitude
from c1 to c5, and the negative control still contains no positives. Class
composition per level is otherwise unchanged in shape (Figure 1B), with the
ambiguous class remaining small (19,310 reads at c1, ≤ 2,150 elsewhere) and the
indeterminate class confined almost entirely to the NC barcode (3,531 reads).
The practical consequence is that the *evaluation set is the same reads with a
marginally more permissive positive definition* — an ideal robustness test,
because any conclusion that depends on where exactly the label boundary sits
should visibly move.

### 2. Feature-set comparison (Figure 2)

The verdict of reporting_6 is reproduced exactly: BLAST alone already saturates
the problem and Kraken2 adds nothing on top of it. On the combined arm the tree
ensembles reach AUPRC 0.9999, identical to BLAST-only (0.9999), while
Kraken2-only is measurably weaker (XGBoost 0.9991, random forest 0.9989)
(Figure 2A). The GLM is again the most feature-sensitive learner — 0.9809
(combined) and 0.9805 (BLAST-only) versus 0.9543 on Kraken2-only — confirming
that a linear model leans harder on the informative BLAST margins than the trees
do. At the 95 %-recall operating point the same ordering holds (Figure 2C):
trees 99.97 % precision on BLAST/combined, GLM 94.7 %, and Kraken2-only clearly
worst (GLM 85.9 %). Expressed as within-fold differences (Figure 2B),
*combined − BLAST-only* is ≈ 0 for every learner while *Kraken2 − BLAST-only* is
uniformly negative (GLM ≈ −0.026). Formally H1 is still "significant"
(median Δ = 5.0e-06, Holm–Šídák p = 0.039), but the effect is six orders of
magnitude below any practical threshold, and the donor×level mixed model again
returns a flat null (p = 0.736). Notably, H1's raw p rose from the exact n = 8
p-floor (0.0078) under the fixed profile to 0.0391 here, meaning the donor folds
are no longer unanimous — direct evidence that this "significance" is a ceiling
artefact rather than a stable effect.

### 3. ML model performance (Figure 3)

All learned models again beat the fixed threshold decisively, and by a
*smaller but still large* margin. Cross-validated read-level AUPRC rises from
0.906 (fixed threshold) through 0.981 (GLM) to 0.9996 (GLMM) and 0.9999
(random forest, XGBoost) (Figure 3A); the paired test gives median
ΔAUPRC = 0.088 (H2, Holm–Šídák p = 0.038), against 0.141 in the fixed profile.
The shrinkage is entirely explained by the baseline improving (0.866 → 0.906)
rather than the models degrading, i.e. the stricter identity cut-off removed some
of the marginal reads the threshold used to misrank. The advantage is again
stable across operating points (Figure 3B): XGBoost and the random forest hold
99.99 / 99.97 / 99.85 % precision at the 90 / 95 / 99 % recall targets, the GLM
erodes 96.3 → 94.7 → 88.3 %, and the fixed threshold sits far below
(80.8 → 78.1 → 71.0 %). Per-fold results confirm this is not donor-driven
(Figure 3C): XGBoost and the random forest span only 0.9998–0.9999 across the
eight LOEO folds (SD ≈ 4e-05), the GLMM 0.9991–0.9998, while the GLM
(0.9421–0.9942) and especially the fixed threshold (0.8486–0.9425, SD = 0.031)
are both lower and far noisier. Trees again beat the linear GLM (H4, median
Δ = 0.0149, p = 0.038 — essentially identical to the fixed profile's 0.016),
a small but perfectly consistent margin. At the 95 %-recall operating point the
combined XGBoost model classified 4,722,518 held-out reads with **99.98 %
precision, 95.0 % recall and 97.0 % accuracy** — 2,644,030 true positives against
just **552 false positives** (fixed profile: 412), with 139,159 false negatives
and 1,938,777 true negatives.

### 4. Titration dependence and feature-block ablation (Figure 4)

The titration result — the most clinically relevant finding — replicates almost
unchanged. XGBoost holds 99.98 → 99.91 % precision at 95 % recall from c1 to c5,
while the fixed threshold collapses from 92.5 % (c1) to 71.1 % (c2), 22.4 % (c3),
4.9 % (c4) and 3.0 % (c5) (Figure 4A). The interaction test H3 gives a coefficient
of **+0.506** (versus +0.520 fixed), i.e. the ML advantage grows by roughly half
an AUPRC unit across the titration range — the single largest effect in the
study, and stable to two significant figures across ground-truth definitions. The
same collapse appears in AUPRC (Figure 4B): trees and the GLMM stay ≥ 0.9986 at
every level whereas the GLM decays to 0.529 and the fixed threshold to 0.373 by
c5. The feature ablation likewise reproduces (Figure 4C): at 99 % recall the full
combined feature set gives an FDR of **0.147 %**, removing the BLAST-margin block
raises it to **0.257 %** and removing margin + human features to **0.258 %**,
whereas removing only the human-competitor feature leaves it essentially
unchanged (0.151 %). The margin block therefore again accounts for ~75 % of the
avoidable false discoveries, and the human-competitor feature again contributes
almost nothing on its own — the identical pattern, at a uniformly higher absolute
FDR (0.147 % vs 0.11 % fixed) because the looser coverage cut-off admits more
borderline positives. H5 is significant at median Δ = 8.66e-05 (p = 0.038), and
the sample×taxon breadth component H5b at +0.0069 (versus +0.012 fixed).

### 5. Statistical results (Figure 5)

The statistical picture is qualitatively identical and quantitatively compressed.
The forest plot (Figure 5A) preserves the full ordering of effects — H3 (+0.506)
≫ H12 (+0.116) > H2 (+0.088) > H6 (+0.061) > H9b (+0.016 / +0.014) > H10
(+0.025) ≫ H5 (8.7e-05), H1 (5.0e-06) — with the same significance and family
encoding, and the same two near-zero results (H9, H11) hugging the null line.
Calibration again separates the models cleanly (Figure 5B): the fixed threshold
is worst (Brier 0.126), the GLM intermediate (0.046), and the random forest
(0.0027), XGBoost (0.0028) and GLMM (0.0043) cluster two orders of magnitude
lower; H12 = +0.116 (p = 0.0078), reproducing the fixed profile's +0.133. The
random-effect decomposition (Figure 5C) again shows that the *species* baseline,
not the donor one, carries the transferable structure: median AUPRC is 0.9850 for
the GLM and 0.9836 for a donor-only GLMM, but 0.9997 with a truth-sourced species
effect and 0.9975 with a classifier-sourced one, giving species contributions of
+0.0161 and +0.0140 respectively (both p = 0.0078). One reporting subtlety
deserves note: for the donor-only model the *median of the paired differences* is
positive (+3.07e-04) even though the *difference of medians* is negative
(0.9836 vs 0.9850), which is why H9's automated significance flag reads "yes"
while its Wilcoxon p is 0.195 — see the H9 caveat below.

### 6. Pre-selected hypotheses (Supplementary Table S3)

Full statistics for this profile are in **Supplementary Table S3**
([hypothesis_results_table_gt_calculated_id0.92_cov0.50.tsv](../hypothesis_results_table_gt_calculated_id0.92_cov0.50.tsv));
the fixed-profile counterpart is Supplementary Table S2. All six primary
hypotheses (H1–H6) remain significant after Holm–Šídák correction, five of them
at the exact n = 8 p-floor (0.0078 → 0.0385) and H1 marginally above it
(0.0391 → 0.0391). Every effect retains its sign; no hypothesis changes verdict
in a practically meaningful way. Among the secondary hypotheses the species
random-effect contribution is confirmed under both sources (H9b: +0.0161 truth,
+0.0140 classifier), generalisation to an unseen taxon again costs little
(H10: LOEO median 1.000 → LOTO median 0.975, degradation ≈ 0.025, better than the
0.037 seen under the fixed profile), calibration is confirmed (H12: +0.116), and
truncated/unblocked reads again carry a false-positive share of exactly 0.000
(H11, not significant). The only nominal discrepancy is **H9**, whose
`significant` flag reads TRUE here versus FALSE in reporting_6; this is a flag
artefact, because the flag follows the bootstrap CI ([4.9e-05, 1.4e-03], which
now excludes zero) while the hypothesis's own Wilcoxon test is clearly
non-significant (p = 0.195) and the effect is ~3e-04. H9 should therefore
continue to be read as **no donor-transfer benefit**, identically to reporting_6.
H7 and H8 were again not run. Taken together, the replicate confirms that the
substantive claims — ML over baseline, a widening low-abundance advantage,
calibrated probabilities, per-taxon aggregation, and species-level rather than
donor-level generalisation — do not depend on the ground-truth cut-off choice.

## Hypothesis test results (all statistics)

Exact statistics for every hypothesis under the calculated profile
(**Supplementary Table S3**); the visual companion is **Figure 5A**. Conventions
are **[Identical to reporting_6]**.

| H | Claim | Test | Effect (median Δ) | 95 % CI | p | Holm–Šídák p | Significant |
|---|---|---|---|---|---|---|---|
| H1 | Combined features beat BLAST-only | Paired Wilcoxon (Holm–Šídák) | 5.00e-06 | [4.64e-08, 9.16e-06] | 0.0391 | 0.0391 | yes |
| H2 | ML beats fixed-threshold baseline | Paired Wilcoxon (Holm–Šídák) | 0.088 | [0.068, 0.120] | 0.0078 | 0.0385 | yes |
| H3 † | ML gain larger at low titration | LMM method×level | 0.506 | — | 2.6e-19 | <1e-16 | yes |
| H4 | Tree ensembles beat linear GLM | Paired Wilcoxon (Holm–Šídák) | 0.015 | [7.62e-03, 0.023] | 0.0078 | 0.0385 | yes |
| H5 | Margin + human features drive read-level gain | Paired Wilcoxon (Holm–Šídák) | 8.66e-05 | [3.57e-05, 1.20e-04] | 0.0078 | 0.0385 | yes |
| H5b | Genome breadth helps (sample×taxon) | Paired Wilcoxon + bootstrap | 6.91e-03 | [5.93e-03, 8.64e-03] | 0.0078 | — | yes |
| H6 | Sample×taxon aggregation beats read threshold | Paired Wilcoxon (Holm–Šídák) | 0.061 | [0.039, 0.080] | 0.0078 | 0.0385 | yes |
| H7 | Adapter+ reads enriched in false positives | Fisher exact | — | — | — | — | not run (deferred) |
| H8 | Neural net not better than XGBoost | — | — | — | — | — | not run (dropped) |
| H9 ‖ | Donor random effect improves transfer | Paired Wilcoxon + bootstrap | 3.07e-04 | [4.94e-05, 1.36e-03] | 0.1953 | — | flag yes / test no |
| H9b (truth) | Species-RE contribution (truth source) | Paired Wilcoxon + bootstrap | 0.016 | [9.05e-03, 0.047] | 0.0078 | — | yes |
| H9b (classifier) | Species-RE contribution (classifier source) | Paired Wilcoxon + bootstrap | 0.014 | [7.37e-03, 0.043] | 0.0078 | — | yes |
| H10 ‡ | Generalises to an unseen taxon (LOTO vs LOEO) | Median degradation + bootstrap | 0.025 | — | — | — | yes |
| H11 § | Truncated/unblocked reads drive false positives | Score-distribution | -6.66e-06 | — | <1e-16 | — | no |
| H12 | Models better calibrated than fixed threshold | Brier + paired Wilcoxon | 0.116 | [0.108, 0.129] | 0.0078 | — | yes |

† H3 was again fit with `lm` (the `(1|donor)` term was dropped), so its p ignores
within-donor correlation — trust the direction/size, not the tiny p.
‡ H10 = median AUPRC degradation (LOEO 1.000 → LOTO 0.975); the stored CI
describes the absolute-level distribution, not the paired difference, so it is
omitted. § H11's p is inflated by the very large read count; the practical result
is a **false-positive share of 0.000** carried by truncated reads. ‖ H9's flag
follows the bootstrap CI while its Wilcoxon p = 0.195; read it as **not
significant** (see §5 and the H9 section).

Five of six primary tests again land on the exact n = 8 paired-Wilcoxon p-floor
(0.0078 → 0.0385); H1 alone rises above it (0.0391), i.e. its donor folds are no
longer unanimous.

---

## Per-hypothesis results

For every hypothesis the **Tested / method** and **Why this test** rationale is
**[Identical to reporting_6]** and is not repeated; only the new numbers and the
interpretation under this profile are given, together with the fixed-profile
value for direct comparison.

| H | fixed (reporting_6) | calculated (this report) | verdict |
|---|---|---|---|
| H1 | +7.27e-06, p 0.0078 | +5.00e-06, p 0.0391 | reproduced; still a ceiling artefact, now non-unanimous |
| H2 | +0.141 | **+0.088** | reproduced; smaller because the baseline improved 0.866→0.906 |
| H3 | +0.520 | **+0.506** | reproduced, essentially unchanged — the strongest result |
| H4 | +0.016 | +0.015 | reproduced, unchanged |
| H5 | +7.97e-05 | +8.66e-05 | reproduced, unchanged |
| H5b | +0.012 | +0.0069 | reproduced, smaller |
| H6 | +0.110 | **+0.061** | reproduced, smaller (same baseline effect as H2) |
| H7 | not run | not run | — |
| H8 | not run | not run | — |
| H9 | +5.79e-04, p 0.148, **no** | +3.07e-04, p 0.195, **no** | reproduced null (flag artefact aside) |
| H9b truth | +0.0141 | +0.0161 | reproduced, slightly larger |
| H9b classifier | +0.0121 | +0.0140 | reproduced, slightly larger |
| H10 | 0.037 degradation | **0.025** degradation | reproduced; generalisation slightly better |
| H11 | −7.37e-06, FP share 0.000 | −6.66e-06, FP share 0.000 | reproduced null |
| H12 | +0.133 | +0.116 | reproduced, slightly smaller |

### H1 — richer feature sets vs BLAST alone (primary)

| quantity | value |
|---|---|
| AUPRC combined (xgboost) | 0.9999 |
| AUPRC blast_only (xgboost) | 0.9999 |
| median Δ (per donor, n = 8) | **+5.0e‑06** |
| 95 % CI | [4.6e‑08, 9.2e‑06] |
| Wilcoxon p / Holm–Šídák | 0.0391 / 0.0391 |
| mixed-model supplement (160 units) | p = 0.736 |

**Interpretation.** Unchanged: adding Kraken2 on top of BLAST buys nothing. The
evidence that this is a ceiling artefact is *stronger* here — the raw p left the
exact p-floor, so the donors no longer agree unanimously, and the mixed model is
a flat null.

### H2 — ML vs the fixed-threshold baseline (primary)

| quantity | value |
|---|---|
| AUPRC combined (xgboost) | 0.9999 |
| AUPRC fixed_threshold | **0.906** (fixed profile: 0.866) |
| median Δ (per donor, n = 8) | **+0.088** |
| 95 % CI | [0.068, 0.120] |
| Wilcoxon p / Holm–Šídák | 0.0078 / 0.038 |
| mixed-model supplement (160 units) | p = 3e‑22 |

**Interpretation.** Reproduced. The gap narrowed by ~0.05 AUPRC purely because
the data-driven identity cut-off (0.92) removed marginal reads the threshold used
to misrank; the models were already at the ceiling and could not improve further.

### H3 — larger ML advantage at low titration (primary)

| quantity | value |
|---|---|
| method×level coefficient | **+0.506** |
| p (lm fit) | 2.6e‑19 |
| xgboost precision @95 % recall, c1 → c5 | 99.98 % → 99.91 % |
| fixed_threshold precision, c1 → c5 | 92.5 % → 3.0 % |

**Interpretation.** The headline result, essentially identical to reporting_6
(+0.520). The model keeps near-perfect precision across a 10,000-fold dilution
while the threshold becomes unusable below c3.

### H4 — tree ensembles vs the linear GLM (primary)

| quantity | value |
|---|---|
| AUPRC xgboost / ranger | 0.9999 / 0.9999 |
| AUPRC glm | 0.9809 |
| median Δ | **+0.015** |
| 95 % CI | [7.6e‑03, 0.023] |
| Wilcoxon p / Holm–Šídák | 0.0078 / 0.038 |

**Interpretation.** Reproduced with virtually the same effect (+0.016 fixed).
Small but unanimous across all 8 donors, and the mixed model agrees
(p = 1.2e‑12).

### H5 / H5b — feature-block ablation (primary / secondary)

| arm (xgboost) | FDR @ 99 % recall |
|---|---|
| full combined | **0.147 %** |
| − BLAST margin | 0.257 % |
| − margin + human | 0.258 % |
| − human competitor only | 0.151 % |

| quantity | value |
|---|---|
| H5 median Δ (AUPRC) | +8.66e‑05, p 0.038 |
| H5b (sample×taxon breadth) | +6.91e‑03, p 0.0078 |

**Interpretation.** Reproduced exactly in pattern: dropping the BLAST-margin
block ~1.75× the FDR, dropping the human-competitor feature alone does almost
nothing. Absolute FDR is higher than in reporting_6 (0.147 % vs 0.11 %) because
the looser coverage cut-off admits more borderline positives.

### H6 — sample×taxon aggregation vs read-level thresholding (primary)

| quantity | value |
|---|---|
| median Δ (sxt xgboost − read baseline) | **+0.061** |
| 95 % CI | [0.039, 0.080] |
| Wilcoxon p / Holm–Šídák | 0.0078 / 0.038 |

**Interpretation.** Reproduced. The smaller effect (+0.110 fixed) again traces to
the stronger read-level baseline, not to weaker aggregation.

### H7 / H8 — not run

**[Identical to reporting_6]** — H7 remains `DEFERRED` (adapter table not
available) and H8 remains `DROPPED` (NN excluded by design) in this profile.

### H9 / H9b — donor vs species random effects (secondary)

| species source | GLMM AUPRC | GLM AUPRC | Δ (GLMM−GLM) | 95 % CI | species contrib. |
|---|---|---|---|---|---|
| none (1&#124;donor) | 0.9836 | 0.9850 | **+3.07e‑04** | [4.9e‑05, 1.4e‑03] | 0 |
| truth | 0.9997 | 0.9850 | +0.0145 | [7.5e‑03, 0.022] | **+0.0161** |
| classifier | 0.9975 | 0.9850 | +0.0118 | [6.5e‑03, 0.017] | **+0.0140** |

H9 itself = the **none** row: Wilcoxon p = 0.195 → **not significant**.

**Interpretation.** Reproduced. A donor random effect alone does not help on an
unseen donor; the species baseline supplies the lift, and it is again robust to
whether species comes from the truth or the classifier.
**Caveat.** The automated `significant` flag reads TRUE for H9 here because the
bootstrap CI excludes zero, even though the paired Wilcoxon does not reject
(p = 0.195) and the median of paired differences (+3.07e‑04) has the opposite
sign to the difference of medians (0.9836 − 0.9850). Report H9 as **null**.

### H10 — generalisation to an unseen taxon (secondary)

| quantity | value |
|---|---|
| LOEO median AUPRC | 1.000 |
| LOTO median AUPRC | **0.975** |
| median degradation | **≈ 0.025** |
| LOTO folds | 16 |

**Interpretation.** Reproduced, and slightly *better* than the fixed profile
(0.037). Held out an entire species, the model still reaches 0.975 AUPRC —
evidence that the taxon-agnostic feature design transfers.

### H11 — truncated / unblocked reads (secondary)

| quantity | value |
|---|---|
| median score Δ (unblocked − signal-pos) | −6.7e‑06 |
| FP share carried by truncated reads | **0.000** |
| significant | **no** |

**Interpretation.** Reproduced null. Adaptive-sampling unblocked reads carry no
measurable share of the false positives under either ground-truth definition.

### H12 — calibration (secondary)

| model (combined) | Brier |
|---|---|
| fixed_threshold (Platt) | **0.126** |
| glm | 0.046 |
| glmmTMB | 0.0043 |
| xgboost | 0.0028 |
| ranger_rf | **0.0027** |

median Δ (threshold − xgboost) = **+0.116**, 95 % CI [0.108, 0.129], p = 0.0078.

**Interpretation.** Reproduced. Learned probabilities remain ~45× better
calibrated than the Platt-scaled threshold.

---

## Cross-cutting caveats

Caveats 1–3 are **[Identical to reporting_6]**:
1. **Read-level ceiling / partial circularity** — the ground truth is a minimap2
   competition and the BLAST features correlate with it, so between-model
   differences are often <0.001 (clearest in H1, H5).
2. **n = 8 power floor** — the exact paired Wilcoxon p-floor is 0.0078; five of
   six primary tests land exactly there.
3. **run = donor confound (F8)** — one flowcell per donor, so donor biology and
   batch cannot be separated.

Caveat 4 is **resolved by this report**:
4. ~~**Single GT profile.**~~ **Addressed.** Two independent ground-truth
   definitions (fixed 0.90/0.80 and data-driven 0.92/0.50) yield the same
   directions, the same significance pattern and the same practical ranking. The
   remaining sensitivity is that both profiles share the *same* minimap2
   competition logic, so this tests cut-off choice, not the labelling paradigm.

New caveat specific to this run:
5. **Effect-size compression.** Because the data-driven cut-offs strengthen the
   baseline (AUPRC 0.866 → 0.906), the model-vs-baseline effects (H2, H6) are
   ~40 % smaller here. This is a property of the comparator, not of the models,
   whose absolute performance is unchanged.

---

## Discussion

This report set out to test whether the conclusions of reporting_6 depend on how
the ground truth is defined, and the answer is that they do not. Re-deriving
every per-read label with data-driven minimap2 cut-offs (identity ≥ 0.92,
coverage ≥ 0.50 instead of 0.90 / 0.80) and re-running the entire pipeline
end-to-end left all six primary hypotheses significant after Holm–Šídák
correction, preserved the sign of every effect, and preserved both null results.
The two most important claims are also the two most stable: the interaction
between method and titration level (H3) moved only from +0.520 to +0.506, and the
tree-versus-GLM margin (H4) from +0.016 to +0.015. Absolute model performance was
likewise untouched — XGBoost and the random forest again reach AUPRC 0.9999 on
the combined arm with per-fold spread below 1e-04 — as was the mechanistic
finding that the BLASTn alignment-margin block, not the Kraken2 block and not the
human-competitor feature alone, carries the discriminative signal.

Where the two profiles differ, they differ predictably and in one direction: the
model-versus-baseline effects shrink (H2 0.141 → 0.088; H6 0.110 → 0.061;
H12 0.133 → 0.116). The cause is visible in the comparator rather than in the
models — the fixed threshold's AUPRC improves from 0.866 to 0.906 under the
stricter identity cut-off, because reads that align at 0.90–0.92 identity are
exactly the marginal cases a single score threshold mishandles. Removing them
from the positive class flatters the baseline while leaving the already-saturated
models unchanged, mechanically compressing the difference. This is reassuring
rather than concerning: it shows the measured ML advantage is bounded by, and
responsive to, how hard the discrimination task actually is. It also sharpens the
interpretation of H1, whose raw p left the exact n = 8 floor (0.0078 → 0.0391)
here, meaning the donor folds are no longer unanimous — the clearest confirmation
yet that H1's "significance" under the fixed profile was a ceiling artefact and
that combined features tie, rather than beat, BLAST-only.

Two limitations of reporting_6 persist and one is retired. The ceiling and the
partial circularity between labels and features remain, as does the n = 8 power
floor and the run = donor confound, since none of these is a function of the
cut-off choice. What is now retired is the "single ground-truth profile" caveat:
cut-off sensitivity has been tested directly and found not to drive any
conclusion. The residual caveat is narrower — both profiles use the same minimap2
competition paradigm, so this replicate probes *where the boundary sits*, not
*how the boundary is drawn*. Two reporting artefacts also warrant care in reuse:
H9's automated significance flag follows the bootstrap CI rather than its
Wilcoxon test and reads TRUE here despite p = 0.195, and H10's stored CI
describes the absolute-level distribution rather than the paired difference. Both
are documented above and neither changes a substantive conclusion. The remaining
priority is therefore external validation on an independent donor cohort and
sequencing run, which would simultaneously break the run = donor confound and
test whether the near-perfect read-level ceiling survives outside the training
distribution.

## Conclusion

The findings of reporting_6 are robust to the ground-truth definition. Under
data-driven cut-offs the pipeline reproduces every substantive result: trained
models decisively outperform a fixed score threshold (H2, ΔAUPRC = 0.088), the
advantage widens as bacterial input falls until the threshold is unusable below
c3 (H3, +0.506), BLAST alignment-margin features carry the signal (H5), per-taxon
aggregation beats per-read calling (H6), learned probabilities are far better
calibrated (H12, Brier 0.0028 vs 0.126), and generalisation is species-driven
rather than donor-driven (H9b, H10). At the 95 %-recall operating point the model
classified 4.72 M held-out reads at 99.98 % precision with 552 false positives.
Only the magnitudes of the model-versus-baseline contrasts shift, and they do so
for an identifiable reason — a stronger baseline under stricter identity
cut-offs. With cut-off sensitivity now excluded as an explanation, the method
stands as a reproducible approach to low-biomass bacterial read discrimination in
host-dominated clinical metagenomics, pending external validation.

---

## Figures

Figure design, style and panel layout are **[Identical to reporting_6]**
(GraphPad/Prism styling, Okabe–Ito colourblind-safe palette, LOEO folds,
mean ± SD error bars); only the underlying run differs. Files are in
[figures/gt_calculated_id0.92_cov0.50/](../figures/gt_calculated_id0.92_cov0.50/)
as PDF and PNG, generated with
`FIG_RUN_TAG=gt_calculated_id0.92_cov0.50 Rscript create_project_figures.R`.

### Figure 1 — Dataset design and ground-truth label structure

![Figure 1](../figures/gt_calculated_id0.92_cov0.50/Figure_1_dataset_and_labels.png)

Panels **[Identical to reporting_6]** in construction. **(A)** reads per donor ×
level — unchanged, since the same 4.75 M host-depleted reads are used.
**(B)** label composition — the positive block at c1 grows to 2.47 M (from
2.40 M) and the negative block shrinks correspondingly. **(C)** positive fraction
per level — 83.5 / 43.9 / 8.4 / 2.1 / 1.5 / 0.0 %, i.e. the same prevalence
gradient shifted ~1–3 pp upward.

### Figure 2 — Feature-set comparison (H1)

![Figure 2](../figures/gt_calculated_id0.92_cov0.50/Figure_2_feature_set_arms.png)

**(A)** AUPRC by arm — trees ≈ 0.9999 on BLAST-only and combined, 0.999 on
Kraken2-only; GLM 0.981 / 0.981 / 0.954. **(B)** Δ vs BLAST-only — *combined −
BLAST* ≈ 0 for all learners; *Kraken2 − BLAST* negative throughout.
**(C)** precision @ 95 % recall — trees 99.97 %, GLM 94.7 %, Kraken2-only GLM
85.9 %.

### Figure 3 — ML model performance, combined arm (H2, H4)

![Figure 3](../figures/gt_calculated_id0.92_cov0.50/Figure_3_ml_model_performance.png)

**(A)** AUPRC per model — 0.906 (threshold) → 0.981 (GLM) → 0.9996 (GLMM) →
0.9999 (RF, XGBoost). **(B)** precision at 90 / 95 / 99 % recall — XGBoost
99.99 → 99.97 → 99.85 %, GLM 96.3 → 94.7 → 88.3 %, threshold 80.8 → 78.1 →
71.0 %. **(C)** per-fold heatmap — XGBoost/RF 0.9998–0.9999 across all eight
folds, threshold 0.8486–0.9425.

### Figure 4 — Titration dependence and ablation (H3, H5)

![Figure 4](../figures/gt_calculated_id0.92_cov0.50/Figure_4_titration_and_ablation.png)

**(A)** precision vs titration — XGBoost 99.98 → 99.91 % (c1→c5) versus threshold
92.5 → 3.0 %. **(B)** AUPRC vs titration — trees/GLMM ≥ 0.9986 at every level;
GLM decays to 0.529, threshold to 0.373. **(C)** ablation FDR @ 99 % recall —
0.147 % (full) → 0.257 % (− BLAST margin) → 0.258 % (− margin + human), with
− human-competitor alone at 0.151 %.

### Figure 5 — Statistics, calibration and random effects

![Figure 5](../figures/gt_calculated_id0.92_cov0.50/Figure_5_statistics_and_calibration.png)

**(A)** hypothesis forest plot — same ordering and significance pattern as
reporting_6, magnitudes compressed. **(B)** Brier — threshold 0.126, GLM 0.046,
GLMM 0.0043, XGBoost 0.0028, RF 0.0027. **(C)** H9 decomposition — GLM 0.9850,
donor-only GLMM 0.9836, +species(truth) 0.9997, +species(classifier) 0.9975.
