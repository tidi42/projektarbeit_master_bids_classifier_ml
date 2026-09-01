# Reporting 8 — Ground-truth comparison: `--gt fixed` vs `--gt calculated`

Head-to-head analysis of the two complete pipeline runs. Reproduced by:

```bash
Rscript compare_gt_profiles.R            # -> results/comparison_fixed_vs_calculated/*.tsv
Rscript create_gt_comparison_figures.R   # -> figures/comparison_fixed_vs_calculated/Figure_C*.{pdf,png}
```

Companion reports: [reporting_6](reporting_6_results_id0.90_cov0.80.md) (fixed
profile), [reporting_7](reporting_7_results_id0.92_cov0.50.md) (calculated
profile). This document contains only the **comparison**.

---

## 1. Introduction — what is actually different between the two runs

### The one thing that changed

Every read in this study is labelled **positive** (a genuine ZymoBIOMICS read) or
**negative** (a false-positive taxonomic hit) by aligning it with minimap2 against
the Zymo reference genomes *and* against the human genome, and asking two
questions:

1. **Does it win the competition?** The read must score better against Zymo than
   against human. *(Unchanged between the two runs.)*
2. **Is the alignment good enough?** The read must clear a minimum **identity**
   (how well the bases match) and a minimum **coverage** (how much of the read
   aligns at all). ← **This is the only thing that differs.**

| | `--gt fixed` | `--gt calculated` |
|---|---|---|
| minimum identity | **0.90** | **0.92** (stricter) |
| minimum coverage | **0.80** | **0.50** (looser) |
| where the numbers come from | chosen a priori, before seeing the data | derived from the data by `prepare_cutoff_sensitivity.R` |
| run tag | `gt_fixed_id0.90_cov0.80` | `gt_calculated_id0.92_cov0.50` |
| wall-clock | 35.9 h | 36.9 h |
| errors | 0 | 0 |

**Everything else is byte-for-byte identical** — the same 8 donors, the same
4,747,846 host-depleted reads, the same BLASTn/Kraken2 feature table, the same
five models, the same nested leave-one-donor-out cross-validation, the same
paired-Wilcoxon + Holm–Šídák statistics, the same seed (1729). Any difference in
the results is therefore attributable to the label definition and to nothing else.
That is what makes this a clean controlled experiment rather than two loosely
related analyses.

### How the data-driven cut-offs were chosen (in plain words)

The risk with a hand-picked threshold is that it is arbitrary. The sensitivity
procedure removes the arbitrariness by using a criterion that is **independent of
the cut-offs being tested**:

1. Take every read whose Zymo alignment beats its human alignment. These are
   *confident-Zymo* reads — we are sure they are microbial, without ever
   consulting identity or coverage. There are **2,562,772** of them.
2. Sweep a grid of (identity × coverage) pairs and, for each pair, record what
   share of that confident set would survive.
3. Keep the **strictest** pair that still retains **≥ 95 %** of them — strict
   enough to exclude junk, permissive enough not to throw away real signal.

The winner is **identity ≥ 0.92, coverage ≥ 0.50** (retains 97.2 %).

> **A finding worth stating up front.** The a-priori pair (0.90 / 0.80) retains
> only **94.6 %** — it *fails* the ≥ 95 % criterion. The original coverage gate of
> 0.80 was too aggressive: minimap2's `map-ont` preset under-aligns the ends of
> ONT reads, so demanding 80 % end-to-end coverage silently discarded genuinely
> microbial reads. The data-driven procedure fixes this by trading a looser
> coverage gate for a stricter identity gate. This is not a cosmetic change to a
> parameter; it is a measurable correction of a labelling bias, which is exactly
> why the comparison below matters.

### What we are testing in this report

> **Are the conclusions of this project a property of the biology and the models,
> or an artefact of where we drew the label boundary?**

Concretely: 13 testable hypotheses were evaluated independently under both label
definitions. We ask whether their **direction**, **significance** and **practical
ranking** agree, and — where numbers do move — whether the movement has an
identifiable mechanical cause.

---

## 2. Results

### 2.1 The ground truth changed by 2.8 % of reads, exactly at the boundary

Joining the two label tables read-by-read (both contain the same 4,747,846
`read_id`s) gives a complete confusion matrix of the label change:

| fixed ↓ / calculated → | positive | negative | ambiguous | indeterminate |
|---|---|---|---|---|
| **positive** | 2,677,448 | **27,548** | 0 | 0 |
| **negative** | **105,741** | 1,911,734 | 0 | 511 |
| **ambiguous** | 0 | 0 | 21,795 | 0 |
| **indeterminate** | 0 | 47 | 0 | 3,022 |

**133,847 reads (2.82 %) changed class**; the net effect is **+78,193 positives**.
The two flows have distinct, interpretable causes:

| transition | reads | identity range | coverage range | cause |
|---|---|---|---|---|
| negative → positive | 105,741 | 0.920 – 1.000 | **0.500 – 0.7999** | coverage gate relaxed 0.80 → 0.50 |
| positive → negative | 27,548 | **0.900 – 0.9200** | 0.800 – 1.000 | identity gate tightened 0.90 → 0.92 |

The identity of the reclassified reads is confined to the exact boundary bands —
no read outside `[0.90, 0.92]` identity lost its positive label, and no read
below 0.50 coverage gained one. The relabelling is therefore **surgical**: it
moves precisely the reads the two rules disagree about and leaves the other
97.2 % of the dataset untouched. The `ambiguous` class (rRNA/repeat-masked reads,
21,795) is bit-identical, as expected — that routing does not depend on the
cut-offs.

Churn is concentrated where the reads are: **88 % of it sits at c1** (93,061
gained + 24,474 lost), falling to a few hundred reads by c4–c5. Crucially, the
**prevalence gradient that defines the task is preserved** — the pooled positive
fraction moves only 0.2–2.3 percentage points per level (c1 82.0 → 84.4 %,
c2 45.5 → 46.8 %, c3 8.1 → 8.4 %, c4 1.7 → 1.9 %, c5 1.1 → 1.3 %, NC 0.0 → 0.0 %).
The negative control contains zero positives under both definitions.

### 2.2 Learned model performance is invariant; only the baseline moves

| combined arm | AUPRC fixed | AUPRC calculated | Δ |
|---|---|---|---|
| **fixed_threshold (baseline)** | 0.8520 | **0.9056** | **+0.0536** |
| glm | 0.9816 | 0.9809 | −0.0007 |
| glmmTMB | 0.9982 | 0.9996 | +0.0014 |
| ranger_rf | 0.99988 | 0.99989 | +0.00001 |
| xgboost | 0.99991 | 0.99991 | 0.00000 |

This is the central mechanical result of the comparison. **Every learned model
reproduces to 3–5 decimal places, while the fixed-threshold comparator gains
0.054 AUPRC.** The per-fold view makes the asymmetry unambiguous: all eight
threshold folds shift by +0.03 to +0.12, whereas every learned model's folds are
pinned at Δ ≈ 0.

The reason is intuitive. A single score cut-off has no way to handle reads that
align at 0.90–0.92 identity — precisely the population the stricter identity gate
removed from the positive class. Deleting its hardest cases flatters the
baseline. The trained models were already at the performance ceiling and had
nothing left to gain, so they do not move.

Calibration behaves identically: the threshold's Brier improves 0.1410 → 0.1257
while the learned models stay put (XGBoost 0.00276 → 0.00284; RF 0.00265 →
0.00274; GLMM 0.00435 → 0.00426).

At the clinical operating point (combined / XGBoost @ 95 % recall):

| | fixed | calculated |
|---|---|---|
| true positives | 2,569,746 | 2,644,030 |
| **false positives** | **412** | **552** |
| false negatives | 135,250 | 139,159 |
| true negatives | 2,017,574 | 1,938,777 |
| precision | 99.984 % | 99.979 % |
| accuracy | 97.13 % | 97.04 % |

### 2.3 The dose-response and the feature dependence replicate

The clinically decisive result — that the model's advantage grows as bacterial
input falls — is reproduced with near-identical shape:

| precision @ 95 % recall | c1 | c2 | c3 | c4 | c5 |
|---|---|---|---|---|---|
| XGBoost, fixed | 99.98 % | 99.99 % | 99.98 % | 99.93 % | 99.97 % |
| XGBoost, calculated | 99.98 % | 99.98 % | 99.97 % | 99.88 % | 99.91 % |
| threshold, fixed | 89.72 % | 68.85 % | 21.44 % | 4.39 % | 2.50 % |
| threshold, calculated | 92.46 % | 71.07 % | 22.43 % | 4.90 % | 2.96 % |
| **gap, fixed** | 10.3 pp | 31.1 pp | 78.5 pp | 95.5 pp | 97.5 pp |
| **gap, calculated** | 7.5 pp | 28.9 pp | 77.5 pp | 95.0 pp | 96.9 pp |

The gap grows ~10-fold (fixed) to ~13-fold (calculated) from c1 to c5 under both
definitions. The H3 interaction coefficient is +0.520 (fixed) vs +0.506
(calculated) — agreement to within 3 %.

Feature-block ablation is likewise reproduced in pattern, at a uniformly higher
absolute level (the looser coverage gate admits more borderline positives):

| XGBoost, FDR @ 99 % recall | fixed | calculated |
|---|---|---|
| full combined | 0.110 % | 0.147 % |
| − BLAST margin | 0.219 % (**2.00×**) | 0.257 % (**1.75×**) |
| − margin + human | 0.216 % (1.97×) | 0.258 % (1.75×) |
| − human competitor only | 0.110 % (1.00×) | 0.151 % (1.03×) |

Under both definitions, removing the BLAST-margin block roughly doubles the false
discovery rate while removing the human-competitor feature alone does nothing.
The mechanistic conclusion — the alignment margin is what carries the signal — is
cut-off independent.

### 2.4 Hypothesis-level concordance: 13/13 sign, 12/13 verdict

| H | fixed | calculated | Δ | sign | verdict |
|---|---|---|---|---|---|
| H1 | 7.27e-06 | 5.00e-06 | −2.3e-06 | ✓ | ✓ sig |
| H2 | 0.1412 | 0.0883 | **−0.0529** | ✓ | ✓ sig |
| H3 | 0.5204 | 0.5058 | −0.0147 | ✓ | ✓ sig |
| H4 | 0.0162 | 0.0149 | −0.0013 | ✓ | ✓ sig |
| H5 | 7.97e-05 | 8.66e-05 | +6.8e-06 | ✓ | ✓ sig |
| H5b | 0.0125 | 0.0069 | −0.0056 | ✓ | ✓ sig |
| H6 | 0.1096 | 0.0608 | **−0.0488** | ✓ | ✓ sig |
| H7 | not run | not run | — | — | — |
| H8 | not run | not run | — | — | — |
| H9 | 5.79e-04 | 3.07e-04 | −2.7e-04 | ✓ | ✗ flag only |
| H9b (truth) | 0.0141 | 0.0161 | +0.0020 | ✓ | ✓ sig |
| H9b (classifier) | 0.0121 | 0.0140 | +0.0019 | ✓ | ✓ sig |
| H10 | 0.0372 | 0.0248 | −0.0124 | ✓ | ✓ sig |
| H11 | −7.37e-06 | −6.66e-06 | +7.1e-07 | ✓ | ✓ n.s. |
| H12 | 0.1332 | 0.1161 | −0.0171 | ✓ | ✓ sig |

**Sign agreement 13/13. Verdict agreement 12/13.** The single disagreement is
**H9**, and it is a reporting artefact, not a scientific one: the automated
`significant` flag follows the bootstrap CI (which excludes zero in the
calculated run, [4.9e-05, 1.4e-03]) while H9's own paired Wilcoxon test is
non-significant under **both** definitions (p = 0.148 fixed, p = 0.195
calculated). The effect is ~3–6e-04 in both. H9 should be read as **null in both
runs**; the discrepancy is a threshold-crossing in an ancillary interval, not a
change of conclusion.

**The one systematic shift, and its exact cause.** Two effects move materially:
H2 (−0.0529) and H6 (−0.0488). Both are AUPRC contrasts measured *against the
fixed-threshold baseline*, and the baseline gained **+0.0536** AUPRC. In other
words, H2 and H6 lose almost precisely what the baseline gained — the shift is
fully accounted for by the comparator, not by any change in the models. Every
other effect moves by < 0.02. H12 (−0.0171) is the same phenomenon in Brier
units: the threshold's Brier improved by 0.0153.

The mixed-model supplement agrees across profiles, including on the negative
result: H1 p = 0.727 (fixed) vs 0.736 (calculated) — a flat null both times,
confirming that H1's formally "significant" paired test is a ceiling artefact.
Consistent with that, H1's raw p rose off the exact n = 8 floor (0.0078 → 0.0391)
in the calculated run, i.e. its donor folds are no longer unanimous.

The species random-effect decomposition also reproduces, and is the only place
where effects grow slightly:

| species source | fixed contribution | calculated contribution |
|---|---|---|
| none (donor only) | 0.0000 | 0.0000 |
| truth | +0.0141 | +0.0161 |
| classifier | +0.0121 | +0.0140 |

Under both definitions a donor random effect alone adds nothing while a species
baseline adds ~0.014, and the lift is insensitive to whether species comes from
the ground truth or the classifier's own guess.

### 2.5 Summary

| question | answer |
|---|---|
| Do any conclusions reverse? | **No.** 13/13 effects keep their sign. |
| Do any verdicts change? | **One (H9)**, and only via a CI flag; its hypothesis test is null in both runs. |
| Do the models perform differently? | **No.** All learned models agree to 3–5 decimals. |
| What did change? | The **baseline** improved by +0.054 AUPRC, shrinking the two baseline-referenced effects (H2, H6) by almost exactly that amount. |
| Which profile is better justified? | **Calculated.** The a-priori cut-offs retained only 94.6 % of confident-Zymo reads, below the 95 % design target; the data-driven pair retains 97.2 %. |

---

## 3. Figures

All figures are in
[figures/comparison_fixed_vs_calculated/](../figures/comparison_fixed_vs_calculated/)
as PDF (vector) and PNG (320 dpi). Style matches the main report figures
(GraphPad/Prism axes, Okabe–Ito colourblind-safe palette). Throughout,
**blue = fixed (0.90/0.80)** and **orange = calculated (0.92/0.50)**.

### Figure C1 — What the data-driven ground truth changed

![Figure C1](../figures/comparison_fixed_vs_calculated/Figure_C1_ground_truth_definition.png)

- **(A) Cut-off sensitivity landscape.** Percentage of the 2,562,772
  confident-Zymo reads retained at each (identity × coverage) pair; darker =
  more retained, values printed in each cell. The two profiles are outlined in
  black. The a-priori pair (0.90/0.80) retains **94.6 %**, *below* the 95 %
  design target, while the selected pair (0.92/0.50) retains **97.2 %**. Reading
  across a row shows that coverage is the dominant lever: at fixed identity,
  tightening coverage from 0.5 to 0.9 costs 5–6 points of retention.
- **(B) Reads relabelled, by titration level.** Diverging bars: green above the
  axis = reads that *gained* a positive label because the coverage gate was
  relaxed (105,741 total); pink below = reads that *lost* it because the identity
  gate was tightened (27,548 total). 88 % of all churn sits at c1, and by c4–c5
  the bars are visually indistinguishable from zero.
- **(C) Positive-read prevalence per level.** The two curves are nearly
  superimposed, collapsing from ~82–84 % at c1 to ~0 % at the no-template
  control. The prevalence gradient — the property that makes this a hard,
  clinically realistic problem — is unaffected by the relabelling.

### Figure C2 — Model performance is invariant to the ground-truth definition

![Figure C2](../figures/comparison_fixed_vs_calculated/Figure_C2_performance_concordance.png)

- **(A) Read-level AUPRC, combined arm.** Dumbbells connect each model's value
  under the two definitions (large blue dot = fixed, small orange dot =
  calculated). Only the fixed-threshold baseline has a visible connector
  (0.852 → 0.906); the four learned models show a single overlapping dot,
  i.e. no measurable change.
- **(B) Per-fold change in AUPRC.** Δ = calculated − fixed for each of the eight
  donor folds, with the median marked. The fixed-threshold folds scatter from
  +0.03 to +0.12; the GLM, GLMM, random-forest and XGBoost folds all sit on the
  zero line. This rules out the possibility that the aggregate agreement in (A)
  hides compensating per-donor differences.
- **(C) Calibration (Brier, lower = better).** The same asymmetry: the threshold
  improves (0.141 → 0.126) while every learned model stays within 0.0001 of its
  original value.

### Figure C3 — Dose-response and feature dependence replicate

![Figure C3](../figures/comparison_fixed_vs_calculated/Figure_C3_titration_and_ablation.png)

- **(A) Precision vs titration under both definitions.** Colour = model, line
  type = ground-truth profile. The four curves collapse into two: XGBoost holds
  ~100 % precision across the whole titration under both definitions, and the two
  threshold curves are nearly superimposed as they collapse from ~90 % to ~3 %.
  The dose-response is a property of the data, not of the labels.
- **(B) ML-over-baseline advantage (H3).** The precision gap (XGBoost − threshold)
  per level, side by side. It rises monotonically from ~8–10 pp at c1 to ~97 pp at
  c5 under both definitions — a ~10–13× growth. The calculated bars are slightly
  shorter at every level, reflecting the improved baseline rather than a weaker
  model.
- **(C) Feature ablation (XGBoost), FDR @ 99 % recall.** Mean ± SD over folds.
  Dropping the BLAST-margin block (alone, or together with the human-competitor
  features) roughly doubles the FDR under both definitions, while dropping the
  human-competitor feature alone changes nothing. The absolute level is uniformly
  higher in the calculated run because the looser coverage gate admits more
  borderline positives.

### Figure C4 — Hypothesis-level agreement

![Figure C4](../figures/comparison_fixed_vs_calculated/Figure_C4_hypothesis_agreement.png)

- **(A) Hypothesis effect sizes under both definitions.** Dumbbell plot on a
  pseudo-log axis (so sub-0.001 effects remain legible). Every hypothesis stays on
  the same side of zero; the visibly long connectors are H2, H6, H9, H5b and H10.
- **(B) Effect-size concordance.** |effect| under one definition against the
  other, log–log, with the identity line dashed. All 13 points lie on or beside
  the line across five orders of magnitude (from H1/H11 at ~1e-05 to H3 at 0.5),
  the compact visual statement of reproducibility.
- **(C) Absolute change in effect size (calculated − fixed).** The diagnostic
  panel. Bars are ordered from most negative; the orange dashed line marks
  **−0.054**, the AUPRC the fixed-threshold baseline gained. The two AUPRC
  contrasts measured against that baseline — **H2 and H6** — fall almost exactly
  on the line, i.e. they lose precisely what the comparator gained. H12 (amber)
  is the same effect in Brier units. Every remaining hypothesis (blue) moves by
  less than 0.02, and H9b_truth/H9b_classifier move slightly *upward*.

---

## 4. Conclusion

The two runs constitute a controlled experiment on the ground-truth definition:
identical reads, features, models, folds and statistics, with only the
identity/coverage rule changed. That change relabelled 2.82 % of reads, surgically
and only at the boundary, and left the prevalence gradient intact.

Under this perturbation **no conclusion of the project reverses**. All 13 testable
hypotheses keep the sign of their effect; 12 of 13 keep their verdict, and the
single exception (H9) is a bootstrap-CI flag crossing zero while its actual
hypothesis test remains non-significant in both runs. Every learned model
reproduces its performance to 3–5 decimal places, per fold as well as in
aggregate. The clinically important findings — the ML-over-threshold advantage,
its growth at low abundance, the dependence on BLAST alignment-margin features,
the calibration advantage, and species-level rather than donor-level
generalisation — are all cut-off independent.

The only systematic movement is a shrinkage of the two effects measured against
the fixed-threshold baseline (H2, H6), and it is quantitatively explained: the
baseline itself gained +0.054 AUPRC because the stricter identity gate removed
the marginal reads it handled worst. This is a property of the comparator, not of
the models, and it makes the measured ML advantage a *conservative* estimate
under the better-justified label definition.

Two practical recommendations follow. First, the **calculated profile should be
treated as the primary analysis**: the a-priori cut-offs retained only 94.6 % of
confidently microbial reads, missing their own 95 % design target, whereas the
data-driven pair retains 97.2 %. Second, the `significant` flag for
bootstrap-CI-based secondary hypotheses should be read alongside its hypothesis
test rather than alone, as H9 demonstrates.

With cut-off sensitivity now eliminated as an explanation, the principal
remaining limitations are those the comparison cannot address — the read-level
performance ceiling, the eight-donor power floor, and the run = donor confound —
all of which point to external validation on an independent cohort as the next
step.
