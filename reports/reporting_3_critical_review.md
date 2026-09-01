# Reporting 3 — Critical review of the design and pipeline

**Scope.** An independent, adversarial review of the whole study design and the
`scripts/` pipeline (stages 00–07 + prep), with a **companion automated test
suite** ([scripts/tests/test_pipeline.R](scripts/tests/test_pipeline.R)) that
checks each section and confirms every finding. Document date: 2026-08-03.

**How this review was produced.** Every stage script was read in full; a test
harness then (a) validated the sound parts as *correctness checks* and (b)
confirmed each weakness as a *review probe*. Current result:

> **Original review: 40 / 40 correctness checks pass · 9 / 9 findings confirmed ·
> synthetic end-to-end (stages 04→07) runs clean.**
>
> **Remediation update (2026-08-05): 74 / 74 correctness checks pass. F1, F2, F4,
> F5, F6, F7, F9 and the minor items M1, M2, M3 are fixed and F3 is mitigated
> (design reasoning in
> [reporting_4](reports/reporting_4_model_selection_and_power.md)); only the
> documented design limit F8, the documented F10, and the minor item M4 remain.**

Run it with `Rscript scripts/tests/test_pipeline.R` (exit code 0 = all
correctness checks pass; findings are printed as `DETECTED`).

---

## 1. Verdict

The pipeline is **well-engineered and conceptually strong**: the ground truth is
genuinely non-circular, features are taxon-agnostic with identity columns
actively stripped, the CV is donor-grouped with **no** leakage in the primary
LOEO scheme, the honest baseline tunes its threshold on training data only, and
the metric/statistics helpers are numerically correct (all verified — §3).

The original review found **three validity-threatening issues** plus several major
and minor ones. **As of 2026-08-05 the code has been remediated: F1, F2, F4, F5,
F6, F7, F9 and the minor items M1, M2, M3 are fixed and F3 is mitigated (reasoning
in [reporting_4](reports/reporting_4_model_selection_and_power.md)); only the
design-limit F8, the documented F10, and the minor item M4 remain.** The `Status`
column marks each; the detail sections carry a status tag.

| # | Severity | Finding | Affects | Status |
|---|----------|---------|---------|--------|
| F1 | **Critical** | LOTO reuses the entire negative set in train **and** test (leakage) | H10 | ✅ fixed |
| F2 | **Critical** | Winner-take-all model/arm selection on the same folds used for testing | H1, H2, H4, H6 | ✅ fixed (inner-CV) |
| F3 | **Critical** | n = 8 gives an exact-Wilcoxon power floor; only unanimous effects are declarable | all primary tests | ⚠ mitigated (mixed model) |
| F4 | Major | H12 compares Brier of probabilities vs a hard-call baseline — does not isolate calibration | H12 | ✅ fixed (Platt) |
| F5 | Major | Sample×taxon count uses one global threshold pooled across incomparable model score scales | sample×taxon count feature | ✅ fixed (per-model) |
| F6 | Major | H5 claims "breadth" drives the gain, but breadth is never ablated | H5 | ✅ fixed (sxt ablation) |
| F9 | Major | Ground-truth identity/coverage cutoffs are fixed and never sensitivity-tested | the entire supervised target | ✅ fixed (cutoff sweep) |
| F10 | Moderate | Poisson "detectability" floor models loading, not sequenceability | sample×taxon truth, H6 | ○ open (documented) |
| F8 | Design | Run/flow-cell is perfectly confounded with donor | external validity | ○ open (design) |
| F7 | Moderate | Barcode-leakage flag is computed (randomly) but never consumed | leakage correction | ✅ resolved (diagnostic) |
| M1, M3 | Minor | Declared-but-unproduced subject props; always-NA LCA-rank columns | feature hygiene | ✅ fixed (dropped) |
| M2 | Minor | `tune_inner` passed HPs via a global `grid_row` superassignment | tuning robustness | ✅ fixed (explicit arg) |
| M4 | Minor | Holm–Šídák assumes independence under correlated tests | polish | ○ open |

The rest of this report details each with **evidence** and a **recommended
solution**. §6 is a prioritised remediation roadmap.

---

## 2. What is done well (validated, keep it)

These are not assumptions — each is asserted by the test suite (§ = test section):

- **Non-circular 3-class ground truth (§D).** Positives require minimap2
  identity ≥ 0.90 **and** coverage ≥ 0.80 **and** beating the human score;
  rRNA/plasmid/low-complexity hits are routed to `ambiguous` and excluded rather
  than forced into the binary. Ground truth is minimap2-only, i.e. independent of
  both classifiers under test.
- **Taxon-agnostic features (§A).** `model_features("combined")` provably
  contains **no** identity columns (`species/genus/top_species/top_genus/staxid/
  k2_taxid`); a `setdiff` belt-and-braces removes them even if a block leaks one.
- **Donor-grouped LOEO with no leakage (§C).** For every LOEO fold, no read is in
  both train and test, each fold tests exactly one donor, and ambiguous/
  indeterminate rows are excluded from both sides.
- **Honest baseline (§, code).** `fit_fixed_threshold` chooses its feature,
  orientation and threshold **on training reads only**, then scores test reads —
  the correct null model for H2.
- **Correct metrics/stats (§B).** AUPRC (average precision), precision/FDR at
  recall, Brier, Holm–Šídák monotonicity, exact paired Wilcoxon (matches
  `stats::wilcox.test`), the Poisson formula $P(\ge1)=1-e^{-E[N]}$, Shannon
  entropy, and PAF-interval breadth/evenness all return the analytically correct
  values.
- **Reproducibility.** Single fixed seed, pinned tool/package versions, frozen
  feature/hypothesis registries, deterministic prep.

---

## 3. Critical findings

### F1 · LOTO leaks the entire negative set into training *and* test — ✅ FIXED

**What.** In `fold_masks()` ([scripts/04_cv_splits.R](scripts/04_cv_splits.R)),
the LOTO branch builds

```r
train = keep & !(is_pos & species == test_species)
test  = keep &  (is_pos & species == test_species) | (keep & label == "negative")
```

so **every negative read is in both train and test**. The function's own comment
says negatives should be "split by donor to avoid donor leakage", but no such
split is implemented. The probe confirms it: **2350 negative rows appear in both
train and test across the 5 synthetic LOTO folds** (in the real data it is the
*whole* negative class, every fold).

**Why it matters.** H10 (generalisation to unseen taxa) is the honest transfer
test. With the negative class memorised at training time, the held-out taxon's
positives are scored against negatives the model has already seen → **LOTO
precision/AUPRC is optimistically biased**, exactly for the hypothesis meant to
be the hardest.

**Recommended solution.** Partition negatives disjointly. The cleanest design is
to **nest LOTO inside the donor structure**: hold out taxon *t*'s positives from
training **and** evaluate on negatives from **held-out donors only** (donor-split
the negatives, as the comment intended). Concretely, give each LOTO fold a
`test_donors` subset and set `test = (held-out-taxon positives ∪ negatives from
test_donors)`, `train = (other-taxon positives ∪ negatives from train_donors)`.
Add the §C disjointness assertion to LOTO so the harness fails if it ever
regresses.

---

### F2 · Post-hoc "winner" selection on the same folds used for testing — ✅ FIXED (inner-CV)

**What.** In [scripts/07_hypothesis_tests.R](scripts/07_hypothesis_tests.R),
`best_family(M,…)` picks the model family with the highest **mean AUPRC across
the eight outer test folds**, `test_H1` picks `best_single` by the same argmax,
and `test_H4` uses `pmax(RF, XGB)` per fold — all using the metric matrix `M`
that the subsequent paired Wilcoxon then tests.

**Why it matters.** Selecting the maximum over noisy estimators and then testing
that maximum is textbook **double dipping / winner's curse**: it inflates the
reported effect and biases p-values anti-conservatively. The probe quantifies the
bias under the null (all families equal): selecting the best-of-4 by fold mean
inflates the reported gap by **E[max − mean] ≈ 0.365 fold-SD** — a spurious
positive effect from noise alone. This touches **H1, H2, H4, H6**.

**Recommended solution.** Remove selection from the inference:

1. **Pre-register a single primary model family** (e.g. XGBoost) for H1/H2/H6 and
   test only it; report the others as secondary/descriptive.
2. If a data-driven choice is required, select it **inside each training fold**
   via the existing inner CV — never from outer-fold metrics — so the outer folds
   remain untouched by selection.
3. For H4, pre-register RF **or** XGB (or test both with the multiplicity already
   in the family) instead of `pmax`.
4. For "combined vs best single arm" (H1), either pre-register the comparison
   arm or use a **max-T / union-intersection permutation test** that accounts for
   the selection.

---

### F3 · The n = 8 design can only declare *unanimous* effects — ⚠ MITIGATED (mixed model)

**What.** The primary tests are exact paired Wilcoxon over 8 LOEO folds. The
smallest attainable two-sided p is $2/2^{8}=0.0078$ (all eight fold differences
share sign). The probe shows the consequence after Holm–Šídák across the 6-test
family:

- all 8 folds agree → p = 0.0078 → adjusted **0.046 < 0.05 (passes)**;
- **one** discordant fold → p = 0.0156 → adjusted **0.090 ≥ 0.05 (fails)**.

**Why it matters.** A primary hypothesis can reach significance **only if the
effect is perfectly consistent across all eight donors**. A single donor going
the other way — entirely plausible with per-donor biology and the run confound
(F8) — makes significance impossible after correction, regardless of effect size.
This is a *structural* power ceiling, not a tuning issue.

**Recommended solution.** Do not rely on the 8-fold Wilcoxon as the sole
inference.

1. **Lead with effect sizes + donor-bootstrap CIs** (already computed) and treat
   p-values as secondary — the honest framing for n = 8.
2. Add a **pre-registered one-sided** test where direction is hypothesised
   (halves the floor to $1/2^{8}=0.0039$).
3. Gain power by modelling at a **finer unit**: a mixed-effects model over
   (donor × titration-level) fold metrics, `metric ~ method + (1|donor)`, uses
   partial pooling and ~40 units instead of 8 while respecting donor grouping —
   the same structure H3 already uses. Consider it the primary inference and keep
   Wilcoxon as a robustness check.
4. If feasible later, **multiplex ≥ 2 donors per flow cell** in any new run to
   raise the effective replication and break F8.

---

## 4. Major findings

### F4 · H12 does not actually test calibration — ✅ FIXED (Platt-calibrated baseline)

`test_H12` compares a model's Brier score against a **hard-call (0/1)** baseline.
Brier mixes calibration and refinement, and against 0/1 labels it is dominated by
the operating threshold, not calibration. The probe makes this concrete: a
correctly-*ranked* but mis-calibrated probabilistic model scores **Brier = 0.202**
while a hard-call baseline scores **0.000** — i.e. the test can *favour the
baseline* on a "models are better calibrated" claim.
**Solution:** compare calibration only among probabilistic models via **reliability
curves / Expected Calibration Error**, and give the fixed-threshold baseline a
*calibrated* probability (Platt or isotonic fit on train) before any Brier
comparison.

### F5 · One global threshold across incomparable score scales — ✅ FIXED (per-arm/model)

`aggregate_sample_taxon()` computes `thr <- quantile(score, 0.75)` pooled over
**all** arms/models/folds/schemes. Model scores live on different scales
(XGBoost/GLM probabilities in [0,1] vs `fixed_threshold` raw bitscore/matches).
The probe: pooling a [0,1] model with a [0,1000] model puts the 0.75 quantile at
**≈ 512**, above which **0 %** of the first model's reads fall — so
`n_reads_above_thr` / `reads_above_thr_per_million` are not comparable across
models. (The headline sample×taxon metric uses `max_score`, so this corrupts the
*count feature*, not the primary ranking — hence Major, not Critical.)
**Solution:** compute the threshold **per (arm, model)** — ideally from the
train-fold score distribution — or drop the raw count in favour of the
scale-free `max_score`/rank aggregate.

### F6 · H5's "breadth" claim is never tested — ✅ FIXED (sample×taxon ablation)

H5 is stated as "margin + human-competitor **+ breadth** drive the gain", but
`cfg$ablation_sets` only ablated `blast_margin` and `human_competitor`;
`genome_breadth`/`coverage_evenness` exist only at the sample×taxon level, no
ablation set removed them, and stage 05 **skips** any ablation that "drops no
read-level features". So the breadth component was **asserted, not measured**
(probe confirmed).
**Fix applied.** `06_evaluate.R` now carries a dedicated sample×taxon ablation
(`sample_taxon_ablation()`): it aggregates the out-of-fold read scores to
(library × taxon) cells, joins genome-breadth/coverage-evenness, then fits the
sample×taxon aggregate score **with vs without** the breadth block under
leave-one-donor-out CV and reports per-fold AUPRC for both (`sxt_full` vs
`sxt_minus_breadth`, registered in `cfg$sxt_ablation_sets`). Stage 07 pairs them as
**H5b** (median ΔAUPRC + donor-bootstrap CI, reported in the secondary family so the
six-test primary cap is preserved). The breadth claim is now **measured**; the
harness asserts both variants are produced and that withholding breadth moves the
sample×taxon AUPRC.

### F9 · Ground-truth cutoffs are fixed and never stress-tested — ✅ FIXED (cutoff sweep)

The whole supervised target hinges on identity ≥ 0.90 / coverage ≥ 0.80, which
are unvalidated for a **hac r10.4.1** basecaller (true-Zymo identity mode
≈ 0.95–0.98): the label set could shift materially with the cutoff.
**Fix applied.** New `scripts/prepare_cutoff_sensitivity.R` (run after stage 01)
collects **confident** Zymo reads at the high-titre level(s) — those whose Zymo
minimap2 hit beats human, a criterion **independent of** the identity/coverage
cutoffs under test — reports their empirical identity/coverage distribution
(`work/cutoff_zymo_distribution.tsv`), and sweeps a grid of (identity × coverage)
cutoffs recording the fraction of that confident set retained
(`work/cutoff_sensitivity.tsv`). It recommends the strictest cutoffs that still
retain ≥ `cfg$params$gt_retain_frac` (0.95). This quantifies how sensitive the
positive label set — the input to *every* endpoint — is to the cutoff, so key
results can be shown stable across a defensible range (C4 elevated from optional).
Method note added to [reporting_1](reports/reporting_1_methodology.md) §7.

### F10 · The Poisson floor models loading, not sequenceability

`expected_present` = (input cells × rel. abundance × recovery ≥ floor). It
captures "≥ 1 genome went into the tube", **not** whether adaptive sampling +
yield actually produced reads. R6 correctly counts undetected-but-expected taxa
as FN, so the *recall* accounting is honest — but `p_detect` should not be read as
a sequencing detection probability.
**Solution:** document it explicitly as a **loading-based lower bound**;
optionally add an observed-depth term (reads/​genome) for a sequencing-aware
detectability estimate, and keep `poisson_p_min` as the sensitivity knob.

---

## 5. Design & minor findings

### F8 · Run confounded with donor (design)
All 8 flow cells carry exactly one donor (probe confirmed), so any per-run
technical effect (chemistry, load, date) is inseparable from donor. LOEO removes
it but generalisation is to "**new donor + new run jointly**". Acknowledged in
`check_run_confounder()`; state the generalisation target explicitly, and add a
shared control library per flow cell (or multiplex donors) in future runs to make
the run effect estimable.

### F7 · Leakage flag computed but never used — ✅ RESOLVED (principled diagnostic)
`correct_barcode_leakage()` sets `suspected_leakage` by **randomly** flagging a
share of c4/c5 positives equal to the negative-barcode positive rate, but 05/06/07
never reference it (probe confirmed). It is both **inert** and, as written, **not
principled** (random flagging).
**Solution:** either remove it, or (a) flag by a real criterion — e.g. a taxon
present in a co-loaded high-mass barcode but absent from the negative — and (b)
consume it in a sensitivity analysis that re-runs the endpoints with flagged
reads down-weighted or excluded.

### Minor
- **M1 — declared-but-unproduced features — ✅ FIXED.** `subject_assembly_level`
  and `subject_is_wgs_draft` were in `cfg$subject_props` but `parse_blast` emits
  only `subject_genome_len` (BLAST `outfmt 6` carries no assembly-level / WGS-draft
  metadata, and no NCBI `assembly_summary` lookup is wired), so they were silently
  dropped by the `intersect` in stage 05. **Removed** from `cfg$subject_props` → the
  declared subject props now match what is produced (`subject_genome_len`).
- **M2 — hidden coupling — ✅ FIXED.** `tune_inner()` used to pass the grid row
  through a global `grid_row <<-` read from the model closure's calling frame.
  It now passes the row **explicitly** as the closure's third argument
  (`predict_fun(a, b, row)`); the `fit_ranger` / `fit_xgboost` closures read
  `row$mtry` / `row$nrounds` etc., and the global superassignment is gone.
- **M3 — inert LCA-rank features — ✅ FIXED.** `lca_rank` (BLAST) / `k2_lca_rank`
  (Kraken2) were always `NA` (no taxonomy tree wired). **Dropped** from
  `cfg$feature_blocks$blast_margin` / `$kraken2` and from `parse_blast` /
  `parse_kraken2`, so the feature set reflects reality; the species/genus grouping
  used by the count/entropy/margin features is unaffected.
- **M4 — multiplicity under dependence.** The 6 primary tests reuse the same
  folds and often the same base model, so they are correlated; **Šídák** assumes
  independence and is mildly anti-conservative here. Plain **Holm** is valid under
  arbitrary dependence — or use a permutation **max-T** for exact FWER control.
- **Negative-class heterogeneity.** Negatives pool residual human + donor bacteria
  + unaligned survivors — fine for the Zymo-vs-not task, but it limits the FP
  interpretation (already documented; keep it explicit).

---

## 6. Recommended remediation roadmap

Priority order (impact × how load-bearing), with the fix and where it lives:

| Order | Finding | Fix | Where | Status |
|-------|---------|-----|-------|--------|
| 1 | **F1** LOTO leakage | Donor-split negatives (nest LOTO in donor grouping); add LOTO disjointness assertion | `04_cv_splits.R` `fold_masks`, test §C | ✅ done |
| 2 | **F2** selection bias | Pre-register one primary model family; select inside inner CV only; drop `pmax` | `07_hypothesis_tests.R` | ✅ done (inner-CV) |
| 3 | **F3** power | Effect-size-first + one-sided pre-registration + mixed-model primary inference | `07_hypothesis_tests.R`, report framing | ⚠ mixed model added; one-sided/framing pending |
| 4 | **F9** cutoff sensitivity | Empirical c1 identity/coverage sweep | `prepare_cutoff_sensitivity.R` + `reporting_1` §7 | ✅ done |
| 5 | **F4** calibration | Reliability/ECE among probabilistic models; calibrate the baseline | `06_evaluate.R`, `07` `test_H12` | ✅ done (Platt) |
| 6 | **F6** breadth ablation | Sample×taxon ablation (`sxt_full` vs `sxt_minus_breadth`) → H5b | `06_evaluate.R` `sample_taxon_ablation`, `07` `test_H5b` | ✅ done |
| 7 | **F5** cross-model threshold | Per-(arm,model) threshold from train scores | `06_evaluate.R` `aggregate_sample_taxon` | ✅ done |
| 8 | **F7** leakage flag | Consume it (sensitivity) or remove it | `02`/`06` | ✅ done (diagnostic) |
| 9 | **M1, M3** feature hygiene | Drop unproduced subject props + always-NA LCA ranks | `00_config.R`, `03_build_features.R` | ✅ done |
| 10 | **M2** tuning robustness | Pass the inner-CV grid row explicitly (drop global `grid_row`) | `05_train_models.R` `tune_inner` | ✅ done |
| 11 | **F10, M4, F8** | Document / tidy as above | various | ○ open |

**Effort note.** F1, F2, F5, F6, F7, M1–M3 are small, localized code changes.
F3, F4, F9 are analysis-design changes (and pre-registration) rather than large
rewrites. None require re-collecting data; F8 is the only one that would benefit
from a future wet-lab change.

---

## 7. Bottom line

The scaffold is sound and, uniquely, **self-checking** — the accompanying suite
locks in the properties that are correct and will fail loudly if a fix regresses
them. **As of 2026-08-05 the inference layer and the ground-truth layer have both
been remediated:** F1 (LOTO leakage), F2 (selection bias), F4 (calibration), F5
(cross-model threshold), F6 (breadth ablation), F7 (leakage flag) and F9
(ground-truth cutoff sensitivity) are **fixed**, and F3 (power) is **mitigated**
with a donor-clustered mixed model — reasoning in
[reporting_4](reports/reporting_4_model_selection_and_power.md), all validated by
the suite (**74/74**). The headline read-level, the sample×taxon / transfer, and
the label-definition axes are now defensible; the feature set reflects reality
(M1, M3 dropped the declared-but-inert columns) and the inner-CV tuner carries no
hidden global state (M2). **Remaining:** the documented design limit F8 (run⟂donor,
needs a wet-lab change), the documented F10 (Poisson floor semantics), and one
minor polish item M4 (multiplicity under dependence).

*Companion artefacts:* [scripts/tests/test_pipeline.R](scripts/tests/test_pipeline.R)
(this review's evidence), [reports/reporting_4_model_selection_and_power.md](reports/reporting_4_model_selection_and_power.md)
(remediation reasoning), [open_items_v3.md](open_items_v3.md) (live tracker),
[reports/reporting_1_methodology.md](reports/reporting_1_methodology.md) (method),
[reports/reporting_2_masking_and_progress.md](reports/reporting_2_masking_and_progress.md)
(masking + scope).
