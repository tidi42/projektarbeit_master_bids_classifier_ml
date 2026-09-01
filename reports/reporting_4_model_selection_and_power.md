# Reporting 4 — Model selection, the n = 8 problem, and the negative control

This report gives the **reasoning** behind three design decisions taken while
remediating the critical review (reporting_3): the human-depleted negative
control (item 2), the data-driven model choice (item 3), and the additional test
for the n = 8 power problem (item 4). It is the design-rationale companion to
reporting_3 (which tracks *what* was fixed) and to
[scripts/tests/test_pipeline.R](scripts/tests/test_pipeline.R) (which validates
the changes: **56/56 correctness checks pass**). Document date: 2026-08-03.

---

## 1. The negative control: is it human-depleted, and does that make sense? (item 2)

**Is it depleted?** **Yes.** Stage 01 processes *every* library identically —
`groundtruth_align` → `deplete_human` → `nonhuman.fastq` → QC/BLAST/Kraken2
([scripts/01_external_tools.R](scripts/01_external_tools.R)) — and the loop runs
over the whole sample sheet, which includes the per-donor **negative (c0)**
barcodes. So the negative goes through the same GRCh38 depletion as the spiked
libraries.

**Is the premise "the negative should only be human" correct?** **Only partly —
it is a useful thing to check, but the negative is not pure human.** The negative
barcode is the donor's *no-spike* clinical sample: the same matrix as the spiked
barcodes minus the ZymoBIOMICS community. A clinical sample is not sterile human
DNA — it carries the donor's **endogenous microbiome + reagent kitome + residual
human**. After depletion the negative therefore *retains* non-human reads, and
those surviving reads are precisely the **realistic false-positive background**.

**Why depleting the negative is correct (reasoning).**

1. **Shared read universe.** H1 (combined vs single arm) and every paired
   comparison must be evaluated on the *same* read universe across libraries. If
   the negative were *not* depleted, it would be scored on a ~99 %-human universe
   while spiked libraries were scored on the non-human universe — confounding the
   comparison. Identical processing is a design requirement.
2. **The task is defined on non-human reads.** The classifier separates true Zymo
   from false positives *among non-human survivors*. The negative's surviving
   non-human reads **are** the hard negatives the model must reject. Leaving human
   in would bury that signal under trivially-removable human.
3. **If the negative really were pure human, depletion is harmless.** It would
   simply leave ~0 reads and contribute nothing. The fact that negatives retain
   reads is *informative*: those reads are endogenous microbiome / kitome, and any
   that map to Zymo are contamination / cross-barcode leakage.

**What we added (also the item-7 decision).** A principled per-run diagnostic,
`estimate_leakage()` → `work/leakage_estimates.tsv`, reporting for each run: the
negative's surviving-read count, its non-Zymo **background** count, its
Zymo-**positive** (contamination) count, and the **leakage upper bound** (the
negative's positive rate). This (a) *verifies* the negative is endogenous /
contamination background rather than a depletion failure, and (b) supplies the
principled cross-barcode-leakage upper bound. Inspect it: a negative with
near-zero survivors means that donor's matrix was effectively pure human
(contributes little, no harm); a negative with many Zymo-mapping reads is genuine
contamination to account for (feeds the deferred H7 and a leakage sensitivity
analysis).

**Bottom line:** depleting the negative is the correct, necessary choice; the
negative is the donor's endogenous + contamination background, not pure human, and
that is exactly what makes it a useful hard-negative / false-positive control.

---

## 2. A data-driven model choice without double-dipping (item 3)

**The problem (F2).** Choosing the "best" model by its performance on the *same*
outer folds that the paired tests then use is a winner's curse: the test suite
shows that selecting the best-of-4 by fold mean inflates the reported gap by
**E[max − mean] ≈ 0.365 fold-SD under the null** — a spurious effect from noise.

**The requirement.** Choose the model data-drivenly, and compare models on
*overall* fit, *taxon-specific* fit, and *concentration-specific* fit.

**Solution — two deliberately separated jobs.**

**(a) Selection for inference — leakage-free (nested CV done correctly).**
Stage 05 now records, for the combined arm × each model × each outer fold, the
mean **inner-CV AUPRC** (`inner_cv_score`, default hyper-parameters) →
`results/inner_cv_scores.tsv`. Stage 07's `select_primary_model()` ranks families
by mean inner-CV AUPRC and picks the primary family for H1/H2/H4/H6; H4 no longer
takes `pmax(RF, XGB)` but selects the tree family the same way.
*Reasoning:* the inner CV lives entirely **inside each outer training set**, so the
outer test folds — the ones the paired tests consume — are never used for
selection. Model selection is thus part of the training procedure and is evaluated
out-of-(outer)-fold, which is the textbook nested-CV cure for F2, while remaining
fully data-driven.

**(b) Comparison for interpretation — a multi-faceted leaderboard.**
Stage 06 `model_comparison()` → `results/model_comparison.tsv` reports, per
(arm, model), AUPRC / precision@recall / FDR at three **facets**:

| facet | grouping | question it answers |
|-------|----------|---------------------|
| `overall` | all reads | which model fits best on average |
| `taxon` | per species | is a model better/worse for a *specific pathogen* |
| `concentration` | per titration level c1–c5 | is a model better/worse at the *clinically critical low abundances* |

pooled across LOEO folds for stable estimates.
*Reasoning:* a single scalar hides that a model can win overall yet lose on a
specific taxon or at low titration (where detection matters most). The three
facets are exactly the "overall / specific taxon / specific concentration" views
requested, and let you see whether the inner-CV-selected family is also *robust*
across taxa and concentrations before committing to it.

**Why the two are kept separate (the key point).** Using the multi-faceted
*outer* leaderboard to pick-then-test on the same folds would re-introduce F2. So
the leaderboard is **descriptive** — it informs interpretation and the
pre-registration of the primary model for a confirmatory run — while the
**inner-CV** signal drives the automatic, leakage-free selection used by the
tests. Concretely: read `model_comparison.tsv` to understand behaviour; the
inference already uses the inner-CV-selected family; if the leaderboard reveals a
taxon- or concentration-specific weakness, report it as a caveat or pre-register a
per-stratum model for the next run.

---

## 3. The n = 8 problem: the most scientific additional test (item 4)

**The problem (F3).** Eight LOEO folds give an exact paired-Wilcoxon two-sided
p-floor of $2/2^{8}=0.0078$; after Holm–Šídák a primary hypothesis can reach
significance only if **all eight donors agree in sign**, and a single discordant
donor ($p=0.0156$) already fails (adj $0.090$). Crucially this is a property of
**any** exact test on 8 exchangeable units — switching to a donor-level
permutation / sign-flip test keeps the same $2^{8}$ granularity and the same
floor. The floor cannot be removed at the donor level; it can only be escaped by
**adding (properly modelled) exchangeable units**.

**The solution — a donor-clustered linear mixed model.** Each donor contributes
**five titration levels**, i.e. genuine within-donor replication. Fitting

$$\text{AUPRC} \sim \text{method} + (1\mid\text{donor})$$

over the per-(donor × titration-level) metric rows uses ≈ 40 units for the
method fixed-effect test, while the $(1\mid\text{donor})$ random intercept absorbs
donor clustering (and, since run ≡ donor, the run confound). This is the standard,
defensible way to gain power with a small number of clusters without
pseudo-replication.

**Why this and not the alternatives.**

- *Donor-level permutation / sign-flip:* still $2^{8}$ granularity → same floor. **Rejected.**
- *Pool all reads, test at read level:* reads within a donor are not independent →
  massive pseudo-replication, invalid p-values. **Rejected.**
- *Treat donor × level as independent rows in a plain test:* ignores clustering →
  anti-conservative. The random intercept is exactly what fixes this. **This is why the mixed model, not a plain lm, is used.**
- *Consistency:* the pre-registered H3 already uses
  `metric ~ method × level + (1|donor)`, so extending the same model to H1/H2/H4
  is coherent with the plan.

**Implementation.** `mixed_effect_test()`
([scripts/utils.R](scripts/utils.R)) + `run_mixed_supplement()`
([scripts/07_hypothesis_tests.R](scripts/07_hypothesis_tests.R)) fit the LMM for
the H1/H2/H4 comparisons on per-(donor × level) AUPRC and write
`results/mixed_model_supplement.tsv` (method estimate, p, `n_units`, engine)
beside the Wilcoxon results. `lmerTest` supplies the Satterthwaite p-value; if it
is unavailable the code falls back to a donor-fixed-effect `lm`.

**The recommended inference layer (framing).** For 8 donors, report all three,
in this order of emphasis:

1. **Effect sizes + donor bootstrap CIs** (already computed) — the most honest
   primary summary at small n.
2. **Donor-clustered mixed model** — the *powered* confirmatory test, using
   within-donor replication.
3. **Exact paired Wilcoxon** — retained as a *conservative* robustness check
   (the unanimity test).

Plus a **pre-registered one-sided** direction wherever a direction is
hypothesised (halves the floor to $1/2^{8}=0.0039$).

**Caveat.** The mixed model treats titration level as replication for H1/H2/H4
(level is itself the *effect* only in H3). Report the H1/H2/H4 models both with and
without a `level` fixed effect as a robustness check, since a strong method × level
interaction (H3) would otherwise leak into the method estimate.

---

## 4. What changed (pointer to reporting_3 §"Fixed")

| Item | Change | File(s) |
|------|--------|---------|
| 1 (F1) | LOTO donors partitioned so no read is shared train/test | `04_cv_splits.R` |
| 2 | Confirmed negative is depleted + not pure human; added leakage diagnostic | `01` (verified), `02_ground_truth_labels.R` |
| 3 (F2) | Inner-CV `select_primary_model` + `model_comparison` (overall/taxon/concentration) | `05`, `06`, `07` |
| 4 (F3) | Donor-clustered `mixed_effect_test` + `run_mixed_supplement` | `utils.R`, `07` |
| 5 (F4) | Platt-calibrated `fixed_threshold`; H12 fair Brier; calibration includes it | `05`, `06`, `07` |
| 6 (F5) | Per-(arm, model) sample×taxon threshold | `06_evaluate.R` |
| 7 (F7) | Random flag removed; principled `estimate_leakage` diagnostic | `02_ground_truth_labels.R` |

Validated by the test suite: **70/70 correctness checks pass**; F1–F7 and F9 now
pass as correctness/mitigation checks (F6 breadth ablation and F9 cutoff sensitivity
were added 2026-08-04 — see [reporting_3](reporting_3_critical_review.md)). Still
open (documented, not code-fixable here): **F8** (run⟂donor, design), **F10**
(Poisson floor semantics), **M1–M4** (minor polish).
