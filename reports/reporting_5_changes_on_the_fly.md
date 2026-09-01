# Reporting 5 — On-the-fly fixes: stage-02/03 crashes + stage-05 speedup

Short log of the problems that stopped or stalled the run and how they were
fixed. The three bugs were **latent** — unreachable until stage 01 first
completed end-to-end (48/48, 24.4 h) — and the fourth was a **performance**
problem that only showed once training (stage 05) ran on real data.
Document date: 2026-08-11 – 08-15.

---

## Bug 1 — `flag_ambiguous_regions()` wrong column names (fatal)

**Problem.** Stage 02 aborted at labelling start with
`Error in eval(...) : object 'tname' not found`. The function referenced bare
`tname`/`tstart`/`tend`, but the per-read label table carries the Zymo best-hit
coords under the `zymo_` prefix (`zymo_tname`/`zymo_tstart`/`zymo_tend`).

**Fix.** Use the `zymo_*` columns; drop reads with no Zymo hit (NA coords, which
cannot fall in a Zymo ambiguous region); map the overlap results back to the
original rows via a `.row` index.
See [scripts/02_ground_truth_labels.R](scripts/02_ground_truth_labels.R#L73).

## Bug 2 — `read_qc.tsv` corrupted by ONT BAM tags

**Problem.** These ONT FASTQs carry tab-separated BAM tags in the header
(`BC:Z:`, `qs:f:`, `du:f:`, `ns:i:`, …). `seqkit fx2tab --name` emitted the whole
header, so `length`/`GC`/`avg.qual` were shifted into trailing columns
(header = 4 cols vs data = 24–26 cols). This would feed corrupt QC features into
stage 03.

**Fix.** Added `--only-id` to the `seqkit fx2tab` call so only the read ID is
kept. The 48 existing `read_qc.tsv` were regenerated from `nonhuman.fastq`
out-of-band (**0/48 header/data mismatches** afterwards).
See [scripts/01_external_tools.R](scripts/01_external_tools.R).

## Bug 3 — stage 03 lost `read_id` (knock-on of the Bug 2 fix)

**Problem.** After the stage-02 fixes the run cleared stage 02 but crashed
entering stage 03 with
`The following columns listed in by are missing from x: [read_id]`. The Bug 2
fix (`--only-id`) changed the `read_qc.tsv` header column from `#name` to `#id`,
but `parse_qc`'s rename map only knew `#name`/`name`, so the QC table kept `#id`
and the merge on `read_id` failed.

**Fix.** Added `#id`/`id` → `read_id` to the rename map (kept `#name`/`name` for
back-compat). See
[scripts/03_build_features.R](scripts/03_build_features.R#L124).

## Perf 4 — stage 05 training ran ~2 cores wide (would take weeks)

**Problem.** Stages 02–04 finished, but after ~25 h stage 05 was still on arm 1
of 3, fold 1 of 8 (only ~16 of ~500 random forests done) — a projected **5+
weeks**. Measurement showed the 64‑core/502 GB box was **97 % idle**, no I/O‑wait,
no swapping: the process held many threads but only **2 were doing work**. Cause:
`fit_ranger`/`fit_xgboost` never passed a thread count, and this ranger/xgboost
build **defaults to 2 threads**, so every fit crawled on 2 cores.

**Actions taken.**
1. **ranger → all cores:** pass `num.threads` to both `ranger()` calls.
2. **xgboost → all cores + `tree_method="hist"`** on both fits (the GPU isn't
   usable from this R build — it silently falls back to CPU — so CPU `hist` gives
   the equivalent speedup). Was pinned to `nthread = 2`.
3. New config param [`train_threads`](scripts/00_config.R) (resolves to 60 here);
   thread count changes **speed only, not the models**.
4. **Trees, middle way** (good science + efficiency): keep the **final** models at
   `num.trees = 500`; cut only the throwaway inner‑CV **tuning** forests
   300 → 100 (they just rank hyperparameters).

**Result.** A forest dropped from ~40 min → ~1 min; ranger now spreads to ~15
cores (~7×). Stage 05 estimate: **weeks → ~15–20 h**. Required killing the run
(no checkpointing) and relaunching `--gt fixed --from 05` (stages 02–04 outputs
reused).

## Bugs 5 & 6 — xgboost & glmmTMB failed on the real run (stage 05)

**Problem.** Once training actually ran (fast ranger got us there), every
**xgboost** fold failed with `argument "y" is missing` and every **glmmTMB**
fold with `re.form must equal NULL, NA, or ~0`. Both are API mismatches with the
installed package versions that the earlier crawl never reached: xgboost 3.x
replaced the old `xgboost(data, …)` with a new sklearn-style `xgboost(x, y, …)`,
and glmmTMB's `predict` does not accept an lme4-style `re.form = ~(1|species)`.
The driver caught each per fold, so the run kept going but produced **no
xgboost/glmmTMB predictions**.

**Fix.** xgboost → the stable low-level `xgb.train(params, data, nrounds)`.
glmmTMB → `re.form = NULL`; with `allow.new.levels = TRUE` the held-out LOEO
donor still collapses to population level while trained species REs are applied,
so the intended "integrate out donor, apply species" semantics is preserved.
See [scripts/05_train_models.R](scripts/05_train_models.R). Both validated on a
real leave-one-donor-out split (80k predictions each, all finite in [0,1]);
relaunched `--from 05`.

---

## Test used

[scripts/tests/test_pipeline.R](scripts/tests/test_pipeline.R) — full suite
re-run after each fix: **95 passed, 0 failed**. Extra checks: the out-of-band
regeneration of the 48 QC files (0/48 mismatches, Bug 2), and an out-of-band
feature build for the first and last libraries (28 cols, `read_id` present,
Bug 3). For Perf 4, an out-of-band timing on the real 4.7M-row feature table
confirmed a 500-tree forest at 60 threads takes ~5.7 min (vs ~40 min on 2
cores), and the live relaunch was verified to spread ranger across ~15 cores.
The run resumed each time (`--gt fixed`, now `--from 05`) and progressed past
the crash/stall points with 0 new errors.

---

## Perf 7 — GLMM training-row cap (stage-05 speed lever, 40 h deadline)

**Problem.** With the trees now fast, the remaining bottleneck was the GLMM. On
the `combined` arm (24 fixed effects + `(1|donor)+(1|species)`) a single
`glmmTMB` fit on the full ~4.16 M training rows takes ~6–7 h, and the GLMM is fit
~48× in stage 05 (8 LOEO folds + the 24-fit H9 species-RE decomposition). Hourly
snapshots showed ~25 of the last ~30 h were spent inside GLMM fits — projecting
stage 05 to **~1–2 weeks**, far past the required ~40 h.

**Change.** Cap the GLMM's *training* rows via a **label-stratified subsample at
the natural class ratio**:
- New config knob
  [`glmm_max_train_rows = 250000`](scripts/00_config.R) (`NA` = no cap).
- New helper
  [`stratified_cap_rows()`](scripts/05_train_models.R) returns row indices to keep
  so `nrow ≤ cap` while keeping the *same sampling fraction per class* (uses
  `sample.int` on each class's index vector — avoids `sample()`'s length-1 trap).
- [`glmm_fit_predict()`](scripts/05_train_models.R) applies it to the training
  frame only when it exceeds the cap.

**Reasoning.**
- A ~24-parameter logistic GLMM with two variance components is estimated to the
  **same accuracy** from a few 100 k rows as from millions — 4 M is statistical
  overkill for that many parameters.
- Subsampling **both classes at the same fraction preserves prevalence**, so the
  fitted intercept stays calibrated to the full test set and **no prior offset**
  is needed. This is deliberately *not* the negative-only downsampling knob
  (`negative_downsample_ratio`), which was rejected here: the negatives are the
  false-positive background the classifier must learn to reject *and* are the
  minority class (2.02 M neg vs 2.70 M pos), so thinning them would both hurt FP
  discrimination and worsen the imbalance.
- The cap is the **GLMM only**. `ranger`/`xgboost`/`glm` keep the full data (they
  are fast once threaded and must see every negative), and prediction is always
  on the full test set — so read-level classification and H1–H12 are unchanged in
  scope. This keeps **all five models and every hypothesis (including H9)**.

**Validation.** On a real held-out-donor `combined`-arm split: full train
4,163,510 rows → capped to exactly **250,000**, prevalence **0.5807 → 0.5807**
(unchanged); the fit ran in **7.9 min** (vs ~6–7 h uncapped) and returned 559,472
finite predictions with **AUPRC 0.9987** — i.e. no measurable loss of GLMM
quality. A standalone sweep gave AUPRC 0.998 @50 k and 0.999 @250 k. Test suite
**112 → 120 passed, 0 failed** (new `stratified_cap_rows` unit tests + a
glmmTMB-with-cap smoke test).

**Result.** Stage 05 projection **~1–2 weeks → ~1–1.5 days**, within the ~40 h
window, with the full model set and all hypotheses intact. Required a kill +
relaunch `--gt fixed --from 05` (stages 01–04 outputs reused); the ~3 days of
in-memory arm-1/2 progress were unavoidably discarded, but that run wrote no
predictions until the very end (~days away) so nothing usable was lost for the
deadline.
