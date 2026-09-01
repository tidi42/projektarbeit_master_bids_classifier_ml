# Zymo-in-human read-classification: ML hypothesis-testing pipeline

Distinguish true ZymoBIOMICS reads from false positives in a ~99 % human,
ONT adaptive-sampling matrix, and **test the pre-registered hypotheses
H1–H12** under a nested cross-validation design. Built from `project_plan.md`
(the hypotheses + guardrails) and `notes.txt` (issues A–L).

> **Things to read first**
> 1. `00_config.R` **SECTION 1** — every hardcoded file path lives here.
> 2. `00_config.R` **SECTION 5** — the *still-open items* that block a full run.
>    `Rscript scripts/run_pipeline.R --status` prints both, plus which stage each blocks.
> 3. [`required_input_files.md`](../required_input_files.md) (repo root) — full inventory of
>    every input file, reference/database, external tool and R package the pipeline needs,
>    with expected format/schema for each.

## Pipeline stages

| Stage | File | Does | Key notes |
|------|------|------|-----------|
| — | `00_config.R` | Hardcoded paths, parameters, feature blocks, hypothesis + **open-items** registry | all |
| — | `utils.R` | AUPRC, precision@recall, FDR, Holm–Šídák, bootstrap CI, Poisson floor, complexity | C, H, J |
| 01 | `01_external_tools.R` | QC + `minimap2` ground truth + **independent** BLASTn & Kraken2 arms | A, D, E, F |
| 02 | `02_ground_truth_labels.R` | 3-class labels (pos/neg/**ambiguous**) + **Poisson floor** + **leakage** correction | A, H, L |
| 03 | `03_build_features.R` | Taxon-**agnostic** feature table (no taxon identity) | B, D, F |
| 04 | `04_cv_splits.R` | Nested **LOEO** (8 folds) + **LOTO** splits, donor-grouped | B, guardrails |
| 05 | `05_train_models.R` | fixed-threshold / GLM / RF / XGBoost / **GLMM**; ablations; **no NN** | I, K, F |
| 06 | `06_evaluate.R` | AUPRC, P@recall=0.80, FDR — stratified by titration, ± truncated reads, read + **sample×taxon** | C, E, G, H12 |
| 07 | `07_hypothesis_tests.R` | H1–H6 paired Wilcoxon + **Holm–Šídák**; H9–H12 with CIs (H7 deferred, H8 dropped) | J, all |

## How the notes are enforced

- **A** ground truth via `minimap2` with a third *ambiguous* class, never circular.
- **B** taxon identity is excluded from the model matrix; `model_features()` strips it; LOTO arm added.
- **C** AUPRC / precision@recall=0.80 / FDR, stratified c1–c5 (aggregate is dominated by c5).
- **D** BLASTn and Kraken2 are scored as independent arms; neither gates the other.
- **E** `end_reason` kept; every metric reported with and without unblocked/truncated reads.
- **F** feature list per note F; residual adapter content is **not** a model feature (only used for H7).
- **G** sample×taxon calling is a co-primary endpoint, not optional.
- **H** Poisson `P(≥1 genome)` floor routes low-probability species×level cells to an *indeterminate* stratum.
- **I** GLMM via `glmmTMB` with `(1|donor)+(1|species)`; optional negative downsampling + logit prior-correction.
- **J** primary family capped at 6 (exact Wilcoxon floor at n=8 is 2/2⁸ = 0.0078); Holm–Šídák within family.
- **K** feed-forward NN dropped (H8 inactive).
- **L** cross-barcode leakage into low-mass c4/c5 estimated from the negative barcode and flagged.

## Run it

```bash
# status only: hardcoded-path gaps + open items + hypothesis list
Rscript scripts/run_pipeline.R --status

# full run -- BOTH ground-truth cutoff profiles, each into its own results folder
Rscript scripts/run_pipeline.R                 # == --gt both (default)

# a single ground-truth profile
Rscript scripts/run_pipeline.R --gt fixed      # run 1: id>=0.90, cov>=0.80
Rscript scripts/run_pipeline.R --gt calculated # run 2: cutoffs from the sensitivity sweep

# resume / subset (applies to whichever profile(s) run)
Rscript scripts/run_pipeline.R --from 03
Rscript scripts/run_pipeline.R --only 05,06,07
```

### Two ground-truth runs (`--gt`)
The whole supervised target hinges on the ground-truth identity/coverage cutoffs,
so the pipeline is run **twice**, once per cutoff **profile** (selected by the
`GT_PROFILE` env var, which `run_pipeline.R --gt` sets for you):

| Profile | Cutoffs | Output folder |
|---------|---------|---------------|
| `fixed` | `gt_min_identity >= 0.90`, `gt_min_coverage >= 0.80` (SECTION 3 of `00_config.R`) | `results/gt_fixed_id0.90_cov0.80/` |
| `calculated` | the pair **recommended** by `prepare_cutoff_sensitivity.R` (`work/cutoff_recommended.tsv`) | `results/gt_calculated_id<id>_cov<cov>/` |

`--gt both` runs `fixed` first (which produces the shared stage-01 alignments),
auto-builds the recommended cutoffs, then runs `calculated` starting at stage 02
(stage 01 is cutoff-independent and shared). Each run's cutoff-**dependent**
outputs — `labels.tsv`, `feature_table`, `cv_splits`, models, metrics, hypotheses
and a `ground_truth_settings.tsv` provenance file — land in that profile's own
`results/<tag>/` folder, so the two runs never overwrite each other. Standalone:
`GT_PROFILE=calculated Rscript scripts/02_ground_truth_labels.R`.

### Dependencies
- **External tools** (PATH): `minimap2`, `samtools`, `blastn`, `kraken2`, optionally `seqkit`.
- **R packages**: `data.table`, `ranger`, `xgboost`, `glmmTMB`, optionally `arrow`,
  `Biostrings`, `lmerTest`. Each stage fails with an actionable message if one is missing.

## Outputs (`results/<ground-truth-tag>/`)
Each ground-truth profile writes into its own `results/gt_<profile>_id<id>_cov<cov>/`:
- `ground_truth_settings.tsv` — the cutoffs that define this run's folder
- `labels.tsv`, `feature_table.parquet`, `cv_splits.rds` — cutoff-dependent intermediates
- `predictions.tsv.gz` — per-read scores for every arm × model × fold
- `metrics_read_level.tsv`, `metrics_sample_taxon.tsv`, `calibration.tsv`
- `hypothesis_tests.tsv` — H1–H12 with effect sizes, CIs, Wilcoxon + Holm–Šídák p-values

Shared, cutoff-**independent** artefacts stay in `work/` (per-library `minimap2`/BLAST/
Kraken2 output, `sample_taxon_coverage.tsv`, `expected_sample_taxon.tsv`, and the
`cutoff_sensitivity.tsv` / `cutoff_recommended.tsv` sweep) and are reused by both runs.

## Status: this is a wired skeleton
The control flow, statistics, and note-driven design decisions are complete and
run end-to-end **once real inputs exist**. Items still open (paths, cutoffs,
certificate of analysis, adapter detection, GLMM stack, etc.) are enumerated in
`00_config.R` SECTION 5 and surfaced by `--status`. Resolve them and flip each
`status` to `"resolved"`.
