## =============================================================================
## 07_hypothesis_tests.R  --  the pre-registered tests, H1-H12
## -----------------------------------------------------------------------------
## Primary family H1-H6: paired exact Wilcoxon across the 8 LOEO folds, effect
## sizes as median fold differences with 95% donor-level bootstrap CIs, and
## Holm-Sidak control WITHIN the family. With n = 8 the exact two-sided floor is
## 2/2^8 = 0.0078, so at most ~6 tests can clear alpha = 0.05 -- which is exactly
## why the family is capped at six. [note J / project_plan]
##
## Secondary H9-H12: reported with CIs, NO alpha spending. H7 (adapter content) is
## DEFERRED to a later experiment [2026-08-03] and H8 (NN) is dropped by note K;
## both emit a provenance row only (not run). See reporting_2. [OI T1 / note K]
## =============================================================================

suppressWarnings(suppressMessages(library(data.table)))
if (!exists("cfg")) {
  .sd <- local({
    a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1])))
    else "/home/tdinse/master_timo/projektarbeit/projektarbeit_ml/scripts"
  })
  source(file.path(.sd, "00_config.R")); source(file.path(.sd, "utils.R"))
}

## ---- metric accessors -------------------------------------------------------
mvec <- function(M, a_arm, a_model, metric = "auprc", a_level = "read",
                 a_stratum = "all", a_scheme = "LOEO", a_trunc = TRUE) {
  d <- M[arm == a_arm & model == a_model & level == a_level &
           stratum == a_stratum & scheme == a_scheme]
  if ("truncated_included" %in% names(d))
    d <- d[is.na(truncated_included) | truncated_included == a_trunc]
  if ("recall_target" %in% names(d))          # collapse the R9 multi-target rows
    d <- d[recall_target == cfg$params$recall_primary]
  setNames(d[[metric]], as.character(d$fold))
}

best_family <- function(M, arm = "combined",
                        fams = c("glm", "ranger_rf", "xgboost", "glmmTMB")) {
  fams <- intersect(fams, unique(M$model))
  means <- vapply(fams, function(f) mean(mvec(M, arm, f), na.rm = TRUE), numeric(1))
  if (all(is.na(means))) return("xgboost")
  names(means)[which.max(means)]
}

## [item 3 / F2] Leakage-free primary-model selection. Ranks families by the mean
## INNER-CV AUPRC written by stage 05 (results/inner_cv_scores.tsv) -- which never
## touches the outer test folds -- instead of the outer metrics `M` (that would be
## the F2 double-dip). Falls back to best_family(M) only if inner scores are
## absent (documented as second-best).
load_inner_scores <- function() {
  p <- file.path(cfg$paths$out_root, "inner_cv_scores.tsv")
  if (file.exists(p)) fread(p) else NULL
}
select_primary_model <- function(M = NULL, arm = "combined",
                                 fams = c("glm", "ranger_rf", "xgboost", "glmmTMB")) {
  ics <- load_inner_scores()
  if (!is.null(ics) && nrow(ics)) {
    a_arm <- arm
    agg <- ics[arm == a_arm & model %in% fams, .(m = mean(inner_auprc, na.rm = TRUE)), by = model]
    if (nrow(agg) && any(is.finite(agg$m))) return(agg$model[which.max(agg$m)])
  }
  if (!is.null(M)) return(best_family(M, arm, fams))
  if ("xgboost" %in% fams) "xgboost" else fams[1]
}

run_paired <- function(x, y, id, family, comparison, note = "") {
  f <- intersect(names(x), names(y))
  x <- x[f]; y <- y[f]; d <- x - y
  wt <- paired_wilcox(x, y, cfg$params$wilcox_exact)
  ci <- donor_bootstrap_ci(d, cfg$params$n_boot, seed = cfg$params$seed)
  data.table(id = id, family = family, comparison = comparison, n_folds = length(f),
             median_diff = unname(ci["est"]), ci_lo = unname(ci["lo"]),
             ci_hi = unname(ci["hi"]), wilcox_p = wt$p.value, note = note)
}

## ---- PRIMARY FAMILY ---------------------------------------------------------
test_H1 <- function(M) {  # combined > best single arm, model family fixed
  fam <- select_primary_model(M, "combined")
  singles <- c("blast_only", "kraken2_only")
  best_single <- singles[which.max(vapply(singles,
                  function(a) mean(mvec(M, a, fam), na.rm = TRUE), numeric(1)))]
  run_paired(mvec(M, "combined", fam), mvec(M, best_single, fam),
             "H1", "primary",
             sprintf("AUPRC combined(%s) - %s(%s)", fam, best_single, fam))
}

test_H2 <- function(M) {  # best ML > fixed-threshold baseline (same arm)
  fam <- select_primary_model(M, "combined")
  run_paired(mvec(M, "combined", fam), mvec(M, "combined", "fixed_threshold"),
             "H2", "primary",
             sprintf("AUPRC combined(%s) - combined(fixed_threshold)", fam))
}

test_H3 <- function(M) {  # ML gain larger at low than high titration (interaction)
  fam <- select_primary_model(M, "combined")
  ## titration order: c1/c2 = HIGHEST concentration (high abundance), c4/c5 = LOWEST
  ## (low abundance). H3 predicts a LARGER ML gain at LOW abundance -> low := c4/c5.
  ## (Same c1/c2 <-> c4/c5 inversion fixed in low_mass_levels, note L; was reversed here.)
  low <- c("c4", "c5"); high <- c("c1", "c2")
  d <- M[level == "read" & scheme == "LOEO" & arm == "combined" &
           model %in% c(fam, "fixed_threshold") &
           stratum %in% c(low, high) & (truncated_included %in% TRUE)]
  if ("recall_target" %in% names(d)) d <- d[recall_target == cfg$params$recall_primary]
  if (!nrow(d)) return(data.table(id = "H3", family = "primary",
        comparison = "interaction method x level", n_folds = 0L,
        median_diff = NA_real_, ci_lo = NA, ci_hi = NA, wilcox_p = NA_real_,
        note = "no stratified metrics"))
  d[, method := ifelse(model == "fixed_threshold", "fixed", "ML")]
  d[, level_group := ifelse(stratum %in% low, "low", "high")]
  p_int <- NA_real_; est_int <- NA_real_; how <- "lm"
  fit <- tryCatch({
    if (requireNamespace("lmerTest", quietly = TRUE)) {
      how <<- "lmerTest"
      lmerTest::lmer(auprc ~ method * level_group + (1 | fold), data = d)
    } else lm(auprc ~ method * level_group, data = d)
  }, error = function(e) NULL)
  if (!is.null(fit)) {
    co <- summary(fit)$coefficients
    ix <- grep(":", rownames(co))
    if (length(ix)) { est_int <- co[ix[1], 1]
      p_int <- co[ix[1], ncol(co)] }
  }
  data.table(id = "H3", family = "primary",
             comparison = "interaction method x level (AUPRC ~ method*level + (1|donor))",
             n_folds = uniqueN(d$fold), median_diff = est_int, ci_lo = NA_real_,
             ci_hi = NA_real_, wilcox_p = p_int, note = paste0("fit=", how))
}

test_H4 <- function(M) {  # tree ensemble > linear GLM (combined arm)
  ## [F2] pick the tree family by inner-CV (leakage-free); do NOT pmax the two
  ## outer-fold vectors (taking the max of two noisy estimators biases upward).
  tree <- select_primary_model(M, "combined", fams = c("ranger_rf", "xgboost"))
  run_paired(mvec(M, "combined", tree), mvec(M, "combined", "glm"), "H4", "primary",
             sprintf("AUPRC %s - GLM (combined); tree family chosen by inner-CV", tree))
}

test_H5 <- function(M) {  # margin+human(+breadth) drive the gain: full vs ablated
  fam <- "xgboost"
  ablated <- "combined_minus_H5_key"
  if (!ablated %in% unique(M$arm)) ablated <- "combined_minus_blast_margin"
  run_paired(mvec(M, "combined", fam), mvec(M, ablated, fam), "H5", "primary",
             sprintf("AUPRC combined(%s) - %s(%s)", fam, ablated, fam),
             note = "read-level margin+human ablation; breadth measured separately as H5b (sample x taxon) [F6]")
}

## [F6] H5 'breadth' component, measured at the sample x taxon level. genome_breadth
## + coverage_evenness cannot be dropped from a read-level arm, so stage 06 fits the
## sample x taxon aggregate score WITH vs WITHOUT them (sxt_full vs sxt_minus_breadth)
## under leave-one-donor-out CV; this pairs their per-fold AUPRC. Reported with a CI in
## the secondary family (no alpha spending) so the primary family stays capped at 6.
test_H5b <- function(M) {
  full <- mvec(M, "combined", "sxt_full",          a_level = "sample_taxon")
  abl  <- mvec(M, "combined", "sxt_minus_breadth", a_level = "sample_taxon")
  common <- intersect(names(full), names(abl))
  if (length(common) < 2L)
    return(data.table(id = "H5b", family = "secondary",
      comparison = "sample x taxon breadth ablation (full - minus breadth)",
      n_folds = length(common), median_diff = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
      wilcox_p = NA_real_,
      note = "breadth/evenness unavailable (no coverage table) -> not measured [F6]"))
  run_paired(full[common], abl[common], "H5b", "secondary",
             "AUPRC sample x taxon full - minus-breadth (genome_breadth+coverage_evenness)",
             note = "H5 'breadth' component, measured by sample x taxon ablation [F6]")
}

test_H6 <- function(M) {  # sample x taxon aggregation > read-level thresholding
  fam <- select_primary_model(M, "combined")
  sxt <- mvec(M, "combined", fam, a_level = "sample_taxon")
  read_thr <- mvec(M, "combined", "fixed_threshold", a_level = "read")
  run_paired(sxt, read_thr, "H6", "primary",
             sprintf("AUPRC sample_taxon(%s) - read fixed-threshold", fam))
}

## ---- SECONDARY (CIs only, no alpha spending) --------------------------------
## H7 DEFERRED [2026-08-03]: the adapter-content enrichment test is a separate
## experiment for later (see reporting_2). It is not run now -- emit a provenance
## row only, mirroring the dropped H8. Re-enable by restoring the Fisher/OR body,
## setting cfg$hypotheses active=TRUE for H7, and producing work/adapter_content.tsv.
test_H7 <- function() data.table(id = "H7", family = "secondary",
  comparison = "adapter+ enrichment in FP vs TP", n_folds = NA_integer_,
  median_diff = NA_real_, ci_lo = NA, ci_hi = NA, wilcox_p = NA_real_,
  note = "DEFERRED -- experiment for later (adapter content); not run [OI T1]")

test_H8 <- function() data.table(id = "H8", family = "secondary",
  comparison = "feed-forward NN vs XGBoost", n_folds = NA_integer_,
  median_diff = NA_real_, ci_lo = NA, ci_hi = NA, wilcox_p = NA_real_,
  note = "DROPPED per note K -- NN excluded from the comparison")

## H9 DECOMPOSED into two transparent terms + a species-SOURCE sensitivity axis.
##   core  : GLMM(1|donor) - GLM  -> does modelling donor variance in TRAINING aid
##           transfer to a held-out donor? (no species term -> the honest question)
##   H9b_* : GLMM(1|donor+1|species[source]) - GLMM(1|donor) -> how much the species
##           baseline adds, for source in {truth = minimap2, classifier = coarse rank}.
## Robustness: the sign of GLMM-GLM is checked across sources {none, truth, classifier};
## if it FLIPS, the species source IS the finding (recorded in the H9 note). [H9]
test_H9 <- function(M) {
  arm     <- cfg$params$h9_arm %||% "combined"
  glm_v   <- mvec(M, arm, "glm")
  donor_v <- mvec(M, arm, "glmmTMB_none")
  if (!length(donor_v))
    return(data.table(id = "H9", family = "secondary",
      comparison = sprintf("GLMM(1|donor) - GLM (%s, held-out donors)", arm),
      n_folds = 0L, median_diff = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
      wilcox_p = NA_real_,
      note = "H9 decomposition unavailable (glmmTMB donor-only variant not produced)"))
  rows <- list(run_paired(donor_v, glm_v, "H9", "secondary",
                 sprintf("AUPRC GLMM(1|donor) - GLM (%s, held-out donors)", arm),
                 note = "honest donor-transfer term; species contribution decomposed as H9b"))
  signs <- c(none = sign(median(donor_v - glm_v[names(donor_v)], na.rm = TRUE)))
  for (lab in setdiff(names(cfg$params$h9_species_sources), "none")) {
    full_v <- mvec(M, arm, paste0("glmmTMB_", lab))
    if (!length(full_v)) next
    dc <- intersect(names(full_v), names(donor_v))
    rows[[length(rows) + 1L]] <- run_paired(full_v[dc], donor_v[dc], paste0("H9b_", lab),
        "secondary",
        sprintf("AUPRC GLMM(1|donor+1|species[%s]) - GLMM(1|donor): species-term contribution", lab),
        note = sprintf("species source = %s", lab))
    gc <- intersect(names(full_v), names(glm_v))
    signs[lab] <- sign(median(full_v[gc] - glm_v[gc], na.rm = TRUE))
  }
  robust <- length(unique(signs[is.finite(signs) & signs != 0])) <= 1L
  rows[[1]]$note <- paste0(rows[[1]]$note,
    sprintf("; GLMM-GLM sign by source {%s} -> %s",
            paste(sprintf("%s:%+d", names(signs), as.integer(signs)), collapse = ", "),
            if (robust) "ROBUST" else "SOURCE-DEPENDENT (flips -> the finding)"))
  rbindlist(rows, fill = TRUE)
}

## Transparent side-by-side table for H9: per species source, the GLMM-vs-GLM gap
## (with donor-level bootstrap CI) and the species-term contribution over the
## donor-only model. -> results/<run>/h9_species_sensitivity.tsv. [H9]
run_h9_sensitivity <- function(M) {
  arm     <- cfg$params$h9_arm %||% "combined"
  glm_v   <- mvec(M, arm, "glm")
  donor_v <- mvec(M, arm, "glmmTMB_none")
  if (!length(donor_v) || !length(glm_v)) return(NULL)
  mk <- function(lab) {
    v <- if (lab == "none") donor_v else mvec(M, arm, paste0("glmmTMB_", lab))
    if (!length(v)) return(NULL)
    gc <- intersect(names(v), names(glm_v)); if (!length(gc)) return(NULL)
    ci <- donor_bootstrap_ci(v[gc] - glm_v[gc], cfg$params$n_boot, seed = cfg$params$seed)
    dc <- intersect(names(v), names(donor_v))
    data.table(species_source = lab,
               glmm_median_auprc = median(v[gc], na.rm = TRUE),
               glm_median_auprc  = median(glm_v[gc], na.rm = TRUE),
               glmm_minus_glm    = unname(ci["est"]),
               ci_lo = unname(ci["lo"]), ci_hi = unname(ci["hi"]),
               species_contribution = if (lab == "none") 0 else
                 median(v[dc] - donor_v[dc], na.rm = TRUE),
               n_folds = length(gc))
  }
  tab <- rbindlist(Filter(Negate(is.null),
           lapply(names(cfg$params$h9_species_sources), mk)), fill = TRUE)
  if (!nrow(tab)) return(NULL)
  fin <- tab$glmm_minus_glm[is.finite(tab$glmm_minus_glm)]
  tab[, glmm_beats_glm := glmm_minus_glm > 0]
  tab[, robust_sign := length(unique(sign(fin[fin != 0]))) <= 1L]
  tab[]
}

test_H10 <- function(M) {  # leave-one-taxon-out degradation vs LOEO
  loeo <- mvec(M, "combined", "xgboost", a_scheme = "LOEO")
  loto <- mvec(M, "combined", "xgboost", a_scheme = "LOTO")
  ci_l <- donor_bootstrap_ci(loeo, cfg$params$n_boot, seed = cfg$params$seed)
  ci_t <- donor_bootstrap_ci(loto, cfg$params$n_boot, seed = cfg$params$seed)
  data.table(id = "H10", family = "secondary",
             comparison = "AUPRC LOEO vs LOTO (median degradation)",
             n_folds = length(loto),
             median_diff = unname(ci_l["est"] - ci_t["est"]),
             ci_lo = unname(ci_t["lo"]), ci_hi = unname(ci_t["hi"]),
             wilcox_p = NA_real_,
             note = sprintf("LOEO med=%.3f, LOTO med=%.3f", ci_l["est"], ci_t["est"]))
}

test_H11 <- function(preds) {  # truncated/unblocked reads differ & drive FP
  d <- preds[arm == "combined" & model == "xgboost" & scheme == "LOEO" &
               !is.na(end_reason_unblock)]
  if (!nrow(d)) return(data.table(id = "H11", family = "secondary",
      comparison = "truncated read score shift & FP share", n_folds = NA_integer_,
      median_diff = NA_real_, ci_lo = NA, ci_hi = NA, wilcox_p = NA_real_,
      note = "BLOCKED: no end_reason (sequencing_summary) [OI 6]"))
  neg <- d[y == 0]
  wt <- suppressWarnings(wilcox.test(score ~ end_reason_unblock, data = neg))
  thr <- d[, stats::quantile(score, 0.75, na.rm = TRUE)]
  fp_share_trunc <- d[score > thr & y == 0, mean(end_reason_unblock == 1)]
  data.table(id = "H11", family = "secondary",
             comparison = "unblocked vs signal_positive negative-read scores",
             n_folds = NA_integer_,
             median_diff = neg[end_reason_unblock == 1, median(score, na.rm = TRUE)] -
                           neg[end_reason_unblock == 0, median(score, na.rm = TRUE)],
             ci_lo = NA_real_, ci_hi = NA_real_, wilcox_p = wt$p.value,
             note = sprintf("FP share truncated = %.3f", fp_share_trunc))
}

test_H12 <- function(cal, preds) {  # calibration: model better calibrated than the Platt-scaled baseline
  ## [item 5 / F4] both the model AND the fixed_threshold baseline are now
  ## CALIBRATED probabilities (stage 05 Platt-scales the baseline on train), so
  ## Brier is a fair like-for-like calibration comparison from the cal table.
  fam <- select_primary_model(NULL, "combined")
  bx <- cal[arm == "combined" & model == fam & scheme == "LOEO", setNames(brier, as.character(fold))]
  bf <- cal[arm == "combined" & model == "fixed_threshold" & scheme == "LOEO", setNames(brier, as.character(fold))]
  common <- intersect(names(bf), names(bx))
  if (!length(common)) return(data.table(id = "H12", family = "secondary",
      comparison = "Brier model vs Platt-scaled fixed threshold", n_folds = NA_integer_,
      median_diff = NA_real_, ci_lo = NA, ci_hi = NA, wilcox_p = NA_real_,
      note = "insufficient calibration data"))
  run_paired(bf[common], bx[common], "H12", "secondary",
             sprintf("Brier fixed_threshold(Platt) - Brier %s (positive = model better calibrated)", fam))
}

## [item 4 / F3] Complementary donor-clustered mixed-model inference for the key
## primary comparisons, using per-(donor x titration) metric rows (~40 units)
## rather than 8, to escape the exact-Wilcoxon sign floor. A (1|donor) random
## intercept keeps donor clustering honest. Written to
## results/mixed_model_supplement.tsv beside the Wilcoxon results.
mixed_pair <- function(M, id, armA, modelA, armB, modelB, label, a_level = "read") {
  d <- M[level == a_level & scheme == "LOEO" & stratum != "all" &
           ((arm == armA & model == modelA) | (arm == armB & model == modelB))]
  if ("recall_target" %in% names(d)) d <- d[recall_target == cfg$params$recall_primary]
  if (!nrow(d)) return(NULL)
  d <- copy(d)
  d[, method := ifelse(arm == armA & model == modelA, "A", "B")]
  d[, donor := fold]
  res <- mixed_effect_test(d[, .(auprc, method, donor)], "auprc")
  data.table(id = id, comparison = label, mixed_estimate = res$estimate,
             mixed_p = res$p.value, n_units = res$n, engine = res$engine)
}

run_mixed_supplement <- function(M) {
  fam  <- select_primary_model(M, "combined")
  tree <- select_primary_model(M, "combined", fams = c("ranger_rf", "xgboost"))
  singles <- c("blast_only", "kraken2_only")
  best_single <- singles[which.max(vapply(singles,
                  function(a) mean(mvec(M, a, fam), na.rm = TRUE), numeric(1)))]
  rbindlist(list(
    mixed_pair(M, "H1", "combined", fam, best_single, fam,
               sprintf("combined(%s) vs %s(%s)", fam, best_single, fam)),
    mixed_pair(M, "H2", "combined", fam, "combined", "fixed_threshold",
               sprintf("combined(%s) vs combined(fixed_threshold)", fam)),
    mixed_pair(M, "H4", "combined", tree, "combined", "glm",
               sprintf("combined(%s) vs combined(glm)", tree))
  ), fill = TRUE)
}

## ---- driver -----------------------------------------------------------------
run_stage07 <- function() {
  cfg_init_dirs()
  rl  <- fread(cfg$paths$metrics_read)
  sxt <- fread(cfg$paths$metrics_sxt)
  M   <- rbindlist(list(rl, sxt), use.names = TRUE, fill = TRUE)
  preds <- fread(file.path(cfg$paths$out_root, "predictions.tsv.gz"))
  cal   <- fread(file.path(cfg$paths$out_root, "calibration.tsv"))

  primary <- rbindlist(list(test_H1(M), test_H2(M), test_H3(M),
                            test_H4(M), test_H5(M), test_H6(M)), fill = TRUE)
  ## Holm-Sidak WITHIN the primary family only. [note J]
  primary[, holm_sidak_p := holm_sidak(wilcox_p)]
  primary[, significant := !is.na(holm_sidak_p) & holm_sidak_p < cfg$params$alpha]

  secondary <- rbindlist(list(test_H7(), test_H8(), test_H9(M),
                              test_H10(M), test_H11(preds), test_H12(cal, preds),
                              test_H5b(M)),
                         fill = TRUE)
  secondary[, `:=`(holm_sidak_p = NA_real_,
                   significant = !is.na(ci_lo) & !is.na(ci_hi) & (ci_lo > 0 | ci_hi < 0))]

  results <- rbindlist(list(primary, secondary), use.names = TRUE, fill = TRUE)
  fwrite(results, cfg$paths$hypotheses_out, sep = "\t")

  ## [item 4] complementary mixed-model inference (n=8 power supplement)
  mixed <- tryCatch(run_mixed_supplement(M), error = function(e) NULL)
  if (!is.null(mixed) && nrow(mixed)) {
    fwrite(mixed, file.path(cfg$paths$out_root, "mixed_model_supplement.tsv"), sep = "\t")
    message("  mixed-model supplement -> ", file.path(cfg$paths$out_root, "mixed_model_supplement.tsv"))
  }

  ## [H9] transparent glmmTMB decomposition + species-source sensitivity table
  h9s <- tryCatch(run_h9_sensitivity(M), error = function(e) NULL)
  if (!is.null(h9s) && nrow(h9s)) {
    fwrite(h9s, file.path(cfg$paths$out_root, "h9_species_sensitivity.tsv"), sep = "\t")
    message("  H9 species-source sensitivity -> ",
            file.path(cfg$paths$out_root, "h9_species_sensitivity.tsv"))
  }

  message("Stage 07 complete -> ", cfg$paths$hypotheses_out)
  print(results[, .(id, family, comparison, median_diff, wilcox_p, holm_sidak_p, significant)])
  invisible(results)
}

if (sys.nframe() == 0L) run_stage07()
