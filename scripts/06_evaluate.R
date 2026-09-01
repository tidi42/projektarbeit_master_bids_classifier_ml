## =============================================================================
## 06_evaluate.R  --  read-level and sample x taxon metrics
## -----------------------------------------------------------------------------
## Turns stage-05 predictions into per-fold metrics. Endpoints [note C]:
##   AUPRC, precision @ recall = 0.80, FDR @ that operating point.
## Reported:
##   * stratified by titration level c1..c5  (aggregate is dominated by c5) [note C]
##   * with AND without truncated/unblocked adaptive-sampling reads          [note E]
##   * at read level AND at sample x taxon level (co-primary)                [note G]
##   * plus Brier score for calibration                                      [H12]
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

## Metrics per group, ONE ROW PER clinical recall target [R9]. AUPRC is
## threshold-free so it is repeated across the target rows; precision & FDR are
## evaluated at each target.
metric_rows <- function(y, score, targets = cfg$params$recall_target) {
  ok <- is.finite(score) & !is.na(y); y <- y[ok]; score <- score[ok]
  np <- sum(y == 1); nn <- sum(y == 0)
  if (np == 0 || nn == 0)
    return(data.table(recall_target = targets, auprc = NA_real_,
                      prec_at_recall = NA_real_, fdr_at_recall = NA_real_,
                      n_pos = np, n_neg = nn))
  ap <- auprc(y, score)
  data.table(recall_target = targets, auprc = ap,
             prec_at_recall = vapply(targets, function(t) precision_at_recall(y, score, t), numeric(1)),
             fdr_at_recall  = vapply(targets, function(t) fdr_at_recall(y, score, t), numeric(1)),
             n_pos = np, n_neg = nn)
}

## Per-group metrics across (arm, model, scheme, fold, stratum, truncated flag).
read_level_metrics <- function(preds) {
  variants <- rbindlist(lapply(cfg$params$report_truncated, function(incl_trunc) {
    d <- if (incl_trunc) preds else preds[is.na(end_reason_unblock) | end_reason_unblock == 0]
    d <- copy(d); d[, truncated_included := incl_trunc]; d
  }))
  ## overall + per-titration strata [note C]
  by_all <- variants[, metric_rows(y, score),
                     by = .(arm, model, scheme, fold, truncated_included)][, stratum := "all"]
  by_lvl <- variants[, metric_rows(y, score),
                     by = .(arm, model, scheme, fold, truncated_included, titration_level)]
  by_lvl[, stratum := titration_level][, titration_level := NULL]
  rbindlist(list(by_all, by_lvl), use.names = TRUE, fill = TRUE)[, level := "read"][]
}

## ---- sample x taxon aggregation (co-primary) [note G] -----------------------
## A (library x species) call. Aggregate read scores into: n reads above a
## train-free quantile threshold, mean/max score. Genome breadth + coverage
## evenness are the strongest real-vs-artefact discriminators [note F] and should
## be joined here from a subject-coverage pass (samtools depth) when available. [OI 11]
##
## DEPTH NORMALISATION: read counts vary a lot across donors/barcodes (10x
## dilution series), so the raw count is also expressed per-million reads of the
## SAME library (n / library_total * 1e6). Score aggregates (mean/max) are
## already depth-robust, and per-fold ranking metrics do not mix libraries.
aggregate_sample_taxon <- function(preds) {
  ## [item 6 / F5] threshold computed PER (arm, model) from that model's pooled
  ## out-of-fold score distribution (many data points), NOT a single global
  ## quantile across incomparable score scales -- so n_reads_above_thr is
  ## comparable across models.
  preds <- copy(preds)
  preds[, thr_am := stats::quantile(score, 0.75, na.rm = TRUE), by = .(arm, model)]
  ## per-library sequencing depth (reads scored in this arm/model/fold slice)
  depth <- preds[, .(library_total = .N),
                 by = .(arm, model, scheme, fold, donor, titration_level)]
  sxt <- preds[!is.na(species), .(
    n_reads_above_thr = sum(score > thr_am, na.rm = TRUE),
    mean_score = mean(score, na.rm = TRUE),
    max_score  = max(score, na.rm = TRUE),
    ## observed presence: >=1 ground-truth read mapped to this taxon
    y_obs = as.integer(any(y == 1))
  ), by = .(arm, model, scheme, fold, donor, titration_level, species)]
  sxt <- merge(sxt, depth, by = c("arm", "model", "scheme", "fold", "donor", "titration_level"),
               all.x = TRUE)
  sxt[, reads_above_thr_per_million := n_reads_above_thr / pmax(library_total, 1L) * 1e6]

  ## sample x taxon TRUTH: prefer the a-priori expectation (concentration x
  ## rel_abundance) from stage 02 -- a taxon is truly present iff it clears the
  ## Poisson floor; a Zymo taxon detected where it is NOT expected (e.g. a
  ## negative barcode) is a leakage/contamination FALSE POSITIVE. Cells that are
  ## spiked but below the floor are 'indeterminate' and dropped. Falls back to
  ## observed presence if the expected table is absent. [note H / note G]
  exp_path <- file.path(cfg$paths$work_dir, "expected_sample_taxon.tsv")
  if (file.exists(exp_path)) {
    ex <- fread(exp_path)[, .(donor, titration_level, species, expected_cells,
                              p_detect, expected_present)]
    sxt <- merge(sxt, ex, by = c("donor", "titration_level", "species"), all.x = TRUE)
    indet <- !is.na(sxt$p_detect) & sxt$p_detect > 0 & sxt$p_detect < cfg$params$poisson_p_min
    sxt <- sxt[!indet]
    sxt[, y := as.integer(expected_present %in% TRUE)]
  } else {
    sxt[, y := y_obs]
  }

  ## optional breadth/evenness join (per library x species) if produced upstream
  cov_path <- file.path(cfg$paths$work_dir, "sample_taxon_coverage.tsv")
  if (file.exists(cov_path)) {
    cov <- fread(cov_path)  # donor, titration_level, species, genome_breadth, coverage_evenness
    sxt <- merge(sxt, cov, by = c("donor", "titration_level", "species"), all.x = TRUE)
  } else {
    message("  sample_taxon_coverage.tsv absent -> breadth/evenness NA; using score aggregates only. [OI 11]")
    sxt[, `:=`(genome_breadth = NA_real_, coverage_evenness = NA_real_)]
  }

  ## R6: honest recall universe -- add expected-present taxa that produced NO
  ## reads as FALSE NEGATIVES (score below every observed call), per (arm, model,
  ## scheme, fold) test library. Without this, a wholly-undetected pathogen never
  ## counts against recall. [note G / H6]
  if (file.exists(exp_path)) {
    ex_pos <- fread(exp_path)[expected_present == TRUE, .(donor, titration_level, species)]
    libs <- unique(preds[, .(arm, model, scheme, fold, donor, titration_level)])
    universe <- libs[ex_pos, on = .(donor, titration_level), allow.cartesian = TRUE, nomatch = 0L]
    missing <- universe[!sxt, on = .(arm, model, scheme, fold, donor, titration_level, species)]
    if (nrow(missing)) {
      fin <- sxt$max_score[is.finite(sxt$max_score)]
      sentinel <- if (length(fin)) min(fin) - 1 else -1
      missing[depth, on = .(arm, model, scheme, fold, donor, titration_level),
              library_total := i.library_total]
      missing[, `:=`(n_reads_above_thr = 0L, mean_score = sentinel, max_score = sentinel,
                     reads_above_thr_per_million = 0, y_obs = 0L, y = 1L,
                     genome_breadth = 0, coverage_evenness = 0)]
      sxt <- rbindlist(list(sxt, missing), use.names = TRUE, fill = TRUE)
    }
  }
  sxt
}

sample_taxon_metrics <- function(preds) {
  sxt <- aggregate_sample_taxon(preds)
  ## score for the call = max read score (swap for a fitted aggregate if desired)
  by_all <- sxt[, metric_rows(y, max_score), by = .(arm, model, scheme, fold)][, stratum := "all"]
  by_lvl <- sxt[, metric_rows(y, max_score), by = .(arm, model, scheme, fold, titration_level)]
  by_lvl[, stratum := titration_level][, titration_level := NULL]
  rbindlist(list(by_all, by_lvl), use.names = TRUE, fill = TRUE)[
    , `:=`(level = "sample_taxon", truncated_included = TRUE)][]
}

## ---- sample x taxon breadth ablation (H5 'breadth' component) [F6] -----------
## H5 claims genome-breadth + coverage-evenness drive the gain, but those blocks
## live ONLY at the sample x taxon level, so NO read-level ablation can remove them
## (stage 05 correctly skips any ablation that drops no read-level feature). Here we
## MEASURE their contribution: aggregate the (arm, model) out-of-fold read scores to
## sample x taxon cells (joining breadth/evenness), then fit the sample x taxon
## aggregate score WITH vs WITHOUT the breadth block under leave-one-DONOR-out CV,
## reporting per-fold AUPRC for both. Emitted as pseudo-models 'sxt_full' and
## 'sxt_minus_breadth' at level='sample_taxon' so stage 07 runs the paired H5-breadth
## test with the same machinery. Leakage-free: the read scores are already OOF and the
## held-out donor never appears in the aggregator's training rows. [F6 / note F,G]
sample_taxon_ablation <- function(preds, arm = "combined", model = NULL) {
  a_arm <- arm
  avail <- unique(preds[arm == a_arm & scheme == "LOEO", model])
  a_model <- if (!is.null(model)) model else {
    pref <- intersect(c("xgboost", "ranger_rf", "glm", "glmmTMB"), avail)
    if (length(pref)) pref[1] else if (length(avail)) avail[1] else NA_character_
  }
  if (is.na(a_model)) return(data.table())
  d <- preds[arm == a_arm & model == a_model & scheme == "LOEO"]
  if (!nrow(d)) return(data.table())
  sxt <- aggregate_sample_taxon(d)
  breadth_feats <- intersect(cfg$sxt_ablation_sets$breadth, names(sxt))
  base_feats    <- intersect(cfg$sxt_score_features, names(sxt))
  ## breadth must be genuinely present (real coverage table joined upstream), not
  ## all-NA / all-zero sentinels -- otherwise the ablation is not measurable.
  bvals <- if (length(breadth_feats)) unlist(sxt[, ..breadth_feats]) else numeric(0)
  if (!length(breadth_feats) || !length(base_feats) ||
      sum(is.finite(bvals) & bvals != 0) < 5L) {
    message("  sample_taxon_ablation: breadth/evenness unavailable (no coverage table) -> skipped. [F6]")
    return(data.table())
  }
  full_feats <- c(base_feats, breadth_feats)

  ## leave-one-donor-out OOF scores from a logistic aggregator over the sxt features
  fit_predict <- function(feats) {
    donors <- sort(unique(sxt$donor))
    rbindlist(lapply(donors, function(dd) {
      tr <- sxt[donor != dd]; te <- sxt[donor == dd]
      if (!nrow(te) || uniqueN(tr$y) < 2L) return(NULL)
      med <- lapply(tr[, ..feats], function(v) { m <- suppressWarnings(median(v, na.rm = TRUE))
        if (is.finite(m)) m else 0 })
      impute <- function(dt) { X <- copy(dt[, ..feats])
        for (f in feats) { v <- X[[f]]; v[!is.finite(v)] <- med[[f]]; X[[f]] <- v }; X }
      Xtr <- impute(tr); Xte <- impute(te)
      fit <- tryCatch(suppressWarnings(glm(y ~ ., data = cbind(y = tr$y, Xtr),
                                           family = binomial())), error = function(e) NULL)
      sc <- if (is.null(fit)) rowSums(scale(as.matrix(Xte)))   # standardized-sum fallback
            else as.numeric(stats::predict(fit, newdata = Xte, type = "response"))
      data.table(fold = te$fold, y = te$y, score = sc)
    }), fill = TRUE)
  }
  mk_rows <- function(feats, mlabel) {
    oof <- fit_predict(feats)
    if (is.null(oof) || !nrow(oof)) return(data.table())
    oof[, metric_rows(y, score), by = .(fold)][
      , `:=`(arm = a_arm, model = mlabel, scheme = "LOEO", stratum = "all",
             level = "sample_taxon", truncated_included = TRUE)][]
  }
  rbindlist(list(mk_rows(full_feats, "sxt_full"),
                 mk_rows(base_feats, "sxt_minus_breadth")), use.names = TRUE, fill = TRUE)
}

## ---- calibration [H12] ------------------------------------------------------
calibration_metrics <- function(preds) {
  ## [item 5] fixed_threshold is now Platt-calibrated in stage 05, so EVERY model
  ## emits a probability in [0,1] and is Brier-comparable on equal footing.
  d <- preds[is.finite(score)]
  d[, .(brier = brier_score(y, pmin(pmax(score, 0), 1)),
        n = .N), by = .(arm, model, scheme, fold)]
}

## ---- descriptive multi-faceted model comparison [item 3] --------------------
## Per (arm, model) AUPRC / precision@recall / FDR at three facets: OVERALL,
## per-SPECIES (taxon-specific fit) and per-TITRATION level (concentration-
## specific fit), pooled across LOEO folds for stable estimates. A leaderboard to
## INTERPRET model behaviour and to pre-register the primary model -- it is NOT
## used to pick-then-test on the same folds (that would reintroduce F2).
model_comparison <- function(preds) {
  d <- preds[scheme == "LOEO" & is.finite(score)]
  if (!nrow(d)) return(data.table())
  overall <- d[, metric_rows(y, score), by = .(arm, model)][
    , `:=`(facet = "overall", key = "all")]
  bytax <- d[!is.na(species), metric_rows(y, score), by = .(arm, model, species)][
    , `:=`(facet = "taxon", key = species)][, species := NULL]
  bylvl <- d[, metric_rows(y, score), by = .(arm, model, titration_level)][
    , `:=`(facet = "concentration", key = titration_level)][, titration_level := NULL]
  rbindlist(list(overall, bytax, bylvl), use.names = TRUE, fill = TRUE)
}

run_stage06 <- function() {
  cfg_init_dirs()
  pred_path <- file.path(cfg$paths$out_root, "predictions.tsv.gz")
  preds <- fread(pred_path)

  rl <- read_level_metrics(preds)
  fwrite(rl, cfg$paths$metrics_read, sep = "\t")

  ## sample x taxon metrics + the H5 breadth ablation (sxt_full / sxt_minus_breadth) [F6]
  sxt <- sample_taxon_metrics(preds)
  abl <- tryCatch(sample_taxon_ablation(preds), error = function(e) {
    message("  sample_taxon_ablation failed: ", conditionMessage(e)); data.table() })
  sxt <- rbindlist(list(sxt, abl), use.names = TRUE, fill = TRUE)
  fwrite(sxt, cfg$paths$metrics_sxt, sep = "\t")

  cal <- calibration_metrics(preds)
  fwrite(cal, file.path(cfg$paths$out_root, "calibration.tsv"), sep = "\t")

  mc <- model_comparison(preds)
  fwrite(mc, file.path(cfg$paths$out_root, "model_comparison.tsv"), sep = "\t")

  message("Stage 06 complete:")
  message("  read-level     -> ", cfg$paths$metrics_read)
  message("  sample x taxon -> ", cfg$paths$metrics_sxt, "  (incl. F6 breadth ablation: sxt_full / sxt_minus_breadth)")
  message("  calibration    -> ", file.path(cfg$paths$out_root, "calibration.tsv"))
  message("  model compare  -> ", file.path(cfg$paths$out_root, "model_comparison.tsv"))
  invisible(list(read = rl, sxt = sxt, cal = cal))
}

if (sys.nframe() == 0L) run_stage06()
