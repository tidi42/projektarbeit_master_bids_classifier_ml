## =============================================================================
## 05_train_models.R  --  baselines + supervised models under nested CV
## -----------------------------------------------------------------------------
## For every (classifier arm x model x outer fold) trains on the fold's training
## reads and scores its held-out reads. Writes one long predictions table that
## stage 06 turns into metrics and stage 07 into hypothesis tests.
##
## Models [note K excludes the feed-forward NN]:
##   fixed_threshold  best single raw score, ORIENTATION + threshold tuned on
##                    TRAIN ONLY  -> the honest baseline for H2 [guardrail]
##   glm              logistic regression (linear)                    -> H4 ref
##   ranger_rf        probability random forest, mtry/node tuned inner-CV
##   xgboost          gradient boosting, depth/eta/nrounds tuned inner-CV
##   glmmTMB          logistic GLMM, (1|donor)+(1|species) random intercepts,
##                    optional negative downsampling + logit prior-correction [note I]
##
## Arms are the INDEPENDENT classifier feature sets (blast_only / kraken2_only /
## combined) [note D]. The H5 ablation drops feature blocks from 'combined'.
## Only OUTER-fold predictions are kept [guardrail].
## =============================================================================

suppressWarnings(suppressMessages(library(data.table)))
.sd <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1])))
  else "/home/tdinse/master_timo/projektarbeit/projektarbeit_ml/scripts"
})
if (!exists("cfg")) { source(file.path(.sd, "00_config.R")); source(file.path(.sd, "utils.R")) }
## stage 05 depends on model_features() (stage 03) and fold_masks() (stage 04)
if (!exists("model_features")) source(file.path(.sd, "03_build_features.R"))
if (!exists("fold_masks"))     source(file.path(.sd, "04_cv_splits.R"))

## ---- data helpers -----------------------------------------------------------
load_feature_table <- function() {
  p <- cfg$paths$feature_table
  if (file.exists(p) && requireNamespace("arrow", quietly = TRUE))
    return(as.data.table(arrow::read_parquet(p)))
  alt <- sub("\\.parquet$", ".tsv.gz", p)
  if (file.exists(alt)) return(fread(alt))
  stop("Feature table not found. Run stage 03 first.", call. = FALSE)
}

titration_num <- function(x) match(x, c("c1", "c2", "c3", "c4", "c5"))

## Median imputation learned on train, reused on test. Tree models can keep NA.
make_imputer <- function(train_X) {
  med <- vapply(train_X, function(v) suppressWarnings(median(v, na.rm = TRUE)), numeric(1))
  med[is.na(med)] <- 0
  med
}
apply_imputer <- function(X, med) {
  for (j in names(X)) X[[j]][is.na(X[[j]])] <- med[[j]]
  X
}

## ---- models: each returns test scores in [higher = more positive] -----------

fit_fixed_threshold <- function(tr, te, features) {
  y <- tr$y
  best <- list(feat = NA, sign = 1, auprc = -Inf)
  for (f in features) {
    v <- tr$X[[f]]
    if (all(is.na(v)) || length(unique(v[!is.na(v)])) < 2) next
    v[is.na(v)] <- median(v, na.rm = TRUE)
    for (s in c(1, -1)) {
      a <- auprc(y, s * v)
      if (!is.na(a) && a > best$auprc) best <- list(feat = f, sign = s, auprc = a)
    }
  }
  if (is.na(best$feat)) return(rep(0.5, nrow(te$X)))
  med    <- median(tr$X[[best$feat]], na.rm = TRUE)
  tr_raw <- best$sign * ifelse(is.na(tr$X[[best$feat]]), med, tr$X[[best$feat]])
  te_raw <- best$sign * ifelse(is.na(te$X[[best$feat]]), med, te$X[[best$feat]])
  ## [item 5 / F4] Platt scaling: fit a logistic map raw-score -> probability on
  ## TRAIN ONLY, so the baseline emits a CALIBRATED probability for a fair Brier
  ## comparison. The map is monotone, so AUPRC / precision@recall are unchanged.
  cal <- tryCatch(stats::glm(y ~ tr_raw, family = binomial()), error = function(e) NULL)
  sc <- if (is.null(cal)) {
    rng <- range(c(tr_raw, te_raw), finite = TRUE)
    (te_raw - rng[1]) / max(rng[2] - rng[1], 1e-9)
  } else {
    as.numeric(stats::predict(cal, newdata = data.frame(tr_raw = te_raw), type = "response"))
  }
  attr(sc, "chosen_feature") <- best$feat
  sc
}

fit_glm <- function(tr, te, features) {
  med <- make_imputer(tr$X[, ..features])
  Xtr <- apply_imputer(copy(tr$X[, ..features]), med)
  Xte <- apply_imputer(copy(te$X[, ..features]), med)
  df <- cbind(y = tr$y, Xtr)
  fit <- suppressWarnings(glm(y ~ ., data = df, family = binomial(), weights = tr$w))
  as.numeric(predict(fit, newdata = Xte, type = "response"))
}

fit_ranger <- function(tr, te, features, inner) {
  require_pkgs("ranger")
  med <- make_imputer(tr$X[, ..features])
  Xtr <- apply_imputer(copy(tr$X[, ..features]), med)
  Xte <- apply_imputer(copy(te$X[, ..features]), med)
  if (is.null(inner) || !length(inner)) {
    best <- list(mtry = max(1L, floor(0.5 * length(features))), min.node.size = 10)  # default HPs (selection pass)
  } else {
    ## small inner-CV grid on mtry (report outer only) [guardrail]
    grid <- expand.grid(mtry = unique(pmax(1, floor(c(0.33, 0.66) * length(features)))),
                        min.node.size = c(5, 20))
    best <- tune_inner(tr, features, inner, function(a, b, row) {
      m <- ranger::ranger(x = apply_imputer(copy(a$X[, ..features]), med), y = factor(a$y),
                          probability = TRUE, num.trees = 100,
                          mtry = row$mtry, min.node.size = row$min.node.size,
                          num.threads = cfg$params$train_threads)
      predict(m, apply_imputer(copy(b$X[, ..features]), med))$predictions[, "1"]
    }, grid)
  }
  m <- ranger::ranger(x = Xtr, y = factor(tr$y), probability = TRUE, num.trees = 500,
                      mtry = best$mtry, min.node.size = best$min.node.size,
                      case.weights = tr$w, num.threads = cfg$params$train_threads)
  predict(m, Xte)$predictions[, "1"]
}

fit_xgboost <- function(tr, te, features, inner) {
  require_pkgs("xgboost")
  Xtr <- as.matrix(tr$X[, ..features]); Xte <- as.matrix(te$X[, ..features])  # NA kept
  spw <- sum(tr$y == 0) / max(sum(tr$y == 1), 1)   # class imbalance weight
  if (is.null(inner) || !length(inner)) {
    best <- list(max_depth = 6, eta = 0.1, nrounds = 300)          # default HPs (selection pass)
  } else {
    grid <- expand.grid(max_depth = c(4, 6), eta = c(0.1, 0.3), nrounds = c(200, 400))
    best <- tune_inner(tr, features, inner, function(a, b, row) {
      da <- xgboost::xgb.DMatrix(as.matrix(a$X[, ..features]), label = a$y, missing = NA)
      m <- xgboost::xgb.train(params = list(objective = "binary:logistic",
                                            max_depth = row$max_depth, eta = row$eta,
                                            scale_pos_weight = spw, tree_method = "hist",
                                            nthread = cfg$params$train_threads),
                              data = da, nrounds = row$nrounds, verbose = 0)
      predict(m, xgboost::xgb.DMatrix(as.matrix(b$X[, ..features]), missing = NA))
    }, grid)
  }
  d <- xgboost::xgb.DMatrix(Xtr, label = tr$y, weight = tr$w, missing = NA)
  m <- xgboost::xgb.train(params = list(objective = "binary:logistic",
                                        max_depth = best$max_depth, eta = best$eta,
                                        scale_pos_weight = spw, tree_method = "hist",
                                        nthread = cfg$params$train_threads),
                          data = d, nrounds = best$nrounds, verbose = 0)
  predict(m, xgboost::xgb.DMatrix(Xte, missing = NA))
}

## Label-stratified row cap. Returns the row indices to KEEP so nrow(kept) <= cap
## while preserving each class's share (same sampling fraction per label) -> the
## class ratio (prevalence) is unchanged, so a GLMM fit on the subset keeps a
## calibrated intercept with NO prior offset. Deterministic (seeded). Uses
## sample.int on the within-class index vector (avoids sample()'s length-1 trap).
## Stage-05 GLMM SPEED lever; ranger/xgboost/glm are NOT subsampled. [perf 2026-08-15]
stratified_cap_rows <- function(y, cap, seed = cfg$params$seed) {
  n <- length(y)
  if (is.null(cap) || is.na(cap) || n <= cap) return(seq_len(n))
  frac <- cap / n
  set.seed(seed)
  keep <- integer(0)
  for (lv in unique(y)) {
    idx <- which(y == lv)
    k   <- round(frac * length(idx))
    keep <- c(keep, if (k >= length(idx)) idx else idx[sample.int(length(idx), k)])
  }
  sort(keep)
}

## Logistic GLMM with random intercepts. [note I / H9]
## glmm_fit_predict() is the general fitter used BOTH as the `glmmTMB` model family
## and by the H9 decomposition. sp_train/sp_test choose the species grouping:
##   NULL              -> (1|donor) only. At a held-out (LOEO) donor the donor RE is
##                        a NEW level, so prediction is population-level -- exactly the
##                        HONEST "does modelling donor variance in training aid transfer?"
##   character vector  -> (1|donor)+(1|species); species can recur, so its RE IS applied
##                        at predict (re.form = ~(1|species)) while donor is integrated out.
## Donor-depth imbalance is absorbed by (1|donor), so balance_unit is NOT applied here.
glmm_fit_predict <- function(tr, te, features, sp_train = NULL, sp_test = NULL) {
  require_pkgs("glmmTMB")
  use_sp <- !is.null(sp_train)
  med  <- make_imputer(tr$X[, ..features])
  d_tr <- cbind(y = tr$y, apply_imputer(copy(tr$X[, ..features]), med), donor = tr$grp$donor)
  if (use_sp) d_tr$species <- as.character(sp_train)
  ## optional negative downsampling + logit prior-correction for scale [note I / OI 16]
  s <- 1
  if (!is.na(cfg$params$negative_downsample_ratio)) {
    pos <- which(d_tr$y == 1); neg <- which(d_tr$y == 0)
    keep_neg <- min(length(neg), round(cfg$params$negative_downsample_ratio * length(pos)))
    s <- keep_neg / length(neg)
    set.seed(cfg$params$seed)
    d_tr <- d_tr[c(pos, sample(neg, keep_neg)), ]
  }
  ## GLMM row cap (SPEED): when the training fold exceeds cfg$params$glmm_max_train_rows,
  ## subsample it label-stratified at the NATURAL class ratio. Prevalence is preserved,
  ## so the intercept stays calibrated and no prior offset is applied (s unchanged). A
  ## ~24-parameter GLMM fits to the same accuracy from a few 100k rows as from millions
  ## (AUPRC 0.999 @250k vs full 4M) for ~50x less time. GLMM-only; trees/glm keep full data.
  if (nrow(d_tr) > (cfg$params$glmm_max_train_rows %||% Inf))
    d_tr <- d_tr[stratified_cap_rows(d_tr$y, cfg$params$glmm_max_train_rows), ]
  re   <- if (use_sp) "+ (1 | donor) + (1 | species)" else "+ (1 | donor)"
  form <- as.formula(paste("y ~", paste(features, collapse = " + "), re))
  fit  <- glmmTMB::glmmTMB(form, data = d_tr, family = binomial())
  d_te <- cbind(apply_imputer(copy(te$X[, ..features]), med), donor = te$grp$donor)
  if (use_sp) d_te$species <- as.character(sp_test)
  ## glmmTMB::predict allows re.form only NULL / NA / ~0 (NOT lme4-style ~(1|species)).
  ## NULL conditions on ALL REs, but the held-out LOEO donor is a NEW level, so
  ## allow.new.levels sets its donor RE to 0 (population) while the trained species
  ## REs ARE applied -- exactly the intended "integrate out donor, apply species".
  reform <- if (use_sp) NULL else NA   # donor unseen at LOEO -> donor RE = 0 (population)
  eta <- predict(fit, newdata = d_te, type = "link", re.form = reform, allow.new.levels = TRUE)
  plogis(eta + log(s))   # prior correction if negatives were downsampled
}

## The `glmmTMB` model family = the full (1|donor)+(1|species truth) model. [note I]
fit_glmmTMB <- function(tr, te, features)
  glmm_fit_predict(tr, te, features, sp_train = tr$grp$species, sp_test = te$grp$species)

## Generic inner-CV grid search over donor-grouped inner folds; returns best row.
## [M2] The grid row is passed EXPLICITLY as the 3rd arg of predict_fun(a, b, row);
## there is no global `grid_row` and no hidden coupling to the closure's frame.
tune_inner <- function(tr, features, inner, predict_fun, grid) {
  scores <- numeric(nrow(grid))
  for (r in seq_len(nrow(grid))) {
    row <- grid[r, , drop = FALSE]
    fold_auprc <- vapply(inner, function(val_donors) {
      is_val <- tr$grp$donor %in% val_donors
      if (!any(is_val) || !any(!is_val)) return(NA_real_)
      a <- list(X = tr$X[!is_val], y = tr$y[!is_val], grp = tr$grp[!is_val])
      b <- list(X = tr$X[is_val],  y = tr$y[is_val],  grp = tr$grp[is_val])
      auprc(b$y, predict_fun(a, b, row))
    }, numeric(1))
    scores[r] <- mean(fold_auprc, na.rm = TRUE)
  }
  grid[which.max(scores), , drop = FALSE]
}

## ---- driver -----------------------------------------------------------------
## Per-read training weights so each donor (or library) contributes EQUALLY to
## the pooled fit despite very different read depths (10x dilution series). This
## is the training-side read-count normalisation; controlled by
## cfg$params$balance_unit = "none" | "library" | "donor".
balance_weights <- function(grp, unit = cfg$params$balance_unit %||% "none") {
  n <- nrow(grp)
  if (is.null(unit) || unit == "none") return(rep(1, n))
  key <- switch(unit,
    donor   = grp$donor,
    library = paste(grp$donor, grp$titration_level),
    stop("unknown balance_unit: ", unit))
  cnt <- stats::ave(seq_len(n), key, FUN = length)  # reads per unit in this pool
  w <- 1 / cnt                                       # equal TOTAL weight per unit
  w * n / sum(w)                                      # normalise to mean weight 1
}

split_xy <- function(dt, features) {
  grp <- dt[, .(donor, species, titration_level)]
  list(X = dt[, ..features],
       y = as.integer(dt$label == "positive"),
       grp = grp,
       w = balance_weights(grp))
}

predict_arm_model <- function(model, tr, te, features, inner) {
  switch(model,
    fixed_threshold = fit_fixed_threshold(tr, te, features),
    glm             = fit_glm(tr, te, features),
    ranger_rf       = fit_ranger(tr, te, features, inner),
    xgboost         = fit_xgboost(tr, te, features, inner),
    glmmTMB         = fit_glmmTMB(tr, te, features),
    stop("unknown model: ", model))
}

## [item 3 / F2] Leakage-free model-selection signal: mean AUPRC of `model` across
## the inner donor folds of the training set, using DEFAULT hyperparameters (no
## nested tuning) so families can be ranked cheaply WITHOUT ever touching outer
## test reads. Stage 07's select_primary_model() aggregates these to pick the
## primary family for H1/H2/H4/H6 instead of the (double-dipping) outer metrics.
inner_cv_score <- function(model, tr, features, inner) {
  if (is.null(inner) || !length(inner)) return(NA_real_)
  sc <- vapply(inner, function(val_donors) {
    is_val <- tr$grp$donor %in% val_donors
    if (!any(is_val) || !any(!is_val)) return(NA_real_)
    a <- list(X = tr$X[!is_val], y = tr$y[!is_val], grp = tr$grp[!is_val], w = tr$w[!is_val])
    b <- list(X = tr$X[is_val],  y = tr$y[is_val],  grp = tr$grp[is_val],  w = tr$w[is_val])
    s <- tryCatch(predict_arm_model(model, a, b, features, inner = NULL),
                  error = function(e) rep(NA_real_, length(b$y)))
    auprc(b$y, s)
  }, numeric(1))
  mean(sc, na.rm = TRUE)
}

## [imputation caveat / F-missing] Read-only diagnostic: how much median-imputation
## the non-tree models incur, and -- critically -- which feature columns are
## ENTIRELY missing within a TRAINING fold (median falls back to 0 -> a constant,
## zero-information column that is silently inert). Lets us judge whether the
## imputation caveat is negligible. (M1/M3 fixed: the formerly-inert
## subject_assembly_level / lca_rank / k2_lca_rank columns are now dropped, so
## remaining NAs are GENUINE per-read missingness, e.g. k2_* on a BLAST-only read.)
## Does NOT touch training; failures here never abort the run.
feature_missingness_report <- function(ft, splits, arms = names(cfg$classifier_arms)) {
  if (is.null(splits$loeo) || !length(splits$loeo)) return(invisible(NULL))
  rows <- list(); i <- 0L
  for (arm in arms) {
    feats <- intersect(model_features(arm), names(ft))
    if (!length(feats)) next
    for (fold in splits$loeo) {
      m   <- fold_masks(ft, fold)
      Xtr <- ft[m$train, ..feats]
      ntr <- nrow(Xtr)
      for (f in feats) {
        n_na <- sum(is.na(Xtr[[f]]))
        i <- i + 1L
        rows[[i]] <- data.table(
          arm = arm, scheme = "LOEO", fold = fold$fold, feature = f,
          n_train = ntr, n_missing = n_na,
          pct_missing = round(100 * n_na / max(ntr, 1L), 2),
          all_na = n_na == ntr)   # -> imputed to constant 0: zero-variance, inert
      }
    }
  }
  if (!length(rows)) return(invisible(NULL))
  rep  <- rbindlist(rows)
  path <- file.path(cfg$paths$out_root, "feature_missingness.tsv")
  fwrite(rep, path, sep = "\t")
  ## columns that are entirely NA in one or more folds (silently dead features)
  allna <- rep[all_na == TRUE, .(folds_all_na = .N), by = feature][order(-folds_all_na)]
  if (nrow(allna))
    message("  WARNING: ", nrow(allna), " feature(s) entirely NA in >=1 fold ",
            "(imputed to constant 0, no signal): ",
            paste(sprintf("%s[%d]", allna$feature, allna$folds_all_na), collapse = ", "))
  worst <- rep[pct_missing > 0, .(mean_pct = round(mean(pct_missing), 1)), by = feature
             ][order(-mean_pct)][seq_len(min(5L, .N))]
  message("  feature missingness -> ", path,
          if (nrow(worst)) paste0(" (top imputed: ",
            paste(sprintf("%s %.1f%%", worst$feature, worst$mean_pct), collapse = ", "), ")")
          else " (no imputation needed)")
  invisible(rep)
}

run_stage05 <- function(models = cfg$params$models,
                        arms = names(cfg$classifier_arms),
                        do_ablation = TRUE, do_loto = TRUE, do_h9 = TRUE) {
  cfg_init_dirs()
  ft <- load_feature_table()
  splits <- readRDS(cfg$paths$cv_splits)
  ## imputation/missingness diagnostic (read-only; never aborts the run)
  tryCatch(feature_missingness_report(ft, splits),
           error = function(e) message("  feature_missingness_report skipped: ", conditionMessage(e)))
  preds <- list(); k <- 0L
  inner_scores <- list(); ki <- 0L

  ## ---- primary: arms x models x LOEO folds --------------------------------
  for (arm in arms) {
    message(sprintf("-- training arm %d/%d: %s", match(arm, arms), length(arms), arm))
    feats <- intersect(model_features(arm), names(ft))
    for (fold in splits$loeo) {
      m <- fold_masks(ft, fold)
      tr <- split_xy(ft[m$train], feats); te_dt <- ft[m$test]; te <- split_xy(te_dt, feats)
      for (model in models) {
        sc <- tryCatch(predict_arm_model(model, tr, te, feats, fold$inner),
                       error = function(e) { message("  ", arm, "/", model,
                         " fold ", fold$fold, " failed: ", conditionMessage(e)); rep(NA_real_, nrow(te_dt)) })
        ## [item 3] record the leakage-free inner-CV score for the combined arm
        if (identical(arm, "combined")) {
          isc <- tryCatch(inner_cv_score(model, tr, feats, fold$inner), error = function(e) NA_real_)
          ki <- ki + 1L
          inner_scores[[ki]] <- data.table(arm = arm, model = model, fold = fold$fold, inner_auprc = isc)
        }
        k <- k + 1L
        preds[[k]] <- data.table(read_id = te_dt$read_id, donor = te_dt$donor,
          titration_level = te_dt$titration_level, species = te_dt$species,
          y = as.integer(te_dt$label == "positive"),
          end_reason_unblock = te_dt$end_reason_unblock,
          arm = arm, model = model, scheme = "LOEO", fold = fold$fold, score = as.numeric(sc))
      }
    }
  }
  if (length(inner_scores))
    fwrite(rbindlist(inner_scores), file.path(cfg$paths$out_root, "inner_cv_scores.tsv"), sep = "\t")

  ## ---- H5 ablation: drop each block set from 'combined' (xgboost) ----------
  if (do_ablation) {
    base <- intersect(model_features("combined"), names(ft))
    for (setname in names(cfg$ablation_sets)) {
      drop <- unlist(cfg$feature_blocks[cfg$ablation_sets[[setname]]])
      feats <- setdiff(base, drop)
      if (length(feats) == length(base)) {
        message("  ablation '", setname, "' drops no read-level features ",
                "(evaluated at sample x taxon level) -- skipping. [OI 11]")
        next
      }
      for (fold in splits$loeo) {
        m <- fold_masks(ft, fold)
        tr <- split_xy(ft[m$train], feats); te_dt <- ft[m$test]; te <- split_xy(te_dt, feats)
        sc <- tryCatch(fit_xgboost(tr, te, feats, fold$inner),
                       error = function(e) rep(NA_real_, nrow(te_dt)))
        k <- k + 1L
        preds[[k]] <- data.table(read_id = te_dt$read_id, donor = te_dt$donor,
          titration_level = te_dt$titration_level, species = te_dt$species,
          y = as.integer(te_dt$label == "positive"),
          end_reason_unblock = te_dt$end_reason_unblock,
          arm = paste0("combined_minus_", setname), model = "xgboost",
          scheme = "LOEO", fold = fold$fold, score = as.numeric(sc))
      }
    }
  }

  ## ---- H10 LOTO: combined arm, best tree family, unseen taxa ---------------
  if (do_loto && length(splits$loto)) {
    feats <- intersect(model_features("combined"), names(ft))
    for (fold in splits$loto) {
      m <- fold_masks(ft, fold)
      tr <- split_xy(ft[m$train], feats); te_dt <- ft[m$test]; te <- split_xy(te_dt, feats)
      sc <- tryCatch(fit_xgboost(tr, te, feats, fold$inner),
                     error = function(e) rep(NA_real_, nrow(te_dt)))
      k <- k + 1L
      preds[[k]] <- data.table(read_id = te_dt$read_id, donor = te_dt$donor,
        titration_level = te_dt$titration_level, species = te_dt$species,
        y = as.integer(te_dt$label == "positive"),
        end_reason_unblock = te_dt$end_reason_unblock,
        arm = "combined", model = "xgboost", scheme = "LOTO",
        fold = fold$fold, score = as.numeric(sc))
    }
  }

  ## ---- H9 decomposition: glmmTMB (1|donor) vs +(1|species), species-source swept ----
  ## Emits per-read scores under model names glmmTMB_<source> on the H9 arm so stage 06
  ## scores them like any model and stage 07 test_H9 pairs them:
  ##   glmmTMB_none        (1|donor) only          -> honest donor-transfer term
  ##   glmmTMB_truth       + (1|species) minimap2   -> species baseline (ground truth)
  ##   glmmTMB_classifier  + (1|species) top_genus  -> species baseline (classifier coarse rank)
  ## The difference (truth/classifier - none) is the species term's contribution, and the
  ## three-way source sweep exposes whether H9's conclusion is robust to it. [H9]
  if (do_h9 && requireNamespace("glmmTMB", quietly = TRUE)) {
    h9arm <- cfg$params$h9_arm %||% "combined"
    if (h9arm %in% arms) {
      feats <- intersect(model_features(h9arm), names(ft))
      srcs  <- cfg$params$h9_species_sources
      message(sprintf("-- H9 decomposition on '%s' (species sources: %s)",
                      h9arm, paste(names(srcs), collapse = ", ")))
      for (fold in splits$loeo) {
        m <- fold_masks(ft, fold)
        tr_dt <- ft[m$train]; te_dt <- ft[m$test]
        tr <- split_xy(tr_dt, feats); te <- split_xy(te_dt, feats)
        for (lab in names(srcs)) {
          col <- srcs[[lab]]
          if (!is.null(col) && !is.na(col) && nzchar(col)) {        # species source requested
            if (!col %in% names(ft)) next                            # column absent -> skip source
            sp_tr <- as.character(tr_dt[[col]]); sp_tr[is.na(sp_tr) | !nzchar(sp_tr)] <- "unknown"
            sp_te <- as.character(te_dt[[col]]); sp_te[is.na(sp_te) | !nzchar(sp_te)] <- "unknown"
            if (uniqueN(sp_tr) < 2L) next                            # single-level RE -> skip
          } else { sp_tr <- NULL; sp_te <- NULL }                    # 'none' -> (1|donor) only
          sc <- tryCatch(glmm_fit_predict(tr, te, feats, sp_tr, sp_te),
                         error = function(e) { message("  H9 glmmTMB[", lab, "] fold ",
                           fold$fold, " failed: ", conditionMessage(e)); rep(NA_real_, nrow(te_dt)) })
          k <- k + 1L
          preds[[k]] <- data.table(read_id = te_dt$read_id, donor = te_dt$donor,
            titration_level = te_dt$titration_level, species = te_dt$species,
            y = as.integer(te_dt$label == "positive"),
            end_reason_unblock = te_dt$end_reason_unblock,
            arm = h9arm, model = paste0("glmmTMB_", lab), scheme = "LOEO",
            fold = fold$fold, score = as.numeric(sc))
        }
      }
    }
  }

  out <- rbindlist(preds, fill = TRUE)
  pred_path <- file.path(cfg$paths$out_root, "predictions.tsv.gz")
  fwrite(out, pred_path, sep = "\t")
  message("Stage 05 complete: ", nrow(out), " predictions -> ", pred_path)
  invisible(out)
}

if (sys.nframe() == 0L) run_stage05()
