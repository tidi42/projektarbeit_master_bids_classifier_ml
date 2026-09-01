#!/usr/bin/env Rscript
# =============================================================================
# compare_gt_profiles.R -- head-to-head comparison of the two ground-truth runs
#
#   A = gt_fixed_id0.90_cov0.80       identity >= 0.90, coverage >= 0.80  (a priori)
#   B = gt_calculated_id0.92_cov0.50  identity >= 0.92, coverage >= 0.50  (data-driven,
#       = strictest pair still retaining >= gt_retain_frac of confident-Zymo reads)
#
# Everything else in the pipeline (reads, features, models, CV, statistics) is
# identical, so every difference below is attributable to the label definition.
#
# Writes tidy comparison tables to results/comparison_fixed_vs_calculated/ and
# prints a console digest. Consumed by create_gt_comparison_figures.R.
# =============================================================================

suppressPackageStartupMessages(library(data.table))

A_TAG <- Sys.getenv("GT_A_TAG", "gt_fixed_id0.90_cov0.80")
B_TAG <- Sys.getenv("GT_B_TAG", "gt_calculated_id0.92_cov0.50")
A_LAB <- "fixed (0.90/0.80)"
B_LAB <- "calculated (0.92/0.50)"
OUT   <- file.path("results", "comparison_fixed_vs_calculated")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

rd  <- function(tag, f) fread(file.path("results", tag, f))
say <- function(...) cat("\n== ", ..., " ==\n", sep = "")
wr  <- function(dt, name) { fwrite(dt, file.path(OUT, name), sep = "\t"); dt }

MODELS <- c("fixed_threshold", "glm", "glmmTMB", "ranger_rf", "xgboost")
ARMS3  <- c("blast_only", "kraken2_only", "combined")

# -----------------------------------------------------------------------------
# 1. Ground-truth profile provenance
# -----------------------------------------------------------------------------
gtA <- rd(A_TAG, "ground_truth_settings.tsv"); gtB <- rd(B_TAG, "ground_truth_settings.tsv")
prof <- rbind(cbind(profile_label = A_LAB, gtA), cbind(profile_label = B_LAB, gtB))
say("ground-truth profiles"); print(prof); wr(prof, "gt_profiles.tsv")

# the sweep that produced the 'calculated' pair (profile-independent, lives in work/)
sens <- if (file.exists("work/cutoff_sensitivity.tsv")) fread("work/cutoff_sensitivity.tsv") else NULL
if (!is.null(sens)) {
  sens[, chosen := gt_min_identity == gtB$gt_min_identity & gt_min_coverage == gtB$gt_min_coverage]
  sens[, apriori := gt_min_identity == gtA$gt_min_identity & gt_min_coverage == gtA$gt_min_coverage]
  wr(sens, "cutoff_sensitivity.tsv")
  say("cutoff sweep: the two profiles in context")
  print(sens[chosen | apriori])
}

# -----------------------------------------------------------------------------
# 2. Per-read label concordance (the actual ground-truth change)
#    zymo_ident / zymo_cov are carried in labels.tsv, so we can show exactly
#    WHERE in (identity, coverage) space the reclassified reads sit.
# -----------------------------------------------------------------------------
cols <- c("read_id", "titration_level", "label", "zymo_ident", "zymo_cov")
la <- rd(A_TAG, "labels.tsv")[, ..cols]; setnames(la, "label", "label_A")
lb <- rd(B_TAG, "labels.tsv")[, c("read_id", "label")];  setnames(lb, "label", "label_B")
lab <- merge(la, lb, by = "read_id")
rm(la, lb); invisible(gc())

xt <- lab[, .N, by = .(label_A, label_B)][order(-N)]
xt[, changed := label_A != label_B]
say("per-read label cross-tab (A = fixed rows, B = calculated cols)")
print(dcast(xt, label_A ~ label_B, value.var = "N", fill = 0L))
wr(xt, "label_crosstab.tsv")

n_tot <- nrow(lab); n_chg <- xt[changed == TRUE, sum(N)]
say(sprintf("reads relabelled: %s / %s (%.3f%%)", format(n_chg, big.mark = ","),
            format(n_tot, big.mark = ","), 100 * n_chg / n_tot))

# where do the flipped reads live in (identity, coverage) space?
flip <- lab[label_A != label_B]
if (nrow(flip)) {
  fl <- flip[, .(n = .N,
                 ident_min = round(min(zymo_ident, na.rm = TRUE), 4),
                 ident_med = round(median(zymo_ident, na.rm = TRUE), 4),
                 ident_max = round(max(zymo_ident, na.rm = TRUE), 4),
                 cov_min   = round(min(zymo_cov,  na.rm = TRUE), 4),
                 cov_med   = round(median(zymo_cov, na.rm = TRUE), 4),
                 cov_max   = round(max(zymo_cov,  na.rm = TRUE), 4)),
             by = .(label_A, label_B)]
  say("(identity, coverage) profile of the relabelled reads"); print(fl)
  wr(fl, "label_flip_profile.tsv")
  wr(flip[, .N, by = .(titration_level, label_A, label_B)][order(titration_level)],
     "label_flip_by_level.tsv")
}

# per-level composition + positive fraction, both profiles
comp <- rbind(
  lab[, .(profile = A_LAB, .N), by = .(titration_level, label = label_A)],
  lab[, .(profile = B_LAB, .N), by = .(titration_level, label = label_B)])
wr(comp, "label_composition.tsv")
posfrac <- comp[, .(pos_pct = round(100 * sum(N[label == "positive"]) / sum(N), 2)),
                by = .(profile, titration_level)]
say("positive-read fraction per level"); print(dcast(posfrac, titration_level ~ profile, value.var = "pos_pct"))
wr(posfrac, "positive_fraction.tsv")
rm(lab, flip); invisible(gc())

# -----------------------------------------------------------------------------
# 3. Read-level performance (LOEO, all reads, truncated included)
# -----------------------------------------------------------------------------
mrl <- function(tag) {
  m <- rd(tag, "metrics_read_level.tsv")
  m[scheme == "LOEO" & stratum == "all" & truncated_included == TRUE]
}
mA <- mrl(A_TAG); mB <- mrl(B_TAG)

perf_one <- function(m, lab) m[model %in% MODELS & arm %in% ARMS3,
  .(auprc = mean(auprc), prec = mean(prec_at_recall), fdr = mean(fdr_at_recall),
    auprc_sd = sd(auprc)), by = .(arm, model, recall_target)][, profile := lab][]
perf <- rbind(perf_one(mA, A_LAB), perf_one(mB, B_LAB))
wr(perf, "read_level_metrics.tsv")

pw <- dcast(perf[recall_target == 0.95], arm + model ~ profile, value.var = c("auprc", "prec"))
setnames(pw, gsub(" .*$", "", names(pw)))
say("AUPRC / precision@95% recall by arm x model, both profiles")
print(pw[order(arm, model)])
wr(pw, "read_level_wide.tsv")

# per-fold, combined arm -- fold-level paired view
fold <- rbind(mA[arm == "combined" & recall_target == 0.95 & model %in% MODELS,
                 .(fold, model, auprc, profile = A_LAB)],
              mB[arm == "combined" & recall_target == 0.95 & model %in% MODELS,
                 .(fold, model, auprc, profile = B_LAB)])
wr(fold, "per_fold_auprc.tsv")
say("per-fold AUPRC spread (combined arm)")
print(fold[, .(min = round(min(auprc), 4), med = round(median(auprc), 4),
               max = round(max(auprc), 4), sd = round(sd(auprc), 5)), by = .(profile, model)])

# feature-block ablation (XGBoost, FDR @ 99% recall)
ABL <- c("combined", "combined_minus_blast_margin", "combined_minus_human_competitor",
         "combined_minus_H5_key")
abl <- rbind(
  mA[model == "xgboost" & recall_target == 0.99 & arm %in% ABL,
     .(fdr_pct = 100 * mean(fdr_at_recall), sd = 100 * sd(fdr_at_recall)), by = arm][, profile := A_LAB][],
  mB[model == "xgboost" & recall_target == 0.99 & arm %in% ABL,
     .(fdr_pct = 100 * mean(fdr_at_recall), sd = 100 * sd(fdr_at_recall)), by = arm][, profile := B_LAB][])
abl[, fold_change := round(fdr_pct / fdr_pct[arm == "combined"], 3), by = profile]
say("feature ablation: FDR @ 99% recall (%)"); print(dcast(abl, arm ~ profile, value.var = "fdr_pct"))
wr(abl, "ablation.tsv")

# -----------------------------------------------------------------------------
# 4. Titration dependence
# -----------------------------------------------------------------------------
titr_one <- function(tag, lab) {
  mc <- rd(tag, "model_comparison.tsv")
  mc[facet == "concentration" & arm == "combined" & recall_target == 0.95 &
       key %chin% paste0("c", 1:5) & model %chin% MODELS,
     .(level = key, model, auprc, prec = prec_at_recall, profile = lab)]
}
titr <- rbind(titr_one(A_TAG, A_LAB), titr_one(B_TAG, B_LAB))
wr(titr, "titration.tsv")
say("precision @95% recall vs titration (combined arm)")
print(dcast(titr[model %chin% c("fixed_threshold", "xgboost")],
            model + level ~ profile, value.var = "prec"))

# the H3 quantity made explicit: model - baseline gap per level
gap <- dcast(titr, profile + level ~ model, value.var = "prec")
gap[, gap_xgb_minus_fixed := xgboost - fixed_threshold]
wr(gap, "titration_gap.tsv")
say("ML-over-baseline precision gap per level"); print(dcast(gap, level ~ profile, value.var = "gap_xgb_minus_fixed"))

# -----------------------------------------------------------------------------
# 5. Calibration
# -----------------------------------------------------------------------------
cal_one <- function(tag, lab) {
  c0 <- rd(tag, "calibration.tsv")
  c0 <- unique(c0[scheme == "LOEO" & arm == "combined" & model %chin% MODELS], by = c("model", "fold"))
  c0[, .(brier = mean(brier), sd = sd(brier)), by = model][, profile := lab][]
}
cal <- rbind(cal_one(A_TAG, A_LAB), cal_one(B_TAG, B_LAB))
say("Brier score (combined arm, lower = better)"); print(dcast(cal, model ~ profile, value.var = "brier"))
wr(cal, "calibration.tsv")

# -----------------------------------------------------------------------------
# 6. Hypothesis tests -- the headline concordance table
# -----------------------------------------------------------------------------
ht_one <- function(tag, lab) rd(tag, "hypothesis_tests.tsv")[
  , .(id, family, effect = median_diff, ci_lo, ci_hi, p = wilcox_p,
      p_adj = holm_sidak_p, sig = as.logical(significant), profile = lab)]
htA <- ht_one(A_TAG, A_LAB); htB <- ht_one(B_TAG, B_LAB)

ORD <- c("H1","H2","H3","H4","H5","H5b","H6","H7","H8",
         "H9","H9b_truth","H9b_classifier","H10","H11","H12")
hyp <- merge(htA[, .(id, family, effect_A = effect, p_A = p, padj_A = p_adj, sig_A = sig)],
             htB[, .(id, effect_B = effect, p_B = p, padj_B = p_adj, sig_B = sig)],
             by = "id", all = TRUE)
hyp <- hyp[match(ORD, id)]
hyp[, `:=`(
  delta        = effect_B - effect_A,
  ratio        = ifelse(is.finite(effect_A) & effect_A != 0, effect_B / effect_A, NA_real_),
  sign_agree   = ifelse(is.na(effect_A) | is.na(effect_B), NA,
                        sign(effect_A) == sign(effect_B)),
  verdict_agree = ifelse(is.na(effect_A) & is.na(effect_B), NA, sig_A == sig_B))]
say("hypothesis-level concordance")
print(hyp[, .(id, family, effect_A = signif(effect_A, 3), effect_B = signif(effect_B, 3),
              ratio = round(ratio, 3), sig_A, sig_B, sign_agree, verdict_agree)])
wr(hyp, "hypothesis_concordance.tsv")

tested <- hyp[!is.na(effect_A) & !is.na(effect_B)]
say(sprintf("tested hypotheses: %d | sign agreement: %d/%d | verdict agreement: %d/%d",
            nrow(tested), sum(tested$sign_agree), nrow(tested),
            sum(tested$verdict_agree), nrow(tested)))

# H9 random-effect decomposition
h9 <- rbind(cbind(profile = A_LAB, rd(A_TAG, "h9_species_sensitivity.tsv")),
            cbind(profile = B_LAB, rd(B_TAG, "h9_species_sensitivity.tsv")))
say("H9 / H9b random-effect decomposition")
print(h9[, .(profile, species_source, glmm = round(glmm_median_auprc, 4),
             glm = round(glm_median_auprc, 4), species_contribution = signif(species_contribution, 3))])
wr(h9, "h9_decomposition.tsv")

# mixed-model supplement (the honest check on the n = 8 primary tests)
mm <- rbind(cbind(profile = A_LAB, rd(A_TAG, "mixed_model_supplement.tsv")),
            cbind(profile = B_LAB, rd(B_TAG, "mixed_model_supplement.tsv")))
say("mixed-model supplement (160 donor x level units)")
print(mm[, .(profile, id, estimate = signif(mixed_estimate, 3), p = signif(mixed_p, 3), n_units)])
wr(mm, "mixed_model_supplement.tsv")

# -----------------------------------------------------------------------------
# 7. Operating-point confusion matrices (combined / XGBoost @ 95% recall)
# -----------------------------------------------------------------------------
conf_one <- function(tag, lab) {
  f <- file.path("results", tag, "read_calls_combined_xgboost.tsv.gz")
  if (!file.exists(f)) return(NULL)
  d <- fread(cmd = sprintf("zcat %s | cut -f7,9", shQuote(f)))   # ground_truth, model_call
  m <- d[, .N, by = .(ground_truth, model_call)]
  tp <- m[ground_truth == "positive" & model_call == "positive", sum(N)]
  fp <- m[ground_truth == "negative" & model_call == "positive", sum(N)]
  fn <- m[ground_truth == "positive" & model_call == "negative", sum(N)]
  tn <- m[ground_truth == "negative" & model_call == "negative", sum(N)]
  data.table(profile = lab, tp, fp, fn, tn, n = tp + fp + fn + tn,
             precision = tp / (tp + fp), recall = tp / (tp + fn),
             accuracy = (tp + tn) / (tp + fp + fn + tn))
}
conf <- rbindlist(list(conf_one(A_TAG, A_LAB), conf_one(B_TAG, B_LAB)))
if (nrow(conf)) {
  say("operating point: combined / XGBoost @ 95% recall"); print(conf)
  wr(conf, "confusion.tsv")
}

cat("\n\nAll comparison tables written to: ", normalizePath(OUT), "\n", sep = "")
