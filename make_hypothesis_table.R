#!/usr/bin/env Rscript
# =============================================================================
# make_hypothesis_table.R -- build the hypothesis RESULTS table (exact stats)
# from results/<tag>/hypothesis_tests.tsv. Writes hypothesis_results_table.tsv
# (project root) and prints a report-ready Markdown table to stdout.
# Env: FIG_RUN_TAG (default gt_fixed_id0.90_cov0.80), HYP_TABLE_OUT (output path).
# =============================================================================
suppressPackageStartupMessages(library(data.table))

tag <- Sys.getenv("FIG_RUN_TAG", "gt_fixed_id0.90_cov0.80")
out <- Sys.getenv("HYP_TABLE_OUT", "hypothesis_results_table.tsv")
ht  <- fread(file.path("results", tag, "hypothesis_tests.tsv"))
setkey(ht, id)

claim <- c(
  H1 = "Combined features beat BLAST-only",
  H2 = "ML beats fixed-threshold baseline",
  H3 = "ML gain larger at low titration",
  H4 = "Tree ensembles beat linear GLM",
  H5 = "Margin + human features drive read-level gain",
  H5b = "Genome breadth helps (sample x taxon)",
  H6 = "Sample x taxon aggregation beats read threshold",
  H7 = "Adapter+ reads enriched in false positives",
  H8 = "Neural net not better than XGBoost",
  H9 = "Donor random effect improves transfer",
  H9b_truth = "Species RE contribution (truth source)",
  H9b_classifier = "Species RE contribution (classifier source)",
  H10 = "Generalises to an unseen taxon (LOTO vs LOEO)",
  H11 = "Truncated/unblocked reads drive false positives",
  H12 = "Models better calibrated than fixed threshold")
test <- c(
  H1 = "Paired Wilcoxon (Holm-Sidak)", H2 = "Paired Wilcoxon (Holm-Sidak)",
  H3 = "LMM method x level (lm fallback)", H4 = "Paired Wilcoxon (Holm-Sidak)",
  H5 = "Paired Wilcoxon (Holm-Sidak)", H5b = "Paired Wilcoxon + bootstrap CI",
  H6 = "Paired Wilcoxon (Holm-Sidak)", H7 = "Fisher exact (deferred)",
  H8 = "-- (dropped)", H9 = "Paired Wilcoxon + bootstrap CI",
  H9b_truth = "Paired Wilcoxon + bootstrap CI",
  H9b_classifier = "Paired Wilcoxon + bootstrap CI",
  H10 = "Median degradation + bootstrap CI",
  H11 = "Score-distribution test", H12 = "Brier + paired Wilcoxon")

fmt  <- function(x) ifelse(is.na(x), "-",
                    ifelse(abs(x) >= 0.01, sprintf("%.3f", x), sprintf("%.2e", x)))
fmtp <- function(x) { x <- suppressWarnings(as.numeric(x))
  ifelse(is.na(x), "-", ifelse(x == 0, "<1e-16",
    ifelse(x < 1e-4, sprintf("%.1e", x), sprintf("%.4f", x)))) }

ord <- c("H1","H2","H3","H4","H5","H5b","H6","H7","H8",
         "H9","H9b_truth","H9b_classifier","H10","H11","H12")

res <- rbindlist(lapply(ord, function(i) {
  r  <- ht[id == i]
  md <- if (nrow(r)) r$median_diff else NA_real_
  lo <- if (nrow(r)) r$ci_lo else NA_real_
  hi <- if (nrow(r)) r$ci_hi else NA_real_
  ci_ok <- is.finite(lo) && is.finite(hi) && lo <= md && md <= hi   # drop H10's mislabeled CI
  data.table(
    hypothesis   = i,
    claim        = claim[[i]],
    family       = if (nrow(r)) r$family else NA_character_,
    test         = test[[i]],
    effect       = fmt(md),
    ci_95        = if (isTRUE(ci_ok)) sprintf("[%s, %s]", fmt(lo), fmt(hi)) else "-",
    p_raw        = fmtp(if (nrow(r)) r$wilcox_p else NA),
    p_holm_sidak = fmtp(if (nrow(r)) r$holm_sidak_p else NA),
    significant  = if (is.na(md)) "not run" else if (isTRUE(as.logical(r$significant))) "yes" else "no",
    note         = if (nrow(r)) r$note else "")
}))

fwrite(res, out, sep = "\t")

## report-ready Markdown (compact column set)
cat("\n<<<MARKDOWN>>>\n")
cat("| H | Claim | Test | Effect (median \u0394) | 95% CI | p | Holm\u2013\u0160\u00edd\u00e1k p | Significant |\n")
cat("|---|---|---|---|---|---|---|---|\n")
for (k in seq_len(nrow(res))) with(res[k],
  cat(sprintf("| %s | %s | %s | %s | %s | %s | %s | %s |\n",
              hypothesis, claim, test, effect, ci_95, p_raw, p_holm_sidak, significant)))
cat("<<<END>>>\n\nwrote ", out, " (", nrow(res), " rows)\n", sep = "")
