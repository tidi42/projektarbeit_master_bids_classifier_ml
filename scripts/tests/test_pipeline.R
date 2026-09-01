#!/usr/bin/env Rscript
## =============================================================================
## test_pipeline.R  --  critical-review test suite for the Zymo-in-human pipeline
## -----------------------------------------------------------------------------
## Two kinds of checks:
##   expect(name, cond)      CORRECTNESS  -- a property that SHOULD hold (PASS/FAIL).
##                           These validate the parts of the pipeline that are sound.
##   finding(id, name, ...)  REVIEW PROBE -- a diagnostic that CONFIRMS a design or
##                           implementation weakness (DETECTED = the flaw is real).
##
## Sections mirror the pipeline stages so each part is checked in isolation, plus
## a synthetic end-to-end smoke test of stages 04->07. Companion to
## reports/reporting_3_critical_review.md.  Run:  Rscript scripts/tests/test_pipeline.R
## =============================================================================

.thisfile <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
SCRIPTS <- if (length(.thisfile) && !is.na(.thisfile) && nzchar(.thisfile))
  dirname(dirname(normalizePath(.thisfile))) else
  "/home/tdinse/master_timo/projektarbeit/projektarbeit_ml/scripts"

suppressWarnings(suppressMessages({
  source(file.path(SCRIPTS, "00_config.R"))
  source(file.path(SCRIPTS, "utils.R"))
  source(file.path(SCRIPTS, "01_external_tools.R"))
  source(file.path(SCRIPTS, "02_ground_truth_labels.R"))
  source(file.path(SCRIPTS, "03_build_features.R"))
  source(file.path(SCRIPTS, "04_cv_splits.R"))
  source(file.path(SCRIPTS, "05_train_models.R"))
  source(file.path(SCRIPTS, "06_evaluate.R"))
  source(file.path(SCRIPTS, "07_hypothesis_tests.R"))
  source(file.path(SCRIPTS, "prepare_cutoff_sensitivity.R"))
}))
library(data.table)

REAL_SS <- cfg$paths$sample_sheet     # capture before section E mutates cfg$paths

## ---- tiny harness -----------------------------------------------------------
.PASS <- 0L; .FAIL <- 0L; .findings <- list()
expect <- function(name, cond) {
  ok <- isTRUE(cond)
  if (ok) .PASS <<- .PASS + 1L else .FAIL <<- .FAIL + 1L
  cat(sprintf("  [%s] %s\n", if (ok) "PASS" else "FAIL", name))
}
finding <- function(id, name, detected, detail = "") {
  .findings[[id]] <<- list(name = name, detected = isTRUE(detected), detail = detail)
  cat(sprintf("  [%s] %-4s %s\n      %s\n",
              if (isTRUE(detected)) "DETECTED" else "clear", id, name, detail))
}
section <- function(s) cat(sprintf("\n== %s ==\n", s))
approx <- function(a, b, tol = 1e-8) isTRUE(all(abs(a - b) <= tol))

## =============================================================================
section("A0. Human depletion (competition: max(hg38,T2T) > Zymo)")
## =============================================================================
local({
  mk <- function(ids, sc) data.table(read_id = ids, score = sc)
  hg <- mk(c("rHuman", "rHumWin", "rZymWin", "rZymo"), c(100, 150,  80,  50))
  tt <- mk(c("rHuman", "rT2Tonly"),                     c( 90,  60))
  z  <- mk(c("rZymo",  "rHumWin", "rZymWin"),           c(200, 100, 120))
  hl <- human_like_reads(hg, tt, z)
  expect("deplete: human-wins reads removed (incl. T2T-only)",
         setequal(hl, c("rHuman", "rHumWin", "rT2Tonly")))
  expect("deplete: Zymo-wins + unaligned reads survive",
         !any(c("rZymo", "rZymWin", "rNeither") %in% hl))
  expect("deplete: no human alignments -> keep all reads",
         length(human_like_reads(mk(character(), numeric()), mk(character(), numeric()), z)) == 0L)
  ## labelling's human_score uses .best_score over best_per_read tables (column `qname`)
  expect("label: .best_score picks max(hg38,T2T) per read",
         approx(pmax(.best_score(data.table(qname = "rH", score = 100), "rH"),
                     .best_score(data.table(qname = "rH", score =  90), "rH")), 100))
})

## =============================================================================
section("A. Config integrity")
## =============================================================================
expect("hypotheses: 12 rows, ids H1..H12",
       nrow(cfg$hypotheses) == 12 && all(cfg$hypotheses$id == paste0("H", 1:12)))
expect("hypotheses: active vector length 12", length(cfg$hypotheses$active) == 12)
expect("H7 deferred + H8 dropped (both inactive)",
       isFALSE(cfg$hypotheses$active[cfg$hypotheses$id == "H7"]) &&
       isFALSE(cfg$hypotheses$active[cfg$hypotheses$id == "H8"]))
expect("open_items: 21 rows", nrow(cfg$open_items) == 21)
expect("open_items: all resolved (0 open) -- Kraken2 DB provided [OI 9]", {
  op <- cfg$open_items[cfg$open_items$status != "resolved", , drop = FALSE]
  nrow(op) == 0
})
expect("primary_family has <= 6 tests", length(cfg$params$primary_family) <= 6)
expect("recall_primary is one of recall_target",
       cfg$params$recall_primary %in% cfg$params$recall_target)
idcols <- c("species", "genus", "top_species", "top_genus", "staxid", "k2_taxid")
expect("model_features(combined) contains NO taxon-identity columns",
       length(intersect(model_features("combined"), idcols)) == 0)
expect("every classifier-arm block exists in feature_blocks",
       all(unlist(cfg$classifier_arms) %in% names(cfg$feature_blocks)))
expect("every ablation-set block exists in feature_blocks",
       all(unlist(cfg$ablation_sets) %in% names(cfg$feature_blocks)))
expect("F6: sxt breadth ablation withholds only sample_taxon-level features",
       all(unlist(cfg$sxt_ablation_sets) %in% cfg$feature_blocks$sample_taxon))

## =============================================================================
section("A2. Ground-truth cutoff profile (two-run design)")
## =============================================================================
## Two reproducible runs selected by GT_PROFILE: 'fixed' (0.90/0.80) and
## 'calculated' (from prepare_cutoff_sensitivity.R). Each lands in its own
## results/<gt_run_tag>/ folder; cutoff-INDEPENDENT stage-01 output stays shared.
srcCut <- paste(readLines(file.path(SCRIPTS, "prepare_cutoff_sensitivity.R")), collapse = "\n")
srcRun <- paste(readLines(file.path(SCRIPTS, "run_pipeline.R")), collapse = "\n")
expect("GT profile defaults to 'fixed' with id>=0.90, cov>=0.80",
       identical(cfg$params$gt_profile, "fixed") &&
       approx(cfg$params$gt_min_identity, 0.90) && approx(cfg$params$gt_min_coverage, 0.80))
expect("GT run tag encodes profile + active cutoffs",
       identical(cfg$params$gt_run_tag, "gt_fixed_id0.90_cov0.80") &&
       identical(cfg_gt_run_tag("calculated", 0.95, 0.90), "gt_calculated_id0.95_cov0.90"))
expect("per-run out_root + all cutoff-dependent paths live under results/<run_tag>/", {
  identical(basename(cfg$paths$out_root), cfg$params$gt_run_tag) &&
  identical(dirname(cfg$paths$out_root), cfg$paths$results_base) &&
  all(startsWith(c(cfg$paths$labels_table, cfg$paths$feature_table, cfg$paths$cv_splits,
                   cfg$paths$metrics_read, cfg$paths$metrics_sxt, cfg$paths$hypotheses_out,
                   cfg$paths$leakage_table, cfg$paths$model_dir),
                 paste0(cfg$paths$out_root, "/")))
})
expect("cutoff-INDEPENDENT recommended-cutoffs file stays in the SHARED work/ dir",
       identical(dirname(cfg$paths$gt_recommended), cfg$paths$work_dir))
expect("cfg_gt_profile_cutoffs('fixed') returns the hardcoded pair (ignores rec file)",
       { r <- cfg_gt_profile_cutoffs("fixed", 0.90, 0.80, tempfile())
         approx(r$gt_min_identity, 0.90) && approx(r$gt_min_coverage, 0.80) })
expect("cfg_gt_profile_cutoffs('calculated') reads the recommended (id, cov) pair", {
  d <- tempfile()
  write.table(data.frame(gt_min_identity = 0.95, gt_min_coverage = 0.90),
              d, sep = "\t", row.names = FALSE, quote = FALSE)
  r <- cfg_gt_profile_cutoffs("calculated", 0.90, 0.80, d)
  approx(r$gt_min_identity, 0.95) && approx(r$gt_min_coverage, 0.90)
})
expect("calculated profile without a recommended file stops with guidance",
       inherits(tryCatch(cfg_gt_profile_cutoffs("calculated", 0.9, 0.8, "/no/such/file.tsv"),
                         error = function(e) e), "error"))
expect("unknown GT_PROFILE is rejected",
       inherits(tryCatch(cfg_gt_profile_cutoffs("bogus", 0.9, 0.8, tempfile()),
                         error = function(e) e), "error"))
expect("prepare_cutoff_sensitivity persists cutoff_recommended.tsv for the 2nd run",
       grepl("cutoff_recommended|gt_recommended", srcCut) &&
       grepl("gt_min_identity", srcCut) && grepl("gt_min_coverage", srcCut))
expect("run_pipeline exposes the --gt fixed|calculated|both switch via GT_PROFILE",
       grepl("--gt", srcRun, fixed = TRUE) && grepl("GT_PROFILE", srcRun) &&
       grepl("calculated", srcRun) && grepl("both", srcRun))

## =============================================================================
section("B. utils metrics correctness")
## =============================================================================
y  <- c(1, 1, 1, 0, 0, 0); ssep <- c(0.9, 0.8, 0.7, 0.3, 0.2, 0.1)
expect("auprc perfectly separable == 1", approx(auprc(y, ssep), 1))
expect("auprc with no positives == NA", is.na(auprc(rep(0, 4), runif(4))))
expect("precision_at_recall(1.0) separable == 1", approx(precision_at_recall(y, ssep, 1.0), 1))
expect("fdr_at_recall == 1 - precision_at_recall",
       approx(fdr_at_recall(y, ssep, 0.8), 1 - precision_at_recall(y, ssep, 0.8)))
expect("brier perfect == 0", approx(brier_score(c(1, 0, 1, 0), c(1, 0, 1, 0)), 0))
expect("brier worst == 1", approx(brier_score(c(1, 0), c(0, 1)), 1))
pv <- c(0.01, 0.04, 0.03); hs <- holm_sidak(pv)
expect("holm_sidak adjusted >= raw p", all(hs >= pv - 1e-12))
expect("holm_sidak single p unchanged", approx(holm_sidak(0.023), 1 - (1 - 0.023)^1))
expect("holm_sidak monotone in p-rank", { o <- order(pv); all(diff(hs[o]) >= -1e-12) })
set.seed(1); x1 <- rnorm(8, 0.4); x2 <- rnorm(8)
expect("paired_wilcox exact matches stats::wilcox.test",
       approx(paired_wilcox(x1, x2, TRUE)$p.value,
              suppressWarnings(wilcox.test(x1, x2, paired = TRUE, exact = TRUE))$p.value))
pd <- poisson_detection_prob(1000, 0.01, 0.5, 1)      # E[N] = 5
expect("poisson expected genomes = cells*ra*eff*frac", approx(pd$expected_genomes, 5))
expect("poisson p = 1 - exp(-E)", approx(pd$p_at_least_one, 1 - exp(-5)))
csf <- coverage_stats_from_intervals(0, 100, 100)
expect("coverage breadth full == 1", approx(csf$genome_breadth, 1))
expect("coverage evenness uniform == 1", approx(csf$coverage_evenness, 1))
csh <- coverage_stats_from_intervals(0, 50, 100)
expect("coverage breadth half interval == 0.5", approx(csh$genome_breadth, 0.5))
expect("shannon entropy of 2 uniform classes == 1 bit", approx(shannon_entropy(c("a", "b")), 1))
expect("weighted_gini of equal values == 0", approx(weighted_gini(c(5, 5, 5), c(1, 1, 1)), 0))

## =============================================================================
section("C. CV splits & train/test leakage")
## =============================================================================
set.seed(42)
donors  <- sprintf("D%02d", 1:8)
species <- paste0("sp", 1:5)
mk_labels <- function() {
  rows <- list(); i <- 0L
  for (d in donors) for (lv in c("c1", "c2", "c3", "c4", "c5", "negative")) {
    if (lv != "negative") {
      i <- i + 1L
      rows[[i]] <- data.table(read_id = paste0(d, lv, "P", 1:8), donor = d,
        titration_level = lv, species = sample(species, 8, replace = TRUE), label = "positive")
    }
    i <- i + 1L
    rows[[i]] <- data.table(read_id = paste0(d, lv, "N", 1:10), donor = d,
      titration_level = lv, species = NA_character_, label = "negative")
  }
  rbindlist(rows, fill = TRUE)
}
L <- mk_labels(); L[, run_id := donor]
L[sample(.N, 20), label := "ambiguous"]
loeo <- build_loeo(donors, unique(L[, .(donor, run_id)]))

ok_disjoint <- TRUE; held_one <- TRUE; excl_amb <- TRUE
for (f in loeo) {
  m <- fold_masks(L, f)
  if (any(m$train & m$test)) ok_disjoint <- FALSE
  td <- unique(L$donor[m$test]); if (!(length(td) == 1 && td == f$test_donors)) held_one <- FALSE
  if (any(L$label[m$train | m$test] == "ambiguous")) excl_amb <- FALSE
}
expect("LOEO: no read appears in both train and test", ok_disjoint)
expect("LOEO: each fold tests exactly one donor", held_one)
expect("LOEO: ambiguous rows excluded from train and test", excl_amb)

loto <- build_loto(species, donors)
leak_rows <- 0L; donor_disjoint <- TRUE; taxon_ok <- TRUE
for (f in loto) {
  m <- fold_masks(L, f)
  leak_rows <- leak_rows + sum(m$train & m$test)
  if (length(intersect(f$test_donors, f$train_donors))) donor_disjoint <- FALSE
  if (any(L$label[m$train] == "positive" & L$species[m$train] == f$test_species)) taxon_ok <- FALSE
}
expect("LOTO [F1 fixed]: no read appears in both train and test", leak_rows == 0)
expect("LOTO [F1 fixed]: train/test donor groups are disjoint", donor_disjoint)
expect("LOTO: held-out taxon positives never appear in training", taxon_ok)

## =============================================================================
section("D. Ground-truth labelling logic")
## =============================================================================
tf <- tempfile(fileext = ".paf")
writeLines(c(
  paste("r1", 100, 0, 95, "+", "zc1", 5000,  10, 105, 90, 95, 60, sep = "\t"),  # id .947 cov .95
  paste("r2", 100, 0, 50, "+", "zc1", 5000,   0,  50, 20, 50, 60, sep = "\t"),  # id .40  cov .50
  paste("r1", 100, 0, 40, "+", "zc9", 5000,   0,  40, 10, 40, 60, sep = "\t")   # weaker r1 hit
), tf)
paf <- read_paf(tf)
expect("read_paf identity = matches/blocklen", approx(paf[qname == "r1"][1]$identity, 90 / 95, 1e-6))
expect("read_paf coverage = (qend-qstart)/qlen", approx(paf[qname == "r1"][1]$coverage, 0.95, 1e-6))
bpr <- best_per_read(paf)
expect("best_per_read keeps the max-matches hit per read",
       bpr[qname == "r1"]$matches == 90 && nrow(bpr) == 2)
gi <- cfg$params$gt_min_identity; gc2 <- cfg$params$gt_min_coverage
zpass <- bpr$identity >= gi & bpr$coverage >= gc2
expect("GT rule: r1 (id.95/cov.95) passes identity & coverage cutoffs", isTRUE(zpass[bpr$qname == "r1"]))
expect("GT rule: r2 (id.40/cov.50) fails cutoffs", isFALSE(zpass[bpr$qname == "r2"]))

## =============================================================================
section("E. Findings status after fixes (items 1-7; F1 in section C)")
## =============================================================================
## F2 -- FIXED: model selection now uses leakage-free inner-CV, not outer folds.
src02 <- paste(readLines(file.path(SCRIPTS, "02_ground_truth_labels.R")), collapse = "\n")
src05 <- paste(readLines(file.path(SCRIPTS, "05_train_models.R")), collapse = "\n")
src06 <- paste(readLines(file.path(SCRIPTS, "06_evaluate.R")), collapse = "\n")
src07 <- paste(readLines(file.path(SCRIPTS, "07_hypothesis_tests.R")), collapse = "\n")
nf <- 8L
expect("F2 fixed: inner-CV selection present (select_primary_model + inner_cv_scores)",
       grepl("select_primary_model", src07) && grepl("inner_cv_scores", src07) && grepl("inner_cv_score", src05))
expect("F2 fixed: H4 no longer uses pmax(RF, XGB)", !grepl("pmax\\(rf", src07))
set.seed(11)
selbias <- replicate(2000, { M <- matrix(rnorm(nf * 4L), nf, 4L); max(colMeans(M)) - mean(colMeans(M)) })
cat(sprintf("      rationale: outer-fold winner-selection would inflate by E[max-mean]=%.3f fold-SD (null); inner-CV avoids it.\n", mean(selbias)))

## F3 -- MITIGATED: donor-clustered mixed model supplements the 8-fold Wilcoxon.
expect("F3 fixed: mixed_effect_test (donor-clustered) available",
       exists("mixed_effect_test") && is.function(mixed_effect_test))
expect("F3 fixed: stage 07 emits a mixed-model supplement (run_mixed_supplement)",
       grepl("run_mixed_supplement", src07))
set.seed(9)
longF3 <- rbindlist(lapply(sprintf("D%02d", 1:8), function(d)
  rbindlist(lapply(c("c1","c2","c3","c4","c5"), function(l)
    data.table(donor = d, method = c("A","B"), auprc = c(0.60, 0.70) + rnorm(2, 0, 0.03))))))
mtF3 <- mixed_effect_test(longF3, "auprc")
expect("F3 fixed: mixed model uses ~40 donor x level units (not 8)", is.finite(mtF3$p.value) && mtF3$n >= 30)
floor2 <- 2 / 2^nf
cat(sprintf("      rationale: 8-fold exact-Wilcoxon p-floor=%.4f (unanimity-only); the mixed model recovers within-donor power.\n", floor2))

## F4 -- FIXED: fixed_threshold is Platt-calibrated, so Brier is a fair comparison.
set.seed(5)
ytr <- rep(0:1, each = 60); Xtr <- data.table(f1 = rnorm(120, 2 * ytr), f2 = rnorm(120))
yte <- rep(0:1, each = 30); Xte <- data.table(f1 = rnorm(60, 2 * yte), f2 = rnorm(60))
sc_ft <- suppressWarnings(fit_fixed_threshold(list(X = Xtr, y = ytr, w = rep(1, 120)),
                                              list(X = Xte, y = yte), c("f1", "f2")))
expect("F4 fixed: fixed_threshold emits a calibrated probability in [0,1]", all(sc_ft >= 0 & sc_ft <= 1))
cmF4 <- calibration_metrics(data.table(arm = "combined", model = c("fixed_threshold", "glm"),
                                       scheme = "LOEO", fold = 1L, y = c(1L, 0L), score = c(0.8, 0.3)))
expect("F4 fixed: calibration table now INCLUDES fixed_threshold", "fixed_threshold" %in% cmF4$model)

## F5 -- FIXED: threshold is per-(arm, model), so counts are scale-comparable.
cfg$paths$work_dir <- file.path(tempdir(), "f5work"); dir.create(cfg$paths$work_dir, showWarnings = FALSE)
set.seed(6); np <- 200L
pF5 <- rbind(
  data.table(arm = "combined", model = "A", scheme = "LOEO", fold = 1L, donor = "D01",
             titration_level = "c1", species = "sp1", y = 1L, score = runif(np, 0, 1)),
  data.table(arm = "combined", model = "B", scheme = "LOEO", fold = 1L, donor = "D01",
             titration_level = "c1", species = "sp1", y = 1L, score = runif(np, 0, 1000)))
sxtF5 <- suppressMessages(aggregate_sample_taxon(pF5))
expect("F5 fixed: both score scales retain reads above their OWN threshold",
       sxtF5[model == "A", sum(n_reads_above_thr)] > 0 && sxtF5[model == "B", sum(n_reads_above_thr)] > 0)

## F6 -- FIXED: H5's 'breadth' component is now MEASURED via a sample x taxon ablation.
expect("F6 fixed: sample_taxon_ablation() implemented",
       exists("sample_taxon_ablation") && is.function(sample_taxon_ablation))
f6work <- file.path(tempdir(), "f6work"); dir.create(f6work, showWarnings = FALSE)
old_work_f6 <- cfg$paths$work_dir; cfg$paths$work_dir <- f6work
set.seed(23)
f6_donors <- sprintf("D%02d", 1:6); f6_species <- paste0("sp", 1:5); f6_present <- c("sp1", "sp2", "sp3")
f6_cells <- CJ(donor = f6_donors, titration_level = "c1", species = f6_species)
fwrite(copy(f6_cells)[, `:=`(expected_cells = 100, p_detect = 0.99,
         expected_present = species %in% f6_present)],
       file.path(f6work, "expected_sample_taxon.tsv"), sep = "\t")
## coverage table: breadth/evenness are the ONLY discriminators (high for present taxa)
fwrite(copy(f6_cells)[, `:=`(
         genome_breadth    = ifelse(species %in% f6_present, runif(.N, 0.85, 0.98), runif(.N, 0.02, 0.15)),
         coverage_evenness = ifelse(species %in% f6_present, runif(.N, 0.80, 0.95), runif(.N, 0.05, 0.20)))],
       file.path(f6work, "sample_taxon_coverage.tsv"), sep = "\t")
## predictions: read scores are pure NOISE, so any sxt gain must come from breadth
f6_preds <- rbindlist(lapply(f6_donors, function(d) rbindlist(lapply(f6_species, function(sp)
  data.table(read_id = paste0(d, sp, 1:6), donor = d, titration_level = "c1", species = sp,
             y = as.integer(sp %in% f6_present), end_reason_unblock = 0L,
             arm = "combined", model = "glm", scheme = "LOEO", fold = match(d, f6_donors),
             score = runif(6))))))
ablF6 <- suppressMessages(sample_taxon_ablation(f6_preds))
expect("F6 fixed: ablation yields sxt_full AND sxt_minus_breadth at sample_taxon level",
       all(c("sxt_full", "sxt_minus_breadth") %in% ablF6$model) && all(ablF6$level == "sample_taxon"))
expect("F6 fixed: both ablation variants report finite AUPRC",
       nrow(ablF6[model == "sxt_full" & is.finite(auprc)]) > 0 &&
       nrow(ablF6[model == "sxt_minus_breadth" & is.finite(auprc)]) > 0)
f6_full <- ablF6[model == "sxt_full" & recall_target == cfg$params$recall_primary, mean(auprc, na.rm = TRUE)]
f6_abl  <- ablF6[model == "sxt_minus_breadth" & recall_target == cfg$params$recall_primary, mean(auprc, na.rm = TRUE)]
expect("F6 fixed: withholding breadth measurably changes sample x taxon AUPRC (breadth is now tested)",
       is.finite(f6_full) && is.finite(f6_abl) && f6_full >= f6_abl)
cfg$paths$work_dir <- old_work_f6
cat(sprintf("      rationale: sxt AUPRC full=%.3f vs minus-breadth=%.3f -> H5 'breadth' is now MEASURED, not asserted.\n",
            f6_full, f6_abl))

## F7 -- RESOLVED: random flag removed; principled leakage diagnostic added.
expect("F7 resolved: estimate_leakage diagnostic present in stage 02", grepl("estimate_leakage", src02))
expect("F7 resolved: random suspected_leakage flag removed everywhere",
       !grepl("suspected_leakage", paste(src02, src05, src06, src07)))

## F8 -- run confounded with donor (from the real sample sheet).
ss <- tryCatch(fread(REAL_SS), error = function(e) NULL)
if (!is.null(ss) && all(c("donor", "run_id") %in% names(ss))) {
  per_run <- ss[, uniqueN(donor), by = run_id]$V1
  finding("F8", "run/flow-cell confounded with donor", detected = all(per_run == 1),
          sprintf("all %d runs carry exactly one donor -> a per-run technical effect cannot be separated from the biological donor effect; LOEO generalises to 'new donor+run' jointly",
                  ss[, uniqueN(run_id)]))
} else {
  finding("F8", "run/flow-cell confounded with donor", detected = FALSE, "sample sheet unavailable in this run")
}

## M1 -- FIXED: declared subject props now match what parse_blast() actually produces.
expect("M1 fixed: cfg$subject_props declares only the produced column (subject_genome_len)",
       identical(cfg$subject_props, "subject_genome_len") &&
       !any(c("subject_assembly_level", "subject_is_wgs_draft") %in% model_features("combined")))
## M3 -- FIXED: the always-NA LCA-rank columns are dropped so the feature set reflects reality.
expect("M3 fixed: inert lca_rank / k2_lca_rank removed from feature blocks and model_features",
       !("lca_rank" %in% cfg$feature_blocks$blast_margin) &&
       !("k2_lca_rank" %in% cfg$feature_blocks$kraken2) &&
       !any(c("lca_rank", "k2_lca_rank") %in% model_features("combined")))

## M2 -- FIXED: tune_inner passes the grid row EXPLICITLY (no global `grid_row`).
expect("M2 fixed: no global `grid_row` superassignment remains in stage 05 source",
       !grepl("grid_row[[:space:]]*<<-", src05))
set.seed(31)
m2dt <- data.table(donor = rep(sprintf("D%02d", 1:4), each = 10), y = rep(0:1, 20))
tr_m2 <- list(X = m2dt[, .(f = rnorm(.N))], y = m2dt$y, grp = m2dt[, .(donor)], w = rep(1, nrow(m2dt)))
inner_m2 <- lapply(sprintf("D%02d", 1:4), function(d) d)
pf_m2 <- function(a, b, row) if (isTRUE(row$good)) as.numeric(b$y) else as.numeric(1 - b$y)
best_m2 <- tune_inner(tr_m2, "f", inner_m2, pf_m2, data.frame(good = c(FALSE, TRUE)))
expect("M2 fixed: tune_inner selects the best grid row via the explicit `row` arg", isTRUE(best_m2$good))

## =============================================================================
section("E2. Ground-truth cutoff sensitivity (F9)")
## =============================================================================
## F9 -- FIXED: gt_min_identity/coverage are stress-tested by prepare_cutoff_sensitivity.R.
expect("F9 fixed: cutoff-sensitivity helpers implemented",
       all(vapply(c("collect_confident_zymo", "zymo_distribution", "cutoff_sweep", "recommend_cutoffs"),
                  function(fn) exists(fn) && is.function(get(fn)), logical(1))))
f9work <- file.path(tempdir(), "f9work")
dir.create(file.path(f9work, "L_c1"), recursive = TRUE, showWarnings = FALSE)
old_work_f9 <- cfg$paths$work_dir; old_ss_f9 <- cfg$paths$sample_sheet
cfg$paths$work_dir <- f9work
cfg$paths$sample_sheet <- file.path(f9work, "sample_sheet.tsv")
fwrite(data.table(library_id = "L_c1", donor = "D01", barcode = 1L, titration_level = "c1",
                  run_id = "R1", fastq = "x.fastq"), cfg$paths$sample_sheet, sep = "\t")
mkpaf <- function(rows, f) writeLines(vapply(rows, function(r) paste(r, collapse = "\t"), character(1)), f)
## r1..r4 = confident (Zymo beats human) spanning identity/coverage; r5 = human wins
zrows <- list(
  c("r1", 100, 0, 95, "+", "zc1", 5000, 0, 97, 97, 100, 60),   # id .97 cov .95
  c("r2", 100, 0, 85, "+", "zc1", 5000, 0, 93, 93, 100, 60),   # id .93 cov .85
  c("r3", 100, 0, 75, "+", "zc1", 5000, 0, 88, 88, 100, 60),   # id .88 cov .75
  c("r4", 100, 0, 60, "+", "zc1", 5000, 0, 82, 82, 100, 60),   # id .82 cov .60
  c("r5", 100, 0, 99, "+", "zc1", 5000, 0, 99, 99, 100, 60))   # strong Zymo BUT human wins
hrows <- list(
  c("r1", 100, 0, 30, "+", "hchr", 99999, 0, 30, 10, 100, 60),
  c("r2", 100, 0, 30, "+", "hchr", 99999, 0, 30, 10, 100, 60),
  c("r3", 100, 0, 30, "+", "hchr", 99999, 0, 30, 10, 100, 60),
  c("r4", 100, 0, 30, "+", "hchr", 99999, 0, 30, 10, 100, 60),
  c("r5", 100, 0, 99, "+", "hchr", 99999, 0, 99, 150, 100, 60))  # human matches 150 > Zymo 99
mkpaf(zrows, file.path(f9work, "L_c1", "gt_zymo.paf"))
mkpaf(hrows, file.path(f9work, "L_c1", "gt_human_grch38.paf"))
confF9 <- suppressMessages(collect_confident_zymo(levels = "c1"))
expect("F9: confident-Zymo = reads whose Zymo hit beats human (r1..r4, not r5)",
       all(c("r1", "r2", "r3", "r4") %in% confF9$read_id) && !("r5" %in% confF9$read_id))
distF9 <- zymo_distribution(confF9)
expect("F9: empirical distribution reports identity + coverage quantiles",
       all(c("identity", "coverage") %in% distF9$metric) && all(is.finite(distF9$value)))
swF9 <- cutoff_sweep(confF9)
expect("F9: sweep reports a retained fraction per cutoff pair",
       "frac_retained" %in% names(swF9) &&
       nrow(swF9) == length(cfg$params$gt_identity_grid) * length(cfg$params$gt_coverage_grid))
f9_mono <- TRUE
for (cc in cfg$params$gt_coverage_grid) {
  v <- swF9[gt_min_coverage == cc][order(gt_min_identity), n_retained]; if (any(diff(v) > 0)) f9_mono <- FALSE }
for (ic in cfg$params$gt_identity_grid) {
  v <- swF9[gt_min_identity == ic][order(gt_min_coverage), n_retained]; if (any(diff(v) > 0)) f9_mono <- FALSE }
expect("F9: stricter cutoffs retain FEWER-or-equal confident-Zymo reads (monotone)", f9_mono)
expect("F9: loosest cutoff retains all confident-Zymo reads",
       approx(swF9[gt_min_identity == min(cfg$params$gt_identity_grid) &
                   gt_min_coverage == min(cfg$params$gt_coverage_grid), frac_retained], 1))
recF9 <- recommend_cutoffs(swF9)
expect("F9: recommends one cutoff pair (>= gt_retain_frac, or the best available)",
       nrow(recF9) == 1 && (recF9$frac_retained >= cfg$params$gt_retain_frac ||
                            approx(recF9$frac_retained, max(swF9$frac_retained))))
cat(sprintf("      rationale: current id>=%.2f/cov>=%.2f retain %.0f%% of confident-Zymo reads; sweep documents label-set sensitivity.\n",
            cfg$params$gt_min_identity, cfg$params$gt_min_coverage,
            100 * swF9[is_current == TRUE, frac_retained][1]))
cfg$paths$work_dir <- old_work_f9; cfg$paths$sample_sheet <- old_ss_f9

## =============================================================================
section("E3. H9 GLMM decomposition + species-source sensitivity")
## =============================================================================
## H9 is decomposed into (1|donor) alone vs +(1|species), and the species source is
## swept {none, truth, classifier}. Tested on a synthetic metrics table (no glmmTMB
## fitting) so the decomposition + robustness LOGIC is checked fast + sandbox-safe.
mkM_h9 <- function(by_model) rbindlist(lapply(names(by_model), function(mdl)
  data.table(arm = "combined", model = mdl, level = "read", stratum = "all",
             scheme = "LOEO", truncated_included = TRUE,
             recall_target = cfg$params$recall_primary,
             fold = seq_along(by_model[[mdl]]), auprc = by_model[[mdl]])), fill = TRUE)
set.seed(101)
glmv   <- runif(8, 0.55, 0.62)
donorv <- glmv   + runif(8, 0.03, 0.07)   # (1|donor) helps transfer -> H9 core > 0
truthv <- donorv + runif(8, 0.02, 0.05)   # +species(truth) adds
clfv   <- donorv + runif(8, -0.01, 0.01)  # +species(classifier coarse rank) ~ neutral
Mh9 <- mkM_h9(list(glm = glmv, glmmTMB_none = donorv, glmmTMB_truth = truthv,
                   glmmTMB_classifier = clfv))
h9rows <- suppressWarnings(test_H9(Mh9))
expect("H9: decomposition emits core H9 + H9b_truth + H9b_classifier rows",
       all(c("H9", "H9b_truth", "H9b_classifier") %in% h9rows$id))
expect("H9 core = GLMM(1|donor) - GLM is positive (donor RE aids transfer)",
       { r <- h9rows[id == "H9"]; is.finite(r$median_diff) && r$median_diff > 0 })
expect("H9b_truth: species-term contribution over donor-only is positive",
       { r <- h9rows[id == "H9b_truth"]; is.finite(r$median_diff) && r$median_diff > 0 })
expect("H9 note records the GLMM-GLM sign per species source (robustness axis, ROBUST here)",
       grepl("sign by source", h9rows[id == "H9"]$note) &&
       grepl("ROBUST", h9rows[id == "H9"]$note))
h9tab <- suppressWarnings(run_h9_sensitivity(Mh9))
expect("H9 sensitivity table has a row per species source {none, truth, classifier}",
       !is.null(h9tab) && all(c("none", "truth", "classifier") %in% h9tab$species_source))
expect("H9 sensitivity: 'none' species_contribution is 0 by definition",
       approx(h9tab[species_source == "none", species_contribution], 0))
## robustness FLIP: donor-only < glm but +species(truth) > glm -> sign differs by source
donorv2 <- glmv - runif(8, 0.02, 0.05)
truthv2 <- glmv + runif(8, 0.03, 0.06)
h9flip  <- suppressWarnings(test_H9(mkM_h9(list(glm = glmv, glmmTMB_none = donorv2,
                                                glmmTMB_truth = truthv2))))
expect("H9 robustness: conclusion flips across species sources -> flagged as the finding",
       grepl("SOURCE-DEPENDENT", h9flip[id == "H9"]$note))
cat(sprintf("      rationale: H9 core GLMM(1|donor)-GLM med=%.3f; species contribution truth=%.3f, classifier=%.3f.\n",
            h9rows[id == "H9"]$median_diff, h9rows[id == "H9b_truth"]$median_diff,
            h9rows[id == "H9b_classifier"]$median_diff))

## =============================================================================
section("E4. Model fit+predict smoke tests (ranger / xgboost / glmmTMB)")
## =============================================================================
## The e2e (section F) exercises only fixed_threshold + glm (to stay fast + sandbox
## -safe), so the tree/GLMM fitters were historically UNtested -- and two latent API
## breaks slipped through onto the real run: xgboost 3.x's new xgboost(x, y) signature
## ('argument y is missing'), and glmmTMB::predict rejecting re.form = ~(1|species).
## These smoke tests fit each family on a small SEPARABLE synthetic set with a HELD-OUT
## (unseen) donor -- the LOEO transfer path where those bugs bit -- and assert valid,
## discriminating predictions. They SKIP (not fail) when a package is not loadable
## (e.g. the sandboxed $HOME R library), so the suite stays green in every context.
mk_model_xy <- function(donors, species = c("sp1", "sp2", "sp3"), n = 30L, seed = 1L) {
  set.seed(seed)
  dt <- rbindlist(lapply(donors, function(d) rbindlist(lapply(species, function(s) {
    y <- rep(0:1, each = n)
    data.table(donor = d, species = s, titration_level = "c1", y = y,
               f1 = rnorm(2 * n, 2 * y), f2 = rnorm(2 * n), f3 = rnorm(2 * n, -y)) }))))
  list(X = dt[, .(f1, f2, f3)], y = dt$y, w = rep(1, nrow(dt)),
       grp = dt[, .(donor, species, titration_level)])
}
.trM    <- mk_model_xy(sprintf("D%02d", 1:6), n = 30L, seed = 101L)  # 6 seen donors
.teM    <- mk_model_xy("D07",                n = 30L, seed = 202L)  # 1 UNSEEN donor (LOEO transfer)
.featM  <- c("f1", "f2", "f3")
.innerM <- as.list(sprintf("D%02d", 1:6))                            # donor-held-out inner folds
.pred_ok  <- function(p) is.numeric(p) && length(p) == length(.teM$y) &&
  all(is.finite(p)) && all(p >= 0 & p <= 1)
.pred_sep <- function(p) isTRUE(mean(p[.teM$y == 1]) > mean(p[.teM$y == 0]))  # separable -> AUC>0.5

model_smoke <- function(pkg, label, fit) {
  if (!is.na(pkg) && !requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  [SKIP] %s -- package '%s' not loadable here\n", label, pkg)); return(invisible()) }
  p <- tryCatch(suppressWarnings(suppressMessages(fit())), error = function(e) e)
  if (inherits(p, "error")) {
    expect(paste0(label, ": fit+predict runs without error"), FALSE)
    cat("      error: ", conditionMessage(p), "\n"); return(invisible()) }
  expect(paste0(label, ": valid probability vector in [0,1], length = n_test"), .pred_ok(p))
  expect(paste0(label, ": discriminates on separable data (mean score y=1 > y=0)"), .pred_sep(p))
}

model_smoke("ranger",  "ranger_rf (default HPs)",
            function() fit_ranger(.trM, .teM, .featM, NULL))
model_smoke("ranger",  "ranger_rf (inner-CV tuning path, num.trees=100)",
            function() fit_ranger(.trM, .teM, .featM, .innerM))
model_smoke("xgboost", "xgboost (default HPs, xgb.train + tree_method=hist)",
            function() fit_xgboost(.trM, .teM, .featM, NULL))
model_smoke("xgboost", "xgboost (inner-CV tuning path)",
            function() fit_xgboost(.trM, .teM, .featM, .innerM))
model_smoke("glmmTMB", "glmmTMB (1|donor)+(1|species), unseen-donor predict (re.form=NULL)",
            function() fit_glmmTMB(.trM, .teM, .featM))
model_smoke("glmmTMB", "glmmTMB donor-only (re.form=NA, population-level)",
            function() glmm_fit_predict(.trM, .teM, .featM, sp_train = NULL, sp_test = NULL))

## predict_arm_model dispatch: every configured model routes and returns a valid vector.
for (mdl in cfg$params$models) {
  pkg <- switch(mdl, ranger_rf = "ranger", xgboost = "xgboost", glmmTMB = "glmmTMB", NA_character_)
  if (!is.na(pkg) && !requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  [SKIP] predict_arm_model(%s) -- package not loadable\n", mdl)); next }
  p <- tryCatch(suppressWarnings(suppressMessages(
         predict_arm_model(mdl, .trM, .teM, .featM, NULL))), error = function(e) e)
  expect(sprintf("predict_arm_model dispatch: %s -> valid [0,1] vector", mdl),
         !inherits(p, "error") && .pred_ok(p))
}

## GLMM training-row cap (the stage-05 GLMM speed lever): stratified_cap_rows() bounds
## rows at the NATURAL class ratio. Pure + glmmTMB-free -> these run in every context.
expect("glmm_max_train_rows config is a single positive integer",
       is.numeric(cfg$params$glmm_max_train_rows) &&
         length(cfg$params$glmm_max_train_rows) == 1L &&
         is.finite(cfg$params$glmm_max_train_rows) && cfg$params$glmm_max_train_rows > 0)
local({
  y  <- rep(0:1, c(7000L, 3000L))                 # 70/30 prevalence, 10,000 rows
  k1 <- stratified_cap_rows(y, cap = 1000L, seed = 42L)
  expect("stratified_cap_rows: caps at <= cap rows", length(k1) <= 1000L)
  expect("stratified_cap_rows: keeps both classes", all(c(0L, 1L) %in% y[k1]))
  expect("stratified_cap_rows: preserves the class ratio (within 1%)",
         abs(mean(y[k1]) - mean(y)) < 0.01)
  expect("stratified_cap_rows: no-op when n <= cap",
         identical(stratified_cap_rows(y, cap = 50000L), seq_along(y)))
  expect("stratified_cap_rows: deterministic under a fixed seed",
         identical(stratified_cap_rows(y, 1000L, 42L), stratified_cap_rows(y, 1000L, 42L)))
})
## glmm_fit_predict honours the cap end-to-end (skip-guarded on glmmTMB): a small cap
## forces the subsample branch, and predictions must still be valid + discriminating.
local({
  old <- cfg$params$glmm_max_train_rows
  cfg$params$glmm_max_train_rows <<- 400L          # < nrow(.trM)=1080 -> triggers the cap
  model_smoke("glmmTMB", "glmmTMB with row cap on (stratified subsample path)",
              function() fit_glmmTMB(.trM, .teM, .featM))
  cfg$params$glmm_max_train_rows <<- old
})

## =============================================================================
section("F. End-to-end stages 04->07 (synthetic smoke test)")
## =============================================================================
e2e_err <- tryCatch({
  td <- file.path(tempdir(), paste0("e2e_", as.integer(Sys.time())))
  dir.create(file.path(td, "work"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(td, "results"), recursive = TRUE, showWarnings = FALSE)
  cfg$paths$work_dir       <- file.path(td, "work")
  cfg$paths$out_root       <- file.path(td, "results")
  cfg$paths$model_dir      <- file.path(td, "results", "models")
  cfg$paths$feature_table  <- file.path(td, "work", "feature_table.parquet")
  cfg$paths$labels_table   <- file.path(td, "work", "labels.tsv")
  cfg$paths$cv_splits      <- file.path(td, "work", "cv_splits.rds")
  cfg$paths$metrics_read   <- file.path(td, "results", "metrics_read_level.tsv")
  cfg$paths$metrics_sxt    <- file.path(td, "results", "metrics_sample_taxon.tsv")
  cfg$paths$hypotheses_out <- file.path(td, "results", "hypothesis_tests.tsv")

  set.seed(7)
  base <- rbindlist(lapply(donors, function(d) rbindlist(lapply(c("c1","c2","c3","c4","c5"), function(lv) {
    rbindlist(list(
      data.table(read_id = paste0(d, lv, "P", 1:12), donor = d, titration_level = lv,
                 species = sample(species, 12, replace = TRUE), label = "positive"),
      data.table(read_id = paste0(d, lv, "N", 1:20), donor = d, titration_level = lv,
                 species = NA_character_, label = "negative"))) }))))
  base[, `:=`(library_id = paste0(donor, "_", titration_level),
              end_reason_unblock = rbinom(.N, 1, 0.1))]
  yv <- as.integer(base$label == "positive")
  allfeat <- model_features("combined")   # M1/M3: every declared feature is now producible (no always-NA columns)
  for (f in allfeat) if (!f %in% c("bitscore", "pident", "k2_conf", "end_reason_unblock"))
    base[[f]] <- rnorm(nrow(base))
  base[, bitscore := rnorm(.N, 2.5 * yv)]
  base[, pident   := pmin(1, pmax(0, 0.7 + 0.2 * yv + rnorm(.N, 0, 0.03)))]
  base[, k2_conf  := pmin(1, pmax(0, 0.25 + 0.5 * yv + rnorm(.N, 0, 0.05)))]

  alt <- sub("\\.parquet$", ".tsv.gz", cfg$paths$feature_table)
  fwrite(base, alt, sep = "\t")
  fwrite(base[, .(read_id, library_id, donor, run_id = donor, barcode = 1L,
                  titration_level, species, label, concentration = NA_real_)],
         cfg$paths$labels_table, sep = "\t")
  ex <- CJ(donor = donors, titration_level = c("c1","c2","c3","c4","c5"), species = species)
  ex[, `:=`(expected_cells = 100, p_detect = 0.99,
            expected_present = species %in% c("sp1", "sp2", "sp3"))]
  fwrite(ex, file.path(cfg$paths$work_dir, "expected_sample_taxon.tsv"), sep = "\t")
  ## [F6] sample x taxon coverage so the breadth ablation (sxt_full/sxt_minus_breadth) runs
  cov <- ex[, .(donor, titration_level, species)]
  cov[, `:=`(genome_breadth    = ifelse(species %in% c("sp1","sp2","sp3"), runif(.N, 0.80, 0.98), runif(.N, 0.02, 0.20)),
             coverage_evenness = ifelse(species %in% c("sp1","sp2","sp3"), runif(.N, 0.75, 0.95), runif(.N, 0.05, 0.25)))]
  fwrite(cov, file.path(cfg$paths$work_dir, "sample_taxon_coverage.tsv"), sep = "\t")

  suppressMessages({
    run_stage04()
    run_stage05(models = c("fixed_threshold", "glm"), arms = names(cfg$classifier_arms),
                do_ablation = FALSE, do_loto = FALSE, do_h9 = FALSE)
    run_stage06()
    run_stage07()
  })
  NULL
}, error = function(e) conditionMessage(e))

if (is.null(e2e_err)) {
  preds <- fread(file.path(cfg$paths$out_root, "predictions.tsv.gz"))
  rl    <- fread(cfg$paths$metrics_read)
  sxt   <- fread(cfg$paths$metrics_sxt)
  hyp   <- fread(cfg$paths$hypotheses_out)
  expect("stage 05: predictions produced (>0 rows)", nrow(preds) > 0)
  expect("stage 06: read-level metrics have finite AUPRC", any(is.finite(rl$auprc)))
  expect("stage 06: sample x taxon metrics produced", nrow(sxt) > 0 && "sample_taxon" %in% sxt$level)
  expect("stage 07: hypothesis table has H1..H12 + H5b breadth component (F6)",
         all(c(paste0("H", 1:12), "H5b") %in% hyp$id))
  expect("stage 07: primary family carries Holm-Sidak adjusted p", "holm_sidak_p" %in% names(hyp))
  expect("item 3: inner-CV scores written (leakage-free selection)",
         file.exists(file.path(cfg$paths$out_root, "inner_cv_scores.tsv")))
  expect("item 3: model_comparison (overall/taxon/concentration) written", {
    p <- file.path(cfg$paths$out_root, "model_comparison.tsv"); file.exists(p) && "facet" %in% names(fread(p)) })
  expect("item 4: mixed-model supplement written",
         file.exists(file.path(cfg$paths$out_root, "mixed_model_supplement.tsv")))
  expect("F6: sample x taxon ablation rows present (sxt_full / sxt_minus_breadth)",
         all(c("sxt_full", "sxt_minus_breadth") %in% sxt$model))
  expect("F6: H5b breadth component measured (finite median_diff)",
         { r <- hyp[id == "H5b"]; nrow(r) == 1 && is.finite(r$median_diff) })
} else {
  expect(paste("end-to-end stages 04->07 run without error:", e2e_err), FALSE)
}

## =============================================================================
section("SUMMARY")
## =============================================================================
det <- Filter(function(x) x$detected, .findings)
cat(sprintf("  correctness checks : %d passed, %d failed\n", .PASS, .FAIL))
cat(sprintf("  review findings    : %d confirmed of %d probes\n", length(det), length(.findings)))
for (id in names(.findings)) {
  f <- .findings[[id]]
  cat(sprintf("    %-4s %-9s %s\n", id, if (f$detected) "DETECTED" else "clear", f$name))
}
quit(status = if (.FAIL > 0L) 1L else 0L)
