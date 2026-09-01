## =============================================================================
## utils.R  --  shared helpers: metrics, statistics, Poisson floor
## -----------------------------------------------------------------------------
## Sourced by the evaluation and hypothesis-testing stages. No side effects.
## =============================================================================

## Lazy package loader: load if present, else stop with an actionable message.
require_pkgs <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(sprintf("Missing R packages: %s\n  install.packages(c(%s))",
                 paste(missing, collapse = ", "),
                 paste(sprintf('\"%s\"', missing), collapse = ", ")),
         call. = FALSE)
  }
  invisible(TRUE)
}

## Per-item text progress bar for long loops (works under Rscript). Returns
## tick(i) to advance and done() to finish; no-op when total <= 0. [progress UX]
make_progress <- function(total, label = "") {
  total <- as.integer(total)
  if (is.na(total) || total <= 0)
    return(list(tick = function(i) invisible(NULL), done = function() invisible(NULL)))
  if (nzchar(label)) cat(sprintf("%s (%d items)\n", label, total))
  pb <- utils::txtProgressBar(min = 0, max = total, style = 3, width = 40)
  list(tick = function(i) utils::setTxtProgressBar(pb, i),
       done = function() { utils::setTxtProgressBar(pb, total); close(pb); cat("\n") })
}

## --- Classification metrics --------------------------------------------------
## All take a binary `labels` vector (1 = positive/true-Zymo, 0 = negative) and a
## numeric `scores` vector (higher = more positive). Reads in the 'ambiguous' or
## 'indeterminate' strata must be removed by the caller BEFORE scoring. [notes A,H]

.roc_pr_grid <- function(labels, scores) {
  stopifnot(length(labels) == length(scores))
  ord <- order(scores, decreasing = TRUE)
  y <- labels[ord]
  P <- sum(labels == 1); N <- sum(labels == 0)
  tp <- cumsum(y == 1)
  fp <- cumsum(y == 0)
  list(tp = tp, fp = fp, P = P, N = N,
       precision = tp / pmax(tp + fp, 1),
       recall    = if (P > 0) tp / P else rep(NA_real_, length(tp)))
}

## Area under the precision-recall curve (average precision). [note C]
auprc <- function(labels, scores) {
  g <- .roc_pr_grid(labels, scores)
  if (g$P == 0) return(NA_real_)
  rec <- c(0, g$recall)
  prec <- c(1, g$precision)
  sum(diff(rec) * prec[-1])  # step-wise average precision
}

## Precision achievable at a target recall: the highest precision among all
## thresholds whose recall >= target (upper envelope). [note C]
precision_at_recall <- function(labels, scores, target_recall = 0.80) {
  g <- .roc_pr_grid(labels, scores)
  if (g$P == 0) return(NA_real_)
  ok <- which(g$recall >= target_recall)
  if (!length(ok)) return(NA_real_)
  max(g$precision[ok])
}

## FDR (= 1 - precision) at the score threshold that first reaches target recall. [note C]
fdr_at_recall <- function(labels, scores, target_recall = 0.80) {
  p <- precision_at_recall(labels, scores, target_recall)
  if (is.na(p)) NA_real_ else 1 - p
}

## Brier score for calibrated probabilities in [0,1]. Lower is better. [H12]
brier_score <- function(labels, probs) mean((probs - labels)^2)

## --- Multiplicity correction -------------------------------------------------
## Step-down Holm-Sidak adjusted p-values. With m<=6 primary tests at n=8 folds
## the exact Wilcoxon floor (2/2^8 = 0.0078) can still clear alpha. [note J]
holm_sidak <- function(pvals) {
  m <- length(pvals)
  o <- order(pvals)
  p_sorted <- pvals[o]
  adj_sorted <- 1 - (1 - p_sorted)^(m - seq_len(m) + 1)
  adj_sorted <- cummax(adj_sorted)          # enforce monotonicity
  adj <- numeric(m); adj[o] <- pmin(adj_sorted, 1)
  adj
}

## --- Paired comparison across the 8 LOEO folds -------------------------------
## Exact paired Wilcoxon signed-rank on per-fold metric differences. [project_plan]
paired_wilcox <- function(x, y, exact = TRUE) {
  d <- x - y
  d <- d[d != 0]
  if (length(d) < 1) return(list(p.value = NA_real_, statistic = NA_real_, n = 0L))
  wt <- suppressWarnings(stats::wilcox.test(x, y, paired = TRUE, exact = exact))
  list(p.value = unname(wt$p.value), statistic = unname(wt$statistic), n = length(d))
}

## Donor-level (cluster) bootstrap CI for the median paired difference.
## Resamples folds/donors with replacement -> percentile CI. [project_plan]
donor_bootstrap_ci <- function(diffs, n_boot = 2000L, conf = 0.95, seed = 1L) {
  diffs <- diffs[is.finite(diffs)]
  if (length(diffs) < 2) return(c(est = median(diffs), lo = NA_real_, hi = NA_real_))
  set.seed(seed)
  boot_med <- replicate(n_boot, median(sample(diffs, replace = TRUE)))
  a <- (1 - conf) / 2
  c(est = median(diffs),
    lo  = unname(quantile(boot_med, a, names = FALSE)),
    hi  = unname(quantile(boot_med, 1 - a, names = FALSE)))
}

## Donor-clustered mixed-effects test for the n = 8 power problem. [item 4 / F3]
## Fits  value ~ method + (1 | donor)  over per-(donor x stratum) metric rows,
## using WITHIN-donor replication (the titration levels) to raise the effective
## sample size above the 8-fold exact-Wilcoxon sign floor, while a (1 | donor)
## random intercept absorbs donor clustering (and the run confound). Returns the
## method fixed-effect estimate + p (positive = the second factor level is
## better). Falls back to a donor fixed-effect lm if lmerTest is unavailable.
mixed_effect_test <- function(long, value = "auprc") {
  ok <- is.finite(long[[value]]) & !is.na(long$method) & !is.na(long$donor)
  long <- long[ok, , drop = FALSE]
  if (nrow(long) < 4 || length(unique(long$method)) < 2 || length(unique(long$donor)) < 3)
    return(list(estimate = NA_real_, p.value = NA_real_, n = nrow(long), engine = "none"))
  long$method <- factor(long$method)
  if (requireNamespace("lmerTest", quietly = TRUE)) {
    fit <- try(lmerTest::lmer(stats::as.formula(sprintf("%s ~ method + (1 | donor)", value)),
                              data = long), silent = TRUE)
    if (!inherits(fit, "try-error")) {
      co <- summary(fit)$coefficients
      ix <- grep("^method", rownames(co))
      if (length(ix))
        return(list(estimate = unname(co[ix[1], 1]), p.value = unname(co[ix[1], ncol(co)]),
                    n = nrow(long), engine = "lmerTest"))
    }
  }
  fit <- stats::lm(stats::as.formula(sprintf("%s ~ method + donor", value)), data = long)
  co <- summary(fit)$coefficients; ix <- grep("^method", rownames(co))
  list(estimate = if (length(ix)) unname(co[ix[1], 1]) else NA_real_,
       p.value  = if (length(ix)) unname(co[ix[1], ncol(co)]) else NA_real_,
       n = nrow(long), engine = "lm")
}

## --- Poisson expected-copy-number floor [note H] -----------------------------
## E[N] = total_cells * rel_abundance * extraction_eff * fraction_loaded
## P(>=1 genome) = 1 - exp(-E[N]).  Species x level below poisson_p_min are
## routed to the 'indeterminate' stratum and excluded from recall accounting.
poisson_detection_prob <- function(total_cells, rel_abundance,
                                    extraction_eff = 1, fraction_loaded = 1) {
  en <- total_cells * rel_abundance * extraction_eff * fraction_loaded
  list(expected_genomes = en, p_at_least_one = 1 - exp(-en))
}

## Given a per-species relative-abundance table + per-level cell counts, return
## the species x level detectability grid and the indeterminate mask. [note H]
build_poisson_floor <- function(coa, cells_per_level, extraction_eff,
                                fraction_loaded, p_min = 0.95) {
  levels <- names(cells_per_level)
  grid <- expand.grid(species = coa$species, level = levels,
                      stringsAsFactors = FALSE)
  grid$rel_abundance <- coa$rel_abundance[match(grid$species, coa$species)]
  grid$total_cells   <- cells_per_level[grid$level]
  pd <- poisson_detection_prob(grid$total_cells, grid$rel_abundance,
                               extraction_eff, fraction_loaded)
  grid$expected_genomes <- pd$expected_genomes
  grid$p_detect         <- pd$p_at_least_one
  grid$indeterminate    <- is.na(grid$p_detect) | grid$p_detect < p_min
  grid
}

## --- Sequence-complexity helpers [note F] ------------------------------------
gc_fraction <- function(seq) {
  s <- toupper(seq)
  (nchar(gsub("[^GC]", "", s))) / pmax(nchar(s), 1)
}

## Longest homopolymer run as a fraction of read length.
homopolymer_frac <- function(seq) {
  runs <- gregexpr("(.)\\1*", toupper(seq), perl = TRUE)[[1]]
  if (runs[1] == -1) return(0)
  max(attr(runs, "match.length")) / max(nchar(seq), 1)
}

## Simple DUST-like low-complexity score from triplet frequencies (higher =
## lower complexity). A drop-in for symmetric DUST; swap for sdust if available. [note F / OI 2]
dust_score <- function(seq) {
  s <- toupper(seq); L <- nchar(s)
  if (L < 3) return(0)
  triplets <- substring(s, 1:(L - 2), 3:L)
  counts <- table(triplets)
  sum(counts * (counts - 1)) / (2 * (L - 2))
}

## Shannon entropy (bits) of a discrete vector -- used for taxonomic entropy of
## the top-N BLAST hits. [note F]
shannon_entropy <- function(labels) {
  if (!length(labels)) return(0)
  p <- prop.table(table(labels))
  -sum(p * log2(p))
}

## Genome breadth (fraction of a subject genome covered) and coverage evenness
## (1 - Gini of per-base depth). Strongest real-vs-artefact discriminator. [note F,G]
coverage_evenness <- function(per_base_depth) {
  d <- sort(per_base_depth)
  n <- length(d)
  if (n == 0 || sum(d) == 0) return(0)
  gini <- (2 * sum(seq_len(n) * d) / (n * sum(d))) - (n + 1) / n
  1 - gini
}

## Weighted Gini via the Lorenz curve (values x with segment weights w). [R3]
weighted_gini <- function(x, w) {
  keep <- w > 0; x <- x[keep]; w <- w[keep]
  if (!length(x)) return(0)
  o <- order(x); x <- x[o]; w <- w[o]
  W <- sum(w); sw <- sum(w * x)
  if (sw == 0) return(0)                       # no coverage anywhere
  px <- c(0, cumsum(w) / W)                    # cumulative population fraction
  vy <- c(0, cumsum(w * x) / sw)               # cumulative value fraction
  B <- sum((px[-1] - px[-length(px)]) * (vy[-1] + vy[-length(vy)]) / 2)  # area under Lorenz
  max(0, min(1, 1 - 2 * B))
}

## Breadth + evenness for ONE reference of length glen from a set of alignment
## intervals [start, end) (equivalent to `samtools depth` but computed directly
## from the minimap2 PAF, no BAM needed). Depth is run-length-encoded via a
## sweep line; evenness = 1 - weighted Gini over per-base depth incl. zero gaps. [R3 / note F,G]
coverage_stats_from_intervals <- function(starts, ends, glen) {
  glen <- as.numeric(glen[1])
  if (is.na(glen) || glen <= 0) return(list(genome_breadth = NA_real_, coverage_evenness = NA_real_))
  s <- pmax(0, pmin(glen, as.numeric(starts)))
  e <- pmax(0, pmin(glen, as.numeric(ends)))
  ok <- e > s; s <- s[ok]; e <- e[ok]
  if (!length(s)) return(list(genome_breadth = 0, coverage_evenness = 0))
  bp <- data.table(pos = c(s, e), delta = c(rep(1L, length(s)), rep(-1L, length(e))))
  bp <- bp[, .(delta = sum(delta)), by = pos][order(pos)]
  depth_after <- cumsum(bp$delta)                     # depth just AFTER each breakpoint
  seg_start <- c(0, bp$pos); seg_end <- c(bp$pos, glen)
  seg_dep   <- c(0, depth_after)                      # depth before first breakpoint = 0
  seg_len   <- seg_end - seg_start
  keep <- seg_len > 0; seg_len <- seg_len[keep]; seg_dep <- seg_dep[keep]
  list(genome_breadth   = sum(seg_len[seg_dep > 0]) / glen,
       coverage_evenness = 1 - weighted_gini(seg_dep, seg_len))
}
