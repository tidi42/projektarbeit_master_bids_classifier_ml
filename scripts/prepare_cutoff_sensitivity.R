## =============================================================================
## prepare_cutoff_sensitivity.R  --  ground-truth cutoff sensitivity analysis [F9]
## -----------------------------------------------------------------------------
## The whole supervised target hinges on two fixed cutoffs:
##   gt_min_identity = 0.90   gt_min_coverage = 0.80
## which were never validated for a hac r10.4.1 basecaller (true-Zymo identity mode
## is ~0.95-0.98). This standalone analysis STRESS-TESTS them so key results can be
## shown stable across a defensible cutoff range (addresses reporting_3 F9).
##
## Method (cutoff-INDEPENDENT confident set):
##   1. Collect CONFIDENT Zymo reads at the high-titre level(s) (cfg$gt_calib_levels):
##      a read whose best Zymo minimap2 hit beats its best human hit
##      (zymo_score > human_score + gt_human_margin) -- the "beats human" rule that
##      defines a candidate positive, WITHOUT using the identity/coverage cutoffs, so
##      the set does not depend on the values under test.
##   2. Report the empirical identity + coverage distribution of that set.
##   3. Sweep a grid of (identity x coverage) cutoffs and report the FRACTION of the
##      confident set retained at each -- i.e. how sensitive the positive label set
##      (the input to every downstream endpoint) is to the cutoff.
##   4. Recommend the strictest cutoffs that still retain >= cfg$gt_retain_frac.
##
## Run AFTER stage 01 (needs work/<lib>/gt_zymo.paf + gt_human_grch38.paf), under
## the default/fixed profile (it reads only the shared work/ PAFs):
##   Rscript scripts/prepare_cutoff_sensitivity.R
## Outputs (all in the SHARED work/ dir -- profile-independent):
##   work/cutoff_zymo_distribution.tsv   identity/coverage quantiles of confident set
##   work/cutoff_sensitivity.tsv         retained fraction per (identity, coverage)
##   work/cutoff_recommended.tsv         the single recommended (identity, coverage)
##                                       pair -- CONSUMED by GT_PROFILE=calculated as
##                                       the second run's ground-truth cutoffs.
## =============================================================================

suppressWarnings(suppressMessages(library(data.table)))
.sd <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1])))
  else "/home/tdinse/master_timo/projektarbeit/projektarbeit_ml/scripts"
})
if (!exists("cfg")) { source(file.path(.sd, "00_config.R")); source(file.path(.sd, "utils.R")) }
## reuse read_paf() / best_per_read() / map_contig_to_species() from stage 02
## (sourcing does NOT auto-run it: the sys.nframe() guard only fires under Rscript).
if (!exists("read_paf")) source(file.path(.sd, "02_ground_truth_labels.R"))

## ---- (1) collect confident Zymo reads at the calibration level(s) -----------
## Confident = best Zymo hit beats best human hit. This is INDEPENDENT of the
## identity/coverage cutoffs under test, so it is a fair yardstick for them.
collect_confident_zymo <- function(levels = cfg$params$gt_calib_levels) {
  ss <- fread(cfg$paths$sample_sheet)
  if ("titration_level" %in% names(ss)) ss <- ss[titration_level %in% levels]
  if (!nrow(ss)) {
    message("  no libraries at calibration level(s): ", paste(levels, collapse = ", "))
    return(data.table())
  }
  rows <- vector("list", nrow(ss))
  for (i in seq_len(nrow(ss))) {
    lib <- as.list(ss[i])
    out_dir <- file.path(cfg$paths$work_dir, lib$library_id)
    zpaf <- file.path(out_dir, "gt_zymo.paf")
    hpaf <- file.path(out_dir, "gt_human_grch38.paf")
    if (!file.exists(zpaf)) next
    z <- best_per_read(read_paf(zpaf))
    if (!nrow(z)) next
    h  <- best_per_read(read_paf(hpaf))
    h2 <- best_per_read(read_paf(file.path(out_dir, "competitor_t2t.paf")))  # T2T
    ## human_score = max(hg38, T2T) -- SAME competition as stage-01 depletion and
    ## stage-02 labelling, so this confident-Zymo yardstick matches the real labels. [R2]
    z[, human_score := pmax(.best_score(h, qname), .best_score(h2, qname))]
    conf <- z[score > human_score + cfg$params$gt_human_margin]
    if (!nrow(conf)) next
    conf[, `:=`(library_id = lib$library_id,
                donor = if (!is.null(lib$donor)) lib$donor else NA_character_,
                titration_level = if (!is.null(lib$titration_level)) lib$titration_level else NA_character_,
                species = map_contig_to_species(tname))]
    rows[[i]] <- conf[, .(library_id, donor, titration_level, species, read_id = qname,
                          identity, coverage, zymo_score = score, human_score)]
  }
  out <- rbindlist(rows, fill = TRUE)
  if (!nrow(out)) message("  found gt_zymo.paf files but no reads beat human at these level(s).")
  out
}

## ---- (2) empirical identity + coverage distribution (long format) -----------
zymo_distribution <- function(conf,
                              probs = c(0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99)) {
  if (!nrow(conf)) return(data.table())
  mk <- function(x, nm) {
    q <- stats::quantile(x, probs, na.rm = TRUE)
    data.table(metric = nm, quantile = paste0("p", probs * 100),
               value = as.numeric(q), mean = mean(x, na.rm = TRUE), n = sum(is.finite(x)))
  }
  rbindlist(list(mk(conf$identity, "identity"), mk(conf$coverage, "coverage")))
}

## ---- (3) sweep the (identity x coverage) cutoff grid ------------------------
## For each pair, how many confident-Zymo reads would still be LABELLED positive.
cutoff_sweep <- function(conf,
                         id_grid = cfg$params$gt_identity_grid,
                         cov_grid = cfg$params$gt_coverage_grid) {
  if (!nrow(conf)) return(data.table())
  ntot <- nrow(conf)
  grid <- CJ(gt_min_identity = id_grid, gt_min_coverage = cov_grid)
  grid[, n_confident := ntot]
  grid[, n_retained := mapply(function(ic, cc)
        sum(conf$identity >= ic & conf$coverage >= cc, na.rm = TRUE),
        gt_min_identity, gt_min_coverage)]
  grid[, frac_retained := n_retained / pmax(n_confident, 1L)]
  grid[, is_current := abs(gt_min_identity - cfg$params$gt_min_identity) < 1e-9 &
                       abs(gt_min_coverage - cfg$params$gt_min_coverage) < 1e-9]
  grid[order(gt_min_identity, gt_min_coverage)]
}

## ---- (4) recommend the strictest cutoffs retaining >= gt_retain_frac --------
recommend_cutoffs <- function(sweep, retain_frac = cfg$params$gt_retain_frac) {
  if (!nrow(sweep)) return(data.table())
  ok <- sweep[frac_retained >= retain_frac]
  if (!nrow(ok)) return(sweep[which.max(frac_retained)][1])
  ok[order(-gt_min_identity, -gt_min_coverage)][1]   # strictest = highest id then cov
}

## ---- driver -----------------------------------------------------------------
run_cutoff_sensitivity <- function() {
  cfg_init_dirs()
  conf <- collect_confident_zymo()
  if (!nrow(conf)) {
    message("Cutoff sensitivity [F9]: no confident Zymo reads. Run stage 01 first ",
            "(needs work/<lib>/gt_zymo.paf + gt_human_grch38.paf) and ensure libraries ",
            "exist at level(s): ", paste(cfg$params$gt_calib_levels, collapse = ", "), ".")
    return(invisible(NULL))
  }
  dist  <- zymo_distribution(conf)
  sweep <- cutoff_sweep(conf)
  rec   <- recommend_cutoffs(sweep)

  dist_path  <- file.path(cfg$paths$work_dir, "cutoff_zymo_distribution.tsv")
  sweep_path <- file.path(cfg$paths$work_dir, "cutoff_sensitivity.tsv")
  rec_path   <- cfg$paths$gt_recommended
  fwrite(dist, dist_path, sep = "\t")
  fwrite(sweep, sweep_path, sep = "\t")
  ## Persist the single recommended (identity, coverage) pair so the SECOND run
  ## (GT_PROFILE=calculated) can pick it up as its ground-truth cutoffs. [F9]
  if (nrow(rec))
    fwrite(rec[, .(gt_min_identity, gt_min_coverage, frac_retained, n_confident, n_retained)],
           rec_path, sep = "\t")

  message(sprintf("Cutoff sensitivity [F9]: %d confident-Zymo reads from %d librar(ies) at level(s) %s.",
                  nrow(conf), uniqueN(conf$library_id), paste(cfg$params$gt_calib_levels, collapse = ", ")))
  message("  identity/coverage quantiles -> ", dist_path)
  message("  cutoff sweep                -> ", sweep_path)
  if (nrow(rec)) message("  recommended cutoffs         -> ", rec_path,
                         "  (consumed by GT_PROFILE=calculated)")
  cur <- sweep[is_current == TRUE]
  if (nrow(cur))
    message(sprintf("  CURRENT cutoffs id>=%.2f cov>=%.2f retain %.1f%% of confident-Zymo reads.",
                    cur$gt_min_identity[1], cur$gt_min_coverage[1], 100 * cur$frac_retained[1]))
  if (nrow(rec))
    message(sprintf("  strictest cutoffs retaining >= %.0f%%: id>=%.2f cov>=%.2f (retain %.1f%%).",
                    100 * cfg$params$gt_retain_frac, rec$gt_min_identity, rec$gt_min_coverage,
                    100 * rec$frac_retained))
  print(dist)
  print(sweep)
  invisible(list(conf = conf, distribution = dist, sweep = sweep, recommended = rec))
}

if (sys.nframe() == 0L) run_cutoff_sensitivity()
