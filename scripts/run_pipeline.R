## =============================================================================
## run_pipeline.R  --  orchestrator for the Zymo-in-human read-classification
##                     ML hypothesis-testing pipeline
## -----------------------------------------------------------------------------
## Usage:
##   Rscript scripts/run_pipeline.R                 # print status + open items, then run all
##   Rscript scripts/run_pipeline.R --status        # only print config status + open items
##   Rscript scripts/run_pipeline.R --from 03       # run from stage 03 onwards
##   Rscript scripts/run_pipeline.R --only 05,06,07 # run just those stages
##
## The pipeline is a straight line:
##   01 external tools -> 02 labels -> 03 features -> 04 CV splits
##      -> 05 train -> 06 evaluate -> 07 hypothesis tests
## =============================================================================

SCRIPTS <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1])))
  else "/home/tdinse/master_timo/projektarbeit/projektarbeit_ml/scripts"
})
THIS_FILE <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) normalizePath(sub("^--file=", "", a[1])) else file.path(SCRIPTS, "run_pipeline.R")
})

## -----------------------------------------------------------------------------
## GROUND-TRUTH PROFILE selection / dispatch. MUST precede sourcing 00_config.R,
## because config reads GT_PROFILE to resolve the cutoffs + the per-run out folder.
##   --gt fixed        run 1 only : gt id>=0.90 cov>=0.80  -> results/gt_fixed_.../
##   --gt calculated   run 2 only : cutoffs from prepare_cutoff_sensitivity.R
##   --gt both         (DEFAULT)   run 1 THEN run 2, each in its own results folder
## The GT_PROFILE env var takes precedence over --gt and pins THIS process to a
## single profile -- that is how the per-profile child processes below run exactly
## one profile each without re-dispatching.
## -----------------------------------------------------------------------------
.rp_args <- commandArgs(TRUE)
.rp_opt  <- function(flag, default = NA_character_) {
  i <- which(.rp_args == flag)
  if (length(i) && i[1] < length(.rp_args)) .rp_args[i[1] + 1L] else default
}
.rp_strip <- function(args, flag) {
  i <- which(args == flag); if (length(i)) args[-c(i[1], i[1] + 1L)] else args
}

if (nzchar(Sys.getenv("GT_PROFILE", ""))) {
  ## launched as a per-profile child (or the user preset the env): run this one inline.
  Sys.setenv(GT_PROFILE = tolower(Sys.getenv("GT_PROFILE")))
} else {
  .rp_gt       <- tolower(.rp_opt("--gt", "both"))
  .rp_profiles <- if (.rp_gt == "both") c("fixed", "calculated") else .rp_gt
  if (!all(.rp_profiles %in% c("fixed", "calculated")))
    stop("--gt must be one of: fixed, calculated, both", call. = FALSE)

  ## Multi-profile request -> launch one clean child Rscript per profile, then exit.
  ## (--status is single-shot: fall through and print status for the fixed profile.)
  if (length(.rp_profiles) > 1L && !("--status" %in% .rp_args)) {
    .rscript  <- file.path(R.home("bin"), "Rscript")
    .pass     <- .rp_strip(.rp_args, "--gt")
    .rec_file <- file.path(dirname(SCRIPTS), "work", "cutoff_recommended.tsv")
    for (i in seq_along(.rp_profiles)) {
      p     <- .rp_profiles[i]
      child <- .pass
      ## stage 01 (alignments/BLAST/Kraken2) is cutoff-INDEPENDENT and shared, so
      ## only the FIRST run needs it; later runs start at stage 02 UNLESS the user
      ## already constrained the stage set with --from/--only.
      if (i > 1L && !any(c("--from", "--only") %in% child)) child <- c(child, "--from", "02")
      ## the 'calculated' run needs the recommended cutoffs -> build them (from the
      ## shared work/ PAFs) if the first run has not already produced them.
      if (identical(p, "calculated") && !file.exists(.rec_file)) {
        cat("\n>> building recommended cutoffs (prepare_cutoff_sensitivity.R) ...\n")
        system2(.rscript, shQuote(file.path(SCRIPTS, "prepare_cutoff_sensitivity.R")),
                env = "GT_PROFILE=fixed")
      }
      cat("\n############################################################\n")
      cat(sprintf("##  GROUND-TRUTH RUN %d/%d  --  GT_PROFILE=%s\n", i, length(.rp_profiles), p))
      cat("############################################################\n")
      st <- system2(.rscript, c(shQuote(THIS_FILE), child), env = paste0("GT_PROFILE=", p))
      if (!identical(as.integer(st), 0L))
        message("   ground-truth run '", p, "' exited with status ", st)
    }
    quit(save = "no")
  }
  ## single profile (or --status): run inline under that profile.
  Sys.setenv(GT_PROFILE = .rp_profiles[1])
}

source(file.path(SCRIPTS, "00_config.R"))
source(file.path(SCRIPTS, "utils.R"))

STAGES <- list(
  "01" = list(file = "01_external_tools.R",      fun = "run_stage01", name = "external tools (QC/minimap2/BLAST/Kraken2)"),
  "02" = list(file = "02_ground_truth_labels.R", fun = "run_stage02", name = "ground-truth labels + Poisson floor"),
  "03" = list(file = "03_build_features.R",      fun = "run_stage03", name = "taxon-agnostic feature table"),
  "04" = list(file = "04_cv_splits.R",           fun = "run_stage04", name = "nested LOEO/LOTO CV splits"),
  "05" = list(file = "05_train_models.R",        fun = "run_stage05", name = "train models (GLM/RF/XGB/GLMM)"),
  "06" = list(file = "06_evaluate.R",            fun = "run_stage06", name = "evaluate (read + sample x taxon)"),
  "07" = list(file = "07_hypothesis_tests.R",    fun = "run_stage07", name = "hypothesis tests H1-H12")
)

## -----------------------------------------------------------------------------
## Status report: hardcoded-path gaps + unresolved OPEN ITEMS. Always printed.
## -----------------------------------------------------------------------------
print_status <- function() {
  cat("\n========================================================================\n")
  cat("  PIPELINE STATUS\n")
  cat("========================================================================\n")

  cat(sprintf("\n-- GROUND-TRUTH PROFILE: %s   (gt_min_identity >= %.2f, gt_min_coverage >= %.2f)\n",
              cfg$params$gt_profile, cfg$params$gt_min_identity, cfg$params$gt_min_coverage))
  cat(sprintf("   this run's output folder : %s\n", cfg$paths$out_root))

  cat("\n-- Hardcoded paths still missing / placeholder (SECTION 1 of 00_config.R):\n")
  gaps <- cfg_missing_paths()
  if (nrow(gaps)) print(gaps, row.names = FALSE) else cat("   (none -- all input paths resolve)\n")

  n_all  <- nrow(cfg$open_items)
  n_open <- sum(cfg$open_items$status != "resolved")
  cat(sprintf("\n-- OPEN ITEMS: %d open, %d resolved (SECTION 5; full triage in open_items_v2.md):\n",
              n_open, n_all - n_open))
  oi <- cfg_open_items()
  if (nrow(oi)) {
    for (i in seq_len(nrow(oi)))
      cat(sprintf("   [%-3s] (%s)  %s\n", oi$id[i], oi$blocks_stage[i], oi$item[i]))
  } else cat("   (none open)\n")

  cat("\n-- Hypotheses under test:\n")
  h <- cfg$hypotheses
  for (i in seq_len(nrow(h)))
    cat(sprintf("   %-9s %s%s\n", h$family[i], h$test[i],
                if (!h$active[i]) "   [INACTIVE]" else ""))
  cat("\n========================================================================\n\n")
}

## -----------------------------------------------------------------------------
## Progress reporting: a run plan, a per-stage progress bar + timing, and a
## completion breakdown so the user sees what is running and what is done.
## -----------------------------------------------------------------------------
fmt_dur <- function(secs) {
  if (!is.finite(secs)) "?"
  else if (secs < 60)   sprintf("%.1fs", secs)
  else if (secs < 3600) sprintf("%.1fm", secs / 60)
  else                  sprintf("%.2fh", secs / 3600)
}
stage_bar <- function(done, total, width = 22L) {
  f <- if (total > 0) round(width * done / total) else 0L
  sprintf("[%s%s] %3d%%", strrep("#", f), strrep(".", width - f),
          if (total > 0) round(100 * done / total) else 0L)
}
print_plan <- function(keys) {
  cat(sprintf("\n-- RUN PLAN: %d stage(s) --------------------------------------------\n", length(keys)))
  for (i in seq_along(keys))
    cat(sprintf("   %d/%d  stage %s  %s\n", i, length(keys), keys[i], STAGES[[keys[i]]]$name))
  cat("--------------------------------------------------------------------\n")
}

run_stage <- function(key, idx, n) {
  st <- STAGES[[key]]
  if (is.null(st)) stop("unknown stage: ", key)
  cat(sprintf("\n%s  (%d/%d)  RUNNING  stage %s : %s\n",
              stage_bar(idx - 1L, n), idx, n, key, st$name))
  t0 <- Sys.time()
  ## source into the global env so each stage's functions persist for later ones
  source(file.path(SCRIPTS, st$file), local = FALSE)
  get(st$fun)()
  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("%s  (%d/%d)  DONE     stage %s : %s  [%s]\n",
              stage_bar(idx, n), idx, n, key, st$name, fmt_dur(dt)))
  dt
}

main <- function() {
  args <- commandArgs(TRUE)
  print_status()
  if ("--status" %in% args) return(invisible())

  keys <- names(STAGES)
  if ("--from" %in% args) {
    from <- args[which(args == "--from") + 1]
    keys <- keys[keys >= from]
  }
  if ("--only" %in% args) {
    only <- strsplit(args[which(args == "--only") + 1], ",")[[1]]
    keys <- intersect(keys, trimws(only))
  }

  cfg_init_dirs()
  print_plan(keys)
  t_all <- Sys.time(); times <- setNames(rep(NA_real_, length(keys)), keys)
  for (i in seq_along(keys)) times[keys[i]] <- run_stage(keys[i], i, length(keys))

  cat("\n==================== PIPELINE COMPLETE ====================\n")
  for (i in seq_along(keys))
    cat(sprintf("   [x] stage %s  %-42s %8s\n", keys[i], STAGES[[keys[i]]]$name, fmt_dur(times[keys[i]])))
  cat(sprintf("   total: %s\n", fmt_dur(as.numeric(difftime(Sys.time(), t_all, units = "secs")))))
  cat("==========================================================\n")
  message("\nPrimary results: ", cfg$paths$hypotheses_out)
}

if (!interactive()) main()
