## =============================================================================
## 04_cv_splits.R  --  nested cross-validation design
## -----------------------------------------------------------------------------
## Two outer schemes, both donor-grouped so a donor's barcodes are NEVER split
## across train/test [guardrail]:
##
##   LOEO  leave-one-donor-out   -> 8 outer folds. The primary scheme.
##   LOTO  leave-one-taxon-out   -> tests generalisation to UNSEEN taxa; the
##                                  honest test of transfer, per note B / H10.
##
## Nested CV: each outer training set is further split into `inner_folds` DONOR
## groups for hyperparameter tuning. Only OUTER-fold performance is reported
## [guardrail]. If donors were spread across flow cells, folds are additionally
## checked against the run confounder. [OI 8]
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

## Assign a set of groups (donors) into k inner folds as evenly as possible.
group_kfold <- function(groups, k, seed = cfg$params$seed) {
  set.seed(seed)
  g <- sample(unique(groups))
  split(g, cut(seq_along(g), breaks = min(k, length(g)), labels = FALSE))
}

## Warn if donor and run are confounded (each donor on its own flow cell). [OI 8]
check_run_confounder <- function(meta) {
  tab <- unique(meta[, .(donor, run_id)])
  per_run_donors <- tab[, .(n = uniqueN(donor)), by = run_id]
  if (all(per_run_donors$n == 1)) {
    message("NOTE: each run has exactly one donor -> run is CONFOUNDED with donor. ",
            "LOEO already removes it, but you cannot separately estimate a run effect. [OI 8]")
  }
}

build_loeo <- function(donors, meta) {
  lapply(seq_along(donors), function(i) {
    test_d <- donors[i]
    train_d <- setdiff(donors, test_d)
    list(fold = i, scheme = "LOEO",
         test_donors = test_d,
         train_donors = train_d,
         inner = group_kfold(train_d, cfg$params$inner_folds, cfg$params$seed + i))
  })
}

build_loto <- function(species, donors) {
  species <- species[!is.na(species) & nzchar(species)]
  donors  <- sort(unique(donors))
  ## [F1] Donors are partitioned into two groups so that NEITHER positives nor
  ## negatives are shared between the LOTO train and test sets. The held-out
  ## taxon's positives are removed from training; at test the model must call them
  ## from taxon-agnostic features only, scored against negatives it has NOT seen.
  ## The test group rotates across folds so both donor groups are exercised.
  grpA <- donors[seq_along(donors) %% 2L == 1L]
  grpB <- setdiff(donors, grpA)
  lapply(seq_along(species), function(i) {
    test_donors <- if (i %% 2L == 1L) grpA else grpB
    if (!length(test_donors) || !length(setdiff(donors, test_donors)))
      test_donors <- donors[1]                         # degenerate guard (<= 1 donor)
    train_donors <- setdiff(donors, test_donors)
    list(fold = i, scheme = "LOTO",
         test_species = species[i],
         train_species = setdiff(species, species[i]),
         test_donors = test_donors, train_donors = train_donors,
         inner = group_kfold(train_donors, cfg$params$inner_folds, cfg$params$seed + 100L + i))
  })
}

run_stage04 <- function() {
  cfg_init_dirs()
  labels <- fread(cfg$paths$labels_table)
  meta <- unique(labels[, .(donor, run_id, titration_level, barcode)])
  check_run_confounder(meta)

  donors  <- intersect(cfg$params$donors, unique(labels$donor))
  species <- sort(unique(labels[label == "positive", species]))

  splits <- list(
    loeo = build_loeo(donors, meta),
    loto = build_loto(species, donors),
    meta = list(donors = donors, species = species,
                inner_folds = cfg$params$inner_folds, seed = cfg$params$seed)
  )
  saveRDS(splits, cfg$paths$cv_splits)
  message(sprintf("Stage 04 complete: %d LOEO folds, %d LOTO folds -> %s",
                  length(splits$loeo), length(splits$loto), cfg$paths$cv_splits))
  invisible(splits)
}

## Materialise a train/test row mask for a given outer fold. Used by stage 05/06.
fold_masks <- function(dt, fold) {
  if (fold$scheme == "LOEO") {
    list(train = dt$donor %in% fold$train_donors & dt$label != "ambiguous" & dt$label != "indeterminate",
         test  = dt$donor %in% fold$test_donors  & dt$label != "ambiguous" & dt$label != "indeterminate")
  } else { # LOTO [F1]: taxon held out of TRAIN positives; donors partitioned so
           # neither positives nor negatives are shared between train and test.
    is_pos <- dt$label == "positive"
    is_neg <- dt$label == "negative"
    in_test  <- dt$donor %in% fold$test_donors
    in_train <- dt$donor %in% fold$train_donors
    list(train = in_train & ((is_pos & dt$species != fold$test_species) | is_neg),
         test  = in_test  & ((is_pos & dt$species == fold$test_species) | is_neg))
  }
}

if (sys.nframe() == 0L) run_stage04()
