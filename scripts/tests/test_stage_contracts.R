#!/usr/bin/env Rscript
## =============================================================================
## test_stage_contracts.R -- isolated leakage/configuration tests for stages 01-07
## =============================================================================

.this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
TESTS_DIR <- if (length(.this_file) && !is.na(.this_file) && nzchar(.this_file))
  dirname(normalizePath(.this_file)) else
  file.path(getwd(), "scripts", "tests")
SCRIPTS_DIR <- dirname(TESTS_DIR)
TEST_PROJECT_ROOT <- dirname(SCRIPTS_DIR)

Sys.setenv(GT_PROFILE = "fixed")
source(file.path(TESTS_DIR, "stage_test_parameters.R"))

suppressWarnings(suppressMessages({
  source(file.path(SCRIPTS_DIR, "00_config.R"))
  source(file.path(SCRIPTS_DIR, "utils.R"))
  source(file.path(SCRIPTS_DIR, "01_external_tools.R"))
  source(file.path(SCRIPTS_DIR, "02_ground_truth_labels.R"))
  source(file.path(SCRIPTS_DIR, "03_build_features.R"))
  source(file.path(SCRIPTS_DIR, "04_cv_splits.R"))
  source(file.path(SCRIPTS_DIR, "05_train_models.R"))
  source(file.path(SCRIPTS_DIR, "06_evaluate.R"))
  source(file.path(SCRIPTS_DIR, "07_hypothesis_tests.R"))
}))
suppressWarnings(suppressMessages(library(data.table)))

PARAMS <- PIPELINE_TEST_PARAMETERS
REGISTRY <- as.data.table(PARAMS$registry)
TEST_RESULTS <- data.table()

deep_copy <- function(x) unserialize(serialize(x, NULL))

values_equal <- function(actual, expected, tolerance = PARAMS$numeric_tolerance) {
  if (is.list(expected)) {
    if (!is.list(actual) || length(actual) != length(expected)) return(FALSE)
    if (!is.null(names(expected)) && !identical(names(actual), names(expected))) return(FALSE)
    return(all(vapply(seq_along(expected), function(i)
      values_equal(actual[[i]], expected[[i]], tolerance), logical(1))))
  }
  if (is.numeric(expected)) {
    if (!is.numeric(actual) || length(actual) != length(expected)) return(FALSE)
    if (!identical(names(actual), names(expected))) return(FALSE)
    both_na <- is.na(actual) & is.na(expected)
    one_na <- xor(is.na(actual), is.na(expected))
    if (any(one_na)) return(FALSE)
    return(all(both_na | abs(actual - expected) <= tolerance))
  }
  identical(actual, expected)
}

display_value <- function(x) {
  if (inherits(x, "error")) return(paste0("ERROR: ", conditionMessage(x)))
  paste(capture.output(dput(x)), collapse = " ")
}

record_error <- function(id, error) {
  reg_index <- base::match(id, REGISTRY[["id"]])
  reg <- REGISTRY[reg_index, ]
  TEST_RESULTS <<- rbind(TEST_RESULTS, data.table(
    id = id, stage = reg$stage, kind = reg$kind, description = reg$description,
    expected = display_value(PARAMS$expectations[[id]]),
    actual = paste0("ERROR: ", conditionMessage(error)), matches_expected = FALSE,
    status = "ERROR", detail = conditionMessage(error)
  ), fill = TRUE)
  cat(sprintf("  [ERROR]    %-34s %s\n", id, conditionMessage(error)))
}

record_test <- function(id, actual, detail = "") {
  if (!(id %in% REGISTRY$id)) stop("Unregistered test id: ", id, call. = FALSE)
  if (id %in% TEST_RESULTS$id) stop("Duplicate test id: ", id, call. = FALSE)
  reg_index <- base::match(id, REGISTRY[["id"]])
  reg <- REGISTRY[reg_index, ]
  value <- tryCatch(force(actual), error = function(e) e)
  if (inherits(value, "error")) {
    record_error(id, value)
    return(invisible(value))
  }
  expected <- PARAMS$expectations[[id]]
  matched <- values_equal(value, expected)
  status <- if (reg$kind == "finding") {
    if (isTRUE(value)) "DETECTED" else "CLEAR"
  } else if (matched) "PASS" else "FAIL"
  TEST_RESULTS <<- rbind(TEST_RESULTS, data.table(
    id = id, stage = reg$stage, kind = reg$kind, description = reg$description,
    expected = display_value(expected), actual = display_value(value),
    matches_expected = matched, status = status, detail = as.character(detail)
  ), fill = TRUE)
  cat(sprintf("  [%-8s] %-34s %s\n", status, id, reg$description))
  invisible(value)
}

run_stage_group <- function(stage, fun) {
  cat(sprintf("\n== Stage %s ==\n", stage))
  stage_error <- tryCatch({ fun(); NULL }, error = function(e) e)
  missing <- setdiff(REGISTRY$id[REGISTRY$stage == stage], TEST_RESULTS$id)
  if (length(missing)) {
    if (is.null(stage_error)) stage_error <- simpleError("test was not executed")
    for (id in missing) record_error(id, stage_error)
  }
}

set_global_cfg <- function(value) assign("cfg", value, envir = .GlobalEnv)

paf_row <- function(qname, qlen = 100, qstart = 0, qend = 95, strand = "+",
                    tname = "zc1", tlen = 5000, tstart = 0, tend = 95,
                    matches = 90, blocklen = 95, mapq = 60, tags = character()) {
  c(qname, qlen, qstart, qend, strand, tname, tlen, tstart, tend,
    matches, blocklen, mapq, tags)
}

write_paf_rows <- function(path, rows) {
  writeLines(vapply(rows, function(x) paste(x, collapse = "\t"), character(1)), path)
  invisible(path)
}

write_blast_rows <- function(path, rows) {
  writeLines(vapply(rows, function(x) paste(x, collapse = "\t"), character(1)), path)
  invisible(path)
}

stage01_tests <- function() {
  old_cfg <- deep_copy(cfg)
  old_run_cmd <- get("run_cmd", envir = .GlobalEnv)
  old_path <- Sys.getenv("PATH")
  on.exit({
    set_global_cfg(old_cfg)
    assign("run_cmd", old_run_cmd, envir = .GlobalEnv)
    Sys.setenv(PATH = old_path)
  }, add = TRUE)

  td <- tempfile("stage01_")
  dir.create(td, recursive = TRUE)
  tc <- deep_copy(old_cfg)
  tc$paths$work_dir <- file.path(td, "work")
  tc$paths$zymo_refs_fasta <- file.path(td, "zymo.fa")
  tc$paths$human_grch38 <- file.path(td, "human.fa")
  tc$paths$human_t2t_chm13 <- file.path(td, "t2t.fa")
  tc$paths$blast_db <- file.path(td, "blastdb")
  tc$paths$kraken2_db <- file.path(td, "kraken_db")
  set_global_cfg(tc)

  ragged <- file.path(td, "ragged.paf")
  write_paf_rows(ragged, list(
    paf_row("r1", matches = 80, tags = "cg:Z:80M"),
    paf_row("r2", matches = 70, tags = c("cg:Z:70M", "NM:i:2")),
    paf_row("r1", matches = 95, blocklen = 100, tags = c("cg:Z:95M", "NM:i:1", "AS:i:95")),
    paf_row("r3", matches = 60)
  ))
  best <- .paf_best(ragged)
  record_test("S01_RAGGED_PAF", list(
    n_best = nrow(best), r1_score = best[read_id == "r1", score]
  ))

  mk <- function(ids, scores) data.table(read_id = ids, score = scores)
  human_ids <- human_like_reads(
    mk(c("r_human", "r_zymo", "r_tie"), c(100, 40, 70)),
    mk("r_t2t", 90),
    mk(c("r_human", "r_zymo", "r_tie"), c(80, 120, 70))
  )
  record_test("S01_HUMAN_COMPETITION", sort(human_ids))

  bad_sheet <- file.path(td, "bad_sample_sheet.tsv")
  fwrite(data.table(library_id = "L1", donor = "D1"), bad_sheet, sep = "\t")
  tc <- deep_copy(cfg); tc$paths$sample_sheet <- bad_sheet; set_global_cfg(tc)
  schema_error <- tryCatch({ read_sample_sheet(); NULL }, error = function(e) e)
  record_test("S01_SAMPLE_SHEET_SCHEMA",
              inherits(schema_error, "error") && grepl("missing columns", conditionMessage(schema_error)))

  fake_bin <- file.path(td, "bin")
  dir.create(fake_bin)
  fake_seqkit <- file.path(fake_bin, "seqkit")
  writeLines(c("#!/bin/sh", "exit 0"), fake_seqkit)
  Sys.chmod(fake_seqkit, mode = "0755")
  Sys.setenv(PATH = paste(fake_bin, old_path, sep = .Platform$path.sep))

  calls <- list()
  spy <- function(bin, args, stdout_file = "") {
    calls[[length(calls) + 1L]] <<- list(bin = as.character(bin),
                                         args = as.character(args),
                                         stdout = as.character(stdout_file))
    invisible(0L)
  }
  assign("run_cmd", spy, envir = .GlobalEnv)
  lib <- list(library_id = "L1", donor = "D1", run_id = "run-1",
              fastq = file.path(td, "raw.fastq"), seq_summary = NA_character_)
  out_dir <- file.path(td, "work", "L1")
  dir.create(out_dir, recursive = TRUE)
  depleted <- file.path(td, "depleted.fastq")
  groundtruth_align(lib, out_dir, "minimap2", "samtools", threads = 7L)
  suppressMessages(qc_library(lib, out_dir, reads_fastq = depleted))
  blast_classify(lib, out_dir, "blastn", "makeblastdb", threads = 7L,
                 top_n = 5L, reads_fastq = depleted)
  kraken2_classify(lib, out_dir, "kraken2", threads = 7L, reads_fastq = depleted)

  call_text <- vapply(calls, function(x) paste(c(x$bin, x$args), collapse = " "), character(1))
  mm_calls <- call_text[grepl("^minimap2 ", call_text)]
  qc_call <- call_text[grepl("seqkit .*fx2tab", call_text)]
  fq2fa_call <- call_text[grepl("seqkit .*fq2fa", call_text)]
  blast_call <- call_text[grepl("^blastn ", call_text)]
  kraken_call <- call_text[grepl("^kraken2 ", call_text)]
  command_ok <- length(mm_calls) == 3L && all(grepl(lib$fastq, mm_calls, fixed = TRUE)) &&
    all(grepl("-t 7", mm_calls, fixed = TRUE)) && length(qc_call) == 1L &&
    grepl("--only-id", qc_call, fixed = TRUE) && grepl(depleted, qc_call, fixed = TRUE) &&
    length(fq2fa_call) == 1L && grepl(depleted, fq2fa_call, fixed = TRUE) &&
    length(blast_call) == 1L && grepl("-task megablast", blast_call, fixed = TRUE) &&
    grepl("-max_hsps 1", blast_call, fixed = TRUE) &&
    grepl("-max_target_seqs 5", blast_call, fixed = TRUE) &&
    grepl("-num_threads 7", blast_call, fixed = TRUE) &&
    length(kraken_call) == 1L && grepl(depleted, kraken_call, fixed = TRUE) &&
    grepl("--memory-mapping", kraken_call, fixed = TRUE) &&
    grepl("--threads 7", kraken_call, fixed = TRUE)
  record_test("S01_COMMAND_WIRING", command_ok)

  incomplete <- file.path(td, "incomplete")
  dir.create(incomplete)
  terminal_files <- c("reads.fasta", "blastn.tsv", "kraken2.out", "kraken2.report")
  for (name in terminal_files) writeLines("x", file.path(incomplete, name))
  file.create(file.path(incomplete, "blastn.tsv"), showWarnings = FALSE)
  record_test("S01_EMPTY_OUTPUT_INCOMPLETE", !lib_complete(incomplete))

  resume_only <- file.path(td, "resume_only")
  dir.create(resume_only)
  for (name in terminal_files) writeLines("x", file.path(resume_only, name))
  record_test("F01_RESUME_OMITS_PAF_QC",
              lib_complete(resume_only) &&
                !file.exists(file.path(resume_only, "gt_zymo.paf")) &&
                !file.exists(file.path(resume_only, "read_qc.tsv")))
}

stage02_tests <- function() {
  old_cfg <- deep_copy(cfg)
  on.exit(set_global_cfg(old_cfg), add = TRUE)

  td <- tempfile("stage02_")
  work <- file.path(td, "work")
  ref <- file.path(td, "ref")
  dir.create(work, recursive = TRUE)
  dir.create(ref, recursive = TRUE)
  tc <- deep_copy(old_cfg)
  tc$paths$project_root <- td
  tc$paths$work_dir <- work
  tc$paths$out_root <- file.path(td, "results")
  tc$paths$leakage_table <- file.path(td, "results", "leakage.tsv")
  tc$params$ambiguous_bed <- file.path(ref, "ambiguous.bed")
  tc$params$gt_min_identity <- 0.90
  tc$params$gt_min_coverage <- 0.80
  tc$params$gt_human_margin <- 0
  set_global_cfg(tc)
  dir.create(cfg$paths$out_root, recursive = TRUE)

  fwrite(data.table(contig = "zc1", species = "Zymo_species"),
         file.path(ref, "zymo_contig2species.tsv"), sep = "\t")
  fwrite(data.table(V1 = "zc1", V2 = 1000L, V3 = 1100L),
         cfg$params$ambiguous_bed, sep = "\t", col.names = FALSE)

  lib_dir <- file.path(work, "L1")
  dir.create(lib_dir)
  survivor_ids <- c("r_positive", "r_human", "r_tie", "r_lowid", "r_lowcov",
                    "r_ambiguous", "r_none")
  qc <- data.table(id = survivor_ids, length = 100L, avg_qual = 20, gc = 50)
  setnames(qc, c("#id", "length", "avg.qual", "GC(%)"))
  fwrite(qc, file.path(lib_dir, "read_qc.tsv"), sep = "\t")
  write_paf_rows(file.path(lib_dir, "gt_zymo.paf"), list(
    paf_row("r_positive", matches = 90),
    paf_row("r_human", matches = 90),
    paf_row("r_tie", matches = 90),
    paf_row("r_lowid", matches = 80, blocklen = 100),
    paf_row("r_lowcov", qend = 50, tend = 50, matches = 48, blocklen = 50),
    paf_row("r_ambiguous", tstart = 1000, tend = 1095, matches = 90),
    paf_row("r_paf_only", matches = 95,
            tags = c("cg:Z:95M", "NM:i:0", "AS:i:95"))
  ))
  write_paf_rows(file.path(lib_dir, "gt_human_grch38.paf"), list(
    paf_row("r_positive", tname = "chr1", matches = 20),
    paf_row("r_human", tname = "chr1", matches = 100, blocklen = 100),
    paf_row("r_tie", tname = "chr1", matches = 90),
    paf_row("r_ambiguous", tname = "chr1", matches = 10)
  ))
  write_paf_rows(file.path(lib_dir, "competitor_t2t.paf"), list(
    paf_row("r_human", tname = "chrT2T", matches = 105, blocklen = 110)
  ))
  lib <- list(library_id = "L1", donor = "D1", titration_level = "c1",
              run_id = "R1", barcode = "BC01", concentration = 10)
  labels <- suppressMessages(label_library(lib))
  expected_labels <- PARAMS$expectations$S02_LABEL_ASSIGNMENT
  observed_labels <- setNames(labels$label, labels$read_id)[names(expected_labels)]
  record_test("S02_LABEL_ASSIGNMENT", observed_labels)
  record_test("S02_QC_DEFINES_UNIVERSE", sort(labels$read_id))

  coa <- file.path(td, "coa.tsv")
  fwrite(data.table(species = c("high", "low"), rel_abundance = c(0.5, 0.0001)),
         coa, sep = "\t")
  tc <- deep_copy(cfg)
  tc$paths$zymo_coa <- coa
  tc$params$extraction_eff <- 1
  tc$params$fraction_loaded <- 1
  tc$params$poisson_p_min <- 0.95
  set_global_cfg(tc)
  poisson_input <- data.table(
    read_id = c("high", "low", "background"),
    species = c("high", "low", "low"),
    label = c("positive", "positive", "negative"),
    concentration = 10, titration_level = "c1"
  )
  poisson_out <- apply_poisson_floor(poisson_input)
  poisson_actual <- setNames(poisson_out$label, poisson_out$read_id)[
    names(PARAMS$expectations$S02_POISSON_ROUTING)]
  record_test("S02_POISSON_ROUTING", poisson_actual)

  sheet <- file.path(td, "sample_sheet.tsv")
  fwrite(data.table(library_id = c("L1", "L2"), donor = c("D1", "D2"),
                    titration_level = c("c1", "c1"), concentration = c(10, 20)),
         sheet, sep = "\t")
  tc <- deep_copy(cfg); tc$paths$sample_sheet <- sheet; set_global_cfg(tc)
  expectation <- expected_sample_taxon()
  record_test("S02_EXPECTATION_GRID", list(
    n_rows = nrow(expectation),
    n_expected_present = sum(expectation$expected_present)
  ))

  leakage_labels <- data.table(
    run_id = c(rep("R1", 5), "R2"),
    titration_level = c(rep("negative", 4), "c4", "c4"),
    label = c("positive", "negative", "negative", "negative", "positive", "positive")
  )
  leakage <- suppressMessages(estimate_leakage(leakage_labels))
  leakage_actual <- setNames(leakage$leakage_upper_bound, leakage$run_id)[c("R1", "R2")]
  record_test("S02_LEAKAGE_ESTIMATE", leakage_actual)

  fallback_dir <- file.path(work, "L_no_qc")
  dir.create(fallback_dir)
  write_paf_rows(file.path(fallback_dir, "gt_zymo.paf"), list(paf_row("r_zymo")))
  write_paf_rows(file.path(fallback_dir, "gt_human_grch38.paf"),
                 list(paf_row("r_human", tname = "chr1")))
  write_paf_rows(file.path(fallback_dir, "competitor_t2t.paf"), list())
  fallback <- suppressMessages(label_library(list(
    library_id = "L_no_qc", donor = "D1", titration_level = "c1",
    run_id = "R1", barcode = "BC02", concentration = 10
  )))
  record_test("F02_MISSING_QC_LOSES_UNALIGNED", !("r_unaligned" %in% fallback$read_id))

  tc <- deep_copy(cfg); tc$params$ambiguous_bed <- file.path(td, "missing.bed"); set_global_cfg(tc)
  mask_input <- data.table(zymo_tname = "zc1", zymo_tstart = 1000L, zymo_tend = 1095L)
  missing_mask <- suppressMessages(flag_ambiguous_regions(mask_input))
  record_test("F02_MISSING_BED_DISABLES_MASK", identical(missing_mask, FALSE))
}

stage03_tests <- function() {
  old_cfg <- deep_copy(cfg)
  on.exit(set_global_cfg(old_cfg), add = TRUE)

  td <- tempfile("stage03_")
  work <- file.path(td, "work")
  dir.create(work, recursive = TRUE)
  tc <- deep_copy(old_cfg)
  tc$paths$project_root <- td
  tc$paths$work_dir <- work
  tc$params$top_n_hits <- 2L
  set_global_cfg(tc)

  blast_path <- file.path(td, "blast.tsv")
  write_blast_rows(blast_path, list(
    c("r1", "s1", 1, 99, 100, 0, 0, 1, 100, 1, 100, "1e-50", 100, 100, 5000),
    c("r1", "s2", 2, 95, 90, 5, 0, 1, 90, 1, 90, "1e-30", 70, 100, 4000)
  ))
  taxres <- function(ids) data.table(species = paste0("sp", ids), genus = paste0("g", ids))
  blast <- parse_blast(blast_path, top_n = 2L, taxres = taxres)
  record_test("S03_BLAST_FEATURES", list(
    margin = blast[read_id == "r1", bitscore_margin_species],
    n_species = blast[read_id == "r1", n_species_topN],
    entropy = blast[read_id == "r1", tax_entropy_topN],
    subject_len = blast[read_id == "r1", subject_genome_len]
  ))

  kraken_path <- file.path(td, "kraken.out")
  writeLines(c("C\tr1\t10\t100\t10:7 20:3", "U\tr2\t0\t100\t0:0"), kraken_path)
  kraken <- parse_kraken2(kraken_path)
  kraken_actual <- c(
    r1_conf = kraken[read_id == "r1", k2_conf],
    r1_distinct = kraken[read_id == "r1", k2_distinct_minimizers],
    r2_conf = kraken[read_id == "r2", k2_conf],
    r2_distinct = kraken[read_id == "r2", k2_distinct_minimizers]
  )
  record_test("S03_KRAKEN_FEATURES", kraken_actual)

  lib_dir <- file.path(work, "L3")
  dir.create(lib_dir)
  write_blast_rows(file.path(lib_dir, "blastn.tsv"), list(
    c("r_blast", "s1", 1, 99, 100, 0, 0, 1, 100, 1, 100, "1e-50", 100, 100, 5000)
  ))
  write_paf_rows(file.path(lib_dir, "competitor_t2t.paf"), list(
    paf_row("r_competitor", tname = "chrT2T", tags = c("cg:Z:90M", "NM:i:1"))
  ))
  writeLines("C\tr_kraken\t10\t100\t10:8 0:2", file.path(lib_dir, "kraken2.out"))
  qc <- data.table(id = "r_qc", length = 101L, avg_qual = 18, gc = 51)
  setnames(qc, c("#id", "length", "avg.qual", "GC(%)"))
  fwrite(qc, file.path(lib_dir, "read_qc.tsv"), sep = "\t")
  union_features <- suppressMessages(features_for_library(list(library_id = "L3"), taxres))
  record_test("S03_OUTER_UNION", sort(union_features$read_id))

  blast_arm <- model_features("blast_only")
  kraken_arm <- model_features("kraken2_only")
  combined_arm <- model_features("combined")
  arm_ok <- !any(cfg$feature_blocks$kraken2 %in% blast_arm) &&
    !any(unlist(cfg$feature_blocks[c("blast_core", "blast_margin", "human_competitor")]) %in% kraken_arm) &&
    setequal(combined_arm, union(blast_arm, kraken_arm))
  record_test("S03_ARM_ISOLATION", arm_ok)

  tc <- deep_copy(cfg)
  tc$feature_blocks$read_qc <- unique(c(tc$feature_blocks$read_qc,
                                        "label", "donor", "species"))
  set_global_cfg(tc)
  injected <- intersect(model_features("combined"), c("label", "donor", "species"))
  record_test("F03_FORBIDDEN_CONFIG_INJECTION", length(injected) > 0L)
}

stage04_tests <- function() {
  old_cfg <- deep_copy(cfg)
  on.exit(set_global_cfg(old_cfg), add = TRUE)
  tc <- deep_copy(old_cfg); tc$params$inner_folds <- 3L; tc$params$seed <- PARAMS$seed
  set_global_cfg(tc)

  donors <- paste0("D", 1:4)
  rows <- rbindlist(lapply(donors, function(d) rbindlist(list(
    data.table(read_id = paste0(d, "_p1"), donor = d, library_id = paste0(d, "_L1"),
               barcode = "BC01", run_id = paste0("R", d), species = "sp1", label = "positive"),
    data.table(read_id = paste0(d, "_p2"), donor = d, library_id = paste0(d, "_L2"),
               barcode = "BC02", run_id = paste0("R", d), species = "sp2", label = "positive"),
    data.table(read_id = paste0(d, "_n"), donor = d, library_id = paste0(d, "_L0"),
               barcode = "BC00", run_id = paste0("R", d), species = NA_character_, label = "negative"),
    data.table(read_id = paste0(d, "_a"), donor = d, library_id = paste0(d, "_L1"),
               barcode = "BC01", run_id = paste0("R", d), species = "sp1", label = "ambiguous"),
    data.table(read_id = paste0(d, "_i"), donor = d, library_id = paste0(d, "_L1"),
               barcode = "BC01", run_id = paste0("R", d), species = "sp1", label = "indeterminate")
  ))))
  loeo <- build_loeo(donors, unique(rows[, .(donor, run_id)]))
  test_counts <- integer(nrow(rows))
  loeo_ok <- TRUE
  for (fold in loeo) {
    masks <- fold_masks(rows, fold)
    test_counts <- test_counts + as.integer(masks$test)
    loeo_ok <- loeo_ok && !any(masks$train & masks$test) &&
      !length(intersect(unique(rows$donor[masks$train]), unique(rows$donor[masks$test]))) &&
      !any(rows$label[masks$train | masks$test] %in% c("ambiguous", "indeterminate")) &&
      setequal(unique(rows$donor[masks$test]), fold$test_donors)
  }
  eligible <- rows$label %in% c("positive", "negative")
  loeo_ok <- loeo_ok && all(test_counts[eligible] == 1L) && all(test_counts[!eligible] == 0L)
  record_test("S04_LOEO_CONTRACT", loeo_ok)

  inner_ok <- all(vapply(loeo, function(fold) {
    inner <- unlist(fold$inner, use.names = FALSE)
    setequal(inner, fold$train_donors) && !anyDuplicated(inner) &&
      !length(intersect(inner, fold$test_donors)) && length(fold$inner) <= cfg$params$inner_folds
  }, logical(1)))
  record_test("S04_INNER_DONOR_CONTRACT", inner_ok)

  loto <- build_loto(c("sp1", "sp2"), donors)
  loto_ok <- all(vapply(loto, function(fold) {
    masks <- fold_masks(rows, fold)
    !any(masks$train & masks$test) &&
      !length(intersect(unique(rows$donor[masks$train]), unique(rows$donor[masks$test]))) &&
      !any(rows$label[masks$train | masks$test] %in% c("ambiguous", "indeterminate")) &&
      !any(rows$label[masks$train] == "positive" & rows$species[masks$train] == fold$test_species) &&
      all(rows$species[masks$test & rows$label == "positive"] == fold$test_species)
  }, logical(1)))
  record_test("S04_LOTO_CONTRACT", loto_ok)
  record_test("S04_DETERMINISTIC_SPLITS",
              identical(loeo, build_loeo(donors, unique(rows[, .(donor, run_id)]))))

  shared <- copy(rows); shared[, run_id := "R_SHARED"]
  split_run <- any(vapply(loeo, function(fold) {
    masks <- fold_masks(shared, fold)
    length(intersect(unique(shared$run_id[masks$train]), unique(shared$run_id[masks$test]))) > 0L
  }, logical(1)))
  record_test("F04_SHARED_RUN_SPLIT", split_run)
}

stage05_tests <- function() {
  old_cfg <- deep_copy(cfg)
  old_predict <- get("predict_arm_model", envir = .GlobalEnv)
  on.exit({
    set_global_cfg(old_cfg)
    assign("predict_arm_model", old_predict, envir = .GlobalEnv)
  }, add = TRUE)

  train_x <- data.table(f = c(1, 2, 3, NA_real_))
  med <- make_imputer(train_x)
  test_x <- apply_imputer(data.table(f = NA_real_), med)
  record_test("S05_TRAIN_ONLY_IMPUTER", c(train_median = med[["f"]],
                                           imputed_test = test_x$f))

  y_train <- c(0L, 0L, 0L, 1L, 1L, 1L)
  tr <- list(X = data.table(signal = c(0.1, 0.2, 0.3, 0.7, 0.8, 0.9),
                            noise = c(1, 0, 1, 0, 1, 0)),
             y = y_train, w = rep(1, length(y_train)))
  te <- list(X = data.table(signal = c(0.15, 0.85), noise = c(0, 1)), y = c(0L, 1L))
  baseline <- suppressWarnings(fit_fixed_threshold(tr, te, c("signal", "noise")))
  record_test("S05_FIXED_BASELINE", list(
    chosen_feature = attr(baseline, "chosen_feature"),
    probabilities_valid = all(is.finite(baseline) & baseline >= 0 & baseline <= 1),
    positive_ranked_higher = baseline[2] > baseline[1]
  ))

  groups <- data.table(donor = c("D1", "D1", rep("D2", 4)),
                       titration_level = "c1")
  weights <- balance_weights(groups, unit = "donor")
  totals <- data.table(donor = groups$donor, weight = weights)[, sum(weight), by = donor]$V1
  record_test("S05_BALANCE_WEIGHTS",
              abs(diff(totals)) < PARAMS$numeric_tolerance && abs(mean(weights) - 1) < PARAMS$numeric_tolerance)

  tune_dt <- data.table(donor = rep(paste0("D", 1:4), each = 4), y = rep(0:1, 8), f = seq_len(16))
  tune_tr <- list(X = tune_dt[, .(f)], y = tune_dt$y,
                  grp = tune_dt[, .(donor)], w = rep(1, nrow(tune_dt)))
  tune_inner_folds <- as.list(paste0("D", 1:4))
  best <- tune_inner(tune_tr, "f", tune_inner_folds,
                     function(a, b, row) if (isTRUE(row$good)) b$y else 1 - b$y,
                     data.frame(good = c(FALSE, TRUE)))
  record_test("S05_INNER_GRID_SELECTION", isTRUE(best$good))

  td <- tempfile("stage05_")
  dir.create(file.path(td, "results"), recursive = TRUE)
  tc <- deep_copy(old_cfg)
  tc$paths$out_root <- file.path(td, "results")
  tc$paths$model_dir <- file.path(td, "results", "models")
  tc$paths$feature_table <- file.path(td, "feature_table.parquet")
  tc$paths$cv_splits <- file.path(td, "cv_splits.rds")
  tc$feature_blocks$test_block <- "bitscore"
  tc$classifier_arms <- list(combined = "test_block")
  tc$subject_props <- character()
  tc$params$inner_folds <- 2L
  set_global_cfg(tc)

  donors <- paste0("D", 1:4)
  ft <- rbindlist(lapply(donors, function(d) data.table(
    read_id = paste0(d, "_", 1:4), donor = d, species = c("sp1", NA, "sp2", NA),
    titration_level = "c1", label = rep(c("positive", "negative"), 2),
    end_reason_unblock = 0L, bitscore = c(0.9, 0.1, 0.8, 0.2)
  )))
  fwrite(ft, sub("\\.parquet$", ".tsv.gz", cfg$paths$feature_table), sep = "\t")
  splits <- list(loeo = build_loeo(donors, unique(ft[, .(donor, run_id = donor)])), loto = list())
  saveRDS(splits, cfg$paths$cv_splits)

  calls <- list()
  recorder <- function(model, tr, te, features, inner) {
    calls[[length(calls) + 1L]] <<- list(train = sort(unique(tr$grp$donor)),
                                         test = sort(unique(te$grp$donor)))
    rep(0.5, length(te$y))
  }
  assign("predict_arm_model", recorder, envir = .GlobalEnv)
  predictions <- suppressMessages(run_stage05(models = "stub", arms = "combined",
                                               do_ablation = FALSE, do_loto = FALSE, do_h9 = FALSE))
  isolation <- length(calls) > 0L && all(vapply(calls, function(x)
    !length(intersect(x$train, x$test)), logical(1)))
  record_test("S05_OOF_DONOR_ISOLATION", isolation)
  coverage <- predictions[, .N, by = read_id]
  record_test("S05_OOF_COVERAGE",
              nrow(coverage) == nrow(ft) && all(coverage$N == 1L) && setequal(coverage$read_id, ft$read_id))

  tc <- deep_copy(cfg); tc$feature_blocks$test_block <- c("bitscore", "missing_configured"); set_global_cfg(tc)
  selected <- intersect(model_features("combined"), names(ft))
  record_test("F05_MISSING_FEATURE_DROPPED",
              "missing_configured" %in% model_features("combined") && !("missing_configured" %in% selected))

  assign("predict_arm_model", function(...) stop("synthetic model failure"), envir = .GlobalEnv)
  failed_run <- tryCatch(
    suppressMessages(run_stage05(models = "broken", arms = "combined",
                                 do_ablation = FALSE, do_loto = FALSE, do_h9 = FALSE)),
    error = function(e) e
  )
  record_test("F05_MODEL_FAILURE_SWALLOWED",
              !inherits(failed_run, "error") && nrow(failed_run) > 0L && all(is.na(failed_run$score)))
}

stage06_tests <- function() {
  old_cfg <- deep_copy(cfg)
  on.exit(set_global_cfg(old_cfg), add = TRUE)

  td <- tempfile("stage06_")
  work <- file.path(td, "work")
  dir.create(work, recursive = TRUE)
  tc <- deep_copy(old_cfg)
  tc$paths$work_dir <- work
  tc$params$recall_target <- c(0.5, 1.0)
  tc$params$report_truncated <- c(TRUE, FALSE)
  tc$params$poisson_p_min <- 0.95
  set_global_cfg(tc)

  hand <- metric_rows(c(1L, 1L, 0L, 0L), c(0.9, 0.8, 0.2, 0.1), targets = c(0.5, 1.0))
  record_test("S06_HAND_CALCULATED_METRICS", list(
    auprc = unique(hand$auprc), precision = unname(hand$prec_at_recall),
    fdr = unname(hand$fdr_at_recall)
  ))

  trunc_preds <- data.table(
    read_id = paste0("r", 1:3), donor = "D1", titration_level = "c1",
    species = c("sp1", NA, NA), y = c(1L, 0L, 0L), score = c(0.9, 0.8, 0.1),
    end_reason_unblock = c(0L, 1L, 0L), arm = "combined", model = "glm",
    scheme = "LOEO", fold = 1L
  )
  trunc_metrics <- read_level_metrics(trunc_preds)
  trunc_actual <- c(
    included = trunc_metrics[stratum == "all" & truncated_included == TRUE & recall_target == 0.5, n_neg][1],
    excluded = trunc_metrics[stratum == "all" & truncated_included == FALSE & recall_target == 0.5, n_neg][1]
  )
  storage.mode(trunc_actual) <- "integer"
  record_test("S06_TRUNCATED_COUNTS", trunc_actual)

  expected <- data.table(
    donor = "D1", titration_level = "c1", species = c("spA", "spB", "spC"),
    expected_cells = c(10, 0, 10), p_detect = c(0.99, 0, 0.99),
    expected_present = c(TRUE, FALSE, TRUE)
  )
  fwrite(expected, file.path(work, "expected_sample_taxon.tsv"), sep = "\t")
  sxt_preds <- data.table(
    read_id = c("rA1", "rA2", "rB1", "rN"), donor = "D1", titration_level = "c1",
    species = c("spA", "spA", "spB", NA_character_), y = c(0L, 0L, 1L, 0L),
    score = c(0.9, 0.8, 0.7, 0.1), end_reason_unblock = 0L,
    arm = "combined", model = "glm", scheme = "LOEO", fold = 1L
  )
  sxt <- suppressMessages(aggregate_sample_taxon(sxt_preds))
  truth <- setNames(sxt$y, sxt$species)[c("spA", "spB", "spC")]
  storage.mode(truth) <- "integer"
  record_test("S06_SAMPLE_TAXON_TRUTH", truth)
  observed_min <- min(sxt$max_score[sxt$species != "spC"], na.rm = TRUE)
  record_test("S06_MISSING_TAXON_IS_FN",
              nrow(sxt[species == "spC" & y == 1L]) == 1L &&
                sxt[species == "spC", max_score] < observed_min)
  record_test("S06_DEPTH_NORMALIZATION",
              sxt[species == "spA", reads_above_thr_per_million])

  no_truth_work <- file.path(td, "no_truth")
  dir.create(no_truth_work)
  tc <- deep_copy(cfg); tc$paths$work_dir <- no_truth_work; set_global_cfg(tc)
  threshold_fixture <- function(other_scores) rbind(
    data.table(read_id = c("a1", "a2"), donor = "A", titration_level = "c1",
               species = "sp", y = 1L, score = c(0.7, 0.8)),
    data.table(read_id = c("b1", "b2"), donor = "B", titration_level = "c1",
               species = "sp", y = 0L, score = other_scores)
  )[, `:=`(arm = "combined", model = "glm", scheme = "LOEO", fold = match(donor, c("A", "B")))]
  count_a <- function(p) suppressMessages(aggregate_sample_taxon(p))[
    donor == "A", n_reads_above_thr]
  pooled_changed <- count_a(threshold_fixture(c(0.1, 0.2))) !=
    count_a(threshold_fixture(c(100, 101)))
  record_test("F06_POOLED_OOF_THRESHOLD", isTRUE(pooled_changed))

  no_species <- data.table(
    read_id = c("p", "n"), donor = "D1", titration_level = "c1",
    species = c("spA", NA_character_), y = c(1L, 0L), score = c(0.9, 0.8),
    arm = "combined", model = "glm", scheme = "LOEO", fold = 1L
  )
  no_species_sxt <- suppressMessages(aggregate_sample_taxon(no_species))
  record_test("F06_NA_SPECIES_NEGATIVE_DROPPED",
              any(no_species$y == 0L) && !any(no_species_sxt$y == 0L))
}

stage07_tests <- function() {
  old_cfg <- deep_copy(cfg)
  on.exit(set_global_cfg(old_cfg), add = TRUE)

  td <- tempfile("stage07_")
  out <- file.path(td, "results")
  dir.create(out, recursive = TRUE)
  tc <- deep_copy(old_cfg)
  tc$paths$out_root <- out
  tc$paths$metrics_read <- file.path(out, "metrics_read_level.tsv")
  tc$paths$metrics_sxt <- file.path(out, "metrics_sample_taxon.tsv")
  tc$paths$hypotheses_out <- file.path(out, "hypothesis_tests.tsv")
  tc$params$n_boot <- 100L
  tc$params$recall_primary <- 0.95
  set_global_cfg(tc)

  paired <- suppressWarnings(run_paired(
    c(`1` = 0.9, `2` = 0.8, `3` = 0.7),
    c(`3` = 0.6, `1` = 0.8, `2` = 0.7),
    "pair", "test", "fold pairing"
  ))
  record_test("S07_FOLD_PAIRING", paired$median_diff)

  inner_scores <- data.table(
    arm = "combined", model = rep(c("glm", "xgboost", "ranger_rf", "glmmTMB"), each = 4),
    fold = rep(1:4, 4), inner_auprc = rep(c(0.90, 0.80, 0.70, 0.60), each = 4)
  )
  fwrite(inner_scores, file.path(out, "inner_cv_scores.tsv"), sep = "\t")
  selection_metrics <- rbind(
    data.table(arm = "combined", model = "glm", fold = 1:4, auprc = 0.4),
    data.table(arm = "combined", model = "xgboost", fold = 1:4, auprc = 0.95)
  )[, `:=`(level = "read", stratum = "all", scheme = "LOEO",
            truncated_included = TRUE, recall_target = 0.95)]
  record_test("S07_INNER_MODEL_SELECTION",
              select_primary_model(selection_metrics, fams = c("glm", "xgboost")))

  h3_rows <- rbindlist(lapply(1:6, function(fold) rbindlist(lapply(
    c("c1", "c2", "c4", "c5"), function(level) data.table(
      arm = "combined", model = c("fixed_threshold", "glm"), fold = fold,
      auprc = c(0.50, if (level %in% c("c4", "c5")) 0.70 else 0.55),
      level = "read", stratum = level, scheme = "LOEO",
      truncated_included = TRUE, recall_target = 0.95
    )
  ))))
  h3 <- suppressWarnings(suppressMessages(test_H3(h3_rows)))
  record_test("S07_H3_LOW_ABUNDANCE_DIRECTION",
              is.finite(h3$median_diff) && h3$median_diff > 0)

  make_metric_fixture <- function(arm, model, base, scheme = "LOEO", level = "read",
                                  stratum = "all", folds = 1:8) {
    data.table(arm = arm, model = model, scheme = scheme, fold = folds,
               level = level, stratum = stratum, truncated_included = TRUE,
               recall_target = 0.95, auprc = base + folds * 0.001,
               prec_at_recall = base, fdr_at_recall = 1 - base,
               n_pos = 10L, n_neg = 10L)
  }
  read_parts <- list(
    make_metric_fixture("combined", "fixed_threshold", 0.50),
    make_metric_fixture("combined", "glm", 0.70),
    make_metric_fixture("combined", "ranger_rf", 0.72),
    make_metric_fixture("combined", "xgboost", 0.75),
    make_metric_fixture("combined", "glmmTMB", 0.73),
    make_metric_fixture("combined", "glmmTMB_none", 0.71),
    make_metric_fixture("combined", "glmmTMB_truth", 0.74),
    make_metric_fixture("combined", "glmmTMB_classifier", 0.73),
    make_metric_fixture("blast_only", "glm", 0.66),
    make_metric_fixture("kraken2_only", "glm", 0.64),
    make_metric_fixture("combined_minus_H5_key", "xgboost", 0.70),
    make_metric_fixture("combined", "xgboost", 0.60, scheme = "LOTO")
  )
  for (level_name in c("c1", "c2", "c4", "c5")) {
    read_parts[[length(read_parts) + 1L]] <- make_metric_fixture(
      "combined", "fixed_threshold", 0.50, stratum = level_name)
    read_parts[[length(read_parts) + 1L]] <- make_metric_fixture(
      "combined", "glm", if (level_name %in% c("c4", "c5")) 0.70 else 0.55,
      stratum = level_name)
    read_parts[[length(read_parts) + 1L]] <- make_metric_fixture(
      "combined", "xgboost", if (level_name %in% c("c4", "c5")) 0.73 else 0.60,
      stratum = level_name)
    read_parts[[length(read_parts) + 1L]] <- make_metric_fixture(
      "blast_only", "glm", 0.52, stratum = level_name)
  }
  metrics_read <- rbindlist(read_parts, fill = TRUE)
  metrics_sxt <- rbindlist(list(
    make_metric_fixture("combined", "glm", 0.80, level = "sample_taxon"),
    make_metric_fixture("combined", "sxt_full", 0.85, level = "sample_taxon"),
    make_metric_fixture("combined", "sxt_minus_breadth", 0.78, level = "sample_taxon")
  ), fill = TRUE)
  fwrite(metrics_read, cfg$paths$metrics_read, sep = "\t")
  fwrite(metrics_sxt, cfg$paths$metrics_sxt, sep = "\t")

  predictions <- rbindlist(lapply(1:8, function(fold) data.table(
    read_id = paste0("f", fold, "_", 1:4), donor = paste0("D", fold),
    titration_level = "c1", species = c("sp1", NA, NA, "sp1"),
    y = c(1L, 0L, 0L, 1L), end_reason_unblock = c(0L, 0L, 1L, 1L),
    arm = "combined", model = "xgboost", scheme = "LOEO", fold = fold,
    score = c(0.9, 0.2, 0.8, 0.85)
  )))
  fwrite(predictions, file.path(out, "predictions.tsv.gz"), sep = "\t")
  calibration <- rbind(
    data.table(arm = "combined", model = "glm", scheme = "LOEO", fold = 1:8, brier = 0.10, n = 20L),
    data.table(arm = "combined", model = "fixed_threshold", scheme = "LOEO", fold = 1:8, brier = 0.20, n = 20L)
  )
  fwrite(calibration, file.path(out, "calibration.tsv"), sep = "\t")
  fwrite(inner_scores, file.path(out, "inner_cv_scores.tsv"), sep = "\t")

  hypothesis_output <- suppressWarnings(suppressMessages(run_stage07()))
  expected_ids <- PARAMS$expectations$S07_REQUIRED_ROWS
  id_counts <- table(hypothesis_output$id)
  ids_ok <- setequal(names(id_counts), expected_ids) && all(id_counts == 1L)
  record_test("S07_REQUIRED_ROWS", if (ids_ok) expected_ids else sort(hypothesis_output$id))
  multiplicity_ok <- all(is.na(hypothesis_output[family == "secondary", holm_sidak_p])) &&
    all(hypothesis_output[!is.na(holm_sidak_p), id] %in% paste0("H", 1:6)) &&
    all(paste0("H", 1:6) %in% hypothesis_output$id)
  record_test("S07_MULTIPLICITY_SCOPE", multiplicity_ok)

  unlink(file.path(out, "inner_cv_scores.tsv"))
  fallback_model <- select_primary_model(selection_metrics, fams = c("glm", "xgboost"))
  record_test("F07_OUTER_SELECTION_FALLBACK", identical(fallback_model, "xgboost"))

  fwrite(inner_scores, file.path(out, "inner_cv_scores.tsv"), sep = "\t")
  arm_fixture <- function(blast_base, kraken_base) rbindlist(list(
    make_metric_fixture("combined", "glm", 0.80),
    make_metric_fixture("blast_only", "glm", blast_base),
    make_metric_fixture("kraken2_only", "glm", kraken_base)
  ))
  comparison_a <- suppressWarnings(test_H1(arm_fixture(0.70, 0.60)))$comparison
  comparison_b <- suppressWarnings(test_H1(arm_fixture(0.50, 0.75)))$comparison
  record_test("F07_SINGLE_ARM_OUTER_SELECTION",
              grepl("blast_only", comparison_a) && grepl("kraken2_only", comparison_b))
}

run_stage_group("01", stage01_tests)
run_stage_group("02", stage02_tests)
run_stage_group("03", stage03_tests)
run_stage_group("04", stage04_tests)
run_stage_group("05", stage05_tests)
run_stage_group("06", stage06_tests)
run_stage_group("07", stage07_tests)

unrecorded <- setdiff(REGISTRY$id, TEST_RESULTS$id)
for (id in unrecorded) record_error(id, simpleError("test was not executed"))
TEST_RESULTS <- TEST_RESULTS[match(REGISTRY$id, id)]

dir.create(dirname(PARAMS$results_file), recursive = TRUE, showWarnings = FALSE)
fwrite(TEST_RESULTS, PARAMS$results_file, sep = "\t")

markdown_escape <- function(x, max_chars = 180L) {
  x <- gsub("[\r\n]+", " ", as.character(x))
  x <- gsub("|", "\\|", x, fixed = TRUE)
  ifelse(nchar(x) > max_chars, paste0(substr(x, 1L, max_chars - 3L), "..."), x)
}

correctness_failures <- TEST_RESULTS[kind == "correctness" & status != "PASS"]
overall <- if (nrow(correctness_failures)) "FAIL" else "PASS"
stage_summary <- TEST_RESULTS[, .(
  correctness_pass = sum(kind == "correctness" & status == "PASS"),
  correctness_fail = sum(kind == "correctness" & status %in% c("FAIL", "ERROR")),
  findings_detected = sum(kind == "finding" & status == "DETECTED"),
  findings_clear = sum(kind == "finding" & status == "CLEAR"),
  errors = sum(status == "ERROR")
), by = stage]

summary_lines <- c(
  "# Pipeline Stage Test Summary",
  "",
  sprintf("- **Run:** %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  sprintf("- **Specification:** `%s`", PARAMS$specification_version),
  sprintf("- **Overall correctness:** **%s**", overall),
  sprintf("- **Correctness:** %d passed, %d failed/error",
          TEST_RESULTS[kind == "correctness" & status == "PASS", .N], nrow(correctness_failures)),
  sprintf("- **Pre-registered findings:** %d detected, %d clear, %d error",
          TEST_RESULTS[kind == "finding" & status == "DETECTED", .N],
          TEST_RESULTS[kind == "finding" & status == "CLEAR", .N],
          TEST_RESULTS[kind == "finding" & status == "ERROR", .N]),
  "- **Command:** `Rscript scripts/tests/test_stage_contracts.R`",
  "- **Production data:** not read or modified; all fixtures were generated under `tempdir()`.",
  "",
  "Correctness checks are blocking. Findings are deliberately non-blocking probes of current",
  "behavior; `DETECTED` means the pre-registered weakness was reproduced, not that the test run failed.",
  "",
  "## Stage Summary",
  "",
  "| Stage | Correctness pass | Correctness fail/error | Findings detected | Findings clear | Errors |",
  "|---:|---:|---:|---:|---:|---:|"
)
summary_lines <- c(summary_lines, apply(stage_summary, 1, function(row) sprintf(
  "| %s | %s | %s | %s | %s | %s |", row[["stage"]], row[["correctness_pass"]],
  row[["correctness_fail"]], row[["findings_detected"]], row[["findings_clear"]], row[["errors"]]
)))
summary_lines <- c(summary_lines,
  "", "## Detailed Results", "",
  "| ID | Stage | Type | Expected | Observed | Status | Test |",
  "|---|---:|---|---|---|---|---|"
)
summary_lines <- c(summary_lines, vapply(seq_len(nrow(TEST_RESULTS)), function(i) {
  row <- TEST_RESULTS[i]
  sprintf("| `%s` | %s | %s | `%s` | `%s` | **%s** | %s |",
          row$id, row$stage, row$kind, markdown_escape(row$expected),
          markdown_escape(row$actual), row$status, markdown_escape(row$description))
}, character(1)))

detected <- TEST_RESULTS[kind == "finding" & status == "DETECTED"]
summary_lines <- c(summary_lines, "", "## Detected Findings", "")
if (nrow(detected)) {
  summary_lines <- c(summary_lines, vapply(seq_len(nrow(detected)), function(i) sprintf(
    "- `%s` (stage %s): %s", detected$id[i], detected$stage[i], detected$description[i]
  ), character(1)))
} else {
  summary_lines <- c(summary_lines, "No pre-registered finding was detected.")
}
summary_lines <- c(summary_lines, "", "## Files", "",
  "- Expected outcomes: `scripts/tests/stage_test_parameters.R`",
  "- Test runner: `scripts/tests/test_stage_contracts.R`",
  "- Machine-readable results: `reports/pipeline_stage_test_results.tsv`"
)
writeLines(summary_lines, PARAMS$summary_file)

cat(sprintf("\nSummary: %s\n", PARAMS$summary_file))
cat(sprintf("Details: %s\n", PARAMS$results_file))
cat(sprintf("Correctness: %d passed, %d failed/error; findings: %d detected, %d clear\n",
            TEST_RESULTS[kind == "correctness" & status == "PASS", .N],
            nrow(correctness_failures),
            TEST_RESULTS[kind == "finding" & status == "DETECTED", .N],
            TEST_RESULTS[kind == "finding" & status == "CLEAR", .N]))

quit(status = if (nrow(correctness_failures)) 1L else 0L)