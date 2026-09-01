## =============================================================================
## 00_config.R  --  Zymo-in-human read-classification ML pipeline
## -----------------------------------------------------------------------------
## Single source of truth. EVERY other script sources this file first.
##
## Review two things before running:
##   SECTION 1 -- HARDCODED FILE PATHS   (edit to match your system)
##   SECTION 5 -- OPEN ITEMS             (decisions that still block a full run)
##
## Language note: the pipeline is implemented in R + shelled-out bioinformatics
## tools (minimap2 / BLASTn / Kraken2). R was chosen because the statistical
## layer (glmmTMB GLMM, exact Wilcoxon, Holm-Sidak, bootstrap CIs) is native
## here. See OPEN ITEM T2 if you would rather drive this from Python/Snakemake.
## =============================================================================

cfg <- list()

## Null-coalesce, used across stages. Defined here so it exists after sourcing config.
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

## =============================================================================
## SECTION 1 -- HARDCODED FILE PATHS        <<< EDIT THESE >>>
## -----------------------------------------------------------------------------
## Any value that is NA_character_ or contains "<FILL_IN>" is a KNOWN GAP and is
## cross-referenced by id in SECTION 5 (OPEN ITEMS).
## =============================================================================

PROJECT_ROOT <- "/home/tdinse/master_timo/projektarbeit/projektarbeit_ml"

cfg$paths <- list(
  project_root = PROJECT_ROOT,

  ## --- inputs: raw sequencing ------------------------------------------------
  ## One basecalled fastq[.gz] per library; 8 donors x 6 barcodes = 48 libraries.
  reads_dir       = file.path(PROJECT_ROOT, "data", "reads"),
  ## Sample sheet, one row per library. Required columns documented in SECTION 2.
  sample_sheet    = file.path(PROJECT_ROOT, "data", "sample_sheet.tsv"),
  ## ONT sequencing_summary.txt per run (end_reason, channel, start_time).  [OI 6]
  seq_summary_dir = file.path(PROJECT_ROOT, "data", "sequencing_summary"),

  ## --- inputs: references & databases ----------------------------------------
  ## Zymo per-species genomes ship as a DIRECTORY (D6331 bundle). zymo_members.fasta
  ## is BUILT from it (each species collapsed to one species-named record) by
  ## scripts/prepare_zymo_genomes.R -- run that once before stage 01. [OI 1]
  zymo_refs_dir   = "/home/tdinse/Downloads/D6331.refseq/genomes",          # 21 member genomes  [OI 1]
  zymo_ssrrna_dir = "/home/tdinse/Downloads/D6331.refseq/ssrRNAs",          # 16S/18S -> ambiguous BED [OI 2]
  zymo_refs_fasta = file.path(PROJECT_ROOT, "ref", "zymo_members.fasta"),    # BUILT from zymo_refs_dir [OI 1]
  human_grch38    = "/home/tdinse/new_hg38/hg38.fa",                        # live depletion + GT human check [OI 5]
  human_t2t_chm13 = "/home/tdinse/new_T2T/T2T.fasta",                       # human-competitor score  [note F]
  ## Kraken2: a CURATED, species-level DB (one reference per species, incl. fungi/
  ## protists) -- NOT Standard (no fungi -> misses Zymo Candida/Saccharomyces) and
  ## NOT core_nt (nt redundancy degrades LCA + per-species k-mer features). The NCBI
  ## reference DB below (k2_NCBI_reference_20251007) satisfies this; it ships as a
  ## .tar.gz and is extracted in place -> the directory holds hash.k2d/opts.k2d/
  ## taxo.k2d + taxonomy. [OI 9 RESOLVED 2026-08-04]
  kraken2_db      = file.path(PROJECT_ROOT, "database_kraken2", "k2_NCBI_reference_20251007"),  # curated species-level DB [OI 9]
  blast_db        = "/home/tdinse/core_nt/core_nt",                         # multi-volume prefix, taxdb present [OI 9]

  ## --- inputs: wet-lab metadata ----------------------------------------------
  ## Lot-specific certificate of analysis: per-strain relative abundance.  [note H / OI 3]
  zymo_coa        = file.path(PROJECT_ROOT, "data", "zymo_coa_lot_LOT_ID.tsv"),

  ## --- outputs ---------------------------------------------------------------
  out_root    = file.path(PROJECT_ROOT, "results"),
  work_dir    = file.path(PROJECT_ROOT, "work"),     # per-read intermediate tables
  reports_dir = file.path(PROJECT_ROOT, "reports")
)

## Base results root. The ACTUAL per-run output dir (cfg$paths$out_root) and every
## cutoff-DEPENDENT derived path (labels / features / cv / models / metrics /
## hypotheses) are finalised in SECTION 3b, once the ground-truth cutoff PROFILE is
## known, so each run lands in its OWN results/<gt_run_tag>/ folder. Expensive,
## cutoff-INDEPENDENT stage-01 alignments + the cutoff-sensitivity tables stay in
## the shared work/ dir (per-library subdirs + sample_taxon_coverage / expected /
## cutoff_* tables) and are reused across both runs.
cfg$paths$results_base <- cfg$paths$out_root

## =============================================================================
## SECTION 2 -- SAMPLE SHEET SCHEMA
## -----------------------------------------------------------------------------
## cfg$paths$sample_sheet is a TSV with (at minimum) these columns:
##   library_id        unique id, e.g. D01_c1
##   donor             one of cfg$params$donors  (the LOEO fold unit)  [guardrail]
##   barcode           ONT barcode id within the run
##   titration_level   one of {negative, c1..c5}
##   concentration     input cells for the library (= cells_per_level[level]); metadata
##   run_id            flow-cell / run id  -- needed to check the run<->donor
##                     confounder and to model cross-barcode leakage  [note L / OI 8]
##   fastq             absolute path to the barcode's reads (fastq[.gz])
##   seq_summary       FOLDER holding that run's ONT reports; stage 01 picks the
##                     sequencing_summary_*<run_id>*.txt inside it        [note E / OI 6]
## =============================================================================

## =============================================================================
## SECTION 3 -- PARAMETERS
## =============================================================================
cfg$params <- list(
  seed = 1729L,

  ## stage-05 model-training threads (ranger + xgboost). Resolve from TRAIN_THREADS,
  ## else PIPE_THREADS, else all detected cores. Affects SPEED ONLY, not the models
  ## (thread count does not bias RF/GBM results). ranger/xgboost otherwise defaulted
  ## to 2 threads here -> stage 05 ran ~2 cores wide on a 64-core box. [perf 2026-08-12]
  train_threads    = local({
    n <- suppressWarnings(as.integer(Sys.getenv("TRAIN_THREADS", Sys.getenv("PIPE_THREADS", ""))))
    if (length(n) != 1L || is.na(n) || n < 1L) max(1L, parallel::detectCores()) else n
  }),

  ## sequencing provenance [R1 / OI 7] -- VERIFIED identical across all 8 runs
  ## (2026-08-03): every run's basecalled fastq carries RG:Z =
  ## <run_id>_dna_r10.4.1_e8.2_400bps_hac@v6.0.0_SQK-NBD114-96, and every demux
  ## path is 'exp_c3_dorado200_demux' -> Dorado 2.0.0 for all 8 flow cells.
  ## HAC (not SUP) simplex -> single-read accuracy & Q-scores reflect the hac model.
  basecaller       = "dorado 2.0.0",
  basecaller_model = "dna_r10.4.1_e8.2_400bps_hac@v6.0.0",
  library_kit      = "SQK-NBD114-96",

  ## --- experimental design ---------------------------------------------------
  ## 8 LOEO folds, derived from the sample sheet so the set always matches the data
  ## (this cohort is D01,D03-D09; D02 absent, D09 present). [guardrail]
  donors           = local({
    ss <- cfg$paths$sample_sheet
    d  <- tryCatch(sort(unique(read.delim(ss, sep = "\t", stringsAsFactors = FALSE)$donor)),
                   error = function(e) NULL)
    if (length(d)) d else sprintf("D%02d", 1:8)
  }),
  titration_levels = c("c1", "c2", "c3", "c4", "c5"),# 'negative' handled separately
  low_mass_levels  = c("c4", "c5"),                  # lowest-input levels, most inflated by cross-barcode leakage [note L]
  barcodes_per_run = 6L,                             # negative + c1..c5

  ## --- ground-truth calling via minimap2 [note A] ----------------------------
  gt_min_identity  = 0.90,   # read must align to a Zymo genome above this identity
  gt_min_coverage  = 0.80,   # fraction of read length that must align
  gt_human_margin  = 0.0,    # positive requires zymo_score > human_score + margin (bits)
  ## offline human depletion before classification (R2), by the SAME human-vs-Zymo
  ## COMPETITION the labeller uses: a read is dropped ('human_like') iff it aligns
  ## BETTER to human than to any Zymo member -- human_score = max(GRCh38, T2T)
  ## residue matches > best Zymo score. NO identity/coverage gate: those gates reject
  ## genuine SHORT human reads (minimap2 map-ont is a long-read preset) and let ~83%
  ## human LEAK into the classifier, whereas the competition removes ~83% here yet
  ## keeps every read a Zymo member wins (only ~0.1% align to both). The identical
  ## human_score is used in stage-01 depletion, stage-02 labelling and
  ## prepare_cutoff_sensitivity.R, so depletion and the positive/negative labels
  ## cannot disagree. [note D/E / R2]
  ## reads hitting rRNA operons / mobile elements / plasmid backbones / low
  ## complexity are routed to the 'ambiguous' class, NOT forced into binary. [note A / OI 2]
  ambiguous_bed    = file.path(PROJECT_ROOT, "ref", "zymo_ambiguous_regions.bed"), # [OI 2]
  ## only region-mask low-complexity blocks >= this many bp (note A); shorter
  ## micro-stretches are captured per-read by the DUST feature and would otherwise
  ## over-flag long reads via any-overlap.
  lowcomplexity_min_len = 300L,

  ## --- ground-truth cutoff SENSITIVITY [F9] ----------------------------------
  ## gt_min_identity / gt_min_coverage above are STRESS-TESTED, not assumed:
  ## scripts/prepare_cutoff_sensitivity.R builds the empirical identity/coverage
  ## distribution of CONFIDENT Zymo reads (Zymo hit beats human) at the high-titre
  ## level(s) below, sweeps a grid of cutoffs, and reports the fraction of that
  ## confident set retained -- so key results can be shown stable across a
  ## defensible cutoff range. Run AFTER stage 01 (needs the gt_zymo/gt_human PAFs).
  gt_calib_levels  = c("c1"),                              # high-titre libs for the empirical distribution
  gt_identity_grid = c(0.80, 0.85, 0.90, 0.92, 0.95, 0.97),
  gt_coverage_grid = c(0.50, 0.60, 0.70, 0.80, 0.90),
  gt_retain_frac   = 0.95,   # recommend the strictest cutoffs still retaining >= this frac of confident-Zymo reads

  ## --- metrics [note C / R9] -------------------------------------------------
  ## Clinical pathogen detection wants HIGH sensitivity (a missed pathogen is worse
  ## than a false positive that triggers a confirmatory test). We report precision +
  ## FDR at several high-recall operating points and run the hypothesis tests at the
  ## primary one. (Supersedes the project_plan 0.80 placeholder.)
  ## Operating point CONFIRMED 2026-08-03 [OI 14]: primary recall = 0.95.
  recall_target    = c(0.90, 0.95, 0.99),
  recall_primary   = 0.95,        # operating point used for the hypothesis tests
  ## FDR is REPORTED at each recall target (stage 06); no single hard acceptance
  ## FDR was mandated by the clinical group. Record the agreed threshold here if
  ## one is set later (e.g. 0.05); NA = report-only. [OI 14]
  target_fdr       = NA_real_,
  strata           = "titration_level",
  report_truncated = c(TRUE, FALSE),  # report with and without unblocked reads [note E]

  ## --- feature engineering [note F] ------------------------------------------
  top_n_hits       = 10L,   # margin / entropy features computed over top-N BLAST hits

  ## --- Poisson expected-copy-number floor [note H] ---------------------------
  poisson_p_min    = 0.95,  # species x level with P(>=1 genome) below this -> indeterminate
  ## total input cells per titration level (c1 worked example = 7.88e5 in note H).
  cells_per_level  = c(c1 = 78800000, c2 = 7880000, c3 = 788000, c4 = 78800, c5 = 7880),  # [OI 3]
  ## Unmeasured -> set to defensible rough values (C3). Only the PRODUCT matters
  ## (combined recovery ~0.5): community-average metagenomic DNA extraction ~50%,
  ## prepared library loaded ~in full. Kept sub-unity so we do not over-claim which
  ## lowest-mass taxa were truly sequenceable; sweep poisson_p_min for sensitivity. [OI 3]
  extraction_eff   = 0.5,
  fraction_loaded  = 1.0,

  ## --- modelling -------------------------------------------------------------
  ## feed-forward NN deliberately excluded from the comparison. [note K]
  models           = c("fixed_threshold", "glm", "ranger_rf", "xgboost", "glmmTMB"),
  inner_folds      = 5L,          # nested-CV inner loop for hyperparameter tuning
  ## Read counts vary a lot across donors/barcodes (10x dilution series), so a
  ## deeply-sequenced library must not dominate the pooled training fit. Weight
  ## each read so every unit contributes equally: "none" | "library" | "donor".
  ## (Sample x taxon COUNT features are additionally depth-normalised per-million
  ## in stage 06; the ranking metrics are per-fold so already depth-robust.)
  balance_unit     = "donor",
  ## NO downsampling of negatives (R10 / overall decision): at pathogen-relevant
  ## abundance the true signal can be ~1-10 reads per MILLION, so dropping reads
  ## would discard the key differentiators. Keep the full read set and fit the GLMM
  ## on all data (glmmTMB; segment by chunk if memory-bound). NA = keep all. [note I]
  negative_downsample_ratio = NA_real_,

  ## GLMM training-row cap (stage-05 SPEED lever) [perf 2026-08-15]. glmmTMB is
  ## single-threaded TMB; on the combined arm (~24 fixed effects + (1|donor)+
  ## (1|species)) a full ~4M-row fit takes ~6-7 h, and the GLMM is fit ~48x in
  ## stage 05 (8 LOEO folds + the 24-fit H9 decomposition) -> a multi-DAY tail that
  ## dominates the run. A ~24-parameter GLMM is estimated to the SAME accuracy from a
  ## few 100k rows as from millions (measured 2026-08-15: combined-arm AUPRC 0.999
  ## @250k vs full 4M; fit 6.9 min vs hours). So cap the GLMM's TRAINING rows by a
  ## label-STRATIFIED subsample at the NATURAL class ratio -- keeping the same
  ## fraction of each class preserves prevalence, so the fitted intercept stays
  ## calibrated to the test set and NO prior offset is needed. This is the GLMM ONLY:
  ## ranger/xgboost/glm still train on the FULL data (they are fast once threaded and
  ## must see every negative to learn the false-positive background), and prediction
  ## is always on the FULL test set. NA = no cap (fit on all rows). [note I]
  glmm_max_train_rows = 250000L,

  ## --- H9: GLMM decomposition + species-source sensitivity [H9] ---------------
  ## H9 is reported as TWO transparent terms: (1|donor) ALONE (does modelling donor
  ## variance in TRAINING aid transfer to a held-out donor?) and (1|donor)+(1|species)
  ## (adds a species baseline). Their difference = the species term's contribution.
  ## The species grouping is a SENSITIVITY axis, run three ways so H9's conclusion can
  ## be shown robust (or not) to WHERE 'species' comes from:
  ##   none        no species term            -> the (1|donor)-only variant
  ##   truth       minimap2 ground-truth species (labels `species`)
  ##   classifier  the aligner's COARSE-rank call (BLAST top-hit genus, `top_genus`)
  ## A source whose column is absent or single-level in a fold is skipped there.
  h9_arm             = "combined",
  h9_species_sources = c(none = NA, truth = "species", classifier = "top_genus"),

  ## --- statistics [note J + project_plan] ------------------------------------
  alpha          = 0.05,
  primary_family = c("H1", "H2", "H3", "H4", "H5", "H6"),  # <=6 -> Holm-Sidak feasible at n=8 [note J]
  n_boot         = 2000L,         # donor-level bootstrap CIs for median differences
  wilcox_exact   = TRUE           # exact Wilcoxon signed-rank across the 8 folds
)

## =============================================================================
## SECTION 3b -- GROUND-TRUTH CUTOFF PROFILE  (two reproducible runs)
## -----------------------------------------------------------------------------
## The supervised target hinges on the ground-truth identity/coverage cutoffs, so
## the pipeline is designed to be run TWICE, once per cutoff PROFILE, selected by
## the GT_PROFILE environment variable (so the choice propagates to run_pipeline.R
## AND to any stage run standalone):
##   GT_PROFILE=fixed       gt_min_identity = 0.90, gt_min_coverage = 0.80 (the
##                          values set in SECTION 3 above) -- the first run.
##   GT_PROFILE=calculated  the cutoffs RECOMMENDED by prepare_cutoff_sensitivity.R
##                          (work/cutoff_recommended.tsv): the strictest (identity,
##                          coverage) pair still retaining >= gt_retain_frac of the
##                          confident-Zymo reads -- the second run.
## Each profile writes into its OWN results/<gt_run_tag>/ folder (the tag encodes
## the profile + the active cutoffs), so the two runs never overwrite each other.
## Run BOTH in one shot with:  Rscript scripts/run_pipeline.R --gt both
## =============================================================================

## Resolve the active (identity, coverage) cutoffs for a profile. Pure + testable:
## given the fixed defaults and the path to the recommended-cutoffs file, it returns
## the pair to use, or stops with actionable guidance when 'calculated' is requested
## before prepare_cutoff_sensitivity.R has produced its recommendation.
cfg_gt_profile_cutoffs <- function(profile, fixed_identity, fixed_coverage, rec_path) {
  profile <- tolower(profile)
  if (identical(profile, "fixed"))
    return(list(gt_min_identity = fixed_identity, gt_min_coverage = fixed_coverage))
  if (identical(profile, "calculated")) {
    if (!file.exists(rec_path))
      stop("GT_PROFILE=calculated needs recommended cutoffs, but\n  ", rec_path,
           "\nis missing. Run stage 01, then:  Rscript scripts/prepare_cutoff_sensitivity.R\n",
           "(or run:  Rscript scripts/run_pipeline.R --gt both, which builds it automatically).",
           call. = FALSE)
    rec <- utils::read.delim(rec_path, sep = "\t", stringsAsFactors = FALSE)
    if (!nrow(rec) || !all(c("gt_min_identity", "gt_min_coverage") %in% names(rec)))
      stop("Recommended-cutoffs file is malformed (need columns gt_min_identity + ",
           "gt_min_coverage):\n  ", rec_path, call. = FALSE)
    return(list(gt_min_identity = as.numeric(rec$gt_min_identity[1]),
                gt_min_coverage = as.numeric(rec$gt_min_coverage[1])))
  }
  stop("Unknown GT_PROFILE='", profile, "'. Use 'fixed' or 'calculated'.", call. = FALSE)
}

## Filesystem-safe run tag encoding the profile + the ACTIVE cutoffs.
cfg_gt_run_tag <- function(profile, identity, coverage)
  sprintf("gt_%s_id%.2f_cov%.2f", tolower(profile), identity, coverage)

cfg$params$gt_profile        <- tolower(Sys.getenv("GT_PROFILE", "fixed"))
cfg$params$gt_fixed_identity <- cfg$params$gt_min_identity   # 0.90, before any override
cfg$params$gt_fixed_coverage <- cfg$params$gt_min_coverage   # 0.80, before any override
## The recommended-cutoffs file is profile-INDEPENDENT (built from the shared work/
## PAFs by prepare_cutoff_sensitivity.R), so it lives in work/, not the run folder.
cfg$paths$gt_recommended     <- file.path(cfg$paths$work_dir, "cutoff_recommended.tsv")

.gtc <- cfg_gt_profile_cutoffs(cfg$params$gt_profile,
                               cfg$params$gt_fixed_identity, cfg$params$gt_fixed_coverage,
                               cfg$paths$gt_recommended)
cfg$params$gt_min_identity <- .gtc$gt_min_identity
cfg$params$gt_min_coverage <- .gtc$gt_min_coverage
cfg$params$gt_run_tag      <- cfg_gt_run_tag(cfg$params$gt_profile,
                                             cfg$params$gt_min_identity, cfg$params$gt_min_coverage)
rm(.gtc)

## Per-run output root + every cutoff-DEPENDENT derived path beneath it.
cfg$paths$out_root       <- file.path(cfg$paths$results_base, cfg$params$gt_run_tag)
cfg$paths$labels_table   <- file.path(cfg$paths$out_root, "labels.tsv")
cfg$paths$feature_table  <- file.path(cfg$paths$out_root, "feature_table.parquet")
cfg$paths$cv_splits      <- file.path(cfg$paths$out_root, "cv_splits.rds")
cfg$paths$leakage_table  <- file.path(cfg$paths$out_root, "leakage_estimates.tsv")
cfg$paths$model_dir      <- file.path(cfg$paths$out_root, "models")
cfg$paths$metrics_read   <- file.path(cfg$paths$out_root, "metrics_read_level.tsv")
cfg$paths$metrics_sxt    <- file.path(cfg$paths$out_root, "metrics_sample_taxon.tsv")
cfg$paths$hypotheses_out <- file.path(cfg$paths$out_root, "hypothesis_tests.tsv")

## =============================================================================
## SECTION 4 -- FEATURE BLOCKS, ARMS & HYPOTHESIS REGISTRY
## =============================================================================

## Taxon-agnostic feature blocks. Used for the H5 ablation and for building the
## classifier arms below. Taxon IDENTITY (species/genus name, taxid) is never a
## feature. [note B]
cfg$feature_blocks <- list(
  blast_core       = c("qlen", "pident", "evalue", "mismatch", "aln_fraction", "query_cov", "bitscore"),
  blast_margin     = c("bitscore_margin_species", "n_species_topN", "n_genera_topN", "tax_entropy_topN"),
  human_competitor = c("human_bitscore", "human_pident", "human_minus_best_margin"),
  kraken2          = c("k2_conf", "k2_kmers_taxon_frac", "k2_distinct_minimizers"),
  complexity       = c("dust_score", "gc", "homopolymer_frac"),
  read_qc          = c("read_len", "mean_q", "end_reason_unblock"),
  ## sample x taxon aggregates -- built in stage 06, co-primary endpoint [note G]
  sample_taxon     = c("genome_breadth", "coverage_evenness", "n_reads_above_thr", "mean_score", "max_score")
)

## Only taxon-agnostic properties of the DB subject are admissible. [note B]
## [M1] Only subject_genome_len is actually produced by parse_blast (BLAST outfmt 6
## = slen). Assembly-level / WGS-draft status would need NCBI assembly metadata that
## is not wired here, so they are intentionally NOT declared (declared == produced).
cfg$subject_props <- c("subject_genome_len")

## The two classifiers are scored INDEPENDENTLY, never sequentially. [note D]
cfg$classifier_arms <- list(
  blast_only   = c("blast_core", "blast_margin", "human_competitor", "complexity", "read_qc"),
  kraken2_only = c("kraken2", "complexity", "read_qc"),
  combined     = c("blast_core", "blast_margin", "human_competitor", "kraken2", "complexity", "read_qc")
)

## H5 ablation: drop each of these feature-block sets in turn from the 'combined'
## arm. `H5_key` is the pre-registered union (alignment-margin + human-competitor)
## used as the single primary H5 test; the single-block sets are for interpretation.
## Genome-breadth lives at the sample x taxon level; it is ablated THERE by
## sample_taxon_ablation() via cfg$sxt_ablation_sets, not here. [note F / H5 / F6]
cfg$ablation_sets <- list(
  H5_key           = c("blast_margin", "human_competitor"),
  blast_margin     = "blast_margin",
  human_competitor = "human_competitor"
)

## [F6] Sample x taxon-level ablation for the H5 'breadth' component. genome_breadth
## and coverage_evenness exist ONLY at the sample x taxon endpoint, so they cannot be
## dropped from the read-level 'combined' arm (stage 05 correctly skips any ablation
## that removes no read-level feature). Instead stage 06 sample_taxon_ablation() fits
## the sample x taxon aggregate score WITH vs WITHOUT these blocks under leave-one-
## donor-out CV, so H5's 'breadth' claim is MEASURED, not asserted. The base sample x
## taxon score uses the scale-free read aggregates below. [F6 / note F,G]
cfg$sxt_score_features <- c("mean_score", "max_score", "reads_above_thr_per_million")
cfg$sxt_ablation_sets  <- list(
  breadth = c("genome_breadth", "coverage_evenness")
)

## Hypothesis registry. `family` = primary (alpha-spent) vs secondary (CIs only).
## `active` = FALSE means it is documented but not run (see notes K, F).
cfg$hypotheses <- data.frame(
  id = paste0("H", 1:12),
  family = c(rep("primary", 6), rep("secondary", 6)),
  ## H7 (adapter content) DEFERRED to later [2026-08-03]; H8 (NN) dropped [note K].
  active = c(rep(TRUE, 6), FALSE, FALSE, TRUE, TRUE, TRUE, TRUE),
  test = c(
    "H1  combined BLASTn+Kraken2 > best single arm     : paired Wilcoxon, model family fixed",
    "H2  best ML > best fixed-threshold baseline        : paired Wilcoxon, thresholds tuned on train only",
    "H3  ML gain larger at low titration than high      : mixed model metric ~ method*level + (1|donor)",
    "H4  tree ensembles (RF/XGB) > linear GLM           : paired Wilcoxon",
    "H5  margin+human-competitor+breadth drive the gain : paired Wilcoxon on feature-block ablation",
    "H6  sample x taxon aggregation > read thresholding  : paired Wilcoxon over donors",
    "H7  adapter-positive reads enriched in FP vs TP     : DEFERRED -- experiment for later (adapter content); see reporting_2 [OI T1]",
    "H8  feed-forward NN not better than XGBoost         : DROPPED per note K (equivalence not run)",
    "H9  GLMM (1|donor) generalises > GLM; species RE decomposed + source-swept (H9b) : paired Wilcoxon, held-out donors",
    "H10 leave-one-taxon-out generalisation             : quantify degradation vs LOEO",
    "H11 truncated/unblocked reads differ & drive FP     : score-distribution test + FP share",
    "H12 models better calibrated than fixed thresholds  : Brier score + calibration curves"
  ),
  stringsAsFactors = FALSE
)

## =============================================================================
## SECTION 5 -- OPEN ITEMS  (still blocking a full, defensible run)
## -----------------------------------------------------------------------------
## Machine-readable so run_pipeline.R can print the unresolved ones before it
## starts. IDs 1-16 mirror open_items.txt; T1-T5 are tensions found while
## reconciling notes.txt with project_plan.md. Set status = "resolved" (and fill
## the relevant path/param above) as you close each one.
## =============================================================================
cfg$open_items <- data.frame(
  id = c(as.character(1:16), paste0("T", 1:5)),
  ## status updated as items were closed across development; run_pipeline.R --status
  ## prints only those != "resolved". Cross-referenced in open_items_v2.md.
  status = c(
    "resolved", "resolved", "resolved", "resolved",   # 1-4
    "resolved", "resolved", "resolved", "resolved",   # 5-8
    "resolved", "resolved", "resolved", "resolved",   # 9-12
    "resolved", "resolved", "resolved", "resolved",   # 13-16
    "resolved", "resolved", "resolved", "resolved", "resolved"  # T1-T5
  ),
  blocks_stage = c(
    "02_labels", "02_labels", "02_labels/05_train", "02_labels",   # 1-4
    "01_tools",  "01_tools",  "03_features", "04_cv",              # 5-8
    "01_tools/03_features", "02_labels", "01_tools", "01_tools",    # 9-12
    "06_eval",   "06_eval",   "04_cv",     "05_train",              # 13-16
    "03_features/07_tests", "all", "02_labels", "07_tests", "02_labels"  # T1-T5
  ),
  item = c(
    "Exact Zymo member assemblies for THIS lot (incl. fungal/archaeal)? -> set zymo_refs_fasta.",
    "Only MOBILE-element masking deferred (ISEScan, optional & least critical); rRNA operons (barrnap), plasmid backbones, and substantial low-complexity (dustmasker >=300 bp) are now masked; cutoffs accepted (C4). [note A / R8]",
    "Lot-specific certificate of analysis (per-strain rel. abundance) + extraction_eff + fraction_loaded -> Poisson floor.",
    "Non-Zymo bacterial reads in spiked samples: negatives, or indeterminate (donor bacteraemia/kitome)?",
    "GRCh38 vs T2T for LIVE adaptive sampling, and which reference for the offline human-competitor score -> human_* paths.",
    "Is sequencing_summary retained for all 48 libraries (end_reason/channel/start_time)? -> seq_summary_dir.",
    "Basecaller version+model (Dorado hac vs sup, simplex/duplex)? Q-score features are not comparable across models.",
    "All 48 libraries on one flow cell/run or several? Determines whether run confounds donor -> CV nesting.",
    "RESOLVED 2026-08-04: provided k2_NCBI_reference_20251007 (NCBI reference DB, one reference per species incl. fungi/protists), extracted to database_kraken2/ and wired to cfg$paths$kraken2_db. Original guidance: use a CURATED species-level DB -- PrackenDB (one genome/species, ideal for the k-mer-per-species features; covers fungi/protists) or PlusPFP fallback. NOT Standard (no fungi -> misses Zymo Candida/Saccharomyces), NOT core_nt (nt redundancy degrades Kraken2 LCA). BLAST keeps core_nt (H7 contamination); minimap2 ground truth is classifier-independent so H1 compares feature sets as deployed. Optional controlled arm: build Kraken2 from core_nt (prepare_kraken2_from_core_nt.sh). NB a species-level DB collapses E. coli to ~species -> consider collapsing the 5 CoA E. coli strains at the sample x taxon endpoint.",
    "No reagent-only extraction blanks exist; per-donor negative (c0) barcodes are the only no-spike controls (donor+reagent background, not a pure kitome) -> reagent kitome not separately characterisable (documented limitation).",
    "DECIDED 2026-08-03: FULL-RUN analysis only, no subsampling -- BLAST/Kraken2 score every non-human read (human-depleted universe, R2). [OI 11]",
    "DECIDED 2026-08-03: no subsampling; full-run analysis throughout (development + final eval). [OI 12]",
    "Primary endpoint at read level, sample x taxon level, or both? (Affects the title.)",
    "CONFIRMED 2026-08-03: clinical operating point = primary recall 0.95 (report 0.90/0.95/0.99). FDR reported at each target; no single hard acceptance FDR mandated (cfg$params$target_fdr=NA hook). [R9/OI14]",
    "Are the 8 donors the only folds, or also hold out an entire flow cell / a taxon set (LOTO)?",
    "GLMM software stack at scale (glmmTMB / Julia MixedModels) and whether negative downsampling + calibration correction is acceptable.",
    "TENSION: note F says ignore residual barcode/adapter as a feature, but H7 tests adapter content. Plan: adapter features computed ONLY for the H7 exploratory test, excluded from all primary models. Confirm.",
    "TENSION: implementation language/stack is R + bash. Confirm, or request Python / Snakemake / Nextflow.",
    "TENSION: cross-barcode leakage correction (note L) needs per-run demultiplex error rate / barcode-assignment config to model leakage into c4/c5.",
    "TENSION: NN (H8) dropped per note K; registry keeps it inactive for provenance. Confirm it stays out of the alpha-spending family.",
    "TENSION: Poisson floor also needs extraction efficiency and fraction-loaded (see item 3); until then species x level strata cannot be finalised."
  ),
  stringsAsFactors = FALSE
)

## =============================================================================
## SECTION 6 -- HELPERS TO VALIDATE / REPORT CONFIG
## =============================================================================

## Return open items that are not yet resolved.
cfg_open_items <- function(config = cfg) {
  config$open_items[config$open_items$status != "resolved", , drop = FALSE]
}

## Flag hardcoded paths that are still placeholders or missing on disk.
cfg_missing_paths <- function(config = cfg) {
  p <- unlist(config$paths)
  is_placeholder <- grepl("<FILL_IN>", p) | is.na(p)
  ## only check *input* paths for existence; outputs (incl. the BUILT
  ## zymo_refs_fasta) are created by the pipeline. reads_dir is omitted because
  ## the sample sheet carries absolute `fastq` paths (validated per-row in stage 01).
  input_keys <- c("sample_sheet", "zymo_refs_dir",
                  "human_grch38", "human_t2t_chm13", "kraken2_db", "blast_db", "zymo_coa")
  exists_input <- function(k, path) {
    if (k == "blast_db") length(Sys.glob(paste0(path, "*"))) > 0L  # prefix, multi-volume
    else file.exists(path)
  }
  missing_on_disk <- vapply(names(p), function(k) {
    k %in% input_keys && !is_placeholder[[k]] && !exists_input(k, p[[k]])
  }, logical(1))
  flagged <- is_placeholder | missing_on_disk
  data.frame(key = names(p)[flagged], path = unname(p[flagged]),
             reason = ifelse(is_placeholder[flagged], "placeholder", "not found"),
             stringsAsFactors = FALSE)
}

## Create output scaffolding.
cfg_init_dirs <- function(config = cfg) {
  for (d in c(config$paths$out_root, config$paths$work_dir, config$paths$reports_dir,
              config$paths$model_dir)) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
  ## Provenance: record the ground-truth cutoffs that DEFINE this run's folder, so
  ## results/<gt_run_tag>/ is self-documenting even if moved or shared.
  gt <- data.frame(gt_profile      = config$params$gt_profile,
                   gt_min_identity = config$params$gt_min_identity,
                   gt_min_coverage = config$params$gt_min_coverage,
                   run_tag         = config$params$gt_run_tag,
                   stringsAsFactors = FALSE)
  try(utils::write.table(gt, file.path(config$paths$out_root, "ground_truth_settings.tsv"),
                         sep = "\t", row.names = FALSE, quote = FALSE), silent = TRUE)
  invisible(TRUE)
}

set.seed(cfg$params$seed)
