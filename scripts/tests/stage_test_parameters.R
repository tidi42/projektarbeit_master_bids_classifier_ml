## =============================================================================
## stage_test_parameters.R -- pre-registered expectations for stages 01-07
## =============================================================================

if (!exists("TEST_PROJECT_ROOT", inherits = FALSE))
  stop("TEST_PROJECT_ROOT must be defined before sourcing stage_test_parameters.R")

PIPELINE_TEST_PARAMETERS <- list(
  specification_version = "2026-09-01",
  seed = 20260901L,
  numeric_tolerance = 1e-8,
  summary_file = file.path(TEST_PROJECT_ROOT, "reports", "pipeline_stage_test_summary.md"),
  results_file = file.path(TEST_PROJECT_ROOT, "reports", "pipeline_stage_test_results.tsv"),
  expectations = list(
    S01_RAGGED_PAF = list(n_best = 3L, r1_score = 95),
    S01_HUMAN_COMPETITION = c("r_human", "r_t2t"),
    S01_SAMPLE_SHEET_SCHEMA = TRUE,
    S01_COMMAND_WIRING = TRUE,
    S01_EMPTY_OUTPUT_INCOMPLETE = TRUE,
    F01_RESUME_OMITS_PAF_QC = TRUE,

    S02_LABEL_ASSIGNMENT = c(
      r_ambiguous = "ambiguous",
      r_human = "negative",
      r_lowcov = "negative",
      r_lowid = "negative",
      r_none = "negative",
      r_positive = "positive",
      r_tie = "negative"
    ),
    S02_QC_DEFINES_UNIVERSE = c(
      "r_ambiguous", "r_human", "r_lowcov", "r_lowid", "r_none",
      "r_positive", "r_tie"
    ),
    S02_POISSON_ROUTING = c(high = "positive", low = "indeterminate", background = "negative"),
    S02_EXPECTATION_GRID = list(n_rows = 4L, n_expected_present = 2L),
    S02_LEAKAGE_ESTIMATE = c(R1 = 0.25, R2 = NA_real_),
    F02_MISSING_QC_LOSES_UNALIGNED = TRUE,
    F02_MISSING_BED_DISABLES_MASK = TRUE,

    S03_BLAST_FEATURES = list(margin = 30, n_species = 2L, entropy = 1, subject_len = 5000),
    S03_KRAKEN_FEATURES = c(r1_conf = 0.7, r1_distinct = 2, r2_conf = 0, r2_distinct = 0),
    S03_OUTER_UNION = c("r_blast", "r_competitor", "r_kraken", "r_qc"),
    S03_ARM_ISOLATION = TRUE,
    F03_FORBIDDEN_CONFIG_INJECTION = TRUE,

    S04_LOEO_CONTRACT = TRUE,
    S04_INNER_DONOR_CONTRACT = TRUE,
    S04_LOTO_CONTRACT = TRUE,
    S04_DETERMINISTIC_SPLITS = TRUE,
    F04_SHARED_RUN_SPLIT = TRUE,

    S05_TRAIN_ONLY_IMPUTER = c(train_median = 2, imputed_test = 2),
    S05_FIXED_BASELINE = list(chosen_feature = "signal", probabilities_valid = TRUE,
                              positive_ranked_higher = TRUE),
    S05_BALANCE_WEIGHTS = TRUE,
    S05_INNER_GRID_SELECTION = TRUE,
    S05_OOF_DONOR_ISOLATION = TRUE,
    S05_OOF_COVERAGE = TRUE,
    F05_MISSING_FEATURE_DROPPED = TRUE,
    F05_MODEL_FAILURE_SWALLOWED = TRUE,

    S06_HAND_CALCULATED_METRICS = list(auprc = 1, precision = c(1, 1), fdr = c(0, 0)),
    S06_TRUNCATED_COUNTS = c(included = 2L, excluded = 1L),
    S06_SAMPLE_TAXON_TRUTH = c(spA = 1L, spB = 0L, spC = 1L),
    S06_MISSING_TAXON_IS_FN = TRUE,
    S06_DEPTH_NORMALIZATION = 250000,
    F06_POOLED_OOF_THRESHOLD = TRUE,
    F06_NA_SPECIES_NEGATIVE_DROPPED = TRUE,

    S07_FOLD_PAIRING = 0.1,
    S07_INNER_MODEL_SELECTION = "glm",
    S07_H3_LOW_ABUNDANCE_DIRECTION = TRUE,
    S07_REQUIRED_ROWS = c(paste0("H", 1:12), "H5b", "H9b_classifier", "H9b_truth"),
    S07_MULTIPLICITY_SCOPE = TRUE,
    F07_OUTER_SELECTION_FALLBACK = TRUE,
    F07_SINGLE_ARM_OUTER_SELECTION = TRUE
  )
)

.test_registry <- data.frame(
  id = names(PIPELINE_TEST_PARAMETERS$expectations),
  stage = sub("^[SF]([0-9]{2}).*$", "\\1", names(PIPELINE_TEST_PARAMETERS$expectations)),
  kind = ifelse(startsWith(names(PIPELINE_TEST_PARAMETERS$expectations), "F"),
                "finding", "correctness"),
  stringsAsFactors = FALSE
)

.test_descriptions <- c(
  S01_RAGGED_PAF = "Ragged minimap2 PAF rows are read completely and best hits are retained",
  S01_HUMAN_COMPETITION = "Only strict human-score winners are depleted",
  S01_SAMPLE_SHEET_SCHEMA = "Missing required sample-sheet columns are rejected",
  S01_COMMAND_WIRING = "Truth, QC, BLAST and Kraken2 commands receive the intended inputs and flags",
  S01_EMPTY_OUTPUT_INCOMPLETE = "Empty terminal artifacts do not satisfy the resume guard",
  F01_RESUME_OMITS_PAF_QC = "Resume guard can declare completion without truth PAF or QC artifacts",
  S02_LABEL_ASSIGNMENT = "Synthetic reads receive the predeclared positive, negative and ambiguous labels",
  S02_QC_DEFINES_UNIVERSE = "The depleted QC read set defines the label universe",
  S02_POISSON_ROUTING = "Only low-detectability positives are routed to indeterminate",
  S02_EXPECTATION_GRID = "Expected sample-by-taxon truth is a complete library-by-taxon grid",
  S02_LEAKAGE_ESTIMATE = "Leakage is estimated independently from each run's negative barcode",
  F02_MISSING_QC_LOSES_UNALIGNED = "The missing-QC fallback omits unaligned surviving reads",
  F02_MISSING_BED_DISABLES_MASK = "A missing ambiguous BED silently disables region masking",
  S03_BLAST_FEATURES = "BLAST top-hit, margin, diversity and subject features are correct",
  S03_KRAKEN_FEATURES = "Kraken2 confidence and minimizer-count features are correct",
  S03_OUTER_UNION = "Independent classifier and QC outputs are outer-joined by read ID",
  S03_ARM_ISOLATION = "Single classifier arms contain no features from the other classifier",
  F03_FORBIDDEN_CONFIG_INJECTION = "Adversarial label/group fields can enter model_features through configuration",
  S04_LOEO_CONTRACT = "LOEO keeps donors disjoint, excludes unusable labels and tests each eligible row once",
  S04_INNER_DONOR_CONTRACT = "Inner validation donors partition only the outer-training donors",
  S04_LOTO_CONTRACT = "LOTO withholds the target taxon and partitions donor backgrounds",
  S04_DETERMINISTIC_SPLITS = "Fold generation is deterministic for a fixed seed",
  F04_SHARED_RUN_SPLIT = "LOEO can split one technical run across training and test when donors share a run",
  S05_TRAIN_ONLY_IMPUTER = "Imputation values are learned from training data only",
  S05_FIXED_BASELINE = "Fixed baseline selects the training signal and emits ordered probabilities",
  S05_BALANCE_WEIGHTS = "Donor balancing gives each donor equal total training weight",
  S05_INNER_GRID_SELECTION = "Inner CV selects the predeclared better grid row",
  S05_OOF_DONOR_ISOLATION = "Every model call has disjoint training and validation/test donors",
  S05_OOF_COVERAGE = "Every eligible read receives exactly one requested outer-fold prediction",
  F05_MISSING_FEATURE_DROPPED = "A configured but absent feature is silently removed before fitting",
  F05_MODEL_FAILURE_SWALLOWED = "A failed model fold is converted to all-NA predictions without failing the stage",
  S06_HAND_CALCULATED_METRICS = "A separable fixture returns the predeclared AUPRC, precision and FDR",
  S06_TRUNCATED_COUNTS = "With/without-unblocked reports contain the expected negative counts",
  S06_SAMPLE_TAXON_TRUTH = "A-priori sample-by-taxon truth overrides observed read labels",
  S06_MISSING_TAXON_IS_FN = "Expected-present taxa with no reads are inserted as false negatives",
  S06_DEPTH_NORMALIZATION = "Read counts are normalized by the matching library depth",
  F06_POOLED_OOF_THRESHOLD = "A held-out donor changes another donor's pooled score threshold",
  F06_NA_SPECIES_NEGATIVE_DROPPED = "Negative reads without a truth species disappear from sample-by-taxon aggregation",
  S07_FOLD_PAIRING = "Paired effects are aligned by fold identifier rather than row order",
  S07_INNER_MODEL_SELECTION = "Inner-CV evidence controls primary model selection despite opposite outer metrics",
  S07_H3_LOW_ABUNDANCE_DIRECTION = "H3 treats c4/c5 as lower abundance than c1/c2",
  S07_REQUIRED_ROWS = "The hypothesis driver emits every required primary, secondary and provenance row",
  S07_MULTIPLICITY_SCOPE = "Holm-Sidak adjustment is confined to H1-H6",
  F07_OUTER_SELECTION_FALLBACK = "Missing inner scores causes model selection on outer evaluation metrics",
  F07_SINGLE_ARM_OUTER_SELECTION = "H1 chooses its single-arm comparator on the same outer metrics it tests"
)

.test_registry$description <- unname(.test_descriptions[.test_registry$id])
PIPELINE_TEST_PARAMETERS$registry <- .test_registry
rm(.test_registry, .test_descriptions)