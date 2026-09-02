## =============================================================================
## 02_ground_truth_labels.R  --  three-class labels + Poisson floor + leakage
## -----------------------------------------------------------------------------
## Consumes minimap2 PAFs from stage 01 and produces the per-read label table.
##
## GROUND TRUTH IS NOT CIRCULAR [note A]:
##   positive     read aligns to a Zymo member genome above identity+coverage
##                cutoffs AND scores better there than against human.
##   negative     read aligns better to human, or aligns well to neither.
##   ambiguous    read's best Zymo hit falls in an rRNA operon / mobile element /
##                plasmid backbone / low-complexity region  -> EXCLUDED, not
##                forced into the binary.
##
## Then two corrections that change which positives 'count':
##   Poisson floor [note H]  species x level with P(>=1 genome) < poisson_p_min
##                           -> 'indeterminate' stratum, excluded from recall.
##   Leakage [note L]        c4/c5 are low-mass and share a flow cell with c1-c3,
##                           so cross-barcode leakage inflates apparent positives.
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

PAF_COLS <- c("qname", "qlen", "qstart", "qend", "strand", "tname",
              "tlen", "tstart", "tend", "matches", "blocklen", "mapq")

read_paf <- function(path) {
  if (!file.exists(path) || file.size(path) == 0)
    return(data.table(qname = character()))
  ## cut to the 12 fixed PAF columns first -- minimap2 -c appends a VARIABLE number
  ## of tag columns, so fread(select=1:12) stops early and truncates. [BUGFIX 2026-08-07]
  paf <- fread(cmd = paste("cut -f1-12", shQuote(path)), header = FALSE, sep = "\t",
               col.names = PAF_COLS)
  paf[, `:=`(identity = matches / pmax(blocklen, 1L),
             coverage = (qend - qstart) / pmax(qlen, 1L),
             score = matches)]  # residue matches used as a strand-agnostic score proxy
  paf
}

## Best alignment per read (max residue matches).
best_per_read <- function(paf) {
  if (!nrow(paf)) return(paf)
  paf[order(-score), .SD[1L], by = qname]
}

## Per-read residue-match score from a best-hit table, aligned to `ids`; 0 where
## the read is absent or the table is empty (robust to a missing PAF, e.g. T2T).
.best_score <- function(best, ids) {
  if (!nrow(best) || !("score" %in% names(best))) return(numeric(length(ids)))
  s <- best$score[match(ids, best$qname)]
  fifelse(is.na(s), 0, s)
}

## Map a Zymo contig name -> member species. Prefers an explicit lookup file;
## otherwise falls back to the contig-name prefix. [OI 1]
map_contig_to_species <- function(tnames) {
  lut_path <- file.path(cfg$paths$project_root, "ref", "zymo_contig2species.tsv")
  if (file.exists(lut_path)) {
    lut <- fread(lut_path)  # columns: contig, species
    return(lut$species[match(tnames, lut$contig)])
  }
  message("  ref/zymo_contig2species.tsv missing -- deriving species from contig prefix. [OI 1]")
  sub("[:_].*$", "", tnames)
}

## Overlap best Zymo hits with the ambiguous-region BED (rRNA/mobile/plasmid). [note A / OI 2]
flag_ambiguous_regions <- function(zbest) {
  bed_path <- cfg$params$ambiguous_bed
  if (!file.exists(bed_path)) {
    message("  ambiguous_bed missing -- region-based ambiguity OFF; relying on complexity only. [OI 2]")
    return(rep(FALSE, nrow(zbest)))
  }
  bed <- fread(bed_path, header = FALSE, select = 1:3,
               col.names = c("tname", "start", "end"))
  setkey(bed, tname, start, end)
  ## zbest is the per-read label table; its Zymo best-hit coords are zymo_* prefixed.
  ## Reads with no Zymo hit (NA coords) cannot fall in a Zymo ambiguous region.
  q <- zbest[, .(tname = zymo_tname, start = zymo_tstart, end = zymo_tend, .row = .I)]
  q <- q[!is.na(tname) & !is.na(start) & !is.na(end)]
  if (!nrow(q)) return(rep(FALSE, nrow(zbest)))
  ov <- foverlaps(q, bed, type = "any", which = TRUE, nomatch = NA)
  flagged <- unique(q$.row[ov$xid[!is.na(ov$yid)]])
  seq_len(nrow(zbest)) %in% flagged
}

## -----------------------------------------------------------------------------
## Per-library label assignment.
## -----------------------------------------------------------------------------
label_library <- function(lib) {
  out_dir <- file.path(cfg$paths$work_dir, lib$library_id)
  z  <- best_per_read(read_paf(file.path(out_dir, "gt_zymo.paf")))
  h  <- best_per_read(read_paf(file.path(out_dir, "gt_human_grch38.paf")))  # GRCh38 [OI 5]
  h2 <- best_per_read(read_paf(file.path(out_dir, "competitor_t2t.paf")))   # T2T -- human_score = max(hg38,T2T), same rule as depletion [R2]

  ## Analysis universe = reads that SURVIVED human depletion (QC ran on the
  ## non-human fastq in stage 01). This makes the negative class the interesting
  ## one -- non-Zymo bacteria, residual human, and unaligned survivors -- not
  ## trivially-removed human. Falls back to the PAF union if QC is unavailable. [R2]
  qc_path <- file.path(out_dir, "read_qc.tsv")
  reads <- character(0)
  if (file.exists(qc_path) && file.size(qc_path) > 0) {
    qc <- fread(qc_path)
    if (nrow(qc)) reads <- sub("\\s.*$", "", as.character(qc[[1]]))
  }
  if (!length(reads)) reads <- unique(c(z$qname, h$qname))
  dt <- data.table(read_id = reads, library_id = lib$library_id,
                   donor = lib$donor, titration_level = lib$titration_level,
                   run_id = lib$run_id, barcode = lib$barcode,
                   concentration = if (!is.null(lib$concentration)) as.numeric(lib$concentration) else NA_real_)

  zi <- z[match(dt$read_id, qname)]
  dt[, `:=`(zymo_score = .best_score(z, read_id),
            zymo_ident = zi$identity, zymo_cov = zi$coverage,
            zymo_tname = zi$tname, zymo_tstart = zi$tstart, zymo_tend = zi$tend,
            ## human score = max(hg38, T2T) -- IDENTICAL to the stage-01 depletion
            ## competition, so 'beats human' here cannot disagree with what was depleted. [R2]
            human_score = pmax(.best_score(h, read_id), .best_score(h2, read_id)))]
  dt[, species := map_contig_to_species(zymo_tname)]

  zymo_pass <- !is.na(dt$zymo_ident) &
    dt$zymo_ident >= cfg$params$gt_min_identity &
    dt$zymo_cov   >= cfg$params$gt_min_coverage
  beats_human <- dt$zymo_score > dt$human_score + cfg$params$gt_human_margin
  in_ambig <- flag_ambiguous_regions(dt)

  dt[, label := "negative"]
  dt[zymo_pass & beats_human, label := "positive"]
  dt[in_ambig, label := "ambiguous"]   # excluded region overrides binary [note A]
  dt[]
}

## -----------------------------------------------------------------------------
## Poisson expected-copy-number floor -> indeterminate stratum. [note H]
## -----------------------------------------------------------------------------
## Ground-truth detectability from wet-lab inputs:
##   E[genomes] = concentration (library input cells, from the SAMPLE SHEET)
##                * rel_abundance (per species, from the CoA)
##                * extraction_eff * fraction_loaded
##   P(>=1 genome) = 1 - exp(-E).  A (library x taxon) below poisson_p_min is not
##   reliably present, so its would-be positives are routed to the 'indeterminate'
##   stratum and excluded from recall accounting. extraction_eff/fraction_loaded
##   default to 1 (an upper bound on detectability) until measured. [OI 3]
apply_poisson_floor <- function(labels) {
  coa_path <- cfg$paths$zymo_coa
  if (grepl("<FILL_IN>", coa_path) || !file.exists(coa_path)) {
    message("Poisson floor SKIPPED -- no CoA (rel_abundance) file. [OI 3]")
    labels[, `:=`(expected_cells = NA_real_, p_detect = NA_real_, indeterminate = FALSE)]
    return(labels)
  }
  eff <- cfg$params$extraction_eff; frac <- cfg$params$fraction_loaded
  if (is.na(eff) || is.na(frac))
    message("Poisson floor: extraction_eff/fraction_loaded unset -> using 1 (upper bound). [OI 3]")
  eff <- if (is.na(eff)) 1 else eff; frac <- if (is.na(frac)) 1 else frac

  coa <- fread(coa_path)  # columns: species, rel_abundance
  ra  <- coa$rel_abundance[match(labels$species, coa$species)]
  ## per-library input cells from the sample-sheet `concentration`; fall back to
  ## cells_per_level[titration_level] if the column is absent.
  cells <- labels$concentration
  if (is.null(cells) || all(is.na(cells)))
    cells <- cfg$params$cells_per_level[labels$titration_level]
  en <- cells * ra * eff * frac
  labels[, `:=`(expected_cells = en,
                p_detect       = 1 - exp(-en),
                indeterminate  = !is.na(en) & (1 - exp(-en)) < cfg$params$poisson_p_min)]
  ## keep the pre-floor call: the negative barcodes have E[N] = 0, so the line below
  ## forces every Zymo-positive control read to 'indeterminate' and the leakage
  ## diagnostic would otherwise be zero by construction. [item 3]
  labels[, label_prefloor := label]
  labels[label == "positive" & indeterminate == TRUE, label := "indeterminate"]
  labels[]
}

## A-priori (library x taxon) ground-truth EXPECTATION, from the sample sheet
## `concentration` x CoA `rel_abundance` alone (no reads needed). This is the
## "total cells x rel_abundance >= ~1" truth: expected_present = TRUE when the
## taxon clears the Poisson floor. Written to work/expected_sample_taxon.tsv and
## used as the sample x taxon truth in stage 06. [note H / note G]
expected_sample_taxon <- function() {
  coa_path <- cfg$paths$zymo_coa
  if (grepl("<FILL_IN>", coa_path) || !file.exists(coa_path)) return(NULL)
  coa <- fread(coa_path)                 # species, rel_abundance
  ss  <- fread(cfg$paths$sample_sheet)   # library_id, donor, titration_level, concentration
  eff <- cfg$params$extraction_eff; frac <- cfg$params$fraction_loaded
  eff <- if (is.na(eff)) 1 else eff; frac <- if (is.na(frac)) 1 else frac
  grid <- CJ(library_id = ss$library_id, species = coa$species)
  grid[ss,  on = "library_id", `:=`(donor = i.donor, titration_level = i.titration_level,
                                    concentration = as.numeric(i.concentration))]
  grid[coa, on = "species", rel_abundance := i.rel_abundance]
  grid[, expected_cells   := concentration * rel_abundance * eff * frac]
  grid[, p_detect         := 1 - exp(-expected_cells)]
  grid[, expected_present := !is.na(p_detect) & p_detect >= cfg$params$poisson_p_min]
  grid[]
}

## -----------------------------------------------------------------------------
## Cross-barcode leakage correction for low-mass c4/c5. [note L / OI T3]
## -----------------------------------------------------------------------------
## Model: within a run, a fraction `leakage_rate` of reads in a barcode actually
## originate from other barcodes (index hopping / demux crosstalk). For low-mass
## c4/c5 this can dominate the true signal. Until the per-run demux error rate is
## known [OI T3], we (a) estimate an upper bound from the NEGATIVE barcode's Zymo
## positives on the same run, and (b) flag suspected leaked positives so they can
## be down-weighted or sensitivity-analysed.
## -----------------------------------------------------------------------------
## Cross-barcode leakage / contamination DIAGNOSTIC. [F7 decision / item 2 / OI T3]
## -----------------------------------------------------------------------------
## The negative (c0) barcode is the donor's no-spike background, run through the
## SAME human depletion as every library, so its surviving reads are endogenous
## microbiome + residual human + reagent kitome -- NOT pure human. Any Zymo-
## positive read there is therefore leakage/contamination, and its positive rate
## is a principled per-run UPPER BOUND on cross-barcode leakage into the low-mass
## c4/c5 barcodes. We REPORT this (plus the negative's read makeup) as a diagnostic
## for a downstream sensitivity analysis rather than RANDOMLY flagging individual
## reads (which was inert and unprincipled). -> work/leakage_estimates.tsv.
estimate_leakage <- function(labels) {
  ## Count on the PRE-floor call. Control barcodes have zero input, so p_detect = 0 and
  ## apply_poisson_floor() reroutes all of their positives to 'indeterminate'; counting
  ## the post-floor label returns 0 for every run regardless of the data. [item 3]
  lab <- if ("label_prefloor" %in% names(labels)) "label_prefloor" else "label"
  runs <- unique(labels$run_id)
  out <- rbindlist(lapply(runs, function(r) {
    neg <- labels[run_id == r & titration_level == "negative"]
    npos <- sum(neg[[lab]] == "positive")
    data.table(run_id = r,
               negative_reads         = nrow(neg),
               negative_background    = sum(neg[[lab]] == "negative"),  # endogenous + residual human
               negative_zymo_positive = npos,                           # contamination / leakage FP
               leakage_upper_bound    = if (nrow(neg)) npos / nrow(neg) else NA_real_,
               ## share of the run's own Zymo signal that reached a foreign barcode
               leakage_vs_positive_pool = local({
                 p <- sum(labels[run_id == r][[lab]] == "positive")
                 if (p) npos / p else NA_real_
               }))
  }), fill = TRUE)
  ## leakage is cutoff-DEPENDENT (built from labels) -> per-run folder. [two-run design]
  cov_path <- cfg$paths$leakage_table %||% file.path(cfg$paths$work_dir, "leakage_estimates.tsv")
  fwrite(out, cov_path, sep = "\t")
  message("  leakage / contamination diagnostic -> ", cov_path)
  print(out)
  invisible(out)
}

## -----------------------------------------------------------------------------
## Sample x taxon coverage: genome breadth + evenness per (library x species) from
## the minimap2 Zymo alignments -- the single strongest real-vs-artefact
## discriminator [note F/G]. Written to work/sample_taxon_coverage.tsv, joined by
## stage 06. Computed directly from the PAF intervals (no samtools/BAM). [R3]
## -----------------------------------------------------------------------------
compute_sample_taxon_coverage <- function() {
  ss <- fread(cfg$paths$sample_sheet)
  rows <- vector("list", nrow(ss))
  for (i in seq_len(nrow(ss))) {
    lib <- as.list(ss[i])
    paf <- read_paf(file.path(cfg$paths$work_dir, lib$library_id, "gt_zymo.paf"))
    if (!nrow(paf)) next
    ## per-contig breadth/evenness (one length each), then aggregate to species;
    ## for collapsed E. coli (5 strain contigs) take the best-covered strain
    ## rather than mixing coordinate spaces. [R3 / R5]
    cov_contig <- paf[, coverage_stats_from_intervals(tstart, tend, tlen), by = .(tname)]
    cov_contig[, species := map_contig_to_species(tname)]
    cov <- cov_contig[, .(genome_breadth = max(genome_breadth, na.rm = TRUE),
                          coverage_evenness = coverage_evenness[which.max(genome_breadth)]),
                      by = species]
    cov[, `:=`(donor = lib$donor, titration_level = lib$titration_level)]
    rows[[i]] <- cov
  }
  out <- rbindlist(rows, fill = TRUE)
  if (!nrow(out)) { message("  no Zymo alignments -> no coverage table."); return(invisible(NULL)) }
  setcolorder(out, c("donor", "titration_level", "species", "genome_breadth", "coverage_evenness"))
  cov_path <- file.path(cfg$paths$work_dir, "sample_taxon_coverage.tsv")
  fwrite(out, cov_path, sep = "\t")
  message("  sample x taxon coverage -> ", cov_path, " (", nrow(out), " rows)")
  invisible(out)
}

## -----------------------------------------------------------------------------
## Driver
## -----------------------------------------------------------------------------
run_stage02 <- function() {
  cfg_init_dirs()
  ss <- fread(cfg$paths$sample_sheet)
  pb <- make_progress(nrow(ss), "-- labelling libraries")
  all_labels <- rbindlist(lapply(seq_len(nrow(ss)), function(i) {
    r <- label_library(as.list(ss[i])); pb$tick(i); r
  }), fill = TRUE)
  pb$done()

  all_labels <- apply_poisson_floor(all_labels)
  estimate_leakage(all_labels)   # per-run leakage/contamination diagnostic [F7 / item 2]

  ## a-priori sample x taxon ground-truth expectation (concentration x rel_abundance)
  est <- expected_sample_taxon()
  if (!is.null(est)) {
    est_path <- file.path(cfg$paths$work_dir, "expected_sample_taxon.tsv")
    fwrite(est, est_path, sep = "\t")
    message("  expected sample x taxon truth -> ", est_path,
            "  (", sum(est$expected_present), "/", nrow(est), " library x taxon expected present)")
  }

  compute_sample_taxon_coverage()   # breadth/evenness per (library x species) [R3]

  fwrite(all_labels, cfg$paths$labels_table, sep = "\t")
  message("Label summary:"); print(all_labels[, .N, by = label])
  message("Stage 02 complete -> ", cfg$paths$labels_table)
  invisible(all_labels)
}

if (sys.nframe() == 0L) run_stage02()
