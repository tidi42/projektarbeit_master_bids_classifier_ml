## =============================================================================
## 03_build_features.R  --  taxon-agnostic per-read feature table
## -----------------------------------------------------------------------------
## Parses stage-01 outputs into one feature row per read.
##
## HARD RULE [note B]: taxon IDENTITY never enters the model matrix. Species/taxid
## are retained ONLY as grouping keys (sample x taxon aggregation, Poisson floor).
## Margin / entropy / count features USE taxonomy to be COMPUTED but are themselves
## taxon-agnostic, and are explicitly endorsed in note F.
##
## Feature blocks are the ones declared in 00_config.R $feature_blocks:
##   blast_core, blast_margin, human_competitor, kraken2, complexity, read_qc
## Residual barcode/adapter content is deliberately NOT a model feature [note F];
## it is computed separately only for the H7 exploratory test (stage 07). [OI T1]
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

BLAST_COLS <- c("qseqid", "sseqid", "staxid", "pident", "length", "mismatch",
                "gapopen", "qstart", "qend", "sstart", "send", "evalue",
                "bitscore", "qlen", "slen")

## Resolve staxid -> (species, genus) taxon groups for the count/entropy/margin
## features. Uses ref/taxid_lineage.tsv if present; otherwise treats staxid itself
## as the species group (fine for counts/entropy). [OI 9]
taxonomy_resolver <- function() {
  lut_path <- file.path(cfg$paths$project_root, "ref", "taxid_lineage.tsv")
  if (file.exists(lut_path)) {
    lut <- fread(lut_path)  # columns: taxid, species, genus[, rank]
    ## NOTE: the resolver's argument is deliberately NOT named `taxid` -- inside
    ## data.table's `lut[...]`, a bare `taxid` symbol resolves to lut's own
    ## `taxid` column, so `match(taxid, taxid)` would silently match lut against
    ## itself instead of against the queried staxids.
    function(query_taxid) lut[match(query_taxid, taxid), .(species, genus)]
  } else {
    message("  ref/taxid_lineage.tsv missing -- using staxid as the species group. [OI 9]")
    function(taxid) data.table(species = as.character(taxid),
                               genus = as.character(taxid))
  }
}

## --- BLAST core + margin features [note F] -----------------------------------
parse_blast <- function(path, top_n = cfg$params$top_n_hits, taxres) {
  if (!file.exists(path) || file.size(path) == 0) return(data.table(read_id = character()))
  b <- fread(path, header = FALSE, col.names = BLAST_COLS)
  lin <- taxres(b$staxid)
  b[, `:=`(species = lin$species, genus = lin$genus,
           aln_fraction = length / pmax(qlen, 1L),
           query_cov = (qend - qstart + 1L) / pmax(qlen, 1L))]
  setorder(b, qseqid, -bitscore)
  b[, hit_rank := seq_len(.N), by = qseqid]
  topn <- b[hit_rank <= top_n]

  ## top-hit core features
  core <- b[hit_rank == 1L, .(
    read_id = qseqid, qlen, pident, evalue, mismatch,
    aln_fraction, query_cov, bitscore,
    top_species = species, top_genus = genus,
    subject_genome_len = slen              # taxon-agnostic subject property [note B]
  )]

  ## margin / diversity features over the top-N hits
  margin <- topn[, {
    best_bs <- bitscore[1L]
    other_sp <- bitscore[species != species[1L]]
    .(bitscore_margin_species = best_bs - (if (length(other_sp)) max(other_sp) else 0),
      n_species_topN = uniqueN(species),
      n_genera_topN  = uniqueN(genus),
      tax_entropy_topN = shannon_entropy(species))
  }, by = .(read_id = qseqid)]

  merge(core, margin, by = "read_id", all = TRUE)
}

## --- Human-competitor score from the T2T-CHM13 PAF [note F] -------------------
parse_competitor <- function(paf_path) {
  if (!file.exists(paf_path) || file.size(paf_path) == 0) return(data.table(read_id = character()))
  ## cut to the 12 fixed PAF columns first (variable tag columns else truncate fread). [BUGFIX 2026-08-07]
  p <- fread(cmd = paste("cut -f1-12", shQuote(paf_path)), header = FALSE, sep = "\t",
             col.names = c("qname","qlen","qstart","qend","strand","tname",
                           "tlen","tstart","tend","matches","blocklen","mapq"))
  p[, `:=`(human_bitscore = matches, human_pident = matches / pmax(blocklen, 1L))]
  p[order(-human_bitscore), .(read_id = qname, human_bitscore, human_pident)][, .SD[1L], by = read_id]
}

## --- Kraken2 internals [note F] ----------------------------------------------
## kraken2 --output line: C/U  read_id  taxid  read_len  "taxid:count taxid:count ..."
parse_kraken2 <- function(out_path) {
  if (!file.exists(out_path) || file.size(out_path) == 0) return(data.table(read_id = character()))
  k <- fread(out_path, header = FALSE, sep = "\t",
             col.names = c("status", "read_id", "k2_taxid", "k2_readlen", "kmer_lca"))
  parse_line <- function(s, assigned_taxid) {
    toks <- strsplit(s, " ", fixed = TRUE)[[1]]
    toks <- toks[nzchar(toks)]
    if (!length(toks)) return(list(total = 0L, taxon = 0L, distinct = 0L))
    parts <- tstrsplit(toks, ":", fixed = TRUE)
    tx <- parts[[1]]; cnt <- as.integer(parts[[2]])
    total <- sum(cnt, na.rm = TRUE)
    taxon <- sum(cnt[tx == as.character(assigned_taxid)], na.rm = TRUE)
    list(total = total, taxon = taxon, distinct = uniqueN(tx[tx != "0"]))
  }
  agg <- k[, {
    pl <- parse_line(kmer_lca, k2_taxid)
    .(k2_conf = if (pl$total > 0) pl$taxon / pl$total else 0,
      k2_kmers_taxon_frac = if (pl$total > 0) pl$taxon / pl$total else 0,
      k2_distinct_minimizers = pl$distinct)   # approx; exact needs --report-minimizer-data [OI 9]
  }, by = read_id]
  agg
}

## --- Read QC + complexity + end_reason [notes E,F] ---------------------------
parse_qc <- function(out_dir) {
  qc_path <- file.path(out_dir, "read_qc.tsv")
  qc <- if (file.exists(qc_path) && file.size(qc_path) > 0)
    fread(qc_path, header = TRUE) else data.table(read_id = character())
  ## normalise seqkit fx2tab column names to canonical ones.
  ## fx2tab --only-id emits "#id"; plain --name emits "#name" (keep both for back-compat).
  ren <- c("#id" = "read_id", "id" = "read_id", "#name" = "read_id", "name" = "read_id",
           "length" = "read_len",
           "avg-qual" = "mean_q", "avg.qual" = "mean_q", "GC" = "gc", "GC(%)" = "gc")
  for (nm in names(ren)) if (nm %in% names(qc)) setnames(qc, nm, ren[[nm]])

  ## end_reason -> unblock flag (adaptive-sampling truncated reads) [note E]
  er_path <- file.path(out_dir, "end_reason.tsv")
  if (file.exists(er_path)) {
    er <- fread(er_path)
    qc <- merge(qc, er[, .(read_id, end_reason)], by = "read_id", all.x = TRUE)
    qc[, end_reason_unblock := as.integer(grepl("unblock", end_reason %||% "", ignore.case = TRUE))]
  } else {
    qc[, end_reason_unblock := NA_integer_]
  }

  ## complexity features require the read sequence; compute if reads.fasta exists.
  fa <- file.path(out_dir, "reads.fasta")
  if (file.exists(fa) && requireNamespace("Biostrings", quietly = TRUE)) {
    seqs <- Biostrings::readDNAStringSet(fa)
    cx <- data.table(read_id = sub("\\s.*$", "", names(seqs)),
                     gc2 = as.numeric(Biostrings::letterFrequency(seqs, "GC", as.prob = TRUE)),
                     dust_score = vapply(as.character(seqs), dust_score, numeric(1)),
                     homopolymer_frac = vapply(as.character(seqs), homopolymer_frac, numeric(1)))
    qc <- merge(qc, cx, by = "read_id", all.x = TRUE)
    if (!"gc" %in% names(qc)) qc[, gc := gc2]
    qc[, gc2 := NULL]
  } else {
    message("  complexity features skipped for ", basename(out_dir),
            " (need reads.fasta + Biostrings). [OI 11]")
    qc[, `:=`(dust_score = NA_real_, homopolymer_frac = NA_real_)]
    if (!"gc" %in% names(qc)) qc[, gc := NA_real_]
  }
  qc
}

## -----------------------------------------------------------------------------
## Assemble per-library feature rows and merge with labels.
## -----------------------------------------------------------------------------
features_for_library <- function(lib, taxres) {
  out_dir <- file.path(cfg$paths$work_dir, lib$library_id)
  blast <- parse_blast(file.path(out_dir, "blastn.tsv"), taxres = taxres)
  comp  <- parse_competitor(file.path(out_dir, "competitor_t2t.paf"))
  k2    <- parse_kraken2(file.path(out_dir, "kraken2.out"))
  qc    <- parse_qc(out_dir)

  Reduce(function(a, b) merge(a, b, by = "read_id", all = TRUE),
         list(blast, comp, k2, qc))[
           , library_id := lib$library_id][]
}

run_stage03 <- function() {
  cfg_init_dirs()
  taxres <- taxonomy_resolver()
  ss <- fread(cfg$paths$sample_sheet)
  pb <- make_progress(nrow(ss), "-- building features per library")
  feats <- rbindlist(lapply(seq_len(nrow(ss)), function(i) {
    r <- features_for_library(as.list(ss[i]), taxres); pb$tick(i); r
  }), fill = TRUE)
  pb$done()

  ## human-competitor derived margin [note F]
  if (!"human_bitscore" %in% names(feats)) feats[, human_bitscore := NA_real_]
  if (!"human_pident"  %in% names(feats)) feats[, human_pident := NA_real_]
  feats[, human_minus_best_margin := fifelse(is.na(human_bitscore), NA_real_, human_bitscore - bitscore)]

  ## join labels/strata from stage 02
  labels <- fread(cfg$paths$labels_table)
  ft <- merge(labels, feats, by = c("read_id", "library_id"), all.x = TRUE)

  ## persist (parquet if arrow present, else compressed tsv)
  if (requireNamespace("arrow", quietly = TRUE)) {
    arrow::write_parquet(ft, cfg$paths$feature_table)
  } else {
    alt <- sub("\\.parquet$", ".tsv.gz", cfg$paths$feature_table)
    fwrite(ft, alt, sep = "\t"); cfg$paths$feature_table <<- alt
    message("  arrow not installed -> wrote ", alt)
  }
  message("Stage 03 complete: ", nrow(ft), " reads x ", ncol(ft), " cols.")
  invisible(ft)
}

## Canonical model-feature selector: taxon-agnostic columns for a classifier arm. [notes B,D]
model_features <- function(arm = "combined") {
  blocks <- cfg$classifier_arms[[arm]]
  feats <- unique(c(unlist(cfg$feature_blocks[blocks]), cfg$subject_props))
  ## belt-and-braces: never let identity columns through
  setdiff(feats, c("species", "genus", "top_species", "top_genus", "staxid", "k2_taxid"))
}

if (sys.nframe() == 0L) run_stage03()
