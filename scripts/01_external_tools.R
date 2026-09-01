## =============================================================================
## 01_external_tools.R  --  QC + ground-truth alignment + the two classifier arms
## -----------------------------------------------------------------------------
## Shells out to minimap2 / samtools / blastn / kraken2 once per library.
## Paths and parameters come from 00_config.R (single source of truth).
##
## KEY DESIGN POINTS
##   * minimap2 vs Zymo members AND vs human -> ground-truth evidence      [note A]
##   * minimap2 vs T2T-CHM13 -> human-competitor score (a covariate)       [note F]
##   * BLASTn and Kraken2 are run INDEPENDENTLY on the same QC'd reads.
##     Neither gates the other; both classify every read.                  [note D]
##   * QC keeps ONT end_reason so unblocked/truncated reads can be
##     flagged and reported separately downstream.                         [note E]
## =============================================================================

suppressWarnings(suppressMessages(library(data.table)))

## --- locate & load config (works whether run via Rscript or source()) --------
if (!exists("cfg")) {
  .sd <- local({
    a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1])))
    else "/home/tdinse/master_timo/projektarbeit/projektarbeit_ml/scripts"
  })
  source(file.path(.sd, "00_config.R"))
}

## Resolve an external binary or stop with guidance.
need_tool <- function(bin) {
  path <- Sys.which(bin)
  if (nzchar(path)) return(unname(path))
  stop(sprintf("External tool '%s' not found on PATH. Install it (conda/mamba) before running stage 01.", bin),
       call. = FALSE)
}

run_cmd <- function(bin, args, stdout_file = "") {
  args <- as.character(args)
  message("+ ", bin, " ", paste(args, collapse = " "),
          if (nzchar(stdout_file)) paste0(" > ", stdout_file) else "")
  ## system2 runs via the shell, so shQuote'd path/format args are honoured;
  ## redirection MUST go through stdout=, not a literal ">" argument.
  status <- system2(bin, args = args, stdout = stdout_file, stderr = "")
  if (is.integer(status) && status != 0L)
    stop(sprintf("Command failed (%s), exit status %d", bin, status), call. = FALSE)
  invisible(status)
}

## -----------------------------------------------------------------------------
## Load the sample sheet and validate the hardcoded inputs first.
## -----------------------------------------------------------------------------
stage01_preflight <- function() {
  gaps <- cfg_missing_paths()
  input_gaps <- gaps[gaps$key %in% c("sample_sheet", "zymo_refs_dir",
                                     "human_grch38", "human_t2t_chm13",
                                     "kraken2_db", "blast_db"), , drop = FALSE]
  if (nrow(input_gaps)) {
    message("Stage 01 blocked -- resolve these hardcoded inputs (see OPEN ITEMS):")
    print(input_gaps, row.names = FALSE)
    stop("Fix the paths in 00_config.R SECTION 1 and re-run.", call. = FALSE)
  }
  ## zymo_members.fasta is BUILT from zymo_refs_dir -- build it first if absent. [OI 1]
  if (!file.exists(cfg$paths$zymo_refs_fasta))
    stop("zymo_members.fasta not built. Run:  Rscript scripts/prepare_zymo_genomes.R\n",
         "then re-run stage 01. [OI 1]", call. = FALSE)
  ## every library's fastq must exist
  ss <- read_sample_sheet()
  miss_fq <- ss$library_id[!file.exists(ss$fastq)]
  if (length(miss_fq)) {
    message("Stage 01 blocked -- fastq not found for: ", paste(miss_fq, collapse = ", "))
    stop("Fix the `fastq` column in the sample sheet.", call. = FALSE)
  }
  invisible(TRUE)
}

read_sample_sheet <- function() {
  ss <- fread(cfg$paths$sample_sheet)
  required <- c("library_id", "donor", "barcode", "titration_level", "run_id", "fastq")
  miss <- setdiff(required, names(ss))
  if (length(miss)) stop("sample_sheet missing columns: ", paste(miss, collapse = ", "))
  ss
}

## -----------------------------------------------------------------------------
## Per-library steps. Each writes to work_dir/<library_id>/.
## -----------------------------------------------------------------------------

## Given the per-run ONT reports FOLDER (or a direct file) + the library's run_id,
## return the matching sequencing_summary_*<run_id>*.txt, or NA. [note E / OI 6]
resolve_seq_summary <- function(dir, run_id) {
  if (is.null(dir) || length(dir) != 1L || is.na(dir) || !nzchar(dir)) return(NA_character_)
  if (!dir.exists(dir)) return(if (file.exists(dir)) dir else NA_character_)  # direct file path ok
  hits <- list.files(dir, pattern = "^sequencing_summary.*\\.txt$", full.names = TRUE)
  if (!length(hits)) return(NA_character_)
  tag <- substr(gsub("-", "", as.character(run_id)), 1, 8)         # UUID prefix in the filename
  m <- hits[grepl(tag, basename(hits), fixed = TRUE)]
  if (length(m)) m[1] else if (length(hits) == 1L) hits[1] else NA_character_
}

## --- offline human depletion (R2 / note D) ----------------------------------
## Best hit per read from a PAF: read_id, score(=residue matches), identity, coverage.
.paf_best <- function(path) {
  if (!file.exists(path) || file.size(path) == 0) return(data.table(read_id = character()))
  ## minimap2 -c PAF has a VARIABLE number of trailing tag columns, so
  ## fread(select=1:12) STOPS EARLY at the first ragged line and silently truncates
  ## the file (was capping depletion at ~19k reads). Cut to the 12 fixed PAF
  ## columns first so EVERY read is parsed. [BUGFIX 2026-08-07]
  p <- fread(cmd = paste("cut -f1-12", shQuote(path)), header = FALSE, sep = "\t",
             col.names = c("qname","qlen","qstart","qend","strand","tname",
                           "tlen","tstart","tend","matches","blocklen","mapq"))
  p[, `:=`(identity = matches / pmax(blocklen, 1L),
           coverage = (qend - qstart) / pmax(qlen, 1L), score = matches)]
  p[order(-score), .(read_id = qname, score, identity, coverage)][, .SD[1L], by = read_id]
}

## Write a fastq with the given read IDs EXCLUDED (seqkit if present, else awk).
filter_fastq_exclude <- function(fastq, ids_file, out) {
  if (nzchar(Sys.which("seqkit"))) {
    run_cmd("seqkit", c("grep", "-v", "-f", shQuote(ids_file), shQuote(fastq)), stdout_file = out)
    return(invisible(out))
  }
  src <- fastq
  if (grepl("\\.gz$", fastq)) { src <- tempfile(fileext = ".fastq")
    run_cmd("zcat", shQuote(fastq), stdout_file = src) }
  awk_prog <- 'NR==FNR{ids[$1]=1;next} (FNR%4==1){h=substr($1,2);sub(/[ \\t].*/,"",h);keep=!(h in ids)} keep'
  run_cmd("awk", c(shQuote(awk_prog), shQuote(ids_file), shQuote(src)), stdout_file = out)
  invisible(out)
}

## Reads that align BETTER to human than to any Zymo member -> 'human_like',
## removed before classification. human_score = max(hg38, T2T) residue matches; a
## read absent from BOTH human PAFs can never be human_like. This is the SAME
## competition the stage-02 labeller uses (positive = Zymo beats human), so
## depletion and the labels are consistent by construction. NO identity/coverage
## gate -- minimap2 map-ont under-aligns SHORT human reads, so a coverage gate let
## ~83% human leak into BLAST/Kraken2. Pure + unit-testable. [note D / R2]
human_like_reads <- function(hg, tt, z) {
  parts <- list()
  if (nrow(hg)) parts[[length(parts) + 1L]] <- hg[, .(read_id, score)]
  if (nrow(tt)) parts[[length(parts) + 1L]] <- tt[, .(read_id, score)]
  if (!length(parts)) return(character(0))       # nothing aligned to human
  hum <- rbindlist(parts)[, .(human_score = max(score)), by = read_id]
  if (nrow(z)) hum[z, on = "read_id", zymo_score := i.score]
  if (!"zymo_score" %in% names(hum)) hum[, zymo_score := 0]
  hum[is.na(zymo_score), zymo_score := 0]
  hum[human_score > zymo_score, read_id]
}

## Split the library into the non-human universe (kept -> classified by BOTH arms)
## and the human_like set (removed). Writes human_ids.txt + nonhuman.fastq and
## returns the non-human fastq path. Adaptive sampling already depleted human live;
## this is a second, COMPETITION-based pass (see human_like_reads). [R2/D/E]
deplete_human <- function(lib, out_dir) {
  hg <- .paf_best(file.path(out_dir, "gt_human_grch38.paf"))   # GRCh38
  tt <- .paf_best(file.path(out_dir, "competitor_t2t.paf"))    # T2T-CHM13
  z  <- .paf_best(file.path(out_dir, "gt_zymo.paf"))
  human_ids <- human_like_reads(hg, tt, z)
  if (!length(human_ids)) return(lib$fastq)      # nothing aligns better to human -> keep all
  ids_file <- file.path(out_dir, "human_ids.txt")
  writeLines(human_ids, ids_file)
  nonhuman <- file.path(out_dir, "nonhuman.fastq")
  filter_fastq_exclude(lib$fastq, ids_file, nonhuman)
  message(sprintf("  human depletion (max(hg38,T2T) > Zymo): removed %d human-like reads -> %s",
                  length(human_ids), basename(nonhuman)))
  nonhuman
}

## (1) QC: length/quality per read + join ONT end_reason. [note E]
##     Uses NanoPlot/seqkit if available; otherwise a lightweight fastq scan.
qc_library <- function(lib, out_dir, reads_fastq = lib$fastq) {
  qc_out <- file.path(out_dir, "read_qc.tsv")
  seqkit <- Sys.which("seqkit")
  if (nzchar(seqkit)) {
    ## --only-id: these ONT fastq headers carry TAB-separated BAM tags (BC:Z:, qs:f:, ...).
    ## Without it, --name emits the whole tab-containing header and pushes length/GC/qual
    ## into trailing columns (header vs data column mismatch). --only-id keeps just the id.
    run_cmd(seqkit, c("fx2tab", "--header-line", "--only-id", "--name", "--length", "--avg-qual", "--gc",
                      shQuote(reads_fastq)), stdout_file = qc_out)
  } else {
    message("  seqkit not found -- emitting QC stub; fill with your QC tool. [OI 7]")
    fwrite(data.table(read_id = character(), read_len = integer(),
                      mean_q = numeric(), gc = numeric()), qc_out, sep = "\t")
  }
  ## Join end_reason/channel/start_time from the run's sequencing_summary. [note E / OI 6]
  ss_file <- resolve_seq_summary(lib$seq_summary, lib$run_id)
  if (!is.na(ss_file)) {
    ss <- fread(ss_file, select = c("read_id", "end_reason", "channel", "start_time"))
    fwrite(ss, file.path(out_dir, "end_reason.tsv"), sep = "\t")
  } else {
    message("  no sequencing_summary for ", lib$library_id,
            " (folder: ", lib$seq_summary, ") -- end_reason unavailable. [OI 6]")
  }
  qc_out
}

## (2) Ground-truth alignments (minimap2). NOT a classifier -- truth only. [note A]
groundtruth_align <- function(lib, out_dir, mm2, samtools, threads = 4L) {
  ## a) reads vs concatenated Zymo member genomes
  zymo_paf <- file.path(out_dir, "gt_zymo.paf")
  run_cmd(mm2, c("-cx", "map-ont", "-t", threads, "--secondary=no",
                 shQuote(cfg$paths$zymo_refs_fasta), shQuote(lib$fastq)), stdout_file = zymo_paf)
  ## b) reads vs human (GRCh38) -- the live-depletion reference. [OI 5]
  human_paf <- file.path(out_dir, "gt_human_grch38.paf")
  run_cmd(mm2, c("-cx", "map-ont", "-t", threads, "--secondary=no",
                 shQuote(cfg$paths$human_grch38), shQuote(lib$fastq)), stdout_file = human_paf)
  ## c) reads vs T2T-CHM13 -- the human-competitor covariate. [note F]
  t2t_paf <- file.path(out_dir, "competitor_t2t.paf")
  run_cmd(mm2, c("-cx", "map-ont", "-t", threads, "--secondary=no",
                 shQuote(cfg$paths$human_t2t_chm13), shQuote(lib$fastq)), stdout_file = t2t_paf)
  c(zymo = zymo_paf, human = human_paf, t2t = t2t_paf)
}

## (3a) Classifier arm A -- BLASTn vs core_nt, INDEPENDENT of Kraken2. [note D]
##      DECIDED 2026-08-03: FULL-RUN analysis only, no subsampling -- BLAST every
##      non-human read (shared human-depleted universe for all arms). [OI 11,12]
blast_classify <- function(lib, out_dir, blastn, makeblastdb, threads = 4L,
                           top_n = cfg$params$top_n_hits, reads_fastq = lib$fastq) {
  reads_fa <- file.path(out_dir, "reads.fasta")
  if (nzchar(Sys.which("seqkit"))) {
    run_cmd("seqkit", c("fq2fa", shQuote(reads_fastq), "-o", shQuote(reads_fa)))
  }
  blast_out <- file.path(out_dir, "blastn.tsv")
  fmt <- paste("6 qseqid sseqid staxid pident length mismatch gapopen",
               "qstart qend sstart send evalue bitscore qlen slen")
  ## -task megablast : made EXPLICIT (it is blastn's default) -- the fast, high-
  ##   identity algorithm appropriate for metagenomic read ID vs core_nt.
  ## -max_hsps 1     : keep only the single best HSP per query-subject pair. Long
  ##   ONT reads tile repetitive core_nt subjects into hundreds of HSPs/read (~726
  ##   observed), which (a) bloats output to ~180 GB/library and (b) lets the top-N
  ##   diversity features (n_species_topN / n_genera_topN / tax_entropy_topN /
  ##   bitscore_margin_species) collapse onto one subject. One best HSP per subject
  ##   is exactly what those top-N features assume -> lossless AND more correct.
  run_cmd(blastn, c("-query", shQuote(reads_fa), "-db", shQuote(cfg$paths$blast_db),
                    "-task", "megablast", "-max_hsps", "1",
                    "-max_target_seqs", top_n, "-num_threads", threads,
                    "-outfmt", shQuote(fmt), "-out", shQuote(blast_out)))
  blast_out
}

## (3b) Classifier arm B -- Kraken2, INDEPENDENT of BLASTn. [note D]
kraken2_classify <- function(lib, out_dir, kraken2, threads = 4L, reads_fastq = lib$fastq) {
  k2_out <- file.path(out_dir, "kraken2.out")
  k2_report <- file.path(out_dir, "kraken2.report")
  ## The DB hash is ~456G and this host has ~502G RAM, yet kraken2 is invoked once
  ## PER library (48x). WITHOUT --memory-mapping each call copies the whole 456G
  ## into process RSS (leaving almost no headroom -> OOM-prone, plus 48 reloads).
  ## --memory-mapping mmaps the hash and shares it via the OS page cache across all
  ## invocations; since DB ~= RAM it stays cached after the first warm-up, so it is
  ## OOM-safe and not slower here.
  run_cmd(kraken2, c("--db", shQuote(cfg$paths$kraken2_db), "--threads", threads,
                     "--memory-mapping",          # mmap 456G hash -> shared page cache across 48 libs (host RAM 502G)
                     "--confidence", "0",         # keep raw; confidence is a FEATURE downstream
                     "--report", shQuote(k2_report),
                     "--output", shQuote(k2_out), shQuote(reads_fastq)))
  c(out = k2_out, report = k2_report)
}

## -----------------------------------------------------------------------------
## Driver
## -----------------------------------------------------------------------------
## `threads` defaults to the PIPE_THREADS env var (falls back to 4 when unset, so
## tests / standalone calls are unchanged). Every shelled-out tool (minimap2,
## blastn, kraken2) inherits it -- set e.g. PIPE_THREADS=32 on a many-core host.

## A library's stage-01 work is COMPLETE when its terminal outputs exist and are
## non-empty. kraken2 is the LAST tool in the per-library sequence (blast -> kraken2),
## so a non-empty kraken2.out/report implies blastn.tsv finished too. This lets a
## re-launched run RESUME -- already-finished libraries are skipped instead of redone.
## Set PIPE_REDO=1 to force every library to be re-processed from scratch.
lib_complete <- function(out_dir) {
  need <- file.path(out_dir, c("reads.fasta", "blastn.tsv", "kraken2.out", "kraken2.report"))
  all(file.exists(need)) && all(file.info(need)$size > 0)
}

run_stage01 <- function(threads = as.integer(Sys.getenv("PIPE_THREADS", "4")), libraries = NULL) {
  stage01_preflight()
  cfg_init_dirs()
  mm2        <- need_tool("minimap2")
  samtools   <- need_tool("samtools")
  blastn     <- need_tool("blastn")
  makeblastdb<- Sys.which("makeblastdb")
  kraken2    <- need_tool("kraken2")

  ss <- read_sample_sheet()
  if (!is.null(libraries)) ss <- ss[library_id %in% libraries]

  for (i in seq_len(nrow(ss))) {
    lib <- as.list(ss[i])
    out_dir <- file.path(cfg$paths$work_dir, lib$library_id)
    if (lib_complete(out_dir) && !nzchar(Sys.getenv("PIPE_REDO", ""))) {
      message("== library ", lib$library_id, " (", i, "/", nrow(ss),
              ") -- already complete, skipping (set PIPE_REDO=1 to force) ==")
      next
    }
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    message("== library ", lib$library_id, " (", i, "/", nrow(ss), ") ==")

    groundtruth_align(lib, out_dir, mm2, samtools, threads)   # PAFs on ALL reads [note A / F]
    reads_fq <- deplete_human(lib, out_dir)                    # non-human universe [R2 / note D]
    qc_library(lib, out_dir, reads_fq)                         # QC on the non-human universe [note E]
    blast_classify(lib, out_dir, blastn, makeblastdb, threads, reads_fastq = reads_fq)  # arm A [note D]
    kraken2_classify(lib, out_dir, kraken2, threads, reads_fastq = reads_fq)            # arm B [note D]
  }
  message("Stage 01 complete -> ", cfg$paths$work_dir)
  invisible(TRUE)
}

if (sys.nframe() == 0L) {
  run_stage01()
}
