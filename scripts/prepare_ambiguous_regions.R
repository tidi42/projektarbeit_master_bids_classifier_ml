## =============================================================================
## prepare_ambiguous_regions.R  --  build ref/zymo_ambiguous_regions.bed from
##                                  the supplier's 16S/18S rRNA references
## -----------------------------------------------------------------------------
## rRNA genes are highly conserved across (and even beyond) the Zymo community,
## so a read aligning there doesn't reliably confirm species identity -- ground
## truth would become circular if these loci counted as ordinary evidence.
## [note A / OI 2]
##
## For each species: take ONE clean copy of its supplier-provided 16S/18S
## sequence (ssrRNAs/*.fasta -- several files ship 3-7 near-identical operon
## copies; a few have a malformed header that swallows the sequence -- the
## first intact copy is enough) and align it back against that SAME species'
## own genome with secondary alignments explicitly ENABLED (`--secondary=yes`,
## high `-N`, relaxed `-p`), so every rRNA operon copy actually present in the
## genome is recovered -- not just minimap2's single best/primary hit. Hits are
## filtered by the same identity/coverage cutoffs used for ground truth
## (cfg$params$gt_min_identity/gt_min_coverage), then overlapping hits are
## merged into non-redundant BED intervals per contig.
##
## Scope: this marks the 16S/18S gene body itself (the only sequence the
## supplier provided), not the full rRNA operon (16S-ITS-23S-5S) -- a complete
## operon-level mask would additionally need something like barrnap.
##
## Idempotent: skips entirely if cfg$params$ambiguous_bed already exists;
## delete the file to force a rebuild.
## =============================================================================

suppressWarnings(suppressMessages(library(data.table)))
.sd <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1])))
  else "/home/tdinse/master_timo/projektarbeit/projektarbeit_ml/scripts"
})
if (!exists("cfg")) source(file.path(.sd, "00_config.R"))
## canonical_species() (for plasmid tnames) comes from prepare_zymo_genomes.R
if (!exists("canonical_species")) source(file.path(.sd, "prepare_zymo_genomes.R"))

need_tool <- function(bin) {
  path <- Sys.which(bin)
  if (nzchar(path)) return(unname(path))
  stop(sprintf("External tool '%s' not found on PATH.", bin), call. = FALSE)
}

## ---- minimal FASTA reader: id + concatenated sequence per record ----------
## A record whose header line swallowed its sequence (missing newline in the
## source file) comes back with nchar == 0 and is dropped by first_valid_query().
read_fasta_records <- function(path) {
  lines <- suppressWarnings(readLines(path))  # source files often lack a final newline
  is_header <- startsWith(lines, ">")
  n <- sum(is_header)
  grp <- cumsum(is_header)
  ids <- sub("^>(\\S+).*$", "\\1", lines[is_header])
  seq_list <- split(lines[!is_header], factor(grp[!is_header], levels = seq_len(n)))
  seqs <- vapply(seq_list, paste, character(1L), collapse = "")
  data.table(id = ids, seq = unname(seqs), nchar = nchar(unname(seqs)))
}

first_valid_query <- function(path, min_len = 200L) {
  recs <- read_fasta_records(path)[nchar >= min_len]
  if (!nrow(recs)) stop("No usable (non-empty) rRNA sequence in ", path, call. = FALSE)
  recs[1L]
}

PAF_COLS <- c("qname", "qlen", "qstart", "qend", "strand", "tname",
              "tlen", "tstart", "tend", "matches", "blocklen", "mapq")

## ---- align one query sequence against one genome, secondaries ENABLED -----
## Secondary alignments must stay on: we feed in one query per species, but a
## genome usually carries several near-identical rRNA copies, and we want ALL
## of their locations, not just minimap2's single best (primary) hit.
align_rrna_hits <- function(query_id, query_seq, genome_path, mm2,
                            min_identity = cfg$params$gt_min_identity,
                            min_coverage = cfg$params$gt_min_coverage) {
  if (!file.exists(genome_path)) stop("genome fasta not found: ", genome_path, call. = FALSE)
  qfa <- tempfile(fileext = ".fasta"); paf <- tempfile(fileext = ".paf")
  on.exit(unlink(c(qfa, paf)), add = TRUE)
  writeLines(c(paste0(">", query_id), query_seq), qfa)

  status <- system2(mm2, c("-cx", "map-ont", "--secondary=yes", "-N", "50", "-p", "0.5",
                          shQuote(genome_path), shQuote(qfa)),
                    stdout = paf, stderr = "")
  if (is.integer(status) && status != 0L)
    stop("minimap2 failed aligning ", query_id, " (exit ", status, ")", call. = FALSE)
  if (!file.exists(paf) || file.size(paf) == 0) return(data.table(tname = character()))

  ## cut to the 12 fixed PAF columns first (variable tag columns else truncate fread). [BUGFIX 2026-08-07]
  paf_dt <- fread(cmd = paste("cut -f1-12", shQuote(paf)), header = FALSE, sep = "\t", col.names = PAF_COLS)
  paf_dt[, `:=`(identity = matches / pmax(blocklen, 1L),
               coverage = (qend - qstart) / pmax(qlen, 1L))]
  paf_dt[identity >= min_identity & coverage >= min_coverage, .(tname, tstart, tend)]
}

## ---- merge overlapping/adjacent BED intervals, per contig -----------------
merge_intervals <- function(bed) {
  if (!nrow(bed)) return(bed)
  setorder(bed, tname, tstart)
  merge_grp <- function(starts, ends) {
    grp <- integer(length(starts)); grp[1] <- 1L; cur_end <- ends[1]
    if (length(starts) > 1) for (i in 2:length(starts)) {
      if (starts[i] > cur_end) { grp[i] <- grp[i - 1L] + 1L; cur_end <- ends[i] }
      else { grp[i] <- grp[i - 1L]; cur_end <- max(cur_end, ends[i]) }
    }
    grp
  }
  bed[, grp := merge_grp(tstart, tend), by = tname]
  out <- bed[, .(tstart = min(tstart), tend = max(tend)), by = .(tname, grp)]
  out[, grp := NULL][]
}

## -----------------------------------------------------------------------------
## Species -> (ssrRNA fasta, genome fasta). ssrRNA queries come from the D6331
## bundle (cfg$paths$zymo_ssrrna_dir). Genomes are the COLLAPSED, species-named
## records built by prepare_zymo_genomes.R into ref/genomes/<species>.fasta, so
## the BED coordinates line up with ref/zymo_members.fasta. [OI 1,2]
## -----------------------------------------------------------------------------
SSRRNA_DIR  <- cfg$paths$zymo_ssrrna_dir
GENOMES_DIR <- file.path(cfg$paths$project_root, "ref", "genomes")

RRNA_SOURCES <- list(
  Akkermansia_muciniphila      = list(rrna = "Akkermansia_muciniphila_16S.fasta",     genome = "Akkermansia_muciniphila.fasta"),
  Bacteroides_fragilis         = list(rrna = "Bacteroides_fragilis_16S.fasta",         genome = "Bacteroides_fragilis.fasta"),
  Bifidobacterium_adolescentis = list(rrna = "Bifidobacterium adolescentis_16S.fasta", genome = "Bifidobacterium_adolescentis.fasta"),
  Candida_albicans             = list(rrna = "Candida_albicans_18S.fasta",
                                      genome = file.path(cfg$paths$project_root, "ref", "genomes", "Candida_albicans.fasta"),
                                      genome_is_absolute = TRUE),
  Clostridioides_difficile     = list(rrna = "Clostridioides_difficile_16S.fasta",     genome = "Clostridioides_difficile.fasta"),
  Clostridium_perfringens      = list(rrna = "Clostridium_perfringens_16S.fasta",      genome = "Clostridium_perfringens.fasta"),
  Enterococcus_faecalis        = list(rrna = "Enterococcus_faecalis_16S.fasta",        genome = "Enterococcus_faecalis.fasta"),
  Escherichia_coli_B1109       = list(rrna = "Escherichia_coli_B-1109_16S.fasta",      genome = "Escherichia_coli_B1109.fasta"),
  Escherichia_coli_B3008       = list(rrna = "Escherichia_coli_B-3008_16S.fasta",      genome = "Escherichia_coli_B3008.fasta"),
  Escherichia_coli_B766        = list(rrna = "Escherichia_coli_B-766_16S.fasta",       genome = "Escherichia_coli_B766.fasta"),
  Escherichia_coli_JM109       = list(rrna = "Escherichia_coli_JM109_16S.fasta",       genome = "Escherichia_coli_JM109.fasta"),
  Escherichia_coli_B2207       = list(rrna = "Escherichia_oli_B-2207_16S.fasta",       genome = "Escherichia_coli_b2207.fasta"),
  Faecalibacterium_prausnitzii = list(rrna = "Faecalibacterium_prausnitzii_16S.fasta", genome = "Faecalibacterium_prausnitzii.fasta"),
  Fusobacterium_nucleatum      = list(rrna = "Fusobacterium nucleatum_16S.fasta",      genome = "Fusobacterium_nucleatum.fasta"),
  Lactobacillus_fermentum      = list(rrna = "Lactobacillus_fermentum_16S.fasta",      genome = "Lactobacillus_fermentum.fasta"),
  Methanobrevibacter_smithii   = list(rrna = "Methanobrevibacter_smithii_16S.fasta",   genome = "Methanobrevibacter_smithii.fasta"),
  Prevotella_corporis          = list(rrna = "Prevotella_corporis.fasta",              genome = "Prevotella_corporis.fasta"),
  Roseburia_hominis            = list(rrna = "Roseburia_hominis_16S.fasta",            genome = "Roseburia_hominis.fasta"),
  Saccharomyces_cerevisiae     = list(rrna = "Saccharomyces_cerevisiae_18S.fasta",     genome = "Saccharomyces_cerevisiae.fasta"),
  Salmonella_enterica          = list(rrna = "Salmonella_enterica_16S.fasta",          genome = "Salmonella_enterica.fasta"),
  Veillonella_rogosae          = list(rrna = "Veillonella_rogosae_16S.fasta",          genome = "Veillonella_rogosae.fasta")
)

## R8: full rRNA OPERONS via barrnap on the collapsed species genomes (16S-23S-5S
## for bacteria/archaea, 18S-28S-5.8S for fungi) -- broader than the 16S/18S gene
## body recovered from the supplier queries. Coordinates are in the same space as
## ref/genomes/<species>.fasta = zymo_members.fasta. Skipped (with a hint) if
## barrnap is not installed. Plasmid backbones & mobile elements are deferred
## (contigs were collapsed with N-spacers, so header-based plasmid masking is lost;
## mobile elements are the least critical category). [R8 / note A]
barrnap_operons <- function() {
  ## resolve how to invoke barrnap: on PATH, else the zymo_prep conda env via
  ## `conda run` (barrnap is a Perl script needing its own env's modules). [R8]
  invoke <- NULL
  if (nzchar(Sys.which("barrnap"))) {
    invoke <- list(bin = "barrnap", pre = character(0))
  } else {
    conda <- Sys.which("conda"); if (!nzchar(conda)) conda <- path.expand("~/miniconda3/condabin/conda")
    if (file.exists(conda) && dir.exists(path.expand("~/miniconda3/envs/zymo_prep")))
      invoke <- list(bin = conda, pre = c("run", "-n", "zymo_prep", "barrnap"))
  }
  if (is.null(invoke)) {
    message("  barrnap not found -> operon masking skipped (16S/18S gene body only). ",
            "Install: conda create -n zymo_prep --override-channels -c bioconda -c conda-forge barrnap. [R8]")
    return(data.table(tname = character(), tstart = integer(), tend = integer()))
  }
  genome_dir <- file.path(cfg$paths$project_root, "ref", "genomes")
  fastas <- list.files(genome_dir, pattern = "\\.fasta$", full.names = TRUE)
  kingdom <- function(sp) if (grepl("Methanobrevibacter", sp)) "arc"
                          else if (grepl("Candida|Saccharomyces", sp)) "euk" else "bac"
  out <- list()
  for (fa in fastas) {
    sp  <- sub("\\.fasta$", "", basename(fa))
    gff <- tempfile(fileext = ".gff")
    system2(invoke$bin, c(invoke$pre, "--kingdom", kingdom(sp), "--quiet", shQuote(fa)),
            stdout = gff, stderr = FALSE)
    if (!file.exists(gff) || file.size(gff) == 0) next
    g <- tryCatch(fread(gff, sep = "\t", header = FALSE, fill = TRUE), error = function(e) NULL)
    unlink(gff)
    if (is.null(g) || !nrow(g)) next
    g <- g[!startsWith(as.character(V1), "#") & V1 != ""]
    if (!nrow(g)) next
    out[[sp]] <- data.table(tname = as.character(g$V1),
                            tstart = pmax(0L, as.integer(g$V4) - 1L),  # GFF 1-based -> BED 0-based
                            tend = as.integer(g$V5))
  }
  rbindlist(out, use.names = TRUE, fill = TRUE)
}

## note A: LOW-COMPLEXITY regions via dustmasker (BLAST+), run on the collapsed
## members so coordinates match zymo_members.fasta. acclist = >SeqId<tab>from<tab>to
## (0-based inclusive) -> BED half-open.
dustmasker_lowcomplexity <- function() {
  dm <- Sys.which("dustmasker"); members <- cfg$paths$zymo_refs_fasta
  if (!nzchar(dm) || !file.exists(members)) {
    if (!nzchar(dm)) message("  dustmasker not found -> low-complexity masking skipped. [note A]")
    return(data.table(tname = character(), tstart = integer(), tend = integer()))
  }
  acc <- tempfile(fileext = ".acclist")
  system2(dm, c("-in", shQuote(members), "-outfmt", "acclist"), stdout = acc, stderr = FALSE)
  if (!file.exists(acc) || file.size(acc) == 0) return(data.table(tname = character()))
  d <- fread(acc, header = FALSE, sep = "\t"); unlink(acc)
  if (!nrow(d)) return(data.table(tname = character()))
  ## keep only SUBSTANTIAL low-complexity blocks; short (~30 bp) micro-stretches
  ## are handled per-read by the DUST feature and would over-flag long reads. [note A/F]
  d <- d[(V3 - V2 + 1L) >= cfg$params$lowcomplexity_min_len]
  if (!nrow(d)) return(data.table(tname = character()))
  data.table(tname = sub("^>", "", as.character(d$V1)),
             tstart = as.integer(d$V2), tend = as.integer(d$V3) + 1L)
}

## note A: PLASMID backbones. The raw D6331 genomes separate plasmid from
## chromosome as distinct contigs; prepare_zymo_genomes.R merges them (in file
## order, with `spacer` N's between contigs) into one species record. Recompute
## each plasmid contig's offset in that collapsed record and mask the whole contig.
plasmid_regions <- function(spacer = 100L) {
  src_dir <- cfg$paths$zymo_refs_dir
  if (!dir.exists(src_dir)) return(data.table(tname = character(), tstart = integer(), tend = integer()))
  out <- list()
  for (fa in list.files(src_dir, pattern = "\\.(fa|fasta|fna)$", full.names = TRUE)) {
    recs <- read_fasta_records(fa)                         # id, seq, nchar in FILE order
    if (!nrow(recs) || !any(grepl("plasmid", recs$id, ignore.case = TRUE))) next
    starts <- cumsum(c(0L, head(recs$nchar, -1L) + spacer)) # start_i = sum(len[<i]) + i*spacer
    is_pl  <- grepl("plasmid", recs$id, ignore.case = TRUE)
    sp <- canonical_species(sub("\\.(fa|fasta|fna)$", "", basename(fa)))
    out[[fa]] <- data.table(tname = sp, tstart = starts[is_pl],
                            tend = starts[is_pl] + recs$nchar[is_pl])
  }
  rbindlist(out, use.names = TRUE, fill = TRUE)
}

run_prepare_ambiguous_regions <- function(sources = RRNA_SOURCES,
                                          out_path = cfg$params$ambiguous_bed) {
  if (file.exists(out_path)) {
    message(out_path, " already exists -- skipping. Delete it to force a rebuild.")
    return(invisible(out_path))
  }
  mm2 <- need_tool("minimap2")
  all_hits <- list()
  for (species in names(sources)) {
    src <- sources[[species]]
    rrna_path   <- file.path(SSRRNA_DIR, src$rrna)
    ## always the collapsed, species-named genome (matches zymo_members.fasta)
    genome_path <- file.path(GENOMES_DIR, paste0(species, ".fasta"))
    message("== ", species, " ==")
    hits <- tryCatch({
      q <- first_valid_query(rrna_path)
      align_rrna_hits(q$id, q$seq, genome_path, mm2)
    }, error = function(e) {
      message("  SKIPPED (", conditionMessage(e), ")")
      data.table(tname = character())
    })
    message(sprintf("  %d hit(s) passing identity/coverage cutoffs", nrow(hits)))
    if (nrow(hits)) all_hits[[species]] <- hits
  }

  ## R8: add full rRNA OPERONS from barrnap (not just the 16S/18S gene body).
  bar <- barrnap_operons()
  if (nrow(bar)) { message("  barrnap: ", nrow(bar), " rRNA feature(s)"); all_hits[["barrnap"]] <- bar }

  ## note A: low-complexity (dustmasker) + plasmid backbones
  lc <- dustmasker_lowcomplexity()
  if (nrow(lc)) { message("  dustmasker: ", nrow(lc), " low-complexity interval(s)"); all_hits[["dust"]] <- lc }
  pl <- plasmid_regions()
  if (nrow(pl)) { message("  plasmids: ", nrow(pl), " plasmid contig(s)"); all_hits[["plasmid"]] <- pl }

  bed <- rbindlist(all_hits, use.names = TRUE, fill = TRUE)
  bed <- merge_intervals(bed)
  setorder(bed, tname, tstart)

  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  fwrite(bed, out_path, sep = "\t", col.names = FALSE)
  message(sprintf("Wrote %d merged rRNA region(s) -> %s", nrow(bed), out_path))
  invisible(out_path)
}

if (sys.nframe() == 0L) run_prepare_ambiguous_regions()
