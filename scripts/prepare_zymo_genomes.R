## =============================================================================
## prepare_zymo_genomes.R  --  collapse fragmented per-species genome FASTAs
## -----------------------------------------------------------------------------
## Some per-species reference genomes ship as many small contigs/unitigs from a
## raw assembler (e.g. Candida_albican.fasta: 13,062 contigs, avg ~1kb) rather
## than curated chromosomes. That blows up ref/zymo_contig2species.tsv to one
## row per contig [OI 1]. merge_fasta_contigs() collapses every sequence in a
## FASTA into a single record, joined by an N-spacer so minimap2 never produces
## a spurious alignment spanning two originally-unrelated contigs.
##
## Idempotent: ensure_merged_fasta() skips the rewrite if `out_path` already
## exists and is already a single-sequence FASTA -- safe to re-run every time.
## =============================================================================

if (!exists("cfg")) {
  .sd <- local({
    a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1])))
    else "/home/tdinse/master_timo/projektarbeit/projektarbeit_ml/scripts"
  })
  source(file.path(.sd, "00_config.R"))
}

## Count sequences (header lines) in a FASTA without a full parse.
n_fasta_sequences <- function(path) {
  if (!file.exists(path)) return(0L)
  sum(startsWith(readLines(path), ">"))
}

## Concatenate every sequence in `in_path` into one record under `header`,
## joined by `spacer_len` N's between originally-separate contigs.
merge_fasta_contigs <- function(in_path, out_path, header, spacer_len = 100L, width = 70L) {
  lines <- readLines(in_path)
  is_header <- startsWith(lines, ">")
  n_seqs <- sum(is_header)
  if (n_seqs == 0L) stop("No FASTA headers found in ", in_path, call. = FALSE)

  seq_group <- cumsum(is_header)[!is_header]
  seqs <- vapply(split(lines[!is_header], seq_group), paste, character(1L), collapse = "")
  merged <- paste(seqs, collapse = strrep("N", spacer_len))

  starts <- seq(1L, nchar(merged), by = width)
  wrapped <- substring(merged, starts, pmin(starts + width - 1L, nchar(merged)))

  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  writeLines(c(paste0(">", header), wrapped), out_path)
  message(sprintf("  merged %d contigs -> 1 sequence (%s), %d bp (incl. %d bp of N-spacers) -> %s",
                  n_seqs, header, nchar(merged), spacer_len * (n_seqs - 1L), out_path))
  invisible(out_path)
}

## Only rebuild if `out_path` isn't already a single-sequence FASTA.
ensure_merged_fasta <- function(in_path, out_path, header, spacer_len = 100L) {
  if (n_fasta_sequences(out_path) == 1L) {
    message("  ", out_path, " already merged (1 sequence) -- skipping.")
    return(invisible(out_path))
  }
  merge_fasta_contigs(in_path, out_path, header = header, spacer_len = spacer_len)
}

## -----------------------------------------------------------------------------
## Driver: build ref/zymo_members.fasta from EVERY genome in cfg$paths$zymo_refs_dir.
## Each species -> one species-named record (multi-contig genomes collapsed with
## N-spacers), all concatenated into the members FASTA. ref/zymo_contig2species.tsv
## is (re)written to match the members exactly, so map_contig_to_species() is a
## trivial identity lookup and the tnames join cleanly to the CoA. [OI 1]
## -----------------------------------------------------------------------------

## D6331 file-stem -> canonical species string used by data/zymo_coa_lot_*.tsv.
canonical_species <- function(stem) {
  s <- sub("^Escherichia_coli_b", "Escherichia_coli_B", stem)  # b2207 -> B2207
  if (s == "Candida_albican") s <- "Candida_albicans"          # supplier typo
  s
}

## Species-level grouping for the sample x taxon / Poisson endpoint (R5): the 5
## near-identical E. coli strains collapse to one species -- they cannot be
## resolved from ONT reads and PrackenDB is species-level. The minimap2 reference
## still keeps the 5 strain records (read-level positive calling); only the
## grouping label collapses.
species_group <- function(s) sub("^Escherichia_coli_.*$", "Escherichia_coli", s)

run_prepare_zymo_genomes <- function(src_dir = cfg$paths$zymo_refs_dir,
                                     members_fasta = cfg$paths$zymo_refs_fasta) {
  if (!dir.exists(src_dir))
    stop("zymo_refs_dir not found: ", src_dir,
         "\n  set cfg$paths$zymo_refs_dir in 00_config.R SECTION 1. [OI 1]", call. = FALSE)
  genome_dir <- file.path(cfg$paths$project_root, "ref", "genomes")
  dir.create(genome_dir, recursive = TRUE, showWarnings = FALSE)

  src <- list.files(src_dir, pattern = "\\.(fa|fasta|fna)$", full.names = TRUE)
  if (!length(src)) stop("no FASTA genomes in ", src_dir, call. = FALSE)

  species_vec <- character(0); per_species <- character(0)
  for (in_path in sort(src)) {
    species  <- canonical_species(sub("\\.(fa|fasta|fna)$", "", basename(in_path)))
    out_path <- file.path(genome_dir, paste0(species, ".fasta"))
    message("== ", species, " ==")
    ensure_merged_fasta(in_path, out_path, header = species)  # 1 species-named record
    species_vec <- c(species_vec, species); per_species <- c(per_species, out_path)
  }

  ## concatenate the one-record-per-species files into the members FASTA
  dir.create(dirname(members_fasta), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(members_fasta)) unlink(members_fasta)
  file.create(members_fasta)
  for (f in per_species) file.append(members_fasta, f)
  message(sprintf("  concatenated %d species -> %s", length(per_species), members_fasta))

  ## (re)write contig2species so tname<->species matches the members exactly.
  ## E. coli strain contigs map to the collapsed 'Escherichia_coli' species. [R5]
  ## (Overwrites any earlier multi-contig map, which would not match these records.)
  c2s <- file.path(cfg$paths$project_root, "ref", "zymo_contig2species.tsv")
  writeLines(c("contig\tspecies", paste(species_vec, species_group(species_vec), sep = "\t")), c2s)
  message("  wrote ", c2s, " (", length(species_vec), " contigs -> ",
          length(unique(species_group(species_vec))), " species)")
  invisible(members_fasta)
}

if (sys.nframe() == 0L) run_prepare_zymo_genomes()
