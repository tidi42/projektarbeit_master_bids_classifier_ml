#!/usr/bin/env bash
# =============================================================================
# prepare_kraken2_from_core_nt.sh     [OPTIONAL / SECONDARY analysis]
#   Build a Kraken2 database from the SAME sequences as the BLAST core_nt db.
#
#   This is NOT the primary Kraken2 database. The PRIMARY Kraken2 arm should use
#   a CURATED, species-level DB (PrackenDB preferred, or PlusPFP) -- see
#   00_config.R kraken2_db and OI 9 -- because nt/core_nt redundancy and
#   contamination degrade Kraken2's LCA resolution and its per-species k-mer-count
#   features (k2_kmers_taxon_frac, k2_distinct_minimizers). Standard
#   is unsuitable here because it lacks fungi/protists (Zymo has Candida/Saccharomyces).
#
#   Use THIS script only for a controlled SECONDARY comparison that isolates the
#   ALGORITHM (alignment vs minimizer-LCA) by giving both arms an identical
#   reference universe. The classifier-independent minimap2 ground truth already
#   keeps H1 fair, so this controlled arm is optional. [note D / H1 / H7 / OI 9]
#
# Feasibility on this host (measured 2026-08-03): core_nt = 277 GB, RAM = 502 GB,
#   64 cores, 11 TB free, blastdbcmd + kraken2-build present. A core_nt-derived
#   Kraken2 db is ~400-450 GB (the prior build was 453 GB) -> fits in RAM.
#
# Steps: (1) stream core_nt -> FASTA with taxid in the header (uses core_nt's own
# taxids via taxdb, so no accession2taxid coverage gaps); (2) download NCBI
# taxonomy; (3) add library; (4) build.
# =============================================================================
set -euo pipefail

# ---- configuration (match 00_config.R paths) --------------------------------
CORE_NT="${CORE_NT:-/home/tdinse/core_nt/core_nt}"     # BLAST db prefix
KDB="${KDB:-/home/tdinse/kraken2_core_nt}"             # kraken2 db to build (= cfg$paths$kraken2_db)
THREADS="${THREADS:-32}"                                # of 64 available
SCRATCH="${SCRATCH:-$KDB/library}"                      # holds the extracted FASTA (~250 GB)
FASTA="$SCRATCH/core_nt.kraken.fna"
# Optional hard cap on db size in BYTES (downsamples minimizers if RAM-limited);
# empty = no cap. e.g. MAX_DB_SIZE=$((400*1024*1024*1024))
MAX_DB_SIZE="${MAX_DB_SIZE:-}"
# NO_MASK=1 disables kraken2-build's low-complexity masking (dustmasker).
NO_MASK="${NO_MASK:-0}"

command -v blastdbcmd   >/dev/null || { echo "blastdbcmd not on PATH"; exit 1; }
command -v kraken2-build >/dev/null || { echo "kraken2-build not on PATH"; exit 1; }
mkdir -p "$SCRATCH"

# ---- 1) core_nt -> kraken-ready FASTA (taxid in header) ---------------------
# blastdbcmd '%T\t%s' = taxid <tab> sequence; reformat to >kraken:taxid|<taxid>.
# Skip records with missing/zero taxid. This uses core_nt's embedded taxids
# (taxdb.btd/bti), so it does not depend on NCBI accession2taxid coverage.
if [[ -s "$FASTA" ]]; then
  echo "[1/4] $FASTA already present -- skipping extraction."
else
  echo "[1/4] extracting core_nt -> $FASTA (this is I/O heavy, ~hours) ..."
  blastdbcmd -db "$CORE_NT" -entry all -outfmt $'%T\t%s' \
    | awk -F '\t' 'NF==2 && $1 ~ /^[0-9]+$/ && $1!="0" && length($2)>0 {print ">kraken:taxid|"$1"\n"$2}' \
    > "$FASTA"
  echo "      wrote $(grep -c '^>' "$FASTA") sequences."
fi

# ---- 2) NCBI taxonomy (needs internet) --------------------------------------
echo "[2/4] downloading NCBI taxonomy into $KDB ..."
kraken2-build --download-taxonomy --db "$KDB" --threads "$THREADS"

# ---- 3) add the core_nt library ---------------------------------------------
echo "[3/4] adding core_nt library ..."
kraken2-build --add-to-library "$FASTA" --db "$KDB" $( [[ "$NO_MASK" == "1" ]] && echo --no-masking )

# ---- 4) build ----------------------------------------------------------------
echo "[4/4] building (peak RAM ~ db size; ~hours) ..."
BUILD_ARGS=(--build --db "$KDB" --threads "$THREADS")
[[ -n "$MAX_DB_SIZE" ]] && BUILD_ARGS+=(--max-db-size "$MAX_DB_SIZE")
[[ "$NO_MASK" == "1" ]] && BUILD_ARGS+=(--no-masking)
kraken2-build "${BUILD_ARGS[@]}"

kraken2-build --clean --db "$KDB"   # remove intermediates (keeps *.k2d)
echo "DONE -> $KDB  (set cfg\$paths\$kraken2_db to this; it already is)"
echo "Verify:  kraken2-inspect --db $KDB | head"
