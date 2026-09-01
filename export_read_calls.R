#!/usr/bin/env Rscript
# =============================================================================
# export_read_calls.R -- per-read table of model call vs ground truth + taxon.
#
# For a chosen arm + model (default: combined arm, XGBoost; LOEO out-of-fold
# predictions) it joins:
#   * predictions.tsv.gz  -> model score + ground-truth label (y) per read
#   * feature_table.parquet -> the BLAST top-hit taxon that was CALLED
# and emits, per read:
#   read_id, donor, titration_level,
#   called_taxon  (BLAST top-hit species) , called_genus (BLAST top-hit genus),
#   true_species  (minimap2 ground-truth species; empty for negatives),
#   ground_truth  (positive / negative)   , model_score,
#   model_call    (positive / negative at the primary recall operating point),
#   correct       (model_call == ground_truth)
#
# The positive/negative CALL needs a threshold: we use the pipeline's primary
# clinical operating point = the score threshold giving CALL_RECALL (default
# 0.95) recall over the positives (95% of true reads scored >= threshold).
#
# Output: results/<tag>/read_calls_<arm>_<model>.tsv.gz  (+ console preview).
# Env: FIG_RUN_TAG (run), CALL_ARM, CALL_MODEL, CALL_RECALL.
# =============================================================================

suppressPackageStartupMessages({ library(data.table); library(arrow); library(dplyr) })

tag    <- Sys.getenv("FIG_RUN_TAG", "gt_fixed_id0.90_cov0.80")
arm    <- Sys.getenv("CALL_ARM", "combined")
model  <- Sys.getenv("CALL_MODEL", "xgboost")
recall <- as.numeric(Sys.getenv("CALL_RECALL", "0.95"))
res    <- file.path("results", tag)
pred_gz <- file.path(res, "predictions.tsv.gz")
stopifnot(file.exists(pred_gz))

## 1) out-of-fold predictions for the chosen arm/model (LOEO scheme) -----------
##    columns: 1 read_id 2 donor 3 titration_level 4 species 5 y 6 unblock
##             7 arm 8 model 9 scheme 10 fold 11 score
cmd <- sprintf("zcat %s | awk -F'\\t' 'NR==1 || ($7==\"%s\" && $8==\"%s\" && $9==\"LOEO\")'",
               shQuote(pred_gz), arm, model)
pr <- fread(cmd = cmd)
if (!nrow(pr)) stop("no predictions for arm='", arm, "' model='", model, "' (LOEO)")

## 2) the BLAST-called taxon per read (feature table; read only what we need) --
ft <- open_dataset(file.path(res, "feature_table.parquet")) |>
  select(read_id, top_species, top_genus, zymo_tname) |>
  collect() |> as.data.table()

d <- merge(pr[, .(read_id, donor, titration_level, true_species = species,
                  y = as.integer(y), model_score = as.numeric(score))],
           ft, by = "read_id", all.x = TRUE)

## 3) primary operating point: threshold at `recall` recall over the positives -
pos_scores <- d$model_score[d$y == 1L]
thr <- as.numeric(stats::quantile(pos_scores, probs = 1 - recall, names = FALSE, na.rm = TRUE))
d[, `:=`(ground_truth = fifelse(y == 1L, "positive", "negative"),
         model_call   = fifelse(model_score >= thr, "positive", "negative"))]
d[, correct := ground_truth == model_call]
for (cc in c("true_species", "top_species", "top_genus", "zymo_tname"))
  d[get(cc) %in% c("", "NA", "\"\""), (cc) := NA_character_]

setnames(d, c("top_species", "top_genus"), c("called_taxon", "called_genus"))
out <- d[, .(read_id, donor, titration_level, called_taxon, called_genus,
             true_species, ground_truth, model_score, model_call, correct)]
setorder(out, donor, titration_level, -model_score)

## 4) write the full per-read table -------------------------------------------
outfile <- file.path(res, sprintf("read_calls_%s_%s.tsv.gz", arm, model))
fwrite(out, outfile, sep = "\t")

## 5) console preview + summary -----------------------------------------------
cat(sprintf("\nrun=%s  arm=%s  model=%s  scheme=LOEO  reads=%s\n",
            tag, arm, model, format(nrow(out), big.mark = ",")))
cat(sprintf("operating point: threshold @ %.0f%% recall = %.4g  (score >= thr -> called positive)\n",
            100 * recall, thr))
cat("\n== confusion matrix (rows = ground truth, cols = model call) ==\n")
print(addmargins(table(ground_truth = out$ground_truth, model_call = out$model_call)))
tp <- out[ground_truth == "positive" & model_call == "positive", .N]
fp <- out[ground_truth == "negative" & model_call == "positive", .N]
fn <- out[ground_truth == "positive" & model_call == "negative", .N]
cat(sprintf("\nprecision=%.4f  recall=%.4f  accuracy=%.4f\n",
            tp / (tp + fp), tp / (tp + fn), mean(out$correct)))

cat("\n== preview: true positives, false positives, false negatives ==\n")
prev <- rbind(head(out[ground_truth == "positive" & model_call == "positive"], 8),
              head(out[ground_truth == "negative" & model_call == "positive"], 4),
              head(out[ground_truth == "positive" & model_call == "negative"], 4))
prev[, model_score := signif(model_score, 4)]
print(prev[, .(read_id, titration_level, called_taxon, true_species,
               ground_truth, model_score, model_call, correct)], nrows = 20)

cat("\n== most-called taxa among reads the model called POSITIVE ==\n")
print(out[model_call == "positive", .N, by = called_taxon][order(-N)][1:15])

cat("\nFull per-read table written to:\n  ", normalizePath(outfile), "\n", sep = "")
