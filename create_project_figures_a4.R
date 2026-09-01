#!/usr/bin/env Rscript
# =============================================================================
# create_project_figures_a4.R -- A4-optimised re-format of create_project_figures.R
#
# Same data, same panels, same Prism/Okabe-Ito style as create_project_figures.R,
# but re-laid out for a DIN-A4 portrait page:
#   * canvas 12.5 x 16.25 in  ->  height : width = 1.3 : 1
#     placed at 16 cm text width this is 16.0 x 20.8 cm, leaving ~3.9 cm
#     (7-8 lines) of caption/body text at the bottom of the page.
#   * every text element, point and line is magnified by S = 1.3 (+30 %),
#     so at the same on-page width the lettering is 30 % larger than before.
#
# Output: ./figures/<run_tag>/Figure_*_format.pdf and _format.png
# Run from the project root:  Rscript create_project_figures_a4.R
#   (set FIG_RUN_TAG to render another run, e.g. gt_calculated_id0.92_cov0.50)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

# -----------------------------------------------------------------------------
# page geometry + global magnification
# -----------------------------------------------------------------------------
S     <- 1.3                    # magnify all text / points / lines by 30 %
BASE  <- 12 * S                 # 15.6 pt base font
FIG_W <- 12.5                   # in  (identical to the old landscape width, so
FIG_H <- FIG_W * 1.3            #      the on-page scale factor is unchanged)
gs    <- function(x) x * S      # geom-level size helper (mm for text, pt-ish else)

run_tag <- Sys.getenv("FIG_RUN_TAG", "gt_fixed_id0.90_cov0.80")
res_dir <- file.path("results", run_tag)
out_dir <- file.path("figures", run_tag)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(res_dir)) stop("results dir not found: ", normalizePath(res_dir, mustWork = FALSE))
message("run_tag  : ", run_tag)
message("results  : ", res_dir)
message("figures  : ", out_dir, "   (", FIG_W, " x ", FIG_H, " in, text x", S, ")")

# -----------------------------------------------------------------------------
# Okabe-Ito colourblind-safe palette + display names
# -----------------------------------------------------------------------------
mdl_name   <- c(fixed_threshold = "Fixed threshold", glm = "Logistic (GLM)",
                glmmTMB = "GLMM", ranger_rf = "Random forest", xgboost = "XGBoost")
mdl_levels <- c("Fixed threshold", "Logistic (GLM)", "GLMM", "Random forest", "XGBoost")
pal_model  <- c("Fixed threshold" = "#999999", "Logistic (GLM)" = "#0072B2",
                "GLMM" = "#CC79A7", "Random forest" = "#009E73", "XGBoost" = "#D55E00")

arm_name   <- c(blast_only = "BLAST only", kraken2_only = "Kraken2 only", combined = "Combined")
arm_levels <- c("BLAST only", "Kraken2 only", "Combined")
pal_arm    <- c("BLAST only" = "#0072B2", "Kraken2 only" = "#E69F00", "Combined" = "#009E73")

pal_class  <- c(positive = "#009E73", negative = "#D55E00",
                ambiguous = "#E69F00", indeterminate = "#999999")

lv_ord <- c("c1", "c2", "c3", "c4", "c5", "negative")
lv_lab <- c("c1", "c2", "c3", "c4", "c5", "NC")

# -----------------------------------------------------------------------------
# GraphPad Prism-like theme, scaled by S
# -----------------------------------------------------------------------------
theme_prism2 <- function(base = BASE) {
  theme_classic(base_size = base) +
    theme(
      axis.line        = element_line(colour = "black", linewidth = 0.8 * S),
      axis.ticks       = element_line(colour = "black", linewidth = 0.8 * S),
      axis.ticks.length = unit(4 * S, "pt"),
      axis.title       = element_text(face = "bold", colour = "black"),
      axis.text        = element_text(colour = "black"),
      plot.title       = element_text(face = "bold", size = base + 1),
      plot.subtitle    = element_text(colour = "grey30", size = base - 1.5),
      plot.tag         = element_text(face = "bold", size = base + 5),
      strip.background = element_blank(),
      strip.text       = element_text(face = "bold", colour = "black"),
      legend.position  = "bottom",
      legend.title     = element_text(face = "bold", size = base * 1.10),
      legend.text      = element_text(size = base * 1.05),   # ggplot default is 0.8 * base
      legend.key       = element_blank(),
      legend.key.spacing.x = unit(6, "pt"),
      legend.margin    = margin(2, 2, 2, 2),
      plot.margin      = margin(9, 12, 9, 9)
    )
}
theme_set(theme_prism2())

# Prism convention: the displayed y-minimum sits exactly on the x-axis (no lower gap);
# the 5 % upper headroom keeps top points / error bars from being clipped.
EXP_Y0 <- expansion(mult = c(0, 0.05))

# figure-level (patchwork) title block
ann_theme <- theme(plot.title    = element_text(face = "bold", size = BASE + 3),
                   plot.subtitle = element_text(colour = "grey30", size = BASE))

# wide colourbar: the enlarged tick labels need room, otherwise they collide
cbar <- function(w = 6) guide_colourbar(theme = theme(
  legend.title.position = "top",                 # beside the bar the outer tick label collides with it
  legend.title      = element_text(hjust = 0.5),
  legend.key.width  = unit(w, "cm"),             # widen when the tick labels are long strings
  legend.key.height = unit(0.55, "cm")))

save_fig <- function(p, name, w = FIG_W, h = FIG_H) {
  stem     <- paste0(name, "_format")
  pdf_path <- file.path(out_dir, paste0(stem, ".pdf"))
  okpdf <- tryCatch({ ggsave(pdf_path, p, width = w, height = h, device = cairo_pdf); TRUE },
                    error = function(e) FALSE)
  if (!okpdf) ggsave(pdf_path, p, width = w, height = h)          # fallback if no cairo
  ggsave(file.path(out_dir, paste0(stem, ".png")), p, width = w, height = h, dpi = 300)
  message("  wrote ", stem, ".pdf / .png")
}

rd <- function(f) { p <- file.path(res_dir, f); if (!file.exists(p)) stop("missing ", p); fread(p) }
fmt_k <- function(x) fifelse(x >= 1e6, paste0(formatC(x / 1e6, format = "f", digits = 1), "M"),
                     fifelse(x >= 1e3, paste0(round(x / 1e3), "k"), as.character(x)))

# Heatmap label colour by WCAG luminance of the underlying fill: white on the dark
# end of viridis, black on the light end (white-on-yellow is unreadable).
fill_text_col <- function(x, limits = range(x, na.rm = TRUE), direction = 1,
                          option = "viridis", trans = identity) {
  lim <- trans(limits)
  t   <- scales::rescale(scales::squish(trans(x), lim), from = lim)
  if (direction < 0) t <- 1 - t
  m   <- grDevices::col2rgb(scales::gradient_n_pal(scales::viridis_pal(option = option)(256))(t)) / 255
  lin <- ifelse(m <= 0.03928, m / 12.92, ((m + 0.055) / 1.055)^2.4)
  L   <- 0.2126 * lin[1, ] + 0.7152 * lin[2, ] + 0.0722 * lin[3, ]
  ifelse(L > 0.36, "black", "white")
}

# =============================================================================
# Figure 1: dataset design + ground-truth label structure
# =============================================================================
fig1_ok <- tryCatch({
  lab <- fread(file.path(res_dir, "labels.tsv"),
               select = c("donor", "titration_level", "label"))
  lab[, titration_level := factor(titration_level, levels = lv_ord, labels = lv_lab)]
  lab <- lab[!is.na(titration_level)]

  # (A) design heatmap: reads per donor x level
  heat <- lab[, .(reads = .N), by = .(donor, titration_level)]
  p1a <- ggplot(heat, aes(titration_level, donor, fill = reads)) +
    geom_tile(colour = "white", linewidth = gs(0.8)) +
    geom_text(aes(label = fmt_k(reads),
                  colour = fill_text_col(reads, direction = -1, trans = log10)),
              size = gs(2.9), fontface = "bold") +
    scale_colour_identity() +
    scale_fill_viridis_c(trans = "log10", labels = fmt_k, direction = -1,
                         breaks = c(2e4, 1e5, 5e5), guide = cbar(9)) +
    labs(title = "Dataset design", subtitle = "Classified (non-human) reads\nper donor x titration level",
         x = "Titration level", y = "Donor", fill = "Reads") +
    theme_prism2()

  # (B) class composition per level (absolute reads, stacked, linear y)
  comp <- lab[, .(reads = .N), by = .(titration_level, label)]
  comp[, label := factor(label, levels = names(pal_class))]
  p1b <- ggplot(comp, aes(titration_level, reads, fill = label)) +
    geom_col(width = 0.78) +
    scale_fill_manual(values = pal_class, drop = FALSE) +
    scale_y_continuous(labels = fmt_k, expand = EXP_Y0) +
    labs(title = "Ground-truth label composition", subtitle = "Reads per class at each level",
         x = "Titration level", y = "Reads", fill = "Class") +
    guides(fill = guide_legend(nrow = 2)) +
    theme_prism2()

  # (C) positive fraction per level (mean +/- SD across donors)
  tot  <- lab[, .(tot = .N), by = .(donor, titration_level)]
  posd <- lab[label == "positive", .(pos = .N), by = .(donor, titration_level)]
  fr   <- merge(tot, posd, by = c("donor", "titration_level"), all.x = TRUE)
  fr[is.na(pos), pos := 0][, frac := 100 * pos / tot]
  fr1c <- fr[, .(m = mean(frac), s = sd(frac)), by = titration_level]
  p1c <- ggplot(fr1c, aes(titration_level, m, group = 1)) +
    geom_errorbar(aes(ymin = pmax(m - s, 0), ymax = m + s), width = 0.12, linewidth = gs(0.7)) +
    geom_line(linewidth = gs(0.9), colour = "#0072B2") +
    geom_point(size = gs(3), colour = "#0072B2") +
    scale_y_continuous(expand = EXP_Y0) +
    labs(title = "Positive-read fraction vs titration", subtitle = "Mean +/- SD across 8 donors",
         x = "Titration level", y = "Positive reads (%)") +
    theme_prism2()

  fig1 <- (p1a | p1b) / p1c +
    plot_layout(heights = c(1.45, 1)) +
    plot_annotation(tag_levels = "A",
      title = "Figure 1. Dataset design and ground-truth label structure",
      subtitle = sprintf("Run %s -- true Zymo reads vs false positives in ~99%% human ONT data", run_tag),
      theme = ann_theme)
  save_fig(fig1, "Figure_1_dataset_and_labels")
  TRUE
}, error = function(e) { message("  Figure 1 skipped: ", conditionMessage(e)); FALSE })

# =============================================================================
# shared: read-level per-fold metrics (LOEO, all reads, unblocked included)
# =============================================================================
mrl <- rd("metrics_read_level.tsv")
mrl <- mrl[scheme == "LOEO" & stratum == "all" & truncated_included == TRUE]
mrl[, model_f := factor(mdl_name[model], levels = mdl_levels)]

# =============================================================================
# Figure 2: feature-set arms (BLAST-only vs Kraken2-only vs combined)
# =============================================================================
arms3 <- c("blast_only", "kraken2_only", "combined")
# all four learned models were fitted on every arm; the fixed threshold is the
# baseline and is compared separately in Figure 3.
learners <- c("glm", "glmmTMB", "ranger_rf", "xgboost")
a2 <- mrl[arm %in% arms3 & model %in% learners & recall_target == 0.95,
          .(arm, model, fold, auprc, prec = prec_at_recall)]
a2[, `:=`(arm_f = factor(arm_name[arm], levels = arm_levels),
          model_f = factor(mdl_name[model], levels = mdl_levels))]

s2 <- a2[, .(auprc = mean(auprc), auprc_sd = sd(auprc),
             prec = mean(prec), prec_sd = sd(prec)), by = .(arm_f, model_f)]

p2a <- ggplot(s2, aes(model_f, auprc, colour = arm_f)) +
  geom_errorbar(aes(ymin = auprc - auprc_sd, ymax = pmin(auprc + auprc_sd, 1)),
                width = 0.18, linewidth = gs(0.7), position = position_dodge(0.6)) +
  geom_point(size = gs(3), position = position_dodge(0.6)) +
  scale_colour_manual(values = pal_arm) +
  scale_y_continuous(expand = EXP_Y0) +
  coord_cartesian(ylim = c(0.90, 1.0)) +
  labs(title = "Read-level AUPRC by feature set", subtitle = "Mean +/- SD over 8 LOEO folds",
       x = NULL, y = "AUPRC", colour = "Feature set") +
  theme_prism2()

# delta AUPRC vs BLAST-only, per fold then summarised
w <- dcast(a2, model + fold ~ arm, value.var = "auprc")
w[, `:=`(d_comb = combined - blast_only, d_krak = kraken2_only - blast_only)]
d2 <- melt(w, id.vars = c("model", "fold"), measure.vars = c("d_comb", "d_krak"),
           variable.name = "contrast", value.name = "value")
d2s <- d2[, .(m = mean(value), s = sd(value)), by = .(model, contrast)]
d2s[, `:=`(model_f = factor(mdl_name[model], levels = mdl_levels),
           contrast = factor(contrast, levels = c("d_comb", "d_krak"),
                             labels = c("Combined - BLAST", "Kraken2 - BLAST")))]
p2b <- ggplot(d2s, aes(model_f, m, fill = contrast)) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = gs(0.5)) +
  geom_col(width = 0.7, position = position_dodge(0.75)) +
  geom_errorbar(aes(ymin = m - s, ymax = m + s), width = 0.2, linewidth = gs(0.6),
                position = position_dodge(0.75)) +
  scale_fill_manual(values = c("Combined - BLAST" = "#009E73", "Kraken2 - BLAST" = "#E69F00")) +
  scale_y_continuous(expand = EXP_Y0) +
  labs(title = "Change in AUPRC vs BLAST-only", subtitle = "Positive = better than BLAST-only",
       x = NULL, y = expression(bold(Delta ~ "AUPRC")), fill = NULL) +
  theme_prism2()

p2c <- ggplot(s2, aes(model_f, prec, colour = arm_f)) +
  geom_errorbar(aes(ymin = prec - prec_sd, ymax = pmin(prec + prec_sd, 1)),
                width = 0.18, linewidth = gs(0.7), position = position_dodge(0.6)) +
  geom_point(size = gs(3), position = position_dodge(0.6)) +
  scale_colour_manual(values = pal_arm) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = EXP_Y0) +
  labs(title = "Precision at 95% recall",
       subtitle = "Clinical operating point; higher = fewer false positives",
       x = NULL, y = "Precision @ 95% recall", colour = "Feature set") +
  theme_prism2()

fig2 <- (p2a / p2b / p2c) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A",
    title = "Figure 2. Feature-set comparison: BLAST-only vs Kraken2-only vs combined",
    subtitle = "Does adding Kraken2 features help beyond BLAST? (H1)",
    theme = ann_theme) &
  theme(legend.position = "bottom")
save_fig(fig2, "Figure_2_feature_set_arms")

# =============================================================================
# Figure 3: ML model performance (combined arm) -- the main result
# =============================================================================
c3 <- mrl[arm == "combined" & model %in% names(mdl_name)]
# (A) AUPRC per model, mean +/- SD across folds
a3 <- c3[recall_target == 0.95, .(m = mean(auprc), s = sd(auprc)), by = model_f]
p3a <- ggplot(a3, aes(m, model_f, colour = model_f)) +
  geom_errorbarh(aes(xmin = m - s, xmax = pmin(m + s, 1)), height = 0.18, linewidth = gs(0.7)) +
  geom_point(size = gs(3.4)) +
  scale_colour_manual(values = pal_model, guide = "none") +
  labs(title = "Cross-validated read-level AUPRC",
       subtitle = "Mean +/- SD over 8 LOEO folds (combined arm)",
       x = "AUPRC", y = NULL) +
  theme_prism2()

# (B) precision vs recall target, per model
b3 <- c3[, .(prec = mean(prec_at_recall)), by = .(model_f, recall_target)]
p3b <- ggplot(b3, aes(recall_target, prec, colour = model_f, group = model_f)) +
  geom_hline(yintercept = 0.90, linetype = "dotted", colour = "grey35", linewidth = gs(0.7)) +
  geom_line(linewidth = gs(0.9)) + geom_point(size = gs(2.8)) +
  scale_colour_manual(values = pal_model) +
  scale_x_continuous(breaks = c(0.90, 0.95, 0.99), labels = percent_format(accuracy = 1)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = EXP_Y0) +
  labs(title = "Precision across recall targets",
       subtitle = "Operating points at 90 / 95 / 99% recall;\ndotted line = 90% precision",
       x = "Recall target", y = "Precision", colour = "Model") +
  guides(colour = guide_legend(nrow = 2)) +
  theme_prism2()

# (C) per-fold AUPRC heatmap (model x donor fold)
# 5 decimals: at 3-4 dp the tree-ensemble cells round to "1.000"/"1.0000",
# which reads as a perfect score; the true values are 0.99985-0.99995
h3 <- c3[recall_target == 0.95, .(model_f, fold, auprc)]
FIG3C_LO <- 0.80                       # stated in the Figure 3 legend; keep the two in sync
p3c <- ggplot(h3, aes(factor(fold), model_f, fill = auprc)) +
  geom_tile(colour = "white", linewidth = gs(0.6)) +
  geom_text(aes(label = sprintf("%.5f", auprc),
                colour = fill_text_col(auprc, limits = c(FIG3C_LO, 1))),
            size = gs(2.4), fontface = "bold") +
  scale_colour_identity() +
  scale_fill_viridis_c(option = "viridis", limits = c(FIG3C_LO, 1), oob = scales::squish,
                       breaks = c(0.80, 0.90, 1.00), guide = cbar()) +
  labs(title = "Per-fold AUPRC (LOEO)", subtitle = "Each column = one held-out donor",
       x = "Held-out donor fold", y = NULL, fill = "AUPRC") +
  theme_prism2()

# model legend is collected inside the A|B row so it sits ABOVE panel C
top3 <- (p3a | p3b) +
  plot_layout(guides = "collect", tag_level = "keep") &
  theme(legend.position = "bottom")
fig3 <- (top3 / p3c) +
  plot_layout(heights = c(1.15, 1)) +
  plot_annotation(tag_levels = "A",
    title = "Figure 3. Machine-learning model performance (combined feature set)",
    subtitle = "ML beats the fixed threshold (H2); trees edge out the GLM (H4)",
    theme = ann_theme)
save_fig(fig3, "Figure_3_ml_model_performance")

# =============================================================================
# Figure 4: performance vs titration (H3) + feature ablation (H5)
# =============================================================================
mc <- rd("model_comparison.tsv")
conc <- mc[facet == "concentration" & arm == "combined" & recall_target == 0.95 &
             key %in% c("c1", "c2", "c3", "c4", "c5") & model %in% names(mdl_name)]
conc[, `:=`(level = factor(key, levels = c("c1", "c2", "c3", "c4", "c5")),
            model_f = factor(mdl_name[model], levels = mdl_levels))]

# Models saturate at the top of both panels, so a plain line plot hides Random
# forest under XGBoost/GLMM. Distinct shapes + a small horizontal offset (stated
# in the legend) keep every series readable without changing any value.
pd4       <- position_dodge(width = 0.30)
shp_model <- c("Fixed threshold" = 15, "Logistic (GLM)" = 16, "GLMM" = 17,
               "Random forest" = 18, "XGBoost" = 4)

# (A) fixed threshold vs XGBoost: precision @95% recall across levels
a4 <- conc[model %in% c("fixed_threshold", "xgboost")]
p4a <- ggplot(a4, aes(level, prec_at_recall, colour = model_f, group = model_f)) +
  geom_line(linewidth = gs(0.9), position = pd4) +
  geom_point(aes(shape = model_f), size = gs(3.2), stroke = gs(1.1), position = pd4) +
  scale_colour_manual(values = pal_model, guide = "none") +   # shared legend comes from panel B
  scale_shape_manual(values = shp_model, guide = "none") +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = EXP_Y0) +
  labs(title = "Precision vs titration: model vs baseline",
       subtitle = "Gap widens as input mass drops\n(c1 high -> c5 low) (H3)",
       x = "Titration level", y = "Precision @ 95% recall", colour = "Model") +
  theme_prism2()

# (B) AUPRC vs level, all models
p4b <- ggplot(conc, aes(level, auprc, colour = model_f, group = model_f)) +
  geom_line(linewidth = gs(0.8), position = pd4) +
  geom_point(aes(shape = model_f), size = gs(2.9), stroke = gs(1.1), position = pd4) +
  scale_colour_manual(values = pal_model) +
  scale_shape_manual(values = shp_model) +
  scale_y_continuous(expand = EXP_Y0) +
  labs(title = "AUPRC vs titration level",
       subtitle = "Combined arm, all models; points offset\nhorizontally to separate overlaps",
       x = "Titration level", y = "AUPRC", colour = "Model", shape = "Model") +
  guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2)) +
  theme_prism2()

# (C) feature ablation (XGBoost): FDR @99% recall, lower = better
abl_arms <- c(combined = "Full (combined)",
              combined_minus_H5_key = "- margin + human",
              combined_minus_blast_margin = "- BLAST margin",
              combined_minus_human_competitor = "- human competitor")
abl <- rd("metrics_read_level.tsv")
abl <- abl[scheme == "LOEO" & stratum == "all" & truncated_included == TRUE &
             model == "xgboost" & recall_target == 0.99 & arm %in% names(abl_arms)]
abl[, arm_f := factor(abl_arms[arm], levels = unname(abl_arms))]
abl4 <- abl[, .(m = mean(fdr_at_recall), s = sd(fdr_at_recall)), by = .(arm, arm_f)]
# fold change vs the full feature set -- resolves which block carries the effect
abl4[, ratio := m / abl4[arm == "combined", m]]
p4c <- ggplot(abl4, aes(arm_f, m)) +
  geom_col(width = 0.68, fill = "#0072B2") +
  geom_errorbar(aes(ymin = pmax(m - s, 0), ymax = m + s), width = 0.2, linewidth = gs(0.6)) +
  # n = 8 folds: show every fold rather than a density/violin, which would be
  # a kernel estimate from 8 points
  geom_point(data = abl, aes(arm_f, fdr_at_recall), inherit.aes = FALSE,
             position = position_jitter(width = 0.13, height = 0, seed = 1729),
             size = gs(1.7), shape = 21, fill = "white", colour = "grey15", stroke = gs(0.6)) +
  geom_text(aes(y = m + s, label = sprintf("%.2f%s", ratio, "\u00d7")),
            vjust = -0.7, size = gs(3.1), fontface = "bold", colour = "grey15") +
  scale_y_continuous(labels = percent_format(accuracy = 0.01),
                     expand = expansion(mult = c(0, 0.16))) +
  labs(title = "Feature ablation (XGBoost)",
       subtitle = "FDR @ 99% recall, mean +/- SD with the 8 LOEO folds shown; labels = fold change vs full set (H5)",
       x = NULL, y = "FDR @ 99% recall") +
  theme_prism2() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

# model legend is collected inside the A|B row so it sits ABOVE panel C
top4 <- (p4a | p4b) +
  plot_layout(guides = "collect", tag_level = "keep") &
  theme(legend.position = "bottom")
fig4 <- (top4 / p4c) +
  plot_layout(heights = c(1.2, 1)) +
  plot_annotation(tag_levels = "A",
    title = "Figure 4. Titration dependence and feature-block ablation",
    subtitle = "ML advantage is largest at low abundance (H3); key BLAST/human features matter (H5)",
    theme = ann_theme)
save_fig(fig4, "Figure_4_titration_and_ablation")

# =============================================================================
# Figure 5: calibration + H9 species decomposition + hypothesis effects
# Panel order is deliberate: the two supporting analyses first, the hypothesis
# forest plot last, so the figure closes on the main statistical result.
# =============================================================================
# (C) forest plot of hypothesis effect sizes -- the closing panel
ht <- rd("hypothesis_tests.tsv")
ord <- c("H1", "H2", "H3", "H4", "H5", "H5b", "H6", "H9",
         "H9b_truth", "H9b_classifier", "H10", "H11", "H12")
ht <- ht[id %in% ord & is.finite(median_diff)]
setorder(ht, median_diff)                          # ascending effect size
ht[, `:=`(sig = as.logical(significant),
          fam = ifelse(family == "primary", "primary", "secondary"))]
lev <- ht$id                                       # explicit bottom->top y-axis order
ht[, valid_ci := is.finite(ci_lo) & is.finite(ci_hi) & ci_lo <= median_diff & median_diff <= ci_hi]
ht[, lab := ifelse(abs(median_diff) >= 0.01, sprintf("%+.3f", median_diff),
                   sprintf("%+.1e", median_diff))]
p5a <- ggplot(ht, aes(median_diff, id)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey45") +
  geom_errorbarh(data = ht[valid_ci == TRUE], aes(xmin = ci_lo, xmax = ci_hi),
                 height = 0.25, linewidth = gs(0.6), colour = "grey30") +
  geom_point(aes(colour = sig, shape = fam), size = gs(3)) +
  geom_text(aes(label = lab), vjust = -0.9, size = gs(2.7), colour = "grey20") +
  scale_colour_manual(values = c(`TRUE` = "#009E73", `FALSE` = "#999999"),
                      labels = c(`TRUE` = "significant", `FALSE` = "n.s."), name = NULL) +
  scale_shape_manual(values = c(primary = 16, secondary = 17), name = NULL) +
  scale_x_continuous(trans = pseudo_log_trans(sigma = 1e-4),
                     breaks = c(0, 1e-3, 1e-2, 1e-1, 0.5)) +
  scale_y_discrete(limits = lev) +
  labs(title = "Hypothesis effect sizes",
       subtitle = "Median difference (95% CI); effects <0.001 sit near 0 -- see labels",
       x = "Median difference (pseudo-log scale)", y = NULL) +
  theme_prism2()

# (A) calibration Brier per model (combined arm, LOEO), lower = better
cal <- rd("calibration.tsv")
cal <- unique(cal[scheme == "LOEO" & arm == "combined" & model %in% names(mdl_name)],
              by = c("model", "fold"))
b5 <- cal[, .(m = mean(brier), s = sd(brier)), by = model]
b5[, model_f := factor(mdl_name[model], levels = rev(mdl_levels))]
p5b <- ggplot(b5, aes(m, model_f, colour = model_f)) +
  geom_errorbarh(aes(xmin = pmax(m - s, 0), xmax = m + s), height = 0.18, linewidth = gs(0.7)) +
  geom_point(size = gs(3.4)) +
  scale_colour_manual(values = pal_model, guide = "none") +
  labs(title = "Calibration (Brier score)", subtitle = "Lower = better calibrated probabilities (H12)",
       x = "Brier score", y = NULL) +
  theme_prism2()

# (B) H9 species random-effect decomposition
h9 <- rd("h9_species_sensitivity.tsv")
src_name <- c(none = "GLMM (donor)", truth = "GLMM (donor+species, truth)",
              classifier = "GLMM (donor+species, classifier)")
h9b <- rbind(
  data.table(bar = "GLM", auprc = h9$glm_median_auprc[1]),
  h9[, .(bar = src_name[species_source], auprc = glmm_median_auprc)]
)
bar_lv <- c("GLM", "GLMM (donor)", "GLMM (donor+species, truth)",
            "GLMM (donor+species, classifier)")
h9b[, bar := factor(bar, levels = bar_lv)]
bar_cols <- c("GLM" = "#999999", "GLMM (donor)" = "#56B4E9",
              "GLMM (donor+species, truth)" = "#0072B2",
              "GLMM (donor+species, classifier)" = "#CC79A7")
p5c <- ggplot(h9b, aes(auprc, bar, fill = bar)) +
  geom_col(width = 0.68) +
  geom_text(aes(label = sprintf("%.4f", auprc)), hjust = 1.1, colour = "white",
            fontface = "bold", size = gs(3)) +
  scale_fill_manual(values = bar_cols, guide = "none") +
  coord_cartesian(xlim = c(0.98, 1.0)) +
  labs(title = "H9: donor vs species random effects",
       subtitle = "A species baseline lifts the GLMM (H9b)",
       x = "Median AUPRC (combined arm)", y = NULL) +
  theme_prism2()

fig5 <- (p5b / p5c / p5a) +
  plot_layout(heights = c(0.85, 1, 1.6)) +
  plot_annotation(tag_levels = "A",
    title = "Figure 5. Probability calibration, random-effect decomposition and hypothesis effects",
    subtitle = "Calibration (H12) and H9 variance components, closing with the hypothesis effect sizes and 95% CI",
    theme = ann_theme)
save_fig(fig5, "Figure_5_statistics_and_calibration")

message("\nDone. Re-formatted figures written to: ", normalizePath(out_dir))
