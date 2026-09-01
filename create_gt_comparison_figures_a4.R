#!/usr/bin/env Rscript
# =============================================================================
# create_gt_comparison_figures_a4.R -- A4-optimised re-format of
# create_gt_comparison_figures.R (fixed 0.90/0.80 vs calculated 0.92/0.50).
#
# Same data and panels, re-laid out for a DIN-A4 portrait page:
#   * canvas 12.5 x 16.25 in  ->  height : width = 1.3 : 1
#     placed at 16 cm text width this is 16.0 x 20.8 cm, leaving ~3.9 cm
#     (7-8 lines) of caption/body text at the bottom of the page.
#   * every text element, point and line is magnified by S = 1.3 (+30 %).
#
# Output: ./figures/comparison_fixed_vs_calculated/Figure_C*_format.pdf / .png
# Run from the project root:  Rscript create_gt_comparison_figures_a4.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(scales)
})

# -----------------------------------------------------------------------------
# page geometry + global magnification
# -----------------------------------------------------------------------------
S     <- 1.3
BASE  <- 12 * S
FIG_W <- 12.5
FIG_H <- FIG_W * 1.3
gs    <- function(x) x * S

CMP     <- file.path("results", "comparison_fixed_vs_calculated")
out_dir <- file.path("figures", "comparison_fixed_vs_calculated")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(CMP)) stop("run compare_gt_profiles.R first: ", CMP, " not found")
message("figures  : ", out_dir, "   (", FIG_W, " x ", FIG_H, " in, text x", S, ")")

A_LAB <- "fixed (0.90/0.80)"; B_LAB <- "calculated (0.92/0.50)"
PROF  <- c(A_LAB, B_LAB)
pal_prof <- setNames(c("#0072B2", "#D55E00"), PROF)

rd <- function(f) { p <- file.path(CMP, f); if (!file.exists(p)) stop("missing ", p); fread(p) }
prof_f <- function(x) factor(x, levels = PROF)

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

mdl_name   <- c(fixed_threshold = "Fixed threshold", glm = "Logistic (GLM)",
                glmmTMB = "GLMM", ranger_rf = "Random forest", xgboost = "XGBoost")
mdl_levels <- unname(mdl_name)

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
      plot.margin      = margin(9, 12, 9, 9))
}
theme_set(theme_prism2())

# Prism convention: the displayed y-minimum sits exactly on the x-axis (no lower gap);
# the 5 % upper headroom keeps top points / error bars from being clipped.
EXP_Y0 <- expansion(mult = c(0, 0.05))

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
  ok <- tryCatch({ ggsave(pdf_path, p, width = w, height = h, device = cairo_pdf); TRUE },
                 error = function(e) FALSE)
  if (!ok) ggsave(pdf_path, p, width = w, height = h)
  ggsave(file.path(out_dir, paste0(stem, ".png")), p, width = w, height = h, dpi = 300)
  message("  wrote ", stem, ".pdf / .png")
}

# =============================================================================
# Figure C1: what the ground-truth change actually did
# =============================================================================
sens <- rd("cutoff_sensitivity.tsv")
sens[, `:=`(idf = factor(gt_min_identity), covf = factor(gt_min_coverage),
            pct = 100 * frac_retained)]
sens[, txt := fill_text_col(pct, direction = -1)]
marks <- sens[chosen | apriori][, tag := fifelse(chosen, "calculated", "fixed")]

p1a <- ggplot(sens, aes(covf, idf, fill = pct)) +
  geom_tile(colour = "white", linewidth = gs(0.8)) +
  geom_text(aes(label = sprintf("%.1f", pct), colour = txt), size = gs(2.8),
            fontface = "bold") +
  geom_tile(data = marks, fill = NA, colour = "black", linewidth = gs(1.4)) +
  geom_text(data = marks, aes(label = tag, colour = txt), nudge_y = -0.29, size = gs(2.7),
            fontface = "bold") +
  scale_colour_identity() +
  scale_fill_viridis_c(option = "viridis", direction = -1,
                       breaks = c(85, 90, 95), guide = cbar(9)) +
  labs(title = "Cut-off sensitivity landscape",
       subtitle = "% of confident-Zymo reads retained;\ntarget >= 95%",
       x = "Minimum coverage", y = "Minimum identity", fill = "% retained") +
  theme_prism2()

flip <- rd("label_flip_by_level.tsv")
flip <- flip[label_A %chin% c("positive", "negative") & label_B %chin% c("positive", "negative")]
flip[, transition := fifelse(label_A == "negative" & label_B == "positive",
                             "negative -> positive (coverage relaxed)",
                             "positive -> negative (identity tightened)")]
flip[, signed := fifelse(grepl("^negative", transition), N, -N)]
flip[, level := factor(titration_level, levels = c("c1", "c2", "c3", "c4", "c5", "negative"),
                       labels = c("c1", "c2", "c3", "c4", "c5", "NC"))]
p1b <- ggplot(flip, aes(level, signed, fill = transition)) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = gs(0.5)) +
  geom_col(width = 0.72) +
  scale_fill_manual(values = c("negative -> positive (coverage relaxed)" = "#009E73",
                               "positive -> negative (identity tightened)" = "#CC79A7")) +
  scale_y_continuous(labels = function(x) format(abs(x), big.mark = ",", scientific = FALSE), expand = EXP_Y0) +
  labs(title = "Reads relabelled by the new cut-offs",
       subtitle = "2.82 % of all reads;\n88 % of the churn sits at c1",
       x = "Titration level", y = "Reads relabelled", fill = NULL) +
  guides(fill = guide_legend(nrow = 2)) +
  theme_prism2()

pf <- rd("positive_fraction.tsv")
pf[, `:=`(profile = prof_f(profile),
          level = factor(titration_level, levels = c("c1", "c2", "c3", "c4", "c5", "negative"),
                         labels = c("c1", "c2", "c3", "c4", "c5", "NC")))]
p1c <- ggplot(pf, aes(level, pos_pct, colour = profile, group = profile)) +
  geom_line(linewidth = gs(0.9)) + geom_point(size = gs(3)) +
  scale_colour_manual(values = pal_prof) +
  scale_y_continuous(expand = EXP_Y0) +
  labs(title = "Positive-read prevalence is preserved",
       subtitle = "Pooled across donors; the gradient shifts up by 0.2-2.3 pp",
       x = "Titration level", y = "Positive reads (%)", colour = "GT profile") +
  theme_prism2()

figC1 <- (p1a | p1b) / p1c +
  plot_layout(heights = c(1.3, 1)) +
  plot_annotation(tag_levels = "A",
    title = "Figure C1. What the data-driven ground truth changed",
    subtitle = "Fixed (identity >= 0.90, coverage >= 0.80) vs calculated (identity >= 0.92, coverage >= 0.50)",
    theme = ann_theme)
save_fig(figC1, "Figure_C1_ground_truth_definition")

# =============================================================================
# Figure C2: performance concordance
# =============================================================================
pw <- rd("read_level_wide.tsv")
setnames(pw, c("auprc_calculated", "auprc_fixed", "prec_calculated", "prec_fixed"),
         c("auprc_B", "auprc_A", "prec_B", "prec_A"), skip_absent = TRUE)
cmb <- pw[arm == "combined"]
cmb[, model_f := factor(mdl_name[model], levels = rev(mdl_levels))]

p2a <- ggplot(cmb, aes(y = model_f)) +
  geom_segment(aes(x = auprc_A, xend = auprc_B, yend = model_f),
               colour = "grey55", linewidth = gs(1.1)) +
  geom_point(aes(x = auprc_A, colour = A_LAB), size = gs(4.6)) +
  geom_point(aes(x = auprc_B, colour = B_LAB), size = gs(2.6)) +
  scale_colour_manual(values = pal_prof, breaks = PROF, name = "GT profile") +
  labs(title = "Read-level AUPRC (combined arm)",
       subtitle = "Only the fixed-threshold baseline moves (+0.054)",
       x = "AUPRC", y = NULL) +
  theme_prism2()

fold <- rd("per_fold_auprc.tsv")
fw <- dcast(fold, fold + model ~ profile, value.var = "auprc")
setnames(fw, c(A_LAB, B_LAB), c("A", "B"))
fw[, `:=`(delta = B - A, model_f = factor(mdl_name[model], levels = rev(mdl_levels)))]
p2b <- ggplot(fw, aes(delta, model_f)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey45") +
  geom_point(size = gs(2.6), alpha = 0.75, colour = "#0072B2") +
  stat_summary(fun = median, geom = "point", shape = 124, size = gs(6), colour = "black") +
  labs(title = "Per-fold change in AUPRC",
       subtitle = "calculated - fixed; one point per donor fold,\nbar = median",
       x = "\u0394 AUPRC (calculated - fixed)", y = NULL) +
  theme_prism2()

cal <- rd("calibration.tsv")
cw <- dcast(cal, model ~ profile, value.var = "brier")
setnames(cw, c(A_LAB, B_LAB), c("A", "B"))
cw[, model_f := factor(mdl_name[model], levels = rev(mdl_levels))]
p2c <- ggplot(cw, aes(y = model_f)) +
  geom_segment(aes(x = A, xend = B, yend = model_f), colour = "grey55", linewidth = gs(1.1)) +
  geom_point(aes(x = A, colour = A_LAB), size = gs(4.6)) +
  geom_point(aes(x = B, colour = B_LAB), size = gs(2.6)) +
  scale_colour_manual(values = pal_prof, breaks = PROF, name = "GT profile") +
  labs(title = "Calibration (Brier, lower = better)",
       subtitle = "Learned models unchanged; the baseline improves",
       x = "Brier score", y = NULL) +
  theme_prism2()

figC2 <- (p2a | p2b) / p2c +
  plot_layout(guides = "collect", heights = c(1.15, 1)) +
  plot_annotation(tag_levels = "A",
    title = "Figure C2. Model performance is invariant to the ground-truth definition",
    subtitle = "Every learned model reproduces; only the fixed-threshold comparator shifts",
    theme = ann_theme) &
  theme(legend.position = "bottom")
save_fig(figC2, "Figure_C2_performance_concordance")

# =============================================================================
# Figure C3: titration + ablation
# =============================================================================
titr <- rd("titration.tsv")
titr[, `:=`(profile = prof_f(profile), level = factor(level, levels = paste0("c", 1:5)))]
sel <- titr[model %chin% c("fixed_threshold", "xgboost")]
sel[, model_f := factor(mdl_name[model], levels = mdl_levels)]
p3a <- ggplot(sel, aes(level, prec, colour = model_f, linetype = profile,
                       group = interaction(model, profile))) +
  geom_line(linewidth = gs(0.9)) + geom_point(size = gs(2.6)) +
  scale_colour_manual(values = c("Fixed threshold" = "#999999", "XGBoost" = "#D55E00"), name = "Model") +
  scale_linetype_manual(values = c("solid", "22"), name = "GT profile") +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = EXP_Y0) +
  labs(title = "Precision vs titration under both definitions",
       subtitle = "Curves overlap: the dose-response is a property\nof the data, not the labels",
       x = "Titration level", y = "Precision @ 95% recall") +
  guides(colour = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2)) +
  theme_prism2()

gap <- rd("titration_gap.tsv")
gap[, `:=`(profile = prof_f(profile), level = factor(level, levels = paste0("c", 1:5)))]
p3b <- ggplot(gap, aes(level, gap_xgb_minus_fixed, fill = profile)) +
  geom_col(width = 0.72, position = position_dodge(0.78)) +
  scale_fill_manual(values = pal_prof, name = "GT profile") +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.08))) +
  labs(title = "ML-over-baseline advantage (H3)",
       subtitle = "Grows ~10-13x from c1 to c5\nunder both definitions",
       x = "Titration level", y = "Precision gap (XGBoost - threshold)") +
  guides(fill = guide_legend(nrow = 2)) +
  theme_prism2()

abl <- rd("ablation.tsv")
abl_lab <- c(combined = "Full (combined)",
             combined_minus_H5_key = "- margin + human",
             combined_minus_blast_margin = "- BLAST margin",
             combined_minus_human_competitor = "- human competitor")
abl[, `:=`(profile = prof_f(profile), arm_f = factor(abl_lab[arm], levels = unname(abl_lab)))]
p3c <- ggplot(abl, aes(arm_f, fdr_pct, fill = profile)) +
  geom_col(width = 0.7, position = position_dodge(0.78)) +
  geom_errorbar(aes(ymin = pmax(fdr_pct - sd, 0), ymax = fdr_pct + sd),
                width = 0.2, linewidth = gs(0.6), position = position_dodge(0.78)) +
  scale_fill_manual(values = pal_prof, name = "GT profile") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(title = "Feature ablation (XGBoost)",
       subtitle = "Dropping the BLAST-margin block ~doubles the FDR under both definitions (H5)",
       x = NULL, y = "FDR @ 99% recall (%)") +
  theme_prism2() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

## guides are deliberately NOT collected here: A carries a colour + linetype
## legend, so a single collected row overflows the page width.
figC3 <- (p3a | p3b) / p3c +
  plot_layout(heights = c(1.2, 1)) +
  plot_annotation(tag_levels = "A",
    title = "Figure C3. Dose-response and feature dependence replicate",
    subtitle = "The low-abundance advantage (H3) and the margin-feature dependence (H5) are cut-off independent",
    theme = ann_theme)
save_fig(figC3, "Figure_C3_titration_and_ablation")

# =============================================================================
# Figure C4: hypothesis-level agreement
# =============================================================================
hyp <- rd("hypothesis_concordance.tsv")
h <- hyp[!is.na(effect_A) & !is.na(effect_B)]
ORD <- rev(c("H1","H2","H3","H4","H5","H5b","H6","H9",
             "H9b_truth","H9b_classifier","H10","H11","H12"))
h[, idf := factor(id, levels = ORD)]
pl <- pseudo_log_trans(sigma = 1e-4)

p4a <- ggplot(h, aes(y = idf)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey45") +
  geom_segment(aes(x = effect_A, xend = effect_B, yend = idf), colour = "grey55", linewidth = gs(1.1)) +
  geom_point(aes(x = effect_A, colour = A_LAB), size = gs(3)) +
  geom_point(aes(x = effect_B, colour = B_LAB), size = gs(3)) +
  scale_colour_manual(values = pal_prof, breaks = PROF, name = "GT profile") +
  scale_x_continuous(trans = pl, breaks = c(0, 1e-3, 1e-2, 1e-1, 0.5)) +
  labs(title = "Hypothesis effect sizes under both definitions",
       subtitle = "Every effect keeps its sign; pseudo-log x-axis",
       x = "Effect (median difference)", y = NULL) +
  guides(colour = guide_legend(nrow = 2)) +
  theme_prism2()

h[, involves_baseline := id %chin% c("H2", "H6", "H12")]
p4b <- ggplot(h, aes(abs(effect_A), abs(effect_B), colour = involves_baseline)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey45") +
  geom_point(size = gs(3)) +
  geom_text(aes(label = id), size = gs(2.7), vjust = -1.2, check_overlap = TRUE,
            show.legend = FALSE) +
  scale_colour_manual(values = c(`FALSE` = "#0072B2", `TRUE` = "#D55E00"),
                      labels = c(`FALSE` = "model vs model", `TRUE` = "vs threshold baseline"),
                      name = NULL) +
  scale_x_log10(expand = expansion(mult = 0.12)) +
  scale_y_log10(expand = expansion(mult = 0.12)) + coord_equal() +
  labs(title = "Effect-size concordance",
       subtitle = "All 13 effects lie on the identity line",
       x = "|effect|, fixed profile", y = "|effect|, calculated profile") +
  guides(colour = guide_legend(nrow = 2)) +
  theme_prism2()

## The only systematic shift: contrasts measured AGAINST the fixed-threshold
## baseline lose exactly what the baseline gained (+0.054 AUPRC).
base_gain <- pw[arm == "combined" & model == "fixed_threshold", auprc_B - auprc_A]
h[, grp := fcase(id %chin% c("H2", "H6"), "AUPRC vs threshold baseline",
                 id == "H12",             "calibration (Brier units)",
                 default =                "model vs model")]
p4c <- ggplot(h, aes(reorder(id, delta), delta, fill = grp)) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = gs(0.5)) +
  geom_hline(yintercept = -base_gain, linetype = "dashed", colour = "#D55E00", linewidth = gs(0.7)) +
  geom_col(width = 0.72) +
  scale_fill_manual(values = c("AUPRC vs threshold baseline" = "#D55E00",
                               "calibration (Brier units)" = "#E69F00",
                               "model vs model" = "#0072B2"), name = NULL) +
  scale_y_continuous(expand = EXP_Y0) +
  labs(title = "Absolute change in effect size (calculated - fixed)",
       subtitle = sprintf(paste("Dashed line = -%.3f, the AUPRC the baseline gained.",
                                "H2 and H6 lose exactly that; all other effects move < 0.02"),
                          base_gain),
       x = NULL, y = expression(bold(Delta ~ "effect"))) +
  guides(fill = guide_legend(nrow = 1)) +
  theme_prism2() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

figC4 <- ((p4a | p4b) + plot_layout(widths = c(1.25, 1))) / p4c +
  plot_layout(heights = c(1.35, 1)) +
  plot_annotation(tag_levels = "A",
    title = "Figure C4. Hypothesis-level agreement between ground-truth definitions",
    subtitle = "13 testable hypotheses: 13/13 sign agreement, 12/13 verdict agreement",
    theme = ann_theme)
save_fig(figC4, "Figure_C4_hypothesis_agreement")

message("\nDone. Re-formatted comparison figures written to: ", normalizePath(out_dir))
