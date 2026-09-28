# --- RNAModBench path bootstrap (added when this file was deposited) ----------
.rb_find <- function(start, marker = "RNAMOD_BENCH_ROOT") {
  d <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(d, marker))) break
    p <- dirname(d)
    if (p == d) break
    d <- p
  }
  d
}
if (is.null(.RB)) {
  .rb_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
  .rb_here <- if (!is.na(.rb_file) && nzchar(.rb_file)) dirname(normalizePath(.rb_file)) else normalizePath(getwd())
  .RB <- Sys.getenv("RNAMODBENCH_ROOT", unset = .rb_find(.rb_here))
  .XB <- Sys.getenv("RNAMODBENCH_LOCAL", unset = file.path(.RB, "_local"))
}
# --------------------------------------------------------------------------- #

#!/usr/bin/env Rscript
# 63 -- Figure S7, full page: metagene rows A-F (62) + quantitative row G-I (61).
#
# Reviewer R3-9 asks *why* the non-m6A tools over-call on unmodified RNA; the
# answer has to be quantitative, not an eyeball comparison of two density curves.
# The bottom row carries that evidence, with every statistic tested against a
# chromosome-stratified permutation null computed at the same sample size (see
# 61_figS8_tables.py):
#
#   G  WT-versus-unmodified-IVT profile divergence per tool.  Points = the nine
#      independent unit pairs (3 WT x 3 IVT), the filled diamond = their mean,
#      grey band = 2.5-97.5% of the null, dashed line = null median.  A mean on
#      top of the band means the two conditions are positionally indistinguishable
#      given the tool's own call numbers.
#   H  Change of the ncRNA-body share between the conditions (WT - IVT) per unit
#      pair plus the majority consensus; the dashed line is no change.
#   I  Distance of each condition's majority profile to its own candidate-universe
#      background (the testable ncRNA territory of that library), with the null of
#      each statistic: it separates "the tool has a positional preference" from
#      "the preference flips with the modification state".
#
# Designed with the nature-figure checklist in mind: conclusion first (G answers
# R3-9, H/I decompose it), one message per panel, no in-panel numbers or
# annotations, no gridlines, Arial throughout, vector PDF + 300 dpi PNG at the
# printed size (180 mm wide, <= 247 mm tall).
#
# Usage
#   conda run -n guitar_asm --no-capture-output Rscript 63_figS8_page.R
#
# Inputs (read-only)
#   figures/figureS8/tables/figS8_panels.rds        (62)
#   .../tables/s8_profile_distance.tsv, s8_null_summary.tsv            (61)
# Outputs
#   figures/FigureS7_rev.{pdf,png}                the page, panels A-I
#   figures/FigureS7_rev_quant.{pdf,png}          QC block of rows G-I
#   figures/FigureS7_print_preview.png            169 mm 300 dpi print check
#   tables/figS8_page_geometry.tsv                page size + font self-check

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(showtext)
})
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(this_file)), "guitar_lib.R"))
arial_setup()

OUTD <- file.path(.RB, "figures/figureS8")
TABD <- file.path(OUTD, "tables")
FIGD <- file.path(OUTD, "figures")

TOOL_ORDER <- c("CHEUI_m5C", "NanoMUD_psi", "NanoMUD_m1psi", "NanoNm",
                "NanoPsu", "NanoSPA_psU")
LABEL_OF <- c(CHEUI_m5C = "CHEUI-m5C", NanoMUD_psi = "NanoMUD-\u03a8",
              NanoMUD_m1psi = "NanoMUD-m1\u03a8", NanoNm = "NanoNm",
              NanoPsu = "NanoPsu", NanoSPA_psU = "NanoSPA-\u03a8")
COLOUR_OF <- c(WT = "#4b81b8", IVT = "#e8a76b")
XLAB <- unname(LABEL_OF[TOOL_ORDER])
FS <- list(axis = 8, ytitle = 8.5, legend = 7)
TAG_PT <- 9                          # bold panel letter (A / B / C / D)
stopifnot(min(unlist(FS)) >= 7)
ROW_W <- 7.087            # 180 mm, same footprint as the Figure 7 / S6 pages
META_H <- 4.00            # the 2 x 3 metagene block (62)
QUANT_H <- 2.70           # the quantitative row (rotated tool labels need room)
#: explicit decimal breaks, so no axis ever shows a 10^k or 1e-03 style tick.
#: three breaks only: five of them sat 7-8 pt apart on the short log axis and the
#: labels visually merged (user report 2026-09-21)
JSD_BREAKS <- c(0.003, 0.03, 0.3)
JSD_LABELS <- function(v) formatC(v, format = "g", digits = 2)

save_pair_cairo <- function(p, out, width, height, dpi = 300) {
  dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
  showtext::showtext_auto(FALSE)
  ggsave(out, plot = p, width = width, height = height, units = "in",
         device = grDevices::cairo_pdf)
  showtext::showtext_auto()
  ggsave(sub("\\.pdf$", ".png", out), plot = p, width = width, height = height,
         units = "in", dpi = dpi)
  showtext::showtext_auto(FALSE)
  message("wrote ", out, " (+png)")
}

read_tsv <- function(f) read.delim(file.path(TABD, f), sep = "\t",
                                   check.names = FALSE, stringsAsFactors = FALSE)
dist <- read_tsv("s8_profile_distance.tsv")
nulls <- read_tsv("s8_null_summary.tsv")
dist$idx <- match(dist$tool, TOOL_ORDER)
nulls$idx <- match(nulls$tool, TOOL_ORDER)
stopifnot(!any(is.na(dist$idx)), !any(is.na(nulls$idx)))

theme_s7 <- function() {
  theme_classic(base_family = "Arial", base_size = FS$axis) +
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(colour = "black", linewidth = 0.5),
          axis.text = element_text(colour = "black", size = FS$axis),
          axis.title = element_text(size = FS$ytitle),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.box = "vertical",
          legend.box.just = "center",
          legend.key.width = grid::unit(11, "pt"),
          legend.key.height = grid::unit(6, "pt"),
          legend.text = element_text(size = FS$legend),
          legend.margin = margin(0, 0, 0, 0),
          legend.box.spacing = grid::unit(1, "pt"),
          plot.margin = margin(2, 3, 1, 2))
}

x_axis <- function() {
  scale_x_continuous(breaks = seq_along(TOOL_ORDER), labels = XLAB,
                     limits = c(0.4, length(TOOL_ORDER) + 0.6), expand = c(0, 0))
}

theme_rot <- function() {
  
  # standing upright; vjust = 1 keeps the slanted text clear of the axis
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
                                   size = FS$legend))
}

y_jsd <- function() {
  # squish (not drop): a null band whose 2.5% quantile sits below the lower
  # limit must still be drawn, otherwise the band would silently disappear
  scale_y_log10(breaks = JSD_BREAKS, labels = JSD_LABELS,
                limits = c(2e-4, 1), oob = scales::squish)
}

#: shared key for the quantitative row: open circle = one independent unit pair,
#: filled diamond = the tested summary, grey swatch = the permutation null band
KEY_COLOUR <- c("unit pair" = "grey30", "mean of pairs" = "#b3413a",
                "majority consensus" = "#b3413a")
KEY_SHAPE <- c("unit pair" = 1, "mean of pairs" = 18,
               "majority consensus" = 18)

#: white hairline around every null band so that bands of neighbouring tools stay
#: visually separated instead of merging into one grey block
band_rect <- function(band, fill, alpha = 0.75, half = 0.38) {
  geom_rect(data = band, inherit.aes = FALSE,
            aes(xmin = idx - half, xmax = idx + half,
                ymin = null_p2_5, ymax = null_p97_5),
            fill = fill, alpha = alpha, colour = "white", linewidth = 0.25)
}

## ---- G: WT vs unmodified IVT profile divergence ------------------------------
panel_g <- function() {
  d <- dist[dist$comparison == "unit_pair", ]
  band <- nulls[nulls$statistic == "jsd_wt_ivt_unit", ]
  mean_by_tool <- data.frame(idx = seq_along(TOOL_ORDER),
                             jsd_mass = vapply(seq_along(TOOL_ORDER),
                               function(i) mean(d$jsd_mass[d$idx == i]), 0.0))
  ggplot() +
    geom_rect(data = band, inherit.aes = FALSE,
              aes(xmin = idx - 0.38, xmax = idx + 0.38, ymin = null_p2_5,
                  ymax = null_p97_5, fill = "null, 2.5\u201397.5%"),
              alpha = 0.85, colour = "white", linewidth = 0.25) +
    geom_segment(data = band, inherit.aes = FALSE,
                 aes(x = idx - 0.38, xend = idx + 0.38,
                     y = null_median, yend = null_median,
                     linetype = "null median"),
                 linewidth = 0.3, colour = "grey35") +
    geom_point(data = d,
               aes(x = idx, y = jsd_mass, colour = "unit pair",
                   shape = "unit pair"),
               size = 0.9, stroke = 0.30,
               position = position_jitter(width = 0.16, height = 0,
                                          seed = 1)) +
    geom_point(data = mean_by_tool,
               aes(x = idx, y = jsd_mass, colour = "mean of pairs",
                   shape = "mean of pairs"), size = 2.0) +
    scale_fill_manual(name = NULL, values = c("null, 2.5\u201397.5%" = "grey85"),
                      breaks = "null, 2.5\u201397.5%") +
    scale_linetype_manual(name = NULL, values = c("null median" = "dashed"),
                          breaks = "null median") +
    scale_colour_manual(name = NULL, values = KEY_COLOUR,
                        breaks = c("unit pair", "mean of pairs")) +
    scale_shape_manual(name = NULL, values = KEY_SHAPE,
                       breaks = c("unit pair", "mean of pairs")) +
    y_jsd() + x_axis() +
    labs(x = NULL, y = "JSD (bits)") +
    theme_s7() + theme_rot() +
    
    ## described in the caption, only the two drawn series keep a key
    guides(colour = guide_legend(order = 1, nrow = 1),
           shape = "none", linetype = "none", fill = "none")
}

## ---- H: change of the ncRNA-body share ---------------------------------------
panel_h <- function() {
  d <- dist[dist$comparison == "unit_pair", ]
  m <- dist[dist$comparison == "majority", ]
  ggplot() +
    geom_hline(aes(yintercept = 0, linetype = "no change"),
               linewidth = 0.3, colour = "grey45") +
    geom_point(data = d,
               aes(x = idx, y = d_share_ncRNA_body, colour = "unit pair",
                   shape = "unit pair"),
               size = 0.9, stroke = 0.30,
               position = position_jitter(width = 0.16, height = 0, seed = 2)) +
    geom_point(data = m,
               aes(x = idx, y = d_share_ncRNA_body,
                   colour = "majority consensus", shape = "majority consensus"),
               size = 2.0) +
    scale_colour_manual(name = NULL, values = KEY_COLOUR,
                        breaks = c("unit pair", "majority consensus")) +
    scale_shape_manual(name = NULL, values = KEY_SHAPE,
                       breaks = c("unit pair", "majority consensus")) +
    scale_linetype_manual(name = NULL, values = c("no change" = "dashed"),
                          breaks = "no change") +
    x_axis() +
    labs(x = NULL, y = "\u0394 ncRNA-body share") +
    theme_s7() + theme_rot() +
    guides(colour = guide_legend(order = 1, nrow = 1),
           shape = "none", linetype = "none")   # "no change" line -> caption
}

## ---- I: distance to the candidate-universe background ------------------------
panel_i <- function() {
  wt <- nulls[nulls$statistic == "jsd_wt_background", ]
  ivt <- nulls[nulls$statistic == "jsd_ivt_background", ]
  obs <- rbind(
    data.frame(idx = wt$idx, observed = wt$observed, role = "WT",
               stringsAsFactors = FALSE),
    data.frame(idx = ivt$idx, observed = ivt$observed, role = "IVT",
               stringsAsFactors = FALSE))
  obs$x <- obs$idx + ifelse(obs$role == "WT", -0.14, 0.14)
  obs$role <- factor(obs$role, levels = c("WT", "IVT"))
  ggplot() +
    band_rect(transform(wt, idx = idx - 0.14), COLOUR_OF[["WT"]], alpha = 0.18,
              half = 0.13) +
    band_rect(transform(ivt, idx = idx + 0.14), COLOUR_OF[["IVT"]], alpha = 0.18,
              half = 0.13) +
    geom_point(data = obs, aes(x = x, y = observed, colour = role,
                               fill = role, shape = role),
               size = 1.6, stroke = 0.45) +
    scale_colour_manual(values = COLOUR_OF,
                        labels = c("WT vs background", "IVT vs background")) +
    scale_fill_manual(values = COLOUR_OF,
                      labels = c("WT vs background", "IVT vs background")) +
    scale_shape_manual(values = c(WT = 21, IVT = 24),
                       labels = c("WT vs background", "IVT vs background")) +
    y_jsd() + x_axis() +
    labs(x = NULL, y = "JSD to background",
         colour = NULL, fill = NULL, shape = NULL) +
    theme_s7() + theme_rot() +
    guides(colour = guide_legend(order = 1, nrow = 2),
           fill = "none", shape = "none")       # one key, not three copies
}

## ---- assemble the page -------------------------------------------------------
main <- function() {
  panels <- readRDS(file.path(TABD, "figS8_panels.rds"))
  stopifnot(length(panels) >= 1)
  g <- panel_g(); h <- panel_h(); i <- panel_i()
  # panel letters: ONE A for the whole 2 x 3 metagene block, then B/C/D for the
  # quantitative row (user decision 2026-09-21).  Tagged per plot on purpose --
  # patchwork's tag_levels would letter all six metagene panels A-F.
  tag_theme <- theme(plot.tag = element_text(family = "Arial", face = "bold",
                                             size = TAG_PT),
                     plot.tag.position = c(0.02, 0.99))
  panels[[1]] <- panels[[1]] + labs(tag = "A") + tag_theme
  g <- g + labs(tag = "B") + tag_theme
  h <- h + labs(tag = "C") + tag_theme
  i <- i + labs(tag = "D") + tag_theme
  plots <- c(panels, list(g, h, i))
  heights <- c(rep(META_H / 2, ceiling(length(panels) / 3)), QUANT_H)
  page <- patchwork::wrap_plots(plots, ncol = 3, heights = heights)
  dir.create(FIGD, showWarnings = FALSE, recursive = TRUE)
  save_pair_cairo(page, file.path(FIGD, "FigureS7_rev.pdf"), ROW_W,
                  META_H + QUANT_H)
  save_pair_cairo(patchwork::wrap_plots(list(g, h, i), ncol = 3),
                  file.path(FIGD, "FigureS7_rev_quant.pdf"), ROW_W, QUANT_H)
  showtext::showtext_auto(TRUE)
  ggsave(file.path(FIGD, "FigureS7_print_preview.png"), plot = page,
         width = 169 / 25.4, height = (META_H + QUANT_H) * 169 / (ROW_W * 25.4),
         units = "in", dpi = 300)
  showtext::showtext_auto(FALSE)
  geom <- data.frame(
    item = c("panels", "ncol", "width_in", "height_in", "width_mm", "height_mm",
             "font_min_pt", "tags"),
    value = c(length(plots), 3, ROW_W, META_H + QUANT_H, ROW_W * 25.4,
              (META_H + QUANT_H) * 25.4, min(unlist(FS)),
              paste(LETTERS[seq_along(plots)], collapse = "")))
  write.table(geom, file.path(TABD, "figS8_page_geometry.tsv"), sep = "\t",
              quote = FALSE, row.names = FALSE)
  message("done: ", length(plots), " panels, ",
          sprintf("%.1f x %.1f mm", ROW_W * 25.4, (META_H + QUANT_H) * 25.4))
}

main()
