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
# Supplementary Figure S9, redrawn for the revision (reviewer R2-2/E1, R3-9,
# R3-8/E8; replicate-structure statement R3-2/E6).
#
# Why this exists.  The published sup9.pdf ("Dorado Other Modifications -
# mRNA Guitar Plots") shows six Dorado other-modification models but was built
# by $RNAMODBENCH_LOCAL/code/code/extracted/$RNAMODBENCH_LOCAL/raw/result_RNA004/scripts/R/
# dorado_other_mods_guitar_wt_ivt.R, which forced every BED line to "+"
# (minus-strand sites read on the plus-strand coordinate), drew one pooled
# curve per condition and carried no numbers at all, so the legend could only
# say "the distribution of m5C and pseU modifications" - a claim the data do
# not support (the unmodified IVT library shows the same profile).
#
# This script redraws the same six panels from the analysis-ready layer
# (harmonisation/callsets, exported as BED by 21b_export_guitar_bed.py: own
# strand, pad 1 bp, Ensembl chromosome spelling, delivered percent_modified
# >= 90 %) and keeps the density kernel of the Bioconductor Guitar package
# itself (samplePoints -> normalize -> .generateDensity_CI) so the curves stay
# comparable with the published SI.
#
#  * one panel per model, WT (blue) versus the unmodified IVT negative control
#    (orange), each panel with its own key; the six Guitar panels form the
#    single block A (the page letters blocks, not individual panels);
#  * RNA004 HeLa WT and IVT are ONE sequencing unit each - there are no
#    replicate curves to draw, and the legend says so (never imply n > 1);
#  * inosine models are deliberately NOT drawn (too few calls; see legend);
#  * the bottom row carries two blocks: B (left) the threshold-resolved
#    false-positive density of the same six models on the unmodified Curlcake
#    control (tables/figS10_curlcake_scan.tsv) and C (right) the score-validity
#    check on the HeLa libraries (tables/figS10_score_validity.tsv) - the
#    "clearly defined FPR metric" reviewer R3-8/E8 asked for.  Numbers stay in
#    the legend and source tables, never inside the figure.
#
# Usage
#   conda run -n guitar_asm --no-capture-output Rscript 23e_figS10_guitar.R
#   ... --merge majority --rt 20 --recompute --page 16x12
#
# Outputs (own tree, never inside $RNAMODBENCH_LOCAL/submission/)
#   figures/figureS10/figures/FigureS9_rev.{pdf,png}
#   figures/figureS10/tables/figS10_density.rds
#   figures/figureS10/tables/figS10_panel_inputs.tsv
#   figures/figureS10/logs/23e_figS10_guitar.log (caller)

suppressPackageStartupMessages({
  library(Guitar)
  library(GenomicFeatures)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(showtext)
})
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(this_file)), "guitar_lib.R"))
arial_setup()

getarg <- function(name, default = NULL) {
  a <- commandArgs(TRUE); i <- which(a == name)
  if (length(i) == 1 && length(a) > i) a[i + 1] else default
}
merge     <- getarg("--merge", "majority")
rt        <- as.integer(getarg("--rt", "20"))
page_arg  <- getarg("--page", "16x12")
recompute <- "--recompute" %in% commandArgs(TRUE)
OUTD <- getarg("--outdir", file.path(.RB, "figures/figureS10"))
TABD <- file.path(OUTD, "tables"); FIGD <- file.path(OUTD, "figures")
dir.create(TABD, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGD, showWarnings = FALSE, recursive = TRUE)

PG_W <- as.numeric(strsplit(page_arg, "x")[[1]][1])
PG_H <- as.numeric(strsplit(page_arg, "x")[[1]][2])
stopifnot(PG_W >= 8, PG_H >= 6)

## ---- export: PDF with the real embedded Arial, PNG as raster preview --------
save_pair_cairo <- function(p, out, width, height, dpi = 300) {
  dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
  showtext::showtext_auto(FALSE)            # PDF: cairo embeds the true Arial
  ggsave(out, plot = p, width = width, height = height, units = "in",
         device = grDevices::cairo_pdf)
  showtext::showtext_auto()                 # PNG: glyph outlines are fine
  ggsave(sub("\\.pdf$", ".png", out), plot = p, width = width, height = height,
         units = "in", dpi = dpi)
  showtext::showtext_auto(FALSE)
  message("wrote ", out, " (+png)")
}

## ---- house style -------------------------------------------------------------
FS <- list(axis = 15, ytitle = 15, title = 19, legend = 13.5, legendtitle = 13,
           struct = 12, tag = 24, axis_small = 12.5)
stopifnot(min(unlist(FS)) >= 7)             # SI floor, house rule
NCOL   <- 3
PANEL_W <- 5.0; PANEL_H <- 4.9

#: Three block letters on the page (user decision 2026-09-21, replacing the
#: single S9 figure number): A = the six Guitar density panels as ONE block,
#: B = the Curlcake threshold scan (bottom left), C = the HeLa score validity
#: (bottom right).  patchwork letters one page element per tag, so each block
#: is wrapped into a single element below - a nested composition would be
#: tagged plot by plot (observed 2026-09-21: a spurious "H" on the right half
#: of the bottom row).
BLOCK_LETTERS <- c("A", "B", "C")
stopifnot(identical(BLOCK_LETTERS, LETTERS[1:3]))
#: Bottom-left half of the page (block B): version hue pair, colour-blind safe
#: and clearly apart from the blue / orange condition palette that the density
#: panels (and block C) use for WT vs unmodified IVT.
B_V500 <- "#6A3D9A"; B_V510 <- "#1B9E77"
#: A = the 3 x 2 Guitar block, the bottom row = B (left) | C (right)
PAGE_DESIGN  <- "AA\nAA\nBC"
PAGE_HEIGHTS <- c(1, 1, 0.75)
MM_PER_PT <- as.numeric(grid::convertUnit(grid::unit(1, "pt"), "mm",
                                          valueOnly = TRUE))

## ---- the six panels, in drawing order (4 pseU models, then the 2 m5C) -------
MODELS <- data.frame(
  mod   = c("Psi", "Psi", "Psi", "Psi", "m5C", "m5C"),
  tool  = c("Dorado_hac@v5.0.0_pseU@v1_otherMod",
            "Dorado_hac@v5.1.0_pseU@v1_otherMod",
            "Dorado_sup@v5.0.0_pseU@v1_otherMod",
            "Dorado_sup@v5.1.0_pseU@v1_otherMod",
            "Dorado_hac@v5.1.0_m5C@v1_otherMod",
            "Dorado_sup@v5.1.0_m5C@v1_otherMod"),
  label = c("hac@v5.0.0_pseU", "hac@v5.1.0_pseU", "sup@v5.0.0_pseU",
            "sup@v5.1.0_pseU", "hac@v5.1.0_m5C", "sup@v5.1.0_m5C"),
  stringsAsFactors = FALSE)
#: every drawn density panel lives inside block A; `slot` = its position in the
#: 3 x 2 grid (R1C1 .. R2C3, left to right, top to bottom).  The page prints no
#: per-panel letter - the letters A, B, C name the three blocks.
MODELS$block <- BLOCK_LETTERS[1]
MODELS$slot  <- sprintf("R%dC%d",
                        (seq_len(nrow(MODELS)) - 1L) %/% NCOL + 1L,
                        (seq_len(nrow(MODELS)) - 1L) %%  NCOL + 1L)
#: never drawn: too few calls for a density curve (exclusion reason in legend,
#: reviewer R1-6); the inosine+m6A models are counted in Fig. S8A/C
EXCLUDED <- c("Dorado_hac@v5.1.0_inosine_m6A_otherMod",
              "Dorado_sup@v5.1.0_inosine_m6A_otherMod")
stopifnot(!any(MODELS$tool %in% EXCLUDED))

CONDS <- c(WT = "RNA004_HeLa_WT", IVT = "RNA004_HeLa_IVT")
COND_LABEL <- c(WT = "WT", IVT = "unmodified IVT")

n_lines <- function(p) if (file.exists(p)) length(readLines(p, warn = FALSE)) else NA_integer_

bed_file <- function(mod, tool, group) {
  file.path(BED, "RNA004", merge, group, mod, paste0(tool, ".bed"))
}

## ---- one row per drawn curve -------------------------------------------------
plan <- do.call(rbind, lapply(seq_len(nrow(MODELS)), function(i) {
  m <- MODELS[i, ]
  do.call(rbind, lapply(names(CONDS), function(cond) {
    f <- bed_file(m$mod, m$tool, CONDS[[cond]])
    n <- n_lines(f)
    if (is.na(n)) stop("missing BED: ", f)
    data.frame(model = m$tool, mod = m$mod, label = m$label,
               block = m$block, slot = m$slot,
               condition = cond, cond_label = COND_LABEL[[cond]],
               group = paste(m$tool, cond, sep = "|"), path = f, n_sites = n,
               stringsAsFactors = FALSE)
  }))
}))
if (any(plan$n_sites < 30))
  message("low-input curves (<30 sites):\n  ",
          paste(plan$label[plan$n_sites < 30], plan$condition[plan$n_sites < 30],
                plan$n_sites[plan$n_sites < 30], sep = " ", collapse = "\n  "))

## ---- Guitar density (cached) -------------------------------------------------
# Guitar's samplePoints dies on the first malformed group; isolate per group so
# one bad curve cannot take the page down (pattern of 23d_figS1_guitar.R).
s9_sites <- function(beds, txtype, gt) {
  sitesGroup <- Guitar:::.getStGroup(stBedFiles = unname(beds),
                                     stGroupName = names(beds))
  relative <- list(); weight <- list(); errs <- character(0)
  for (gname in names(sitesGroup)) {
    ok <- tryCatch({
      sp <- Guitar::samplePoints(sitesGroup[gname], stSampleNum = 3,
                                 stAmblguity = 5, pltTxType = txtype,
                                 stSampleModle = "Equidistance",
                                 mapFilterTranscript = TRUE, gt)
      nz <- Guitar::normalize(sp, gt, txtype, 1, 1)
      relative[[gname]] <- nz[[1]]; weight[[gname]] <- nz[[2]]
      TRUE
    }, error = function(e) {
      errs <<- c(errs, sprintf("%s: %s", gname, conditionMessage(e)))
      FALSE
    })
    if (!ok) message("!! samplePoints failed: ", gname, " -- curve skipped")
  }
  list(relative = relative, weight = weight, errors = errs)
}

cache <- file.path(TABD, "figS10_density.rds")
if (file.exists(cache) && !recompute) {
  message("using cached density: ", basename(cache))
  res <- readRDS(cache)
  stopifnot(identical(sort(res$plan$group), sort(plan$group)),
            res$merge == merge)
  plan$note <- res$plan$note[match(plan$group, res$plan$group)]
  plan$note[is.na(plan$note)] <- ""
} else {
  gt   <- guitar_txdb("Human", "mrna")          # cached; never rebuilt
  beds <- as.list(plan$path); names(beds) <- plan$group
  message(sprintf("computing %d curves over %d models (human mrna Guitar)",
                  length(beds), nrow(MODELS)))
  t0 <- Sys.time()
  s  <- s9_sites(beds, "mrna", gt)
  message(sprintf("site sampling %.1f s",
                  as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  if (length(s$errors)) {
    failed <- sub(":.*$", "", s$errors)
    plan$note <- ifelse(plan$group %in% failed, "samplePoints-failed", "")
    s$relative <- s$relative[!names(s$relative) %in% failed]
    s$weight   <- s$weight[!names(s$weight) %in% failed]
  } else {
    plan$note <- ""
  }
  # CI_ResamplingTime is irrelevant while enableCI = FALSE (probe in
  # figures/figureS1/logs/figS1_probe.log: density identical between 20 and 1000)
  dens <- Guitar:::.generateDensity_CI(s$relative, s$weight,
                                       CI_ResamplingTime = rt, adjust = 1,
                                       enableCI = FALSE)
  res <- list(dens = dens, plan = plan,
              componentWidth = gt$mrna$componentWidthAverage_pct,
              merge = merge, rt = rt)
  saveRDS(res, cache)
  message("cached -> ", cache)
}

write.table(plan[, c("block", "slot", "label", "model", "mod", "condition",
                     "cond_label", "n_sites", "note", "path")],
            file.path(TABD, "figS10_panel_inputs.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

## ---- transcript schematic under every panel (same as 23d/S1) -----------------
s9_pos <- function(peak) {
  pos <- Guitar:::.generate_pos_para(peak)
  pos$fig_bottom    <- -0.12 * peak
  pos$rna_comp_text <- -0.060 * peak
  pos
}
add_structure <- function(p, comp_width, pos) {
  lab_map   <- c(promoter = "1kb", utr5 = "5'UTR", cds = "CDS", utr3 = "3'UTR",
                 tail = "1kb")
  height_map <- c(promoter = 0.002, utr5 = 0.005, cds = 0.010, utr3 = 0.005,
                  tail = 0.002)
  alpha_map  <- c(promoter = 0.99, utr5 = 0.99, cds = 0.272, utr3 = 0.99,
                  tail = 0.99)
  nm  <- names(comp_width)
  end <- cumsum(as.numeric(comp_width)); sta <- c(0, end[-length(end)]) + 0.001
  mid <- (sta + end) / 2
  pk  <- pos$fig_top / 1.05
  h   <- unname(height_map[nm]) * pk
  xlab <- pmin(pmax(mid, 0.05), 0.95)
  for (i in seq_along(nm)) {
    p <- p + annotate("rect", xmin = sta[i], xmax = end[i],
                      ymin = pos$rna_lgd_bl - h[i], ymax = pos$rna_lgd_bl + h[i],
                      fill = "grey35", colour = NA, alpha = unname(alpha_map[nm[i]]))
    p <- p + annotate("text", x = xlab[i], y = pos$rna_comp_text,
                      label = unname(lab_map[nm[i]]), family = "Arial",
                      size = FS$struct * MM_PER_PT)
  }
  if (length(nm) > 1) {
    vp <- data.frame(x = end[-length(end)])
    p <- p + geom_segment(data = vp, inherit.aes = FALSE,
                          aes(x = x, xend = x, y = pos$fig_top,
                              yend = pos$rna_lgd_bl),
                          linetype = "dotted", linewidth = 0.25, colour = "black")
  }
  p
}

## ---- one panel per model -----------------------------------------------------
draw_panel <- function(i) {
  m  <- MODELS[i, ]
  gr <- paste(m$tool, names(CONDS), sep = "|")
  d  <- res$dens[res$dens$group %in% gr, ]
  meta <- plan[, c("group", "condition", "cond_label", "n_sites")]
  d  <- merge(d, meta, by = "group", all.x = TRUE)
  if (!nrow(d)) { message("!! no density for panel ", m$label); return(NULL) }
  d$condition <- factor(d$cond_label, levels = unname(COND_LABEL))
  cols <- setNames(c(WT_COLOR, CMP_COLOR), unname(COND_LABEL))
  peak <- max(d$density, na.rm = TRUE)
  pos  <- s9_pos(peak)

  p <- ggplot(d, aes(x = x, colour = condition, fill = condition)) +
    geom_ribbon(aes(ymin = 0, ymax = density, group = group), alpha = 0.20,
                colour = NA) +
    geom_line(aes(y = density, group = group), linewidth = 0.75) +
    scale_colour_manual(values = cols, limits = names(cols)) +
    scale_fill_manual(values = cols, limits = names(cols)) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(pos$fig_bottom, pos$fig_top), expand = c(0, 0)) +
    labs(x = NULL, y = "Density", title = m$label) +
    guides(colour = guide_legend(override.aes = list(linewidth = 0.8, alpha = 1)),
           fill = "none") +
    theme_classic(base_family = "Arial", base_size = FS$axis) +
    theme(panel.grid = element_blank(),
          panel.background = element_rect(fill = "white", colour = NA),
          axis.line = element_line(colour = "black", linewidth = 0.3),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),          # region labels are the x axis
          axis.title.y = element_text(size = FS$ytitle, colour = "black"),
          axis.text.y = element_text(size = FS$axis, colour = "black"),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.title = element_blank(),
          legend.text = element_text(size = FS$legend),
          legend.key.width = unit(18, "pt"),
          legend.key.height = unit(7, "pt"),
          legend.margin = margin(0, 0, 0, 0, "pt"),
          legend.box.margin = margin(-2, 0, 0, 0, "pt"),
          plot.title = element_text(size = FS$title, face = "bold",
                                    hjust = 0.5, colour = "black"),
          plot.margin = margin(5, 6, 4, 6, "pt"),
          text = element_text(family = "Arial"))
  add_structure(p, res$componentWidth, pos)
}

## ---- block B (bottom left): false-positive density vs modified-read ---------
# Threshold scan on the unmodified Curlcake control; the numbers come from
# tables/figS10_curlcake_scan.tsv (written by 60_figS10_tables.py out of the
# analysis-ready call-set layer, never recomputed here).  One facet per
# modification type, because the candidate-universe normalisation differs
# between the Psi and the m5C models (4,928 vs 3,963 candidates), so a shared
# axis would mislead.
draw_block_b <- function() {
  f <- file.path(TABD, "figS10_curlcake_scan.tsv")
  if (!file.exists(f)) stop("missing ", f, " -- run 60_figS10_tables.py first")
  d <- read.delim(f, stringsAsFactors = FALSE)
  stopifnot(all(d$label %in% MODELS$label), length(unique(d$threshold_pct)) >= 5)
  d$mod        <- factor(d$mod, levels = c("Psi", "m5C"))
  d$label      <- factor(d$label, levels = MODELS$label)
  d$basecaller <- factor(d$basecaller, levels = c("hac", "sup"))
  d$version    <- factor(d$version, levels = c("v5.0.0", "v5.1.0"))
  FLOOR <- 0.5                       # below one call per 10 kb (0.987)
  d$y    <- pmax(d$fp_per_10kb, FLOOR)
  d$zero <- d$fp_per_10kb == 0
  ggplot(d, aes(x = threshold_pct, y = y, colour = version, linetype = basecaller,
                shape = version, group = label)) +
    geom_line(linewidth = 0.9) +
    geom_point(data = subset(d, !zero), size = 2.4, stroke = 0.6) +
    geom_point(data = subset(d, zero), size = 2.4, shape = 1, stroke = 0.8,
               show.legend = FALSE) +                     # hollow = zero calls
    facet_wrap(~ mod, nrow = 1,
               labeller = as_labeller(c(Psi = "pseU models", m5C = "m5C models"))) +
    scale_y_log10(limits = c(FLOOR * 0.7, NA)) +
    scale_x_continuous(breaks = c(10, 30, 50, 70, 90)) +
    scale_shape_manual(values = c(16, 17), labels = c("v5.0.0", "v5.1.0")) +
    scale_colour_manual(values = c(v5.0.0 = B_V500, v5.1.0 = B_V510),
                        labels = c("v5.0.0", "v5.1.0")) +
    scale_linetype_manual(values = c(hac = "solid", sup = "dashed"),
                          labels = c("hac", "sup")) +
    labs(x = "modified-read threshold (%)", y = "FP per 10 kb") +
    theme_classic(base_family = "Arial", base_size = FS$axis) +
    theme(panel.grid = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = FS$axis, colour = "black"),
          axis.line = element_line(colour = "black", linewidth = 0.3),
          axis.text = element_text(size = FS$axis, colour = "black"),
          axis.title = element_text(size = FS$axis, colour = "black"),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.title = element_blank(),
          legend.text = element_text(size = FS$legend),
          legend.key.width = unit(16, "pt"),
          legend.key.height = unit(6, "pt"),
          legend.margin = margin(0, 0, 0, 0, "pt"),
          legend.box.margin = margin(-2, 0, 0, 0, "pt"),
          plot.margin = margin(4, 6, 4, 6, "pt"),
          text = element_text(family = "Arial"))
}

#' Block C (bottom right): does the model's own reported modified fraction
#' separate the wild-type library from the unmodified control?  Empirical
#' cumulative distributions of the per-call scores (WT blue, unmodified IVT
#' orange), one facet per model; the Mann-Whitney AUC lives in the legend /
#' source table.
draw_block_c <- function() {
  f <- file.path(TABD, "figS10_score_validity.tsv")
  if (!file.exists(f)) stop("missing ", f, " -- run 60_figS10_tables.py first")
  d <- read.delim(f, stringsAsFactors = FALSE)
  d$label     <- factor(d$label, levels = MODELS$label)   # drawn order
  d$condition <- factor(d$condition, levels = c("WT", "IVT"))
  cols <- setNames(c(WT_COLOR, CMP_COLOR), c("WT", "IVT"))
  ggplot(d, aes(x = score, colour = condition)) +
    stat_ecdf(geom = "step", linewidth = 0.8, pad = FALSE) +
    facet_wrap(~ label, ncol = 3) +
    scale_colour_manual(values = cols, labels = c("WT", "unmodified IVT")) +
    scale_x_continuous(limits = c(0.90, 1.00), breaks = c(0.90, 0.95, 1.00)) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    labs(x = "reported modified fraction", y = "cumulative fraction") +
    theme_classic(base_family = "Arial", base_size = FS$axis_small) +
    theme(panel.grid = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = FS$axis_small, colour = "black"),
          axis.line = element_line(colour = "black", linewidth = 0.3),
          axis.text = element_text(size = FS$axis_small, colour = "black"),
          axis.title = element_text(size = FS$axis_small, colour = "black"),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.title = element_blank(),
          legend.text = element_text(size = FS$legend),
          legend.key.width = unit(16, "pt"),
          legend.key.height = unit(6, "pt"),
          legend.margin = margin(0, 0, 0, 0, "pt"),
          legend.box.margin = margin(-2, 0, 0, 0, "pt"),
          plot.margin = margin(4, 6, 4, 6, "pt"),
          text = element_text(family = "Arial"))
}

panels <- lapply(seq_len(nrow(MODELS)), draw_panel)
if (any(vapply(panels, is.null, logical(1))))
  stop("panels missing: ",
       paste(MODELS$label[vapply(panels, is.null, logical(1))], collapse = ", "))
# The three taggable blocks.  wrap_elements() turns each one into a SINGLE page
# element, which is what makes plot_annotation() print one letter per block:
# A = the six Guitar density panels (3 x 2, internally untagged), B = the
# Curlcake threshold scan (bottom left), C = the HeLa score validity (bottom
# right).  A nested patchwork would be tagged plot by plot instead (observed
# 2026-09-21: a spurious "H" on the right half of the bottom row).
blk_a <- wrap_elements(wrap_plots(panels, design = "ABC\nDEF"))
blk_b <- wrap_elements(draw_block_b())          # bottom left: threshold scan
blk_c <- wrap_elements(draw_block_c())          # bottom right: score validity

page <- wrap_plots(list(blk_a, blk_b, blk_c), design = PAGE_DESIGN) +
  plot_layout(heights = PAGE_HEIGHTS) +
  plot_annotation(tag_levels = "A") &           # A, B, C in reading order
  theme(plot.tag = element_text(family = "Arial", face = "bold", size = FS$tag),
        plot.tag.position = c(0.010, 0.985))

save_pair_cairo(page, file.path(FIGD, "FigureS9_rev.pdf"), width = PG_W,
                height = PG_H)
message("page: ", PG_W, " x ", PG_H, " in (", PG_W * 72, " x ", PG_H * 72,
        " pt), block letters ", paste(BLOCK_LETTERS, collapse = "/"),
        " and no figure number: ", nrow(MODELS),
        " metagene panels in block A + bottom row B (threshold scan) | ",
        "C (score validity)")
message("done ", format(Sys.time()))
