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
# Figure 8 panel E -- Guitar metagene band of the three Dorado m6A models.
#
# Built exactly like block A of the rebuilt Supplementary Figure S9
# (01_code/code/sites_v2/scripts/23e_figS9_guitar.R): the density kernel of the
# Bioconductor *Guitar* package (samplePoints -> normalize ->
# .generateDensity_CI) on the majority-consensus call set across the technical
# replicates (guitar_metagene_replicates/bed/RNA004/majority/...), one panel per
# model, WT versus the unmodified IVT control, every panel carrying its own key
# and Guitar's native transcript schematic (1 kb - 5'UTR - CDS - 3'UTR - 1 kb
# grey bars with dotted component separators).
#
# Differences from the S9 script: one row of three panels printed at the final
# Figure 8 size (2.28 x 2.05 in each), the block letter E, model names at the
# house rule weight (plain, not the heavy black facet title of the retired
# draft), and the numbers stay outside the figure (n sites -> legend/table).
#
# Usage:
#   conda run -n guitar_asm --no-capture-output Rscript 69_fig8e_guitar.R [--recompute]
#
# Outputs (this directory only):
#   figures/panels/fig8E_guitar.pdf (+ .png, 300 dpi)   the assembled piece
#   tables/fig8e_guitar_density.rds                      cached Guitar densities
#   tables/fig8e_guitar_panel_inputs.tsv                 one row per drawn curve
#   logs/69_fig8e_guitar.log                             caller-side log

suppressPackageStartupMessages({
  library(Guitar)
  library(GenomicFeatures)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(showtext)
})
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(this_file)), "../../../01_code/code/sites_v2/scripts/guitar_lib.R"))
arial_setup()

OUTD <- file.path(.RB, "figures/figure8")
TABD <- file.path(OUTD, "tables")
FIGD <- file.path(OUTD, "figures", "panels")
dir.create(TABD, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGD, showWarnings = FALSE, recursive = TRUE)
recompute <- "--recompute" %in% commandArgs(TRUE)

MERGE <- "majority"                       # same source as Figure S9 block A

## ---- house style (fonts >= 7 pt; only these constants set the type) ----------
FS <- list(axis = 8.5, ytitle = 9.0, title = 9.0, legend = 8.0, struct = 7.5)
E_TITLE_FACE <- "plain"                   # user decision: no heavy black titles
E_TITLE_PT   <- 8.5                       # model name size (regular weight)
stopifnot(min(unlist(FS)) >= 7)
PANEL_W <- 2.28; PANEL_H <- 2.05          # inches per panel / for the whole band
NCOL <- 3L
MM_PER_PT <- as.numeric(grid::convertUnit(grid::unit(1, "pt"), "mm", valueOnly = TRUE))

MODELS <- data.frame(
  tool  = c("Dorado_hac@v5.0.0_m6A@v1",
            "Dorado_hac@v5.1.0_inosine_m6A",
            "Dorado_sup@v5.0.0_m6A@v1"),
  mod   = c("m6A", "m6A", "m6A"),
  label = c("hac@v5.0.0_m6A", "hac@v5.1.0_inosine+m6A", "sup@v5.0.0_m6A"),
  stringsAsFactors = FALSE)

CONDS <- c(WT = "RNA004_HeLa_WT", IVT = "RNA004_HeLa_IVT")
COND_LABEL <- c(WT = "WT", IVT = "unmodified IVT")

n_lines <- function(p) if (file.exists(p)) length(readLines(p, warn = FALSE)) else NA_integer_

bed_file <- function(mod, tool, group)
  file.path(BED, "RNA004", MERGE, group, mod, paste0(tool, ".bed"))

## ---- one row per drawn curve -------------------------------------------------
plan <- do.call(rbind, lapply(seq_len(nrow(MODELS)), function(i) {
  m <- MODELS[i, ]
  do.call(rbind, lapply(names(CONDS), function(cond) {
    f <- bed_file(m$mod, m$tool, CONDS[[cond]])
    n <- n_lines(f)
    if (is.na(n)) stop("missing BED: ", f)
    data.frame(model = m$tool, mod = m$mod, label = m$label,
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
# Guitar's samplePoints dies on the first malformed group: isolate per group so
# one bad curve cannot take the band down (pattern of 23e_figS9_guitar.R).
e_sites <- function(beds, txtype, gt) {
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

cache <- file.path(TABD, "fig8e_guitar_density.rds")
if (file.exists(cache) && !recompute) {
  message("using cached density: ", basename(cache))
  res <- readRDS(cache)
  stopifnot(identical(sort(res$plan$group), sort(plan$group)), res$merge == MERGE)
} else {
  gt   <- guitar_txdb("Human", "mrna")            # cached; never rebuilt here
  beds <- as.list(plan$path); names(beds) <- plan$group
  message(sprintf("computing %d curves over %d models (human mrna Guitar)",
                  length(beds), nrow(MODELS)))
  t0 <- Sys.time()
  s  <- e_sites(beds, "mrna", gt)
  message(sprintf("site sampling %.1f s",
                  as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  if (length(s$errors)) message("failed groups: ", paste(s$errors, collapse = "; "))
  dens <- Guitar:::.generateDensity_CI(s$relative, s$weight,
                                       CI_ResamplingTime = 20, adjust = 1,
                                       enableCI = FALSE)
  res <- list(dens = dens, plan = plan,
              componentWidth = gt$mrna$componentWidthAverage_pct,
              merge = MERGE)
  saveRDS(res, cache)
  message("cached -> ", cache)
}
write.table(plan[, c("label", "model", "mod", "condition", "cond_label",
                     "n_sites", "path")],
            file.path(TABD, "fig8e_guitar_panel_inputs.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

## ---- transcript schematic under every panel (as in S9 / S1) ------------------
e_pos <- function(peak) {
  pos <- Guitar:::.generate_pos_para(peak)
  pos$fig_bottom    <- -0.12 * peak
  pos$rna_comp_text <- -0.060 * peak
  pos
}
add_structure <- function(p, comp_width, pos) {
  # The panels are 2.28 in wide, so the two 1 kb flanks are narrower than their
  # own label: they stay as grey bars + dotted boundaries and are described in
  # the legend, while the three body segments carry the text.
  lab_map    <- c(promoter = "1kb", utr5 = "5'UTR", cds = "CDS", utr3 = "3'UTR",
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
  # All five components carry their label (1 kb - 5'UTR - CDS - 3'UTR - 1 kb, as
  # in Fig. S9).  The average component widths are 0.218/0.066/0.241/0.257/0.218,
  # so the midpoints are 0.109/0.252/0.405/0.653/0.891: at 7.5 pt the five labels
  # stay clear of each other in a 2.28 in panel (checked on the rendered PDF).
  xlab <- pmin(pmax(mid, 0.045), 0.955)
  for (i in seq_along(nm)) {
    p <- p + annotate("rect", xmin = sta[i], xmax = end[i],
                      ymin = pos$rna_lgd_bl - h[i], ymax = pos$rna_lgd_bl + h[i],
                      fill = "grey35", colour = NA, alpha = unname(alpha_map[nm[i]]))
    p <- p + annotate("text", x = xlab[i], y = pos$rna_comp_text,
                      label = unname(lab_map[nm[i]]), family = "Arial",
                      size = FS$struct * MM_PER_PT, colour = "black")
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
  pos  <- e_pos(peak)

  p <- ggplot(d, aes(x = x, colour = condition, fill = condition)) +
    geom_ribbon(aes(ymin = 0, ymax = density, group = group), alpha = 0.20,
                colour = NA) +
    geom_line(aes(y = density, group = group), linewidth = 0.55) +
    scale_colour_manual(values = cols, limits = names(cols)) +
    scale_fill_manual(values = cols, limits = names(cols)) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(pos$fig_bottom, pos$fig_top * 1.18),
                       expand = c(0, 0)) +
    labs(x = NULL, y = "Density") +
    guides(colour = guide_legend(override.aes = list(linewidth = 0.6, alpha = 1)),
           fill = "none") +
    theme_classic(base_family = "Arial", base_size = FS$axis) +
    theme(panel.grid = element_blank(),
          panel.background = element_rect(fill = "white", colour = NA),
          axis.line = element_line(colour = "black", linewidth = 0.4),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),          # region labels are the x axis
          axis.title.y = element_text(size = FS$ytitle, colour = "black"),
          axis.text.y = element_text(size = FS$axis, colour = "black"),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.title = element_blank(),
          legend.text = element_text(size = FS$legend, colour = "black"),
          legend.key.width = unit(14, "pt"),
          legend.key.height = unit(7, "pt"),
          legend.margin = margin(0, 0, 0, 0, "pt"),
          legend.box.margin = margin(-1, 0, 0, 0, "pt"),
          plot.title = element_text(size = FS$title, face = E_TITLE_FACE,
                                    hjust = 0.5, colour = "black"),
          plot.margin = margin(5, 6, 3, 5, "pt"),
          text = element_text(family = "Arial"))
  # The model name is a plain text layer, never the plot title (ggplot would draw
  # the title in its own heavy style, which the user rejected).
  p <- add_structure(p, res$componentWidth, pos)
  p + annotate("text", x = 0.5, y = pos$fig_top * 1.10, label = m$label,
               family = "Arial", fontface = E_TITLE_FACE,
               size = E_TITLE_PT * MM_PER_PT, colour = "black")
}

## ---- export: real Arial in the vector file, outlines fine in the raster ------
# patchwork's own page tag carries no font family, and cairo would then fall
# back to Noto Sans CJK on this machine (checked: the embedded font list of the
# first build held NotoSansCJKkr for the single tag glyph).  Stamping Arial on
# every grob before rendering keeps the PDF Arial-only.
force_family <- function(x, family = "Arial") {
  if (inherits(x$gp, "gpar") && length(x$gp)) x$gp$fontfamily <- family
  if (!is.null(x$children))
    x$children <- lapply(x$children, force_family, family = family)
  if (!is.null(x$grobs))
    x$grobs <- lapply(x$grobs, force_family, family = family)
  x
}

draw_page <- function(p, tag) {
  # patchworkGrob() measures its text on the *open* device, so it is built after
  # the cairo device is up (otherwise R warns that Arial is missing from the
  # PostScript font database and lays the page out with fallback metrics).
  grid::grid.draw(patchwork::patchworkGrob(p))
  if (!is.null(tag))
    grid::grid.text(tag, x = grid::unit(0.045, "npc"), y = grid::unit(0.985, "npc"),
                    just = c("left", "top"),
                    gp = grid::gpar(fontfamily = "Arial", fontface = "bold",
                                    fontsize = 15))
}

save_pair_cairo <- function(p, out, width, height, dpi = 300, tag = NULL) {
  dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
  showtext::showtext_auto(FALSE)            # PDF: cairo embeds the true Arial
  grDevices::cairo_pdf(out, width = width, height = height, family = "Arial")
  grid::grid.newpage(); draw_page(p, tag); grDevices::dev.off()
  showtext::showtext_auto()                 # PNG: glyph outlines are fine
  grDevices::png(sub("\\.pdf$", ".png", out), width = width, height = height,
                 units = "in", res = dpi, type = "cairo")
  grid::grid.newpage(); draw_page(p, tag); grDevices::dev.off()
  showtext::showtext_auto(FALSE)
  message("wrote ", out, " (+png)")
}

panels <- lapply(seq_len(nrow(MODELS)), draw_panel)
if (any(vapply(panels, is.null, logical(1))))
  stop("panels missing: ",
       paste(MODELS$label[vapply(panels, is.null, logical(1))], collapse = ", "))

# The block letter is drawn by draw_page() with an explicit Arial/bold gpar:
# patchwork's own tag carries no font family and fell back to Noto Sans CJK.
# (The bare plot_annotation() is what turns the element into a patchwork object
# that patchworkGrob() accepts.)
band <- patchwork::wrap_elements(wrap_plots(panels, ncol = NCOL)) +
  patchwork::plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0)))
save_pair_cairo(band, file.path(FIGD, "fig8E_guitar.pdf"),
                PANEL_W * NCOL, PANEL_H, tag = "E")

geom <- data.frame(
  item = c("panels", "width_in", "height_in", "merge", "font_min_pt",
           "n_sites_WT_hac5.0", "n_sites_IVT_hac5.0"),
  value = c(nrow(MODELS), PANEL_W * NCOL, PANEL_H, MERGE, min(unlist(FS)),
            plan$n_sites[plan$group == paste(MODELS$tool[1], "WT", sep = "|")],
            plan$n_sites[plan$group == paste(MODELS$tool[1], "IVT", sep = "|")]))
write.table(geom, file.path(TABD, "fig8e_guitar_geometry.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
message("done ", format(Sys.time()))
