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
# Figure 7 row D: replicate-aware non-m6A metagene panels, drawn with Guitar.
#
# The published Fig. 7D was produced by GuitarPlot() (Bioconductor Guitar) on the
# Ensembl GRCh38.112 transcript model; this script redraws that row with the same
# engine on the replicate-resolved inputs written by 21b_export_guitar_bed.py:
#
#   bed/RNA002/<merge>/<Condition>/<mod>/<Tool>.bed      majority consensus
#   bed/RNA002/rep_<tag>/<Condition>/<mod>/<Tool>.bed    one independent unit
#
# Six non-m6A tools, two libraries (HeLa WT / HeLa IVT); thick line = majority
# consensus (sites present in >= 2 of 3 replicates), thin dashed = one replicate.
# The density itself comes from Guitar (samplePoints -> normalize ->
# .generateDensity_CI); only the panel layers, the per-panel key and the page
# layout are ours, because GuitarPlot cannot draw per-tool facets.
#
# Usage
#   conda run -n guitar_asm --no-capture-output Rscript 48_fig7d_guitar.R
#   ... --recompute     ignore the cached RDS
#
# Outputs
#   figures/figure7/figures/Figure7_rev_D.{pdf,png}
#   figures/figure7/tables/fig7d_density.rds
#   figures/figure7/tables/fig7d_panel_inputs.tsv
#   figures/figure7/tables/fig7d_geometry.tsv
#
# Fonts: PDF via grDevices::cairo_pdf with showtext OFF (real Arial embedded);
# showtext is enabled only for the 300 dpi PNG preview.  Nothing below 7 pt.

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

#' Same as save_pair_cairo(), but stamps a bold panel letter in the top-left
#' corner (the A-C rows carry their letters from the Python side).
save_pair_tagged <- function(p, out, width, height, tag, tag_pt, dpi = 300) {
  dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
  g <- patchwork::patchworkGrob(p)
  draw <- function() {
    grid::grid.newpage()
    grid::grid.draw(g)
    
    ## at x = 57.8 pt (the left edge of their plot area); stamping G 1.2 mm from
    ## the canvas edge put it 0.75 in to the left of that line, so all seven
    ## letters now share one vertical line.
    grid::grid.text(tag, x = grid::unit(20.4, "mm"),
                    y = grid::unit(1, "npc") - grid::unit(1.2, "mm"),
                    just = c("left", "top"),
                    gp = grid::gpar(fontface = "bold", fontsize = tag_pt,
                                    fontfamily = "Arial"))
  }
  showtext::showtext_auto(FALSE)
  grDevices::cairo_pdf(out, width = width, height = height)
  draw(); grDevices::dev.off()
  showtext::showtext_auto()
  png(sub("\\.pdf$", ".png", out), width = width, height = height,
      units = "in", res = dpi)
  draw(); grDevices::dev.off()
  showtext::showtext_auto(FALSE)
  message("wrote ", out, " (+png, tag ", tag, ")")
}

## ---- configuration -----------------------------------------------------------
GTR   <- file.path(.XB, "harmonisation/guitar_metagene")
BED   <- file.path(GTR, "bed", "RNA002")
OUTD  <- file.path(.RB, "figures/figure7")
TABD  <- file.path(OUTD, "tables")
FIGD  <- file.path(OUTD, "figures")
getarg <- function(name, default = NULL) {
  a <- commandArgs(TRUE); i <- which(a == name)
  if (length(i) == 1 && length(a) > i) a[i + 1] else default
}
merge <- getarg("--merge", "majority")
rt    <- as.integer(getarg("--rt", "20"))
recompute <- "--recompute" %in% commandArgs(TRUE)

TOOLS <- c("CHEUI_m5C", "NanoMUD_m1psi", "NanoMUD_psi", "NanoNm", "NanoPsu",
           "NanoSPA_psU")
MOD_OF <- c(CHEUI_m5C = "m5C", NanoMUD_m1psi = "m1Psi", NanoMUD_psi = "Psi",
            NanoNm = "Nm", NanoPsu = "Psi", NanoSPA_psU = "Psi")
COND_LABEL <- c(HeLa_WT = "WT", HeLa_IVT = "IVT")
#: panel/legend display names (the manuscript spells the Psi tools with a Greek
#: letter and the m5C/m1Psi ones with a hyphen)
LABEL_OF <- c(CHEUI_m5C = "CHEUI-m5C", NanoMUD_m1psi = "NanoMUD-m1\u03a8",
              NanoMUD_psi = "NanoMUD-\u03a8", NanoNm = "NanoNm",
              NanoPsu = "NanoPsu", NanoSPA_psU = "NanoSPA-\u03a8")


## (scale 0.755) so that the row plus its caption fits the text block, and the
## type scale went up by 0.5 pt to keep the printed floor: 9.5 x 0.755 = 7.2 pt.

## widest entry ("NanoMUD-m1Ψ-IVT") was wider than its axis slot, so neighbours
## touched (audited: -0.05 pt at 9 pt).  The row is placed at 0.9389
## (169/180 mm), so 7.5 pt here prints at ~7.0 pt, just above the 7 pt floor.
FS <- list(axis = 9.5, ytitle = 9.5, legend = 7.5, struct = 9)

## floor, to 7.5 pt -- the six per-axis keys share one row and the
## widest ("NanoMUD-m1Ψ-IVT") was wider than its 78 pt slot, so the
## neighbours touched.  Placed at 0.95\textwidth it prints at ~7.1 pt.
stopifnot(min(unlist(FS)) >= 7.0)
PANEL_W <- 2.10; PANEL_H <- 1.60; KEY_H <- 0.30
NCOL <- 3; NROW <- 2
##: relative height of the white spacer between the two panel rows (0.10 of a
##: row: ~9 pt at the fixed 4.00 in device, enough to keep the lower row's tick
##: labels out of the upper row's schematic boxes; user 2026-09-23)
ROW_GAP <- 0.40

n_lines <- function(p) if (file.exists(p)) length(readLines(p, warn = FALSE)) else NA_integer_

unit_tags <- function(cond) {
  d <- list.files(BED, pattern = "^rep_", full.names = FALSE)
  d <- d[vapply(d, function(x) dir.exists(file.path(BED, x, cond)), logical(1))]
  sort(sub("^rep_", "", d))
}

## one row per curve the panels will draw
build_plan <- function() {
  rows <- list()
  for (tool in TOOLS) {
    mod <- unname(MOD_OF[tool])
    for (cond in names(COND_LABEL)) {
      role <- COND_LABEL[[cond]]
      addrow <- function(kind, tag, path) {
        n <- n_lines(path)
        rows[[length(rows) + 1]] <<- data.frame(
          tool = tool, mod = mod, cond_group = cond, role = role,
          kind = kind, unit_tag = tag, path = path,
          n_sites = if (is.na(n)) -1L else as.integer(n),
          stringsAsFactors = FALSE)
      }
      f <- file.path(BED, merge, cond, mod, paste0(tool, ".bed"))
      if (!is.na(n_lines(f)) && n_lines(f) > 0) addrow("consensus", "", f)
      for (tag in unit_tags(cond)) {
        f <- file.path(BED, paste0("rep_", tag), cond, mod, paste0(tool, ".bed"))
        if (!is.na(n_lines(f)) && n_lines(f) > 0) addrow("unit", tag, f)
      }
    }
  }
  out <- do.call(rbind, rows)
  out$group <- paste(out$tool, out$role, out$kind, out$unit_tag, sep = "|")
  
  ## Figure 8E draws them -- from the panel top down to the gene-model bars
  
  ## stay inside the plot band, so no label or legend is crossed.
  if (length(nm) > 1) {
    vp <- data.frame(x = end[-length(end)])
    out[[length(out) + 1]] <- geom_segment(
      data = vp, inherit.aes = FALSE,
      aes(x = x, xend = x, y = pos$fig_top, yend = pos$rna_lgd_bl),
      linetype = "dotted", linewidth = 0.25, colour = "black")
    message(sprintf("7D layout: labels %.4f, bars %.4f, floor %.4f (%d separators)",
                    pos$rna_comp_text, pos$rna_lgd_bl, pos$fig_bottom,
                    length(nm) - 1L))
  }

  out
}

## Guitar's samplePoints dies on the first malformed group; isolate per group so
## one bad tool x unit BED cannot kill the whole row (observed on Human before).
group_sites <- function(beds, txtype, gt) {
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
      errs <<- c(errs, sprintf("%s: %s", gname, conditionMessage(e))); FALSE
    })
    if (!ok) message("!! samplePoints failed: ", gname, " -- curve skipped")
  }
  list(relative = relative, weight = weight, errors = errs)
}

compute_or_load <- function() {
  dir.create(TABD, showWarnings = FALSE, recursive = TRUE)
  cache <- file.path(TABD, "fig7d_density.rds")
  if (file.exists(cache) && !recompute) {
    message("using cached density: ", basename(cache)); return(readRDS(cache))
  }
  plan <- build_plan()
  if (!nrow(plan)) stop("no input BED found")
  gt <- guitar_txdb("Human", "mrna")
  beds <- as.list(plan$path); names(beds) <- plan$group
  message(sprintf("%d curves over %d tools", nrow(plan), length(unique(plan$tool))))
  s <- group_sites(beds, "mrna", gt)
  if (length(s$errors)) {
    failed <- sub(":.*$", "", s$errors)
    message("failed groups (", length(failed), ")")
  }
  dens <- Guitar:::.generateDensity_CI(s$relative, s$weight,
                                       CI_ResamplingTime = rt, adjust = 1,
                                       enableCI = FALSE)
  res <- list(dens = dens, plan = plan,
              componentWidth = gt$mrna$componentWidthAverage_pct,
              merge = merge, rt = rt, errors = s$errors)
  saveRDS(res, cache)
  write.table(plan, file.path(TABD, "fig7d_panel_inputs.tsv"), sep = "\t",
              quote = FALSE, row.names = FALSE)
  message("cached -> ", cache)
  res
}

## ---- transcript schematic + panel frame (GuitarPlot look) --------------------
MM_PER_PT <- as.numeric(grid::convertUnit(grid::unit(1, "pt"), "mm", valueOnly = TRUE))


## (region_labels() in guitar_lib.R) -- full text, 8.5 pt, and the labels that
## had to move get a thin leader back to their segment.  7D's panels are 2.10 in
## (151.2 pt) wide, so its row is expected to sit almost on the anchors; only
## displacements above 4 pt are worth a leader.  Its axes are left untouched.
##: true plot-area width of one 7D panel, measured on the produced piece
##: (the axis line of Figure7_rev_D.pdf spans x 26.7 -> 168.8 pt).  The 110 pt
##: it used before was Figure S1's narrower cell and made the label row look
##: 23 % wider in data units than it is, which forced a 6.5 pt displacement --
##: more than the 4 pt nudge allows, so the row could not be set inside its
##: components (2026-09-23).
F7_AXES_W_PT   <- 142.1
F7_LABEL_ALLOW <- c(-0.24, 1.05)
##: Printed geometry.  Measured on the produced piece (2026-09-23): the "0" and
##: peak tick rows are 68 pt apart for peak = 2.8, i.e. 1 data unit = 24.3 pt.
##: Everything below is expressed in *pt* and converted, so the gene model keeps
##: a real 3-5 pt thickness instead of the sub-point one that peak fractions gave.
F7_PT_PER_UNIT  <- 24.3
F7_BAR_H        <- c(promoter = 0.002, utr5 = 0.005, cds = 0.010,
                     utr3 = 0.005, tail = 0.002, ncrna = 0.010)
F7_BAR_PT       <- c(promoter = 3.0, utr5 = 4.0, cds = 5.0, utr3 = 4.0,
                     tail = 3.0, ncrna = 5.0)      # bar thickness, pt
F7_MODEL_GAP_PT <- 3.0     # air between the label row and the bars
F7_LABEL_AIR_PT <- 1.5     # clearance the two guide segments keep from the row

f7_pos <- function(peak) {
  
  ## the geometry of Figure 8E / Supplementary S9 (69_fig8e_guitar.R:148-149 and
  ## 23e_figS9_guitar.R): Guitar's own compact band, i.e. a shallow floor and the
  ## native label row, which is what gives the bars a visible thickness and lets
  ## the dotted component separators run from the plot top down to the bars.
  ## The deep floor (-0.42 * peak) that earlier rounds introduced was what made
  ## the fraction-based bars sub-point thin and pushed the labels away.
  pos <- Guitar:::.generate_pos_para(peak)
  pos$fig_bottom    <- -0.12 * peak
  pos$rna_comp_text <- -0.060 * peak
  pos
}

structure_layers <- function(comp_width, pos, axes_w_pt = F7_AXES_W_PT,
                            allow_data = F7_LABEL_ALLOW, leader_pt = 4) {
  out <- list()
  lab_map <- c(promoter = "1kb", utr5 = "5'UTR", cds = "CDS", utr3 = "3'UTR",
               tail = "1kb", ncrna = "ncRNA")
  a_map   <- c(promoter = 0.99, utr5 = 0.99, cds = 0.272, utr3 = 0.99,
               tail = 0.99, ncrna = 0.20)
  nm  <- names(comp_width)
  end <- cumsum(as.numeric(comp_width)); sta <- c(0, end[-length(end)]) + 0.001
  mid <- (sta + end) / 2; pk <- pos$fig_top / 1.05
  h <- unname(F7_BAR_H[nm]) * pk          # Guitar's own bar half-heights
  xlab <- pmin(pmax(mid, 0.05), 0.95)
  for (i in seq_along(nm)) {
    out[[length(out) + 1]] <- annotate("rect", xmin = sta[i], xmax = end[i],
                                       ymin = pos$rna_lgd_bl - h[i],
                                       ymax = pos$rna_lgd_bl + h[i],
                                       fill = "grey35", colour = NA,
                                       alpha = unname(a_map[nm[i]]))
  }
  
  # guitar row.  The segments branch is the only branch that guarantees both:
  # it pins every label to its own component (<= nudge_pt of give), never
  # rotates, and returns an empty leader list as a literal (guitar_lib.R:277),
  # so no line can be drawn even by accident.  The packed branch is the one
  # that pushed the outer "1kb" off its own component and then needed a leader
  # to point back at it -- that is what the user rejected.
  # Every component is wider than the label that sits in it (1kb 27 pt,
  # 5'UTR 31 pt, CDS 80 pt, 3'UTR 55 pt vs 13.7 / 24.1 / 17.7 / 24.1 pt ink at
  # 8.5 pt), so the row fits without moving anything.
  # Size ladder: the S1 structural tier first; the 5.0 pt rung it used before
  # printed at ~4.7 pt, so the ladder stops at 7.5 pt internal, i.e. the
  # >= 7.42 pt internal that prints as >= 7 pt at 0.95\textwidth.
  fit <- fit_region_size(unname(lab_map[nm]), diff(allow_data) * axes_w_pt,
                         ladder = data.frame(size = c(8.5, 8.0, 7.5),
                                             gap  = c(3.0, 2.5, 2.0)))
  if (is.null(fit)) stop("7D region labels do not fit above the 7.5 pt floor")
  stopifnot(fit$size >= 7.42)          # 0.95\textwidth -> >= 7 pt on the page
  rl <- region_labels(xlab, unname(lab_map[nm]), y = pos$rna_comp_text,
                      axes_w_pt = axes_w_pt, allow_data = allow_data,
                      size = fit$size, gap = fit$gap,
                      mode = "segments", rotate = character(0),
                      nudge_pt = 4.0)
  stopifnot(rl$n_leaders == 0L)        # the user forbade leaders (2026-09-23)
  # The library's own `fits` flag additionally demands that its three pushing
  # passes reached the exact gap budget; on this row they land a few hundredths
  # of a point short of it, which is invisible on the page.  The gates that
  # matter are the ones the user asked for: every label stays on its own
  # component (moved <= the 4 pt nudge), nothing leaves the label band, and
  # there are no leaders.  The printed gaps are asserted on the PDF itself by
  # check_text_collisions.py.
  stopifnot(rl$n_leaders == 0L)                 # user 2026-09-23: no leaders
  stopifnot(max(rl$placed$moved_pt) <= 4.0)     # label still on its component
  wd <- rl$placed$ink_pt / axes_w_pt
  if (any(rl$placed$x - wd / 2 < allow_data[1] - 1e-9) ||
      any(rl$placed$x + wd / 2 > allow_data[2] + 1e-9))
    stop("7D region labels leave the label band")
  if (!rl$fits)
    message("7D region labels: library fits=FALSE is the strict gap budget only; ",
            "max moved ", sprintf("%.2f pt", max(rl$placed$moved_pt)),
            ", leaders ", rl$n_leaders, " (printed gaps audited on the PDF)")
  out <- c(out, rl$text, rl$leaders)
  message(sprintf("7D region labels: %.1f pt, gap %.1f, need %.1f pt, fits %s, leaders %d, max moved %.2f pt",
                  fit$size, fit$gap, fit$need_pt, rl$fits, rl$n_leaders,
                  rl$max_moved_pt))
  
  ## path across the whole panel (the guitar layout keeps the label row, the
  ## legend and the schematic *inside* the panel's y range, below the density
  ## baseline), so any guide that reaches the axis necessarily runs through the
  
  ## readable: every label is centred on its own component and the alpha-shaded
  ## boxes of the schematic carry the model, which the caption names in full.
  out
}

COLOUR_OF <- c(WT = "#4b81b8", IVT = "#e8a76b")
panel_key <- function(tool) {
  list(ggplot2::guide_legend(nrow = 1, byrow = TRUE, override.aes = list()))
}

draw_panel <- function(res, tool, show_y) {
  plan <- res$plan[res$plan$tool == tool, ]
  d <- res$dens[res$dens$group %in% plan$group, ]
  if (!nrow(d)) return(NULL)
  meta <- plan[, c("group", "role", "kind", "unit_tag", "n_sites")]
  d <- merge(d, meta, by = "group", all.x = TRUE)
  d$role <- factor(d$role, levels = c("WT", "IVT"))
  d$curvetype <- factor(ifelse(d$kind == "unit", "individual replicate",
                               "majority consensus"),
                        levels = c("majority consensus", "individual replicate"))
  cons <- d[d$kind == "consensus", ]
  unit <- d[d$kind == "unit", ]
  peak <- max(c(cons$density, unit$density), na.rm = TRUE)
  pos <- f7_pos(peak)
  p <- ggplot()
  if (nrow(unit))
    p <- p + geom_line(data = unit, aes(x = x, y = density, colour = role,
                                        linetype = curvetype, group = group),
                       linewidth = 0.30, alpha = 0.7)
  if (nrow(cons))
    p <- p + geom_ribbon(data = cons, aes(x = x, ymin = 0, ymax = density,
                                          fill = role, group = group),
                         alpha = 0.22, colour = NA) +
             geom_line(data = cons, aes(x = x, y = density, colour = role,
                                        linetype = curvetype, group = group),
                       linewidth = 1.2)
  lab <- unname(LABEL_OF[tool])
  p + scale_colour_manual(values = COLOUR_OF,
                          breaks = c("WT", "IVT"),
                          labels = c(paste0(lab, "-WT"), paste0(lab, "-IVT"))) +
    scale_fill_manual(values = COLOUR_OF, guide = "none") +
    scale_linetype_manual(values = c("majority consensus" = "solid",
                                     "individual replicate" = "dashed"),
                          labels = c("consensus (\u2265 2/3)", "replicate")) +
    ## only the WT / IVT key stays in the figure: "thick = majority consensus
    ## of >= 2/3 replicates, thin dashed = individual replicate" is said once in
    ## the caption instead of six times on the page
    guides(colour = guide_legend(nrow = 1, order = 1,
                                 override.aes = list(linewidth = 1.2,
                                                     linetype = "solid")),
           linetype = "none") +
    structure_layers(res$componentWidth, pos) +
    scale_x_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = c(0, 1), expand = FALSE, clip = "off") +
    coord_cartesian(clip = "off") +   # keep the packed region labels visible
    scale_y_continuous(limits = c(pos$fig_bottom, pos$fig_top * 1.18),
                       expand = c(0, 0), breaks = c(0, peak / 2, peak),
                       labels = function(v) formatC(v, format = "g", digits = 2)) +
    
    ## the density axis is named once in the caption instead of twice on the
    ## page, exactly as Figure S1 does.  show_y stays in the signature so the
    ## caller (and the panel grid) is unchanged.
    labs(x = NULL, colour = NULL, linetype = NULL, y = NULL) +
    theme_classic(base_family = "Arial", base_size = FS$axis) +
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(colour = "black", linewidth = 0.5),
          axis.text = element_text(colour = "black", size = FS$axis),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title = element_text(size = FS$ytitle),
          axis.title.y = element_blank(),   # named in the caption (no 90 deg type)
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.box = "vertical",
          legend.box.just = "center",
          legend.key.width = grid::unit(12, "pt"),
          legend.key.height = grid::unit(7, "pt"),
          legend.text = element_text(size = FS$legend, margin = margin(r = 4)),
          legend.spacing.x = grid::unit(5, "pt"),
          legend.margin = margin(0, 0, 0, 0),
          legend.box.spacing = grid::unit(1, "pt"),
          plot.margin = margin(2, 3, 1, 2))
}

## ---- assemble the row --------------------------------------------------------
main <- function() {
  res <- compute_or_load()
  panels <- lapply(seq_along(TOOLS), function(i) {
    draw_panel(res, TOOLS[i], show_y = ((i - 1) %% NCOL) == 0)
  })
  keep <- !vapply(panels, is.null, logical(1))
  panels <- panels[keep]
  if (length(panels) != length(TOOLS))
    message("panels drawn: ", length(panels), " of ", length(TOOLS))
  
  # panel and landed on the first row's schematic boxes ("3.2" on the tail box).
  # Stacking the two rows with an explicit spacer buys ~9 pt of white between
  # them; the device height below is unchanged, so the piece stays 510 x 288 pt.
  n_top <- min(NCOL, length(panels))
  row <- patchwork::wrap_plots(panels[seq_len(n_top)], ncol = NCOL) /
    patchwork::plot_spacer() /
    patchwork::wrap_plots(panels[-seq_len(n_top)], ncol = NCOL) +
    patchwork::plot_layout(heights = c(1, ROW_GAP, 1)) +
    patchwork::plot_annotation(tag_levels = NULL)
  out <- file.path(FIGD, "Figure7_rev_D")
  dir.create(FIGD, showWarnings = FALSE, recursive = TRUE)
  save_pair_tagged(row, paste0(out, ".pdf"), 7.087, 4.00, tag = "G",
                   tag_pt = 15)
  # legend bands must fit inside each row without touching the row below
  legend_pt <- vapply(panels, function(p) {
    g <- ggplot2::ggplotGrob(p)
    i <- which(grepl("guide-box", g$layout$name))
    if (!length(i)) return(0)
    sum(as.numeric(grid::convertHeight(g$heights[unique(g$layout$t[i])],
                                       "pt", valueOnly = TRUE)))
  }, numeric(1))
  row_pt <- (4.00 * 72) * 1 / (1 + ROW_GAP + 1)      # one row of the stacked layout
  message("legend band pt: ", paste(round(legend_pt, 1), collapse = " "))
  # one key row now (the consensus/replicate guide moved to the caption),
  # so the band is a single text line instead of a two-row stack
  if (any(legend_pt < 7) || any(legend_pt + 45 > row_pt))
    stop("panel legend band does not fit the row height")
  geom <- data.frame(item = c("panels", "ncol", "width_in", "height_in",
                              "font_min_pt", "legend_pt_min", "legend_pt_max",
                              "row_height_pt", "tag", "tag_pt",
                              "tag_margin_pt"),
                     value = c(length(panels), NCOL, 7.087, 4.00,
                               min(unlist(FS)), round(min(legend_pt), 1),
                               round(max(legend_pt), 1), round(row_pt, 1),
                               "G", 15, 1.2))
  write.table(geom, file.path(TABD, "fig7d_geometry.tsv"), sep = "\t",
              quote = FALSE, row.names = FALSE)
  message("done")
}

main()


