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
# Figure S8, metagene rows (A-F): replicate-aware non-m6A ncRNA panels.
#
# The published Fig. S8 is a 2 x 3 grid of Guitar density profiles, one panel per
# non-m6A tool, WT versus the unmodified IVT negative control, on the ncRNA axis
# "1kb | ncRNA | 1kb".  Its inputs carried no replicate structure (the legacy
# Guitar BEDs were an undocumented union of the three HeLa replicates) and the
# panel keys still used the legacy tool names -- reviewer points R3-2/E6 and
# R1-7.  This script redraws the grid with the same engine on the
# replicate-resolved inputs written by 21b_export_guitar_bed.py:
#
#   bed/RNA002/<merge>/<Condition>/<mod>/<Tool>.bed      majority consensus (>= 2/3)
#   bed/RNA002/rep_<tag>/<Condition>/<mod>/<Tool>.bed    one independent unit
#
# Six non-m6A tools, two libraries (HeLa WT / HeLa IVT); thick filled line =
# majority consensus, thin dashed = one independent sequencing unit.  The density
# itself comes from Guitar (samplePoints -> normalize -> .generateDensity_CI) on
# the Ensembl GRCh38.112 ncRNA model (never GENCODE); only the panel layers, the
# per-panel key and the page layout are ours, because GuitarPlot cannot draw
# per-tool facets or distinguish consensus from replicate curves.
#
# Usage
#   conda run -n guitar_asm --no-capture-output Rscript 62_figS8_guitar.R
#   ... --recompute     ignore the cached RDS
#   ... --rt 20         CI resampling (switched off, kept for provenance)
#
# Outputs (figures/figureS8/)
#   figures/FigureS8_rev_metagene.{pdf,png}   rows A-F (QC block, no tags)
#   tables/figS8_density.rds                  cached Guitar density + plan
#   tables/figS8_panels.rds                   ggplot objects for 63_figS8_page.R
#   tables/figS8_panel_inputs.tsv             one row per drawn BED (n_sites)
#   tables/figS8_guitar_sites.tsv             Guitar's per-site relative coordinate
#                                             + weight (cross-check vs 61)
#   tables/figS8_guitar_geometry.tsv          component widths, panel size, fonts
#
# Fonts: PDF via grDevices::cairo_pdf with showtext OFF (real Arial embedded);
# showtext is switched on only for the 300 dpi PNG preview.  Nothing below 7 pt.

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

## ---- configuration -----------------------------------------------------------
OUTD <- file.path(.RB, "figures/figureS8")
TABD <- file.path(OUTD, "tables")
FIGD <- file.path(OUTD, "figures")
getarg <- function(name, default = NULL) {
  a <- commandArgs(TRUE); i <- which(a == name)
  if (length(i) == 1 && length(a) > i) a[i + 1] else default
}
merge <- getarg("--merge", "majority")
rt    <- as.integer(getarg("--rt", "20"))
recompute <- "--recompute" %in% commandArgs(TRUE)

TXTYPE <- "ncrna"
PLATFORM <- "RNA002"
BED <- file.path(GTR, "bed", PLATFORM)

#: tool order and classes resolved for R3-9 / the revised Figure 7
TOOLS <- c("CHEUI_m5C", "NanoMUD_psi", "NanoMUD_m1psi", "NanoNm", "NanoPsu",
           "NanoSPA_psU")
MOD_OF <- c(CHEUI_m5C = "m5C", NanoMUD_psi = "Psi", NanoMUD_m1psi = "m1Psi",
            NanoNm = "Nm", NanoPsu = "Psi", NanoSPA_psU = "Psi")
COND_LABEL <- c(HeLa_WT = "WT", HeLa_IVT = "IVT")
#: panel/legend display names (Greek letter for the Psi tools, hyphen for m5C/m1Psi)
LABEL_OF <- c(CHEUI_m5C = "CHEUI-m5C", NanoMUD_psi = "NanoMUD-\u03a8",
              NanoMUD_m1psi = "NanoMUD-m1\u03a8", NanoNm = "NanoNm",
              NanoPsu = "NanoPsu", NanoSPA_psU = "NanoSPA-\u03a8")

FS <- list(axis = 8, ytitle = 8.5, legend = 7, struct = 7.5)
stopifnot(min(unlist(FS)) >= 7)
NCOL <- 3
#: the metagene grid keeps the physical footprint of the Figure 7 row D block
#: (180 mm x 101.6 mm) so the two figures look like siblings on the SI page
ROW_W <- 7.087; ROW_H <- 4.00

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
          tool = tool, display = unname(LABEL_OF[tool]), mod = mod,
          cond_group = cond, role = role, kind = kind, unit_tag = tag,
          path = path, n_sites = if (is.na(n)) -1L else as.integer(n),
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
  out
}

## Guitar's samplePoints dies on the first malformed group; isolate per group so
## one bad tool x unit BED cannot kill the whole grid.
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
  cache <- file.path(TABD, "figS8_density.rds")
  if (file.exists(cache) && !recompute) {
    message("using cached density: ", basename(cache)); return(readRDS(cache))
  }
  plan <- build_plan()
  if (!nrow(plan)) stop("no input BED found")
  gt <- guitar_txdb("Human", TXTYPE)
  beds <- as.list(plan$path); names(beds) <- plan$group
  message(sprintf("%d curves over %d tools (ncRNA model)", nrow(plan),
                  length(unique(plan$tool))))
  s <- group_sites(beds, TXTYPE, gt)
  if (length(s$errors)) {
    message("failed groups (", length(s$errors), ")")
  }
  dens <- Guitar:::.generateDensity_CI(s$relative, s$weight,
                                       CI_ResamplingTime = rt, adjust = 1,
                                       enableCI = FALSE)
  res <- list(dens = dens, plan = plan,
              componentWidth = gt[[TXTYPE]]$componentWidthAverage_pct,
              merge = merge, rt = rt, errors = s$errors)
  saveRDS(res, cache)
  write.table(plan, file.path(TABD, "figS8_panel_inputs.tsv"), sep = "\t",
              quote = FALSE, row.names = FALSE)
  # Guitar's own per-site coordinates, for the cross-check against 61 (Python)
  site_rows <- do.call(rbind, lapply(names(s$relative), function(g) data.frame(
    group = g, relative = as.numeric(s$relative[[g]]),
    weight = as.numeric(s$weight[[g]]), stringsAsFactors = FALSE)))
  write.table(site_rows, file.path(TABD, "figS8_guitar_sites.tsv"), sep = "\t",
              quote = FALSE, row.names = FALSE)
  message("cached -> ", cache)
  res
}

## ---- transcript schematic + panel frame (GuitarPlot look, see 48) ------------
MM_PER_PT <- as.numeric(grid::convertUnit(grid::unit(1, "pt"), "mm", valueOnly = TRUE))

f7_pos <- function(peak) {
  pos <- Guitar:::.generate_pos_para(peak)
  pos$fig_bottom    <- -0.16 * peak
  pos$rna_comp_text <- -0.085 * peak
  pos
}

structure_layers <- function(comp_width, pos) {
  out <- list()
  lab_map <- c(promoter = "1kb", utr5 = "5'UTR", cds = "CDS", utr3 = "3'UTR",
               tail = "1kb", ncrna = "ncRNA")
  h_map   <- c(promoter = 0.002, utr5 = 0.005, cds = 0.010, utr3 = 0.005,
               tail = 0.002, ncrna = 0.010)
  a_map   <- c(promoter = 0.99, utr5 = 0.99, cds = 0.272, utr3 = 0.99,
               tail = 0.99, ncrna = 0.20)
  nm  <- names(comp_width)
  end <- cumsum(as.numeric(comp_width)); sta <- c(0, end[-length(end)]) + 0.001
  mid <- (sta + end) / 2; pk <- pos$fig_top / 1.05
  h <- unname(h_map[nm]) * pk
  xlab <- pmin(pmax(mid, 0.05), 0.95)
  for (i in seq_along(nm)) {
    out[[length(out) + 1]] <- annotate("rect", xmin = sta[i], xmax = end[i],
                                       ymin = pos$rna_lgd_bl - h[i],
                                       ymax = pos$rna_lgd_bl + h[i],
                                       fill = "grey35", colour = NA,
                                       alpha = unname(a_map[nm[i]]))
    out[[length(out) + 1]] <- annotate("text", x = xlab[i],
                                       y = pos$rna_comp_text,
                                       label = unname(lab_map[nm[i]]),
                                       family = "Arial",
                                       size = FS$struct * MM_PER_PT)
  }
  if (length(nm) > 1)
    out[[length(out) + 1]] <- geom_segment(
      data = data.frame(x = end[-length(end)]), inherit.aes = FALSE,
      aes(x = x, xend = x, y = pos$fig_top, yend = pos$rna_lgd_bl),
      linetype = "dotted", linewidth = 0.25, colour = "black")
  out
}

COLOUR_OF <- c(WT = "#4b81b8", IVT = "#e8a76b")

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
  if (!is.finite(peak) || peak <= 0) return(NULL)
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
    # the two condition keys of the longest label (NanoMUD-m1Psi-IVT) do not fit
    # side by side inside a 60 mm panel, so the condition legend uses two rows
    guides(colour = guide_legend(nrow = 2, order = 1,
                                 override.aes = list(linewidth = 1.2,
                                                     linetype = "solid")),
           
           ## stated once in the caption instead of on every panel
           linetype = "none") +
    structure_layers(res$componentWidth, pos) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(pos$fig_bottom, peak * 1.05),
                       expand = c(0, 0), breaks = c(0, peak / 2, peak),
                       labels = function(v) formatC(v, format = "g", digits = 2)) +
    labs(x = NULL, colour = NULL, linetype = NULL,
         y = if (show_y) "Density" else NULL) +
    theme_classic(base_family = "Arial", base_size = FS$axis) +
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(colour = "black", linewidth = 0.5),
          axis.text = element_text(colour = "black", size = FS$axis),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title = element_text(size = FS$ytitle),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.box = "vertical",
          legend.box.just = "center",
          legend.key.width = grid::unit(12, "pt"),
          legend.key.height = grid::unit(7, "pt"),
          legend.spacing.x = grid::unit(2, "pt"),
          legend.text = element_text(size = FS$legend),
          legend.margin = margin(0, 0, 0, 0),
          legend.box.spacing = grid::unit(1, "pt"),
          plot.margin = margin(2, 3, 1, 2))
}

## ---- assemble the metagene block --------------------------------------------
main <- function() {
  res <- compute_or_load()
  panels <- lapply(seq_along(TOOLS), function(i) {
    draw_panel(res, TOOLS[i], show_y = ((i - 1) %% NCOL) == 0)
  })
  names(panels) <- TOOLS
  keep <- !vapply(panels, is.null, logical(1))
  if (any(!keep)) message("panels drawn: ", sum(keep), " of ", length(TOOLS),
                          " (missing: ", paste(TOOLS[!keep], collapse = ", "), ")")
  panels <- panels[keep]
  saveRDS(panels, file.path(TABD, "figS8_panels.rds"))
  row <- patchwork::wrap_plots(panels, ncol = NCOL)
  width <- ROW_W; height <- ROW_H
  save_pair_cairo(row, file.path(FIGD, "FigureS8_rev_metagene.pdf"), width, height)
  geom <- data.frame(
    item = c("panels", "ncol", "width_in", "height_in", "font_min_pt",
             "component_width_promoter", "component_width_ncrna",
             "component_width_tail", "merge", "rt", "failed_groups"),
    value = c(length(panels), NCOL, width, height, min(unlist(FS)),
              as.numeric(res$componentWidth["promoter"]),
              as.numeric(res$componentWidth["ncrna"]),
              as.numeric(res$componentWidth["tail"]), res$merge, res$rt,
              length(res$errors)))
  write.table(geom, file.path(TABD, "figS8_guitar_geometry.tsv"), sep = "\t",
              quote = FALSE, row.names = FALSE)
  message("done")
}

main()
