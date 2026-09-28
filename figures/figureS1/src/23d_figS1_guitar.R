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
# Supplementary Figure S1, redrawn replicate-aware.
#
# Why this exists.  The published sup1.pdf ($RNAMODBENCH_LOCAL/submission/02_AS_working_copy_and_revisions/sup/)
# is a 3 species x 10 tools metagene grid, but its inputs carried no replicate
# structure: per $RNAMODBENCH_LOCAL/REPLICATE_AWARE_FIGURES.md the Arabidopsis
# panels were rep3 only, the Mouse panels the mES_WT study only, and the HeLa
# panels an undocumented union of the three replicates (strand forced to "+").
# This script redraws the same grid from
#
#   bed/RNA002/<merge>/<Condition>/m6A/<Tool>.bed        majority consensus
#   bed/RNA002/rep_<tag>/<Condition>/m6A/<Tool>.bed      one independent unit
#
# written by 21b_export_guitar_bed.py, and draws, inside every tool panel, the
# consensus curve of each library (thick, filled) plus one thin curve per
# replicate / study, so the replicate structure R3-2 / E6 asks for is visible in
# the figure itself.
#
# The density is computed by the Bioconductor Guitar package itself
# (samplePoints -> normalize -> .generateDensity_CI); only the panel layers and
# the page layout are ours, because Guitar cannot draw per-tool facets or
# distinguish consensus from replicate curves.
#
# Usage
#   conda run -n guitar_asm --no-capture-output Rscript 23d_figS1_guitar.R
#   ... --species Mouse            # one species only (density is cached)
#   ... --recompute                # ignore the cached RDS and recompute
#   ... --min-sites 10 --rt 20
#
# Outputs (kept in their own tree, NOT mixed with guitar_metagene_replicates/)
#   figures/figureS1/
#     figures/FigureS1_rev.{pdf,png}            the stitched page (~17.6x17.2 in,
#                                               natural-size panels, NOT A4)
#     figures/FigureS1_rev_<Species>.{pdf,png}  per-species block (QC)
#     tables/figS1_density_<Species>.rds        cached density + inputs
#     tables/figS1_panel_inputs.tsv             one row per drawn curve
#     tables/figS1_geometry.tsv                 layout / font self-check
#     logs/                                     run logs
#
# Fonts: the PDF is written with grDevices::cairo_pdf while showtext is OFF, so
# the real system Arial is EMBEDDED (subsetted TrueType, pdffonts-verifiable);
# showtext is switched back on only for the 300 dpi PNG preview, where glyph
# outlines are fine.  Nothing below 7 pt.

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
arial_setup()   # registers Arial for showtext (PNG path); toggled off for PDFs

## ---- export: PDF with the real embedded Arial, PNG as a raster preview ------
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

GTR    <- file.path(.XB, "harmonisation/guitar_metagene")
BED    <- file.path(GTR, "bed", "RNA002")   # read-only inputs from the other tree
OUTD   <- file.path(.RB, "figures/figureS1")
TABD   <- file.path(OUTD, "tables")
FIGD   <- file.path(OUTD, "figures")
LOGD   <- file.path(OUTD, "logs")
REG    <- file.path(.RB, "metadata/sample_registry.csv")

# --- house style: English, Arial, no gridlines --------------------------------
# ONE font profile, one panel canvas: the standalone files are the page panels,
# stitched 1:1 at their natural size.
# 2026-09-23: the page used to be 17.55 x 17.16 in, which the SI layout
# (\includegraphics[width=\textwidth,height=0.80\textheight,keepaspectratio])
# squeezed by 0.399 -- every 8 pt label then printed at 3.2 pt, far below the
# 7 pt floor.  The grid is therefore drawn at the size it is printed: five
# columns of 1.33 in inside the 7.0 in text width, three blocks down 8.0 in.
FS <- list(axis = 8, ytitle = 8, tool = 11, species = 13, legend = 8.5, struct = 8.5)
stopifnot(min(unlist(FS)) >= 8)
PANEL_W <- 1.33; PANEL_H <- 1.02
NCOL    <- 5
BLOCK_TITLE_H  <- 0.25   # species title row, inches
BLOCK_LEGEND_H <- 0.32   # collected legend strip of each species block
# the page is exactly the sum of its parts (title + 2 panel rows + legend) x 3
PG_W <- round(0.15 + NCOL * PANEL_W + (NCOL - 1) * 0.04, 2)          # 6.96 in
PG_H <- round(0.20 + 3 * (BLOCK_TITLE_H + 2 * PANEL_H + BLOCK_LEGEND_H), 2)  # 8.03 in
stopifnot(PG_W <= 7.0, PG_H <= 8.35)   # the SI text area, so nothing is scaled

## ---- CLI ----
getarg <- function(name, default = NULL) {
  a <- commandArgs(TRUE); i <- which(a == name)
  if (length(i) == 1 && length(a) > i) a[i + 1] else default
}
species_arg <- strsplit(getarg("--species", "Arabidopsis,Mouse,Human"), ",")[[1]]
merge       <- getarg("--merge", "majority")
min_sites   <- as.integer(getarg("--min-sites", "10"))
rt          <- as.integer(getarg("--rt", "20"))
recompute   <- "--recompute" %in% commandArgs(TRUE)   # flag, may be the last arg
mode        <- getarg("--mode", "all")                # panels | assemble | all
only_panel  <- getarg("--only-panel", NULL)           # e.g. Mouse_Nanocompore
gate        <- as.integer(getarg("--gate", "100"))    # consensus drawn only if
                                                        # n_sites >= gate
stopifnot(mode %in% c("panels", "assemble", "all"))

## ---- the ten tools of the published grid, in the order sup1.pdf uses them ----
TOOLS <- c("CHEUI_m6A", "DRUMMER", "ELIGOS2_diff", "ELIGOS2_solo", "EpiNano_Error",
           "MINES", "Nanocompore", "NanoSPA_m6A", "xPore", "yanocomp")
# 2026-09-21: labels follow the manuscript nomenclature (Table S11); the earlier
# "ELIGOS_diff / ELIGOS_solo" spellings are dropped, yanocomp stays capitalised
LABEL_OF <- c(yanocomp = "Yanocomp")
tool_label <- function(t) ifelse(t %in% names(LABEL_OF), unname(LABEL_OF[t]), t)

# Differential tools (ELIGOS2_diff, DRUMMER) attribute significant sites to the
# modified (WT) side by construction; their unmodified (pert: IVT / KO / KD) side
# is nearly empty and is NOT drawn in Figure S1 -- only the WT curve is shown
# (see FigS1_legends.md).  This keeps the 10-tool grid intact while removing the
# misleading near-empty orange curves.
DIFF_TOOLS <- c("ELIGOS2_diff", "DRUMMER")

SPECIES <- list(
  Arabidopsis = list(order = 1, block = "Arabidopsis",
                     conds = c(Arabidopsis_WT = "WT", Arabidopsis_KD = "fip37 KD")),
  Mouse       = list(order = 2, block = "Mouse",
                     conds = c(Mouse_WT = "WT", Mouse_KO = "Mettl3 KO")),
  Human       = list(order = 3, block = "Human",
                     conds = c(HeLa_WT = "WT", HeLa_IVT = "IVT")))

n_lines <- function(p) if (file.exists(p)) length(readLines(p, warn = FALSE)) else NA_integer_

## ---- the independent units of a group (rep_rep1..3 / rep_studyA,B) ---------
unit_tags <- function(cond) {
  d <- list.files(BED, pattern = "^rep_", full.names = FALSE)
  d <- d[vapply(d, function(x) dir.exists(file.path(BED, x, cond)), logical(1))]
  sort(sub("^rep_", "", d))
}
unit_samples <- function(cond, tag) {
  if (!file.exists(REG)) return("")
  reg <- read.delim(REG, stringsAsFactors = FALSE)
  s <- unique(reg$sample[reg$dataset_group == cond & reg$replicate_tag == tag])
  paste(s, collapse = ",")
}

unit_label <- function(sp) if (sp == "Mouse") "individual study" else "individual replicate"

## ---- one row per curve that the panels will draw ----------------------------
build_plan <- function(sp) {
  cfg    <- SPECIES[[sp]]
  conds  <- cfg$conds                       # names = bed group, values = label
  tags   <- unit_tags(names(conds)[1])
  rows   <- list()
  for (tool in TOOLS) {
    for (i in seq_along(conds)) {
      cond  <- names(conds)[i]
      role  <- if (i == 1) "WT" else "pert"
      clab  <- conds[[i]]
      addrow <- function(kind, tag, path, note) {
        n <- n_lines(path)
        rows[[length(rows) + 1]] <<- data.frame(
          species = sp, tool = tool, toollab = tool_label(tool), cond_group = cond,
          role = role, cond_label = clab, kind = kind, unit_tag = tag,
          unit_sample = if (kind == "unit") unit_samples(cond, tag) else "",
          path = path, n_sites = if (is.na(n)) -1L else as.integer(n),
          note = note, stringsAsFactors = FALSE)
      }
      # majority consensus curve
      f <- file.path(BED, merge, cond, "m6A", paste0(tool, ".bed"))
      n <- n_lines(f)
      if (!is.na(n) && n > 0)
        addrow("consensus", "", f,
               if (n < min_sites) paste0("below min-sites(", min_sites, ")") else "")
      # one thin curve per independent unit
      for (tag in tags) {
        f <- file.path(BED, paste0("rep_", tag), cond, "m6A", paste0(tool, ".bed"))
        n <- n_lines(f)
        if (!is.na(n) && n > 0)
          addrow("unit", tag, f,
                 if (n < min_sites) paste0("below min-sites(", min_sites, ")") else "")
      }
    }
  }
  out <- do.call(rbind, rows)
  out$group <- paste(out$tool, out$role, out$kind, out$unit_tag, sep = "|")
  out
}

## ---- guitar_sites with per-group error isolation -----------------------------
# Guitar's samplePoints dies on the first malformed group (observed on Human:
# "Componet_pct[x, ] : subscript out of bounds" for a few tool x unit BEDs).
# Isolating per group keeps the other 79 curves and makes the failures
# auditable instead of fatal.
s1_sites <- function(beds, txtype, gt) {
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
      relative[[gname]] <- nz[[1]]
      weight[[gname]]   <- nz[[2]]
      TRUE
    }, error = function(e) {
      errs <<- c(errs, sprintf("%s: %s", gname, conditionMessage(e)))
      FALSE
    })
    if (!ok) message("!! samplePoints failed: ", gname, " -- curve skipped")
  }
  if (length(errs))
    message("failed groups (", length(errs), "):\n  ",
            paste(errs, collapse = "\n  "))
  list(relative = relative, weight = weight, errors = errs)
}

## ---- density, cached per species -------------------------------------------
compute_or_load <- function(sp) {
  cache <- file.path(TABD, sprintf("figS1_density_%s.rds", sp))
  if (file.exists(cache) && !recompute) {
    message("using cached density: ", basename(cache))
    return(readRDS(cache))
  }
  plan <- build_plan(sp)
  if (!nrow(plan)) stop("no input BED for ", sp)
  gt   <- guitar_txdb(sp, "mrna")
  beds <- as.list(plan$path); names(beds) <- plan$group
  message(sprintf("%s: %d curves over %d tools (%s units)",
                  sp, nrow(plan), length(unique(plan$tool)),
                  paste(unique(plan$unit_tag[plan$kind == "unit"]), collapse = ",")))
  t0  <- Sys.time()
  s   <- s1_sites(beds, "mrna", gt)
  message(sprintf("%s: site sampling %.1fs", sp,
                  as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  if (length(s$errors)) {
    failed <- sub(":.*$", "", s$errors)
    # keep the row (auditability), just flag it; draw_panel finds no density
    # for it and silently omits the curve
    plan$note[plan$group %in% failed] <- paste(
      plan$note[plan$group %in% failed], "samplePoints-failed", sep = "; ")
    beds <- beds[!names(beds) %in% failed]
  }
  # CI_ResamplingTime is irrelevant while enableCI = FALSE (density identical to
  # the last digit between 20 and 1000): see figures/figureS1/logs/figS1_probe.log
  dens <- Guitar:::.generateDensity_CI(s$relative, s$weight,
                                       CI_ResamplingTime = rt, adjust = 1,
                                       enableCI = FALSE)
  res <- list(dens = dens, plan = plan,
              componentWidth = gt$mrna$componentWidthAverage_pct,
              species = sp, merge = merge, min_sites = min_sites, rt = rt)
  saveRDS(res, cache)
  message("cached -> ", cache)
  res
}

## ---- the transcript schematic under every panel -----------------------------
# Mirrors Guitar:::.RNAPlotStructure (same rectangles, alphas, dotted component
# boundaries and labels) with two deliberate differences: the label size is
# pinned in pt so it stays >= 7 pt in these small panels, and the labels run on
# TWO levels (5'UTR / 3'UTR one line lower, with a short connector) so that
# 1kb / 5'UTR / CDS / 3'UTR / 1kb never collide inside a ~34 mm panel.
MM_PER_PT <- as.numeric(grid::convertUnit(grid::unit(1, "pt"), "mm",
                                          valueOnly = TRUE))

## panel cell, so they are packed with a minimum gap and the ones that moved get
## a leader line back to their segment (region_labels() in guitar_lib.R).  The
## row sits at -0.40 * peak: that is below the band of the "0" y tick, which is
## what forced the two-row layout before (the tick numbers and the labels shared
## a band and could not both keep their space).  fig_bottom is unchanged, so the
## curves keep the exact geometry they were audited with.
## Label-row geometry measured on the delivered page (pdftotext -bbox):
## plot area 50.3 pt wide (delta x of the two clamped "1kb" anchors / 0.90),
## panel cell 95.76 pt, column pitch 98.64 pt.
S1_AXES_W_PT   <- 50.3   # reference value; each panel measures its own width
#: per-panel geometry, written to tables/figS1_panel_geometry.tsv at the end so
#: the audit can check the label row against the real plot area
#: narrowest plot area measured on this layout (2026-09-23); the tier is derived
#: from it so every panel keeps the row inside its own plot area
S1_PLOT_MIN_PT <- 72.0
GEO <- new.env(parent = emptyenv())
GEO$rows <- list()
S1_LABEL_ALLOW <- c(-0.50, 1.40)   # data units: cell edge to the next cell's edge

s1_pos <- function(peak) {
  pos <- Guitar:::.generate_pos_para(peak)
  
  # title (1.7 pt of air at -0.40).  The row moves down to -0.48 and the range
  # grows with it so the row keeps ~5 pt above and ~2 pt below.
  
  ## sits 3-4 pt below it, and the panel floor comes up with it (page height is
  ## fixed by PG_H, so this only changes the data-to-pt mapping).
  pos$fig_bottom     <- -0.30 * peak
  pos$rna_comp_text  <- -0.20 * peak
  pos
}
add_structure <- function(p, comp_width, pos, labels = TRUE,
                          label_cfg = NULL) {
  lab_map   <- c(promoter = "1kb", utr5 = "5'UTR", cds = "CDS", utr3 = "3'UTR",
                 tail = "1kb", ncrna = "ncRNA", rna = "RNA")
  height_map <- c(promoter = 0.002, utr5 = 0.005, cds = 0.010, utr3 = 0.005,
                  tail = 0.002, ncrna = 0.010, rna = 0.010)
  alpha_map  <- c(promoter = 0.99, utr5 = 0.99, cds = 0.272, utr3 = 0.99,
                  tail = 0.99, ncrna = 0.20, rna = 0.20)
  nm  <- names(comp_width)
  end <- cumsum(as.numeric(comp_width))
  sta <- c(0, end[-length(end)]) + 0.001
  mid <- (sta + end) / 2
  pk  <- pos$fig_top / 1.05                          # the peak the pos was built on
  h    <- unname(height_map[nm]) * pk                # = fraction * peak
  # never let a region label sit on the panel edge (the flanks especially):
  # hjust = 0.5 text anchored inside [0.05, 0.95] keeps "1kb" readable
  xlab <- pmin(pmax(mid, 0.05), 0.95)
  for (i in seq_along(nm)) {
    p <- p + annotate("rect", xmin = sta[i], xmax = end[i],
                      ymin = pos$rna_lgd_bl - h[i], ymax = pos$rna_lgd_bl + h[i],
                      fill = "grey35", colour = NA, alpha = unname(alpha_map[nm[i]]))
  }
  if (labels && !is.null(label_cfg)) {
    # one tier for the whole figure, derived from the first measured plot area:
    # the row has to fit inside the plot area (not just the panel cell) with a
    # 0.6 pt gap between labels and 1 pt of margin on each side
    if (is.null(GEO$tier)) {
      win <- min(label_cfg$axes_w_pt, S1_PLOT_MIN_PT)
      # margin 3.5 pt: the rendered row measures ~8% wider than the packed
      # estimate (measured 2026-09-23), so the margin absorbs that slack
      GEO$tier <- region_size_for(win, unname(lab_map[nm]),
                                  gap_pt = 3.0, margin_pt = 5.0)
      message(sprintf("tier from the %.1f pt narrowest plot area: %.1f pt, gap %.1f pt, row %.1f pt",
                      win, GEO$tier$size, GEO$tier$gap, GEO$tier$row_pt))
    }
    fit <- fit_region_size(unname(lab_map[nm]), label_cfg$allow_w_pt,
                           ladder = data.frame(size = GEO$tier$size,
                                               gap = GEO$tier$gap))
    if (is.null(fit)) {
      stop(sprintf(paste0("region labels need more than %.1f pt at the 7 pt ",
                          "floor; widen S1_LABEL_ALLOW or shorten the labels"),
                   label_cfg$allow_w_pt))
    }
    rl <- region_labels(xlab, unname(lab_map[nm]), y = pos$rna_comp_text,
                        axes_w_pt = label_cfg$axes_w_pt,
                        allow_data = label_cfg$allow_data,
                        size = fit$size, gap = fit$gap,
                        anchor_y = pos$rna_lgd_bl - h,
                        leader_drop = label_cfg$leader_drop,
                        mode = "segments", rotate = character(0),
                        nudge_pt = 4.0)   # <=4 pt: at 6 pt the row fills 68 of the 72 pt
                        # plot area, so neighbours need a little more room;
                        # the shift is invisible at this scale
    p <- p + rl$text + rl$leaders
    attr(p, "region_fit") <- list(size = fit$size, gap = fit$gap,
                                  need_pt = fit$need_pt, fits = rl$fits,
                                  leaders = rl$n_leaders,
                                  max_moved_pt = rl$max_moved_pt)
    message(sprintf(paste0("%s region labels: %.1f pt, gap %.1f pt, %.1f of ",
                           "%.1f pt used, fits %s, leaders %d, max moved %.2f pt"),
                    label_cfg$tag, fit$size, fit$gap, fit$need_pt,
                    label_cfg$allow_w_pt, rl$fits, rl$n_leaders, rl$max_moved_pt))
  }
  if (length(nm) > 1) {
    vp <- data.frame(x = end[-length(end)])
    p <- p + geom_segment(data = vp, inherit.aes = FALSE,
                          
                          ## area only -- they used to run down to the schematic
                          ## band, i.e. through the label row.
                          aes(x = x, xend = x, y = pos$fig_top,
                              yend = 0),
                          linetype = "dotted", linewidth = 0.25, colour = "black")
  }
  p
}

## ---- one panel per tool ------------------------------------------------------
# draw_panel() is THE panel: the standalone files and the stitched page use the
# same function, same slot size, same fonts.  `legend = FALSE` (panels mode)
# drops the in-panel legend — on the page the legend appears once per species
# block (guides="collect"); in standalone files the color code lives in the
# caption.  Truncation fixes baked in: no numeric x ticks (region labels are
# the x axis), generous top/left margins for the 8 pt tool title and the
# rotated "Density".
draw_panel <- function(res, tool, show_x, legend = TRUE) {
  sp   <- res$species
  plan <- res$plan[res$plan$tool == tool, ]
  d    <- res$dens[res$dens$group %in% plan$group, ]
  if (!nrow(d)) return(NULL)
  meta <- plan[, c("group", "role", "kind", "unit_tag", "cond_label", "cond_group",
                   "n_sites")]
  d    <- merge(d, meta, by = "group", all.x = TRUE)
  d$condition <- d$cond_label
  d$curvetype <- ifelse(d$kind == "consensus", "majority consensus", unit_label(sp))
  # differential tools: the unmodified (pert) side is suppressed entirely only
  # when it is genuinely near-empty (its consensus is below the gate).  The
  # modified (WT) side, and any substantial pert side (e.g. Mouse KO, which
  # carries a large consensus), are always drawn -- this removes the
  # broken-looking empty IVT/KD curves without discarding real signal
  # (FigS1_legends.md).
  if (tool %in% DIFF_TOOLS) {
    pert_cons_n <- plan$n_sites[plan$role == "pert" & plan$kind == "consensus"]
    if (length(pert_cons_n) >= 1 && isTRUE(pert_cons_n[1] < gate)) {
      plan <- plan[plan$role == "WT", , drop = FALSE]
      d    <- d[d$group %in% plan$group, , drop = FALSE]
    }
  }
  # consensus gating (user 2026-09-20): a majority consensus built from fewer
  # than `gate` sites is NOT drawn as a thick filled curve — the per-unit thin
  # lines still are.  Keeps e.g. ELIGOS2_diff HeLa_IVT (24 sites) honest.
  gate_keep <- plan$cond_group[plan$kind == "consensus" & plan$n_sites >= gate]
  cons  <- d[d$kind == "consensus" & d$cond_group %in% gate_keep, ]
  units <- d[d$kind == "unit", ]

  # colour only the conditions actually present (differential tools lose the pert
  # colour because their unmodified side is suppressed above).
  cols_all <- setNames(c(WT_COLOR, CMP_COLOR), c("WT", SPECIES[[sp]]$conds[[2]]))
  cols <- cols_all[names(cols_all) %in% unique(plan$cond_label)]
  lin  <- setNames(c("solid", "dashed"), c("majority consensus", unit_label(sp)))

  peak <- max(d$density, na.rm = TRUE)
  pos  <- s1_pos(peak)

  p <- ggplot(d, aes(x = x))
  if (nrow(units))
    p <- p + geom_line(data = units, aes(y = density, group = group,
                                         colour = condition, linetype = curvetype),
                       linewidth = 0.30, alpha = 0.75)
  if (nrow(cons)) {
    p <- p + geom_ribbon(data = cons, aes(ymin = 0, ymax = density, group = group,
                                          fill = condition),
                         alpha = 0.20, colour = NA)
    p <- p + geom_line(data = cons, aes(y = density, group = group,
                                        colour = condition, linetype = curvetype),
                       linewidth = 0.75)
  }
  p <- p +
    
    ## colour scale (full palette + both condition limits).  The panels
    ## differed before (each kept only the conditions it drew), so
    ## plot_layout(guides = "collect") could not merge the colour keys and
    ## the strip carried 2-3 copies of them; the linetype scale, which was
    ## already identical everywhere, merged fine -- that is the proof.
    scale_colour_manual(values = cols_all, limits = names(cols_all)) +
    scale_fill_manual(values = cols_all, limits = names(cols_all)) +

    scale_linetype_manual(values = lin, limits = names(lin)) +
    # the packed label row deliberately runs past the panel range (five labels
    # are wider than the plot area).  The panel range must therefore come from
    # coord_cartesian(): limits = c(0, 1) on the *scale* censors every value
    # outside it to NA, which silently deleted the two outer "1kb" labels from
    # the PDF (found 2026-09-23: no "kb" token was extractable at all).
    scale_x_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = c(0, 1), expand = FALSE, clip = "off") +
    # no ticks below zero: that band is the single label row's own space now
    scale_y_continuous(limits = c(pos$fig_bottom, pos$fig_top), expand = c(0, 0),
                       breaks = function(lim) {
                         b <- scales::breaks_pretty(4)(lim); b[b >= 0]
                       }) +
    labs(x = NULL, y = NULL, title = tool_label(tool)) +
    
    ## "none", so plot_layout(guides = "collect") has exactly one colour key and
    ## one linetype key to collect.  (Two earlier attempts failed: identical scale
    ## limits were still not merged, and switching legend.position per panel was
    ## overridden by the block-level `&` theme and changed the panel areas.)
    guides(colour = if (legend) guide_legend(override.aes = list(linewidth = 0.8,
                                                                 linetype = "solid",
                                                                 alpha = 1)) else "none",
           fill = "none",
           linetype = if (legend) guide_legend(override.aes = list(colour = "grey35",
                                                                   linewidth = c(0.8, 0.30))) else "none") +
    theme_classic(base_family = "Arial", base_size = FS$axis) +
    theme(panel.grid = element_blank(),
          panel.background = element_rect(fill = "white", colour = NA),
          axis.line = element_line(colour = "black", linewidth = 0.3),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),       # region labels are the x axis;
                                               # numeric 0..1 ticks clipped right
          # the density axis is stated once per block now: dropping the rotated
          # per-panel title frees ~10 pt of width for the label row (2026-09-23)
          axis.title.y = element_blank(),
          axis.text.y = element_text(size = FS$axis, colour = "black"),
          legend.position = if (legend) "bottom" else "none",
          legend.title = element_blank(),
          legend.text = element_text(size = FS$legend),
          legend.key.width = unit(10, "pt"),
          legend.key.height = unit(4, "pt"),
          legend.margin = margin(0, 0, 0, 0, "pt"),
          legend.box.margin = margin(0, 0, 0, 0, "pt"),
          plot.title = element_text(size = FS$tool, face = "bold",
                                    hjust = 0.5, colour = "black"),
          plot.margin = margin(5, 6, 4, 6, "pt"),
          text = element_text(family = "Arial"))
  # the packed label row has to fit the *plot area*, and that width differs per
  # panel (the y tick text is not the same width everywhere).  Read it from the
  # built gtable instead of trusting S1_AXES_W_PT, which is only the reference
  # value measured on the delivered page.
  g <- ggplot2::ggplotGrob(p)
  # the panel column carries a null unit, which convertWidth() reports as 0;
  # everything else (axis tick text, the axis title, the plot margins) resolves
  # to real points, so the plot area is what is left of the sub-plot's width
  fixed_pt <- sum(grid::convertWidth(g$widths, "pt", valueOnly = TRUE))
  axes_w_pt <- PANEL_W * 72 - fixed_pt
  if (!is.finite(axes_w_pt) || axes_w_pt < 20) {
    warning(sprintf("could not measure the plot area of %s (%.1f pt); falling back to %.1f pt",
                    tool_label(tool), axes_w_pt, S1_AXES_W_PT), call. = FALSE)
    axes_w_pt <- S1_AXES_W_PT
  }
  if (abs(axes_w_pt - S1_AXES_W_PT) > 0.15 * S1_AXES_W_PT) {
    message(sprintf("panel %s: plot area %.1f pt (reference %.1f pt)",
                    tool_label(tool), axes_w_pt, S1_AXES_W_PT))
  }
  GEO$rows[[length(GEO$rows) + 1]] <-
    data.frame(species = sp, tool = tool, axes_w_pt = round(axes_w_pt, 2),
               stringsAsFactors = FALSE)
  add_structure(p, res$componentWidth, pos, labels = show_x,
                label_cfg = if (show_x) list(tag = tool_label(tool),
                                             axes_w_pt = axes_w_pt,
                                             allow_data = S1_LABEL_ALLOW,
                                             allow_w_pt = diff(S1_LABEL_ALLOW) *
                                               S1_AXES_W_PT,
                                             leader_drop = -0.13 * peak) else NULL)
}

## ---- one block per species: 2 rows x 5 columns of tool panels ---------------
draw_block <- function(sp, res) {
  tools  <- TOOLS
  nrow_b <- ceiling(length(tools) / NCOL)
  panels <- list()
  for (i in seq_along(tools)) {
    # every panel carries its own region labels: "some panels have no x axis"
    # was the reading when only the bottom row of a block was labelled
    show_x <- TRUE
    
    ## TRUE on every panel made plot_layout(guides = "collect") stack a colour
    ## key per panel (their scales differed by limits), so the strip carried
    ## 2-3 copies of WT/condition.
    p <- draw_panel(res, tools[i], show_x, legend = (length(panels) == 0))
    if (is.null(p)) { message("!! no curve for ", sp, "/", tools[i]); next }
    panels[[length(panels) + 1]] <- p
  }
  if (!length(panels)) return(NULL)
  message(sprintf("%s: %d panels", sp, length(panels)))
  # Species title as a REAL plot row, stacked ABOVE the panel grid with the
  # patchwork "/" operator.  Two earlier attempts failed and are documented
  # here so they are not retried:
  #   * plot_annotation(title=) -- silently dropped when the block is nested
  #     into the page patchwork (titles vanished, 2026-09-20);
  #   * wrap_plots(11 plots) + plot_layout(design="ttttt\nabcde\nfghij") --
  #     the pre-filled 3-row grid (row 3 = the 11th plot, Yanocomp) fights the
  #     design remapping, scattering panels and pushing Yanocomp to the page
  #     top (2026-09-20).  The "/" stack is the verified-safe structure: the
  #     outer rows are one full-width title and one inner 2x5 grid whose
  #     TOOLS order and collected legend are untouched by nesting.
  if (length(panels) == length(TOOLS)) {
    title_p <- ggplot() +
      annotate("text", x = 0.5, y = 0.5,
               
               #: ("Arabidopsis — metagene density"); a colon is the plain-human
               #: punctuation here and keeps the same meaning
               label = paste0(SPECIES[[sp]]$block, ": metagene density"),
               family = "Arial", fontface = "bold",
               size = FS$species * MM_PER_PT) +
      theme_void() + theme(text = element_text(family = "Arial"))
    
    ## panel and thus defeated the per-panel switch; the single panel that keeps
    ## its key already places it at the top (draw_panel), so the strip stays with
    ## the block title and no entry is duplicated.
    inner <- wrap_plots(panels, ncol = NCOL) + plot_layout(guides = "collect") &
      theme(legend.position = "top")   # key with the block title, not under it
    # wrap_elements(full=): patchwork otherwise aligns the stacked rows by
    # their PANEL areas, which pulled the title ~1.2 in off the page centre
    # (observed x_center 546 pt vs page 632 pt); full= spans the real width.
    blk <- wrap_elements(full = title_p) / inner +
      plot_layout(heights = c(BLOCK_TITLE_H, 2 * PANEL_H + BLOCK_LEGEND_H))
  } else {
    message("!! ", sp, ": ", length(panels), "/", length(TOOLS),
            " panels -- drawing without the title row")
    blk <- wrap_plots(panels, ncol = NCOL) + plot_layout(guides = "collect")
  }
  blk
}

## ---- standalone panel files: one figure per species x tool -------------------
# Same function/size/fonts as the page panels (slot size), legend omitted —
# the colour code is documented in FigS1_legends.md.  `--only-panel
# Mouse_Nanocompore` iterates a single panel.
panel_targets <- function(sp) {
  if (is.null(only_panel)) return(TOOLS)
  parts <- strsplit(only_panel, "_", fixed = TRUE)[[1]]
  if (parts[1] != sp) return(character(0))
  paste(parts[-1], collapse = "_")           # tool names contain underscores
}
save_standalone_panels <- function(sp, res) {
  tools <- panel_targets(sp)
  if (!length(tools)) return(invisible())
  dir.create(file.path(FIGD, "panels"), showWarnings = FALSE, recursive = TRUE)
  for (tool in tools) {
    p <- draw_panel(res, tool, show_x = TRUE, legend = FALSE)
    if (is.null(p)) { message("!! no curve for ", sp, "/", tool); next }
    save_pair_cairo(p,
                    file.path(FIGD, "panels",
                              sprintf("FigureS1_rev_%s_%s.pdf", sp, tool)),
                    width = PANEL_W, height = PANEL_H)
  }
}

## ---- run --------------------------------------------------------------------
species_arg <- intersect(species_arg, names(SPECIES))
species_arg <- species_arg[order(sapply(species_arg, function(s) SPECIES[[s]]$order))]
if (!is.null(only_panel)) {                   # single-panel iteration: panels only
  sp0 <- strsplit(only_panel, "_", fixed = TRUE)[[1]][1]
  species_arg <- intersect(sp0, names(SPECIES))
  mode <- "panels"
}
dir.create(TABD, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGD, showWarnings = FALSE, recursive = TRUE)

res_list <- list(); plans <- list()
for (sp in species_arg) {
  message("\n=== ", sp, " ", format(Sys.time()), " ===")
  res <- compute_or_load(sp)
  res_list[[sp]] <- res
  plans[[sp]]    <- res$plan
}

if (mode %in% c("panels", "all"))
  for (sp in species_arg) save_standalone_panels(sp, res_list[[sp]])

blocks <- list()
if (mode %in% c("assemble", "all")) {
  for (sp in species_arg) {
    blk <- draw_block(sp, res_list[[sp]])
    if (is.null(blk)) next
    blocks[[sp]] <- blk
    save_pair_cairo(blk, file.path(FIGD, sprintf("FigureS1_rev_%s.pdf", sp)),
                    width = PG_W,
                    height = BLOCK_TITLE_H + 2 * PANEL_H + BLOCK_LEGEND_H + 0.06)
  }
  if (length(blocks) == length(species_arg)) {
    page <- wrap_plots(blocks[species_arg], ncol = 1)
    save_pair_cairo(page, file.path(FIGD, "FigureS1_rev.pdf"),
                    width = PG_W, height = PG_H)
    message("wrote FigureS1_rev.{pdf,png}")
  } else {
    message("page not assembled (blocks: ",
            paste(names(blocks), collapse = ","), ")")
  }
}

## ---- self-check tables ------------------------------------------------------
if (length(plans)) {
  inp <- do.call(rbind, plans)
  inp$panel_index <- match(inp$tool, TOOLS)
  inp$panel_row   <- ceiling(inp$panel_index / NCOL)
  inp$panel_col   <- (inp$panel_index - 1) %% NCOL + 1
  write.table(inp, file.path(TABD, "figS1_panel_inputs.tsv"), sep = "\t",
              quote = FALSE, row.names = FALSE)
  # consensus curves suppressed for differential tools (their unmodified side is
  # never drawn) are recorded separately from the merely gated ones, so the table
  # stays honest about what is actually shown in the figure.
  suppressed_diff <- inp$kind == "consensus" & inp$tool %in% DIFF_TOOLS &
                    inp$role == "pert" & inp$n_sites < gate
  gated_df <- inp[inp$kind == "consensus" & inp$n_sites < gate & !suppressed_diff,
                  c("species", "tool", "cond_group", "cond_label", "n_sites")]
  gated_df$status <- "gated"
  supp_df <- inp[suppressed_diff,
                 c("species", "tool", "cond_group", "cond_label", "n_sites")]
  supp_df$status <- "suppressed_diff"
  gated <- rbind(gated_df, supp_df)
  write.table(gated, file.path(TABD, "figS1_gated_consensus.tsv"), sep = "\t",
              quote = FALSE, row.names = FALSE)
  geo <- data.frame(
    item = c("page_in", "blocks", "panels_per_block", "ncol", "tools_drawn",
             "min_font_pt", "font_axis_ytitle_tool_species_legend_struct_pt",
             "panel_size_in", "panel_legend", "consensus_gate",
             "gated_consensus_n", "font_backend", "label_layout",
             "merge", "min_sites", "CI_ResamplingTime", "consensus_rule"),
    value = c(sprintf("%.2f x %.2f (natural-size stitch, not A4)", PG_W, PG_H),
              as.character(length(blocks)),
              as.character(length(TOOLS)), as.character(NCOL),
              as.character(length(unique(inp$tool))),
              as.character(min(unlist(FS))),
              paste(FS$axis, FS$ytitle, FS$tool, FS$species, FS$legend,
                    FS$struct, sep = "/"),
              sprintf("%.1f x %.1f", PANEL_W, PANEL_H),
              "none in panels (block-level collected legend)",
              as.character(gate),
              as.character(nrow(gated_df)),
              "cairo_pdf with embedded Arial (showtext off); PNG via showtext",
              "region labels on one line under the panel",
              merge, as.character(min_sites),
              as.character(rt),
              "sites supported by strictly more than half of the group's independent units"))
  write.table(geo, file.path(TABD, "figS1_geometry.tsv"), sep = "\t",
              quote = FALSE, row.names = FALSE)
  low <- inp[inp$note != "", c("species", "tool", "cond_group", "kind",
                               "unit_tag", "n_sites", "note")]
  if (nrow(low)) message("low-input curves:\n", paste(apply(low, 1, paste, collapse = " "),
                                                      collapse = "\n"))
  if (nrow(gated_df)) message("gated consensus curves (n < ", gate, ", thin lines only):\n",
                           paste(apply(gated_df, 1, paste, collapse = " "),
                                 collapse = "\n"))
  if (nrow(supp_df)) message("suppressed consensus curves (differential-tool unmodified side, not drawn):\n",
                           paste(apply(supp_df, 1, paste, collapse = " "),
                                 collapse = "\n"))
  if (length(GEO$rows)) {
    write.table(do.call(rbind, GEO$rows),
                file.path(TABD, "figS1_panel_geometry.tsv"), sep = "\t",
                row.names = FALSE, quote = FALSE)
    w <- vapply(GEO$rows, function(d) d$axes_w_pt, numeric(1))
    message(sprintf("plot area per panel: min %.1f, median %.1f, max %.1f pt (%d panels)",
                    min(w), stats::median(w), max(w), length(w)))
  }
  message("wrote figS1_panel_inputs.tsv / figS1_gated_consensus.tsv / figS1_geometry.tsv / figS1_panel_geometry.tsv")
}
message("done ", format(Sys.time()))
