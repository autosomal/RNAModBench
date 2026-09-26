#!/usr/bin/env Rscript
# Redraw every GUITAR panel of the manuscript, on replicate-merged inputs.
#
# Layouts follow the published figures so they can be dropped straight in:
#
#   fig3c  Fig. 3C  : one panel per species, m6A sites pooled over tools,
#                     WT vs fip37 KD / Mettl3 KO / IVT
#   fig3d  Fig. 3D  : one panel per species, the three tools the paper showed
#                     (DENA, m6Anet, Nanom6A), WT vs perturbed
#   sup1   Fig. S1  : per-tool metagene grid, three species rows
#   sup7   Fig. S7  : non-m6A chemistries in HeLa, ncRNA metagene per tool
#   sup9   Fig. S9  : RNA004 Dorado other-modification mRNA metagene, WT vs IVT
#
# Inputs are the BEDs of 21b_export_guitar_bed.py; the merge rule is chosen with
# --merge (majority = strictly more than half of the group's independent units).
#
#   Rscript 23c_guitar_grids.R --layout fig3c --merge majority
#   Rscript 23c_guitar_grids.R --layout sup1  --merge majority
#   Rscript 23c_guitar_grids.R --layout sup7
#   Rscript 23c_guitar_grids.R --layout sup9

suppressPackageStartupMessages(library(patchwork))
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(this_file)), "guitar_lib.R"))
arial_setup()
getarg <- function(name, default = NULL) {
  i <- which(commandArgs(TRUE) == name)
  v <- commandArgs(TRUE)
  if (length(i) == 1 && length(v) > i) v[i + 1] else default
}
layout <- getarg("--layout", "fig3c")
merge   <- getarg("--merge", "majority")
FIGDIR  <- file.path(GTR, "figures")

SP <- list(Arabidopsis = c("Arabidopsis_WT", "Arabidopsis_KD", "KD"),
           Mouse       = c("Mouse_WT", "Mouse_KO", "KO"),
           Human       = c("HeLa_WT", "HeLa_IVT", "IVT"))
PLOT_TOOLS <- c("DENA", "m6Anet", "Nanom6A")          # the Fig. 3D trio

# ---------------------------------------------------------------- layouts ----
panels_species <- function(txtype, mod, tools, per_tool, platform = "RNA002") {
  plots <- list()
  for (species in names(SP)) {
    conds <- SP[[species]][1:2]
    if (per_tool) {
      for (tool in tools) {
        beds <- list()
        for (cond in conds) {
          b <- bed_for(platform, merge, cond, mod, tool, 30)
          if (!is.null(b)) names(b) <- cond else { beds <- NULL; break }
          beds[[cond]] <- b[[1]]
        }
        if (is.null(beds) || !length(beds)) next
        plots[[paste(species, tool, sep = " / ")]] <- guitar_metagene(
          beds, txtype, species, title = sprintf("%s  %s", species, tool),
          colors = c(WT_COLOR, CMP_COLOR), area = FALSE)
      }
    } else {
      beds <- list()
      for (cond in conds) {
        b <- bed_for(platform, merge, cond, mod, NULL, 30)
        if (!is.null(b)) beds[[cond]] <- b[[1]]
      }
      if (length(beds) < 2) next
      plots[[species]] <- guitar_metagene(beds, txtype, species,
                                          title = species,
                                          colors = c(WT_COLOR, CMP_COLOR))
    }
  }
  plots
}

fig_from <- function(plots, ncol, out, panel_h) {
  if (!length(plots)) { message("no panels for ", out); return(invisible()) }
  p <- wrap_plots(plots, ncol = ncol)
  save_pair(p, out, width = 4.1 * ncol, height = panel_h * ceiling(length(plots) / ncol))
}

if (layout == "fig3c") {
  fig_from(panels_species("mrna", "m6A", NULL, FALSE), 3,
           file.path(FIGDIR, sprintf("GuitarR_fig3c_%s_mrna.pdf", merge)), 4.4)
  fig_from(panels_species("ncrna", "m6A", NULL, FALSE), 3,
           file.path(FIGDIR, sprintf("GuitarR_fig3c_%s_ncrna.pdf", merge)), 4.4)
} else if (layout == "fig3d") {
  # the published panel D shows three tools per species: draw one panel per
  # species with the six tool x library curves together
  plots <- list()
  for (species in names(SP)) {
    conds <- SP[[species]][1:2]
    beds <- list()
    for (tool in PLOT_TOOLS) for (cond in conds) {
      b <- bed_for("RNA002", merge, cond, "m6A", tool, 30)
      if (!is.null(b)) beds[[paste(tool, cond, sep = " | ")]] <- b[[1]]
    }
    if (length(beds) < 2) next
    cols <- rep(c(WT_COLOR, CMP_COLOR), length.out = length(beds))
    names(cols) <- names(beds)
    plots[[species]] <- guitar_metagene(beds, "mrna", species, title = species,
                                        colors = cols, area = FALSE)
  }
  fig_from(plots, 3, file.path(FIGDIR, sprintf("GuitarR_fig3d_%s_mrna.pdf", merge)), 4.6)
} else if (layout == "sup1") {
  for (species in names(SP)) {
    conds <- SP[[species]][1:2]
    tools <- sort(sub("\\.bed$", "", basename(list.files(
      file.path(BED, "RNA002", merge, conds[1], "m6A"), pattern = "\\.bed$"))))
    plots <- list()
    for (tool in tools) {
      beds <- list()
      for (cond in conds) {
        b <- bed_for("RNA002", merge, cond, "m6A", tool, 30)
        if (!is.null(b)) beds[[cond]] <- b[[1]]
      }
      if (length(beds) < 2) next
      plots[[tool]] <- guitar_metagene(beds, "mrna", species, title = tool,
                                       colors = c(WT_COLOR, CMP_COLOR))
    }
    fig_from(plots, 4, file.path(FIGDIR,
           sprintf("GuitarR_sup1_%s_%s_mrna.pdf", species, merge)), 3.3)
  }
} else if (layout == "sup7") {
  # Fig. S7: non-m6A chemistries in HeLa, ncRNA metagene, WT vs IVT
  plots <- list()
  for (mod in c("Psi", "m1Psi", "m5C", "Nm")) {
    tools <- sub("\\.bed$", "", sort(list.files(
      file.path(BED, "RNA002", merge, "HeLa_WT", mod), pattern = "\\.bed$")))
    for (tool in tools) {
      beds <- list()
      for (cond in c("HeLa_WT", "HeLa_IVT")) {
        b <- bed_for("RNA002", merge, cond, mod, tool, 30)
        if (!is.null(b)) beds[[cond]] <- b[[1]]
      }
      if (length(beds) < 2) next
      plots[[paste(mod, tool, sep = " / ")]] <- guitar_metagene(
        beds, "ncrna", "Human", title = sprintf("%s  %s", mod, tool),
        colors = c(WT_COLOR, CMP_COLOR))
    }
  }
  fig_from(plots, 3, file.path(FIGDIR, sprintf("GuitarR_sup7_%s_ncrna.pdf", merge)), 3.4)
} else if (layout == "sup9") {
  # Fig. S9: RNA004 Dorado other-modification mRNA metagene, WT vs IVT
  plots <- list()
  for (mod in c("Psi", "m5C", "inosine")) {
    tools <- sub("\\.bed$", "", sort(list.files(
      file.path(BED, "RNA004", merge, "RNA004_HeLa_WT", mod), pattern = "\\.bed$")))
    for (tool in tools) {
      beds <- list()
      for (cond in c("RNA004_HeLa_WT", "RNA004_HeLa_IVT")) {
        b <- bed_for("RNA004", merge, cond, mod, tool, 30)
        if (!is.null(b)) beds[[sub("RNA004_", "", cond)]] <- b[[1]]
      }
      if (length(beds) < 2) next
      plots[[paste(mod, tool, sep = " / ")]] <- guitar_metagene(
        beds, "mrna", "Human", title = sprintf("%s  %s", mod, tool),
        colors = c(WT_COLOR, CMP_COLOR))
    }
  }
  fig_from(plots, 3, file.path(FIGDIR, sprintf("GuitarR_sup9_%s_mrna.pdf", merge)), 3.6)
} else stop("unknown --layout: ", layout)
