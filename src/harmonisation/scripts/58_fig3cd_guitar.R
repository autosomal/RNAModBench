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
# Figure 3C / 3D, replicate-aware rebuild (R3-2 / E6, R1-4 / E5, R3-6).
#
# Why this exists.  The published metagene panels come from
# tool_scripts/Guitar/{Arabidopsis,Mouse,Hela}.r, which forced `strand <- "+"` on
# every call and drew per-tool BEDs that carried no replicate structure
# (Arabidopsis = rep3 only, mouse = the mES_WT study only, HeLa = an
# undocumented union).  $RNAMODBENCH_LOCAL/REPLICATE_AWARE_FIGURES.md
# ("What was actually wrong with the published GUITAR figure") documents all
# three, and the original submission's figure code/Figure3/Guitar_*.r is the legacy
# code path being replaced.
#
# This script redraws the two manuscript panels on the replicate-aware BED tree
# written by 21b_export_guitar_bed.py (read-only), and makes the replicate
# structure visible inside the panel: every library is drawn as a thick
# majority-consensus curve plus one thin curve per independent unit (dashed),
# so the reader sees the spread the reviewer asked for (R3-2 / E6).
# The density itself is computed by the Bioconductor Guitar package
# (samplePoints -> normalize -> .generateDensity_CI) -- the same code path as
# the manuscript's own Guitar_*.r scripts, not a re-implementation.
#
#   panel C  one panel per species, m6A sites pooled over the library's tools,
#            WT vs fip37 KD / Mettl3 KO / IVT
#   panel D  one panel per species, the three tools the paper showed
#            (DENA, m6Anet, Nanom6A), WT vs treated
#
# Consensus = sites supported by strictly more than half of the group's
# independent units (bed/RNA002/majority*/... written by 21b; for mouse
# "majority" means the two independent studies agree -- never a replicate mean).
# Panel C also needs a *pooled per-unit* input, which 21b does not write; it is
# rebuilt here as the union of that unit's per-tool BEDs into
# figures/figure3/bed/ (this tree only), and the consensus definition actually
# drawn is audited against that rebuild in
# tables/fig3c_consensus_reconciliation.tsv.
#
# Usage
#   conda run -n guitar_asm --no-capture-output Rscript 58_fig3cd_guitar.R
#   ... --recompute            # ignore the cached density RDS
#   ... --panel C|D|both
#   ... --gate 50              # minimum consensus sites for a thick curve
#
# Outputs -> figures/figure3/
#   figures/Fig3C_metagene_species.{pdf,png}
#   figures/Fig3D_metagene_tools.{pdf,png}
#   tables/fig3c_density.rds, fig3d_density.rds
#   tables/fig3c_guitar_inputs.tsv, fig3d_guitar_inputs.tsv
#   tables/fig3c_consensus_reconciliation.tsv
#
# Fonts: PDF written with grDevices::cairo_pdf while showtext is OFF, so the
# real system Arial is embedded (pdffonts-verifiable); showtext is switched on
# only for the 300 dpi PNG preview.  Nothing below 7.2 pt.

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

GTR   <- file.path(.XB, "harmonisation/guitar_metagene")
BED   <- file.path(GTR, "bed", "RNA002")      # read-only inputs from the audited tree
OUTD  <- file.path(.RB, "figures/figure3")
TABD  <- file.path(OUTD, "tables"); FIGD <- file.path(OUTD, "figures")
MYBED <- file.path(OUTD, "bed", "RNA002")     # pooled per-unit BEDs (ours only)

getarg <- function(name, default = NULL) {
  a <- commandArgs(TRUE); i <- which(a == name)
  if (length(i) == 1 && length(a) > i) a[i + 1] else default
}
recompute <- "--recompute" %in% commandArgs(TRUE)
panel_arg <- getarg("--panel", "both")
merge     <- getarg("--merge", "majority")
GATE      <- as.integer(getarg("--gate", "50"))

## ---- house style: final printed size, Arial, no gridlines --------------------
MM_PER_PT <- as.numeric(grid::convertUnit(grid::unit(1, "pt"), "mm", valueOnly = TRUE))
W <- 6.66                                  # 0.95 x \textwidth (Wiley USG layout)
# row heights chosen so the assembled A-D page matches the printed height of the
# published Figure 3 (~7.8 in at this width); LEGEND_H covers the legend row
PANEL_H <- 1.41; LETTER_H <- 0.13
FS <- list(axis = 7.4, ytitle = 7.4, title = 8.6, legend = 7.2, struct = 7.2)
stopifnot(min(unlist(FS)) >= 7.2)

SPECIES <- list(
  Arabidopsis = list(conds = c(Arabidopsis_WT = "WT", Arabidopsis_KD = "KD"), order = 1),
  Mouse       = list(conds = c(Mouse_WT = "WT", Mouse_KO = "KO"), order = 2),
  Human       = list(conds = c(HeLa_WT = "WT", HeLa_IVT = "IVT"), order = 3))
TOOLS_D <- c("DENA", "m6Anet", "Nanom6A")          # the published Fig. 3D trio
# published panel-D grammar: WT warm, treated cool, one colour per tool x condition
# layout v2: panel D encodes condition by colour (WT / treated) and tool by
# linetype; the tool key is carried by the figure-wide legend band
TOOL_LINETYPES <- c(DENA = "solid", m6Anet = "dashed", Nanom6A = "dotted")
WT_COLOR <- "#4b81b8"; CMP_COLOR <- "#e8a76b"      # panel C = the manuscript palette

n_lines <- function(p) if (file.exists(p)) length(readLines(p, warn = FALSE)) else NA_integer_

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

unit_tags <- function(cond) {
  d <- list.files(BED, pattern = "^rep_", full.names = FALSE)
  d <- d[vapply(d, function(x) dir.exists(file.path(BED, x, cond)), logical(1))]
  sort(sub("^rep_", "", d))
}

## ---- panel C: pooled per-unit BED = union of that unit's per-tool BEDs -------
pooled_unit_bed <- function(cond, tag) {
  out <- file.path(MYBED, paste0("rep_", tag), cond, "m6A", "ALL.bed")
  if (file.exists(out)) return(out)
  d <- file.path(BED, paste0("rep_", tag), cond, "m6A")
  files <- list.files(d, pattern = "\\.bed$", full.names = TRUE)
  if (!length(files)) return(NULL)
  lines <- unique(unlist(lapply(files, readLines, warn = FALSE), use.names = FALSE))
  lines <- lines[nzchar(lines)]
  if (!length(lines)) return(NULL)
  dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
  writeLines(lines, out)
  out
}

## ---- the planned curves of one panel -----------------------------------------
build_plan <- function(panel) {
  rows <- list()
  for (sp in names(SPECIES)) {
    conds <- SPECIES[[sp]]$conds
    addrow <- function(kind, unit_tag, path, tool, cond, role, n) {
      rows[[length(rows) + 1]] <<- data.frame(
        panel = panel, species = sp, tool = tool, cond_group = cond, role = role,
        kind = kind, unit_tag = unit_tag, path = path,
        n_sites = if (is.na(n)) -1L else as.integer(n), stringsAsFactors = FALSE)
    }
    if (panel == "C") {
      for (i in seq_along(conds)) {
        cond <- names(conds)[i]; role <- if (i == 1) "WT" else "treated"
        f <- file.path(BED, paste0(merge, "_pooled"), cond, "m6A", "ALL.bed")
        addrow("consensus", "", f, "ALL", cond, role, n_lines(f))
        for (tag in unit_tags(cond)) {
          f <- pooled_unit_bed(cond, tag)
          if (is.null(f)) next
          addrow("unit", tag, f, "ALL", cond, role, n_lines(f))
        }
      }
    } else {
      for (tool in TOOLS_D) for (i in seq_along(conds)) {
        cond <- names(conds)[i]; role <- if (i == 1) "WT" else "treated"
        f <- file.path(BED, merge, cond, "m6A", paste0(tool, ".bed"))
        addrow("consensus", "", f, tool, cond, role, n_lines(f))
        for (tag in unit_tags(cond)) {
          f <- file.path(BED, paste0("rep_", tag), cond, "m6A", paste0(tool, ".bed"))
          if (!file.exists(f)) next
          addrow("unit", tag, f, tool, cond, role, n_lines(f))
        }
      }
    }
  }
  out <- do.call(rbind, rows)
  out$group <- paste(out$panel, out$species, out$tool, out$role, out$kind,
                     out$unit_tag, sep = "|")
  out
}

## ---- Guitar density per species, with per-group error isolation --------------
guitar_density <- function(plan, txtype = "mrna") {
  dens_all <- list()
  for (sp in unique(plan$species)) {
    sub <- plan[plan$species == sp, ]
    gt  <- guitar_txdb(sp, txtype)
    beds <- as.list(sub$path); names(beds) <- sub$group
    sitesGroup <- Guitar:::.getStGroup(stBedFiles = unname(beds), stGroupName = names(beds))
    relative <- list(); weight <- list(); errs <- character(0)
    for (gname in names(sitesGroup)) {
      ok <- tryCatch({
        sp2 <- Guitar::samplePoints(sitesGroup[gname], stSampleNum = 3, stAmblguity = 5,
                                    pltTxType = txtype, stSampleModle = "Equidistance",
                                    mapFilterTranscript = TRUE, gt)
        nz <- Guitar::normalize(sp2, gt, txtype, 1, 1)
        relative[[gname]] <- nz[[1]]; weight[[gname]] <- nz[[2]]; TRUE
      }, error = function(e) {
        errs <<- c(errs, sprintf("%s: %s", gname, conditionMessage(e))); FALSE
      })
      if (!ok) message("!! samplePoints failed: ", gname)
    }
    if (length(errs)) message("failed groups (", length(errs), ")")
    if (!length(relative)) { message("!! no curves for ", sp); next }
    d <- Guitar:::.generateDensity_CI(relative, weight, CI_ResamplingTime = 20,
                                      adjust = 1, enableCI = FALSE)
    d$species <- sp
    d$componentWidth <- list(gt[[txtype]]$componentWidthAverage_pct)
    dens_all[[sp]] <- d
  }
  dens_all
}

compute_or_load <- function(panel) {
  cache <- file.path(TABD, sprintf("fig3%s_density.rds", tolower(panel)))
  if (file.exists(cache) && !recompute) {
    message("using cached density: ", basename(cache)); return(readRDS(cache))
  }
  plan <- build_plan(panel)
  plan <- plan[file.exists(plan$path), ]
  plan <- plan[!(plan$kind == "consensus" & (is.na(plan$n_sites) | plan$n_sites < GATE)), ]
  if (!nrow(plan)) stop("no input BED for panel ", panel)
  dens <- guitar_density(plan)
  res <- list(dens = dens, plan = plan, panel = panel, gate = GATE, merge = merge)
  dir.create(TABD, showWarnings = FALSE, recursive = TRUE)
  saveRDS(res, cache); message("cached -> ", cache)
  res
}

## ---- transcript schematic (compact: 7.2 pt labels in a ~2 in panel) ----------
add_structure <- function(p, comp_width, pos) {
  map_lab   <- c(promoter = "1kb", utr5 = "5'UTR", cds = "CDS", utr3 = "3'UTR",
                 tail = "1kb", ncrna = "ncRNA", rna = "RNA")
  map_h     <- c(promoter = 0.002, utr5 = 0.005, cds = 0.010, utr3 = 0.005,
                 tail = 0.002, ncrna = 0.010, rna = 0.010)
  map_alpha <- c(promoter = 0.99, utr5 = 0.99, cds = 0.272, utr3 = 0.99,
                 tail = 0.99, ncrna = 0.20, rna = 0.20)
  nm <- names(comp_width); end <- cumsum(as.numeric(comp_width))
  sta <- c(0, end[-length(end)]) + 0.001; mid <- (sta + end) / 2
  pk <- pos$fig_top / 1.05
  h <- unname(map_h[nm]) * pk
  xlab <- pmin(pmax(mid, 0.045), 0.955)
  for (i in seq_along(nm)) {
    p <- p + annotate("rect", xmin = sta[i], xmax = end[i],
                      ymin = pos$rna_lgd_bl - h[i], ymax = pos$rna_lgd_bl + h[i],
                      fill = "grey35", colour = NA, alpha = unname(map_alpha[nm[i]]))
    p <- p + annotate("text", x = xlab[i], y = pos$rna_comp_text,
                      label = unname(map_lab[nm[i]]), family = "Arial",
                      size = FS$struct * MM_PER_PT)
  }
  p
}

## ---- one metagene panel ------------------------------------------------------
draw_panel <- function(dens, plan, title, cols, col_labels, lin_lab,
                       tool_linetypes = NULL) {
  meta <- plan[, c("group", "role", "kind", "unit_tag", "cond_group", "tool", "n_sites")]
  # componentWidth is a list column -- keep it off the merged frame
  d <- merge(dens[setdiff(names(dens), "componentWidth")], meta, by = "group",
             all.x = TRUE)
  d$ckey <- if (is.null(tool_linetypes) && !all(d$tool == "ALL"))
               paste(d$tool, d$role, sep = " | ")
             else ifelse(d$role == "WT", "WT", "treated")
  d$curvetype <- if (!is.null(tool_linetypes)) d$tool
                 else ifelse(d$kind == "consensus", lin_lab[1], lin_lab[2])
  peak <- max(d$density, na.rm = TRUE)
  pos <- Guitar:::.generate_pos_para(peak)
  pos$fig_bottom <- -0.18 * peak
  pos$rna_comp_text <- -0.082 * peak
  cons <- d[d$kind == "consensus", ]; units <- d[d$kind == "unit", ]

  p <- ggplot(d, aes(x = x))
  if (nrow(units))
    p <- p + geom_line(data = units, aes(y = density, group = group, colour = ckey,
                                         linetype = curvetype),
                       linewidth = 0.28, alpha = 0.38)
  if (nrow(cons)) {
    p <- p + geom_ribbon(data = cons, aes(ymin = 0, ymax = density, group = group,
                                          fill = ckey), alpha = 0.16, colour = NA)
    p <- p + geom_line(data = cons, aes(y = density, group = group, colour = ckey,
                                        linetype = curvetype), linewidth = 0.8)
  }
  p <- p +
    scale_colour_manual(values = cols, limits = names(cols), labels = col_labels) +
    scale_fill_manual(values = cols, limits = names(cols)) +
    scale_linetype_manual(values = if (is.null(tool_linetypes))
                                      setNames(c("solid", "dashed"), lin_lab)
                                    else tool_linetypes,
                          limits = if (is.null(tool_linetypes)) lin_lab
                                   else names(tool_linetypes)) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(pos$fig_bottom, pos$fig_top), expand = c(0, 0)) +
    labs(x = NULL, y = "Density", title = title) +
    guides(colour = "none", fill = "none", linetype = "none", linewidth = "none") +
    theme_classic(base_family = "Arial", base_size = FS$axis) +
    theme(panel.grid = element_blank(),
          panel.background = element_rect(fill = "white", colour = NA),
          axis.line = element_line(colour = "black", linewidth = 0.3),
          axis.ticks.x = element_blank(), axis.text.x = element_blank(),
          axis.title.y = element_text(size = FS$ytitle, colour = "black"),
          axis.text.y = element_text(size = FS$axis, colour = "black"),
          legend.position = "none",
          legend.title = element_blank(),
          legend.text = element_text(size = FS$legend),
          legend.key.width = unit(9, "pt"), legend.key.height = unit(4, "pt"),
          legend.margin = margin(0, 0, 0, 0, "pt"),
          legend.box.margin = margin(0, 0, 0, 0, "pt"),
          legend.box = "horizontal",
          plot.title = element_text(size = FS$title, face = "bold", hjust = 0.5),
          plot.margin = margin(3, 5, 2, 3, "pt"),
          text = element_text(family = "Arial"))
  add_structure(p, dens$componentWidth[[1]], pos)
}

## ---- one panel row -----------------------------------------------------------
draw_row <- function(res, letter) {
  panel <- res$panel
  panels <- list()
  # legend labels are identical across species on purpose: patchwork can only
  # collect-merge guides whose labels match, and three unmerged legend blocks do
  # not fit the 6.66 in canvas (they were clipped).  The condition name travels
  # in the panel title instead, exactly as in panel B.
  lin_lab <- c("majority consensus", "individual replicate or study")
  for (sp in names(SPECIES)) {
    dens <- res$dens[[sp]]; plan <- res$plan[res$plan$species == sp, ]
    if (is.null(dens) || !nrow(plan)) { message("!! no density for ", sp); next }
    # species name only: the condition mapping lives in the caption (layout v2)
    title <- sp
    if (panel == "C") {
      cols <- c(WT = unname(WT_COLOR), treated = unname(CMP_COLOR))
      labs <- c(WT = "WT", treated = "treated")
      panels[[sp]] <- draw_panel(dens, plan, title, cols, labs, lin_lab)
    } else {
      # panel D: colour = condition, linetype = tool -- the tool key travels in
      # the single legend band under the figure (61_fig3_legend_band.py)
      cols <- c(WT = unname(WT_COLOR), treated = unname(CMP_COLOR))
      labs <- c(WT = "WT", treated = "treated")
      panels[[sp]] <- draw_panel(dens, plan, title, cols, labs, lin_lab,
                                 tool_linetypes = TOOL_LINETYPES)
    }
  }
  if (!length(panels)) stop("no panels drawn for ", panel)
  inner <- wrap_plots(panels, ncol = 3) +
    plot_layout() &
    theme(legend.position = "none")
  letter_p <- ggplot() +
    annotate("text", x = 0, y = 0.5, label = letter, family = "Arial",
             fontface = "bold", size = FS$title * MM_PER_PT, hjust = 0, vjust = 0.5) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme_void()
  wrap_elements(full = letter_p) / inner +
    plot_layout(heights = c(LETTER_H, PANEL_H))
}

## ---- panel-C audit: is the drawn consensus what the units support? -----------
reconcile_pooled <- function(plan) {
  f <- function(cond, tag, role) {
    p <- plan[plan$tool == "ALL" & plan$cond_group == cond & plan$role == role, ]
    list(consensus = p$n_sites[p$kind == "consensus"],
         units = setNames(p$path[p$kind == "unit"],
                          p$unit_tag[p$kind == "unit"]))
  }
  rows <- list()
  for (sp in names(SPECIES)) {
    conds <- SPECIES[[sp]]$conds
    for (i in seq_along(conds)) {
      cond <- names(conds)[i]; role <- if (i == 1) "WT" else "treated"
      x <- f(cond, NULL, role)
      if (!length(x$units)) next
      keys <- lapply(names(x$units), function(tag) {
        l <- readLines(x$units[[tag]], warn = FALSE)
        sub("^([^\t]+\t[^\t]+)\t.*$", "\\1", l)
      })
      names(keys) <- names(x$units)
      # support of every site across the drawn per-unit pooled curves
      tab <- table(unlist(keys, use.names = FALSE))
      n_units <- length(keys)
      rows[[length(rows) + 1]] <- data.frame(
        species = sp, cond_group = cond, role = role,
        n_units_drawn = n_units,
        quorum = n_units %/% 2 + 1,
        n_sites_drawn_consensus = x$consensus,
        n_sites_majority_of_drawn_units = sum(tab >= (n_units %/% 2 + 1)),
        n_sites_union_of_drawn_units = length(tab),
        consensus_definition = paste("21b", merge, "_pooled ALL.bed",
                                     "(union over tools of each tool's majority set)"),
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}

## ---- run --------------------------------------------------------------------
dir.create(TABD, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGD, showWarnings = FALSE, recursive = TRUE)

panels_to_do <- if (panel_arg == "both") c("C", "D") else toupper(panel_arg)
res_list <- list()
for (panel in panels_to_do) {
  res <- compute_or_load(panel)
  res_list[[panel]] <- res
  row_letter <- panel
  p <- draw_row(res, row_letter)
  save_pair_cairo(p, file.path(FIGD, sprintf("Fig3%s_metagene_%s.pdf", row_letter,
                                             if (panel == "C") "species" else "tools")),
                  width = W, height = LETTER_H + PANEL_H)
  plan <- res$plan
  write.table(plan[, c("panel", "species", "tool", "cond_group", "role", "kind",
                       "unit_tag", "n_sites", "path")],
              file.path(TABD, sprintf("fig3%s_guitar_inputs.tsv", tolower(panel))),
              sep = "\t", quote = FALSE, row.names = FALSE)
  message(sprintf("panel %s: %d curves planned", panel, nrow(plan)))
}

if ("C" %in% names(res_list)) {
  rec <- reconcile_pooled(res_list$C$plan)
  write.table(rec, file.path(TABD, "fig3c_consensus_reconciliation.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  message("wrote fig3c_consensus_reconciliation.tsv")
}
message("done ", format(Sys.time()))
