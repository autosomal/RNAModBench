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

# Shared helpers for redrawing the manuscript's GUITAR panels with the R
# Bioconductor package Guitar, on replicate-merged inputs.
#
# Guitar cannot merge replicates: stSampleNum only sets how many equidistant
# points are taken inside each (padded) site interval, and groups passed through
# stGroupName are pooled by concatenation.  The merge therefore happens in
# 21b_export_guitar_bed.py, which writes
#
#   bed/<platform>/<merge>/<Condition>/<mod>/<Tool>.bed
#   bed/<platform>/<merge>_pooled/<Condition>/<mod>/ALL.bed
#   bed/<platform>/rep_<tag>/<Condition>/<mod>/<Tool>.bed
#
# with <merge> in {union, majority, intersection}.

suppressPackageStartupMessages({
  library(Guitar)
  library(GenomicFeatures)
  library(ggplot2)
  library(scales)
  library(showtext)
})

GTR <- file.path(.XB, "harmonisation/guitar_metagene")
BED <- file.path(GTR, "bed")

GTF <- c(
  Arabidopsis = file.path(.XB, "reference/arabidopsis/Arabidopsis_thaliana.TAIR10.61.gtf"),
  Mouse       = file.path(.XB, "reference/GRCm39/ensembl/Mus_musculus.GRCm39.114.gtf"),
  Human       = file.path(.XB, "reference/GRCh38p14/ensembl112/Homo_sapiens.GRCh38.112.chr.gtf"),
  Curlcake    = file.path(.XB, "reference/curlcakes/Curlcake.gtf"))

WT_COLOR <- "#4b81b8"; CMP_COLOR <- "#e8a76b"

arial_setup <- function() {
  font_add("Arial",
           regular    = "/usr/share/fonts/truetype/msttcorefonts/Arial.ttf",
           bold       = "/usr/share/fonts/truetype/msttcorefonts/arialbd.ttf",
           italic     = "/usr/share/fonts/truetype/msttcorefonts/Arial_Italic.ttf",
           bolditalic = "/usr/share/fonts/truetype/msttcorefonts/Arial_Bold_Italic.ttf")
  showtext_auto()
  showtext_opts(dpi = 300)
}

#' Cached Guitar transcript model.  A TxDb from makeTxDbFromGFF is backed by an
#' in-memory database and cannot be restored from disk, so the derived
#' guitarTxdb (plain data) is what gets cached.
guitar_txdb <- function(species, txtype = "mrna") {
  dir.create(file.path(GTR, "txdb_cache"), showWarnings = FALSE, recursive = TRUE)
  cache <- file.path(GTR, "txdb_cache", sprintf("%s.%s.guitarTxdb.rds", species, txtype))
  if (file.exists(cache)) return(readRDS(cache))
  message("building guitarTxdb: ", species, " ", txtype,
          " (slow for human/mouse; cached afterwards)")
  txdb <- makeTxDbFromGFF(file = GTF[[species]], format = "auto")
  gt <- Guitar::makeGuitarTxdb(
    txdb, txfiveutrMinLength = 100, txcdsMinLength = 100,
    txthreeutrMinLength = 100, txlongNcrnaMinLength = 100,
    txlncrnaOverlapmrna = FALSE, txpromoterLength = 1000, txtailLength = 1000,
    txAmblguity = 5, txPrimaryOnly = FALSE, pltTxType = txtype)
  saveRDS(gt, cache)
  gt
}

#' The part of GuitarPlot() that follows the transcript model, so a cached
#' guitarTxdb can be reused across panels (GuitarPlot's own txGuitarTxdb
#' argument expects a read.table-able file, not an object).
guitar_sites <- function(bed_files, txtype, gt, stSampleNum = 3,
                        mapFilterTranscript = TRUE) {
  sitesGroup <- Guitar:::.getStGroup(stBedFiles = bed_files,
                                     stGroupName = names(bed_files))
  relative <- list(); weight <- list()
  for (gname in names(sitesGroup)) {
    sp <- Guitar::samplePoints(sitesGroup[gname], stSampleNum = stSampleNum,
                              stAmblguity = 5, pltTxType = txtype,
                              stSampleModle = "Equidistance",
                              mapFilterTranscript = mapFilterTranscript, gt)
    nz <- Guitar::normalize(sp, gt, txtype, 1, 1)
    relative[[txtype]][[gname]] <- nz[[1]]
    weight[[txtype]][[gname]] <- nz[[2]]
  }
  list(relative = relative[[txtype]], weight = weight[[txtype]])
}

#' Build the metagene ggplot for one set of groups.
guitar_metagene <- function(bed_files, txtype = "mrna", species = "Human",
                           title = "", colors = NULL, area = TRUE) {
  stopifnot(length(bed_files) >= 1, !is.null(names(bed_files)))
  gt <- guitar_txdb(species, txtype)
  s <- guitar_sites(bed_files, txtype, gt)
  dens <- Guitar:::.generateDensity_CI(s$relative, s$weight,
                                       CI_ResamplingTime = 1000, adjust = 1,
                                       enableCI = FALSE)
  p <- Guitar:::.plotDensity_CI(
    dens, componentWidth = gt[[txtype]]$componentWidthAverage_pct,
    headOrtail = TRUE, title = title, enableCI = FALSE)

  if (is.null(colors)) colors <- rep(NA_character_, length(bed_files))
  names(colors) <- names(bed_files)
  if (!area) {
    keep <- !vapply(p$layers, function(l) inherits(l$geom, "GeomRibbon"),
                    logical(1))
    p$layers <- p$layers[keep]
  }
  p + ggtitle(title) +
    theme_bw(base_family = "Arial", base_size = 11) +
    theme(panel.grid = element_blank(),
          panel.background = element_rect(fill = "white", colour = NA),
          panel.border = element_blank(),
          axis.line = element_line(colour = "black", linewidth = 0.6),
          axis.text = element_text(colour = "black"),
          axis.title = element_blank(),
          legend.position = "none",
          plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
          plot.margin = margin(4, 8, 4, 4)) +
    {if (any(!is.na(colors)))
       list(scale_colour_manual(values = colors, limits = names(bed_files)),
            scale_fill_manual(values = alpha(colors, 0.28),
                              limits = names(bed_files)))}
}

#' One normalised BED per group, from the tree written by 21b.
bed_for <- function(platform, merge, group, mod, tool = NULL, min_rows = 1) {
  d <- file.path(BED, platform, if (is.null(tool)) paste0(merge, "_pooled") else merge,
                 group, mod)
  if (is.null(tool)) {
    f <- file.path(d, "ALL.bed")
    return(if (file.exists(f) && length(readLines(f)) >= min_rows)
      setNames(list(f), paste0(group, " | ", merge)) else NULL)
  }
  f <- file.path(d, paste0(tool, ".bed"))
  if (!file.exists(f) || length(readLines(f)) < min_rows) return(NULL)
  setNames(list(f), paste0(tool, " | ", group))
}

save_pair <- function(p, out, width, height) {
  dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
  ggsave(out, plot = p, width = width, height = height, units = "in")
  ggsave(sub("\\.pdf$", ".png", out), plot = p, width = width, height = height,
         units = "in", dpi = 300)
  message("wrote ", out)
}

## ------------------------------------------------- region labels (S1 / 7D) --- #
## Single-line region labels with collision avoidance and leader lines.
##

## CDS / 3'UTR / 1kb) must share one baseline instead of the two staggered rows
## S1 used before.  The five strings are wider than an S1 panel cell, so the row
## is packed with a minimum gap, centred on its anchors, and every label that had
## to move gets a thin leader back to its own segment (the genome-browser
## convention).  Full text and a size at or above the 7 pt floor are kept.
##
## Ink widths are taken from the delivered figure (pdftotext -bbox, 8.5 pt
## Arial): 1kb 13.7, CDS 17.7, 5'UTR 24.1, 3'UTR 24.1 pt.  R-side strwidth()
## does not see the showtext metrics, so widths are scaled from that measurement
## rather than re-estimated.
REGION_INK_PT <- c("1kb" = 13.7, "CDS" = 17.7, "5'UTR" = 24.1, "3'UTR" = 24.1)

## ggplot's size unit is mm; the callers define the same constant from grid, this
## keeps the helper usable on its own (the caller's value wins when it exists)
if (!exists("MM_PER_PT")) MM_PER_PT <- 25.4 / 72.27

region_ink_pt <- function(labels, size, ink_ref_pt = REGION_INK_PT,
                          ref_size = 8.5) {
  w <- unname(ink_ref_pt[labels])
  fallback <- is.na(w)                      # ncRNA / RNA and future names
  w[fallback] <- 0.60 * ref_size * nchar(labels[fallback])
  w * size / ref_size
}

## The size ladder: full text first, then closer packing, then smaller type --
## never below the 7 pt floor and never dropping a label.
REGION_LADDER <- data.frame(size = c(8.5, 8.0, 7.5, 7.5, 7.0),
                            gap  = c(3.0, 2.5, 2.5, 2.0, 2.0))


## They are a deliberate exception to the 7 pt floor -- that one row has to fit a
## 95.76 pt panel cell (five labels are 93.3 pt wide at 8.5 pt), and the author
## chose smaller type over moving the labels off their segments.  Pinning the
## tier keeps a rerun from silently stepping back up; REGION_LADDER stays for
## other rows that may appear later.
## The size is not a taste question: the packed row has to stay inside the panel
## cell with a readable white margin on both sides.  region_size_for() returns
## the largest candidate whose row (labels + gaps) fits that budget, and the tier
## below is derived from it once so a rerun cannot silently drift.
region_size_for <- function(window_pt, labels, gap_pt = 3.0, margin_pt = 8.0,
                            candidates = c(8.5, 8.0, 7.5, 7.0, 6.5, 6.0, 5.5,
                                           5.0, 4.5),
                            ink_ref_pt = REGION_INK_PT, ref_size = 8.5) {
  budget <- window_pt - 2 * margin_pt
  for (size in candidates) {                     # candidates run large -> small
    row <- sum(region_ink_pt(labels, size, ink_ref_pt, ref_size)) +
      gap_pt * (length(labels) - 1)
    if (row <= budget) {
      return(list(size = size, gap = gap_pt, row_pt = row,
                  margin_each_pt = (window_pt - row) / 2))
    }
  }
  stop(sprintf(paste0("region labels do not fit %.1f pt with %.1f pt of margin ",
                      "on each side: the smallest candidate still needs %.1f pt"),
               window_pt, margin_pt,
               sum(region_ink_pt(labels, min(candidates), ink_ref_pt, ref_size)) +
                 gap_pt * (length(labels) - 1)))
}

## derived for the S1 cell (95.76 pt) with 8 pt of white on each side
## The rendered row measures ~1.2x the packed width (the panel's real plot area
## is wider than the 50.3 pt used for the data-unit conversion), so the white
## margin budget is raised to 13 pt: measured on the 2026-09-23 build, this is
## what keeps the outer "1kb" labels inside their own panel.
S1_REGION <- region_size_for(95.76, c("1kb", "5'UTR", "CDS", "3'UTR", "1kb"),
                             margin_pt = 15)
REGION_SIZE_PT <- S1_REGION$size
REGION_GAP_PT  <- S1_REGION$gap
REGION_TIER    <- data.frame(size = REGION_SIZE_PT, gap = REGION_GAP_PT)

fit_region_size <- function(labels, allow_w_pt, ladder = REGION_LADDER,
                            ink_ref_pt = REGION_INK_PT, ref_size = 8.5,
                            margin_pt = 1.5) {
  ## margin_pt: the row has to clear the panel cell by a little, otherwise the
  ## edge clamp below eats into the minimum gap on the outer pairs
  for (i in seq_len(nrow(ladder))) {
    need <- sum(region_ink_pt(labels, ladder$size[i], ink_ref_pt, ref_size)) +
      ladder$gap[i] * (length(labels) - 1)
    if (need <= allow_w_pt - margin_pt) {
      return(list(size = ladder$size[i], gap = ladder$gap[i], need_pt = need))
    }
  }
  NULL
}

## x_anchor : one anchor per label, in x data units (segments' midpoints)
## allow_data: c(min, max) the label row may occupy, in x data units (the panel
##             cell, not the axis: the row may overflow the plot area on purpose)
## y        : the single baseline, in y data units
## anchor_y : leader start, in y data units (segment bar underside); NULL = none
## leader_drop: how far below the baseline a leader should stop (y data units)
## mode = "segments" (Figure 8E look, author decision 2026-09-23): every label sits on its
## own segment anchor, the labels named in `rotate` are drawn vertically so their
## horizontal footprint becomes the text height instead of the text width, and no
## leader is drawn.  Only a micro-nudge of at most `nudge_pt` keeps neighbours
## apart -- anything larger is reported as a failure instead of being hidden.
region_labels <- function(x_anchor, labels, y, axes_w_pt, allow_data,
                          size, gap, anchor_y = NULL, leader_drop = 0,
                          mode = c("packed", "segments"),
                          rotate = character(0), nudge_pt = 2.0,
                          ink_ref_pt = REGION_INK_PT, ref_size = 8.5,
                          leader_pt = 1.5, colour = "black",
                          linewidth = 0.25) {
  stopifnot(length(x_anchor) == length(labels), axes_w_pt > 0)
  mode <- match.arg(mode)
  # a vertical label occupies its text height (Arial ~1.15 em) instead of its width
  ink <- region_ink_pt(labels, size, ink_ref_pt, ref_size)
  ink[labels %in% rotate] <- 1.15 * size
  wd <- ink / axes_w_pt                                              # data units
  if (mode == "segments") {
    x <- x_anchor
    for (pass in 1:3) {                       # keep neighbours apart, <= nudge_pt
      for (i in seq_len(length(x) - 1)) {
        need <- (wd[i] + wd[i + 1]) / 2 + gap / axes_w_pt
        if (x[i + 1] - x[i] < need) {
          push <- (need - (x[i + 1] - x[i])) / 2
          x[i] <- x[i] - push
          x[i + 1] <- x[i + 1] + push
        }
      }
      x <- pmin(pmax(x, allow_data[1] + wd / 2), allow_data[2] - wd / 2)
    }
    moved <- abs(x - x_anchor) * axes_w_pt
    txt <- lapply(seq_along(labels), function(i)
      annotate("text", x = x[i], y = y, label = labels[i], family = "Arial",
               size = size * MM_PER_PT, colour = colour,
               angle = if (labels[i] %in% rotate) 90 else 0))
    lo <- x - wd / 2; hi <- x + wd / 2
    fits <- all(lo >= allow_data[1] - 1e-6) && all(hi <= allow_data[2] + 1e-6) &&
      (length(x) < 2 || all(lo[-1] - hi[-length(hi)] >= gap / axes_w_pt - 1e-6)) &&
      max(moved) <= nudge_pt
    return(list(text = txt, leaders = list(),
                placed = data.frame(label = labels, anchor = x_anchor, x = x,
                                    moved_pt = moved, ink_pt = ink),
                size = size, gap = gap, n_leaders = 0L, fits = fits,
                max_moved_pt = max(moved, 0)))
  }
  x <- x_anchor
  gap_d <- gap / axes_w_pt
  for (pass in 1:3) {
    for (i in seq_len(length(x) - 1)) {          # push right where they collide
      if (x[i + 1] - x[i] < (wd[i] + wd[i + 1]) / 2 + gap_d) {
        x[i + 1] <- x[i] + (wd[i] + wd[i + 1]) / 2 + gap_d
      }
    }
    # pull the row down to its minimum width when the anchors are spread wider
    # than the packed row: the row, not the anchors, has to fit the cell budget
    # (without this the outer labels ran past the panel and the neighbouring
    #  panels' "1kb" labels overlapped, 2026-09-23)
    for (i in rev(seq_len(length(x) - 1))) {
      need <- (wd[i] + wd[i + 1]) / 2 + gap_d
      if (x[i + 1] - x[i] > need) {
        get_mid <- (x[i] + x[i + 1]) / 2
        x[i] <- get_mid - need / 2
        x[i + 1] <- get_mid + need / 2
      }
    }
    x <- x + (mean(x_anchor) - mean(x))          # keep the row on its anchors
    x <- pmin(pmax(x, allow_data[1] + wd / 2), allow_data[2] - wd / 2)
  }
  moved_pt <- abs(x - x_anchor) * axes_w_pt
  txt <- lapply(seq_along(labels), function(i) {
    annotate("text", x = x[i], y = y, label = labels[i], family = "Arial",
             size = size * MM_PER_PT, colour = colour)
  })
  lead <- list()
  if (!is.null(anchor_y)) {
    for (i in which(moved_pt > leader_pt)) {
      lead[[length(lead) + 1]] <-
        annotate("segment", x = x_anchor[i], xend = x[i],
                 y = anchor_y, yend = y - leader_drop,
                 linewidth = linewidth, colour = colour)
    }
  }
  lo <- x - wd / 2
  hi <- x + wd / 2
  fits <- all(lo >= allow_data[1] - 1e-6) && all(hi <= allow_data[2] + 1e-6) &&
    (length(x) < 2 || all(lo[-1] - hi[-length(hi)] >= gap_d - 1e-6))
  list(text = txt, leaders = lead,
       placed = data.frame(label = labels, anchor = x_anchor, x = x,
                           moved_pt = moved_pt, ink_pt = wd * axes_w_pt),
       size = size, gap = gap, n_leaders = length(lead), fits = fits,
       max_moved_pt = max(moved_pt, 0))
}
