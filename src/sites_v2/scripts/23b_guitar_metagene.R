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
# GUITAR metagene, redrawn with the same R package the manuscript used
# (Bioconductor Guitar), on replicate-MERGED inputs.
#
# Why this exists.  The published panels were produced from
# code/guitar/<Condition>_clean/<Tool>.bed.  Comparing those files with the
# rebuilt per-replicate call sets shows they were not replicate-aware:
#   Arabidopsis  rep3 only          (EpiNano 5,833 = rep3 5,833; DRUMMER 339 = 339)
#   Mouse        mES_WT study only  (EpiNano 33,073 = 33,073); the 2nd study unused
#   HeLa         a union of 3 replicates (CHEUI 65,550 vs 68,597 summed)
# Guitar's `stSampleNum` is not a replicate parameter either: samplePoints() does
# `stSampleNum <- 2*stSampleNum-1` and takes that many EQUIDISTANT points inside
# every site interval, i.e. it only spreads a padded site across its own width.
# With enableCI = FALSE in the original scripts, the figure carried no replicate
# information at all.
#
# Every Guitar variant is drawn in ONE GuitarPlot call by passing several groups
# (building the transcript model dominates the runtime), so a species costs a
# handful of calls rather than dozens.
#
# Inputs come from 21b_export_guitar_bed.py:
#   bed/<merge>/<Condition>/<Tool>.bed          one file per tool
#   bed/<merge>_pooled/<Condition>/ALL.bed      all tools of a condition pooled
#   bed/rep_<tag>/<Condition>/<Tool>.bed        a single replicate
#
# Usage
#   Rscript 23b_guitar_metagene.R --species Arabidopsis --mode condition --merge majority
#   Rscript 23b_guitar_metagene.R --species Mouse --mode condition --merge majority,union,intersection
#   Rscript 23b_guitar_metagene.R --species Human --mode pertool --merge majority
#   Rscript 23b_guitar_metagene.R --species Human --mode replicates --merge majority
# Options: --txtype mrna|ncrna  --min-sites N  --out FILE.pdf  --keep-bed DIR

suppressPackageStartupMessages({
  library(Guitar)
  library(GenomicFeatures)
  library(rtracklayer)
  library(ggplot2)
  library(scales)
  library(showtext)
})

## ---- house style: all English, Arial, no in-panel gridlines ---------------
font_add("Arial",
         regular    = "/usr/share/fonts/truetype/msttcorefonts/Arial.ttf",
         bold       = "/usr/share/fonts/truetype/msttcorefonts/arialbd.ttf",
         italic     = "/usr/share/fonts/truetype/msttcorefonts/Arial_Italic.ttf",
         bolditalic = "/usr/share/fonts/truetype/msttcorefonts/Arial_Bold_Italic.ttf")
showtext_auto()
showtext_opts(dpi = 300)

GTR <- file.path(.XB, "sites_v2/guitar_metagene")
GTF <- c(
  Arabidopsis = file.path(.XB, "reference/arabidopsis/Arabidopsis_thaliana.TAIR10.61.gtf"),
  Mouse       = file.path(.XB, "reference/GRCm39/ensembl/Mus_musculus.GRCm39.114.gtf"),
  Human       = file.path(.XB, "reference/GRCh38p14/ensembl112/Homo_sapiens.GRCh38.112.chr.gtf"))
CONDS <- list(Arabidopsis = c("Arabidopsis_WT", "Arabidopsis_KD"),
              Mouse       = c("Mouse_WT", "Mouse_KO"),
              Human       = c("HeLa_WT", "HeLa_IVT"))
WT_COLOR <- "#4b81b8"; CMP_COLOR <- "#e8a76b"

## ---- CLI ------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
getarg <- function(name, default = NULL) {
  i <- which(args == name)
  if (length(i) >= 1 && length(args) > i[1]) args[i[1] + 1] else default
}
species <- getarg("--species", "Arabidopsis")
merges  <- strsplit(getarg("--merge", "majority"), ",")[[1]]
mode    <- getarg("--mode", "condition")     # condition | pertool | replicates
txtype  <- getarg("--txtype", "mrna")        # mrna | ncrna
minsit  <- as.integer(getarg("--min-sites", "100"))
keepbed <- getarg("--keep-bed")
slug <- paste(merges, collapse = "+")
out <- getarg("--out", file.path(GTR, "figures",
            sprintf("GuitarR_%s_%s_%s_%s.pdf", mode, slug, species, txtype)))
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)

## ---- BED handling: reproduce the original cleaning --------------------------
# The published scripts kept 6 columns, forced strand "+", set score 0 and
# widened every interval by 1 bp on each side.  21b already pads, so only the
# strand/score normalisation is repeated here (Guitar ignores strand anyway).
normalise_bed <- function(path) {
  d <- read.table(path, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                  quote = "")
  d <- d[, 1:6, drop = FALSE]
  colnames(d) <- c("chr", "start", "end", "name", "score", "strand")
  d$score <- 0L
  d$strand <- "+"
  d$name <- "."
  tmp <- tempfile(fileext = ".bed")
  write.table(d, tmp, sep = "\t", quote = FALSE, row.names = FALSE,
              col.names = FALSE)
  tmp
}

concat_bed <- function(paths) {
  tmp <- tempfile(fileext = ".bed")
  con <- file(tmp, "w")
  for (p in paths) writeLines(readLines(p), con)
  close(con)
  tmp
}

bed_in <- function(dir, min_rows = minsit) {
  if (!dir.exists(dir)) return(character(0))
  files <- sort(list.files(dir, pattern = "\\.bed$", full.names = TRUE))
  n <- vapply(files, function(f) length(readLines(f)), integer(1))
  files[n >= min_rows]
}

## ---- assemble the groups Guitar will plot ---------------------------------
conds <- CONDS[[species]]
if (is.null(conds)) stop("unknown species: ", species)
bed_files <- character(0); bed_names <- character(0)
group_kind <- character(0)          # "wt" / "cmp", for colour and linetype

add_group <- function(label, path, kind) {
  bed_files[[label]] <<- normalise_bed(path)
  bed_names <<- c(bed_names, label)
  group_kind <<- c(group_kind, kind)
}
kind_of <- function(cond) if (cond == conds[1]) "wt" else "cmp"

if (mode == "condition") {
  # one curve per library (all its tools pooled), as in the published panel C;
  # passing several merge rules draws the merge sensitivity in the same axes
  for (mg in merges) for (cond in conds) {
    f <- bed_in(file.path(GTR, "bed", paste0(mg, "_pooled"), cond), 1)
    if (!length(f)) { message("no pooled BED: ", cond, "/", mg); next }
    label <- if (length(merges) == 1) cond else paste(cond, mg, sep = " | ")
    add_group(label, f[1], kind_of(cond))
  }
} else if (mode == "pertool") {
  # one curve per tool x library, as in the published panel D
  for (mg in merges) for (cond in conds) {
    for (f in bed_in(file.path(GTR, "bed", mg, cond))) {
      add_group(paste0(tools::file_path_sans_ext(basename(f)), " | ", cond),
                f, kind_of(cond))
    }
  }
} else if (mode == "replicates") {
  # one curve per biological replicate: how much does the library move between
  # replicates, before any consensus is taken
  dirs <- sort(basename(list.dirs(file.path(GTR, "bed"), recursive = FALSE)))
  reps <- dirs[startsWith(dirs, "rep_")]
  for (r in reps) for (cond in conds) {
    files <- bed_in(file.path(GTR, "bed", r, cond))
    if (!length(files)) next
    add_group(paste0(cond, " | ", sub("^rep_", "", r)), concat_bed(files),
              kind_of(cond))
  }
} else stop("unknown --mode: ", mode)

if (!length(bed_files)) stop("no input BED files were assembled")
message(sprintf("%s / %s / %s / %s: %d group(s)", species, mode, slug, txtype,
                length(bed_names)))

## ---- Guitar transcript model (cached; a TxDb cannot be restored from disk,
#      but the guitarTxdb it is turned into is plain data, so cache that) ------
cache <- file.path(GTR, "txdb_cache", sprintf("%s.%s.guitarTxdb.rds", species, txtype))
dir.create(dirname(cache), showWarnings = FALSE, recursive = TRUE)
if (file.exists(cache)) {
  guitarTxdb <- readRDS(cache)
  message("loaded cached guitarTxdb: ", basename(cache))
} else {
  message("building guitarTxdb from ", basename(GTF[species]),
          " (slow for human; cached afterwards)")
  txdb <- makeTxDbFromGFF(file = GTF[species], format = "auto")
  guitarTxdb <- Guitar::makeGuitarTxdb(
    txdb, txfiveutrMinLength = 100, txcdsMinLength = 100,
    txthreeutrMinLength = 100, txlongNcrnaMinLength = 100,
    txlncrnaOverlapmrna = FALSE, txpromoterLength = 1000, txtailLength = 1000,
    txAmblguity = 5, txPrimaryOnly = FALSE, pltTxType = txtype)
  saveRDS(guitarTxdb, cache)
}

## ---- the part of GuitarPlot() that follows the transcript model ------------
# (GuitarPlot's own txGuitarTxdb argument expects a read.table-able file, so the
#  cached object is fed to the same internal steps directly.)
stSampleNum <- 3; stAmblguity <- 5; overlapIndex <- 1; siteLengthIndex <- 1
mapFilterTranscript <- TRUE; adjust <- 1; enableCI <- FALSE

sitesGroup <- Guitar:::.getStGroup(stBedFiles = unname(bed_files),
                                  stGroupName = bed_names)
sitesPointsNormlize <- list(); sitesPointsRelative <- list(); pointWeight <- list()
for (i in seq_along(sitesGroup)) {
  gname <- names(sitesGroup)[[i]]
  message("sampling sites for ", gname)
  sitesPoints <- Guitar::samplePoints(sitesGroup[i], stSampleNum = stSampleNum,
                                      stAmblguity = stAmblguity,
                                      pltTxType = txtype,
                                      stSampleModle = "Equidistance",
                                      mapFilterTranscript = mapFilterTranscript,
                                      guitarTxdb)
  nz <- Guitar::normalize(sitesPoints, guitarTxdb, txtype, overlapIndex,
                          siteLengthIndex)
  sitesPointsRelative[[txtype]][[gname]] <- nz[[1]]
  pointWeight[[txtype]][[gname]] <- nz[[2]]
}
densityCI <- Guitar:::.generateDensity_CI(sitesPointsRelative[[txtype]],
                                         pointWeight[[txtype]],
                                         CI_ResamplingTime = 1000,
                                         adjust = adjust, enableCI = enableCI)
p <- Guitar:::.plotDensity_CI(
  densityCI,
  componentWidth = guitarTxdb[[txtype]]$componentWidthAverage_pct,
  headOrtail = TRUE,
  title = paste("Distribution on", if (txtype == "mrna") "mRNA" else "ncRNA"),
  enableCI = enableCI)

## ---- house style ----------------------------------------------------------
# Guitar's own palette is turquoise/pink; the manuscript used blue for the
# reference library and orange for the perturbed one, so colour by library and,
# when several tools or merge rules share a colour, separate them by linetype.
cols <- ifelse(group_kind == "wt", WT_COLOR, CMP_COLOR)
names(cols) <- bed_names
lins <- if (length(bed_names) > length(conds))
  rep(c("solid", "dashed", "dotted", "longdash", "dotdash", "twodash"),
      length.out = length(bed_names)) else rep("solid", length(bed_names))
names(lins) <- bed_names

# the area fill is for the two-library panels only; with three or more groups
# Guitar's own ribbons hide the lines, so drop that layer entirely
if (length(bed_names) > 2) {
  keep <- !vapply(p$layers, function(l) inherits(l$geom, "GeomRibbon"), logical(1))
  p$layers <- p$layers[keep]
}
fill_alpha <- 0.28
# legend width: 26 keys do not fit under a 7.6 in panel
wide <- max(0, length(bed_names) - 2)
WID <- 7.6 + 0.28 * wide
p <- p + theme_bw(base_family = "Arial", base_size = 13) +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "white", colour = NA),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.8),
        axis.text = element_text(colour = "black"),
        axis.title.x = element_blank(),
        legend.position = "bottom", legend.title = element_blank(),
        legend.key = element_rect(fill = NA, colour = NA),
        legend.text = element_text(size = 10),
        plot.title = element_text(face = "bold", hjust = 0.5)) +
  scale_colour_manual(values = cols, limits = bed_names) +
  scale_fill_manual(values = alpha(cols, fill_alpha), limits = bed_names) +
  scale_linetype_manual(values = lins, limits = bed_names) +
  guides(fill = "none",
         colour = guide_legend(nrow = ceiling(length(bed_names) / 6)))

ggsave(out, plot = p, width = WID, height = 5.8, units = "in")
ggsave(sub("\\.pdf$", ".png", out), plot = p, width = WID, height = 5.8,
       units = "in", dpi = 300)
message("wrote ", out)

if (!is.null(keepbed)) {
  dest <- file.path(keepbed, species, mode, slug, txtype)
  dir.create(dest, showWarnings = FALSE, recursive = TRUE)
  for (nm in names(bed_files))
    file.copy(bed_files[[nm]], file.path(dest, paste0(gsub("[ /]", "_", nm), ".bed")),
              overwrite = TRUE)
  message("kept the assembled BEDs under ", dest)
}
