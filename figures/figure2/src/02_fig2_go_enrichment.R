#!/usr/bin/env Rscript
# 02 -- Figure 2 panel C: species-native GO:BP enrichment (tables only).
#
# Species-native gene sets replace the published Enrichr run, which used a
# human GO library (GO_Biological_Process_2021) for Arabidopsis and mouse as
# well.  Each species is now enriched against its own annotation:
#
#   Arabidopsis  org.At.tair.db  keyType TAIR    (AT1G01010 style locus ids)
#   Mouse        org.Mm.eg.db    keyType ENSEMBL (ENSMUSG... from GRCm39.114)
#   Human        org.Hs.eg.db    keyType ENSEMBL (ENSG... from Ensembl 112)
#
# Foreground / background come from 01_fig2_panel_inputs.py:
#   foreground = genes carrying a high-confidence site (support >= 5 in a unit,
#                majority consensus over the species' independent units)
#   background = genes covering the shared measurable universe of the same
#                units (the testable gene space, not the whole genome)
#
# No graphics device is opened here: the panel is drawn by 03_fig2_figure.py
# (matplotlib), which is the single drawing backend of this revision.
#
# Usage
# -----
# conda run -n enrich_r --no-capture-output Rscript 02_fig2_go_enrichment.R

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(AnnotationDbi)
})

here <- normalizePath(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
tab_dir <- file.path(dirname(here), "tables")

SPECIES <- list(
  Arabidopsis = list(orgdb = "org.At.tair.db", keytype = "TAIR"),
  Mouse       = list(orgdb = "org.Mm.eg.db",   keytype = "ENSEMBL"),
  Human       = list(orgdb = "org.Hs.eg.db",   keytype = "ENSEMBL")
)

for (sp in names(SPECIES)) {
  spec <- SPECIES[[sp]]
  f_in <- file.path(tab_dir, sprintf("fig2c_genes_%s.tsv", sp))
  f_out <- file.path(tab_dir, sprintf("fig2c_gobp_%s.tsv", sp))
  if (!file.exists(f_in)) {
    message(sprintf("[skip] %s: %s not found", sp, f_in))
    next
  }
  if (!requireNamespace(spec$orgdb, quietly = TRUE)) {
    stop(sprintf("OrgDb '%s' is not installed in this environment; install it with ",
                 spec$orgdb),
         "conda install -n enrich_r -c conda-forge -c bioconda bioconductor-",
         gsub("\\.", "-", spec$orgdb), " (no cross-species fallback by design)")
  }
  orgdb <- getExportedValue(spec$orgdb, spec$orgdb)

  d <- read.delim(f_in, stringsAsFactors = FALSE)
  fg <- unique(d$gene_id[d$role == "foreground"])
  bg <- unique(d$gene_id[d$role == "background"])
  fg <- fg[nzchar(fg)]
  bg <- bg[nzchar(bg)]

  # keep only ids the OrgDb can map, so that the universe of the test matches
  # the background actually used by enrichGO
  bg_m <- bg[!is.na(AnnotationDbi::mapIds(orgdb, keys = bg, column = "GO",
                                          keytype = spec$keytype, multiVals = "first"))]
  fg_m <- fg[!is.na(AnnotationDbi::mapIds(orgdb, keys = fg, column = "GO",
                                          keytype = spec$keytype, multiVals = "first"))]

  res <- clusterProfiler::enrichGO(
    gene = fg_m, universe = bg_m, OrgDb = orgdb, keyType = spec$keytype,
    ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
    minGSSize = 10, maxGSSize = 500, readable = FALSE)

  out <- as.data.frame(res)
  out$species <- sp
  out$orgdb <- spec$orgdb
  out$keytype <- spec$keytype
  out$n_foreground_input <- length(fg)
  out$n_foreground_mapped <- length(fg_m)
  out$n_background_input <- length(bg)
  out$n_background_mapped <- length(bg_m)
  if (nrow(out)) {
    out <- out[order(out$p.adjust, out$pvalue), ]
    out$rank <- seq_len(nrow(out))
  }
  write.table(out, f_out, sep = "\t", quote = FALSE, row.names = FALSE)
  message(sprintf("%s: %d fg (%d mapped) / %d bg (%d mapped) -> %d BP terms -> %s",
                  sp, length(fg), length(fg_m), length(bg), length(bg_m),
                  nrow(out), basename(f_out)))
}
