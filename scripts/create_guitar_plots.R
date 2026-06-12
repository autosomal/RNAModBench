#!/usr/bin/env Rscript
#' Create Guitar plots for RNA modification visualization.
#'
#' 本脚本读取 RNAModBench 中 `results/<Tool>/<sample>/*_processed.txt`
#' 的标准 7 列 TSV 文件（Chr / Start / End / Status / Prob / Strand / mod_ratio），
#' 将修饰位点投影到 gene body（5'UTR / CDS / 3'UTR）并生成 Guitar 图。
#'
#' 依赖包（请提前安装）：
#'   - GenomicRanges
#'   - rtracklayer
#'   - Guitar（来自 Bioconductor）
#'   - ggplotify
#'   - ggplot2
#'
#' 用法：
#'   Rscript scripts/create_guitar_plots.R \
#'       --input_dir results/summary \
#'       --gtf reference/genes.gtf \
#'       --output_dir results/summary/guitar
#'
#' 参数：
#'   --input_dir    包含 *_processed.txt 的目录（递归搜索）
#'   --gtf          Ensembl/Gencode 风格的 GTF（需包含 transcript_id）
#'   --output_dir   输出目录（自动创建）
#'
#' 输出：
#'   - guitar_plots_mrna.png       : mRNA 上的修饰分布
#'   - guitar_plots_all_transcripts.png : 所有转录本的修饰分布
#'   - guitar_plots.log            : sessionInfo 与运行信息

# ------------- 解析命令行参数 -------------
suppressMessages(library("argparse", quietly = TRUE))
parser <- ArgumentParser(description = "Generate Guitar plots from RNAModBench output")
parser$add_argument("--input_dir",  required = TRUE, help = "目录，包含 *_processed.txt 文件")
parser$add_argument("--gtf",        required = TRUE, help = "参考 GTF 文件")
parser$add_argument("--output_dir", required = TRUE, help = "输出目录")
args   <- parser$parse_args()

# ------------- 加载包 -------------
suppressMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(Guitar)
  library(ggplotify)
  library(ggplot2)
  library(data.table)
})

# ------------- 创建输出目录 -------------
dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------- 记录运行环境 -------------
log_file <- file.path(args$output_dir, "guitar_plots.log")
sink(log_file, split = TRUE)
cat("RNAModBench Guitar plots\n")
cat("input_dir :", args$input_dir,  "\n")
cat("gtf       :", args$gtf,        "\n")
cat("output_dir:", args$output_dir, "\n")
cat("Timestamp :", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("===== sessionInfo =====\n")
print(sessionInfo())
cat("\n")
sink()

# Function to process BED files
clean_bed_files <- function(input_dir, output_dir) {
  # Create output directory
  dir.create(output_dir, showWarnings = FALSE)
  
  # Find all bed/txt files
  files <- list.files(input_dir, pattern = "\\.(bed|txt)$", full.names = TRUE)
  
  # Process each file
  for (f in files) {
    # Read file (allowing more than 6 columns)
    dat <- read.table(f, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    
    # Keep only first 6 columns
    dat <- dat[, 1:6, drop = FALSE]
    
    # Set column names
    colnames(dat) <- c("chr", "start", "end", "name", "score", "strand")
    
    # Set strand to "+" for all entries
    dat$strand <- "+"
    
    # Write cleaned file
    out_file <- file.path(output_dir, basename(f))
    write.table(dat, file = out_file, sep = "\t", quote = FALSE, 
                row.names = FALSE, col.names = FALSE)
  }
  
  return(output_dir)
}

# Function to extend BED regions
extend_bed_regions <- function(bed_dir) {
  files <- list.files(bed_dir, pattern = "\\.bed$", full.names = TRUE)
  
  for (f in files) {
    # Import as GRanges
    gr <- import(f, format = "BED")
    
    # Extend regions by 1 bp on each side
    start(gr) <- start(gr) - 1L
    end(gr) <- end(gr) + 1L
    
    # Write back
    export(gr, f, format = "BED")
  }
  
  return(bed_dir)
}

# Function to create Guitar plots
create_guitar_plots <- function(txdb, bed_dir, sample_name) {
  # Get BED files
  bed_files <- list.files(bed_dir, pattern = "\\.bed$", full.names = TRUE)
  
  if (length(bed_files) == 0) {
    stop("No BED files found in ", bed_dir)
  }
  
  # Create output directory
  output_dir <- file.path("results", "guitar_plots", sample_name)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Create plots for mRNA
  cat("Creating Guitar plots for mRNA...\n")
  guitar_plot_mrna <- GuitarPlot(
    txTxdb = txdb,
    stBedFiles = bed_files,
    headOrtail = TRUE,
    enableCI = FALSE,
    pltTxType = c("mrna"),
    stSampleNum = 3,
    stGroupName = tools::file_path_sans_ext(basename(bed_files))
  )
  
  # Save mRNA plot
  ggsave(file.path(output_dir, "guitar_plot_mrna.png"), 
         guitar_plot_mrna, width = 12, height = 8, dpi = 300)
  
  # Create plots for ncRNA
  cat("Creating Guitar plots for ncRNA...\n")
  guitar_plot_ncrna <- GuitarPlot(
    txTxdb = txdb,
    stBedFiles = bed_files,
    headOrtail = TRUE,
    enableCI = FALSE,
    pltTxType = c("ncrna"),
    stSampleNum = 3,
    stGroupName = tools::file_path_sans_ext(basename(bed_files))
  )
  
  # Save ncRNA plot
  ggsave(file.path(output_dir, "guitar_plot_ncrna.png"), 
         guitar_plot_ncrna, width = 12, height = 8, dpi = 300)
  
  cat("Guitar plots saved to:", output_dir, "\n")
}

# Main function
main <- function(args) {
  # Parse arguments
  gtf_file <- args[1]
  input_dir <- args[2]
  sample_name <- args[3]
  
  # Create TxDb from GTF
  cat("Creating TxDb from GTF:", gtf_file, "\n")
  txdb <- makeTxDbFromGFF(file = gtf_file, format = "auto")
  
  # Process BED files
  cat("Processing BED files in:", input_dir, "\n")
  clean_dir <- file.path("temp", "clean_bed", sample_name)
  cleaned_bed_dir <- clean_bed_files(input_dir, clean_dir)
  
  # Extend BED regions
  extended_bed_dir <- extend_bed_regions(cleaned_bed_dir)
  
  # Create Guitar plots
  create_guitar_plots(txdb, extended_bed_dir, sample_name)
  
  # Clean up temp files
  unlink("temp", recursive = TRUE)
  
  cat("Done!\n")
}

# Run the script
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 3) {
    cat("Usage: Rscript create_guitar_plots.R <gtf_file> <input_dir> <sample_name>\n")
    quit(status = 1)
  }
  main(args)
}