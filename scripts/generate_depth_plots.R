#!/usr/bin/env Rscript
#' Generate depth coverage plots for nanopore sequencing data.
#'
#' 本脚本读取 `samtools depth` 的输出（3 列：Chr / Position / Depth），
#' 为每个样本生成可视化覆盖度图，用于评估 genome / transcriptome
#' 的整体覆盖分布。
#'
#' 依赖包：
#'   - ggplot2
#'   - data.table
#'   - argparse（命令行参数）
#'
#' 用法：
#'   Rscript scripts/generate_depth_plots.R \
#'       --depth_file results/qc/sample1_genome.depth \
#'       --sample_name sample1 \
#'       --output_file results/qc/sample1_depth.png
#'
#' 参数：
#'   --depth_file   `samtools depth` 输出的 3 列文件
#'   --sample_name  展示在 plot 标题中的样本名
#'   --output_file  PNG 输出路径

# ------------- 解析命令行参数 -------------
suppressMessages(library("argparse", quietly = TRUE))
parser <- ArgumentParser(description = "Generate coverage depth plots")
parser$add_argument("--depth_file",  required = TRUE, help = "samtools depth 输出文件")
parser$add_argument("--sample_name", required = TRUE, help = "样本名，展示在标题中")
parser$add_argument("--output_file", required = TRUE, help = "PNG 输出路径")
args   <- parser$parse_args()

# ------------- 加载包 -------------
suppressMessages({
  library(ggplot2)
  library(data.table)
})

# ------------- 记录运行环境（附带 sessionInfo） -------------
log_file <- sub("\\.png$", ".log", args$output_file)
sink(log_file, split = TRUE)
cat("RNAModBench depth plot\n")
cat("depth_file :", args$depth_file,  "\n")
cat("sample     :", args$sample_name, "\n")
cat("output_file:", args$output_file, "\n")
cat("Timestamp  :", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("===== sessionInfo =====\n")
print(sessionInfo())
cat("\n")
sink()

# Function to read and process depth data
read_depth_data <- function(depth_file) {
  # Read depth data (assuming format: chr, pos, depth)
  depth_data <- fread(depth_file, header=FALSE, sep="\t")
  colnames(depth_data) <- c("Chr", "Position", "Depth")
  
  # Add facet information for plotting
  depth_data$Facet <- depth_data$Chr
  
  return(depth_data)
}

# Function to create depth plot
create_depth_plot <- function(depth_data, output_file, sample_name) {
  # Convert column names for plotting
  colnames(depth_data) <- c("Facet", "X", "Y")
  
  # Create the plot
  p <- ggplot(depth_data, aes(x = X, y = Y, color = Facet)) +
    geom_line() +
    geom_point() +  # Optional: add data points
    facet_wrap(~ Facet, ncol = 1) +  # Facet by chromosome
    labs(title = paste("Coverage Depth -", sample_name), 
         x = "Position", 
         y = "Depth") +
    theme_bw() +  # Black and white theme
    theme(legend.position = "none")  # Remove legend for cleaner plot
  
  # Save the plot
  ggsave(output_file, plot = p, width = 12, height = 8, dpi = 300)
  
  cat("Depth plot saved to:", output_file, "\n")
}

# Main function
main <- function(args) {
  # Parse arguments
  depth_file <- args[1]
  output_file <- args[2]
  sample_name <- ifelse(length(args) > 2, args[3], "Sample")
  
  # Read and process data
  cat("Reading depth data from:", depth_file, "\n")
  depth_data <- read_depth_data(depth_file)
  
  # Create plot
  cat("Creating depth plot...\n")
  create_depth_plot(depth_data, output_file, sample_name)
  
  cat("Done!\n")
}

# Run the script
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 2) {
    cat("Usage: Rscript generate_depth_plots.R <depth_file> <output_file> [sample_name]\n")
    quit(status = 1)
  }
  main(args)
}