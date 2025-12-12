#!/usr/bin/env Rscript
"""
Generate depth coverage plots for nanopore sequencing data.
This script creates coverage depth visualization plots.
"""

library(ggplot2)
library(data.table)

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