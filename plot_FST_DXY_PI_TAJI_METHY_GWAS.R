#!/usr/bin/env Rscript

# Genomic Data Visualization Script
# This script creates multi-panel plots for genomic analyses including:
# FST, Tajima's D, xpEHH, nucleotide diversity (pi), DXY, GWAS, and optional methylation data

# Load necessary libraries
library(dplyr)
library(rtracklayer)

#=============================================================================
# USER CONFIGURATION SECTION - MODIFY THESE PATHS AND PARAMETERS
#=============================================================================

# Required file paths (set these to your actual file locations)
GWAS_file <- "path/to/your/gwas_results.assoc.txt"  # GWAS association results
DXY_file <- "path/to/your/fst_dxy_data.csv"        # FST and DXY data (CSV format)
xpEHH_file <- "path/to/your/xpehh_results.tsv"     # xpEHH results
Tajima_file <- "path/to/your/tajima_d_results.Tajima.D"  # Tajima's D results
PI_pop1_file <- "path/to/your/population1_pi.windowed.pi"  # Nucleotide diversity for population 1
PI_pop2_file <- "path/to/your/population2_pi.windowed.pi"  # Nucleotide diversity for population 2
GTF_file <- "path/to/your/genome_annotation.gtf"   # Gene annotation file

# Optional methylation file (set to NULL if not available)
Methylation_file <- "path/to/your/methylation_data.txt"  # Set to NULL if no methylation data
# Methylation_file <- NULL  # Uncomment this line if you don't have methylation data

# Output file prefix
output_prefix <- "genomic_analysis"

# Plot region parameters (modify these for your region of interest)
focal_chromosome <- "chr1"        # Target chromosome (e.g., "chr1", "chr23")
window_start <- 1000000          # Start position of the region to plot
window_end <- 2000000            # End position of the region to plot
highlight_start <- 1200000       # Start of highlight region (optional)
highlight_end <- 1800000         # End of highlight region (optional)

# Population labels for the legend
population1_name <- "Population1"  # Name for first population (e.g., "Littoral", "Control")
population2_name <- "Population2"  # Name for second population (e.g., "Benthic", "Treatment")

#=============================================================================
# SCRIPT EXECUTION - NO MODIFICATIONS NEEDED BELOW THIS LINE
#=============================================================================

# Check if methylation file is provided and exists
have_methylation <- !is.null(Methylation_file) && file.exists(Methylation_file)

# Color scheme for plots
if (have_methylation) {
  colors <- c("royalblue", "royalblue", "royalblue", "goldenrod2", "steelblue2", "black", "blue", "purple")
  legend_labels <- c("FST", "TajimaD", "xpEHH", paste("PI", population1_name), paste("PI", population2_name), "DXY", "GWAS", "Methylation")
} else {
  colors <- c("royalblue", "royalblue", "royalblue", "goldenrod2", "steelblue2", "black", "blue")
  legend_labels <- c("FST", "TajimaD", "xpEHH", paste("PI", population1_name), paste("PI", population2_name), "DXY", "GWAS")
}

# Function to safely read files with error handling
safe_read <- function(file_path, read_function, ...) {
  if (!file.exists(file_path)) {
    stop(paste("File not found:", file_path))
  }
  tryCatch({
    read_function(file_path, ...)
  }, error = function(e) {
    stop(paste("Error reading file", file_path, ":", e$message))
  })
}

# Load and process data files
cat("Loading data files...\n")

# Load FST/DXY data
dxy_fst <- safe_read(DXY_file, read.table, header = TRUE, sep = ",") %>%
  filter(scaffold == focal_chromosome, mid >= window_start, mid <= window_end)

x1 <- dxy_fst %>% transmute(CHR = scaffold, POSITION = mid, value = Fst_benthic_littoral)
x6 <- dxy_fst %>% transmute(CHR = scaffold, POSITION = mid, value = dxy_benthic_littoral)

# Load Tajima's D data
x2 <- safe_read(Tajima_file, read.table, header = TRUE) %>%
  transmute(CHR = gsub("^chr", "", CHROM), POSITION = BIN_START + 5000, value = TajimaD) %>%
  filter(CHR == gsub("^chr", "", focal_chromosome), POSITION >= window_start, POSITION <= window_end)

# Load xpEHH data
x3 <- safe_read(xpEHH_file, read.table, header = TRUE) %>%
  transmute(CHR = gsub("^chr", "", CHR), POSITION = POSITION, value = LOGPVALUE) %>%
  filter(CHR == gsub("^chr", "", focal_chromosome), POSITION >= window_start, POSITION <= window_end)

# Load nucleotide diversity (pi) for both populations
x4_pop1 <- safe_read(PI_pop1_file, read.table, header = TRUE) %>%
  transmute(CHR = gsub("^chr", "", CHROM), POSITION = (BIN_START + BIN_END) / 2, value = PI) %>%
  filter(CHR == gsub("^chr", "", focal_chromosome), POSITION >= window_start, POSITION <= window_end)

x4_pop2 <- safe_read(PI_pop2_file, read.table, header = TRUE) %>%
  transmute(CHR = gsub("^chr", "", CHROM), POSITION = (BIN_START + BIN_END) / 2, value = PI) %>%
  filter(CHR == gsub("^chr", "", focal_chromosome), POSITION >= window_start, POSITION <= window_end)

# Load GWAS results
GWAS_data <- safe_read(GWAS_file, read.table, header = TRUE) %>%
  mutate(chr = ifelse(grepl("^chr", chr), chr, paste0("chr", chr))) %>%
  filter(chr == focal_chromosome, ps >= window_start, ps <= window_end)
x7 <- GWAS_data %>% transmute(CHR = chr, POSITION = ps, value = -log10(p_score))

# Load methylation data (if available)
if (have_methylation) {
  cat("Loading methylation data...\n")
  Meth_data <- safe_read(Methylation_file, read.table, header = TRUE) %>%
    filter(chr == focal_chromosome, start >= window_start, end_benthic <= window_end)
  x8 <- Meth_data %>% transmute(CHR = chr, POSITION = (start + end_benthic) / 2, value = Difference)
}

# Load gene annotations
cat("Loading gene annotations...\n")
gtf_data <- rtracklayer::import(GTF_file)
gtf_df <- as.data.frame(gtf_data)
gtf_subset <- gtf_df %>%
  filter(seqnames == gsub("chr", "", focal_chromosome), start >= window_start, end <= window_end)

# Process gene information for plotting
genes <- gtf_subset %>%
  filter(type == "gene") %>%
  select(gene_id, start, end) %>%
  arrange(start)

# Assign plotting levels to avoid gene overlap
margin <- 1000
genes$level <- 1
genes$t_start <- NA
genes$t_end <- NA
gtf_subset$level <- NA

for (i in seq_len(nrow(genes))) {
  transcripts <- gtf_subset %>% 
    filter(gene_id == genes$gene_id[i] & type == "transcript")
  
  if (nrow(transcripts) > 0) {
    genes$t_start[i] <- min(transcripts$start, na.rm = TRUE)
    genes$t_end[i] <- max(transcripts$end, na.rm = TRUE)
  } else {
    genes$t_start[i] <- genes$start[i]
    genes$t_end[i] <- genes$end[i]
  }
  
  if (i > 1 && genes$t_start[i] <= genes$t_end[i - 1] + margin) {
    genes$level[i] <- genes$level[i - 1] + 1
    if (genes$level[i] > 3) genes$level[i] <- 1
  }
  
  gtf_subset$level[gtf_subset$gene_id == genes$gene_id[i]] <- genes$level[i]
}

max_levels <- max(genes$level, na.rm = TRUE)

# Plotting functions
add_highlight <- function() {
  rect(highlight_start, par("usr")[3], highlight_end, par("usr")[4], 
       col = rgb(0.7, 0.7, 0.7, 0.3), border = NA)
}

plot_genomic_data <- function() {
  # Plot FST
  plot(x1$POSITION, x1$value, type = "l", col = colors[1], lwd = 2,
       xlab = "", ylab = "FST", xlim = c(window_start, window_end), 
       ylim = range(x1$value, na.rm = TRUE), xaxt = "n", las = 1)
  add_highlight()
  lines(x1$POSITION, x1$value, col = colors[1], lwd = 2)
  
  # Plot Tajima's D
  plot(x2$POSITION, x2$value, type = "l", col = colors[2], lwd = 2,
       xlab = "", ylab = "Tajima's D", xlim = c(window_start, window_end),
       ylim = range(x2$value, na.rm = TRUE), xaxt = "n", las = 1)
  add_highlight()
  lines(x2$POSITION, x2$value, col = colors[2], lwd = 2)
  
  # Plot xpEHH
  plot(x3$POSITION, x3$value, type = "l", col = colors[3], lwd = 2,
       xlab = "", ylab = "xpEHH", xlim = c(window_start, window_end),
       ylim = range(x3$value, na.rm = TRUE), xaxt = "n", las = 1)
  add_highlight()
  lines(x3$POSITION, x3$value, col = colors[3], lwd = 2)
  
  # Plot nucleotide diversity (pi) for both populations
  ylims <- range(c(x4_pop1$value, x4_pop2$value), na.rm = TRUE)
  plot(x4_pop1$POSITION, x4_pop1$value, type = "l", col = colors[4], lwd = 2,
       xlab = "", ylab = "Nucleotide Diversity (π)", xlim = c(window_start, window_end),
       ylim = ylims, xaxt = "n", las = 1)
  add_highlight()
  lines(x4_pop1$POSITION, x4_pop1$value, col = colors[4], lwd = 2)
  lines(x4_pop2$POSITION, x4_pop2$value, col = colors[5], lwd = 2)
  legend("topright", legend = c(population1_name, population2_name), 
         col = colors[4:5], lty = 1, lwd = 2, cex = 0.8)
  
  # Plot DXY
  plot(x6$POSITION, x6$value, type = "l", col = colors[6], lwd = 2,
       xlab = "", ylab = "DXY", xlim = c(window_start, window_end),
       ylim = range(x6$value, na.rm = TRUE), xaxt = "n", las = 1)
  add_highlight()
  lines(x6$POSITION, x6$value, col = colors[6], lwd = 2)
  
  # Plot GWAS results
  plot(x7$POSITION, x7$value, type = "n", xlab = "", ylab = "-log10(p)",
       xlim = c(window_start, window_end), ylim = range(x7$value, na.rm = TRUE),
       xaxt = "n", las = 1)
  add_highlight()
  points(x7$POSITION, x7$value, col = colors[7], pch = 16, cex = 0.8)
  
  # Plot methylation data (if available)
  if (have_methylation) {
    plot(x8$POSITION, x8$value, type = "l", col = colors[8], lwd = 2,
         xlab = "", ylab = "Methylation Diff", xlim = c(window_start, window_end),
         ylim = range(x8$value, na.rm = TRUE), xaxt = "n", las = 1)
    add_highlight()
    lines(x8$POSITION, x8$value, col = colors[8], lwd = 2)
  }
}

plot_gene <- function(gene_id, gtf_data, levels) {
  gene_data <- gtf_data[gtf_data$gene_id == gene_id,]
  if (nrow(gene_data) == 0) return()
  
  gene_level <- gene_data[1, "level"]
  level_width <- 1 / levels
  middle <- level_width * (levels - gene_level) + 0.5 * level_width
  top <- middle + 0.3 * level_width
  bottom <- middle - 0.3 * level_width
  
  # Draw gene structure
  exons <- gene_data[gene_data$type == "exon",]
  segments(min(gene_data$start), middle, max(gene_data$end), middle, col = "red", lwd = 3)
  
  if (nrow(exons) > 0) {
    apply(exons, 1, function(exon) {
      rect(exon["start"], bottom, exon["end"], top, col = "red", border = NA)
    })
  }
  
  # Add gene label
  gene_name <- ifelse(!is.na(gene_data$gene_name[1]), gene_data$gene_name[1], gene_id)
  if (gene_data$strand[1] == "+") {
    label <- paste0(gene_name, " >")
  } else {
    label <- paste0("< ", gene_name)
  }
  text(x = mean(c(min(gene_data$start), max(gene_data$end))), y = top + 0.02, 
       labels = label, col = "red", cex = 0.8)
}

plot_gene_annotations <- function(gene_list, gtf_data, levels) {
  plot(0, type = 'n', axes = FALSE, xlab = "", ylab = "", frame.plot = TRUE,
       xlim = c(window_start, window_end), ylim = c(0, 1.15))
  rect(highlight_start, 0, highlight_end, 1.15, 
       col = rgb(0.7, 0.7, 0.7, 0.3), border = NA)
  axis(1, xpd = TRUE)
  
  for (gene_id in gene_list$gene_id) {
    plot_gene(gene_id, gtf_data, levels)
  }
}

# Create the final plot
cat("Creating visualization...\n")

# Calculate plot dimensions
annotation_height <- 0.25 + 0.75 * max_levels
num_panels <- ifelse(have_methylation, 8, 7)  # 7 data panels + 1 annotation panel (or 8 if methylation)
png_height <- 210 * (num_panels + annotation_height)

# Set up the plot layout
png(filename = paste0(output_prefix, "_genomic_visualization.png"), 
    width = 3000, height = png_height)

if (have_methylation) {
  layout(matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 0, 9), ncol = 1),
         widths = rep(17, 9),
         heights = c(rep(2, 8), 0.5, annotation_height))
} else {
  layout(matrix(c(1, 2, 3, 4, 5, 6, 7, 0, 8), ncol = 1),
         widths = rep(17, 8),
         heights = c(rep(2, 7), 0.5, annotation_height))
}

par(mar = c(0, 5, 2, 1), oma = c(5, 4, 2, 1), cex = 1.8)

# Generate all plots
plot_genomic_data()
plot_gene_annotations(genes, gtf_subset, max_levels)

# Add main x-axis label
mtext(text = paste0("Genomic Position on ", focal_chromosome),
      side = 1, line = 1.5, outer = TRUE, cex = 2.2)

dev.off()

cat("Plot saved as:", paste0(output_prefix, "_genomic_visualization.png\n"))
cat("Analysis complete!\n")

if (have_methylation) {
  cat("Methylation data was included in the analysis.\n")
} else {
  cat("No methylation data was provided - analysis completed without methylation panel.\n")
}
