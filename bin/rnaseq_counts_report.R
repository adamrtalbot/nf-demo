#!/usr/bin/env Rscript

# Script: rnaseq_counts_report.R
# Purpose: Parse RNA-seq count data and generate MultiQC-compatible report

# Load required libraries
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("jsonlite", quietly = TRUE)) install.packages("jsonlite")
library(ggplot2)
library(jsonlite)

# Function to analyze RNA-seq count data
analyze_rnaseq_counts <- function(counts_file, file_prefix = "rnaseq") {
    # Extract directory from prefix if it contains a path
    output_dir <- dirname(file_prefix)
    if (output_dir != ".") {
        # Create output directory if it doesn't exist and isn't just current dir
        if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    }

    # Extract base prefix (filename part only)
    base_prefix <- basename(file_prefix)

    # Replace any non-alphanumeric characters with underscore for safety
    safe_prefix <- gsub("[^a-zA-Z0-9]", "_", base_prefix)

    # Create parent section ID for grouping all sections from this sample
    parent_id <- paste0("sample_", safe_prefix)
    parent_name <- paste0("Sample: ", base_prefix)
    parent_description <- paste0("RNA-seq analysis results for: ", base_prefix)

    # Read count data
    message("Reading count data from: ", counts_file)
    counts <- read.csv(counts_file, row.names = 1, check.names = FALSE)

    # Basic stats
    total_genes <- nrow(counts)
    total_samples <- ncol(counts)
    message("Found ", total_genes, " genes across ", total_samples, " samples")

    # Calculate library sizes
    lib_sizes <- colSums(counts)
    detected_genes <- colSums(counts > 0)
    percent_detected <- 100 * detected_genes / total_genes

    # Create summary dataframe
    sample_summary <- data.frame(
        Sample = colnames(counts),
        TotalReads = lib_sizes,
        DetectedGenes = detected_genes,
        PercentDetected = percent_detected
    )

    # Generate MultiQC-compatible general stats - ensure correct format for merging
    general_stats_file <- paste0(file_prefix, "_general_stats_mqc.tsv")

    # Create header with MultiQC configuration for general stats
    # The key is to use the standard format without customizing section or parent ID
    # MultiQC will automatically merge general stats into the main table
    header <- c(
        "# plot_type: 'generalstats'",
        "# headers:",
        "# - TotalReads:",
        "#     title: 'Total Counts'",
        "#     description: 'Total read counts per sample'",
        "#     format: '{:,.0f}'",
        "#     min: 0",
        "# - DetectedGenes:",
        "#     title: 'Detected Genes'",
        "#     description: 'Number of genes with at least one read'",
        "#     format: '{:,.0f}'",
        "#     min: 0",
        "# - PercentDetected:",
        "#     title: '% Detected'",
        "#     description: 'Percentage of genes detected'",
        "#     suffix: '%'",
        "#     max: 100",
        "#     min: 0",
        "#     scale: 'RdYlGn'"
    )
    # Write the header first
    # writeLines(header, general_stats_file)

    # Then write the data with column names
    # suppressWarnings(
    #     write.table(sample_summary, general_stats_file,
    #         append = TRUE, sep = "\t",
    #         quote = FALSE, row.names = FALSE, col.names = TRUE
    #     )
    # )

    # Create distribution plot
    # Generate MultiQC-compatible boxplot from count distribution
    boxplot_file <- paste0(file_prefix, "_count_distribution_mqc.txt")

    # Write header with configuration - using parent grouping
    header <- c(
        paste0("# id: '", safe_prefix, "_count_distribution'"),
        paste0("# parent_id: '", parent_id, "'"),
        paste0("# parent_name: '", parent_name, "'"),
        paste0("# parent_description: '", parent_description, "'"),
        paste0("# section_name: '", base_prefix, " Count Distribution'"),
        "# description: 'Distribution of log10-transformed read counts per sample'",
        "# plot_type: 'boxplot'",
        "# pconfig:",
        paste0("#     id: '", safe_prefix, "_counts_boxplot'"),
        "#     title: 'Read Count Distribution'",
        "#     xlab: 'Samples'",
        "#     ylab: 'log10(counts+1)'",
        "#     logswitch: false"
    )
    writeLines(header, boxplot_file)

    # Format the data for boxplot - create a complete data frame first
    boxplot_data <- data.frame()

    for (sample in colnames(counts)) {
        sample_data <- log10(counts[, sample] + 1)
        # Format as one number per line with sample name as first column
        sample_df <- data.frame(Sample = sample, log10_counts = sample_data)
        boxplot_data <- rbind(boxplot_data, sample_df)
    }

    # Write the complete data frame at once
    suppressWarnings(
        write.table(boxplot_data, boxplot_file,
            append = TRUE,
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE
        )
    )

    # Perform PCA analysis
    message("Performing PCA analysis")
    # Filter for genes with reasonable expression to improve PCA
    keep_genes <- rowSums(counts >= 10) >= 3 # At least 3 samples with 10+ counts
    if (sum(keep_genes) < 10) {
        # If too few genes pass the filter, just take top expressed genes
        gene_means <- rowMeans(counts)
        keep_genes <- gene_means >= sort(gene_means, decreasing = TRUE)[min(500, length(gene_means))]
    }

    filtered_counts <- counts[keep_genes, ]
    message("Using ", nrow(filtered_counts), " genes for PCA")

    # Log transform data for PCA
    log_counts <- log2(filtered_counts + 1)

    # Transpose for PCA (samples as rows)
    pca_data <- t(log_counts)

    # Calculate PCA
    pca_result <- prcomp(pca_data, scale = TRUE, center = TRUE)

    # Extract variance explained
    var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100

    # Format PCA results for MultiQC
    # Create JSON file for scatter plot (MultiQC compatible)
    pca_file <- paste0(file_prefix, "_pca_plot_mqc.json")

    # Extract first two PCs
    pca_df <- data.frame(
        PC1 = pca_result$x[, 1],
        PC2 = pca_result$x[, 2],
        row.names = rownames(pca_data)
    )

    # Create JSON structure with parent grouping
    pca_json <- list(
        id = paste0(safe_prefix, "_pca"),
        parent_id = parent_id,
        parent_name = parent_name,
        parent_description = parent_description,
        section_name = paste0(base_prefix, " PCA Analysis"),
        description = "Principal Component Analysis of RNA-seq samples based on gene expression profiles.",
        plot_type = "scatter",
        pconfig = list(
            id = paste0(safe_prefix, "_pca_scatter"),
            title = "PCA Plot",
            xlab = sprintf("PC1 (%.1f%%)", var_explained[1]),
            ylab = sprintf("PC2 (%.1f%%)", var_explained[2]),
            dataLabels = list(format = "{point.name}")
        ),
        data = list()
    )

    # Add data points for each sample
    for (sample in rownames(pca_df)) {
        pca_json$data[[sample]] <- list(
            x = pca_df[sample, "PC1"],
            y = pca_df[sample, "PC2"]
        )
    }

    # Write JSON to file
    write(toJSON(pca_json, auto_unbox = TRUE, pretty = TRUE), file = pca_file)

    # Also create a heatmap of sample correlations
    corr_matrix <- cor(log_counts)

    # Write correlation matrix to TSV file for MultiQC heatmap
    corr_file <- paste0(file_prefix, "_sample_correlation_mqc.tsv")

    # Add MultiQC configuration headers as comments with parent grouping
    corr_header <- c(
        paste0("# id: '", safe_prefix, "_correlation'"),
        paste0("# parent_id: '", parent_id, "'"),
        paste0("# parent_name: '", parent_name, "'"),
        paste0("# parent_description: '", parent_description, "'"),
        paste0("# section_name: '", base_prefix, " Sample Correlation'"),
        "# description: 'Correlation matrix showing similarity between samples based on gene expression profiles'",
        "# plot_type: 'heatmap'",
        "# pconfig:",
        paste0("#     id: '", safe_prefix, "_correlation_heatmap'"),
        "#     title: 'Sample Correlation Heatmap'",
        "#     xlab: 'Samples'",
        "#     ylab: 'Samples'",
        "#     min: 0",
        "#     max: 1",
        "#     square: true",
        "#     colstops: [[0, '#d53e4f'], [0.5, '#ffffbf'], [1, '#5e8ca5']]",
        "#     tt_decimals: 2"
    )

    # Write the header first
    writeLines(corr_header, corr_file)

    # Then write the matrix with explicit control over column names
    # Add column names manually to first line to avoid warning
    col_names <- paste(c("", colnames(corr_matrix)), collapse = "\t")
    write(col_names, corr_file, append = TRUE)

    # Then write the data without column names
    suppressWarnings(
        write.table(corr_matrix, corr_file,
            sep = "\t", quote = FALSE,
            append = TRUE, row.names = TRUE, col.names = FALSE
        )
    )

    # Return output files
    message("Generated MultiQC-compatible report files with prefix: ", file_prefix)
    message("  - ", general_stats_file)
    message("  - ", boxplot_file)
    message("  - ", pca_file)
    message("  - ", corr_file)
}

# Main execution
# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if required arguments are provided
if (length(args) < 1) {
    cat("Usage: Rscript rnaseq_counts_report.R counts_file.csv [file_prefix]\n")
    cat("  counts_file.csv: Input counts file in CSV format\n")
    cat("  file_prefix: Prefix for output files [default: counts_file basename without extension]\n")
    cat("\nExample: Rscript rnaseq_counts_report.R sample1_counts.csv results/sample1\n")
    cat("This will create files like results/sample1_general_stats_mqc.txt\n")
    cat("\nOutput files will be MultiQC-compatible and can be included in a MultiQC report.\n")
    quit(status = 1)
}

# Parse arguments
counts_file <- args[1]

# Check if file exists
if (!file.exists(counts_file)) {
    stop("Error: Counts file does not exist: ", counts_file)
}

# If prefix not provided, use the counts filename (without extension)
if (length(args) >= 2) {
    file_prefix <- args[2]
} else {
    # Extract basename and remove extension
    file_prefix <- tools::file_path_sans_ext(basename(counts_file))
}

# Run analysis
analyze_rnaseq_counts(counts_file, file_prefix)
