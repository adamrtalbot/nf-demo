#!/usr/bin/env Rscript

# Script: generate_counts.R
# Purpose: Generate fake RNA-seq count data compatible with rnaseq_counts_report.R

# Function to create a fake RNA-seq count dataset
generate_rnaseq_counts <- function(num_genes = 20000,
                                   num_samples = 12,
                                   groups = c("Control", "Treatment"),
                                   prefix = "experiment",
                                   random_seed = 42) {
    # Set seed for reproducibility
    set.seed(random_seed)

    # Define output files based on prefix
    output_file <- paste0(prefix, "_counts.csv")
    metadata_file <- paste0(prefix, "_metadata.csv")

    # Create more descriptive sample names that include the group and prefix
    samples_per_group <- num_samples / length(groups)
    sample_names <- c()

    for (group in groups) {
        group_samples <- paste0(prefix, "_", group, "_Sample_", 1:samples_per_group)
        sample_names <- c(sample_names, group_samples)
    }

    # Create gene names (Ensembl-like IDs)
    gene_names <- paste0("ENSG", formatC(1:num_genes, width = 11, format = "d", flag = "0"))

    # Create expression profiles with realistic properties
    # Log-normal distribution for mean expression levels
    mean_expression <- exp(rnorm(num_genes, mean = 3, sd = 2))

    # Create count matrix
    counts <- matrix(0, nrow = num_genes, ncol = num_samples)

    # Assign samples to groups
    sample_groups <- rep(groups, each = num_samples / length(groups))

    # Generate group-specific effects (treatment effect)
    is_de <- rbinom(num_genes, 1, prob = 0.1) == 1 # 10% of genes are differentially expressed
    effect_size <- rnorm(sum(is_de), mean = 0, sd = 2) # effect sizes for DE genes

    # Add biological variation between samples
    for (i in 1:num_samples) {
        # Add some sample-specific variation (library size differences)
        sample_factor <- runif(1, min = 0.7, max = 1.3)

        # Add group-specific effects
        group_factor <- 1
        if (sample_groups[i] == "Treatment" && any(is_de)) {
            # Apply effect sizes to the differentially expressed genes
            group_effects <- rep(1, num_genes)
            group_effects[is_de] <- 2^effect_size # log2 fold change
            mean_expression_adjusted <- mean_expression * group_effects
        } else {
            mean_expression_adjusted <- mean_expression
        }

        # Generate counts from negative binomial (typical for RNA-seq)
        counts[, i] <- rnbinom(num_genes,
            mu = mean_expression_adjusted * sample_factor,
            size = 5 # dispersion parameter (smaller = more overdispersion)
        )
    }

    # Add some genes with very low or zero counts
    low_expr_genes <- sample(1:num_genes, size = num_genes * 0.3)
    counts[low_expr_genes, ] <- rpois(length(low_expr_genes) * num_samples, lambda = 0.5)

    # Add some genes with no expression
    zero_expr_genes <- sample(1:num_genes, size = num_genes * 0.05)
    counts[zero_expr_genes, ] <- 0

    # Create dataframe
    count_df <- as.data.frame(counts)
    colnames(count_df) <- sample_names
    rownames(count_df) <- gene_names

    # Ensure data directory exists
    dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)

    # Write to CSV
    write.csv(count_df, file = output_file)

    # Create metadata file with more information
    metadata <- data.frame(
        sample_id = sample_names,
        condition = sample_groups,
        batch = sample(1:3, num_samples, replace = TRUE) # random batch assignment
    )

    write.csv(metadata, file = metadata_file, row.names = FALSE)

    message("Created fake RNA-seq count file: ", output_file)
    message("Created sample metadata file: ", metadata_file)
}

# Display usage information
show_usage <- function() {
    cat("
Usage: Rscript generate_counts.R [options]

Generate fake RNA-seq count data compatible with rnaseq_counts_report.R

Options:
  -g, --genes NUM        Number of genes to generate [default: 20000]
  -s, --samples NUM      Number of samples to generate [default: 12]
  -c, --conditions STR   Comma-separated list of condition names [default: Control,Treatment]
  -r, --randomseed NUM   Random seed for reproducibility [default: 42]
  -p, --prefix STR       Prefix for sample names and output files [default: 'experiment']
  -h, --help             Show this help message and exit

Examples:
  Rscript generate_counts.R -p Experiment1
  Rscript generate_counts.R -g 15000 -s 6 -c Control,Treatment,Mutant -p Exp2

Report bugs to: your-email@example.com
")
}

# Parse command line arguments
if (!interactive()) {
    args <- commandArgs(trailingOnly = TRUE)

    # Set defaults
    num_genes <- 20000
    num_samples <- 12
    conditions <- "Control,Treatment"
    random_seed <- 42
    prefix <- "experiment"

    # Show help if no arguments or --help
    if (length(args) == 0 || "--help" %in% args || "-h" %in% args) {
        show_usage()
        quit(save = "no", status = 0)
    }

    # Parse arguments
    i <- 1
    while (i <= length(args)) {
        if (args[i] == "-g" || args[i] == "--genes") {
            if (i + 1 <= length(args)) {
                num_genes <- as.numeric(args[i + 1])
                i <- i + 2
            } else {
                stop("Missing value for ", args[i])
            }
        } else if (args[i] == "-s" || args[i] == "--samples") {
            if (i + 1 <= length(args)) {
                num_samples <- as.numeric(args[i + 1])
                i <- i + 2
            } else {
                stop("Missing value for ", args[i])
            }
        } else if (args[i] == "-c" || args[i] == "--conditions") {
            if (i + 1 <= length(args)) {
                conditions <- args[i + 1]
                i <- i + 2
            } else {
                stop("Missing value for ", args[i])
            }
        } else if (args[i] == "-r" || args[i] == "--randomseed") {
            if (i + 1 <= length(args)) {
                random_seed <- as.numeric(args[i + 1])
                i <- i + 2
            } else {
                stop("Missing value for ", args[i])
            }
        } else if (args[i] == "-p" || args[i] == "--prefix") {
            if (i + 1 <= length(args)) {
                prefix <- args[i + 1]
                i <- i + 2
            } else {
                stop("Missing value for ", args[i])
            }
        } else {
            stop("Unknown argument: ", args[i])
        }
    }

    # Process conditions
    groups <- strsplit(conditions, ",")[[1]]

    # Validate inputs
    if (num_genes <= 0) stop("Number of genes must be positive")
    if (num_samples <= 0) stop("Number of samples must be positive")
    if (length(groups) == 0) stop("At least one condition must be specified")
    if (num_samples %% length(groups) != 0) {
        warning(
            "Number of samples (", num_samples, ") is not divisible by number of conditions (",
            length(groups), "). Sample groups will be unbalanced."
        )
    }

    # Generate data - pass the random seed to the function
    message("Generating ", num_genes, " genes across ", num_samples, " samples in ", length(groups), " conditions with prefix '", prefix, "'")
    generate_rnaseq_counts(num_genes, num_samples, groups, prefix, random_seed)
} else {
    # If running in interactive mode
    message("Run this script directly with arguments or call the generate_rnaseq_counts() function")
}

# Interactive usage example:
# source("generate_counts.R")
# generate_rnaseq_counts(15000, 8, c("Control", "Treatment", "Mutant"), "Experiment1", 123)
