#!/bin/bash

# Generate counts and metadata files
# mkdir -p data results
# Rscript utils/generate_counts.R -g 10000 -s 6 -c A,B,C -r 123 -p Exp1
# Rscript utils/generate_counts.R -g 10000 -s 6 -c A,B,C -r 321 -p Exp2
# Rscript utils/generate_counts.R -g 10000 -s 6 -c A,B,C -r 213 -p Exp3

# Run the report
Rscript rnaseq_counts_report.R data/Exp1_counts.csv results/Exp1
Rscript rnaseq_counts_report.R data/Exp2_counts.csv results/Exp2
Rscript rnaseq_counts_report.R data/Exp3_counts.csv results/Exp3

# Generate the MultiQC report
multiqc -f results/