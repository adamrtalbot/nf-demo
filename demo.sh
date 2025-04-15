#!/bin/bash

rm -f results/*
rm -f multiqc_report.html
rm -f multiqc_data/*

# Create necessary directories
mkdir -p data results

# Generate counts and metadata files with different seeds
Rscript data/generate_counts.R -g 10000 -s 6 -c A,B,C -r 123 -p data/Exp1
Rscript data/generate_counts.R -g 10000 -s 6 -c A,B,C -r 321 -p data/Exp2
Rscript data/generate_counts.R -g 10000 -s 6 -c A,B,C -r 213 -p data/Exp3

# Run the report

Rscript rnaseq_counts_report.R data/Exp1_counts.csv results/Exp1
Rscript rnaseq_counts_report.R data/Exp2_counts.csv results/Exp2
Rscript rnaseq_counts_report.R data/Exp3_counts.csv results/Exp3

# Generate the MultiQC report
multiqc -f results/