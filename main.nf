#!/usr/bin/env nextflow

// Define parameters with defaults
params.counts_dir = "data"
params.output_dir = "results"
params.pattern = "*_counts.csv"

// Define the process to analyze RNA-seq counts and generate MultiQC-compatible reports
process analyzeRnaSeqCounts {
    tag "${sample_id}"

    conda "conda-forge::r-ggplot2 conda-forge::r-jsonlite"
    
    input:
    path counts_file
    
    output:
    path "${sample_id}_count_distribution_mqc.txt", emit: distribution
    path "${sample_id}_pca_plot_mqc.json"         , emit: mqc
    path "${sample_id}_sample_correlation_mqc.tsv", emit: correlation
    
    script:
    sample_id = counts_file.simpleName
    """
    rnaseq_counts_report.R ${counts_file} ${sample_id}
    """
}

// Execute MultiQC to aggregate all reports
process runMultiQC {
    publishDir "${params.output_dir}", mode: 'copy'

    conda "bioconda::multiqc"
    
    input:
    path '*'
    
    output:
    path "multiqc_report.html"
    path "multiqc_data"
    
    script:
    """
    multiqc -f .
    """
}

// Define the workflow
workflow {
    log.info """
    ========================
    RNA-SEQ QC PIPELINE
    ========================
    
    Input parameters:
      counts_dir: ${params.counts_dir}
      output_dir: ${params.output_dir}
      pattern:    ${params.pattern}
    """

    // Define the input channel for count files
    counts_files = Channel
        .fromPath("${params.counts_dir}/${params.pattern}")

    // Run the RNA-seq analysis
    analyzeResults = analyzeRnaSeqCounts(counts_files)
    
    multiqcFiles = Channel.empty()
                    .mix(analyzeResults.distribution)
                    .mix(analyzeResults.mqc)
                    .mix(analyzeResults.correlation)
                    .collect()
    
    // Collect all outputs and run MultiQC
    runMultiQC(multiqcFiles)
}