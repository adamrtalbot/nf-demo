# nf-demo

A Nextflow pipeline for demonstrating the value of Nextflow to R users.

## Overview

This pipeline has three variations on 3 branches:

1. `bash`: An R script manually called via a `bash` script
2. `nextflow`: A Nextflow pipeline
3. `platform`: A Nextflow pipeline configured with additional features to make it fully compatible with Seqera Platform

This pipeline processes RNA-Seq count data, performs quality control, and generates visualizations including PCA plots, count distributions, and sample correlations but this is not it's primary purpose.

## Requirements

- Nextflow
- Conda or Docker or Singularity

OR:

- R with packages:
  - ggplot2
  - jsonlite
- MultiQC

## Usage
Run the pipeline with:

```bash
nextflow run adamrtalbot/nf-demo --input 'path/to/*.csv' --output_dir results
```

## Input

The pipeline expects CSV files with the following input:

```csv
,Sample1,Sample2,Sample3,Sample4,Sample5,Sample6
ENSG00000000001,1,1,0,1,0,0
ENSG00000000002,558,1568,445,513,374,781
```

## Output

The pipeline generates a MultiQC document.

