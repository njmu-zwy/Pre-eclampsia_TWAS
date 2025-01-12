# Overview
This GitHub repository hosts code related to coloc (colocalization analysis) and fusion (fusion gene - associated analysis). These codes are designed to assist researchers in bioinformatics, genetics, and related fields with specific data analysis tasks, enabling them to uncover potential associations and variant information at the gene level.
## 01_coloc_new_gene.r
### Purpose: Likely used to initiate colocalization analysis for novel genes. It probably loads relevant datasets, sets up initial parameters, and starts the first steps of analyzing whether new genes are colocalized with certain genetic traits or markers.
### Usage: Given its name, it might be the starting point for a new colocalization study. In an R environment, you can run it using source("01_coloc_new_gene.r") after ensuring all dependencies are met.
## 01_coloc_split_gene.r
### Purpose: This script may be responsible for splitting gene data, perhaps to break down large genomic regions into more manageable subsets for colocalization analysis. This can simplify the computational process and improve the accuracy of the coloc results.
### Usage: Execute it in an R session with source("01_coloc_split_gene.r"). It may require pre - loaded gene data in a specific format.
## 02_get_COLOC_run.sh
### Purpose: As a shell script, it's designed to trigger the actual running of the colocalization analysis pipeline. It might call the R scripts related to coloc, pass necessary parameters, and manage the overall execution flow, potentially handling things like resource allocation and job scheduling on a Unix - like system.
### Usage: In a terminal with appropriate permissions, run bash 02_get_COLOC_run.sh. Make sure R is properly installed and configured on the system, and any required input data is in place.
## 02_get_FUSION_run.sh
### Purpose: Similar to the previous shell script but focused on the fusion gene analysis pipeline. It will start the series of steps to detect and analyze fusion genes, which could involve running multiple R scripts, data pre - processing, and result compilation.
### Usage: Run it via bash 02_get_FUSION_run.sh in a Unix - based terminal. Ensure that the underlying fusion - related R scripts and data sources are accessible.
