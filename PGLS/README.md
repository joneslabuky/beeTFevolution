# Snakemake workflow for PGLS: testing TF motif associations with sociality
Snakemake is used to manage this workflow and the files in this directory control the pipeline to run phylogenetic least squares (PGLS) models to test associations between TF motif presence in promoter regions and social behavior across bee species. 

This pipeline uses motif count tables from upstream workflows, species tree from OrthoFinder, and a metadata table defining sociality for each species. 

**Config file: `config_PGLS.yaml`**

This configuration file defines all paths used in the PGLS workflow. 

**Snakefile: `Snakefile_PGLS`**

This workflow runs the `run_PGLS.R` script for each motif table ultimately producing one result file per motif. This workflow implements these major steps: 

Reads motif count tables, filters orthogroups with sufficient data, and runs script.

**Script: `run_PGLS.R`**

This R script contails all model logic to load phylogeny and species sociality metadata, MAD MeanAD scale and min-max normalize data, split data into two comparisons for PGLS1 (social vs ancestrally solitary species) and PGLS2 (social vs secondarily solitary species), and applies multiple lambda bounds to optimize model convergence.
