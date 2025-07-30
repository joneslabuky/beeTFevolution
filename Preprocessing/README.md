# Snakemake workflow for preprocessing

Snakemake is used to manage this workflow and the two files control the pipeline to prepare orthogroup and transcription factor data across multiple species for downstream comparative analyses.

**Config file: `config.yaml`**

This configuration file defines all paths, file naming, and containerized tool versions for the workflow.

`paths:` defines locations and file patterns to be used for: 

Input data (genomes, annotations, peptide FASTAS, and chromosome lengths), intermediate outputs (longest isoforms, BED files, Bowtie indices, PWMSCan outputs), PWMScan integration (motif tag generation and alignment results), final tables (per-species motif count and score summarries).  

File patterns use placeholders `{sample}` and `{motif}`, which are then populated dynamically by Snakemaking through the provided sample and motif lists. 

`containers`: specifies versions of all tools run using Docker

**Snakefile: `Snakefile`**

This workflow implements these major steps in preparation for downstream analysis: 

Extracts longest isoforms per gene using AGAT, translates and cleans peptides to use as inputs for OrthoFinder, runs OrthoFinder to infer orthogroups across species, generates bowtie indices, scans promoter regions for TF motifs with PWMscan and genome-specific background compositions, aligns motif tag lists to promoter regions using Bowtie, and generates motif count and score tables for each species. 
