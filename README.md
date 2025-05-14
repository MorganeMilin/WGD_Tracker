[![DOI](https://zenodo.org/badge/855434198.svg)](https://doi.org/10.5281/zenodo.14672010) 

# WGD_Tracker
Created by: MILIN Morgane

Contributors: BURBAN Ewen, AINOUCHE Malika and SALMON Armel

Date: 2024-12-04

## Table of Contents
- [Introduction](#introduction)
- [Installation](#installation)
- [Usage](#usage)
- [Citation](#citation)

### Introduction
WGD_Tracker is a tool designed to ease (intra- and inter-) genomic comparisons for detecting whole genome duplication events and to provide WGD dating by Ks analyses. This tool is fully customizable, allowing users to select stringent or flexible parameters, detailed in this document, depending on their specific goals.

This tool includes five pipelines:
* The RBBH Pipeline identifies homologous gene pairs via reciprocal BLAST best hit (RBBH) analysis from a BLAST output file. Because polyploid genomes contain multiple copies of the same genes, this pipeline can also search for reciprocal blast best hits (RBH) and not just the best hits in order to highlight all putative duplicated gene copies.
* The Ks Pipeline allows the calculation of synonymous substitution rates (Nei & Gojobori model) between gene pairs.
* The Synteny Pipeline accurately identifies syntenic blocks even when dispersed.
* The Dotplot Pipeline generates graphical plots from the outputs of the three previous Pipelines (RBBH; Ks; or Synteny).
* Karyotype Pipeline generates a graphical representation of conserved syntenic blocks between two genomes as karyotypes, for each chromosome across both genomes.

### Installation 
WGD_Tracker does not require any installation; a simple git clone is enough. However, WGD_Tracker does require dependencies: Python3, Java, PAML, MACSE, R, Singularity and Parallel. However, if you only want to use some of the tool's pipelines, not all dependencies are necessary for them to work properly. The following table shows the dependencies required for each pipeline:

However, WGD_Tracker does **require dependencies**:
- RBBH Pipeline: Python and Parallel
- Ks Pipeline: Python, Java, PAML, MACSE, R, Singularity and Parallel
- Synteny Pipeline: Python
- Dotplot Pipeline: Python

Please find the recommended dependencies version, which was used when developing the tool: Python v3.9, Parallel 20190122, Java v1.8.0, PAML v4.9, MACSE v2.05, R v4.1.0, Singularity v3.8.0

Retrieve WGD_Tracker

	git clone https://github.com/MorganeMilin/WGD_Tracker.git 

Dependencies installation (Please note that you must have access to conda)

	cd WGD_Tracker/WGD_Tracker/
 	conda env create -p ./dependencies_conda_env -f ./wgd_tracker_dependencies.yml
	conda activate dependencies_conda_env
 	R
  	install.packages("tinytex")
   	tinytex::install_tinytex()
 	q()
   	conda deactivate

### Usage
Please consult the Manual.pdf for instructions on how to format your working directory folder and input files. The document also contains all the information you need to create the configuration file required to run WGD_Tracker. Additionally, it provides further clarification on the various output files.

#### To run the RBBH Pipeline:
	sbatch --cpus-per-task=2 "{PATH}"/WGD_Tracker/RBBH_Pipeline.txt "{PATH}"/file.config

#### To run the Ks Pipeline:
 	sbatch --cpus-per-task=10 "{PATH}"/WGD_Tracker/Ks_Pipeline.txt "{PATH}"/file.config

#### To run the Synteny Pipeline:
  	sbatch  "{PATH}"/WGD_Tracker/Synteny_Pipeline.txt "{PATH}"/file.config

#### To run the Dotplot Pipeline:
   	sbatch  "{PATH}"/WGD_Tracker/Dotplot_Pipeline.txt "{PATH}"/file.config
    
#### To run the Karyotype Pipeline:
   	sbatch  "{PATH}"/WGD_Tracker/Karyotype_Pipeline.txt "{PATH}"/file.config

### Citation

Milin, M. WGD_Tracker (v1.0). Zenodo https://doi.org/10.5281/zenodo.14672010 (2025).
