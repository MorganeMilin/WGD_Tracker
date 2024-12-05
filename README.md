# WGD_Tracker
Created by: Morgane MILIN

Date: 2024-12-04

## Table of Contents
- [Introduction](#introduction)
- [Installation](#installation)
- [Usage](#usage)
- [Citation](#citation)

### Introduction
WGD_Tracker is a tool designed to ease (intra- and inter-)genomic comparisons for detecting whole genome duplication events and to provide dating estimations by Ks analyses. This tool is fully customizable, allowing users to select stringent or flexible parameters depending on their specific analytical objectives.

This tool includes five pipelines:
* The RBBH Pipeline identifies homologous gene pairs via reciprocal BLAST best hit (RBBH) analysis from a BLAST output file. Because polyploid genomes can contain multiple copies of the same genes, this pipeline can also search for reciprocal blast best hits (RBH) and not just the best hits in order to highlight all putative duplicated gene copies.
* The Ks Pipeline allows the calculation of synonymous substitution rates between gene pairs.
* The Synteny Pipeline accurately identifies synteny blocks even when they are fragmented.
* The Dotplot Pipeline generates graphical plots from the outputs of the three previous Pipelines (RBBH; Ks; or Synteny).
* The Karyotype Pipeline generates graphical representation of the two analyzed genomes karyotypes, illustrating the detected synteny blocks for each chromosome across both genomes.

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

Conda environment creation: (Please note that you must have access to conda)

	conda env create -p ./dependencies_conda_env -f wgd_tracker_dependencies.yml$ cd ./WGD_Tracker


Image sif for R dependencies: (Please note that you must have access to singularity and the generatedd .sif image must be located in the WGD_Tracker folder)

	cd WGD_Tracker 
	singularity build ./rmarkdown.sif ./wgd_tracker_R.def

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

### Citations

