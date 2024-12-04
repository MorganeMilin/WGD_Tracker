# WGD_Tracker
Created by: Morgane MILIN

Date: 2024-12-04

## Table of Contents
- [Introduction](#introduction)
- [Installation](#installation)
- [Usage](#usage)
- [Citation](#citation)

### Introduction
WGD_Tracker is a tool designed to ease (intra- and inter-)genomic comparisons for detecting whole genome duplication events and providing datation estimations by Ks analyses. This tool is fully customizable, allowing users to select stringent or flexible parameters depending on their specific analytical objectives.

This tool includes four pipelines:
* The RBBH Pipeline identifies homologous gene pairs via reciprocal BLAST best hit (RBBH) analysis from a BLAST output file. Because polyploid genomes can contain multiple copies of the same genes, this pipeline can also search for reciprocal blast best hits (RBH) and not just the best hits in order to highlight all putative duplicated gene copies.
* The Ks Pipeline allows the calculation of synonymous substitution rates between gene pairs.
* The Synteny Pipeline accurately identifies synteny blocks even when they are fragmented.
* The Dotplot Pipeline generates graphical plots from the outputs of the three previous Pipelines (RBBH; Ks; or Synteny).

### Installation 
WGD_Tracker does not require any installation, a simple git clone is enough:

	git clone https://github.com/MorganeMilin/WGD_Tracker.git 

However, WGD_Tracker does **require dependencies**:
- RBBH Pipeline: Python and Parallel
- Ks Pipeline: Python, Java, PAML, MACSE, R, Singularity and Parallel
- Synteny Pipeline: Python
- Dotplot Pipeline: Python

Please find the recommended dependencies version, which was used when developing the tool: Python v3.9, Parallel 20190122, Java v1.8.0, PAML v4.9, MACSE v2.05, R v4.1.0, Singularity v3.8.0

### Dependency installation
Create a working directory (ex: your_wd) containing the WGD_Tracker directory (with all scripts) + your data folder (with data files with the same file id but the following different extensions .gff; .fasta; .cds; .masked; .blast)

	cd your_wd/

Connect to a cluster node: 

	srun --pty bash
 
conda environment creation: 

 	. /local/env/envconda3.sh
	conda env create -p ./dependencies_conda_env -f wgd_tracker_dependencies.yml
 	exit

Image sif for R dependencies within the WGD_Tracker directory:

	srun --pty bash
 	cd WGD_Tracker/
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

