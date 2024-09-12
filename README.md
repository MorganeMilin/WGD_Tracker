# WGD_Tracker
Created by: Morgane MILIN

Date: 2024-09-11

## Table of Contents
- [Introduction](#introduction)
- [Installation](#installation)
- [Usage](#usage)
- [Citation](#citation)

### Introduction
WGD_Tracker is a tool designed to facilitate genome comparisons and detect and date whole genome duplication events. This tool is fully customizable, allowing users to select stringent or flexible parameters depending on their specific analytical objectives.
This tool comprises four pipelines:
* The RBBH Pipeline identifies homologous gene pairs via reciprocal BLAST best hit (RBBH) analysis on a BLAST output file. However, polyploid genomes can contain multiple copies of the same genes, so WGD_Tracker can also search for reciprocal blast best hits (RBH) and not just the best hits in order to highlight all duplicated gene copies.
* Ks Pipeline allows the calculation of the synonymous substitution rate between gene pairs.
* Synteny Pipeline accurately identifies synteny blocks even when they are fragmented.
* Dotplot Pipeline generates graphical plots.

### Installation 
WGD_Tracker does not require any installation, a simple git clone is enough:

	git clone https://github.com/MorganeMilin/WGD_Tracker.git 

However, WGD_Tracker does **require dependencies**:
- RBBH Pipeline: Python and Parallel
- Ks Pipeline: Python, Java, PAML, MACSE, R, Singularity and Parallel
- Synteny Pipeline: Python
- Dotplot Pipeline: Python
Please find the recommended dependencies version, which was used when developing the tool: Python v3.9, Parallel 20190122, Java v1.8.0, PAML v4.9, MACSE v2.05, R v4.1.0, Singularity v3.8.0

### Usage
Please consult the Manual.pdf for instructions on how to format your working directory folder and input files. The document also contains all the information you need to create the configuration file required to run WGD_Tracker. Additionally, it provides further clarification on the various output files.

#### To run the RBBH Pipeline:
	bash --cpus-per-task=2 "{PATH}"/WGD_Tracker/RBBH_Pipeline.txt "{PATH}"/file.config

#### To run the Ks Pipeline:
 	bash --cpus-per-task=2 "{PATH}"/WGD_Tracker/Ks_Pipeline.txt "{PATH}"/file.config

#### To run the Synteny Pipeline:
  	bash "{PATH}"/WGD_Tracker/Synteny_Pipeline.txt "{PATH}"/file.config

#### To run the Dotplot Pipeline:
   	bash "{PATH}"/WGD_Tracker/Dotplot_Pipeline.txt "{PATH}"/file.config

### Citations

