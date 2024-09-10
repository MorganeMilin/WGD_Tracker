# WGD_Tracker
Created by: Morgane MILIN
Date: 2024-09-11

## Table of Contents
- [Introduction](#introduction)
- [Installation](#installation)
- [Inputs](#inputs)
- [Outputs](#outputs)
- [Usage](#usage)
- [Citations](#citations)

### Introduction
WGD_Tracker is a tool for performing genomic comparisons. This tool includes four pipelines: (1) RBBH_Pipeline, which identifies homologous gene pairs, (2) Ks_Pipeline, which calculates the synonymous substitution rate between gene pairs, (3) Synteny_Pipeline, which accurately identifies syntenic blocks even when fragmented, and (4) Dotplot_Pipeline, which generates graphical plots. This tool was primarily designed to detect and date whole genome duplication events. This tool is completely customizable, allowing for stringent or flexible parameters depending on your analysis goals.

### Installation 
WGD_Tracker does not require any installation, a simple git clone is enough, but it does **require dependencies**:
- RBBH_Pipeline: Python3 and Parallel
- Ks_Pipeline: Python3, Java, PAML, MACSE, R, Singularity and Parallel
- Synteny_Pipeline: Python3
- Dotplot_Pipeline: Python3

### Inputs
  #### RBBH_Pipeline
Blast output file 

  #### Ks_Pipeline
the output file from the RBBH_Pipeline analysis 

  #### Synteny_Pipeline
the output files from the RBBH_Pipeline and Ks_Pipeline analysis

  #### Dotplot_Pipeline
Any outfile from RBBH, Ks or Synteny Pipelines + another configuration file for the axes

### Outputs

### Usage

### Citations
