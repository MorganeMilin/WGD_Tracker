# WGD_Tracker
Created by: Morgane MILIN
Date: 2024-09-11

## Table of Contents

## Introduction
WGD_Tracker is a tool for performing genomic comparisons. This tool includes four pipelines: (1) RBBH_Pipeline, which identifies homologous gene pairs, (2) Ks_Pipeline, which calculates the synonymous substitution rate between gene pairs, (3) Synteny_Pipeline, which accurately identifies syntenic blocks even when fragmented, and (4) Dotplot_Pipeline, which generates graphical plots. This tool was primarily designed to detect and date whole genome duplication events. This tool is completely customizable, allowing for stringent or flexible parameters depending on your analysis goals.

## Installation 
WGD_Tracker does not require any installation, a simple git clone is enough, but it does require dependencies: Python3, Java, PAML, MACSE, R, Singularity and Parallel. However, if you only want to use some of the tool's pipelines, not all dependencies are necessary for them to work properly. 
- [RBBH_Pipeline](#RBBH_Pipeline)
- [RBBH_Pipeline]
- (#RBBH_Pipeline)
""Required dependencies:""
- RBBH_Pipeline: Python3 and Parallel
- Ks_Pipeline: Python3, Java, PAML, MACSE, R, Singularity and Parallel
- Synteny_Pipeline: Python3
- Dotplot_Pipeline: Python3
