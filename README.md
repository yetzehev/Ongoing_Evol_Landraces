# Ongoing_Evol_Landraces
This repository contains the data and scripts used to  run the analysis included in the paper "Contemporary evolution of maize landraces and their wild relatives influenced by gene flow with modern maize varieties"

## Prerequisites

VCFtools 0.1.15

PLINK v1.9

R version 3.6.1

Admixture 1.3

FastTree

angsd version: 0.921-9-ga30d0b4

custom scripts available at the github repository:

- https://github.com/owensgl/reformat

- https://github.com/owensgl/pop_gen


### R packages:
ggplot2

SNPRelate

permute

lattice

vegan

RColorBrewer

### Directories
The repository is divided in the figures included in the paper "Contemporary evolution of maize landraces and their wild relatives influenced by gene flow with modern maize varieties" https://doi.org/10.1073/pnas.1817664116

1) SNP_Calling
This directory contains the scripts to perform the sequence alignment and single nucleotide polymorphism calling

2) data available at (https://doi.org/10.17605/OSF.IO/PQVT4)
- bam_files: This directory contains the bam files (*.bam, *.bai) to run the gene flow analysis and to compute Tajima's D and the genetic diversity estimators: Theta and Pi.
- Zea_mays.AGPv4.2x_0.8_0.01.NewID.Rojas2019.recode.vcf: The vcf  used to compute Fst between samples collected at different periods
- PCA_Rojas_Romay_CMLs_LRs_MVs: The plink file use to perform the PCAs with the [Romay et al., 2013 data set](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-6-r55).It includes, the landraces and modern varieties genotypes in B73_AGPv2.

2) Ancestry_analysis
Contains the commands used to run admixture.

3) Gene_flow_test_ABBA-BABA
Contains the commands description used to run gene flow models, and examples of metadata files. All the meta data files can be generated with the information available at meta/metadata_385_taxa_Rojas_etal_2019.txt

4) FST_index_across_genome
Contains the executable files to compute Fst per site and window, and one markdown file with examples to run the analysis

5) PCA
Contains the html files thath show how to perform the PCA analysis with the files saved at data directory

6) Nucleotide_diversity
Conatains a text file that explains the steps to compute pi, theta and Tajima's D with the bam files provided.



### Contact
Idalia Rojas

icrojasb@gmail.com

icrojasb@iecologia.unam.mx
