# Ongoing_Evol_Landraces
This repository contain the data and scripts used to  run the analysis included in the paper "Contemporary evolution of maize landraces and their wild relatives influenced by gene flow with modern maize varieties"

##Prerequisites

VCFtools 0.1.15
PLINK v1.9
R version 3.6.1
Admixture 1.3
FastTree
custom scripts available

###R packages:
ggplot2
SNPRelate

### Directories
The repository is divided in the figure included in the paper "Ongoing evolution of maize (Zea mays L.) landraces and their wild relatives: Gene flow with modern maize varieties"

1) SNP_Calling
This directory contains the scripts to perform the sequence alignment and varian calling

2) data
- bam_files: This directory contains the bam files (*.bam, *.bai) to run the gene flow analysis and to compute Tajima's D and the genetic diversity estimators: Theta and Pi.
- Zea_mays.AGPv4.2x_0.8_0.01.NewID.Rojas2019.recode.vcf: The vcf  used to compute Fst between samples collected at different periods
- PCA_Rojas_Romay_CMLs_LRs_MVs: The plink file use to perform the PCAs with the [Romay et al., 2013 data set](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-6-r55).It includes, the landraces and modern varieties genotypes in B73_AGPv2.

2) Ancestry_analysis_Admixture
3) Gene_flow_test_ABBA-BABA
4) FST_index_across_genome
5) PCA
6) Nucleotide_diversity

Each directory contains subdirectories
- bin
- meta

###Contact
Idalia Rojas
icrojasb@gmail.com
