# Genomic polymorphisms in chr 10 underlie variation in gene expression in people with cancer

The aim of this study is to demonstrate by means of an expression
quantitative trait loci (eQTL) analysis an association between SNPs in de genes of chromosome 10, and a difference in gene expression between cancer and non-cancer patients.
During the eQTL analysis, gene expression is compared by regression to SNPs in the genes of chromosome 10. This allows the effect of these SNPs on gene expression to be determined. For this research, a pipeline was created with which the analysis was performed.

## Clone the project
First, start by cloning this project using the following command:
git clone 

## Installation
Make sure all the packages from the requirements.txt have been
installed. STAR, GATK and Snakemake are all tools which needs to be installed on the server/PC as well. 

## SnakeMake
The snakemake is functional at the moment, but can only be used on the server from the HAN. This because the files that are being used are way too large to be put on this Git. Therefore only a couple files can be found on this Git. If you want to access all the other documents, you should login into the server of the HAN. 
Remember to activate the conda environment by executing the following command: conda activate gatk. Otherwise the Snakefile wont be able to run. 

## Additional information
The data that has been used to run the snakemake is the sample from one patient. This does not change the functionality of the script, because the fact the Snakefile can be used for multiple samples with a little change.
