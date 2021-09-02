# Annotate CNVRs with ENSEMBL gene IDs 

#################################################
### Libraries ###################################
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)


#################################################
### STEP 1: Annotate CNVRs ######################

# Define variables; "GWAS_CNVR_del_ANNOVAR_input.txt" from GWAS/07_CNVR/continuous/deletion_only/extract_CNVR_del.R contains CNVR boundaries in ANNOVAR input format
input <- "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/09_CNVR/continuous/white_british/deletion_only/data/final/GWAS_CNVR_del_ANNOVAR_input.txt"
output <- "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/10_ANNOVAR/continuous/white_british/deletion_only/data/final/gene_annotation_ENSEMBL/by_phenotype_CNVR/GWAS_CNVR_del_ANNOVAR_ENSEMBL"

# Run ANNOVAR in terminal
system(paste("/home/cauwerx/scratch/cauwerx/softwares/ANNOVAR-20191024/annotate_variation.pl",
			 "--geneanno -dbtype ensGene",
			 "-out", output, 
		     "-build hg19", 
			 input,
			 "/home/cauwerx/scratch/cauwerx/softwares/ANNOVAR-20191024/humandb/"))	
unlink(paste0(output, ".log"))
