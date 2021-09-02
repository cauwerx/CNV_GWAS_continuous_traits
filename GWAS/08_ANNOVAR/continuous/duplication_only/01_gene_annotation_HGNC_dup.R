# Annotate CNVRs with HGNC gene IDs 

#################################################
### Libraries ###################################
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)


#################################################
### STEP 1: Annotate CNVRs ######################

# Define variables; "GWAS_CNVR_dup_ANNOVAR_input.txt" from GWAS/07_CNVR/continuous/duplication_only/extract_CNVR_dup.R contains CNVR boundaries in ANNOVAR input format
input <- "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/09_CNVR/continuous/white_british/duplication_only/data/final/GWAS_CNVR_dup_ANNOVAR_input.txt"
output <- "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/10_ANNOVAR/continuous/white_british/duplication_only/data/final/gene_annotation_HGNC/by_phenotype_CNVR/GWAS_CNVR_dup_ANNOVAR_HGNC"

# Run ANNOVAR in terminal
system(paste("/home/cauwerx/scratch/cauwerx/softwares/ANNOVAR-20191024/annotate_variation.pl",
			 "--geneanno -dbtype refGene",
			 "-out", output, 
		     "-build hg19", 
			 input,
			 "/home/cauwerx/scratch/cauwerx/softwares/ANNOVAR-20191024/humandb/"))	
unlink(paste0(output, ".log"))
