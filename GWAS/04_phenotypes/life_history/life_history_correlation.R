# Calculate the correlation across life history phenotypes

#################################################
### Libraries ###################################
library(data.table)
library(dplyr)


#################################################
### STEP 1: Load Data ###########################

# This files are the output of GWAS/04_phenotypes/life_history/life_history_extraction.R
all <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/05_phenotypes/life_history/white_british/data/final/raw/life_history_WB_raw_All.txt", header = T))


#################################################
### STEP 2: Calculate the correlations ##########

cor_all <- as.data.frame(cor(all[, -c(1)], use = "pairwise.complete.obs", method = "pearson"))


#################################################
### Save ########################################

fwrite(cor_all, "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/05_phenotypes/life_history/white_british/data/final/raw/correlation/cor_LH_WB_raw_All.txt", col.names = T, row.names = T, quote = F, sep = "\t")
