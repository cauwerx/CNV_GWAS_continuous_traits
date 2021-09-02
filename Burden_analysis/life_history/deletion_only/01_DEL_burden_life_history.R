# Life history traits association with deletion burden

#################################################
### Libraries ###################################
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(matrixStats)


#################################################
### STEP 1: Load data ###########################

# Load deletion burden; from Burden_analysis/CNV_burden/03_DEL_burden_calculation.R 
cnv_burden <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/05_phenotypes/CNV_burden/white_british/data/final/DEL_burden.txt.gz", header = T)) 
print(paste0("Number of samples with DEL burden: ", nrow(cnv_burden)))

# Load INT covariate-corrected life history phenotypes; from GWAS/04_phenotypes/life_history/life_history_extraction.R
pheno <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/05_phenotypes/life_history/white_british/data/final/INT_age_age2_batch_PCs/life_history_WB_INT_age_age2_sex_batch_PCs_All.txt", header = T))
print(paste0("Number of samples with life history phenotype: ", nrow(pheno)))


#################################################
### STEP 2: CNV burden analysis #################

# Create empty dataframe to store results
results <- data.frame()

# Loop over phenotypes; burdens
for(p in names(pheno)[2:ncol(pheno)]) {

	# Define the phenotype
	print(paste0("Testing for association between ", p, " and deletion burden"))

	# Make a temporary datframe 
	df <- left_join(cnv_burden, pheno[, c("IID", p)], by = "IID")
	df <- df[which(!is.na(df[, p])), ]
	print(paste0("Number of individuals analyzed: ", nrow(df)))

	# Loop over burdens (Mb and genes)
	for(b in 3:(ncol(df)-1)) {
		print(paste0("#### ", names(df)[b], " ####"))

		# Fit a linear regression
		fit <- lm(df[, p] ~ df[, b])
		print(summary(fit))

		# Fill in table
		results_temp <- data.frame()
		results_temp[1, "PHENO"] <- p 
		results_temp[1, "BURDEN"] <- names(df)[b]
		results_temp[1, "N"] <- nrow(df)
		results_temp[1, "BETA"] <- as.numeric(coef(summary(fit))[2,1])
		results_temp[1, "SE"] <- as.numeric(coef(summary(fit))[2,2])
		results_temp[1, "T"] <- as.numeric(coef(summary(fit))[2,3])
		results_temp[1, "P"] <- as.numeric(coef(summary(fit))[2,4])
		results <-rbind(results, results_temp)	 	
	}
}


#################################################
### Save ########################################

fwrite(results, "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/12_burden_analysis/life_history/white_british/deletion_only/data/final/summary_stat_DEL_burden_LH.txt", col.names = T, row.names = F, quote = F, sep = "\t")
