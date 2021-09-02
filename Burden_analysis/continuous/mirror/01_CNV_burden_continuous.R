# Continuous trait association with CNV burden - no correction for CNV-GWAS

#################################################
### Libraries ###################################
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)


#################################################
### STEP 1: Load data ###########################

# Load CNV burden; from Burden_analysis/CNV_burden/01_CNV_burden_calculation.R 
cnv_burden <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/05_phenotypes/CNV_burden/white_british/data/final/CNV_burden.txt.gz", header = T)) 
print(paste0("Number of samples with CNV burden: ", nrow(cnv_burden)))

# Load INT covariate-corrected continuous phenotypes - sex-agnostic; from GWAS/04_phenotypes/continuous/phenotype_extraction.R
pheno <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/05_phenotypes/continuous/white_british/data/final/INT_age_age2_batch_PCs/pheno_continuous_WB_INT_age_age2_sex_batch_PCs_All.txt", header = T))
print(paste0("Number of phenotypes: ", ncol(pheno)-1))
print(paste0("Number of samples with phenotype: ", nrow(pheno)))

# Load INT covariate-corrected continuous phenotypes - male-specific; from GWAS/04_phenotypes/continuous/phenotype_extraction.R
pheno_M <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/05_phenotypes/continuous/white_british/data/final/INT_age_age2_batch_PCs/pheno_continuous_WB_INT_age_age2_batch_PCs_M.txt", header = T))
print(paste0("Number of phenotypes (males): ", ncol(pheno_M)-1))
print(paste0("Number of male samples with phenotype: ", nrow(pheno_M)))

# Load INT covariate-corrected continuous phenotypes - female-specific; from GWAS/04_phenotypes/continuous/phenotype_extraction.R
pheno_F <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/05_phenotypes/continuous/white_british/data/final/INT_age_age2_batch_PCs/pheno_continuous_WB_INT_age_age2_batch_PCs_F.txt", header = T))
print(paste0("Number of phenotypes (female): ", ncol(pheno_F)-1))
print(paste0("Number of female samples with phenotype: ", nrow(pheno_F)))

# Combine to get a list of all assessed phenotypes
pheno_names <- unique(c(colnames(pheno)[-1], colnames(pheno_M)[-1], colnames(pheno_F)[-1]))
print(paste0("Number of analyzed phenotype (total): ", length(pheno_names)))


#################################################
### STEP 2: CNV burden analysis #################

# Create empty dataframe to store results
results <- data.frame()

# Loop over phenotypes; sexes; burdens
for(p in pheno_names) {

	# Define the phenotypes
	print(paste0("Testing for association between ", p, " and CNV burden"))

	# Loop over sexes (1 = All, 2 = Males, 3 = Females)
	for (s in 1:3) {

		# Stop the loop for non-compatible sex-phenotype combinations
		if (s == 1 & p %in% c("balding", "facial_hair", "menarche", "birth_weight_first_child", "menopause")) {
			print(paste0("--> Skip: ", p, " - All")); next}
		if (s == 2 & !(p %in% c("balding", "facial_hair"))) {
			print(paste0("--> Skip: ", p, " - Males")); next}
		if (s == 3 & !(p %in% c("menarche", "birth_weight_first_child", "menopause"))) {
			print(paste0("--> Skip: ", p, " - Females")); next}

		# Merge burden and appropriate phenotype data
		if (s == 1) {
			df <- na.omit(left_join(cnv_burden, pheno[, c("IID", p)], by = "IID"))
			sex <- "All"; print(sex)}
		if (s == 2) {
			df <- na.omit(left_join(cnv_burden, pheno_M[, c("IID", p)], by = "IID"))
			sex <- "males"; print(sex)}
		if (s == 3) {
			df <- na.omit(left_join(cnv_burden, pheno_F[, c("IID", p)], by = "IID"))
			sex <- "females"; print(sex)}
	
		# Loop over burden types (Mb or genes)
		for(b in 3:(ncol(df)-1)) {
		
			# Print analyzed burden
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
}
		

#################################################
### STEP 3: Detect significant associations #####

#  Calculate significance threshold (FWER 5%: 57 continuous + 6 life history traits)
thr <- 0.05/63

# Print output
print(paste0("There are ", nrow(results[which(results[,7] <= thr), ])," significant associations across ", length(unique(results[which(results[,7] <= thr), "PHENO"]))," phenotypes"))
print(unique(results[which(results[,7] <= thr), "PHENO"]))


#################################################
### Save ########################################

fwrite(results, "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/12_burden_analysis/continuous/white_british/mirror/data/final/SumStat/SumStat_CNV_burden_continuous.txt", col.names = T, row.names = F, quote = F, sep = "\t")
