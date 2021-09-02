# Continuous trait association with duplication burden - WITH correction for CNV-GWAS

#################################################
### Libraries ###################################
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)


#################################################
### STEP 1: Load data ###########################

# Load duplication burden and select the Mb burden; from Burden_analysis/CNV_burden/02_DUP_burden_calculation.R
cnv_burden <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/05_phenotypes/CNV_burden/white_british/data/final/DUP_burden.txt.gz", header = T, select = c(1,2,4))) 
print(paste0("Number of samples with duplication burden: ", nrow(cnv_burden)))

# Load INT covariate-corrected continuous phenotypes - sex-agnostic; from GWAS/04_phenotypes/continuous/phenotype_extraction.R
pheno <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/05_phenotypes/continuous/white_british/data/final/INT_age_age2_batch_PCs/pheno_continuous_WB_INT_age_age2_sex_batch_PCs_All.txt", header = T))

# Load and add INT covariate-corrected continuous phenotypes - male-specific; from GWAS/04_phenotypes/continuous/phenotype_extraction.R
pheno_M <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/05_phenotypes/continuous/white_british/data/final/INT_age_age2_batch_PCs/pheno_continuous_WB_INT_age_age2_batch_PCs_M.txt", header = T))
pheno <- left_join(pheno, pheno_M[, c("IID", "balding", "facial_hair")], by = "IID"); rm(pheno_M)

# Load and add INT covariate-corrected continuous phenotypes - female-specific; from GWAS/04_phenotypes/continuous/phenotype_extraction.R
pheno_F <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/05_phenotypes/continuous/white_british/data/final/INT_age_age2_batch_PCs/pheno_continuous_WB_INT_age_age2_batch_PCs_F.txt", header = T))
pheno <- left_join(pheno, pheno_F[, c("IID", "menarche", "birth_weight_first_child", "menopause")], by = "IID"); rm(pheno_F)
print(paste0("Number of phenotype: ", ncol(pheno)-1))
print(paste0("Number of samples with phenotype: ", nrow(pheno)))

# Load CNVR covariates and set deletion to 0; from Burden_analysis/CNVR_covariates/CNVR_covariates.R 
cov <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/04_covariates/data/final/covariates_CNVR.txt.gz", header = T))
print(paste0("Number of samples with covariates: ", nrow(cov)))
cov[cov == -1] <- 0

# Load duplication CNVRs & remodel; from GWAS/07_CNVR/continuous/duplication_only/extract_CNVR_dup.R
cnvr <- as.data.table(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/09_CNVR/continuous/white_british/duplication_only/data/final/GWAS_CNVR_dup.txt", header = T, select = c(2,9,10,7), col.names = c("CHR", "START", "END", "PHENO"))) 
cnvr <- cnvr[which(cnvr$CHR != "X"), ] 									# Retain autosomal CNVs
cnvr$CHR <- as.numeric(cnvr$CHR)
cnvr <- cnvr[!duplicated(cnvr), ] 										# Remove duplicates
cnvr$ID <- paste0("chr", cnvr$CHR, "_", cnvr$START, "_", cnvr$END)		# Add unique ID
print(paste0("Number of unique autosomal duplication CNVRs: ", nrow(cnvr)))

# Load all CNVs and retain high confidence autosomal duplications; "ukb_cnv_global.gz" is the final PennCNV output, with all CNV calls in a linear format
cnvs <- as.data.table(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/cnv_calls/final/ukb_cnv_global.gz", header = T, select = c(1:5, 18), col.names = c("IID", "CHR", "START", "END", "CN", "QS")))
cnvs$IID <- sub("_.*", "", cnvs$IID)
cnvs$IID <- as.numeric(cnvs$IID)
cnvs <- cnvs[!cnvs$CHR %in% c(23, 24), ]								# Retain autosomal CNVs
cnvs <- cnvs[which(cnvs$QS >= 0.5), ] 									# Filter for quality score
print(paste0("Number of retained autosomal duplications: ", nrow(cnvs)))


#################################################
### STEP 2: Overlap CNVs with CNVRs #############

# Set key
setkey(cnvr, CHR, START, END)

# Overlap between duplications and duplication CNVRs
overlap <- na.omit(as.data.frame(foverlaps(cnvs, cnvr, maxgap = 0L, minoverlap = 1L, type = "any", mult = "all")))
overlap <- overlap[, c(5:9)]
overlap <- overlap[!duplicated(overlap), ]
print(paste0("Number of overlaps: ", nrow(overlap)))


#################################################
### STEP 3: CNV burden analysis #################

# Create empty dataframe to store results
results <- data.frame()

# Loop over phenotypes; burdens
for(p in names(pheno)[2:ncol(pheno)]) {

	# Define the phenotype
	print(paste0("Testing for association between ", p, " and duplication burden"))

	# Identify relevant covariates and correct both phenotype and burden for them
	if (nrow(cnvr[which(cnvr$PHENO == p), ]) > 0) { 
		print(paste0("Covariates: ", cnvr[which(cnvr$PHENO == p), "ID"]$ID))

		# 1) Correct the phenotype
		df_cov <- cov[, colnames(cov) %in% c("IID", cnvr[which(cnvr$PHENO == p), "ID"]$ID)]
		df_cov <- na.omit(left_join(df_cov, pheno[, c("IID", p)], by = "IID"))
		colnames(df_cov)[ncol(df_cov)] <- "phenotype"
		print(paste0("Number of individuals with CNVR-corrected ", p, ": ", nrow(df_cov)))
		df_cov[,"phenotype_CNVR"] <- residuals(lm(phenotype ~ . , data = df_cov[, -c(1)], na.action = na.exclude))

		# 2) Correct the burden 
		df_mb <- overlap[overlap$ID %in% cnvr[which(cnvr$PHENO == p), "ID"]$ID, ]
		print(paste0("Number of overlapping CNVs with covariate CNVRs: ", nrow(df_mb)))
		df_mb <- aggregate(TO_RM ~ IID, df_mb, sum)
		df_mb <- left_join(cnv_burden, df_mb, by = "IID")
		df_mb[is.na(df_mb$TO_RM), "TO_RM"] <- 0
		df_mb$BURDEN_Mb_COR <- df_mb$BURDEN_Mb - df_mb$TO_RM
		
		# Make a temporary dataframe
		df <- na.omit(left_join(df_mb[, c(1,2,5)], df_cov[, c("IID", "phenotype_CNVR")], by = "IID"))
		colnames(df)[ncol(df)] <- p
		print(paste0("Number of individuals analyzed: ", nrow(df)))
	
	# If no CNVR covariates are present (= no CNV-GWAS signal for that phenotype)
	} else {
		print(paste0("No covariates for ", p))
	
		# Make a temporary dataframe
		df <- na.omit(left_join(cnv_burden, pheno[, c("IID", p)], by = "IID"))
		colnames(df)[3] <- "BURDEN_Mb_COR"
		print(paste0("Number of individuals analyzed: ", nrow(df)))		
	}

	# Fitt a linear regression between phenotype and burden
	fit <- lm(df[, p] ~ df[, "BURDEN_Mb_COR"])
	print(summary(fit))

	# Fill in table
	results_temp <- data.frame()
	results_temp[1, "PHENO"] <- p 
	results_temp[1, "BURDEN"] <- "BURDEN_Mb_COR"
	results_temp[1, "N"] <- nrow(df)
	results_temp[1, "BETA"] <- as.numeric(coef(summary(fit))[2,1])
	results_temp[1, "SE"] <- as.numeric(coef(summary(fit))[2,2])
	results_temp[1, "T"] <- as.numeric(coef(summary(fit))[2,3])
	results_temp[1, "P"] <- as.numeric(coef(summary(fit))[2,4])
	results <-rbind(results, results_temp)	 			
}


#################################################
### STEP 4: Detect significant associations #####

#  Calculate significance threshold (FWER 5%: 57 continuous + 6 life history traits)
thr <- 0.05/63

# Print results
print(paste0("There are ", nrow(results[which(results[,7] <= thr), ])," significant associations across ", length(unique(results[which(results[,7] <= thr), "PHENO"]))," phenotypes"))


#################################################
### Save ########################################

fwrite(results, "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/12_burden_analysis/continuous/white_british/duplication_only/data/final/SumStat_CNVR_Mb/SumStat_DUP_burden_CNVR_Mb_continuous_All.txt", col.names = T, row.names = F, quote = F, sep = "\t")
