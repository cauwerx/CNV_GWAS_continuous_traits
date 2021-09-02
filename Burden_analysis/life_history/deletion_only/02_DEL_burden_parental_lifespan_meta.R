# Meta-analyse maternal and paternal lifespan into parental lifespan
# Burden analysis on parental lifespan 

#################################################
### Libraries ###################################
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(matrixStats)


#################################################
### STEP 1: Correct raw life history traits #####

# Load RAW life history phenotypes; from GWAS/04_phenotypes/life_history/life_history_extraction.R
pheno <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/05_phenotypes/life_history/white_british/data/final/raw/life_history_WB_raw_All.txt", header = T, select = c(1,4,5)))

# Load standard covariates and select unrelated white British individuals; from GWAS/03_covariates/covariate_extraction.R
cov <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/04_covariates/data/final/covariates.txt", header = T, select = c(1:4, 6:46)))
colnames(cov)[1] <- "IID"
cov <- cov[cov$IID %in% pheno$IID, ]

# Create an empty dataframe
pheno_cor <- data.frame(IID = pheno$IID)

# Loop over phenotypes
for (p in 2:ncol(pheno)) {
	
	# Define the phenotype
	phenotype <- colnames(pheno)[p]
	print(paste0("Phenotype: ", phenotype))

	# Regress out covariates
	temp <- right_join(cov, pheno[, c(1,p)], by = "IID")
	print(paste0("Number of covariates: ", ncol(temp)-2))
	colnames(temp)[ncol(temp)] <- "PHENO"
	pheno_cor[, phenotype] <- residuals(lm(PHENO ~ . , data = temp[, -c(1)], na.action = na.exclude))
}; rm (p, phenotype, temp)


#################################################
### STEP 2: CNV burden analysis ##################

# Load deletion burden; from Burden_analysis/CNV_burden/03_DEL_burden_calculation.R 
cnv_burden <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/05_phenotypes/CNV_burden/white_british/data/final/DEL_burden.txt.gz", header = T)) 
print(paste0("Number of samples with DEL burden: ", nrow(cnv_burden)))

# Create empty dataframe to store results 
results <- data.frame()

# Loop over phenotypes; burdens
for(p in names(pheno_cor)[2:ncol(pheno_cor)]) {

	# Define the phenotype
	print(paste0("Testing for association between ", p, " and deletion burden"))

	# Make a temporary dataframe
	df <- left_join(cnv_burden, pheno_cor[, c("IID", p)], by = "IID")
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

# Save
fwrite(results, "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/12_burden_analysis/life_history/white_british/deletion_only/data/final/meta_analysis/association_original_scale/summary_stat_DEL_burden_parental_lifespan_original_scale.txt", col.names = T, row.names = F, quote = F, sep = "\t")


#################################################
### STEP 3: Meta-analysis #######################

# Create a GWAMA input file for "mother lifespan" (i.e mother's age at death)
gwama_mother_input <- data.frame(MARKERNAME = results[which(results$PHENO == "lifespan_mother"), "BURDEN"], EA = "A", NEA = "T", BETA = results[which(results$PHENO == "lifespan_mother"), "BETA"], SE = results[which(results$PHENO == "lifespan_mother"), "SE"], N = results[which(results$PHENO == "lifespan_mother"), "N"]) 
gwama_mother_file <- "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/12_burden_analysis/life_history/white_british/deletion_only/data/temp/GWAMA_input/mother_lifespan_GWAMA_input.txt"
fwrite(gwama_mother_input, gwama_mother_file, col.names = T, row.names = F, quote = F, sep = "\t") 
 
# Create a GWAMA input file for "paternal lifespan" (i.e father's age at death)
gwama_father_input <- data.frame(MARKERNAME = results[which(results$PHENO == "lifespan_father"), "BURDEN"], EA = "A", NEA = "T", BETA = results[which(results$PHENO == "lifespan_father"), "BETA"], SE = results[which(results$PHENO == "lifespan_father"), "SE"], N = results[which(results$PHENO == "lifespan_father"), "N"]) 
gwama_father_file <- "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/12_burden_analysis/life_history/white_british/deletion_only/data/temp/GWAMA_input/father_lifespan_GWAMA_input.txt"
fwrite(gwama_father_input, gwama_father_file, col.names = T, row.names = F, quote = F, sep = "\t") 
 
# Create a GWAMA file list
gwama_file_list <- "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/12_burden_analysis/life_history/white_british/deletion_only/data/temp/GWAMA_input/GWAMA_input_file_list.txt" 
fwrite(data.frame(FILE = c(gwama_mother_file, gwama_father_file)), gwama_file_list, col.names = F, row.names = F, quote = F, sep = "\t")

# Determine the output path
gwama_output  <- "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/12_burden_analysis/life_history/white_british/deletion_only/data/final/meta_analysis/GWAMA/DEL_burden_parental_lifespan_original_meta"

#Run GWAMA 
system(paste("/home/cauwerx/scratch/cauwerx/softwares/GWAMA-2.2.2/GWAMA", 
				"--filelist", gwama_file_list,
				"--output", gwama_output, 
				"--quantitative" ))
