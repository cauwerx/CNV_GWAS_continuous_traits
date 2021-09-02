# Extract phenotype data and do the inverse normal transformation and covariate correction

################################################
### Libraries ##################################
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(matrixStats)


################################################
### STEP 1: Load Main Data #####################

# This corresponds to a .csv file with two columns: field ID and short name for the corresponding phenotype
pheno_list <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/05_phenotypes/continuous/white_british/data/raw/continuous_traits.csv", header = T))


#################################################
### STEP 2: Order raw phenotype files ###########

# As phenotypic data is spread throughout multiple file downloaded from the UKBB portal, we start by listing all phenotype files available and order them
raw_pheno_file <- file.info(list.files("/data/FAC/FBM/DBC/zkutalik/default_sensitive/UKBB/pheno", pattern = "^ukb", full.names = T, recursive = F))
raw_pheno_file <- data.frame(File = rownames(raw_pheno_file), Date = raw_pheno_file[,4])
raw_pheno_file <- raw_pheno_file[order(raw_pheno_file$Date, decreasing = T), ]
print(paste0("There are ", nrow(raw_pheno_file), " raw phenotype files"))


#################################################
### STEP 3: Identify most recent location #######
# As mentioned in STEP 2, there are multiple phenotype files, with some phenotypes being in several files
# For the purpose of the analysis, we use phenotype information originating from the most recent phenotype file

pheno_list$File <- NA
pheno_list$Col <- NA

# Loop over all phenotypes to detect the most recent location
for (p in 1:nrow(pheno_list)) {

	# Define the phenotype
	pheno <- pheno_list[p, "Pheno"] 
	ID <- pheno_list[p, "FieldID"] 
	print(paste0("Extracting most recent file for ", pheno))
	counter <- 1

	# Determine the most recent location --> loop through raw files
	while (is.na(pheno_list[p, "File"])) {
	
		header <- fread(as.character(raw_pheno_file[counter, "File"]), nrow = 0)
		col <- grep(paste0("^", ID, "-"), names(header))
		if (length(col) > 0) {
			pheno_list[p, "File"] <- as.character(raw_pheno_file[counter, "File"])
			print(paste0("Most recent location: ", sub(".*/", "", pheno_list[p, "File"])))
			pheno_list[p, "Col"] <- paste(col, collapse = "_")}
		counter <- counter +1
	}
} 
rm(p, pheno, ID, counter, header, col)
fwrite(pheno_list, "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/05_phenotypes/continuous/white_british/data/final/pheno_list.txt", col.names = T, row.names = F, quote = F, sep = "\t")


#################################################
### STEP 4: Extract phenotypes and average ######

# This file contains all UKBB sample eids, except for the redacted ones
phenotypes <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/general_data/UKBB/eid_no_redacted3.txt.gz", header = T, select = c(2)), drop = F)

# Loop over all phenotypes to extract the data and calculate the average over different measurement instances
for(p in 1:nrow(pheno_list)) {
	
	# Define the phenotype
	pheno <- pheno_list[p, "Pheno"] 
	ID <- pheno_list[p, "FieldID"] 
	file <- pheno_list[p, "File"] 
	col <- pheno_list[p, "Col"] 
	print(paste0("Extracting ", pheno, " from ", sub(".*/", "", file)))

	# Extract columns corresponding to the defined phenotype
	temp_pheno <- as.data.frame(fread(file, header = T, select = c(1, as.numeric(str_split(col, pattern = "_")[[1]]))))

	# Average (if more than one column present)
	if (ncol(temp_pheno) > 2) {
	temp_pheno <- data.frame(eid = temp_pheno[, 1], pheno = rowMeans(temp_pheno[, 2:ncol(temp_pheno)], na.rm = T))}

	# Merge
	colnames(temp_pheno) <- c("eid", pheno)
	phenotypes <- full_join(phenotypes, temp_pheno, by = "eid")

}
rm(pheno, ID, file, col, temp_pheno)
print(paste0("Dimensions of phenotype table: ", ncol(phenotypes), " pheno x ", nrow(phenotypes), " eids"))


#################################################
### STEP 5: Composite traits ####################

# Grip strength: mean left and right
phenotypes$GS <- (phenotypes$grip_strength_left + phenotypes$grip_strength_right)/2
phenotypes <- phenotypes[, !names(phenotypes) %in% c("grip_strength_left", "grip_strength_right")]

# WHR: waist circumference divided by hip circumference
phenotypes$WHR <- phenotypes$waist_circumference/phenotypes$hip_circumference
phenotypes <- phenotypes[, !names(phenotypes) %in% c("waist_circumference", "hip_circumference")]

# WHRadjBMI: waist circumference divided by hip circumference adjusted for BMI 
# --> Done after selecting samples

print(paste0("Dimensions of phenotype table after adding composite traits: ", ncol(phenotypes), " pheno x ", nrow(phenotypes), " eids"))


#################################################
### STEP 6: Replace entries by NA ###############

# -1: do not know; -2: only had twins; -3: prefer not to answer
phenotypes[phenotypes < 0] <- NA


#################################################
### STEP 7: Select samples & add WHRadjBMI ######
# Retain unrelated, white British samples selected during GWAS/01_samples/sample_filtering.R

# All individuals - exclude sex-specific traits
eid_all <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/01_samples/white_british/data/final/samples_white_british_All.txt", header = F, col.names = "eid"))
pheno_all <- phenotypes[phenotypes$eid %in% eid_all$eid, ]
pheno_all <- pheno_all[, !names(pheno_all) %in% c("facial_hair", "balding", "menarche", "menopause","birth_weight_first_child")]
pheno_all$WHRadjBMI <- residuals(lm(pheno_all$WHR ~ pheno_all$BMI, na.action = na.exclude))
colnames(pheno_all)[1] <- "IID"
print(paste0("Dimensions of phenotype table for all unrelated white British individuals: ", ncol(pheno_all), " pheno x ", nrow(pheno_all), " eids"))
fwrite(pheno_all, "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/05_phenotypes/continuous/white_british/data/final/raw/pheno_continuous_WB_raw_All.txt", col.names = T, row.names = F, quote = F, sep = "\t", na = NA)

# Males - only male-specific traits
eid_M <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/01_samples/white_british/data/final/samples_white_british_M.txt", header = F, col.names = "eid"))
pheno_M <- phenotypes[phenotypes$eid %in% eid_M$eid, ]
pheno_M <- pheno_M[, !names(pheno_M) %in% c("menarche", "menopause","birth_weight_first_child")]
pheno_M$WHRadjBMI <- residuals(lm(pheno_M$WHR ~ pheno_M$BMI, na.action = na.exclude))
colnames(pheno_M)[1] <- "IID"
print(paste0("Dimensions of phenotype table for unrelated white British males: ", ncol(pheno_M), " pheno x ", nrow(pheno_M), " males"))
fwrite(pheno_M, "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/05_phenotypes/continuous/white_british/data/final/raw/pheno_continuous_WB_raw_M.txt", col.names = T, row.names = F, quote = F, sep = "\t", na = NA)

# Females - only female-specific traits
eid_F <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/01_samples/white_british/data/final/samples_white_british_F.txt", header = F, col.names = "eid"))
pheno_F <- phenotypes[phenotypes$eid %in% eid_F$eid, ]
pheno_F <- pheno_F[, !names(pheno_F) %in% c("facial_hair", "balding")]
pheno_F$WHRadjBMI <- residuals(lm(pheno_F$WHR ~ pheno_F$BMI, na.action = na.exclude))
colnames(pheno_F)[1] <- "IID"
print(paste0("Dimensions of phenotype table for unrelated white British females: ", ncol(pheno_F), " pheno x ", nrow(pheno_F), " females"))
fwrite(pheno_F, "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/05_phenotypes/continuous/white_british/data/final/raw/pheno_continuous_WB_raw_F.txt", col.names = T, row.names = F, quote = F, sep = "\t", na = NA)


#################################################
### STEP 8: Summary #############################

# All individuals - no sex-specific traits
summary_pheno_all <- data.frame(Pheno = colnames(pheno_all)[2:ncol(pheno_all)],
								Num = colSums(!is.na(pheno_all[, c(2:ncol(pheno_all))])),
								NumNA = colSums(is.na(pheno_all[, c(2:ncol(pheno_all))])),
								Mean = colMeans(pheno_all[, c(2:ncol(pheno_all))], na.rm = T),
								Median = colMedians(as.matrix(pheno_all[, c(2:ncol(pheno_all))]), na.rm = T),
								SD = colSds(as.matrix(pheno_all[, c(2:ncol(pheno_all))]), na.rm = T),
								Min = colMins(as.matrix(pheno_all[, c(2:ncol(pheno_all))]), na.rm = T),
								Max = colMaxs(as.matrix(pheno_all[, c(2:ncol(pheno_all))]), na.rm = T))
fwrite(summary_pheno_all, "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/05_phenotypes/continuous/white_british/data/final/raw/summary_pheno_continuous_WB_raw_All.txt", col.names = T, row.names = F, quote = F, sep = "\t", na = NA)
print("Summary all: done!")

# Males - no female-specific traits
summary_pheno_M <- data.frame(Pheno = colnames(pheno_M)[2:ncol(pheno_M)],
							  Num = colSums(!is.na(pheno_M[, c(2:ncol(pheno_M))])),
							  NumNA = colSums(is.na(pheno_M[, c(2:ncol(pheno_M))])),
							  Mean = colMeans(pheno_M[, c(2:ncol(pheno_M))], na.rm = T),
							  Median = colMedians(as.matrix(pheno_M[, c(2:ncol(pheno_M))]), na.rm = T),
							  SD = colSds(as.matrix(pheno_M[, c(2:ncol(pheno_M))]), na.rm = T),
							  Min = colMins(as.matrix(pheno_M[, c(2:ncol(pheno_M))]), na.rm = T),
							  Max = colMaxs(as.matrix(pheno_M[, c(2:ncol(pheno_M))]), na.rm = T))
fwrite(summary_pheno_M, "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/05_phenotypes/continuous/white_british/data/final/raw/summary_pheno_continuous_WB_raw_M.txt", col.names = T, row.names = F, quote = F, sep = "\t", na = NA)
print("Summary males: done!")

# Females - no male-specific traits
summary_pheno_F <- data.frame(Pheno = colnames(pheno_F)[2:ncol(pheno_F)],
							  Num = colSums(!is.na(pheno_F[, c(2:ncol(pheno_F))])),
							  NumNA = colSums(is.na(pheno_F[, c(2:ncol(pheno_F))])),
							  Mean = colMeans(pheno_F[, c(2:ncol(pheno_F))], na.rm = T),
							  Median = colMedians(as.matrix(pheno_F[, c(2:ncol(pheno_F))]), na.rm = T),
							  SD = colSds(as.matrix(pheno_F[, c(2:ncol(pheno_F))]), na.rm = T),
							  Min = colMins(as.matrix(pheno_F[, c(2:ncol(pheno_F))]), na.rm = T),
							  Max = colMaxs(as.matrix(pheno_F[, c(2:ncol(pheno_F))]), na.rm = T))
fwrite(summary_pheno_F, "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/05_phenotypes/continuous/white_british/data/final/raw/summary_pheno_continuous_WB_raw_F.txt", col.names = T, row.names = F, quote = F, sep = "\t", na = NA)
print("Summary females: done!")


#################################################
### STEP 9: Inverse Normal Transformation #######

#  All individuals
pheno_int_all <- data.frame(IID = pheno_all$IID)
for (p in 2:ncol(pheno_all)) {
	pheno <- colnames(pheno_all)[p]
	print(paste0("INT - All: ", pheno))
	pheno_int_all[, pheno] <- qnorm((rank(pheno_all[, p], na.last = "keep")-0.5)/sum(!is.na(pheno_all[, p])))
}; rm (p, pheno)

fwrite(pheno_int_all, "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/05_phenotypes/continuous/white_british/data/final/INT/pheno_continuous_WB_INT_All.txt", col.names = T, row.names = F, quote = F, sep = "\t", na = NA)


# Males
pheno_int_M <- data.frame(IID = pheno_M$IID)
for (p in 2:ncol(pheno_M)) {
	pheno <- colnames(pheno_M)[p]
	print(paste0("INT - Males: ", pheno))
	pheno_int_M[, pheno] <- qnorm((rank(pheno_M[, p], na.last = "keep")-0.5)/sum(!is.na(pheno_M[, p])))
}; rm (p, pheno)

fwrite(pheno_int_M, "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/05_phenotypes/continuous/white_british/data/final/INT/pheno_continuous_WB_INT_M.txt", col.names = T, row.names = F, quote = F, sep = "\t", na = NA)


# Females
pheno_int_F <- data.frame(IID = pheno_F$IID)
for (p in 2:ncol(pheno_F)) {
	pheno <- colnames(pheno_F)[p]
	print(paste0("INT - Females: ", pheno))
	pheno_int_F[, pheno] <- qnorm((rank(pheno_F[, p], na.last = "keep")-0.5)/sum(!is.na(pheno_F[, p])))
}; rm (p, pheno)

fwrite(pheno_int_F, "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/05_phenotypes/continuous/white_british/data/final/INT/pheno_continuous_WB_INT_F.txt", col.names = T, row.names = F, quote = F, sep = "\t", na = NA)


#################################################
### STEP 10: Correct for covariates #############

# Read in covariates from GWAS/covariates/covariate_extraction.R
cov <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/04_covariates/data/final/covariates.txt", header = T, select = c(1:4, 6:46)))
colnames(cov)[1] <- "IID"

# Select subsets (retained samples; sex is not used when assessing sex-specific traits)
cov_all <- cov[cov$IID %in% eid_all$eid, ] # 331523 eids
print(paste0("Covariates all: ", ncol(cov_all)-1))
print(colnames(cov_all[-1]))
cov_M <- cov[cov$IID %in% eid_M$eid, -c(4)] # 152968 eids
print(paste0("Covariates males: ", ncol(cov_M)-1))  
print(colnames(cov_M[-1]))
cov_F <- cov[cov$IID %in% eid_F$eid, -c(4)] # 178555 eids
print(paste0("Covariates females: ", ncol(cov_F)-1))
print(colnames(cov_F[-1]))


# Regress out covariates - All individuals
pheno_int_cor_all <- data.frame(IID = pheno_int_all$IID)
for (p in 2:ncol(pheno_int_all)) {
	# Define phenotype
	pheno <- colnames(pheno_int_all)[p]
	print(paste0("Covariates - All: ", pheno))
	# Create a temporary dataframe
	temp <- right_join(cov_all, pheno_int_all[, c(1,p)], by = "IID")
	colnames(temp)[ncol(temp)] <- "pheno"
	# Regress out covariates
	pheno_int_cor_all[, pheno] <- residuals(lm(pheno ~ . , data = temp[, -c(1)], na.action = na.exclude))
}; rm (p, pheno, temp)

fwrite(pheno_int_cor_all, "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/05_phenotypes/continuous/white_british/data/final/INT_age_age2_batch_PCs/pheno_continuous_WB_INT_age_age2_sex_batch_PCs_All.txt", col.names = T, row.names = F, quote = F, sep = "\t", na = NA)


# Regress out covariates - Males
pheno_int_cor_M <- data.frame(IID = pheno_int_M$IID)
for (p in 2:ncol(pheno_int_M)) {
	# Define phenotype
	pheno <- colnames(pheno_int_M)[p]
	print(paste0("Covariates - Males: ", pheno))
	# Create a temporary dataframe
	temp <- right_join(cov_M, pheno_int_M[, c(1,p)], by = "IID")
	colnames(temp)[ncol(temp)] <- "pheno"
	# Regress out covariates
	pheno_int_cor_M[, pheno] <- residuals(lm(pheno ~ . , data = temp[, -c(1)], na.action = na.exclude))
}; rm (p, pheno, temp)

fwrite(pheno_int_cor_M, "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/05_phenotypes/continuous/white_british/data/final/INT_age_age2_batch_PCs/pheno_continuous_WB_INT_age_age2_batch_PCs_M.txt", col.names = T, row.names = F, quote = F, sep = "\t", na = NA)


# Regress out covariates - Females
pheno_int_cor_F <- data.frame(IID = pheno_int_F$IID)
for (p in 2:ncol(pheno_int_F)) {
	# Define phenotype
	pheno <- colnames(pheno_int_F)[p]
	print(paste0("Covariates - Females: ", pheno))
	# Create a temporary dataframe
	temp <- right_join(cov_F, pheno_int_F[, c(1,p)], by = "IID")
	colnames(temp)[ncol(temp)] <- "pheno"
	# Regress out covariates
	pheno_int_cor_F[, pheno] <- residuals(lm(pheno ~ . , data = temp[, -c(1)], na.action = na.exclude))
}; rm (p, pheno, temp)

fwrite(pheno_int_cor_F, "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/05_phenotypes/continuous/white_british/data/final/INT_age_age2_batch_PCs/pheno_continuous_WB_INT_age_age2_batch_PCs_F.txt", col.names = T, row.names = F, quote = F, sep = "\t", na = NA)



##########################################################
### STEP 11: Correct for covariates on non-INT data ######

# Regress out covariates - All
pheno_cor_all <- data.frame(IID = pheno_all$IID)
for (p in 2:ncol(pheno_all)) {
	# Define phenotype
	pheno <- colnames(pheno_all)[p]
	print(paste0("Covariates (no INT) - All: ", pheno))
	# Create a temporary dataframe
	temp <- right_join(cov_all, pheno_all[, c(1,p)], by = "IID")
	colnames(temp)[ncol(temp)] <- "pheno"
	# Regress out covariates
	pheno_cor_all[, pheno] <- residuals(lm(pheno ~ . , data = temp[, -c(1)], na.action = na.exclude))
}; rm (p, pheno, temp)

fwrite(pheno_cor_all, "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/05_phenotypes/continuous/white_british/data/final/age_age2_batch_PCs/pheno_continuous_WB_age_age2_sex_batch_PCs_All.txt", col.names = T, row.names = F, quote = F, sep = "\t", na = NA)

# Regress out covariates - Males
pheno_cor_M <- data.frame(IID = pheno_M$IID)
for (p in 2:ncol(pheno_M)) {
	# Define phenotype
	pheno <- colnames(pheno_M)[p]
	print(paste0("Covariates (no INT) - Males: ", pheno))
	# Create a temporary dataframe
	temp <- right_join(cov_M, pheno_M[, c(1,p)], by = "IID")
	colnames(temp)[ncol(temp)] <- "pheno"
	# Regress out covariates
	pheno_cor_M[, pheno] <- residuals(lm(pheno ~ . , data = temp[, -c(1)], na.action = na.exclude))
}; rm (p, pheno, temp)

fwrite(pheno_cor_M, "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/05_phenotypes/continuous/white_british/data/final/age_age2_batch_PCs/pheno_continuous_WB_age_age2_batch_PCs_M.txt", col.names = T, row.names = F, quote = F, sep = "\t", na = NA)


# Regress out covariates - Females
pheno_cor_F <- data.frame(IID = pheno_F$IID)
for (p in 2:ncol(pheno_F)) {
	# Define phenotype
	pheno <- colnames(pheno_F)[p]
	print(paste0("Covariates - Females (no INT): ", pheno))
	# Create a temporary dataframe
	temp <- right_join(cov_F, pheno_F[, c(1,p)], by = "IID")
	colnames(temp)[ncol(temp)] <- "pheno"
	# Regress out covariates
	pheno_cor_F[, pheno] <- residuals(lm(pheno ~ . , data = temp[, -c(1)], na.action = na.exclude))
}; rm (p, pheno, temp)

fwrite(pheno_cor_F, "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/05_phenotypes/continuous/white_british/data/final/age_age2_batch_PCs/pheno_continuous_WB_age_age2_batch_PCs_F.txt", col.names = T, row.names = F, quote = F, sep = "\t", na = NA)
