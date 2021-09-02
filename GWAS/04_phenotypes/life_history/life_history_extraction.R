# Extract life historyphenotype data and do the inverse normal transformation and covariate correction

################################################
### Libraries #################################
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(matrixStats)


#################################################
### STEP 1: Load Main Data ######################

# This corresponds to a .csv file with two columns: field ID and short name for the corresponding phenotype
pheno_list <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/05_phenotypes/life_history/white_british/data/raw/LH_traits.csv", header = T))


#################################################
### STEP 2: Order raw phenotype files ###########

# As phenotypic data is spread throughout multiple file downloaded from the UKBB portal, we start by listing all phenotype files available and order them
raw_pheno_file <- file.info(list.files("/data/FAC/FBM/DBC/zkutalik/default_sensitive/UKBB/pheno", pattern = "^ukb.*csv", full.names = T, recursive = F))
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
fwrite(pheno_list, "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/05_phenotypes/life_history/white_british/data/final/pheno_list.txt", col.names = T, row.names = F, quote = F, sep = "\t")


#################################################
### STEP 4: Extract phenotypes and average ######

# This file contains all UKBB sample eids, except for the redacted ones
phenotypes <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/general_data/UKBB/eid_no_redacted3.txt.gz", header = T, select = c(2)), drop = F)

# Loop over all phenotypes to extract the data
for(p in 1:nrow(pheno_list)) {
	
	# Define the phenotype
	pheno <- pheno_list[p, "Pheno"] 
	ID <- pheno_list[p, "FieldID"] 
	file <- pheno_list[p, "File"] 
	col <- pheno_list[p, "Col"] 
	print(paste0("Extracting ", pheno, " from ", sub(".*/", "", file)))

	# Extract columns corresponding to the defined phenotype
	temp_pheno <- as.data.frame(fread(file, header = T, select = c(1, as.numeric(str_split(col, pattern = "_")[[1]]))))

	# IF: phenotype is a parental lifespan or income --> Set negatives to NA; mean
	if(pheno == "lifespan_mother" | pheno == "lifespan_father" | pheno == "income") {
	temp_pheno[temp_pheno < 0] <- NA
	temp_pheno <- data.frame(eid = temp_pheno[, 1], pheno = rowMeans(temp_pheno[, 2:ncol(temp_pheno)], na.rm = T))	

	# IF: phenotype is EA --> Set negative values -3 (prefer not to answer) and -1 (do not know) to NA and set -2 (never went to school) to 5 (min value);mean 
	} else if (pheno == "EA") {
	temp_pheno[temp_pheno == -1 | temp_pheno == -3] <- NA
	temp_pheno[temp_pheno == -2] <- 5
	temp_pheno <- data.frame(eid = temp_pheno[, 1], pheno = rowMeans(temp_pheno[, 2:ncol(temp_pheno)], na.rm = T)) 

	# IF: none of the above --> Average, if more than one column present
	} else if (ncol(temp_pheno) > 2) {
	print("Averaging over instances")
	temp_pheno <- data.frame(eid = temp_pheno[, 1], pheno = rowMeans(temp_pheno[, 2:ncol(temp_pheno)], na.rm = T))}

	# Merge
	colnames(temp_pheno) <- c("eid", pheno)
	phenotypes <- full_join(phenotypes, temp_pheno, by = "eid")

}
rm(pheno, ID, file, col, temp_pheno)
print(paste0("Dimensions of phenotype table: ", ncol(phenotypes), " pheno x ", nrow(phenotypes), " eids"))


#################################################
### STEP 5: Replace entries by NA ###############

phenotypes[phenotypes < 0] <- NA


#################################################
### STEP 6: Select samples ######################

# Retain unrelated, white British samples selected during GWAS/01_samples/sample_filtering.R

# All retained individuals 
eid_all <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/01_samples/white_british/data/final/samples_white_british_All.txt", header = F, col.names = "eid"))
pheno_all <- phenotypes[phenotypes$eid %in% eid_all$eid, ]
colnames(pheno_all)[1] <- "IID"
print(paste0("Dimensions of phenotype table for all unrelated white British individuals: ", ncol(pheno_all)-1, " pheno x ", nrow(pheno_all), " eids"))
fwrite(pheno_all, "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/05_phenotypes/life_history/white_british/data/final/raw/life_history_WB_raw_All.txt", col.names = T, row.names = F, quote = F, sep = "\t", na = NA)


#################################################
### STEP 7:Summary ##############################

# All retained individuals 
summary_pheno_all <- data.frame(Pheno = colnames(pheno_all)[2:ncol(pheno_all)],
								Num = colSums(!is.na(pheno_all[, c(2:ncol(pheno_all))])),
								NumNA = colSums(is.na(pheno_all[, c(2:ncol(pheno_all))])),
								Mean = colMeans(pheno_all[, c(2:ncol(pheno_all))], na.rm = T),
								Median = colMedians(as.matrix(pheno_all[, c(2:ncol(pheno_all))]), na.rm = T),
								SD = colSds(as.matrix(pheno_all[, c(2:ncol(pheno_all))]), na.rm = T),
								Min = colMins(as.matrix(pheno_all[, c(2:ncol(pheno_all))]), na.rm = T),
								Max = colMaxs(as.matrix(pheno_all[, c(2:ncol(pheno_all))]), na.rm = T))
fwrite(summary_pheno_all, "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/05_phenotypes/life_history/white_british/data/final/raw/summary_life_history_WB_raw_All.txt", col.names = T, row.names = F, quote = F, sep = "\t", na = NA)
print("Summary all: done!")


#################################################
### STEP 8: Inverse Normal Transformation #######

# All retained individuals 
pheno_int_all <- data.frame(IID = pheno_all$IID)
for (p in 2:ncol(pheno_all)) {
	pheno <- colnames(pheno_all)[p]
	print(paste0("INT - All: ", pheno))
	pheno_int_all[, pheno] <- qnorm((rank(pheno_all[, p], na.last = "keep")-0.5)/sum(!is.na(pheno_all[, p])))
}; rm (p, pheno)

fwrite(pheno_int_all, "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/05_phenotypes/life_history/white_british/data/final/INT/life_history_WB_INT_All.txt", col.names = T, row.names = F, quote = F, sep = "\t", na = NA)


#################################################
### STEP 9: Correct for covariates ##############

# Read in covariates from GWAS/covariates/covariate_extraction.R
cov <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/04_covariates/data/final/covariates.txt", header = T, select = c(1:4, 6:46)))
colnames(cov)[1] <- "IID"

# Select subsets (retained samples)
cov_all <- cov[cov$IID %in% eid_all$eid, ] # 331523 eids
print(paste0("Covariates all: ", ncol(cov_all)-1))
print(colnames(cov_all[-1]))

# Regress out covariates - All retained individuals
pheno_int_cor_all <- data.frame(IID = pheno_int_all$IID)
for (p in 2:ncol(pheno_int_all)) {
	
	# Define phenotype
	pheno <- colnames(pheno_int_all)[p]
	print(paste0("Covariates - All: ", pheno))
	
	# Do not correct for age and age^2 when the phenotype is age
	if (pheno == "age") {
		temp <- right_join(cov_all[, -c(2,3)], pheno_int_all[, c(1,p)], by = "IID")
		print(paste0("Number of covariates: ", ncol(temp)-2))
	} else {
		temp <- right_join(cov_all, pheno_int_all[, c(1,p)], by = "IID")
		print(paste0("Number of covariates: ", ncol(temp)-2))}

	# Regress out covariates
	colnames(temp)[ncol(temp)] <- "pheno"
	pheno_int_cor_all[, pheno] <- residuals(lm(pheno ~ . , data = temp[, -c(1)], na.action = na.exclude))
}; rm (p, pheno, temp)

fwrite(pheno_int_cor_all, "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/05_phenotypes/life_history/white_british/data/final/INT_age_age2_batch_PCs/life_history_WB_INT_age_age2_sex_batch_PCs_All.txt", col.names = T, row.names = F, quote = F, sep = "\t", na = NA)
