# Harmonize GWAS effects to A1 = "T"

#################################################
### Libraries ###################################
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)


#################################################
### External arguments ##########################

args <- commandArgs(trailingOnly = TRUE)
sex <- args[1]
print(paste0("Correcting files in ", sex))


#################################################
### STEP 1: List files to correct ###############

# These are the output files from the GWAS
list_files <- list.files(paste0("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/06_gwas/continuous/white_british/duplication_only/data/temp/", sex,"/associations"), pattern = ".glm.linear", full.names = T, recursive = F)
print(paste0("Number of files to be analyzed: ", length(list_files)))


#################################################
### STEP 2: Loop over files & correct ###########

# Loop over GWAS files
for (f in list_files) {

	# Load the file
	df <- as.data.frame(fread(f, header = T)) 

	# Correct the BETA, CI, A1, REF, and ALT
	df[which(df$A1 == "A"), c(9, 11, 12, 13)] <- -1 * df[which(df$A1 == "A"), c(9, 11, 12, 13)]
	df[which(df$A1 == "A"), "REF"] <- "A"
	df[which(df$A1 == "A"), "ALT"] <- "T"
	df[which(df$A1 == "A"), "A1"] <- "T"
	colnames(df)[c(11,12)] <- c("U95", "L95")
	df <- df[, c(1:10, 12, 11, 13:15)]
	
	# Save
	fwrite(df, f, col.names = T, row.names = F, quote = F, sep = "\t", na = NA)

}
