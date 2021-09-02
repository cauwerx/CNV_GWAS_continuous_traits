# Merge files by phenotype and extract significant associations 

#################################################
### Libraries ###################################
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(gtools)


#################################################
### External arguments ##########################

args <- commandArgs(trailingOnly = TRUE)
sex <- args[1]
print(paste0("Correcting files in ", sex))


# ************************************************ GWAS RESULTS ************************************************ #
print("Merging GWAS data")
# Define threshold for significance (FWER 5%)
thr <- 0.05/11804


##########################################
### STEP 1: List phenotypes ##############

# These are the corrected output files from the GWAS
list_files <- list.files(paste0("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/06_gwas/life_history/white_british/mirror/data/temp/", sex,"/associations"), pattern = ".glm.linear", full.names = T, recursive = F)
pheno_list <- unique(gsub(".*\\.", "", sub(".glm.linear", "", list_files)))
print(paste0("There are ", length(pheno_list), " unique phenotypes in ", sex))
rm(list_files)


##########################################
### STEP 2: Merge and extract ############

# Create an empty dataframe
sig_all <- data.frame()

# Loop over phenotypes
for (p in pheno_list) {

	# Define the phenotype 
	print(paste0("Starting analyzing ", p))

	# List files, merge, and save
	merged_gwas <- data.frame()
	files <- mixedsort(list.files(paste0("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/06_gwas/life_history/white_british/mirror/data/temp/", sex,"/associations"), pattern = paste0("\\.", p, ".glm.linear"), full.names = T, recursive = F))
	print(paste0("Number of files to merge: ", length(files)))
	for (f in files) {
		df <- as.data.frame(fread(f, header = T, select = c(1:6, 8:15)), stringsAsFactors = F)
		colnames(df)[1] <- "CHR"
		merged_gwas <- rbind(merged_gwas, df)}
	fwrite(merged_gwas, paste0("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/06_gwas/life_history/white_british/mirror/data/final/", sex,"/gwas_", p, "_", sex, "_mirror.txt.gz"), col.names = T, row.names = T, quote = F, sep = "\t")

	# Select and save significant hits
	sig <- merged_gwas[which(merged_gwas$P <= thr), ]
	if (nrow(sig) >= 1) {
		sig$PHENO <- p
		sig_all <- rbind(sig_all, sig)
		fwrite(sig, paste0("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/06_gwas/life_history/white_british/mirror/data/final/", sex,"/significant_associations/significant_", p, "_", sex, "_mirror.txt.gz"), col.names = T, row.names = T, quote = F, sep = "\t")
		print(paste0(nrow(sig), " significant associations for ", p))
	} else {
		print("No significant association")}

}

print(paste0("Total number of significant associations in ", sex,": ", nrow(sig_all)))
fwrite(sig_all, paste0("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/06_gwas/life_history/white_british/mirror/data/final/", sex,"/significant_associations/00_significant_", sex, "_mirror.txt.gz"), col.names = T, row.names = F, quote = F, sep = "\t")
rm(p, f, files, merged_gwas, sig, sig_all)
