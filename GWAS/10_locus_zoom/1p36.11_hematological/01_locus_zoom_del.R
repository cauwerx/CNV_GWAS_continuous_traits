# Locus zoom between the 1p36.11 locus and associated hematological traits

#################################################
### Libraries ###################################
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)


#################################################
### External arguments ##########################

args <- commandArgs(trailingOnly = TRUE)
chr <- args[1]
start <- args[2]
end <- args[3]
pheno <- gsub("-", " ", args[4])
print(paste0("Locus zoom plot for: chr", chr, ":", start, "-", end, " - ", pheno))


#################################################
### STEP 1: Temporary missingness file ##########

# Here we will perform a GWAS with all probe in the region of interest, only excluding high missingness probes
hm <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/02_probes/white_british/stringent/frequency/data/final/high_missingness/genotype_missingness_0.95.txt"))
fwrite(data.frame(hm$SNP), "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/99_candidates/1p36.11_hematological/data/temp/locus_zoom/high_missingness.txt", col.names = F, row.names = F, quote = F, sep = "\t")


#################################################
### STEP 2: Prepare PLINK input #################

# PLINK file set - deletion-only model
plink_file <- paste0("/data/FAC/FBM/DBC/zkutalik/default_sensitive/cauwerx/CNV_GWA/deletion_only/ukb_cnv_bed/ukb_cnv_chr", chr)

# INT covariate-corrected phenotypes (from GWAS/04_phenotypes/continuous/phenotype_extraction.R)
pheno_file <- "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/05_phenotypes/continuous/white_british/data/final/INT_age_age2_batch_PCs/pheno_continuous_WB_INT_age_age2_sex_batch_PCs_All.txt"

# Unrelated, white British samples (from GWAS/01_samples/sample_filtering.R)
samples_file <- "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/01_samples/white_british/data/final/samples_white_british_All.txt"

# Probes to exclude due to high missingness (see STEP 1)
HM_file <- "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/99_candidates/1p36.11_hematological/data/temp/locus_zoom/high_missingness.txt"

# Output file
output_file <- paste0("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/99_candidates/1p36.11_hematological/data/temp/locus_zoom/del_All_chr", chr, ":", start, ":", end)


#################################################
### STEP 3: Perform local association study #####

# Run PLINK v2
system(paste("/home/cauwerx/scratch/cauwerx/softwares/plink2/plink2",
 			 "--bfile", plink_file,
			 "--pheno",  pheno_file, "--no-psam-pheno", "--pheno-name", pheno,
			 "--glm omit-ref no-x-sex hide-covar allow-no-covars --ci 0.95", 
			 "--keep", samples_file,
			 "--exclude",  HM_file, "--chr", chr, "--from-bp", start, "--to-bp", end,
			 "--out", output_file))
unlink(HM_file)
unlink(paste0(output_file, ".log"))


#################################################
### STEP 4: Correct association direction #######

# Loop over phenotypes tested previously
for(p in unlist(strsplit(pheno, " "))) {

	# Load association file
	df <- as.data.frame(fread(paste0(output_file, ".", p, ".glm.linear"), header = T))

	# Correct the BETA, CI, A1, REF, and ALT
	df[which(df$A1 == "A"), c(9, 11, 12, 13)] <- -1 * df[which(df$A1 == "A"), c(9, 11, 12, 13)]
	df[which(df$A1 == "A"), "REF"] <- "A"
	df[which(df$A1 == "A"), "ALT"] <- "T"
	df[which(df$A1 == "A"), "A1"] <- "T"
	colnames(df)[c(11,12)] <- c("U95", "L95")
	df <- df[, c(1:10, 12, 11, 13:15)]
	
	# Save corrected version and delete the old file
	fwrite(df, paste0("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/99_candidates/1p36.11_hematological/data/final/locus_zoom/del_All_chr", chr, ":", start, "-", end, ".", p, ".glm.linear"), col.names = T, row.names = F, quote = F, sep = "\t", na = NA)
	unlink(paste0(output_file, ".", p, ".glm.linear"))
}

unlink(paste0(output_file, "/*"), expand = T)




