# Perform a step-wise conditional analysis for the sex-agnostic continuous phenotypes

#################################################
### Libraries ###################################
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(gtools)


#################################################
### SET UP WHILE LOOP ###########################

# Set a parameter to count the number of iterations
round <- 1

# Load and count the number of significant hits from the GWAS
n <- nrow(as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/06_gwas/continuous/white_british/duplication_only/data/final/All/significant_associations/00_significant_All_dup.txt.gz", header = T, select = c(1:3, 8:9, 13:15))))

# Start a while loop
while (n > 0) {
	print(paste0("STARTING ROUND ", round))

#################################################
### STEP 1: Select top probes ###################

# Access the significant data; For round = 1 the GWAS output; Else the output of the previous SCA iteration
if (round == 1) {
		new_sig <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/06_gwas/continuous/white_british/duplication_only/data/final/All/significant_associations/00_significant_All_dup.txt.gz", header = T, select = c(1:3, 8:9, 13:15)))
	} else {
		new_sig <- as.data.frame(fread(paste0("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/08_SCA/continuous/white_british/duplication_only/data/temp/significant_associations/All/sig_assoc_dup_SCA_round", round-1, ".txt.gz"), header = T, select = c(1:3, 8:9, 13:15)))}

# Extract top GW-significant variant/phenotype
new_sig <- new_sig[order(new_sig$P), ]
new_sig <- new_sig[!duplicated(new_sig[, "PHENO"]), ]
new_sig$ROUND <- round
print(paste0("Number of phenotypes with at least one significant association added at round ", round,": ", nrow(new_sig)))


#################################################
### STEP 2: Fill summary file ###################

# Create (round = 1) or load (round > 1) and file the summary file
if (round == 1) {
		sum_file <- new_sig
	} else {
		sum_file <- as.data.frame(fread(paste0("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/08_SCA/continuous/white_british/duplication_only/data/temp/summary/All/summary_dup_SCA_round", round-1, ".txt"), header = T))
		sum_file <- rbind(sum_file, new_sig)}

print(paste0("Total number of independent signals after round ",round,": ", nrow(sum_file)))
fwrite(sum_file, paste0("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/08_SCA/continuous/white_british/duplication_only/data/temp/summary/All/summary_dup_SCA_round", round, ".txt"), col.names = T, row.names = F, quote = F, sep = "\t", na = NA)


#################################################
### STEP 3: Probe CNV profile ###################

# Detect the probes whose CNV profile is to be extracted
probes <- sum_file[order(sum_file$ROUND), c(1,3,8,9)]
probes <- probes[!duplicated(probes[, "ID"]), ]
probes <- probes[which(probes$ROUND == round), ]
print(paste0("Number of probes to prepare: ", nrow(probes)))

# Loop over probes to generate their CNV profile
if(nrow(probes) > 0) {
for(i in 1:nrow(probes)) {

	# Define the probe
	rs <- as.character(probes[i, "ID"])
	chr <- as.character(probes[i, "CHR"])
	print(paste0("Starting probe: ", rs, " (chr", chr, ")"))

	# Write file
	fwrite(data.frame(rs), paste0("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/08_SCA/continuous/white_british/duplication_only/data/temp/cnv_profiles/All/", rs), col.names = F, row.names = F, quote = F, sep = "\t")

	# Extract profile with PLINK
	input <- paste0("/data/FAC/FBM/DBC/zkutalik/default_sensitive/cauwerx/CNV_GWA/duplication_only/ukb_cnv_bed/ukb_cnv_chr", chr) # PLINK file set for the duplication-only model
	output <- paste0("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/08_SCA/continuous/white_british/duplication_only/data/temp/cnv_profiles/All/", rs)
	system(paste("/home/cauwerx/scratch/cauwerx/softwares/plink/plink",
				 "--bfile", input, 
				 "--keep /home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/01_samples/white_british/data/final/samples_white_british_plink1.9_All.txt",
				 "--extract", output,
				 "--recode tab --out", output))

	# Correct the .ped in R
 	ped <- as.data.frame(fread(paste0(output, ".ped"), header = F, select = c(2,7), col.names = c("IID", "CNV")))
	ped[which(ped$CNV == "A A"), "CNV"] <- NA
	ped[which(ped$CNV == "A T" | ped$CNV == "T A"), "CNV"] <- 0
	ped[which(ped$CNV == "T T"), "CNV"] <- 1
	colnames(ped)[2] <- rs
	fwrite(ped, paste0(output, "_profile"), col.names = T, row.names = F, quote = F, sep = "\t")

	# Delete the rs file and the .ped/.map
	unlink(output); unlink(paste0(output, ".log")); unlink(paste0(output, ".map")); unlink(paste0(output, ".ped"))

}
}
rm(i, rs, chr, input, output, ped)


#################################################
### STEP 4: Prepare phenotypes ##################
print("Starting to regress out the covariates")

# Read in INT phenotypes for all and select those with significant association; From GWAS/04_phenotypes/continuous/phenotype_extraction.R
pheno <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/05_phenotypes/continuous/white_british/data/final/INT/pheno_continuous_WB_INT_All.txt", header = T))
pheno <- pheno[ , colnames(pheno) %in% c("IID", unique(new_sig$PHENO))]

# Read in covariates; From GWAS/03_covariates/covariate_extraction.R
cov <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/04_covariates/data/final/covariates.txt", header = T, select = c(1:4, 6:46)))
colnames(cov)[1] <- "IID"
cov <- cov[cov$IID %in% pheno$IID, ]

# Correct the phenotypes for covariates
pheno_cor <- data.frame(IID = pheno$IID)
for (i in 2:ncol(pheno)) {

	# Define phenotype and variants to correct for
	p <- colnames(pheno)[i]
	rs_list <- sum_file[which(sum_file$PHENO == p), ]

	# Add SCA variants to the covariate file
	temp <- cov
	for (rs in rs_list$ID) {
		rs_profile <- as.data.frame(fread(paste0("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/08_SCA/continuous/white_british/duplication_only/data/temp/cnv_profiles/All/", rs, "_profile"), header = T))
		temp <- left_join(temp, rs_profile, by = "IID")}

	# Add the phenotype of interest to the temp file
	temp <- right_join(temp, pheno[ , names(pheno) %in% c("IID", p)], by = "IID")
	colnames(temp)[ncol(temp)] <- "PHENO"

	# Regress out covariates
	pheno_cor[, p] <- residuals(lm(PHENO ~ . , data = temp[, -c(1)], na.action = na.exclude))
	
}

fwrite(pheno_cor, paste0("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/08_SCA/continuous/white_british/duplication_only/data/temp/phenotypes/All/pheno_cor_round", round, ".txt"), col.names = T, row.names = F, quote = F, sep = "\t", na = NA)
rm(i, p, rs_list, rs, rs_profile, temp, pheno_cor, new_sig, sum_file)


#################################################
### STEP 5: GWAS ################################

# Note: GWAS will be performed by chromosome, recursively
print("Strating GWAS")
system(paste0("mkdir /home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/08_SCA/continuous/white_british/duplication_only/data/temp/associations/All/round", round))

# List chromosomes
chromosomes <- c(1:22, "X", "XY")

# Loop over chromosomes
for(chr in chromosomes) {

	# Define variables
	input <- paste0("/data/FAC/FBM/DBC/zkutalik/default_sensitive/cauwerx/CNV_GWA/duplication_only/ukb_cnv_bed/ukb_cnv_chr", chr) # PLINK binary file set for duplication-only model
	pheno_file <- paste0("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/08_SCA/continuous/white_british/duplication_only/data/temp/phenotypes/All/pheno_cor_round", round, ".txt") # INT, covariate-corrected phenotype file generate above
	sample_file <- paste0("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/01_samples/white_british/data/final/samples_white_british_All.txt") # Unrelated white British individuals from GWAS/01_samples/sample_filtering.R
	probe_file <- paste0("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/02_probes/white_british/stringent/duplication_only/data/final/pruning/r_0.9999/probes_WB_All_prune_0.9999_dup_chr", chr,".prune.in") # Filtered and pruned (r^2 at 0.9999) from GWAS/02_probes/duplication_only/probe_pruning_dup.R
	output <- paste0("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/08_SCA/continuous/white_british/duplication_only/data/temp/associations/All/round", round,"/gwas_dup_SCA_round", round,"_chr", chr)

	# GWAS in PLINK v2
	system(paste("/home/cauwerx/scratch/cauwerx/softwares/plink2/plink2",  	
				  "--bfile", input, 
				  "--fam /data/FAC/FBM/DBC/zkutalik/default_sensitive/cauwerx/CNV_GWA/mirror/ukb_cnv_bed/ukb_cnv_chrX_onlyF.fam", 
				  "--pheno",  pheno_file, "--no-psam-pheno", 
				  "--glm omit-ref no-x-sex hide-covar allow-no-covars --ci 0.95",
				  "--keep", sample_file, 
				  "--extract", probe_file,
				  "--out", output))

}
rm(chromosomes, chr, input, pheno_file, sample_file, probe_file, output)


#################################################
### STEP 6: Correct GWAS & ID significants ######

# Threshold for significance (FWER 5%)
thr <- 0.05/11804

# List files to corrrect (harmonization of effect allele to "T")
list_gwas <- list.files(paste0("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/08_SCA/continuous/white_british/duplication_only/data/temp/associations/All/round", round), pattern = ".glm.linear", full.names = T, recursive = F)

# Create an empty dataframe
sig_gwas <- data.frame()

# Loop over files to correct the effect allele and extract significant signals
for (f in list_gwas) {

	# Load the file
	df <- as.data.frame(fread(f, header = T))
	p <- sub(".*\\.", "",sub(".glm.linear", "", f)) 

	# Correct the BETA, CI, A1, REF, and ALT
	df[which(df$A1 == "A"), c(9, 11, 12, 13)] <- -1 * df[which(df$A1 == "A"), c(9, 11, 12, 13)]
	df[which(df$A1 == "A"), "REF"] <- "A"
	df[which(df$A1 == "A"), "ALT"] <- "T"
	df[which(df$A1 == "A"), "A1"] <- "T"
	colnames(df)[c(11,12)] <- c("U95", "L95")
	df <- df[, c(1:10, 12, 11, 13:15)]
	
	# Save corrected version under the same name for follow-up analyses
	fwrite(df, f, col.names = T, row.names = F, quote = F, sep = "\t", na = NA)

	# Extract significant signals
	df <- df[, c(1:6, 8:15)]
	colnames(df)[1] <- "CHR"
	df$PHENO <- p
	temp_sig <- df[which(df$P <= thr), ]
	sig_gwas <- rbind(sig_gwas, temp_sig)

}

print(paste0("Number of significant associations after round ", round, ": ", nrow(sig_gwas)))
fwrite(sig_gwas, paste0("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/08_SCA/continuous/white_british/duplication_only/data/temp/significant_associations/All/sig_assoc_dup_SCA_round", round, ".txt.gz"), col.names = T, row.names = F, quote = F, sep = "\t", na = NA)
rm(f, df, temp_sig)


#################################################
### STEP 7: Update loop variables ###############

round <- round + 1
n <- nrow(sig_gwas)
rm(sig_gwas)

}
