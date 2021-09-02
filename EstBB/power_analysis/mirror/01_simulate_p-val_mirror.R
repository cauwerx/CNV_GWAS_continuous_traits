# Simulate 10'000 p-values for each trait-CNVR association
# Use effect size from the UKBB 
# Use CNV frequency and sample sizes from teh EStBB

#################################################
### Libraries ###################################
library(data.table)
library(dplyr)
library(purrr)


#################################################
### External arguments ##########################

args <- commandArgs(trailingOnly = TRUE)
signal <- args[1]
print(paste0("Simulating mirror signal ", signal))


#################################################
### STEP 1: Load files ##########################

# Load replicated mirror trait-CNVR association characteristics; "GWAS_CNVR_MIR.withEstResults.withCounts.txt" contains for each signal, in the following order:
# probe ID, chromosome, position in the UKBB
# effect size, SE, and p-value in the UKBB (see GWAS/05_gwas/continuous/mirror)
# phenotype name, phenotype field ID
# start position, end position, length, and number of probes of the associated CNVR as defined in the UKBB (see GWAS/07_CNVR/continuous/mirror)
# probe ID, chromosome, position of matching signal in the EstBB (provided by our collaborrators at University of Tartu)
# sample size for the matching phenotype in the EstBB (provided by our collaborrators at University of Tartu)
# effect size, SE, and p-value in the UKBB (provided by our collaborrators at University of Tartu)
# number of probes covering the CNVR in the EstBB (provided by our collaborrators at University of Tartu)
# number of deletion carriers, duplication carriers, and copy-neutral individuals at that probe in the EstBB (provided by our collaborrators at University of Tartu)
associations <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/13_EstBB/raw/GWAS_CNVR_MIR.withEstResults.withCounts.txt"))

# Load CNV frequency data; "CNV_probe_full_frequency_summary.RDS" contains for each genotyped probe in the EstBB: chromosome, probe ID, number of CNV, deletion, duplication, copy-neutral, genotype missing individuals, as well as the CNV, deletion, and duplication frequencies 
frequencies <- readRDS("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/13_EstBB/raw/CNV_probe_full_frequency_summary.RDS")

# Load phenotype data; "continuous_traits_clean.txt" contains the phenotype name, number of individuals in the EstBB for which this phenotype is available, as well as the median, mean, and SD of the measurements for this trait
phenotypes <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/13_EstBB/raw/continuous_traits_clean.txt"))


#################################################
### STEP 2: Extract information #################

# Information on the trait-CNVR association to be replicated
pheno <- associations[signal, "PHENO"]; print(paste0("Phenotype: ", pheno))
chr <- associations[signal, "CHR"]; print(paste0("Chromosome: ", chr))
pos_uk <- associations[signal, "POS"]; print(paste0("Position (UKB): ", pos_uk))
pos_est <- associations[signal, "POS_Est"]; print(paste0("Position (EstBB): ", pos_est))
rs_uk <- associations[signal, "ID"]; print(paste0("Rs (UKB): ", rs_uk))
rs_est <- associations[signal, "MARKER_Est"]; print(paste0("Rs (EstBB): ", rs_est))
b_uk <- associations[signal, "BETA"]; print(paste0("Beta (UKB): ", b_uk))

# Information on the CNV frequency in the EstBB of the assessed region
freq_cnv_est <- frequencies[which(frequencies$probe == rs_est), "FreqCNV"]; print(paste0("CNV frequency (Est BB): ", freq_cnv_est))
freq_dup_est <- frequencies[which(frequencies$probe == rs_est), "FreqDUP"]; print(paste0("DUP frequency (Est BB): ", freq_dup_est))
freq_del_est <- frequencies[which(frequencies$probe == rs_est), "FreqDEL"]; print(paste0("DEL frequency (Est BB): ", freq_del_est))

# Information on the sample size in the EstBB of the assessed trait
N_est <- phenotypes[which(phenotypes$Pheno == pheno), "N_Est"]; print(paste0("N (EstBB): ", N_est))


#################################################
### STEP 3: Simulation ##########################

# Number of simulations
n_sim <- 10000

# Skip the analysis if there is no available data
if (is.na(N_est) | is.na(rs_est)) {print(paste0("No data available for ", pheno, " - chr", chr, ":", pos_uk))
	
	# If data is available
	} else {
		print(paste0("Performing ", n_sim, " simulations"))

		# Create empty vector to store generated p-value
		p_val <- vector()

		# Simulation
		for (i in 1:n_sim) {

			# Create a temporary empty matrix
			df <- as.data.frame(matrix(nrow = N_est, ncol = 3, dimnames = list(c(), c("CNV", "error", "Y"))))
		
			# Simulate CNVs: CNV ~ Uniform(n, min = 0, max = 1) --> attribute -1/0/1
			df[, "CNV"] <- runif(N_est, min = 0, max = 1)
			df[which(df$CNV >= (1-freq_dup_est)), "CNV"] <- 1
			df[which(df$CNV <= freq_del_est), "CNV"] <- -1
			df[which(df$CNV < (1-freq_dup_est) & df$CNV > freq_del_est), "CNV"] <- 0
				
			# Simulate the error: error ~ N(0, sig2), with sig^2 = SD^2 -var(CNV)*beta^2 
			sig2 <- 1 - var(df$CNV)*b_uk^2
			df[, "error"] <- rnorm(N_est, mean = 0, sd = sqrt(sig2))

			# Simulate phenotypes (Y): Y = CNV * beta + error
			df[, "Y"] <- df[, "CNV"] * b_uk + df[, "error"]

			# Perform a linear regression Y ~ CNV and retain the p-value
			df <- as.data.frame(df)
			if (nrow(df[which(df$CNV != 0), ]) > 0) {		
				p_val[i] <- summary(lm(Y ~ CNV, data = df))$coefficients[2,4]
			} else { 
				p_val[i] <- NA}

	}
	
	# Save generated p-values
	fwrite(data.frame(p_val = p_val), paste0("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/13_EstBB/power_analysis/mirror/data/temp/p_values/", pheno, ".", chr, ":", pos_uk, ".mirror.txt"), col.names = F, row.names = F, quote = F, sep = "\t")

}
