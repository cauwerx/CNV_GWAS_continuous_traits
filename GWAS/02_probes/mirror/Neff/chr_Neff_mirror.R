# Calculate the chromosome-wide number of effective test, based on the method of Gao et al., 2008 (DOI: 10.1002/gepi.20310) 

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
print(paste0("Analyzing chromosome ", chr))


#################################################
### STEP 1: Clean the .ped file #################

### Load the core files #########################

# All samples; This corresponds to the .ped file matching the binary PLINK file set described in the methods of the paper (Table 1)
all_samples <- as.data.frame(fread(paste0("/data/FAC/FBM/DBC/zkutalik/default_sensitive/cauwerx/CNV_GWA/mirror/ukb_cnv_ped/ukb_cnv_chr", chr,".ped"), header = F, select = 2, col.names = "eid"))
print(paste0("There are ", nrow(all_samples), " samples in the .ped file"))

# Samples to keep; "samples_white_british_All.txt" are the results from GWAS/01_samples/sample_filtering.R
eids_to_keep <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/01_samples/white_british/data/final/samples_white_british_All.txt", header = F, col.names = "eid"))
print(paste0("There are ", nrow(eids_to_keep), " samples to keep"))

# All probes; This corresponds to the .bim of the PLINK file set described in the methods of the paper (Table 1)
all_probes <- as.data.frame(fread(paste0("/data/FAC/FBM/DBC/zkutalik/default_sensitive/cauwerx/CNV_GWA/mirror/ukb_cnv_bed/ukb_cnv_chr", chr,".bim"), header = F, select = c(1,2), col.names = c("chr","rs"))) 
print(paste0("There are ", nrow(all_probes), " probes in the .bim file"))

# Probes to keep; Results from GWAS/02_probes/mirror/probe_pruning_mirror.R at r^2 threshold of 0.9999
probes_to_keep <- as.data.frame(fread(paste0("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/02_probes/white_british/stringent/mirror/data/final/pruning/r_0.9999/probes_WB_All_prune_0.9999_mirror_chr", chr,".prune.in"), header = F, col.names = "rs"))
print(paste0("There are ", nrow(probes_to_keep), " probes to keep"))


#### Split .ped in chuncks of 20'000 columns ##

# This corresponds to the .ped file matching the binary PLINK file set described in the methods of the paper (Table 1)
ped_header <- length(fread(paste0("/data/FAC/FBM/DBC/zkutalik/default_sensitive/cauwerx/CNV_GWA/mirror/ukb_cnv_ped/ukb_cnv_chr", chr, ".ped"), nrow = 1))
ped_split <- split(7:ped_header, ceiling(seq_along(7:ped_header)/20000))


### Loop to load .ped by chuncks & clean it #######
cnv <- matrix()
counter <- 1

for (i in 1:length(ped_split)) {

	# Parameters
	start <- min(ped_split[[i]])
	end <- max(ped_split[[i]])
	NumProbes <- length(ped_split[[i]])/2

	# Load .ped chunk
	ped_chunck <- as.data.frame(fread(paste0("/data/FAC/FBM/DBC/zkutalik/default_sensitive/cauwerx/CNV_GWA/mirror/ukb_cnv_ped/ukb_cnv_chr", chr,".ped"), header = F, select = start:end))
	
	# Clean samples by retaining unrelated samples of white Brritish ancestry (rows)
	rownames(ped_chunck) <- all_samples$eid
	ped_chunck <- ped_chunck[rownames(ped_chunck) %in% eids_to_keep$eid, ]

	# Clean variants by retaining filtered and prruned probes (columns)
	cnv_chunck <- data.frame(mapply(paste0, ped_chunck[][c(T, F)], ped_chunck[][c(F,T)]), stringsAsFactors = F)
	rownames(cnv_chunck) <- rownames(ped_chunck); rm (ped_chunck)
	colnames(cnv_chunck) <- all_probes[counter:(counter + NumProbes - 1), "rs"] 
	cnv_chunck <- cnv_chunck[, colnames(cnv_chunck) %in% probes_to_keep$rs]

	# Merge
	if(i == 1) {cnv <- cnv_chunck} else {cnv <- cbind(cnv, cnv_chunck)}

	# Update counter
	counter <- counter + NumProbes
}

print(paste0("Final dimensions for the cleaned cnv file: ", nrow(cnv), " (samples; row) x ", ncol(cnv), " (probes; column)"))


### Recode & make numeric matrix ################

# AA -> -1 (= deletion); AT -> 0 (= copy-neutral); TT -> 1 (= duplication)
cnv[cnv == "AA"] <- -1
cnv[cnv == "AT" | cnv == "TA"] <- 0
cnv[cnv == "TT"] <- 1
cnv[cnv == "00"] <- NA
cnv <- mutate_all(cnv, function(x) as.numeric(as.character(x)))

# Tranform in matrix
cnv <- as.matrix(cnv)


#################################################
### STEP 2: Calculate the correlation matrix ###

# Calculate the correlation matrix
cnv_cor <- cor(cnv, use = 'pairwise.complete.obs')
rm(cnv)

# Set missing values to 0 to allow svd()
cnv_cor[which(is.na(cnv_cor))] <- 0


#################################################
### STEP 3: Calculate the eigenvalues ###########

# svd() computes the singular-value decomposition of a recatngular matrix, with $d being a vector containing the singular values of the decomposition
cnv_EV <- svd(cnv_cor)$d 


#################################################
### STEP 4: Calculate Neff ######################

# Neff is defined as the number of eigenvalues required to explain 99.5% of the variation from the CNV data 
sum_EV <- 0
count <- 0

while(sum_EV/sum(cnv_EV) < 0.995) {
	count <- count + 1
	sum_EV <- sum_EV + cnv_EV[count]
}

print(paste0("Neff on chromosome ", chr, ": ", count))


#################################################
### Save ########################################
fwrite(data.frame(count), paste0("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/02_probes/white_british/stringent/mirror/data/final/Neff/Neff_chr", chr, ".txt"), col.names = F, row.names = F, quote = F, sep = "\t")

