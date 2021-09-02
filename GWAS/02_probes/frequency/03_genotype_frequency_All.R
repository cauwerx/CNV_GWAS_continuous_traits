# Calculate the CNV, duplication, and deletion frequency of all probes with genotype missingness < 5%

###############################################
### Libraries #################################
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)


#################################################
### STEP 1: Load genotype counts ################

# List all genotype count files from GWAS/02/frequency/01_genotype_count.sh
geno_count_list <- list.files("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/02_probes/white_british/stringent/frequency/data/temp/genotype_count/All", pattern = "frqx.gz", full.names = T, recursive = F)
print(paste0("Number of files to be merged: ", length(geno_count_list)))

# Merge files from different chromosomes
geno_count <- data.frame()
for (f in geno_count_list) {	
	df <- fread(f, header = T)
	geno_count <- rbind(geno_count, df)}
print(paste0("Number of probes assessed: ", nrow(geno_count)))


#################################################
### STEP 2: High-missingness probes #############

# Load high-missingness probes from GWAS/02/frequency/02_high_missingness_0.95.R
miss_probes <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/02_probes/white_british/stringent/frequency/data/final/high_missingness/genotype_missingness_0.95.txt", header = T, select = 2))
print(paste0("Number of probes excluded due to high missingness: ", nrow(miss_probes)))

# Exclude high (5%) missingness probes
geno_count <- geno_count[!geno_count$SNP %in% miss_probes$SNP, ]
print(paste0("Number of probes remaining after excluding high missingness probes: ", nrow(geno_count)))


#####################################################################
### STEP 3: Convert genotype counts to number of CNV carriers #######

# Set up a dataframe
geno_freq <- data.frame(Chr = geno_count$CHR, SNP = geno_count$SNP, NumCNV = NA, NumDup = NA, NumDel = NA, NumCN = NA, NumMiss = NA)

# Number of duplication carriers (NumDup; > 2 copies)
geno_freq[which(geno_count$A1 == "A"), "NumDup"] <- geno_count[which(geno_count$A1 == "A"), "C(HOM A2)"] 
geno_freq[which(geno_count$A1 == "T"), "NumDup"] <- geno_count[which(geno_count$A1 == "T"), "C(HOM A1)"] 

# Number of deletion carriers (NumDel; < 2 copies)
geno_freq[which(geno_count$A1 == "A"), "NumDel"] <- geno_count[which(geno_count$A1 == "A"), "C(HOM A1)"] 
geno_freq[which(geno_count$A1 == "T"), "NumDel"] <- geno_count[which(geno_count$A1 == "T"), "C(HOM A2)"]

# Number of CNV carriers (NumCNV; != 2 copies)
geno_freq$NumCNV <- geno_freq$NumDup + geno_freq$NumDel

# Number of copy-neutral carriers (NumCN, = 2 copies)
geno_freq[, "NumCN"] <- geno_count[, "C(HET)"]

# Number of missing genotypes (NumMiss)
geno_freq[, "NumMiss"] <- geno_count[, "C(MISSING)"]


#################################################
### STEP 4: Genotype frequencies [%] ############

geno_freq$FreqCNV <- 100 * (geno_freq$NumCNV/(geno_freq$NumCNV + geno_freq$NumCN))
geno_freq$FreqDup <- 100 * (geno_freq$NumDup/(geno_freq$NumCNV + geno_freq$NumCN))
geno_freq$FreqDel <- 100 * (geno_freq$NumDel/(geno_freq$NumCNV + geno_freq$NumCN))


#################################################
### STEP 5: Correct sex chromosome labels #######

geno_freq[which(geno_freq$Chr == 23), "Chr"] <- "X"
geno_freq[which(geno_freq$Chr == 25), "Chr"] <- "XY"


#################################################
### STEP 6: Print summary #######################

# CNVs
print("Number of probes with CNV frequency >= 0.005%:") 
print(as.data.frame(geno_freq[which(geno_freq$FreqCNV >= 0.005), ] %>% count(Chr)))
print(paste0("TOTAL: ", nrow(geno_freq[which(geno_freq$FreqCNV >= 0.005), ])))

# Duplications
print("Number of probes with duplication frequency >= 0.005%:") 
print(as.data.frame(geno_freq[which(geno_freq$FreqDup >= 0.005), ] %>% count(Chr)))
print(paste0("TOTAL: ", nrow(geno_freq[which(geno_freq$FreqDup >= 0.005), ])))

# Deletions
print("Number of probes with deletion frequency >= 0.005%:") 
print(as.data.frame(geno_freq[which(geno_freq$FreqDel >= 0.005), ] %>% count(Chr)))
print(paste0("TOTAL: ", nrow(geno_freq[which(geno_freq$FreqDel >= 0.005), ])))


#################################################
### Save ########################################

fwrite(geno_freq, "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/02_probes/white_british/stringent/frequency/data/final/genotype_frequency/genotype_freq_WB_All.txt.gz", col.names = T, row.names = F, quote = F, sep = "\t" )
