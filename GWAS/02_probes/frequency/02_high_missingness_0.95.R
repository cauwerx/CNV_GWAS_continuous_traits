# Probes with high genotype missingness (>5%) will be excluded from the GWAS
# Note that a threshold of 1%, 5%, and 10% gives similar results.

###############################################
### Libraries #################################
library(data.table)
library(dplyr)


#################################################
### STEP 1: Load & merge data ###################

# List all genotype count files from GWAS/02/frequency/01_genotype_count.sh
geno_count_list <- list.files("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/02_probes/white_british/stringent/frequency/data/temp/genotype_count/All", pattern = "frqx.gz", full.names = T, recursive = F)

# Merge files from different chromosomes
geno_count <- data.frame()
for (f in geno_count_list) {	
	df <- fread(f, header = T)
	geno_count <- rbind(geno_count, df)}
print(paste0("Number of probes assessed: ", nrow(geno_count)))


#################################################
### STEP 2: Detect high missingness probes ######

# Calculate the corresponding number of individuals (unrelated, white British) for a missingness threshold at 5% 
Num <- nrow(as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/01_samples/white_british/data/final/samples_white_british_All.txt", header = F)))
missingness_thr <- Num*0.05

# Extract probes with high missingness
missingness_probes <- geno_count[which(geno_count$`C(MISSING)` > missingness_thr), ]
print(paste0("Number of probes with >5% missing genotypes in all WB individuals: ", nrow(missingness_probes)))
print(as.data.frame(missingness_probes %>% count(CHR)))


#################################################
### Save ########################################

fwrite(missingness_probes[,c(1,2)], "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/02_probes/white_british/stringent/frequency/data/final/high_missingness/genotype_missingness_0.95.txt", col.names = T, row.names = F, quote = F, sep = "\t")
