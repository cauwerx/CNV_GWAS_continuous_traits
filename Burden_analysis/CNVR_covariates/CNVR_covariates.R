# Identify samples with a CNV overlapping a CNVR identified through CNV-GWAS
# This information will be used to correct the burden analysis for modifier CNVRs

#################################################
### Libraries ###################################
library(data.table)
library(dplyr)


#################################################
### STEP 1: Load CNVR & remodel #################

# Load all CNVRs for the mirror model
mirror <- as.data.table(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/09_CNVR/continuous/white_british/mirror/data/final/GWAS_CNVR_mirror.txt", header = T, select = c(2,9,10), col.names = c("CHR", "START", "END"))) # "PHENO" is col7
print(paste0("Number of loaded mirror CNVR: ", nrow(mirror)))
mirror <- mirror[which(mirror$CHR != "X"), ] 							# Retain autosomal CNVRs
mirror$CHR <- as.numeric(mirror$CHR)
mirror <- mirror[!duplicated(mirror), ] 								# Remove duplicates
mirror$ID <- paste(mirror$CHR, mirror$START, mirror$END, sep = "_")		# Add unique ID
mirror <- mirror[order(mirror$CHR, mirror$START, mirror$END), ]
print(paste0("Number of unique autosomal mirror CNVRs: ", nrow(mirror)))

# Load all CNVRs for the duplication-only model
dup <- as.data.table(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/09_CNVR/continuous/white_british/duplication_only/data/final/GWAS_CNVR_dup.txt", header = T, select = c(2,9,10), col.names = c("CHR", "START", "END")))  # "PHENO" is col7
print(paste0("Number of loaded duplication CNVR: ", nrow(dup)))
dup <- dup[which(dup$CHR != "X"), ] 									# Retain autosomal CNVRs
dup$CHR <- as.numeric(dup$CHR)
dup <- dup[!duplicated(dup), ] 											# Remove duplicates
dup$ID <- paste(dup$CHR, dup$START, dup$END, sep = "_")					# Add unique ID
dup <- dup[order(dup$CHR, dup$START, dup$END), ]
print(paste0("Number of unique autosomal duplication CNVRs: ", nrow(dup)))

# Load all CNVRs for the deletion-only model
del <- as.data.table(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/09_CNVR/continuous/white_british/deletion_only/data/final/GWAS_CNVR_del.txt", header = T, select = c(2,9,10), col.names = c("CHR", "START", "END")))  # "PHENO" is col7
print(paste0("Number of loaded mirror CNVR: ", nrow(del)))
del <- del[which(del$CHR != "X"), ] 									# Retain autosomal CNVRs
del$CHR <- as.numeric(del$CHR)
del <- del[!duplicated(del), ] 											# Remove duplicates
del$ID <- paste(del$CHR, del$START, del$END, sep = "_")					# Add unique ID
del <- del[order(del$CHR, del$START, del$END), ]
print(paste0("Number of unique autosomal deletions CNVRs: ", nrow(del)))

# Merge all CNVRs across the 3 association models 
CNVR <- rbind(mirror, dup, del)
CNVR <- CNVR[!duplicated(CNVR), ]
print(paste0("Total number of CNVR for covariate file: ", nrow(CNVR)))


#################################################
### STEP 2: Load CNVs & remodel #################

# Load all caled CNVs; "ukb_cnv_global.gz" is the final PennCNV output, with all CNV calls in a linear format
cnvs <- as.data.table(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/cnv_calls/final/ukb_cnv_global.gz", header = T, select = c(1:5, 18), col.names = c("IID", "CHR", "START", "END", "CN", "QS")))
cnvs$IID <- sub("_.*", "", cnvs$IID)									# Clean the sample eid
cnvs <- cnvs[!cnvs$CHR %in% c(23, 24), ]								# Retain autosomal CNVs
cnvs <- cnvs[which(abs(cnvs$QS) >= 0.5), ] 								# Filter for quality score (retain high confidence CNVs)
print(paste0("Number of retained autosomal CNVs: ", nrow(cnvs)))


#################################################
### STEP 3: Calculate overlaps ##################

#### Set Key ####################################
setkey(CNVR, CHR, START, END)

#### Overlap between called CNVs and CNVRs ######
overlap <- na.omit(as.data.frame(foverlaps(cnvs, CNVR, maxgap = 0L, minoverlap = 1L, type = "any", mult = "all")))
print(paste0("Number of overlaps: ", nrow(overlap)))


#################################################
### STEP 4: Create covariate files ##############

# Load unrelated white British samples; # "samples_white_british_All.txt" are the results from GWAS/01_samples/sample_filtering.R 
eids <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/01_samples/white_british/data/final/samples_white_british_All.txt", header = F, col.names = "IID")) 

# Encode the covariate table as: -1 (deletion), 0 (copy-neutral), 1 (duplication)
sum_CNVR <- data.frame(CNVR = paste0("chr", CNVR$ID))
cov_cnvr <- data.frame(IID = eids$IID)
cov_cnvr[, paste0("chr", CNVR$ID)] <- 0
print(paste0("Dimensions of CNVR covariate file: ", nrow(cov_cnvr), " x ", ncol(cov_cnvr)))
for(cnvr in CNVR$ID) {
	# Detect CNV carriers	
	dup_carrier <- as.numeric(overlap[which(overlap$ID == cnvr & overlap$CN > 2), "IID"])
	del_carrier <- as.numeric(overlap[which(overlap$ID == cnvr & overlap$CN < 2), "IID"])
	# Add them to the covariate table 	
	cov_cnvr[cov_cnvr$IID %in% dup_carrier, paste0("chr",cnvr)] <- 1
	cov_cnvr[cov_cnvr$IID %in% del_carrier, paste0("chr",cnvr)] <- -1
	# Fill summary table
	sum_CNVR[which(sum_CNVR$CNVR == paste0("chr",cnvr)), "DUP"] <- length(cov_cnvr[which(cov_cnvr[, paste0("chr",cnvr)] == 1), paste0("chr",cnvr)])
	sum_CNVR[which(sum_CNVR$CNVR == paste0("chr",cnvr)), "DEL"] <- length(cov_cnvr[which(cov_cnvr[, paste0("chr",cnvr)] == -1), paste0("chr",cnvr)])
	sum_CNVR[which(sum_CNVR$CNVR == paste0("chr",cnvr)), "CNV"] <- length(cov_cnvr[which(cov_cnvr[, paste0("chr",cnvr)] != 0), paste0("chr",cnvr)])
}


#################################################
### Save ########################################

# Covariate matrix
fwrite(cov_cnvr, "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/04_covariates/data/final/covariates_CNVR.txt.gz", col.names = T, row.names = F, quote = F, sep = "\t")

# Summary table
fwrite(sum_CNVR, "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/04_covariates/data/final/summary/CNVR_individual_count.txt", col.names = T, row.names = F, quote = F, sep = "\t")
