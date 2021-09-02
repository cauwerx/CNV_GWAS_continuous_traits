# Retrieve samples used for the CNV-GWAS analysis 

################################################
### Libraries ##################################
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)


################################################
### Load Main Data #############################

# Sample QC; This corresponds to the "ukb_sqc_v2.txt" file available from the UKBB portal
sample_qc <- as.data.frame(fread("/data/FAC/FBM/DBC/zkutalik/default_sensitive/UKBB/geno/ukb_sqc_v2.txt", header = F, select = c(3:68)))
colnames(sample_qc) <- c("array", "batch", "plate", "well", "cluster_CR", "dQC", "dna_concentration", "submitted_gender", "inferred_gender", "X_intensity", "Y_intensity", "submitted_plate", "submitted_well", "missing_rate", "heterozygosity", "heterozygosity_pc_corrected", "heterozygosity_missing_outlier", "PSCA", "in_kinship", "excluded_kinship_inference", "excess_relatives", "white_british", "pca_calculation", paste0("PC", seq(1,40)), "phasing_autosome", "phasing_X", "phasing_Y")

# Sample eid; This corresponds to a list of all sample eids (with sex information), retrieved from a ".fam" file available from the UKBB portal
sample_eid <- as.data.frame(fread("/data/FAC/FBM/DBC/zkutalik/default_sensitive/UKBB/plink/_001_ukb_cal_chr1_v2.fam", header = F, select = c(1,5), col.names = c("eid", "sex")))


#################################################
### Merge Main Data #############################

df <- cbind(sample_eid, sample_qc)
print(paste0("Start: ", nrow(df), " individuals"))


#################################################
### Sample Filtering ############################

#################################################
### STEP 1: Exclude related samples (pca_calculation = 1)
df <- df[which(df$pca_calculation == 1), ]
print(paste0("STEP 1: Exclude related samples: ", nrow(df), " individuals"))


#################################################
### STEP 2: Exclude non-white, non-British ancestry samples (white_british = 1)
df <- df[which(df$white_british == 1), ]  
print(paste0("STEP 2: Exclude non-white, non-British ancestry samples: ", nrow(df), " individuals"))


#################################################
### STEP 3: Exclude redacted samples & sample for which CNVs were not called (eid = 6026310)
df <- df[-which(df$eid < 0), ]  
df <- df[-which(df$eid == 6026310), ]
print(paste0("STEP 3: Exclude redacted samples: ", nrow(df), " individuals"))


#################################################
### STEP 4: Exclude retracted samples; This file was available from the UKBB portal 
retracted <- as.data.frame(fread("/data/FAC/FBM/DBC/zkutalik/default_sensitive/UKBB/org/w16389_20200820.csv", header = F, col.names = "eid"))
df <- df[!df$eid %in% retracted$eid, ]
print(paste0("STEP 4: Exclude retracted samples: ", nrow(df), " individuals"))


#################################################
### STEP 5: Exclude genotype plate outliers (samples genotyped on plates with a mean CNV count per sample > 100); Sample eids are contained in "plate_outliers.txt"
plate_outlier <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/01_samples/white_british/data/raw/plate_outliers.txt", header = F, select = c(1), col.names = "eid"))
df <- df[!df$eid %in% plate_outlier$eid, ]
print(paste0("STEP 5: Exclude genotype plate outliers: ", nrow(df), " individuals"))


#################################################
### STEP 6: Exclude extreme CNV profiles (samples with >200 CNVs or or 1 CNV > 10Mb); "ukb_cnv.gz" is the final PennCNV output, with all CNV calls in a linear format

# CNVs > 10Mb
cnv <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/cnv_calls/PennCNV/final/ukb_cnv.gz", header = F, select = c(3,5), col.names = c("Length", "File")))
cnv$Length <- as.numeric(gsub("length=|,", "", cnv$Length))
length_outlier <- cnv[which(cnv$Length > 10000000), "File"]
print(paste0("There are ", length(unique(length_outlier)), " samples with at least 1 CNV > 10Mb"))

# > 200 CNVs
cnv_sqc <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/cnv_calls/PennCNV/final/ukb_cnv_sqc.gz", header = T, select = c(1,10)))
num_outlier <- cnv_sqc[which(cnv_sqc$NumCNV > 200), "File"]
print(paste0("There are ", length(unique(num_outlier)), " samples with > 200 CNVs"))

# Overlap of both criteria
cnv_outlier <- as.numeric(gsub("_.*", "", union(length_outlier, num_outlier)))
print(paste0("There are ", length(unique(cnv_outlier)), " CNV outlier samples"))

# Exclude individuals
df <- df[!df$eid %in% cnv_outlier, ]
print(paste0("STEP 6: Exclude extreme CNV profiles (>200 CNVs or or 1 CNV >10Mb): ", nrow(df), " individuals"))


#################################################
### STEP 7: Exclude samples with non-matching submitted vs. inferred sex
df <- df[which(df$submitted_gender == df$inferred_gender), ]
print(paste0("STEP 7: Exclude sex mismatches: ", nrow(df), " individuals"))


#################################################
## STEP 8: Exclude samples with slef-reported blood malignancies or ICD10 blood malignancy diagnosis
print("STEP 8: Exclude samples with blood malignancies")

# UKBB field 20001 (self-reported cancers), codes: 1047 (lymphoma), 1048 (leukemia), 1050 (multiple myeloma), 1051 (myelofibrosis), 1052 (Hodgkin's lymphoma), 1053 (Non-Hodgkin's lymphoma), 1055 (chronic lymphocytic), 1056 (chronic myeloid), 1058 (other hematological malignancies)
SRBM_code <- "1047|1048|1050|1051|1052|1053|1055|1056|1058"

# All ICD10 diagnoses mapping to the PheCode "cancer of lymphatic and hematopoietic tissue"; "phecode_ICD10_091220.csv" was downoladed from https://phewascatalog.org/phecodes_icd10 on 09/12/2020 
phecode <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/general_data/PheCode/phecode_ICD10_091220.csv", header = T, select = c(1,6)))
DBM_code <- phecode[which(phecode$Excl_Phenotypes == "cancer of lymphatic and hematopoietic tissue"), "ICD10"]
DBM_code<- paste(sub("\\.", "", DBM_code[grep("\\.", DBM_code)]), collapse = "|")

# Select appropriated columns; "ukb44073.csv" contains phenotypic data for field 20001 (self-reported cancer) and 41270 (ICD10 diagnoses) available from the UKBB portal
header_pheno <- fread("/data/FAC/FBM/DBC/zkutalik/default_sensitive/UKBB/pheno/ukb44073.csv", nrow = 0)
col_SRBM <- grep("20001-", names(header_pheno))
col_DBM <- grep("41270-", names(header_pheno))

# Self-reported blood malignancy
SRBM <- as.data.frame(fread("/data/FAC/FBM/DBC/zkutalik/default_sensitive/UKBB/pheno/ukb44073.csv", select = c(1, col_SRBM), header = T, colClass = "character"))
SRBM[SRBM == ""] <- NA
SRBM <- unite(SRBM, code, names(SRBM[2:ncol(SRBM)]), sep = "_", remove = T, na.rm = T)
SRBM <- SRBM %>% filter(str_detect(code, SRBM_code))
print(paste0("Number of samples with self-reported blood malignancy: ", nrow(SRBM)))

# ICD10 diagnosis
DBM <- as.data.frame(fread("/data/FAC/FBM/DBC/zkutalik/default_sensitive/UKBB/pheno/ukb44073.csv", select = c(1, col_DBM), header = T))
DBM[DBM == ""] <- NA
DBM <- unite(DBM, code, names(DBM[2:ncol(DBM)]), sep = "_", remove = T, na.rm = T)
DBM <- DBM %>% filter(str_detect(code, DBM_code))
print(paste0("Number of samples ICD10 diagnosed blood malignancy: ", nrow(DBM)))

# Total number of blood malignancies
BM <- union(SRBM$eid, DBM$eid)
print(paste0("Total number of samples with blood malignancy: ", length(BM)))

# Exclude individuals
df <- df[!df$eid %in% BM, ]
print(paste0("STEP 9: Exclude samples with blood malignancies: ", nrow(df), " individuals"))



#################################################
### Save Data ###################################

# For use with PLINK v2.0 - All
fwrite(df[, "eid", drop = F], "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/01_samples/white_british/data/final/samples_white_british_All.txt", col.names = F, row.names = F, quote = F, sep = "\t")
# Males only
fwrite(df[which(df$sex == 1), "eid", drop = F], "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/01_samples/white_british/data/final/samples_white_british_M.txt", col.names = F, row.names = F, quote = F, sep = "\t")
# Females only
fwrite(df[which(df$sex == 2), "eid", drop = F], "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/01_samples/white_british/data/final/samples_white_british_F.txt", col.names = F, row.names = F, quote = F, sep = "\t")

# For use with PLINK v1.9 - All
df_plink1.9 <- data.frame(FID = 0, IID = df$eid, sex = df$sex) 
fwrite(df_plink1.9[, c(1,2), drop = F], "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/01_samples/white_british/data/final/samples_white_british_plink1.9_All.txt", col.names = F, row.names = F, quote = F, sep = "\t")
# Males only
fwrite(df_plink1.9[which(df_plink1.9$sex == 1), c(1,2), drop = F], "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/01_samples/white_british/data/final/samples_white_british_plink1.9_M.txt", col.names = F, row.names = F, quote = F, sep = "\t")
# Females only
fwrite(df_plink1.9[which(df_plink1.9$sex == 2), c(1,2), drop = F], "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/01_samples/white_british/data/final/samples_white_british_plink1.9_F.txt", col.names = F, row.names = F, quote = F, sep = "\t")



