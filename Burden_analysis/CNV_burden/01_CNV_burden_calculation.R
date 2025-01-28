# Calculate each individual's CNV burden

#################################################
### Libraries ###################################
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)


#################################################
### STEP 1: Load data ###########################

# Load all CNVs; "ukb_cnv_global.gz" is the final PennCNV output, with all CNV calls in a linear format
cnvs_all <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/cnv_calls/final/ukb_cnv_global.gz", header = T))
cnvs_all$Sample_Name <- sub("_.*", "", cnvs_all$Sample_Name)
print(paste0("Number of loaded CNVs: ", nrow(cnvs_all)))

# Load sex info; "eid_no_redacted3.txt.gz" contains all samples eids with sex information, except for redacted samples
eids <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/general_data/UKBB/eid_no_redacted3.txt.gz", header = T, select = c(2,5), col.names = c("IID", "SEX")))

# All unrelated, white British samples from GWAS/01_samples/sample_filtering.R
cnv_burden <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/01_samples/white_british/data/final/samples_white_british_All.txt", header = F, col.names = "IID"))
cnv_burden <- left_join(cnv_burden, eids, by = "IID")


#################################################
### STEP 2: Exclude non-autosomal CNVs ##########

cnvs_all <- cnvs_all[!cnvs_all$Chromosome %in% c(23, 24),]


#################################################
### STEP 3: Calculate the burden ################

#################################################
### BURDEN 1: Mb for |QS| > 0.5 #################

# Filter CNVs with |QS| > 0.5
cnvs <- cnvs_all[which(abs(cnvs_all$Quality_Score) >= 0.5), ]
print(paste0("Number of CNVs with |QS| >= 0.5: ", nrow(cnvs)))

# Calculate the sum of affected Mb
cnv_burden_Mb <- data.frame(IID = as.numeric(cnvs$Sample_Name), BURDEN_Mb = cnvs$Length_bp/1000000)
cnv_burden_Mb <- aggregate(BURDEN_Mb ~ IID, cnv_burden_Mb, sum)
cnv_burden <- left_join(cnv_burden, cnv_burden_Mb, by = "IID")


#################################################
### BURDEN 2: Gene number for |QS| > 0.5 ########

# Save a temporary file compatible with ANNOVAR format
cnvs <- data.frame(CHR = cnvs$Chromosome, START = cnvs$Start_Position_bp, END = cnvs$End_Position_bp, REF = 0, OBS = "-", SAMPLE = cnvs$Sample_Name, LENGTH = cnvs$Length_bp, CN = cnvs$Copy_Number)
fwrite(cnvs, "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/05_phenotypes/CNV_burden/white_british/data/temp/CNV_ANNOVAR_input.txt", col.names = F, row.names = F, quote = F, sep = "\t")

# Run ANNOVAR
input <- "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/05_phenotypes/CNV_burden/white_british/data/temp/CNV_ANNOVAR_input.txt"
output <- "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/05_phenotypes/CNV_burden/white_british/data/temp/CNV_ANNOVAR_output"
system(paste("/home/cauwerx/scratch/cauwerx/softwares/ANNOVAR-20191024/annotate_variation.pl",
			 "--geneanno -dbtype refGene",
			 "-out", output, 
		     "-build hg19", 
			 input,
			 "/home/cauwerx/scratch/cauwerx/softwares/ANNOVAR-20191024/humandb/"))	
unlink(paste0(output, ".log"))

# Load output file
cnvs <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/05_phenotypes/CNV_burden/white_british/data/temp/CNV_ANNOVAR_output.variant_function", header = F, select = c(1:5, 8:10), col.names = c("REGION", "GENE", "CHR", "START", "END", "IID", "LENGTH", "CN")))
print(paste0("Number of annotated CNVs: ", nrow(cnvs)))

# Filter for CNVs affecting the genes (stringent)
regions <- c("exonic", "splicing", "ncRNA", "UTR5", "UTR3")
cnvs <- cnvs[grep(paste(regions, collapse = "|"), cnvs$REGION), ]
print(paste0("Remaining CNVs after filtering for gene impact (stringent): ", nrow(cnvs)))

# Group by SAMPLE
cnv_burden_gene <- aggregate(GENE ~ IID, cnvs[, c("IID", "GENE")], paste, collapse = ",")

# Count the number of unique affected genes
cnv_burden_gene$BURDEN_GENES <- lapply(lapply(strsplit(cnv_burden_gene$GENE, ","), unique), length)
cnv_burden_gene$BURDEN_GENES <- as.numeric(cnv_burden_gene$BURDEN_GENES)

# Merge to main file
cnv_burden <- left_join(cnv_burden, cnv_burden_gene[, c("IID", "BURDEN_GENES")], by = "IID")


#################################################
### STEP 4: Correlation #########################

# Set missing values to 0
cnv_burden[is.na(cnv_burden)] <- 0  
cnv_burden[which(cnv_burden$BURDEN_GENES == "NULL"), "BURDEN_GENES"] <- 0 

# print correlation among Burden measures
print(cor(cnv_burden[, -1]))


#################################################
### STEP 5: Sex differences #####################

sex_diff <- data.frame(BURDEN = c("BURDEN_Mb", "BURDEN_GENES"))

# BURDEN_Mb
t_Mb <- wilcox.test(cnv_burden[which(cnv_burden$SEX == 1), "BURDEN_Mb"], cnv_burden[which(cnv_burden$SEX == 2), "BURDEN_Mb"], alternative = "two.sided", paired = F, conf.int = T)
sex_diff[2, "DIFF"] <- t_Mb$estimate[1]
sex_diff[2, "T"] <- t_Mb$statistic
sex_diff[2, "P"] <- t_Mb$p.value

# BURDEN_GENES
t_GENES <- wilcox.test(cnv_burden[which(cnv_burden$SEX == 1), "BURDEN_GENES"], cnv_burden[which(cnv_burden$SEX == 2), "BURDEN_GENES"], alternative = "two.sided", paired = F, conf.int = T)
sex_diff[8, "DIFF"] <- t_GENES$estimate[1]
sex_diff[8, "T"] <- t_GENES$statistic
sex_diff[8, "P"] <- t_GENES$p.value


#################################################
### Save ########################################

# Burden
fwrite(cnv_burden, "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/05_phenotypes/CNV_burden/white_british/data/final/CNV_burden.txt.gz", col.names = T, row.names = F, quote = F, sep = "\t")

# Sex differences
fwrite(sex_diff, "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/05_phenotypes/CNV_burden/white_british/data/final/sex_difference/CNV_burden_sex_diff.txt", col.names = T, row.names = F, quote = F, sep = "\t")
