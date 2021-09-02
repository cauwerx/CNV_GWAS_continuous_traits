# Merging of CNV calling by PennCNV on the autosomes and chrX

########################################################
# Libraries
########################################################
library(dplyr)
library(data.table)
options(scipen = 999)


########################################################
##### Arguments & QS ###################################

args <- commandArgs(trailingOnly = TRUE)
path_to_PennCNV <- args[1] 	# Path to PennCNV data
batch <- args[2] 			# Batch


########################################################
##### STEP 1: Load data ################################

# Load merged CNVs for autosomes (from CNV_calling/scripts_for_pipeline/05_1_PennCNV_array.sh) and chrX (from CNV_calling/scripts_for_pipeline/05_2_PennCNV_ChrX_array.sh)
merged_cnv_A <- fread(paste0(path_to_PennCNV, "/", batch, "_merged_rawcnv.txt" ), header = F, fill = T)
merged_cnv_X <- fread(paste0(path_to_PennCNV, "/ChrX/", batch, "_chrX_filtered_rawcnv.txt"), header = F, fill = T)
merged_cnv_X$V1 <- gsub("chrX:", "chr24:", merged_cnv_X$V1)

# Load sample QC for autosomes (from CNV_calling/scripts_for_pipeline/05_1_PennCNV_array.sh) and chrX (from CNV_calling/scripts_for_pipeline/05_2_PennCNV_ChrX_array.sh)
sqc_A <- fread(paste0(path_to_PennCNV, "/", batch, "_qc_summary.txt"), header = T, fill = T)
sqc_X <- fread(paste0(path_to_PennCNV, "/ChrX/", batch, "_chrX_qc_summary.txt"), header = T, fill = T)


########################################################
##### STEP 2: Merge CNV files ##########################

rawcnv_AX <- rbind(merged_cnv_A, merged_cnv_X)


########################################################
##### STEP 3: Sum up CNVs from autoosomes/ChrX #########

# Order files 
sqc_A <- sqc_A[order(sqc_A$File), ]
sqc_X <- sqc_X[order(sqc_X$File), ]

# Merge NumCNV
if (identical(sqc_A$File, sqc_X$File)) {
	print("All equal, add up NumCNV")
	sqc_AX <- sqc_A
	sqc_AX$NumCNV <- sqc_A$NumCNV + sqc_X$NumCNV} 


########################################################
##### Save #############################################

write.table(rawcnv_AX, paste0(path_to_PennCNV, "/", batch, "_CNVs_AX.txt"), col.names = F, row.names = F, quote = F, sep = '\t')
write.table(sqc_AX, paste0(path_to_PennCNV, "/", batch, "_SQC_AX.txt"), col.names = T, row.names = F, quote = F, sep = '\t')
print(paste0("Autosome - ChrX file for ", batch, " were successfully merged"))
