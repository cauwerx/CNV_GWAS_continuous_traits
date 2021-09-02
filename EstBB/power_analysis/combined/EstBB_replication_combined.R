# Combine results of power analysis to the most significant signal in the UKBB

#################################################
### Libraries ###################################
library(dplyr)
library(tidyr)
library(data.table)


#################################################
### STEP 1: Load files ##########################

# Load the top model for each signal in the UKBB; "UKBB_top_assoc.txt" contains in the following order:
# phenotype, chromosome, start and end position of the CNVR, all models through which the association was found to be significant
# most significant model and the associated p-value, probe ID, sign (direction of the association), cytogenic band
top_uk <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/13_EstBB/power_analysis/combined/data/raw/UKBB_top_assoc.txt"))

# Load data from the power analysis in the EstBB for mirror association; from EstBB/power_analysis/mirror/02_power_calculation_mirror.R
mirror <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/13_EstBB/power_analysis/mirror/data/final/power_EstBB_mirror.txt"))
mirror$MODEL <- "M"

# Load data from the power analysis in the EstBB for duplication-only association; from EstBB/power_analysis/duplication_only/02_power_calculation_dup.R
dup <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/13_EstBB/power_analysis/duplication_only/data/final/power_EstBB_dup.txt"))
dup$MODEL <- "DUP"

# Load data from the power analysis in the EstBB for deletion-only association; from EstBB/power_analysis/deletion_only/02_power_calculation_del.R
del <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/13_EstBB/power_analysis/deletion_only/data/final/power_EstBB_del.txt"))
del$MODEL <- "DEL"


#################################################
### STEP 2: Select top model ####################

# Create an empty dataframe
df <- data.frame()

# Loop over 131 idependent signals in the UKBB
for (i in 1:nrow(top_uk)) {
  
  # Define the top association in the UKBB
  pheno <- top_uk[i, "PHENO"]
  chr <- top_uk[i, "CHR"]
  type <- top_uk[i, "TOP_TYPE"]
  rs <- top_uk[i, "TOP_ID"]
  print(paste0("UKBB Signal: chr", chr, ":", rs, " - ", pheno, " (", type, ")"))
  
  # Identify if there are matches in the EstBB replication data
  if (type == "M") {
    df_temp <- mirror[which(mirror$PHENO == pheno & mirror$CHR == chr & mirror$ID_UK == rs), ]
    if(nrow(df_temp) == 1) {df <- rbind(df, df_temp[, -c(19)])}
    
   } else if (type == "DUP" | type == "M-DUP") {
    df_temp <- dup[which(dup$PHENO == pheno & dup$CHR == chr & dup$ID_UK == rs), ]
    if(nrow(df_temp) == 1) {df <- rbind(df, df_temp[, -c(19)])}
    
  } else if (type == "DEL") {
    df_temp <- del[which(del$PHENO == pheno & del$CHR == chr & del$ID_UK == rs), ]
    if(nrow(df_temp) == 1) {df <- rbind(df, df_temp[, -c(19)])}
  }
}

# Merge the selected signals from the EstBB to the top UKBB file
df$CHR <- as.character(df$CHR)
df <- right_join(top_uk[, c(1:6, 8:10)], df[, c(1:9, 11:19)], by = c("PHENO", "CHR", "TOP_ID" = "ID_UK"))
df <- df[, c(1:4, 9, 5:6, 8, 7, 10:24)]
colnames(df)[9] <- "ID_UK"


######################################################
### STEP 3: Adjust p-val for effect size direction ###

# If effect sizes agree in directiion, P_COR = p(old)/2
df[which(df$BETA_UK < 0 & df$BETA_Est < 0 | df$BETA_UK > 0 & df$BETA_Est > 0), "P_COR"] <- df[which(df$BETA_UK < 0 & df$BETA_Est < 0 | df$BETA_UK > 0 & df$BETA_Est > 0), "P_Est"]/2

# If effect sizes DO NOT agree in direction, P_COR = 1 - p(old)/2
df[which(df$BETA_UK < 0 & df$BETA_Est > 0 | df$BETA_UK > 0 & df$BETA_Est < 0), "P_COR"] <- 1-(df[which(df$BETA_UK < 0 & df$BETA_Est > 0 | df$BETA_UK > 0 & df$BETA_Est < 0), "P_Est"]/2)

# Indicate nominally and MTC surviving significant associations
df$SIG <- 0
df[which(df$P_COR <= 0.05), "SIG"] <- 1 
df[which(df$P_COR <= 0.05/nrow(df)), "SIG"] <- 2 
df$SIG <- factor(df$SIG)

# Calculate enrichment for nominally significant signals
ExpNSig <- 0.05*nrow(df) # 3.05
NumNSig <- nrow(df[which(df$SIG == 1 | df$SIG == 2), ]) # 22
FoldEnrichNSig <- NumNSig/ExpNSig # 7.2
p_binom_EnrichNSig <- binom.test(NumNSig, nrow(df), p = 0.05) # 7.835e-14


######################################################
### Save #############################################

fwrite(df, "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/13_EstBB/power_analysis/combined/data/final/UKBB_EstBB_top_assoc.txt", col.names = F, row.names = F, quote = F, sep = "\t")
