# Power analysis for the replication study in the EstBB using simulated p-values

#################################################
### Libraries ###################################
library(data.table)
library(dplyr)
library(tibble)


#################################################
### STEP 1: Load file ###########################

# Load replicated duplication-only trait-CNVR association characteristics; "GWAS_CNVR_DUP.withEstResults.withCounts.txt" columns are explained in EstBB/power_analysis/duplication_only/01_simulate_p-val_dup.R
df <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/13_EstBB/raw/GWAS_CNVR_DUP.withEstResults.withCounts.txt", 
						  select = c(1:7, 13, 15:19, 21:23), 
						  col.names = c("ID_UK", "CHR", "POS_UK", "BETA_UK", "SE_UK", "P_UK", "PHENO", "ID_Est", "POS_Est", "N_Est", "BETA_Est", "SE_Est", "P_Est", "NumDEL_Est", "NumDUP_Est", "NumCN_Est")))
df <- df[!is.na(df$BETA_Est), ]


#################################################
### STEP 2: Define significance thresholds ######

# Threshold 1: model-specific MTC
thr1 <- 0.05/nrow(df)
print(paste0("Model-specific MTC p-value threshold (", nrow(df)," tests): ", thr1))

# Threshold 2: top significant model for a given CNVR-trait pair in the UKBB (61 independnet pairs were replicated)
thr2 <- 0.05/61 
print(paste0("Global MTC p-value threshold (", nrow(df)," tests): ", thr2))


#################################################
### STEP 3: Differences in size effect ##########

# Calculate t-stat and corresponding p-value for the difference in effect size between the two cohorts
df <- add_column(df, BETA_T = (df$BETA_UK - df$BETA_Est)/sqrt((df$SE_UK)^2 + (df$SE_Est)^2), .after = "P_Est")
df <- add_column(df, BETA_P = 2*pnorm(-abs(df$BETA_T), mean = 0, sd = 1), .after = "BETA_T")


#################################################
### STEP 4: Power analysis ######################

# Calculate power
for (i in 1:nrow(df)) {

	# Load p-values from simulations
	p_val <- na.omit(as.data.frame(fread(paste0("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/13_EstBB/power_analysis/duplication_only/data/temp/p_values/", df[i, "PHENO"],".", df[i, "CHR"],":", df[i, "POS_UK"],".dup.txt"))))

	# Calculate power as percentage of non-NA p-values <= thr
	df[i, paste0("POWER_0.05/", nrow(df))] <- length(p_val[which(p_val$V1 <= thr1), ])/nrow(p_val)
	df[i, "POWER_0.05/61"] <- length(p_val[which(p_val$V1 <= thr2), ])/nrow(p_val)

}


#################################################
### Save ########################################

fwrite(df, "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/13_EstBB/power_analysis/duplication_only/data/final/power_EstBB_dup.txt", col.names = T, row.names = F, quote = F, sep = "\t")
