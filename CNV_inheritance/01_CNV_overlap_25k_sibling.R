# Calculate the fraction of inherited high confidence CNVs across sibling pairs

#################################################
### Libraries #################################
library(data.table)
library(dplyr)
library(tidyr)


#################################################
### STEP 1: Load data ###########################

# Load unrelated white British individuals; from GWAS/01_samples/sample_filtering.R
eids_unrel <- as.data.frame(fread("//home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/01_samples/white_british/data/final/samples_white_british_All.txt", col.names = "eids"))

# Load kingship data; This corresponds to the table availble on the UKBB portal
kin <- as.data.frame(fread("/data/FAC/FBM/DBC/zkutalik/default_sensitive/UKBB/geno/ukb1638_rel_s488366.dat")) # 107'162 eids

## Load all CNVs; "ukb_cnv_global.gz" is the final PennCNV output, with all CNV calls in a linear format
cnvs <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/cnv_calls/final/ukb_cnv_global.gz"))
cnvs$Sample_Name <- as.numeric(sub("_.*", "", cnvs$Sample_Name))


#################################################
### STEP 2: Detect sibling pairs ################

# Select siblings: Kinship > 0.2 & IBS0 > 0.0012 (as described in Bycroft et al., 2018)
kin <- kin[which(kin$Kinship <= 0.3 & kin$Kinship >= 0.2 & kin$IBS0 >= 0.0012), ] # 22'450

# Determine which individual was in the main analysis
kin[kin$ID1 %in% eids_unrel$eids, "Main1"] <- "ID1"
kin[kin$ID2 %in% eids_unrel$eids, "Main2"] <- "ID2"
kin <- unite(kin, "Main", Main1:Main2, sep = "-", remove = T, na.rm = T) 
kin <- kin[kin$Main %in% c("ID1", "ID2"), ]
# There are 16'179 pairs with 1 sib in the main analysis


#################################################
### STEP 3: Calculate number of shared CNVs #####

# Select high-confidence CNVs
cnvs <- cnvs[which(abs(cnvs$Quality_Score) >= 0.5), ]

# Define the minimal overlap size
OverlapFilter <- 25000

# Set counters to 0
NumOverlap <- 0 # Total number of overlaps
NumOverlapLost <- 0 # Number of overlaps lost due to minimal overlap size filter

# Loop over sibling pairs
for (i in 1:nrow(kin)) { 

	#### DEFINE PAIRS ################

	# Define the pairs: Main (sample in the main GWAS/burden analysis) & Sib (Sibling of the main sample)
	Main <- kin[i, "Main"]
	if (Main == "ID1") {Sib <- "ID2"} 
	if (Main == "ID2") {Sib <- "ID1"}

	# Define the identifiers
	ID_Main <- kin[i, Main]
	ID_Sib <- kin[i, Sib]

	# Create CNV tables for both individuals by CNV type
	CNV_Main <- as.data.table(cnvs[which(cnvs$Sample_Name == ID_Main), c(2,3,4,5)])
	CNV_Main$TYPE <- "DUP"
	CNV_Main[which(CNV_Main$Copy_Number < 2), "TYPE"] <- "DEL"
	CNV_Sib <- as.data.table(cnvs[which(cnvs$Sample_Name == ID_Sib), c(2,3,4,5)])
	CNV_Sib$TYPE <- "DUP"
	CNV_Sib[which(CNV_Sib$Copy_Number < 2), "TYPE"] <- "DEL"

	# Calculate the number of CNVs per pair
	kin[i, "NumCNV_Main"] <- nrow(CNV_Main)
	kin[i, "NumCNV_Sib"] <- nrow(CNV_Sib)

	
	#### OVERLAP ######################	

	# Situation 1: No CNV in any of the individuals in the pair --> NA
	if (kin[i, "NumCNV_Main"] == 0 & kin[i, "NumCNV_Sib"] == 0) {
			kin[i, "NumOverlap"] <- 0
			kin[i, "FractionOverlap_Main"] <- kin[i, "NumOverlap"]/kin[i, "NumCNV_Main"]
			kin[i, "FractionOverlap_Sib"] <- kin[i, "NumOverlap"]/kin[i, "NumCNV_Sib"]
	}

	# Situation 2: Only one individual has CNV(s) --> 0 if CNV(s) in the Main, NA if CNV(s) in the Sib
	if (kin[i, "NumCNV_Main"] > 0 & kin[i, "NumCNV_Sib"] == 0 | kin[i, "NumCNV_Main"] == 0 & kin[i, "NumCNV_Sib"] > 0) {
			kin[i, "NumOverlap"] <- 0
			kin[i, "FractionOverlap_Main"] <- kin[i, "NumOverlap"]/kin[i, "NumCNV_Main"]
			kin[i, "FractionOverlap_Sib"] <- kin[i, "NumOverlap"]/kin[i, "NumCNV_Sib"]
	}

	# Situation 3: Both individuals have CNVs --> Calculate overlap
	if (kin[i, "NumCNV_Main"] > 0 & kin[i, "NumCNV_Sib"] > 0) {

		# Setkeys
		setkey(CNV_Main, Chromosome, TYPE, Start_Position_bp, End_Position_bp)
		setkey(CNV_Sib, Chromosome, TYPE, Start_Position_bp, End_Position_bp)

		# Calculate overlap (min 1 bp overlap)
		overlap <- na.omit(foverlaps(CNV_Sib, CNV_Main, type = "any", maxgap = 0L, minoverlap = 1L))
		ol_original <- nrow(overlap)
		NumOverlap <- NumOverlap + ol_original

		# Calculate overlap length and discard if the overlap is < 25kb
		for (j in 1:nrow(overlap)){
			val <- sort(c(as.numeric(overlap[j, 3]), as.numeric(overlap[j, 4]), as.numeric(overlap[j, 6]), as.numeric(overlap[j, 7])))
			overlap[j, "Overlap_Length"] <- val[3]-val[2]
		}
		overlap <- overlap[which(overlap$Overlap_Length >= OverlapFilter), ]
		ol_filter <- nrow(overlap)
		NumOverlapLost <- NumOverlapLost + (ol_original - ol_filter)

		# Fill in the table with the number of overlapping CNVs (NumOverlap), the fraction of overlap from the perspective of the main indiividual (FractionOverlap_Main)) and sibling (FractionOverlap_Sib)
		kin[i, "NumOverlap"] <- nrow(overlap)
		kin[i, "FractionOverlap_Main"] <- kin[i, "NumOverlap"]/kin[i, "NumCNV_Main"]
		kin[i, "FractionOverlap_Sib"] <- kin[i, "NumOverlap"]/kin[i, "NumCNV_Sib"]
	}

}


#################################################
### STEP 4: Calculate shared fraction of CNVs ###

# Print result summary
print(paste0("Unique number of pairs: ", nrow(kin)))
print(paste0("Unique number of individuals: ", length(unique(c(kin$ID1, kin$ID2)))))
print(paste0("Total number of CNVs (Main): ", sum(kin$NumCNV_Main, na.rm = T), " ;Mean number of CNVs (Main): ", mean(kin$NumCNV_Main, na.rm = T)))
print(paste0("Total number of CNVs (Sib): ", sum(kin$NumCNV_Sib, na.rm = T), " ;Mean number of CNVs (Sib): ", mean(kin$NumCNV_Sib, na.rm = T)))
print(paste0("Total number of overlaps: ", sum(kin$NumOverlap, na.rm = T), " ;Mean number of overlaps: ", mean(kin$NumOverlap, na.rm = T)))
print(paste0(NumOverlap, " overlaps were detected, ", NumOverlapLost, " of which were lost as the overlap was smaller than ", OverlapFilter, " kb"))
print(paste0("Mean fraction of overlap (Main): ", mean(kin$FractionOverlap_Main, na.rm = T), " ;Median fraction of overlap (Main): ", median(kin$FractionOverlap_Main, na.rm = T)))
print(paste0("Mean fraction of overlap (Sib): ", mean(kin$FractionOverlap_Sib, na.rm = T), " ;Median fraction of overlap (Sib): ", median(kin$FractionOverlap_Sib, na.rm = T)))


#################################################
### Save ########################################

fwrite(kin, "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/frequency_analysis/02_CNV_inheritance/white_british/data/final/shared_CNVs/shared_CNV_siblings_overlap_25k.txt", col.names = T, row.names = F, quote = F, sep = "\t")
