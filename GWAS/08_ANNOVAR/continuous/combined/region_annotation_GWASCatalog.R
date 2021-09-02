# Combine GWAS Catalog CNVR annottations from different association models

################################################
### Libraries ##################################
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)


#################################################
### STEP 1: Load CNVRs and annotations ##########

# Load CNVR annotations for the mirror model; from GWAS/08_ANNOVAR/continuous/mirror/02_region_annotation_GWASCatalog_mirror.R
df <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/10_ANNOVAR/continuous/white_british/mirror/data/final/region_annotation_GWASCatalog/GWASCatalog_matching_pheno_annotation_mirror_CNVR.txt", header = T, select = c(1:3, 6, 8)))

# Load CNVR annotations for the duplication-only model; from GWAS/08_ANNOVAR/continuous/duplication_only/02_region_annotation_GWASCatalog_dup.R
dup <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/10_ANNOVAR/continuous/white_british/duplication_only/data/final/region_annotation_GWASCatalog/GWASCatalog_matching_pheno_annotation_dup_CNVR.txt", header = T, select = c(1:3, 6, 8)))

# Load CNVR annotations for the deletion-only model; from GWAS/08_ANNOVAR/continuous/deletion_only/02_region_annotation_GWASCatalog_del.R
del <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/10_ANNOVAR/continuous/white_british/deletion_only/data/final/region_annotation_GWASCatalog/GWASCatalog_matching_pheno_annotation_del_CNVR.txt", header = T, select = c(1:3, 6, 8)))


#################################################
### STEP 2: Extract data ########################

# Extract CNVR with GWAS hits for each association model
df <- df[which(df$GWAS == 1), ]
dup <- dup[which(dup$GWAS == 1), ]
del <- del[which(del$GWAS == 1), ]

# Add model information to each dataframe
df$MODEL <- "M"
dup$MODEL <- "DUP"
del$MODEL <- "DEL"


#################################################
### STEP 3: Add DUP data to DF (mirror) #########

# Loop over duplication-only CNVRs
for (i in 1:nrow(dup)) {

	# Define the CNVR	
	chr <- dup[i, "CHR"]
	pheno <- dup[i, "PHENO"]
	start <- dup[i, "START"]
	end <- dup[i, "END"]
	print(paste0("Duplication signal: chr", chr, ":", start, "-", end, " - ", pheno))

	# Look for matches from the mirror model
	match <- df[which(df$CHR == chr & df$PHENO == pheno), ]

	# If no match --> Add row with this CNVR to the dataframe
	if (nrow(match) == 0) {
		df <- rbind(df, dup[i,])
		print("New signal! Add to the list")

	# Else, check for region overlap
	} else {
		id_match <- vector()
		for (j in 1:nrow(match)) {
			if (match[j, "START"] >= start & match[j, "END"] <= end | match[j, "START"] <= start & match[j, "END"] >= end | match[j, "START"] <= end & match[j, "END"] >= end | match[j, "START"] <= start & match[j, "END"] >= start) {id_match <- j}
		}
		# If no overlap --> Add row with this CNVR to the dataframe
		if(length(id_match) == 0) {
			df <- rbind(df, dup[i,])
			print("New signal! Add to the list")		
		
		# If overlap --> Remove the old row and add a modified row combining information from both association models
		} else {
			# Remove old row
			df <- df[-which(df$CHR == chr & df$PHENO == pheno & df$START == match[j, "START"] & df$END == match[j, "END"]), ]
			# Add a new one
			df[nrow(df) +1, "CHR"] <- chr
			df[nrow(df), "START"] <- min(match[j, "START"], start)
			df[nrow(df), "END"] <- max(match[j, "END"], end)
			df[nrow(df), "PHENO"] <- pheno
			df[nrow(df), "GWAS"] <- 1
			df[nrow(df), "MODEL"] <- paste0(match[j, "MODEL"], "-", dup[i, "MODEL"])
			print("Modify existing signal")

		}
	}
}

print(paste0("There is a total of ", nrow(df), " signals after addition of dup signals"))


#################################################
### STEP 4: Add DEL data to DF (mirror + DUP) ###

# Loop over deletion-only CNVRs
for (i in 1:nrow(del)) {

	# Define the CNVR	
	chr <- del[i, "CHR"]
	pheno <- del[i, "PHENO"]
	start <- del[i, "START"]
	end <- del[i, "END"]
	print(paste0("Deletion signal: chr", chr, ":", start, "-", end, " - ", pheno))

	# Look for matches from the mirror/duplication-only models
	match <- df[which(df$CHR == chr & df$PHENO == pheno), ]

	# If no match --> Add row with this CNVR to the dataframe
	if (nrow(match) == 0) {
		df <- rbind(df, del[i,])
		print("New signal! Add to the list")

	# Else, check for region overlap
	} else {
		id_match <- vector()
		for (j in 1:nrow(match)) {
			if (match[j, "START"] >= start & match[j, "END"] <= end | match[j, "START"] <= start & match[j, "END"] >= end | match[j, "START"] <= end & match[j, "END"] >= end | match[j, "START"] <= start & match[j, "END"] >= start) {id_match <- j}
		}
		# If no overlap --> Add row with this CNVR to the dataframe
		if(length(id_match) == 0) {
			df <- rbind(df, del[i,])
			print("New signal! Add to the list")		
		
		# If overlap --> Remove the old row and add a modified row combining information from all association models
		} else {
			# Remove old row
			df <- df[-which(df$CHR == chr & df$PHENO == pheno & df$START == match[j, "START"] & df$END == match[j, "END"]), ]
			# Add a new one
			df[nrow(df) +1, "CHR"] <- chr
			df[nrow(df), "START"] <- min(match[j, "START"], start)
			df[nrow(df), "END"] <- max(match[j, "END"], end)
			df[nrow(df), "PHENO"] <- pheno
			df[nrow(df), "GWAS"] <- 1
			df[nrow(df), "MODEL"] <- paste0(match[j, "MODEL"], "-", del[i, "MODEL"])
			print("Modify existing signal")

		}
	}
}

print(paste0("There is a total of ", nrow(df), " signals after addition of dup & del signals"))


#################################################
### Save ########################################

fwrite(df, "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/10_ANNOVAR/continuous/white_british/combined/data/final/region_annotation_GWASCatalog/combined_CNVR_phenotype_GWASCatalog.txt", col.names = T, row.names = F, quote = F, sep = "\t")
