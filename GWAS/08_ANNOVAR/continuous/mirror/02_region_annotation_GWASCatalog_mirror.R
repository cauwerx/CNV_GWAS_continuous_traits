# Annotate CNVRs with GWAS Catalog entries

#################################################
### Libraries ###################################
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)


#################################################
### STEP 1: Load CNVRs ##########################

# Load CNVRs; "GWAS_CNVR_mirror_ANNOVAR_input.txt" from GWAS/07_CNVR/continuous/mirror/extract_CNVR_mirror.R contains CNVR boundaries in ANNOVAR input format
cnvr <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/09_CNVR/continuous/white_british/mirror/data/final/GWAS_CNVR_mirror_ANNOVAR_input.txt", select = c(1:7), col.names = c("CHR", "START", "END", "REF", "OBS", "PHENO", "ID")))

# This file contains the phenotype name used in the previous analyses as well as all GWAS Catalog synonyms for that traits (as in Table S1) separated by a comma
pheno_synonym <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/10_ANNOVAR/continuous/raw/continuous_traits_GWAS_Catalog.txt", sep = "\t", header = T, col.names = c("PHENO", "SYN")))


#################################################
### STEP 2: Loop over CNVRs and annotate ########

# Counter for the number of matches
counter <- 0 

# Loop over CNVRs
for (i in 1:nrow(cnvr)) {

	# Define the region and synonyms
	chr <- cnvr[i, "CHR"]
	start <- cnvr[i, "START"]
	end <- cnvr[i, "END"]
	pheno <- cnvr[i, "PHENO"]
	synonyms <- unlist(strsplit(as.character(pheno_synonym[which(pheno_synonym$PHENO == pheno), "SYN"]), ","))
	print(paste0("Analyzing chr", chr, ":", start, "-", end, " - ", pheno))
	print(paste0(synonyms))

	# Define paths
	input <- paste0("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/10_ANNOVAR/continuous/white_british/mirror/data/temp/GWASCatalog/chr", chr, ":", start, "-", end, "_", pheno, ".input.txt")
	output <- paste0("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/10_ANNOVAR/continuous/white_british/mirror/data/temp/GWASCatalog/chr", chr, ":", start, "-", end, "_", pheno)

	# Save temporary ANNNOVAR input
	fwrite(cnvr[i, ], input, col.names = F, row.names = F, quote = F, sep = "\t")

	# Run region annotation for GWAS catalog in ANNOVAR
	system(paste("/home/cauwerx/scratch/cauwerx/softwares/ANNOVAR-20191024/annotate_variation.pl",
			 "--regionanno -dbtype gwasCatalog",
			 "-out", output, 
		     "-build hg19", 
			 input,
			 "/home/cauwerx/scratch/cauwerx/softwares/ANNOVAR-20191024/humandb/"))	
	unlink(paste0(output, ".log"))
	unlink(paste0(input))

	# Load ANNOVAR annotation & filter for correspnding terms
	annot <- as.data.frame(fread(paste0(output,".hg19_gwasCatalog"), sep = "\t", select = c(2:5,8), col.names = c("GWAS", "CHR", "START", "END", "PHENO")))
	if (nrow(annot) == 0) {
		cnvr[i, "GWAS"] <- 0
		print("No GWAS hits in the region")
	} else {	
		annot <- sub("Name=", "", annot$GWAS)
		annot <- unlist(strsplit(annot, ","))
		print(paste0("Number of unique traits with GWAS hits in the region: ", length(annot)))
		match_annot <- intersect(tolower(synonyms), tolower(annot))
		if (length(match_annot) == 0) {
			cnvr[i, "GWAS"] <- 0
			print(paste0("No GWAS hits for ", pheno, " in the region"))
		} else {
			cnvr[i, "GWAS"] <- 1
			print(paste0("GWAS hits for ", pheno," in the region"))
			counter <- counter + 1
		}
	}
	rm(input, output, annot, annot_match)
	
}


#################################################
### Save ########################################

# Print global result summary
print(paste0("Total number of trait-CNVR with a match reported in GWAS catalog: ", counter, "/", nrow(cnvr)))

# Save
fwrite(cnvr, "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/10_ANNOVAR/continuous/white_british/mirror/data/final/region_annotation_GWASCatalog/GWASCatalog_matching_pheno_annotation_mirror_CNVR.txt", col.names = T, row.names = F, quote = F, sep = "\t")
