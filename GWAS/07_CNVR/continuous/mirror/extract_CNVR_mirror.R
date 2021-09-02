# Determine CNVR for each lead signal

#################################################
### Libraries ###################################
library(data.table)
library(dplyr)


#################################################
### STEP 1: Load data ###########################

# Load lead signals from SCA for all; from GWAS/06_SCA_continuous/mirror/SCA_mirror_All.R
files_c <- list.files("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/08_SCA/continuous/white_british/mirror/data/temp/summary/All/", recursive = F, full.names = T)
df_c <- as.data.frame(fread(files_c[length(files_c)], header = T, select = c(1:6,8)))

# Load lead signals from SCA for males; from GWAS/06_SCA_continuous/mirror/SCA_mirror_M.R
files_M <- list.files("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/08_SCA/continuous/white_british/mirror/data/temp/summary/M/", recursive = F, full.names = T)
df_M <- as.data.frame(fread(files_M[length(files_M)], header = T, select = c(1:6,8)))
df_M <- df_M[df_M$PHENO %in% c("balding", "facial_hair"),  ]

# Load lead signals from SCA for females; from GWAS/06_SCA_continuous/mirror/SCA_mirror_F.R
files_F <- list.files("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/08_SCA/continuous/white_british/mirror/data/temp/summary/F/", recursive = F, full.names = T)
df_F <- as.data.frame(fread(files_F[length(files_F)], header = T, select = c(1:6,8)))
df_F <- df_F[df_F$PHENO %in% c("menopause", "menarche", "birth_weight_first_child"),  ]

# Merge signals
df <- rbind(df_c, df_M, df_F)
df <- df[, c(3,1:2,4:7)]
print(paste0("Number of independent signals: ", nrow(df)))


#################################################
### STEP 2: Add field ID ########################

# Add field ID to each CNV signal
field_id <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/05_phenotypes/continuous/white_british/data/raw/continuous_traits.csv", header = T, col.names = c("FIELD_ID", "PHENO")))
df <- left_join(df, field_id, by = "PHENO")


#################################################
### STEP 3: Calculate the CNVR ##################

# Unique combination of chromosome-lead variant (remove duplicates)
unique_combo <- df[!duplicated(df[, c(1,2)]), c(1,2)]
print(paste0("Number of probes to analyze: ", nrow(unique_combo)))

# Load all probes and correct sex chromosome annotation; "probes.txt.gz" contains probe rs number, chromosome, and genomic position
probes <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/general_data/UKBB/probes.txt.gz", header = T, col.names = c("ID", "CHR", "POS")))
probes$CHR <- sub("23", "X", probes$CHR)

# Loop over each probe: identify adjacent probes in LD
for(i in 1:nrow(unique_combo)) {

	# Define settings
	rs <- unique_combo[i, 1]
	chr <- unique_combo[i, 2]
	print(paste0("Starting ", rs, " (chr", chr, ") - number #", i))

	# Save a temporary file with rs number
	rs_file <- paste0("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/09_CNVR/continuous/white_british/mirror/data/temp/", rs)
	fwrite(data.frame(ID = rs), rs_file, col.names = F, row.names = F, quote = F, sep = "\t")
	
	# Run PLINK: detect probes in LD
	input <- paste0("/data/FAC/FBM/DBC/zkutalik/default_sensitive/cauwerx/CNV_GWA/mirror/ukb_cnv_bed/ukb_cnv_chr", chr)
	system(paste("/home/cauwerx/scratch/cauwerx/softwares/plink/plink",
				 "--bfile", input,
				 "--keep /home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/01_samples/white_british/data/final/samples_white_british_plink1.9_All.txt",
				 "--show-tags", rs_file,
				 "--tag-kb 3000 --tag-r2 0.5",
				 "--out", rs_file))
 	unlink(rs_file); unlink(paste0(rs_file, ".log"))

	# Detect the boundaries of the CNVR
	cnvr <- as.data.frame(fread(paste0(rs_file, ".tags"), header = F, col.names = "ID"))
	cnvr <- left_join(cnvr, probes, by = "ID")
	df[which(df$ID == rs), "START_CNVR"] <- min(cnvr$POS)
	df[which(df$ID == rs), "END_CNVR"] <- max(cnvr$POS)
	df[which(df$ID == rs), "LENGTH_CNVR"] <- max(cnvr$POS)-min(cnvr$POS)
	df[which(df$ID == rs), "NUM_PROBES_CNVR"] <- nrow(cnvr)
	print(paste0("Length of associated CNVR: ", max(cnvr$POS)-min(cnvr$POS), " kb (", nrow(cnvr), " probes)"))

}


#################################################
### STEP 4: Final remodeling and save ###########

# Order by p-value
df <- df[order(df$P), ]

# Replace chromosome "XY" by "X"
df$CHR <- sub("XY", "X", df$CHR)

# Save
fwrite(df, "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/09_CNVR/continuous/white_british/mirror/data/final/GWAS_CNVR_mirror.txt", col.names = T, row.names = F, quote = F, sep = "\t")


#################################################
### STEP 5: ANNOVAR input & Save ################

# Order for ANNOVAR input 
df_ANNOVAR <- data.frame(CHR = df$CHR, START = df$START_CNVR, END = df$END_CNVR, REF = 0, OBS = "-", PHENO = df$PHENO, ID = df$ID, BETA = df$BETA, P = df$P)

# Save
fwrite(df_ANNOVAR, "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/09_CNVR/continuous/white_british/mirror/data/final/GWAS_CNVR_mirror_ANNOVAR_input.txt", col.names = F, row.names = F, quote = F, sep = "\t")
