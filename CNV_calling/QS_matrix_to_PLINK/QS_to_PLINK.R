########################################################
# Description
########################################################
# Date: August 16th 2022
# Author: Chiara Auwerx
# Publication: Auwerx et al. 2022 AJHG  (https://www.cell.com/ajhg/pdf/S0002-9297(22)00061-1.pdf). 
# Description: Example script to transform QS matrices (1 x file/chr/batch) into PLINK file sets (1 x file set/chr).
#  Notes:
# - This is an example script on simulated data that does not aim at recreating the true CNV landscape (STEP 0). It should be adapted for usage on real data. 
# - STEP 1 should be applied to each QS matrix, that is 1 x file per UKBB genotyping batch per chromosome (here exemplified for 1 file).
# - STEP 2 merges pre-.ped files from different batches but same chromosome; It results in 1 x .ped/.map file set per chromosome (here exemplified for 2 files from the same chromosome).
# - STEP 3 transforms chromosome-wise .ped/.map file sets into a .bed/.bim/.fam file sets (here exemplified for one .ped/.map file set).


########################################################
# Libraries & directories
########################################################
setwd("/Users/chiara/Desktop/QS_to_PLINK/") # Adapt path
library(dplyr)
library(data.table)
library(gtools)


########################################################
# STEP 0: Generate data - (if not using real data)
########################################################

# Define parameters (can be tuned)
n <- 1000 # n  corresponds to both the number of samples with at least one CNV and the number of probes.
n_no_cnv <- round(n/2, 0) # Number of samples with no CNVs (set as 1/3 of the number of samples with at least one CNV).
max_num_cnv <- 2 # Maximum number of CNVs per individual.
max_cnv_length <- 10 # Maximum number of probes covered by a CNV.
n_cnv_outlier <- round(n/20, 0) # Number of individuals that ought to be excluded (e.g. CNV outliers; set as 1/20 of the number with at least one CNV).


# Fake QS matrix of size probe x sample (n x n) with random QS scores from a uniform distribution (!Not meant to mimic real data structure well!).
# On real data, would be the outputed QS matrices from previous steps.
qs_matrix_1 <- matrix(data = 0, nrow = n, ncol = n)
for (i in 1:n){
  
  #  Number of CNVs (between 1 and 2).
  num_cnv <- sample(1:max_num_cnv, 1)
  for(j in 1:num_cnv) {
   
    # Determine the characteristics of the CNV (start position; length of the CNV in number of probes; QS of the CNV). 
    start <- sample(1:(n-max_cnv_length), 1)
    length <- sample(1:max_cnv_length, 1)
    QS <- runif(1, min = -1, max = 1)
    
    # Modify the empty n x n matrix with random QS scores.
    qs_matrix_1[start:(start+length), i] <- QS
  }
}; rm(i, j, num_cnv, start, length, QS)
rownames(qs_matrix_1) <- paste0("probe_", 1:n) 
colnames(qs_matrix_1) <- paste0("eid_", sample(1:(n + n_no_cnv), n, replace = F)) # The columns will be ordered randomly from n + n_no_cnv possible samples to mimic individuals with no CNVs that are not included in the QS matrix (Macé et al; Auwerx et al.). 


# Fake list will all individuals (n + n_no_cnv). On real data, this would correspond to the list of UKBB individuals (e.g. from UKBB .fam files).
all_samples <- paste0("eid_", 1:(n + n_no_cnv))


# Fake list with CNV outliers. On real data, this is generated from various filtering steps.
cnv_outlier_samples <- paste0("eid_", sample(1:(n + n_no_cnv), n_cnv_outlier, replace = F))


########################################################
# STEP 1: Transform the QS matrix into a .ped file
########################################################

# 1.1 Add samples with no CNVs.
samples_no_cnv <- setdiff(all_samples, colnames(qs_matrix_1))
samples_no_cnv_qs <- matrix(data = 0, nrow = n, ncol = length(samples_no_cnv))
colnames(samples_no_cnv_qs) <- samples_no_cnv
qs_matrix_1 <- cbind(qs_matrix_1, samples_no_cnv_qs); rm(samples_no_cnv_qs)
qs_matrix_1 <- qs_matrix_1[, all_samples]


# 1.2 Remove samples that are non-conform/outliers/not desired for follow-up analyses.
qs_matrix_1 <- qs_matrix_1[ , !colnames(qs_matrix_1) %in% cnv_outlier_samples] 


# 1.3 Apply QS filter (|QS| > 0.5) and encode CNVs as -1 (= deletion)/0 (= copy-neutral)/1 (= duplication). 
# Note that other QS filtering values can be used (e.g. |QS| > 0.1 for a less stringent filtering approach). 
qs_matrix_1[qs_matrix_1 <= -0.5] <- -1
qs_matrix_1[qs_matrix_1 > -0.5 & qs_matrix_1 < 0.5] <- 0
qs_matrix_1[qs_matrix_1 >= 0.5] <- 1


# 1.4 Remove probes that are not in a CNV state in any individual of this batch to reduce dimension of the data.  
ncnv_probes <- rownames(qs_matrix_1[which(rowSums(qs_matrix_1 == 0) == ncol(qs_matrix_1)), ])
qs_matrix_1 <- qs_matrix_1[!rownames(qs_matrix_1) %in% ncnv_probes, ]


# 1.5 Recode the CNV data for PLINK (here according to a mirror model; Other encodings are described in Auwerx et al., 2022 (see Table 1)).
qs_matrix_1[qs_matrix_1 == -1] <- "A A"
qs_matrix_1[qs_matrix_1 == 0] <- "A T"
qs_matrix_1[qs_matrix_1 == 1] <- "T T"


# 1.6 Transpose the matrix.
qs_matrix_1 <- t(qs_matrix_1)


# 1.7 Transform the matrix into a .ped file by adding information columns (Note: Sex information can be derived from UKBB files. Here it is generated randomly).
plink_ped_1 <- data.frame(FID = 0, IID = rownames(qs_matrix_1), PID  = 0, MID = 0, Sex = sample(1:2, nrow(qs_matrix_1), replace = T), Pheno = 0)
plink_ped_1 <- cbind(plink_ped_1, as.data.frame(qs_matrix_1))


########################################################
# STEP 2: Merge .ped files from different batches
########################################################

# 2.0 For the purpose of this example, we create a second .ped file "plink_ped_2" that does not contain the same probes than "plink_ped_1" to mimic situations where different probes are excluded at step 1.4.   
plink_ped_2 <- plink_ped_1
colnames(plink_ped_2)[sample(7:ncol(plink_ped_2), round(length(ncnv_probes)/2, 0), replace = F)] <- sample(ncnv_probes, round(length(ncnv_probes)/2, 0), replace = F)
plink_ped_2$IID <- paste0("eid_", (nrow(plink_ped_2)+1):(2*nrow(plink_ped_2)))


# 2.1 Merge .ped files recursively across batches (here, only for two) and set missing values to "A T" (= copy-neutral according to the mirror model).
plink_ped <- as.data.frame(rbind(setDT(plink_ped_1), setDT(plink_ped_2), fill = T))
plink_ped[is.na(plink_ped)] <- "A T"


# 2.2 Order probes, generate a .map file, and save the file. Note that on real data, this step will use probe information available for genotyped markers (e.g. from UKBB .bim file) that overlap the ones retained in "plink_ped".
probe_order <- mixedsort(colnames(plink_ped)[7:ncol(plink_ped)])
plink_map <- data.frame(CHR = 1, SNP = probe_order, GD = 0, POS = 1:length(probe_order))
fwrite(plink_map, "./data/fake_QS_to_PLINK.map", col.names = F, row.names = F, quote = F, sep = "\t")


# 2.3 Order variants according to the .map file, split each variant column into two columns (to mimic a .ped file), and save the file.
plink_ped <- plink_ped[, c(colnames(plink_ped)[1:6], probe_order)]
plink_ped <- cbind(plink_ped[,1:6], as.data.frame(unlist(lapply(plink_ped[, c(7:ncol(plink_ped))], data.table::tstrsplit, " "), recursive = FALSE)))
fwrite(plink_ped, "./data/fake_QS_to_PLINK.ped", col.names = F, row.names = F, quote = F, sep = "\t")


########################################################
# STEP 3: Transform .ped/.map into .bed/.bim/.fam files
########################################################

system("/Users/chiara/softwares/plink/plink --file ./data/fake_QS_to_PLINK --make-bed --out ./data/fake_QS_to_PLINK") # Adapt path to plink
