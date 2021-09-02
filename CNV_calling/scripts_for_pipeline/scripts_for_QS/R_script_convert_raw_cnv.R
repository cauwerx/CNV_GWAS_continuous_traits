# Comment: This script takes as input the merged, raw CNV data and QC summary from PennCNV. 
# It combines information about the CNVs and eids into a single cnv_data_global.rdata, excluding eids with > 200 CNVs or a single CNV < 1e7 
# It additionally produces a log file indicating the number of excluded eids and CNVs and the total number of remaining deletions and duplications
# 
# Author: Aurelien Macé (Macé et al., 2016; Macé et al., 2017)
# Modified by: Chiara Auwerx (27/04/2020)
###############################################################################


################################################################################################################################################################################################################################
#  Clean workspace
################################################################################################################################################################################################################################
rm(list = ls())
print('R - start to convert raw CNV data')


################################################################################################################################################################################################################################
# Load libraries
################################################################################################################################################################################################################################
library(stringr)
library(parallel)


################################################################################################################################################################################################################################
#  Load the directory path information
################################################################################################################################################################################################################################
args <- commandArgs(TRUE)
cnv_data_path <- args[1]            # Autosomal + chrX merged, raw, CNV file from PennCNV ./clean_cnv.pl (e.g. "test_merged_rawcnv.txt"); from CNV_calling/scripts_for_pipeline/scripts_for_QS/R_script_merge_autosome_ChrX.R
sample_data_path <- args[2]         # Autosomal + chrX QS summary file from PennCNV ./filter_cnv.pl (e.g. "test_qc_summary.txt"); from CNV_calling/scripts_for_pipeline/scripts_for_QS/R_script_merge_autosome_ChrX.R
tmp_dir_path <- args[3]             # Path to store results ("./QS")
nb_cores  <- as.numeric(args[4])    # Number of cores used (e.g. 1)


################################################################################################################################################################################################################################
#  Load the CNV data
################################################################################################################################################################################################################################
cnv_data <- as.matrix(read.table(file = cnv_data_path, sep = '\t'))
sample_data <- read.table(file = sample_data_path, sep = '\t', header = TRUE,  stringsAsFactors = F)


################################################################################################################################################################################################################################
#  Convert the CNV data in a R format
################################################################################################################################################################################################################################
cnv_data   <- t(apply(cnv_data, 1, FUN =  function(CNV.line){
    #browser()
    CNV.lines.tp   <- CNV.line
    CNV.line       <- strsplit(CNV.lines.tp, ' ') 
    snp.idx        <- grep('numsnp', CNV.line)
    cn.idx         <- grep('state', CNV.line)
    id.idx         <- grep('.txt', CNV.line) 
    conf.idx       <- grep('conf', CNV.line)
    tmp            <- str_extract_all(CNV.line[1], "[0-9]+")
    tmp2           <- str_extract_all(CNV.line[cn.idx], "[0-9]+")
    tmp3           <- strsplit(CNV.line[[id.idx]], ".txt") 
    tmp4           <- str_extract_all(CNV.line[conf.idx], "[0-9]+")
    tmp5           <- str_extract_all(CNV.line[snp.idx], "[0-9]+")
    tmp6           <- str_extract_all(CNV.line[3], "[0-9]+")
    chr            <- as.integer(tmp[[1]][1])
    start          <- as.integer(tmp[[1]][2])
    end            <- as.integer(tmp[[1]][3])
    if(length(cn.idx) > 0){cn <- as.integer(tmp2[[1]][2])}else{cn  <- NA}
    if(length(conf.idx) > 0){conf <- as.integer(tmp4[[1]][1]) + as.integer(tmp4[[1]][2]) / 1000}else{conf  <- 0; print('No conf score')}
    
    id <- tmp3[[1]][length(tmp3[[1]])]

    if(length(snp.idx) > 0){numSNP  <- tmp5[[1]][1]}else{numSNP  <- NA}
    
    l.val   <- end - start + 1
    return(c('Chromosome' = chr, 'Start_Position_bp' = start, 'End_Position_bp' = end, 'Copy_Number' = cn, 'Max_Log_BF' = conf, 'Sample_Name' = id, 'No_Probes' = numSNP, 'Length_bp' = l.val))
}))

cnv_data    <- cnv_data[,c('Sample_Name', 'Chromosome', 'Start_Position_bp', 'End_Position_bp', 'Copy_Number', 'Max_Log_BF', 'Length_bp', 'No_Probes')]
cnv_data    <- as.data.frame(cnv_data, stringAsFactors = FALSE)

cnv_data$Sample_Name        <- as.character(cnv_data$Sample_Name)
cnv_data$Chromosome         <- as.numeric(as.character(cnv_data$Chromosome))
cnv_data$Start_Position_bp  <- as.numeric(as.character(cnv_data$Start_Position_bp))
cnv_data$End_Position_bp    <- as.numeric(as.character(cnv_data$End_Position_bp))
cnv_data$Copy_Number        <- as.numeric(as.character(cnv_data$Copy_Number))
cnv_data$Max_Log_BF         <- as.numeric(as.character(cnv_data$Max_Log_BF))
cnv_data$Length_bp          <- as.numeric(as.character(cnv_data$Length_bp))
cnv_data$No_Probes          <- as.numeric(as.character(cnv_data$No_Probes))

# remove CNV with NA confidence score
if(length(which(is.na(cnv_data$Max_Log_BF))) > 0){
    cnv_data<- cnv_data[-which(is.na(cnv_data$Max_Log_BF)),]
}


################################################################################################################################################################################################################################
#  Clean the Sample ID for the sample data
################################################################################################################################################################################################################################
sample_data$File <- sapply(as.character(sample_data$File), FUN = function(data){
                            val <- strsplit(data, '.txt')
                            return(val[[1]][length(val[[1]])])
                           })


################################################################################################################################################################################################################################
#  Merge CNV and sample data
################################################################################################################################################################################################################################
cnv_data_global <- merge(cnv_data, sample_data, by.x = 'Sample_Name', by.y = 'File')


################################################################################################################################################################################################################################
#  Count number of samples and CNV before cleaning
################################################################################################################################################################################################################################
log_data  <- c(nb_sample_before = length(unique(cnv_data_global$Sample_Name)), nb_cnv_before = nrow(cnv_data_global))


################################################################################################################################################################################################################################
# Keep the CNVs with NumCNV < 200 and Length < 1e7
################################################################################################################################################################################################################################
cnv_data_global <- subset(cnv_data_global, subset = NumCNV <= 200)
eid_to_rm <- unique(cnv_data_global[which(cnv_data_global$Length_bp >= 1e7), "Sample_Name"])
cnv_data_global <- cnv_data_global[!cnv_data_global$Sample_Name %in% eid_to_rm ,]


################################################################################################################################################################################################################################
#  Count number of samples and CNV after cleaning
################################################################################################################################################################################################################################
log_data <- c(log_data, nb_sample_after = length(unique(cnv_data_global$Sample_Name)), nb_cnv_after = nrow(cnv_data_global), nb_deletion = length(which(cnv_data_global$Copy_Number < 2)), nb_duplication = length(which(cnv_data_global$Copy_Number > 2)), eid_to_rm = eid_to_rm)


################################################################################################################################################################################################################################
#  Save the data
################################################################################################################################################################################################################################
save(cnv_data_global, file = paste(tmp_dir_path, 'cnv_data_global.rdata', sep = '_'))
write.table(log_data, file = paste(tmp_dir_path, 'log_info_raw.txt', sep = '_'), quote = FALSE, row.names = TRUE, col.names = FALSE)


################################################################################################################################################################################################################################
#  Write text for the command windows
################################################################################################################################################################################################################################
print('R - finish to convert raw CNV data')
