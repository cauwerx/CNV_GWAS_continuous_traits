# Comment: This script takes as input cnv_data_global.rds from 02_CNV_QS.sh. 
# It converts the CNV level QS to probe level QS, generating a QS matrix. 
#
# Author: Aurelien Macé (Macé et al., 2016; Macé et al., 2017)
# Modified by: Chiara Auwerx (28/04/2020) --> modify read in/out function
###############################################################################


################################################################################################################################################################################################################################
#  Clean workspace
################################################################################################################################################################################################################################
rm(list = ls())
print('R - start to calculate the Quality Score')


################################################################################################################################################################################################################################
# Load libraries
################################################################################################################################################################################################################################
library(parallel)
library(plyr)


################################################################################################################################################################################################################################
# Define function
################################################################################################################################################################################################################################
probe_level_general <- function(dataset, probe_info, chr_list, sample_name, parameter_list, save_dir, nb_cores){
        probe_mat_cn    <- mclapply(chr_list, FUN = function(chr, d_dataset, d_probe_info, d_sample_name, d_parameter_list, d_save_dir){
        sub_probe           <- subset(d_probe_info, subset = Chr == chr)
        sub_probe           <- sub_probe[sort(sub_probe$Position, index.return = TRUE)$ix,]
        sub_dataset         <- subset(d_dataset , subset = Chromosome == chr, select = c('Sample_Name', 'Start_Position_bp', 'End_Position_bp', d_parameter_list))
        
        probe_level_para    <- matrix(0, nrow = nrow(sub_probe), ncol = length(d_sample_name), dimnames = list(sub_probe$Position, d_sample_name))

        probe_vector        <- sub_probe$Position
        
        data_mat_list_tp  <- ddply(sub_dataset, .(Sample_Name), .fun = function(data, d_probe_vector){
            probe_res           <- rep(0, length(probe_vector))
            names(probe_res)    <- probe_vector
            if(nrow(data) > 0){
                for(cnv in 1:nrow(data)){
                    idx_probe <- which(probe_vector >= data$Start_Position_bp[cnv] &probe_vector <= data$End_Position_bp[cnv])
                    if(length(idx_probe) > 0){
                        probe_res[idx_probe] <- data[cnv, d_parameter_list]
                    }
                }
            }
            return(probe_res)
        }, d_probe_vector = probe_vector)
        probe_level_para[, data_mat_list_tp$Sample_Name] <- t(data_mat_list_tp[, -1])

        save(probe_level_para, file = paste0(d_save_dir, "_probe_QS_chr", chr, ".rdata"))
        rm(probe_level_para)
        return(1)
    }, d_dataset = dataset, d_probe_info = probe_info, d_sample_name = sample_name, d_parameter_list = parameter_list, d_save_dir = save_dir, mc.preschedule = FALSE, mc.cores = nb_cores)
}

################################################################################################################################################################################################################################
#  Load the directory path information
################################################################################################################################################################################################################################
args <- commandArgs(TRUE)
tmp_dir_path <- args[1]             # Temporary directory
pfb_path <- args[2]                 # Path to PFB files
nb_cores <- as.numeric(args[3])     # Number of cores used


################################################################################################################################################################################################################################
#  Load the CNV data and create separate variable for deletions and duplications
################################################################################################################################################################################################################################
# "cnv_data_global.rdata" is the output of CNV_calling/scripts_for_pipeline/scripts_for_QS/R_script_calculate_quality_score.R with QS appended
load(file = paste(tmp_dir_path, 'cnv_data_global.rdata', sep = '_'))

    
################################################################################################################################################################################################################################
#  Load PFB files
################################################################################################################################################################################################################################
probe_info <- read.table(file = pfb_path, sep = '\t', header = TRUE)
probe_info$Chr <- sub("X", 24, probe_info$Chr)                        # Convert ChrX to Chr24 --> allow it to become numeric
probe_info$Chr  <- as.numeric(as.character(probe_info$Chr))


################################################################################################################################################################################################################################
#  Get the list of unique samples
################################################################################################################################################################################################################################
sample_name <- unique(cnv_data_global$Sample_Name)


################################################################################################################################################################################################################################
#  Calculate probe level parameters with the above described function
################################################################################################################################################################################################################################
res_cal <- probe_level_general(dataset = cnv_data_global, probe_info = probe_info, chr_list = 1:24, sample_name = sample_name, parameter_list = c('Quality_Score'), save_dir = tmp_dir_path, nb_cores = nb_cores) 


################################################################################################################################################################################################################################
#  Write text for the command windows
################################################################################################################################################################################################################################
print('R - finish to create the probe level matrice with the Quality Score')
