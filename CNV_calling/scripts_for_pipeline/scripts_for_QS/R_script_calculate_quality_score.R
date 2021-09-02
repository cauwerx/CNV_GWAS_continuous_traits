# Comment: This script takes as input cnv_data_global.rdata from 01_PennCNV_to_Rdata.sh. 
# It calculates a quality score using fixed parameters coeffecicients (Mace et al., 2016) for each CNV and appends it to cnv_data_global.rds
# In addition, it provides a QS summary file and histogram of the QS distribution
#
# Author: Aurelien Macé (Macé et al., 2016; Macé et al., 2017)
# Modified by: Chiara Auwerx (27/04/2020)
###############################################################################


################################################################################################################################################################################################################################
#  Clean workspace
################################################################################################################################################################################################################################
rm(list = ls())
print('R - start to calculate the Quality Score')


################################################################################################################################################################################################################################
# Define function
################################################################################################################################################################################################################################
calculate_QS <- function(cnv_data, model_data){
    QS_val <- rep(model_data['(Intercept)', 'Estimate'], nrow(cnv_data))
    for(para in 2:nrow(model_data)){
        para_name   <- rownames(model_data)[para]
        QS_val      <- QS_val + model_data[para_name, 'Estimate'] * scale(cnv_data[, para_name])
    }
    QS_val  <- 1 / (1 + exp(-QS_val))
    return(QS_val)
}

################################################################################################################################################################################################################################
#  Load the directory path information
################################################################################################################################################################################################################################
args <- commandArgs(TRUE)
tmp_dir_path <- args[1]     # Path to store results ("./QS")


################################################################################################################################################################################################################################
#  Load the CNV data and create separate variable for the deletions and duplications
################################################################################################################################################################################################################################
# "cnv_data_global.rdata" is the output of CNV_calling/scripts_for_pipeline/scripts_for_QS/R_script_convert_raw_cnv.R
load(file = paste(tmp_dir_path, 'cnv_data_global.rdata', sep = '_'))
cnv_data_global_del <- subset(cnv_data_global, subset = Copy_Number < 2)
cnv_data_global_dup <- subset(cnv_data_global, subset = Copy_Number > 2)


################################################################################################################################################################################################################################
#  define the QS coefficients for the deletion and duplication models
################################################################################################################################################################################################################################
QS_coefficient_duplication              <- cbind(c(-0.81529185, 3.90248885, -0.90184846, 0.34376765, -0.06351849, -0.01497046, -0.30665878, -0.10354156, -0.44829901, 0.06039380, 0.03645913), rep(0, 11))
rownames(QS_coefficient_duplication)    <- c('(Intercept)', 'Max_Log_BF', 'No_Probes', 'LRR_SD', 'BAF_drift', 'WF', 'LRR_mean', 'NumCNV', 'BAF_SD', 'Length_bp', 'BAF_mean')
colnames(QS_coefficient_duplication)    <- c('Estimate', 'Pr(>|z|)')

QS_coefficient_deletion             <- cbind(c(-1.75648377, 4.64898485, -2.50150285, -0.47552224, 0.28876320, 0.10205302, 0.14363692, 0.02959571, 0.00000000, -0.00000000), rep(0, 10))
rownames(QS_coefficient_deletion)   <- c('(Intercept)', 'Max_Log_BF', 'No_Probes', 'NumCNV', 'LRR_SD', 'Length_bp', 'LRR_mean', 'WF', 'BAF_drift', 'BAF_SD')
colnames(QS_coefficient_deletion)   <- c('Estimate', 'Pr(>|z|)')


################################################################################################################################################################################################################################
#  Calculate the QS for deletions and duplications
################################################################################################################################################################################################################################
QS_deletion <- calculate_QS(cnv_data = cnv_data_global_del, model_data = QS_coefficient_deletion)
QS_duplication <- calculate_QS(cnv_data = cnv_data_global_dup, model_data = QS_coefficient_duplication)


################################################################################################################################################################################################################################
#  Calculate the QS and append it to cnv_data_global
################################################################################################################################################################################################################################
QS <- matrix(NA, nrow = nrow(cnv_data_global), ncol = 1)
QS[which(cnv_data_global$Copy_Number < 2),] <- QS_deletion * -1
QS[which(cnv_data_global$Copy_Number > 2),] <- QS_duplication
colnames(QS) <- c('Quality_Score')
cnv_data_global  <- cbind(cnv_data_global, QS)


################################################################################################################################################################################################################################
#  Save the new CNV dataset with the QS 
################################################################################################################################################################################################################################

# batch_cnv_data_global.rdata
save(cnv_data_global, file = paste(tmp_dir_path, 'cnv_data_global.rdata', sep = '_'))
write.table(cnv_data_global, file = paste(tmp_dir_path, 'cnv_data_global.txt', sep = '_'), row.names= F, col.names = T, sep = "\t", quote = F)

# Summary file
write.table(summary(cnv_data_global)[,-c(1:4)], paste(tmp_dir_path, 'log_cnv_summary.txt', sep = '_'), col.names = T, sep = "\t", quote = F)

# QS distribution info
QS_histogram <- hist(cnv_data_global$Quality_Score, breaks = 20, plot = FALSE)
save(QS_histogram, file = paste(tmp_dir_path, 'log_QS_hist.rdata', sep = '_'))


################################################################################################################################################################################################################################
#  Write text for the command windows
################################################################################################################################################################################################################################
print('R - finish to calculate the Quality Score')
