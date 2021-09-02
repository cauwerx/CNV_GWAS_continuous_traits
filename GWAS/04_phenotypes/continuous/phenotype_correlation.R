# Calculate the correlation across continuous phenotypes

#################################################
### Libraries ###################################
library(data.table)
library(dplyr)


#################################################
### Load Data ###################################

# These files are the output of GWAS/04_phenotypes/continuous/phenotype_extraction.R
all <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/05_phenotypes/continuous/white_british/data/final/raw/pheno_continuous_WB_raw_All.txt", header = T))
males <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/05_phenotypes/continuous/white_british/data/final/raw/pheno_continuous_WB_raw_M.txt", header = T))
females <- as.data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/05_phenotypes/continuous/white_british/data/final/raw/pheno_continuous_WB_raw_F.txt", header = T))


#################################################
### Complete with sex-specific traits ###########

all <- left_join(all, males[, c("IID", "facial_hair", "balding"), ], by = "IID")
all <- left_join(all, females[, c("IID", "menarche", "birth_weight_first_child", "menopause"), ], by = "IID")


#################################################
### Calculate the correlations ##################

cor_all <- as.data.frame(cor(all[, -c(1)], use = "pairwise.complete.obs", method = "pearson"))


#################################################
### Save ########################################

fwrite(cor_all, "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/05_phenotypes/continuous/white_british/data/final/raw/correlation/cor_pheno_continuous_WB_raw_All.txt", col.names = T, row.names = T, quote = F, sep = "\t")
