#!/bin/bash

#################
#   RUN INFO    #
#################

#SBATCH --job-name=gwas_del      # Give your job a name to recognize it in queue overview
#SBATCH --nodes=1          			# Define number of nodes required
#SBATCH --time=0-00:20:00   		# Define how long the job will run (max. is 7 days, default is 1h)
#SBATCH --partition=normal  		# Define the partition on which teh job shall run. May be omitted
#SBATCH --mem=5GB           		# Memory required per node

# Parallelize by chromsome
#SBATCH --array=1-24 
#SBATCH --output=/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/06_gwas/life_history/white_british/deletion_only/data/log/01_gwas_PD_clumping_del-%A_%a.out


#################
#   VARIABLES   #
#################

echo "SLURM_JOB_ID: " $SLURM_JOB_ID
echo "SLURM_ARRAY_ID: " $SLURM_ARRAY_ID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID


#################
#   JOB INFO    #
#################

chromosomes=$(echo 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 XY X)
chr=$(echo ${chromosomes} | cut -f ${SLURM_ARRAY_TASK_ID} -d ' ')


##################
#  Part I: GWAS  #
##################

## DEFINE VARIABLES #######################

# PLINK file set - deletion-only model
input=$(echo "/data/FAC/FBM/DBC/zkutalik/default_sensitive/cauwerx/CNV_GWA/deletion_only/ukb_cnv_bed/ukb_cnv_chr${chr}")
# For analysis in All and Males, a fake .fam is used, with only female individuals, to avoid het. haploid problems on chromosome X in males

# INT covariate-corrected phenotypes (from GWAS/04_phenotypes/life_history/life_history_extraction.R)
pheno_M=$(echo "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/05_phenotypes/life_history/white_british/data/final/INT_age_age2_batch_PCs/life_history_WB_INT_age_age2_batch_PCs_M.txt")
pheno_F=$(echo "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/05_phenotypes/life_history/white_british/data/final/INT_age_age2_batch_PCs/life_history_WB_INT_age_age2_batch_PCs_F.txt")
pheno_all=$(echo "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/05_phenotypes/life_history/white_british/data/final/INT_age_age2_batch_PCs/life_history_WB_INT_age_age2_sex_batch_PCs_All.txt")

# Unrelated, white British samples (from GWAS/01_samples/sample_filtering.R)
samples_M=$(echo "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/01_samples/white_british/data/final/samples_white_british_M.txt")
samples_F=$(echo "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/01_samples/white_british/data/final/samples_white_british_F.txt")
samples_all=$(echo "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/01_samples/white_british/data/final/samples_white_british_All.txt")

# Probes (from GWAS/02_probes/deletion_only/probe_pruning_del.R; r^2 at 0.9999)
probes=$(echo "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/02_probes/white_british/stringent/deletion_only/data/final/pruning/r_0.9999/probes_WB_All_prune_0.9999_del_chr${chr}.prune.in")

# Output
output_gwas_M=$(echo "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/06_gwas/life_history/white_british/deletion_only/data/temp/M/associations/gwas_PD_del_M_chr${chr}")
output_gwas_F=$(echo "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/06_gwas/life_history/white_british/deletion_only/data/temp/F/associations/gwas_PD_del_F_chr${chr}")
output_gwas_all=$(echo "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/06_gwas/life_history/white_british/deletion_only/data/temp/All/associations/gwas_PD_del_All_chr${chr}")


## JOB: GWAS MALES ########################
echo "Starting GWAS males"
/home/cauwerx/scratch/cauwerx/softwares/plink2/plink2 	--bfile ${input} --fam /data/FAC/FBM/DBC/zkutalik/default_sensitive/cauwerx/CNV_GWA/mirror/ukb_cnv_bed/ukb_cnv_chrX_onlyF.fam \
														--pheno ${pheno_M} --no-psam-pheno \
														--glm omit-ref no-x-sex hide-covar allow-no-covars --ci 0.95 \
														--keep ${samples_M} \
														--extract ${probes} \
														--out ${output_gwas_M}

## JOB: GWAS FEMALES ######################
echo "Starting GWAS females"
/home/cauwerx/scratch/cauwerx/softwares/plink2/plink2 	--bfile ${input} \
														--pheno ${pheno_F} --no-psam-pheno \
														--glm omit-ref no-x-sex hide-covar allow-no-covars --ci 0.95 \
														--keep ${samples_F} \
														--extract ${probes} \
														--out ${output_gwas_F}


## JOB: GWAS ALL ##########################
echo "Starting GWAS all"
/home/cauwerx/scratch/cauwerx/softwares/plink2/plink2 	--bfile ${input} --fam /data/FAC/FBM/DBC/zkutalik/default_sensitive/cauwerx/CNV_GWA/mirror/ukb_cnv_bed/ukb_cnv_chrX_onlyF.fam \
														--pheno ${pheno_all} --no-psam-pheno \
														--glm omit-ref no-x-sex hide-covar allow-no-covars --ci 0.95 \
														--keep ${samples_all} \
														--extract ${probes} \
														--out ${output_gwas_all}

wait

echo "Job done"
