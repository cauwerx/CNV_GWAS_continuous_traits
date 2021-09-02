#!/bin/bash

#################
#   RUN INFO    #
#################

#SBATCH --job-name=gwas_correct_dup  # Give your job a name to recognize it in queue overview
#SBATCH --nodes=1          				# Define number of nodes required
#SBATCH --time=0-01:00:00   			# Define how long the job will run (max. is 7 days, default is 1h)
#SBATCH --partition=normal  			# Define the partition on which teh job shall run. May be omitted
#SBATCH --mem=1GB           			# Memory required per node

# Parallelize over sex group
#SBATCH --array=1-3 
#SBATCH --output=/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/06_gwas/life_history/white_british/duplication_only/data/log/02_gwas_PD_correct_dup-%A_%a.out


#################
#   VARIABLES   #
#################

echo "SLURM_JOB_ID: " $SLURM_JOB_ID
echo "SLURM_ARRAY_ID: " $SLURM_ARRAY_ID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID


#################
#   JOB INFO    #
#################

sex=$(echo M F All)
s=$(echo ${sex} | cut -f ${SLURM_ARRAY_TASK_ID} -d ' ')


##################
#    RUN JOB     #
##################

Rscript /home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/06_gwas/life_history/white_british/duplication_only/script/02_correct_gwas_LH_dup.R ${s}

echo "Job Done!"
