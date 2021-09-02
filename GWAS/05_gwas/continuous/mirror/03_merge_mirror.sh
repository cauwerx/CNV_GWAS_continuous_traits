#!/bin/bash

#################
#   RUN INFO    #
#################

#SBATCH --job-name=gwas_merge_mirror  # Give your job a name to recognize it in queue overview
#SBATCH --nodes=1          				# Define number of nodes required
#SBATCH --time=0-00:30:00   			# Define how long the job will run (max. is 7 days, default is 1h)
#SBATCH --partition=normal  			# Define the partition on which teh job shall run. May be omitted
#SBATCH --mem=1GB           			# Memory required per node

# Parallelize over sex group
#SBATCH --array=1-3 
#SBATCH --output=/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/06_gwas/continuous/white_british/mirror/data/log/03_gwas_merge_mirror-%A_%a.out


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

Rscript /home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/06_gwas/continuous/white_british/mirror/script/03_merge_mirror.R ${s}

echo "Job Done!"
