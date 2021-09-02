#!/bin/bash

#################
#   RUN INFO    #
#################

#SBATCH --job-name=power_EstBB_del	    # Give your job a name to recognize it in queue overview
#SBATCH --nodes=1          				# Define number of nodes required
#SBATCH --time=0-00:15:00   			# Define how long the job will run (max. is 7 days, default is 1h)
#SBATCH --partition=normal  			# Define the partition on which teh job shall run. May be omitted
#SBATCH --mem=1GB           			# Memory required per node

# Parallelize by replicated CNVR
#SBATCH --array=1-71 
#SBATCH --output=/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/13_EstBB/power_analysis/deletion_only/data/log/01_simmulate_p-val_del-%A_%a.out


#################
#   VARIABLES   #
#################

echo "SLURM_JOB_ID: " $SLURM_JOB_ID
echo "SLURM_ARRAY_ID: " $SLURM_ARRAY_ID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID


##################
#    RUN JOB     #
##################

Rscript /home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/13_EstBB/power_analysis/deletion_only/script/01_simulate_p-val_del.R ${SLURM_ARRAY_TASK_ID}

echo "Job Done!"
