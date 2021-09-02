#!/bin/bash

#################
#   RUN INFO    #
#################

#SBATCH --job-name=power_EstBB_combo    # Give your job a name to recognize it in queue overview
#SBATCH --nodes=1          				# Define number of nodes required
#SBATCH --time=0-00:25:00   			# Define how long the job will run (max. is 7 days, default is 1h)
#SBATCH --partition=normal  			# Define the partition on which teh job shall run. May be omitted
#SBATCH --mem=1GB           			# Memory required per node
#SBATCH --output=/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/13_EstBB/power_analysis/combined/data/log/EstBB_replication_combined-%A_%a.out



##################
#    RUN JOB     #
##################

Rscript /home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/13_EstBB/power_analysis/mirror/script/01_simulate_p-val_mirror.R ${SLURM_ARRAY_TASK_ID}

echo "Job Done!"
