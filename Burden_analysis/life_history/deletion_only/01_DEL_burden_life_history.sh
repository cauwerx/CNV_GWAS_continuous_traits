#!/bin/bash

#################
#   RUN INFO    #
#################

#SBATCH --job-name=LH_DEL_burden   	    # Give your job a name to recognize it in queue overview
#SBATCH --nodes=1          				# Define number of nodes required
#SBATCH --time=0-00:10:00   			# Define how long the job will run (max. is 7 days, default is 1h)
#SBATCH --partition=normal  			# Define the partition on which teh job shall run. May be omitted
#SBATCH --mem=1GB           			# Memory required per node
#SBATCH --output=/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/12_burden_analysis/life_history/white_british/deletion_only/data/log/01_DEL_burden_LH-%j.out

#################
#   JOB INFO    #
#################

Rscript /home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/12_burden_analysis/life_history/white_british/deletion_only/script/01_DEL_burden_life_history.R

echo "Job done"
