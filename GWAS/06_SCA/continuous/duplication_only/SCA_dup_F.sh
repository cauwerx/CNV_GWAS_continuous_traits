#!/bin/bash

#################
#   RUN INFO    #
#################

#SBATCH --job-name=SCA_dup_F        # Give your job a name to recognize it in queue overview
#SBATCH --nodes=1          			# Define number of nodes required
#SBATCH --time=0-03:00:00   		# Define how long the job will run (max. is 7 days, default is 1h)
#SBATCH --partition=normal  		# Define the partition on which teh job shall run. May be omitted
#SBATCH --mem=5GB           		# Memory required per node
#SBATCH --output=/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/08_SCA/continuous/white_british/duplication_only/data/log/SCA_dup_F-%j.out

#################
#   JOB INFO    #
#################

Rscript /home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/08_SCA/continuous/white_british/duplication_only/script/SCA_dup_F.R

echo "Job done"
