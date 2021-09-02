#!/bin/bash

#################
#   RUN INFO    #
#################

#SBATCH --job-name=pheno_LH   	# Give your job a name to recognize it in queue overview
#SBATCH --nodes=1          			# Define number of nodes required
#SBATCH --time=0-00:30:00   		# Define how long the job will run (max. is 7 days, default is 1h)
#SBATCH --partition=normal  		# Define the partition on which teh job shall run. May be omitted
#SBATCH --mem=10GB           		# Memory required per node
#SBATCH --output=/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/05_phenotypes/life_history/white_british/data/log/life_history_extraction-%j.out

#################
#   JOB INFO    #
#################

Rscript /home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/05_phenotypes/life_history/white_british/script/life_history_extraction.R

echo "Job done"

