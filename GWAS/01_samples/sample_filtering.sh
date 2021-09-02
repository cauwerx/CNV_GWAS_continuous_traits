#!/bin/bash

#################
#   RUN INFO    #
#################

#SBATCH --job-name=sample_filter   	# Give your job a name to recognize it in queue overview
#SBATCH --nodes=1          			# Define number of nodes required
#SBATCH --time=0-00:15:00   		# Define how long the job will run (max. is 7 days, default is 1h)
#SBATCH --partition=normal  		# Define the partition on which teh job shall run. May be omitted
#SBATCH --mem=15GB           		# Memory required per node

#SBATCH --output=/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/01_samples/white_british/data/log/sample_filtering-%j.out

#################
#   JOB INFO    #
#################

Rscript /home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/01_samples/white_british/script/sample_filtering.R

echo "Job done"

