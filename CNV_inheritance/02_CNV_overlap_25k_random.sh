#!/bin/bash

#################
#   RUN INFO    #
#################

#SBATCH --job-name=CNV_random   	# Give your job a name to recognize it in queue overview
#SBATCH --nodes=1          			# Define number of nodes required
#SBATCH --time=0-30:00:00   		# Define how long the job will run (max. is 7 days, default is 1h)
#SBATCH --partition=normal  		# Define the partition on which teh job shall run. May be omitted
#SBATCH --mem=2GB           		# Memory required per node
#SBATCH --output=/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/frequency_analysis/02_CNV_inheritance/white_british/data/log/02_CNV_overlap_25k_random-%j.out

#################
#   JOB INFO    #
#################

Rscript /home/cauwerx/scratch/cauwerx/projects/cnv_gwas/frequency_analysis/02_CNV_inheritance/white_british/script/02_CNV_overlap_25k_random.R

echo "Job done"
