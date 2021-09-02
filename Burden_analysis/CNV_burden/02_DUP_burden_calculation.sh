#!/bin/bash

#################
#   RUN INFO    #
#################

#SBATCH --job-name=pheno_DUP_burden   	# Give your job a name to recognize it in queue overview
#SBATCH --nodes=1          				# Define number of nodes required
#SBATCH --time=0-00:20:00   			# Define how long the job will run (max. is 7 days, default is 1h)
#SBATCH --partition=normal  			# Define the partition on which teh job shall run. May be omitted
#SBATCH --mem=5GB           			# Memory required per node
#SBATCH --output=/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/05_phenotypes/CNV_burden/white_british/data/log/02_DUP_burden_calculation-%j.out

#################
#   JOB INFO    #
#################

Rscript /home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/05_phenotypes/CNV_burden/white_british/script/02_DUP_burden_calculation.R

echo "Job done"

