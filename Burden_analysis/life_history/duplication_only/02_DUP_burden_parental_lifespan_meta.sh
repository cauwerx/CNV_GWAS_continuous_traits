#!/bin/bash

#################
#   RUN INFO    #
#################

#SBATCH --job-name=GWAMA_dup_PLS        # Give your job a name to recognize it in queue overview
#SBATCH --nodes=1          				# Define number of nodes required
#SBATCH --time=0-00:15:00   			# Define how long the job will run (max. is 7 days, default is 1h)
#SBATCH --partition=normal  			# Define the partition on which teh job shall run. May be omitted
#SBATCH --mem=2GB           			# Memory required per node
#SBATCH --output=/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/12_burden_analysis/life_history/white_british/duplication_only/data/log/02_DUP_burden_parental_lifespan_meta-%j.out

#################
#   JOB INFO    #
#################

Rscript /home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/12_burden_analysis/life_history/white_british/duplication_only/script/02_DUP_burden_parental_lifespan_meta.R

echo "Job done"
