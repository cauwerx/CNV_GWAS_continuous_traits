#!/bin/bash

#################
#   RUN INFO    #
#################

#SBATCH --job-name=geno_miss    	# Give your job a name to recognize it in queue overview
#SBATCH --nodes=1          			# Define number of nodes required
#SBATCH --time=0-00:10:00   		# Define how long the job will run (max. is 7 days, default is 1h)
#SBATCH --partition=normal  		# Define the partition on which teh job shall run. May be omitted
#SBATCH --mem=1GB           		# Memory required per node
#SBATCH --output=/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/02_probes/white_british/stringent/frequency/data/log/02_high_missingness_0.95_All-%j.out


#################
#   JOB INFO    #
#################

Rscript /home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/02_probes/white_british/stringent/frequency/script/02_high_missingness_0.95.R 

echo "Job done"

