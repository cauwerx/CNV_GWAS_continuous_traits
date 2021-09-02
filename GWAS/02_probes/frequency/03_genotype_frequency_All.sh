#!/bin/bash

#################
#   RUN INFO    #
#################

#SBATCH --job-name=geno_freq    	# Give your job a name to recognize it in queue overview
#SBATCH --nodes=1          			# Define number of nodes required
#SBATCH --time=0-00:20:00   		# Define how long the job will run (max. is 7 days, default is 1h)
#SBATCH --partition=normal  		# Define the partition on which teh job shall run. May be omitted
#SBATCH --mem=1GB           		# Memory required per node
#SBATCH --output=/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/02_probes/white_british/stringent/frequency/data/log/03_genotype_frequency_All-%j.out


#################
#   JOB INFO    #
#################

Rscript /home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/02_probes/white_british/stringent/frequency/script/03_genotype_frequency_All.R 

echo "Job done"
