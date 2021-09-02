#!/bin/bash

#################
#   RUN INFO    #
#################

#SBATCH --job-name=CNVR        # Give your job a name to recognize it in queue overview
#SBATCH --nodes=1          		 # Define number of nodes required
#SBATCH --time=0-02:00:00        # Define how long the job will run (max. is 7 days, default is 1h)
#SBATCH --partition=normal  	 # Define the partition on which teh job shall run. May be omitted
#SBATCH --mem=1GB           	 # Memory required per node
#SBATCH --output=/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/09_CNVR/continuous/white_british/deletion_only/data/log/01_GWAS_CNVR_del-%j.out


##################
#    RUN JOB     #
##################

Rscript /home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/09_CNVR/continuous/white_british/deletion_only/script/prepare_input_del.R

echo "Job Done!"
