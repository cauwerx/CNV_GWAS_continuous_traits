#!/bin/bash

#################
#   RUN INFO    #
#################

#SBATCH --job-name=ANNOVAR_GWAS_mirror	# Give your job a name to recognize it in queue overview
#SBATCH --nodes=1          				# Define number of nodes required
#SBATCH --time=0-00:30:00   			# Define how long the job will run (max. is 7 days, default is 1h)
#SBATCH --partition=normal  			# Define the partition on which teh job shall run. May be omitted
#SBATCH --mem=1GB           			# Memory required per node
#SBATCH --output=/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/10_ANNOVAR/continuous/white_british/mirror/data/log/02_region_annotation_GWASCatalog_mirror-%j.out


##################
#    RUN JOB     #
##################

Rscript /home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/10_ANNOVAR/continuous/white_british/mirror/script/02_region_annotation_GWASCatalog_mirror.R

echo "Job Done!"
