#!/bin/bash

#################
#   RUN INFO    #
#################

#SBATCH --job-name=locus_zoom           # Give your job a name to recognize it in queue overview
#SBATCH --nodes=1          				# Define number of nodes required
#SBATCH --time=0-00:10:00   			# Define how long the job will run (max. is 7 days, default is 1h)
#SBATCH --partition=normal  			# Define the partition on which teh job shall run. May be omitted
#SBATCH --mem=2GB           			# Memory required per node
#SBATCH --output=/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/99_candidates/1p36.11_hematological/data/log/01_locus_zoom-%j.out


##################
#   VARIABLES    #
##################

chr=1                                                   # Chromosome
start=25000000                                          # Start position 
end=26000000                                            # End posittion
pheno=$(echo reticulocyte_count-platelet_count-HbA1c)   # Associated phenotypes


##################
#    RUN JOB     #
##################

Rscript /home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/99_candidates/1p36.11_hematological/script/01_locus_zoom_del.R ${chr} ${start} ${end} ${pheno}

echo "Job Done!"
