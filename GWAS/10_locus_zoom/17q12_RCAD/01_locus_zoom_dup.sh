#!/bin/bash

#################
#   RUN INFO    #
#################

#SBATCH --job-name=locus_zoom			# Give your job a name to recognize it in queue overview
#SBATCH --nodes=1          				# Define number of nodes required
#SBATCH --time=0-00:10:00   			# Define how long the job will run (max. is 7 days, default is 1h)
#SBATCH --partition=normal  			# Define the partition on which teh job shall run. May be omitted
#SBATCH --mem=2GB           			# Memory required per node
#SBATCH --output=/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/99_candidates/17q12_RCAD/data/log/01_locus_zoom-%j.out


##################
#   VARIABLES    #
##################

chr=17                                  # Chromosome
start=34000000                          # Start position
end=37000000                            # End position
pheno=$(echo CRP-Cr-cystatinC-urea)     # Associated phenotypes


##################
#    RUN JOB     #
##################

Rscript /home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/99_candidates/17q12_RCAD/script/01_locus_zoom_dup.R ${chr} ${start} ${end} ${pheno}

echo "Job Done!"
