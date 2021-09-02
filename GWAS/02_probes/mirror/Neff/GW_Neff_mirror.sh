#!/bin/bash

#################
#   RUN INFO    #
#################

#SBATCH --job-name=GW_Neff   	    # Give your job a name to recognize it in queue overview
#SBATCH --nodes=1          			# Define number of nodes required
#SBATCH --time=0-00:02:00   		# Define how long the job will run (max. is 7 days, default is 1h)
#SBATCH --partition=normal  		# Define the partition on which teh job shall run. May be omitted
#SBATCH --mem=1GB           		# Memory required per node

#SBATCH --output=/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/02_probes/white_british/stringent/mirror/data/log/Neff/GW_Neff_mirror-%A_%a.out


#################
#   JOB INFO    #
#################

cd /home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/02_probes/white_british/stringent/mirror/data/final/Neff
paste Neff_chr* | awk '{sum=0; for(i=1; i<=NF; i++) sum += $i; print sum}' > /home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/02_probes/white_british/stringent/mirror/data/final/Neff/Neff_GW_mirror.txt

echo "Job done"

