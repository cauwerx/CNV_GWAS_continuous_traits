#!/bin/bash

#################
#   RUN INFO    #
#################

#SBATCH --job-name=geno_count   	# Give your job a name to recognize it in queue overview
#SBATCH --nodes=1          			# Define number of nodes required
#SBATCH --time=0-00:20:00   		# Define how long the job will run (max. is 7 days, default is 1h)
#SBATCH --partition=normal  		# Define the partition on which teh job shall run. May be omitted
#SBATCH --mem=1GB           		# Memory required per node

# Parallelize by chromosome

#SBATCH --array=1-24 
#SBATCH --output=/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/02_probes/white_british/stringent/frequency/data/log/01_genotype_count-%A_%a.out


#################
#   INTERNALS   #
#################

echo "SLURM_JOB_ID: " $SLURM_JOB_ID
echo "SLURM_ARRAY_ID: " $SLURM_ARRAY_ID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID


#################
#   VARIABLES   #
#################

# Define chromosome
chromosomes=$(echo 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 XY X)
chr=$(echo ${chromosomes} | cut -f ${SLURM_ARRAY_TASK_ID} -d ' ')

# Prepare input 
input=$(echo "/data/FAC/FBM/DBC/zkutalik/default_sensitive/cauwerx/CNV_GWA/mirror/ukb_cnv_bed/ukb_cnv_chr${chr}") # This corresponds to the .bim of the PLINK file set described in the methods of the paper (Table 1)
samples_all=$(echo "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/01_samples/white_british/data/final/samples_white_british_plink1.9_All.txt") # "samples_white_british_plink1.9_All.txt" are the results from GWAS/01_samples/sample_filtering.R
output_all=$(echo "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/02_probes/white_british/stringent/frequency/data/temp/genotype_count/All/genotype_count_chr${chr}_WB_all")


#################
#   JOB INFO    #
#################

# A fake .fam is used with only female individuals to avoid het. haploid prroblems on chromosome X in males 
# PLINK v1.9 
/home/cauwerx/scratch/cauwerx/softwares/plink/plink --bfile ${input} --fam /data/FAC/FBM/DBC/zkutalik/default_sensitive/cauwerx/CNV_GWA/mirror/ukb_cnv_bed/ukb_cnv_chrX_onlyF.fam \
													--keep ${samples_all} \
													--freqx gz \
													--out ${output_all}

# Clean
rm /home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/02_probes/white_british/stringent/frequency/data/temp/genotype_count/All/genotype_count_chr${chr}_WB_all.log
rm /home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/02_probes/white_british/stringent/frequency/data/temp/genotype_count/All/genotype_count_chr${chr}_WB_all.hh

echo "Job done"
