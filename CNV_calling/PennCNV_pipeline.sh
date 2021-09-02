#!/bin/bash

#################
#   RUN INFO    #
#################

##### Job ##############
#SBATCH --job-name=PennCNV_pipeline         	# Give your job a name to recognize it in queue overview
#SBATCH --nodes=1                           	# Define number of nodes required
#SBATCH --time=0-05:00:00                   	# Define how long the job will run (max. is 7 days, default 1h)
#SBATCH --partition=normal                  	# Define the partition on which the job runs. May be omitted
#SBATCH --mem=5MB                              	# Memory required per node
#SBATCH --output=UKBB_500K_CNV/log-PennCNV_pipeline_Part6-%j.out   	# Name of log file


###########################################
#      STEP 1 : split BAF/LRR files       #
###########################################
# On the UKBB portal, there is one BAF and one LRR file per chromosome, which contain values for all individuals in the biobank.
# The first step splits these files by genotyping batch into smaller files.

# BAF
jid1_1=$(sbatch --parsable --output=UKBB_500K_CNV/intensity_files/log/01_1_split_baf_by_batch/log-01_1_split_baf_by_batch-%A-%a.out UKBB_500K_CNV/scripts/01_1_split_baf_by_batch_array.sh)
echo "Job ID BAF split: ${jid1_1}"
	
# LRR
jid1_2=$(sbatch --parsable --output=UKBB_500K_CNV/intensity_files/log/01_2_split_l2r_by_batch/log-01_2_split_l2r_by_batch-%A-%a.out UKBB_500K_CNV/scripts/01_2_split_l2r_by_batch_array.sh)
echo "Job ID LRR split: ${jid1_2}"	
	

###########################################
#      STEP 2 : Merge chromosomes         #
###########################################
# For each genotyping batch, chromosme-wide BAF and LRR files are merged into single genome-wide BAF and LRR files.

# BAF
jid2_1=$(sbatch --parsable --dependency=afterok:${jid1_1} --output=UKBB_500K_CNV/intensity_files/log/02_1_merge_baf/log-02_1_merge_baf-%A-%a.out UKBB_500K_CNV/scripts/02_1_merge_baf_array.sh)
echo "Job ID BAF chr merge: ${jid2_1}"		

# LRR
jid2_2=$(sbatch --parsable --dependency=afterok:${jid1_2} --output=UKBB_500K_CNV/intensity_files/log/02_2_merge_l2r/log-02_2_merge_l2r-%A-%a.out UKBB_500K_CNV/scripts/02_2_merge_l2r_array.sh)
echo "Job ID LRR chr merge: ${jid2_2}"


###########################################
#      STEP 3 : generate PFB files        #
###########################################
# Generate batch-wise population frequency of B allele files

# PFB files
jid3=$(sbatch --parsable --dependency=afterok:${jid2_1} --output=UKBB_500K_CNV/pfb_files/log/log-03_PFB-%A-%a.out UKBB_500K_CNV/scripts/03_PFB_array.sh)
echo "Job ID PFB generation: ${jid3}"	


###########################################
#     STEP 4 : split intensity files      #
###########################################
# Use output of STEP 1 and 2 to create individual, genome-wide intensity files (= BAF + LRR)

# Individual intensity files
jid4=$(sbatch --parsable --dependency=afterok:${jid2_1}:${jid2_2} --output=UKBB_500K_CNV/intensity_files/log/04_split_eid/log-04_split_eid-%A-%a.out UKBB_500K_CNV/scripts/04_split_eid_array.sh) 
echo "Job ID eid split: ${jid4}"
	

###########################################
#           STEP 5 : PennCNV              # 
###########################################
# Run PennCNV on individual intensity files to call CNVs on autosomes and chromosome X

# PennCNV - Autosomes
jid5_1=$(sbatch --parsable --dependency=afterok:${jid3}:${jid4} --output=UKBB_500K_CNV/PennCNV/log/Autosomes/log-05_1_PennCNV-%A-%a.out UKBB_500K_CNV/scripts/05_1_PennCNV_array.sh) 
echo "Job ID PennCNV autosomes: ${jid5_1}"	

# PennCNV - ChrX
jid5_2=$(sbatch --parsable --dependency=afterok:${jid5_1} --output=UKBB_500K_CNV/PennCNV/log/ChrX/log-05_2_PennCNV_ChrX-%A-%a.out UKBB_500K_CNV/scripts/05_2_PennCNV_ChrX_array.sh) 
echo "Job ID PennCNV ChrX: ${jid5_2}"


###########################################
#        STEP 6 : Quality Score           #
###########################################
# CNV quality score pipeline

# Quality score
jid6=$(sbatch --parsable --dependency=afterok:${jid5_1} --output=UKBB_500K_CNV/QS/log/log-06_QS-%A-%a.out UKBB_500K_CNV/scripts/06_QS_array.sh)
echo "Job ID QS: ${jid6}"

echo "Job done"    
