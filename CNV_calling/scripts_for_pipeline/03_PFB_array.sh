#!/bin/bash

#################
#   RUN INFO    #
#################

##### Job ##############
#SBATCH --job-name=PFB                      # Give your job a name to recognize it in queue overview
#SBATCH --nodes=1                           # Define number of nodes required
#SBATCH --time=0-02:00:00                   # Define how long the job will run (max. is 7 days, default 1h)
#SBATCH --partition normal                  # Define the partition on which the job runs. May be omitted
#SBATCH --mem=1GB                           # Memory required per node

# Parallize by genotyping batch
#SBATCH --array=1-106 						


###################
# ARRAY VARIABLES # 
###################

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID


####################
# DEFINE VARIABLES # 
####################

# Path to file with probe informtaion; This file contains probe ID, chromosome, and genomic position
path_to_probe_info=$(echo "/home/cauwerx/scratch/cauwerx/projects/CNV_GWA/UKBB_500K_CNV/probe_info/probes.txt") 

# Main array variable: batch BAF file from CNV_calling/scripts_for_pipeline/02_1_merge_baf_arrray.sh 
cd /home/cauwerx/scratch/cauwerx/projects/CNV_GWA/UKBB_500K_CNV/intensity_files/data/
baf=$(find $(pwd) -name "*baf.txt"| sort --version-sort | sed -n "${SLURM_ARRAY_TASK_ID}"p) # $(pwd) gives full path
echo "Used BAF file: ${baf}"

# Input variable: batch
batch=$(echo "${baf}" | cut -d "/" -f11)
echo "Generate PFB for: ${batch}"

# Output variables: temporary and final batch-wise PFB file
temp_pfb_output=$(echo "/home/cauwerx/scratch/cauwerx/projects/CNV_GWA/UKBB_500K_CNV/pfb_files/data/${batch}/${batch}_pfb_temp.txt")
pfb_output=$(echo "/home/cauwerx/scratch/cauwerx/projects/CNV_GWA/UKBB_500K_CNV/pfb_files/data/${batch}/${batch}_pfb.txt")
	

################
#     JOB      #
################

# Calculate PFB; extract last line (mean); Name the column (PFB); temporarily save
awk '{sum=cnt=0; for (i=1;i<=NF;i++) if ($i != "NA") { sum+=$i; cnt++ } print $0, (cnt ? sum/cnt : "NA") }' ${baf} | awk '{print $NF}' | sed '1 i\PFB' > ${temp_pfb_output}

# Check that dimensions match
lines_probe=$(wc -l < "${path_to_probe_info}")
echo "Number of lines in probe file: ${lines_probe}"
lines_temp_pfb=$(wc -l < "${temp_pfb_output}")
echo "Number of lines in temp. PFB file: ${lines_temp_pfb}"

# Combine the probe information file and the PFB file
paste ${path_to_probe_info} ${temp_pfb_output} | grep -vw NA > ${pfb_output} 

# Check that dimensions match
lines_pfb=$(wc -l < "${pfb_output}")
echo "Number of lines in final PFB file: ${lines_pfb}"


################
#   CLEANING   #
################

# Remove temporary data
rm ${temp_pfb_output}

echo "Job done"
