#!/bin/bash

#################
#   RUN INFO    #
#################

##### Job ##############
#SBATCH --job-name=merge_l2r                # Give your job a name to recognize it in queue overview
#SBATCH --nodes=1                           # Define number of nodes required
#SBATCH --time=0-05:00:00                   # Define how long the job will run (max. is 7 days, default 1h)
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

# Main array variable
cd /home/cauwerx/scratch/cauwerx/projects/CNV_GWA/UKBB_500K_CNV/intensity_files/data/
dir=$(ls -d $PWD/* | sort --version-sort | sed -n "${SLURM_ARRAY_TASK_ID}"p) # -d $PWD/* gives the full path
cd ${dir}

# Input variables: batch and LRR files (from CNV_calling/scripts_for_pipeline/01_2_split_l2r_by_batch_array.sh)
batch=$(echo "${dir}" | cut -d "/" -f11)
echo "Merging LRR batch: ${batch}"

l2r=$(ls | grep "${batch}_l2r" | sort --version-sort)
echo "Files merged: ${l2r}"

# Output variables: concatenated LRR file 
l2r_output=$(echo "${batch}_l2r.txt")


################
#     JOB      #
################

# Concatenate LRR files
cat ${l2r} > ${l2r_output}


################
#   CLEANING   #
################

# Remove temporary data
rm ${l2r}

echo "Job done"
