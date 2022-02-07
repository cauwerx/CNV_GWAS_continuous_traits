#!/bin/bash

#################
#   RUN INFO    #
#################

##### Job ##############
#SBATCH --job-name=QS                       # Give your job a name to recognize it in queue overview
#SBATCH --nodes=1                           # Define number of nodes required
#SBATCH --time=0-20:00:00                   # Define how long the job will run (max. is 7 days, default 1h)
#SBATCH --partition normal                  # Define the partition on which the job runs. May be omitted
#SBATCH --mem=70GB                          # Memory required per node

# Parallize by genotype batch
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

# Fixed variables: number of cores
nb_core=1

# Main variable
cd /home/cauwerx/scratch/cauwerx/projects/CNV_GWA/UKBB_500K_CNV/QS/data
dir=$(ls -d $PWD/* | sort --version-sort | sed -n "${SLURM_ARRAY_TASK_ID}"p) # -d $PWD/* gives the full path

# Input variables: batch, PennCNV outputs, PFB files
batch=$(echo "${dir}" | cut -d "/" -f11)
echo "Batch: ${batch}"

path_to_PennCNV=$(echo "/home/cauwerx/scratch/cauwerx/projects/CNV_GWA/UKBB_500K_CNV/PennCNV/data/${batch}")
echo "Path to PennCNV results: ${path_to_PennCNV}"

merged_rawcnv=$(echo "/home/cauwerx/scratch/cauwerx/projects/CNV_GWA/UKBB_500K_CNV/PennCNV/data/${batch}/${batch}_CNVs_AX.txt")
echo "Merged CNV input file: ${merged_rawcnv}"

qc_summary=$(echo "/home/cauwerx/scratch/cauwerx/projects/CNV_GWA/UKBB_500K_CNV/PennCNV/data/${batch}/${batch}_SQC_AX.txt")
echo "PennCNV QS summary input file: ${qc_summary}"

path_to_pfb=$(echo "/home/cauwerx/scratch/cauwerx/projects/CNV_GWA/UKBB_500K_CNV/pfb_files/data/${batch}/${batch}_pfb.txt") 
echo "PFB file: ${path_to_pfb}"

# Output variable
output_path=$(echo "/home/cauwerx/scratch/cauwerx/projects/CNV_GWA/UKBB_500K_CNV/QS/data/${batch}/${batch}")


################
#     JOB      #
################
# NOTE: scripts run here are at CNV_calling/scripts_for_pipeline/scripts_for_QS/...

# Merge PennCNV results from autosomes and ChrX
echo "Starting to combine autosomes/chrX CNVs: ${batch}"
Rscript /home/cauwerx/scratch/cauwerx/softwares/CNV_QS/R/R_script_merge_autosomes_chrX.R ${path_to_PennCNV} ${batch} 

# Change to appropriated directory
cd ${dir} 

# Quality score pipeline: Convert PennCNV data to R format
echo "Starting to convert: ${batch}"
Rscript /home/cauwerx/scratch/cauwerx/softwares/CNV_QS/R/R_script_convert_raw_cnv.R ${merged_rawcnv} ${qc_summary} ${output_path} ${nb_core} 

# Quality score pipeline: Calculate CNV quality score
echo "Starting CNV merge: ${batch}"
Rscript /home/cauwerx/scratch/cauwerx/softwares/CNV_QS/R/R_script_calculate_quality_score.R ${output_path} 

# Quality score pipeline: Calculate probe level QS
echo "Starting CNV merge: ${batch}"
Rscript /home/cauwerx/scratch/cauwerx/softwares/CNV_QS/R/R_script_create_probe_level_data.R ${output_path} ${path_to_pfb} ${nb_core} 

echo "Job done"
