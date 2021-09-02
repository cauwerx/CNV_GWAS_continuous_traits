#!/bin/bash

#################
#   RUN INFO    #
#################

##### Job ##############
#SBATCH --job-name=PennCNV                  # Give your job a name to recognize it in queue overview
#SBATCH --nodes=1                           # Define number of nodes required
#SBATCH --time=7-00:00:00                   # Define how long the job will run (max. is 7 days, default 1h)
#SBATCH --partition=normal                  # Define the partition on which the job runs. May be omitted
#SBATCH --mem=3GB                           # Memory required per node

# Parallize by genotyping batch	
#SBATCH --array=1-106						


###################
# ARRAY VARIABLES # 
###################

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID


####################
#     OLD PERL     #
####################

# Load old Perl required by PennCNV
export PATH=$HOME/scratch/cauwerx/softwares/Perl-5.14.2/bin/:$PATH


####################
# DEFINE VARIABLES # 
####################

# Fixed varaibles: path to the HMM file and GC model used by PennCNV (downloaded, see methods)
path_to_hmm=/home/cauwerx/scratch/cauwerx/projects/CNV_GWA/UKBB_500K_CNV/hmm/affygw6.hmm
path_to_gc=/home/cauwerx/scratch/cauwerx/projects/CNV_GWA/UKBB_500K_CNV/GCmodel/GCmodel.txt

# Main variable
cd /home/cauwerx/scratch/cauwerx/projects/CNV_GWA/UKBB_500K_CNV/PennCNV/data
dir=$(ls -d $PWD/* | sort --version-sort | sed -n "${SLURM_ARRAY_TASK_ID}"p) # -d $PWD/* gives the full path

# Input variables: batch, path to and list of intensity files (from CNV_calling/scripts_for_pipeline/04_split_eid_arrray.sh), path to PFB file (from CNV_calling/scripts_for_pipeline/03_PFB_arrray.sh)
batch=$(echo "${dir}" | cut -d "/" -f11)
echo "Batch: ${batch}"

intensity_dir=$(echo "/home/cauwerx/scratch/cauwerx/projects/CNV_GWA/UKBB_500K_CNV/intensity_files/data/${batch}/eid_intensity_files/")
echo "Individual intensity files directory: ${intensity_dir}"

intensity_list=$(echo "/home/cauwerx/scratch/cauwerx/projects/CNV_GWA/UKBB_500K_CNV/intensity_files/data/${batch}/${batch}_intensities.txt")
echo "Individual intensity files list: ${intensity_list}"

path_to_pfb=$(echo "/home/cauwerx/scratch/cauwerx/projects/CNV_GWA/UKBB_500K_CNV/pfb_files/data/${batch}/${batch}_pfb.txt")
echo "PFB file: ${path_to_pfb}"

# Output variables: log file, raw CNV calls, merged CNV calls, CNV call QC, QC summary
log_file=$(echo "${batch}_log.txt")
output_file1=$(echo "${batch}_rawcnv.txt")
output_file2=$(echo "${batch}_merged_rawcnv.txt")
output_file3=$(echo "${batch}_qc_cnv.txt")
qc_summary=$(echo "${batch}_qc_summary.txt")


################
#     JOB      #
################

# Change to appropriated directory
cd ${dir}

# Run PennCNV to detect autosomal CNVs
echo "Starting CNV detection: ${batch}"
/scratch/beegfs/chuv_kuta/cauwerx/softwares/PennCNV-1.0.5/detect_cnv.pl -test -confidence --lastchr 23 -directory ${intensity_dir} -listfile ${intensity_list} -hmmfile ${path_to_hmm} -pfbfile ${path_to_pfb} -gcmodelfile ${path_to_gc} -logfile ${log_file} -output ${output_file1} 

# Run PennCNV to merge adjacent autosomal CNVs
echo "Starting CNV merge: ${batch}"
/scratch/beegfs/chuv_kuta/cauwerx/softwares/PennCNV-1.0.5/clean_cnv.pl combineseg ${output_file1} -signalfile ${path_to_pfb} -output ${output_file2} 

# Run PennCNV to perform QC on autosomal CNVs (default: -qclrrsd 0.3 -qcbafdrift 0.01 -qcwf 0.05)
echo "Starting CNV QC: ${batch}"
/scratch/beegfs/chuv_kuta/cauwerx/softwares/PennCNV-1.0.5/filter_cnv.pl ${output_file2} --chroms 1-23 -qclogfile ${log_file} -qcsumout ${qc_summary} -output ${output_file3} 


echo "Job done"
