#!/bin/bash

#################
#   RUN INFO    #
#################

##### Job ##############
#SBATCH --job-name=split_l2r                # Give your job a name to recognize it in queue overview
#SBATCH --nodes=1                           # Define number of nodes required
#SBATCH --time=5-00:00:00                   # Define how long the job will run (max. is 7 days, default 1h)
#SBATCH --partition=normal                  # Define the partition on which the job runs. May be omitted
#SBATCH --mem=1GB                           # Memory required per node

# Parallize by chromosome
#SBATCH --array=1-24						


###################
# ARRAY VARIABLES # 
###################

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

chr_list=$(echo 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 XY X)
chr=$(echo ${chr_list} | cut -d" " -f${SLURM_ARRAY_TASK_ID})


################
#  VARIABLES   #
################

# Path to the zipped LRR files; These are available from the UKBB portal
path_to_Zl2r="/data/chuv/sgg/data/UKBB/l2r/_001_ukb_l2r_chr${chr}_v2.txt.gz" 
echo "Used zipped LRR file: ${path_to_Zl2r}"

# Path to a file that contains all UKBB sample eids, except for the redacted ones
batch_info=/home/cauwerx/scratch/cauwerx/projects/CNV_GWA/UKBB_500K_CNV/eid_info/batch_info_no_redacted3.txt 


################
#     JOB      #
################

# Go to directory containing the intensity files
cd /scratch/beegfs/chuv_kuta/cauwerx/projects/CNV_GWA/UKBB_500K_CNV/intensity_files/data/

# Unzip the LRR file; exclude the 11 redacted3 samples; rename chrXY to chr23
if [ "${chr}" == "XY" ]
then
	unzip_chr=$(echo "${path_to_Zl2r:38:-12}23_temp.txt")
	zcat ${path_to_Zl2r} | awk '{$57327=$61154=$68759=$127547=$249705=$254293=$308492=$387485=$409933=$425326=$451667="";gsub(/ +/," ")}1' > ${unzip_chr}
else
	unzip_chr=$(echo "${path_to_Zl2r:38:-10}_temp.txt")
	zcat ${path_to_Zl2r} | awk '{$57327=$61154=$68759=$127547=$249705=$254293=$308492=$387485=$409933=$425326=$451667="";gsub(/ +/," ")}1' > ${unzip_chr}
fi

# Extract intensity data per genotyping batch (106 batches)
for ((i=2; i<=107; i++)); do 

	# Set parameters for batch i
	batch=$(awk 'NR=='"$i"'{print $1}' ${batch_info})
	start=$(awk 'NR=='"$i"'{print $3}' ${batch_info})
	end=$(awk 'NR=='"$i"'{print $4}' ${batch_info})
	output_file=$(echo "${batch}/${batch}_${unzip_chr::-9}.txt") 

	# Extract data frrom the main intensity file
	echo "Start extracting: ${batch}"
	awk -v b=${start} -v e=${end} '{for (c=b;c<e;c++) {printf("%s ", $c)} print $e}' ${unzip_chr} | tr ' ' '\t' > ${output_file};
done


################
#   CLEANING   #
################

# Delete the unzipped chromosome
rm ${unzip_chr}

echo "Job done"
