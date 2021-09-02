#!/bin/bash

#################
#   RUN INFO    #
#################

##### Job ##############
#SBATCH --job-name=split_eid                # Give your job a name to recognize it in queue overview
#SBATCH --nodes=1                           # Define number of nodes required
#SBATCH --time=7-00:00:00                   # Define how long the job will run (max. is 7 days, default 1h)
#SBATCH --partition=normal                  # Define the partition on which the job runs. May be omitted
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
#     OLD PERL     #
####################

# Load old Perl required by PennCNV
export PATH=$HOME/scratch/cauwerx/softwares/Perl-5.14.2/bin/:$PATH


####################
# DEFINE VARIABLES # 
####################

# Main variable
cd /home/cauwerx/scratch/cauwerx/projects/CNV_GWA/UKBB_500K_CNV/intensity_files/data/
dir=$(ls -d $PWD/* | sort --version-sort | sed -n "${SLURM_ARRAY_TASK_ID}"p) # -d $PWD/* gives the full path
cd ${dir}


# Input variables: batch, probe IDs, sample IDS, batch BAF files (from CNV_calling/scripts_for_pipeline/02_1_merge_baf_arrray.sh), batch LRR files (from CNV_calling/scripts_for_pipeline/02_2_merge_l2r_arrray.sh)
batch=$(echo "${dir}" | cut -d "/" -f11)
echo "Batch: ${batch}"

probe_names=$(echo "/home/cauwerx/scratch/cauwerx/projects/CNV_GWA/UKBB_500K_CNV/intensity_files/data/${batch}/probe_names.txt")
cut -d$'\t' -f1 "/home/cauwerx/scratch/cauwerx/projects/CNV_GWA/UKBB_500K_CNV/probe_info/probes.txt" | sed '1 s/SNP/Name/' > ${probe_names} 

EIDLIST=$(grep -w "${batch}" /scratch/beegfs/chuv_kuta/cauwerx/projects/CNV_GWA/UKBB_500K_CNV/eid_info/batch_eid_no_redacted3.txt | cut -d$'\t' -f4) 
nb_eid=$(wc -w <<< "${EIDLIST}")
echo "Number of eids: ${nb_eid}"

baf_file=$(echo "${dir}/${batch}_baf.txt") 
echo "BAF file used: ${baf_file}"

l2r_file=$(echo "${dir}/${batch}_l2r.txt")
echo "LRR file used: ${l2r_file}"


################
#     JOB      #
################

# Go to the directory in which individual intensity files will be saved
cd ${dir}/eid_intensity_files

# Split batch BAF into individual BAF files with kcolumn.pl (from PennCNV) in batches of 1000 samples
echo "Start BAF split: 1-1000"
/home/cauwerx/scratch/cauwerx/softwares/PennCNV-1.0.5/kcolumn.pl $baf_file split 1 -heading 0 -start_split 1 -end_split 1000 -tab -out baf

echo "Start BAF split: 1001-2000"
/home/cauwerx/scratch/cauwerx/softwares/PennCNV-1.0.5/kcolumn.pl $baf_file split 1 -heading 0 -start_split 1001 -end_split 2000 -tab -out baf

echo "Start BAF split: 2001-3000"
/home/cauwerx/scratch/cauwerx/softwares/PennCNV-1.0.5/kcolumn.pl $baf_file split 1 -heading 0 -start_split 2001 -end_split 3000 -tab -out baf

echo "Start BAF split: 3001-4000"
/home/cauwerx/scratch/cauwerx/softwares/PennCNV-1.0.5/kcolumn.pl $baf_file split 1 -heading 0 -start_split 3001 -end_split 4000 -tab -out baf

echo "Start BAF split: 4001-end"
/home/cauwerx/scratch/cauwerx/softwares/PennCNV-1.0.5/kcolumn.pl $baf_file split 1 -heading 0 -start_split 4001 -tab -out baf


# Split batch LRR into individual LRR files with kcolumn.pl (from PennCNV) in batches of 1000 samples
echo "Start LRR split: 1-1000"
/home/cauwerx/scratch/cauwerx/softwares/PennCNV-1.0.5/kcolumn.pl $l2r_file split 1 -heading 0 -start_split 1 -end_split 1000 -tab -out l2r

echo "Start LRR split: 1001-2000"
/home/cauwerx/scratch/cauwerx/softwares/PennCNV-1.0.5/kcolumn.pl $l2r_file split 1 -heading 0 -start_split 1001 -end_split 2000 -tab -out l2r

echo "Start LRR split: 2001-3000"
/home/cauwerx/scratch/cauwerx/softwares/PennCNV-1.0.5/kcolumn.pl $l2r_file split 1 -heading 0 -start_split 2001 -end_split 3000 -tab -out l2r

echo "Start LRR split: 3001-4000"
/home/cauwerx/scratch/cauwerx/softwares/PennCNV-1.0.5/kcolumn.pl $l2r_file split 1 -heading 0 -start_split 3001 -end_split 4000 -tab -out l2r

echo "Start LRR split: 4001-end"
/home/cauwerx/scratch/cauwerx/softwares/PennCNV-1.0.5/kcolumn.pl $l2r_file split 1 -heading 0 -start_split 4001 -tab -out l2r


# Combine individual BAF and individual LRR files into an individual intensity file
# Loop over samples
for ((i=1; i<=${nb_eid}; i++)); do 

	# Prepare output file
	eid=$(echo ${EIDLIST} | cut -d" " -f"$i")
	output_file=$(echo "${eid}.txt") 

	# Define the BAF and LRR files to merge
	baf_temp=$(echo "baf.split${i}") 
	l2r_temp=$(echo "l2r.split${i}")
	
	# Add header to the BAF and LRR files
	sed -i '1 i\BAF' ${baf_temp} 
	sed -i '1 i\LRR' ${l2r_temp}

	# Control dimensions
	echo "Start Merging: ${eid}"
	lines_probe=$(wc -l < "${probe_names}")
	echo "Number of lines in probe file: ${lines_probe}"
	lines_baf=$(wc -l < "${baf_temp}")
	echo "Number of lines in temp. baf: ${lines_baf}"
	lines_l2r=$(wc -l < "${l2r_temp}")
	echo "Number of lines in temp. l2r: ${lines_l2r}"

	# Merging files and convert header to PennCNV format
	paste ${probe_names} ${baf_temp} ${l2r_temp} | sed "1 s/BAF/${eid}.B Allele Freq/" | sed "1 s/LRR/${eid}.Log R Ratio/" | grep -vw NA > ${output_file}
		
	# Delete temporary files
	rm ${baf_temp}
 	rm ${l2r_temp}

done;

# Write a file with the path of all the individual intensity files for eids in that given batch 
eid_intesnity_list_file=$(echo "${dir}/${batch}_intensities.txt")
ls > ${eid_intesnity_list_file}


################
#   CLEANING   #
################

# Delete the probe_names file
rm ${probe_names}

echo "Job done"
