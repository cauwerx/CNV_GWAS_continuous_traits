# Create a summary of the effect of CNV frequency filter and pruning at different thresholds

#################################################
### Libraries ###################################
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)


#################################################
### Create an empty dataframe ###################

probe_sum <- data.frame(Chr = c(seq(1:22), "X", "XY", "All"))


#################################################
### STEP 1: Total number of probes ##############

# Load all probe identities and correct chromosome name for the PAR; "probes.txt.gz" contains the rs numberr, chromosome, and position for each genotyped probe
all_probes <- data.frame(fread("/home/cauwerx/scratch/cauwerx/projects/general_data/UKBB/probes.txt.gz", header = T))
all_probes$Chr <- gsub(23, "XY", all_probes$Chr)

# Count the number of probes per chromosome
all_probes_ct <- all_probes %>% count(Chr)

# Merge to the empty dataframe
probe_sum <- left_join(probe_sum, all_probes_ct, by = "Chr")
colnames(probe_sum)[2] <- "NumProbes"


#################################################
### STEP 2: CNV, pruned & Neff probes ###########

# CNV frequency > 0.005% (QS > 0.5) and missing in < 5% samples

for (c in as.character(probe_sum$Chr[1:24])) {

	# CNV frequency > 0.005% (QS > 0.5) and missing in < 5% samples
	probe_sum[which(probe_sum$Chr == c), "NumCNVProbes"] <- nrow(as.data.frame(fread(paste0("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/02_probes/white_british/stringent/mirror/data/temp/probes_WB_All_mirror_0.005_chr", c,".txt"), header = F))) # Temporary output from GWAS/02_probes/mirror/probe_pruning_mirror.R

	# Retained after pruning (r < 0.9)
	probe_sum[which(probe_sum$Chr == c), "NumPruneProbes_09"] <- as.numeric(nrow(as.data.frame(fread(paste0("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/02_probes/white_british/stringent/mirror/data/final/pruning/r_0.9/probes_WB_All_prune_0.9_mirror_chr", c,".prune.in"), header = F)))) # Final output from GWAS/02_probes/mirror/probe_pruning_mirror.R

	# Retained after pruning (r < 0.99)
	probe_sum[which(probe_sum$Chr == c), "NumPruneProbes_099"] <- as.numeric(nrow(as.data.frame(fread(paste0("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/02_probes/white_british/stringent/mirror/data/final/pruning/r_0.99/probes_WB_All_prune_0.99_mirror_chr", c,".prune.in"), header = F)))) # Final output from GWAS/02_probes/mirror/probe_pruning_mirror.R

	# Retained after pruning (r < 0.999)
	probe_sum[which(probe_sum$Chr == c), "NumPruneProbes_0999"] <- as.numeric(nrow(as.data.frame(fread(paste0("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/02_probes/white_british/stringent/mirror/data/final/pruning/r_0.999/probes_WB_All_prune_0.999_mirror_chr", c,".prune.in"), header = F)))) # Final output from GWAS/02_probes/mirror/probe_pruning_mirror.R

	# Retained after pruning (r < 0.9999)
	probe_sum[which(probe_sum$Chr == c), "NumPruneProbes_09999"] <- as.numeric(nrow(as.data.frame(fread(paste0("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/02_probes/white_british/stringent/mirror/data/final/pruning/r_0.9999/probes_WB_All_prune_0.9999_mirror_chr", c,".prune.in"), header = F)))) # Final output from GWAS/02_probes/mirror/probe_pruning_mirror.R

	# Chromosome Neff
	probe_sum[which(probe_sum$Chr == c), "Neff"] <- as.numeric(as.data.frame(fread(paste0("/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/02_probes/white_british/stringent/mirror/data/final/Neff/Neff_chr", c,".txt"), header = F))[1,1]) # Output from GWAS/02_probes/mirror/Neff/chr_Neff_mirror.R

}

#################################################
### STEP 4: Column Sum ##########################

for (c in 2:8) {
	probe_sum[25, c] <- sum(probe_sum[1:24, c])
}


#################################################
### Save ########################################

fwrite(probe_sum, "/home/cauwerx/scratch/cauwerx/projects/cnv_gwas/gwas/02_probes/white_british/stringent/mirror/data/final/global_summary.txt", col.names = T, row.names = F, quote = F, sep = "\t")

