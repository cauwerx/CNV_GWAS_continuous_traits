# CNV-GWAS on continuous traits
Code repository for *"The individual and global impact of copy number variants on complex human traits"*

https://doi.org/10.1101/2021.08.10.21261839

**Contact:** Chiara Auwerx (chiara.auwerx -at- unil.ch) or Zoltan Kutalik (zoltan.kutalik -at- unil.ch).


## Description of content: 

**Burden_analysis:** Contains the scripts to calculate each individual's autosomal CNV burden (`CNV_burden`), determine which individual carries a CNV overlapping CNV regions (CNVRs) identified through CNV-GWAS (`CNVR_covariates`). Resulting data is used to perform CNV burden analyses on continuous traits (`continuous`; both with and without correction for modifier CNVs) and life history traits (`life_history`). 
  

**CNV_calling:** Contains the pipeline for raw BAF and LRR data processing, CNV calling with PennCNV, quality score (QS) calculation, and transformation of QSs to the probe level. 


**CNV_inheritance:** Contains the scripts to calculate the fraction of CNVs individuals in our main analysis share with their UKBB siblings or random individuals.


**EstBB:** Contains the scripts to perform simulations aiming at establishing the power of our replication study of the CNV-GWAS in the Estonian Biobank. Results from the distinct association models are combined (`combined`).  


**GWAS:** Contains script for sample selection (`01_samples`), probe selection (`02_probes`), covariate (`03_covariates`) and phenotype (`04_phenotypes`) extraction and correction, GWAS with PLINK (`05_gwas`), identification of independent signals through stepwise conditional analysis (`06_SCA`), CNVR definition (`07_CNVR`) and annotation (`08_ANNOVAR`), and locus zoom on specific candidate CNVRs (`10_locus_zoom`).



## Data availability: 

CNV frequency and summary statistics available at:
