#=============================================================================#
# posthoc helper file for asp_ref
#
# set variables
# call scripts to generate components for rmarkdown report components
# (make sure you can call this script externally for tmux)
#=============================================================================#
library(tidyverse)
library(data.table)
library(stargazer)
library(lmtest)
library(broom)
library(figifs)
library(kableExtra)
rm(list = ls())

# variables and paths
exposure = 'calcium_dietqc2'
hrc_version = 'v2.3'
output_dir <- paste0("/media/work/gwis/posthoc/", exposure, "/")

covariates_set1 <- c("age_ref_imp", "sex", "energytot_imp", "pc1", "pc2", "pc3", "study_gxe")
covariates_list <- list(covariates_set1)
mod <- 'age_ref_imp+sex+energytot_imp+pc1+pc2+pc3+study_gxe'

#-----------------------------------------------#
# Compile suggestive and significant hits
# extract dosages from binarydosage
# needs 'exposure' and 'hrc_version' defined globally
#-----------------------------------------------#
source("/home/rak/Dropbox/FIGI/FIGI_code/results/posthoc/posthoc_01_combine_results.R")

#-----------------------------------------------#
# INPUT DATA
# rename dose variables to simply chr_bp_ref_alt
# keep original dosage + probabilities
#-----------------------------------------------#
figi <- posthoc_input(exposure, hrc_version)

# # also flip exposure for exposure
# figi <- figi %>% 
#   mutate({{exposure}} := 3 - .data[[exposure]]) %>% 
#   inner_join(geno, 'vcfid')

# GWIS findings
# filter out case only results. biased because strong population E|G association
snps <- readRDS(paste0("/media/work/gwis/posthoc/gwis_sig_results_input_", exposure, ".rds")) %>% 
  mutate(snps = paste0('chr', gsub("\\:", "\\_", SNP))) %>% 
  filter(method != paste0("chiSqCase_", exposure))


# D|G main effects table - check if risk alleles need to be recoded
lm_out <- data.frame(snps = snps$snps, coeffs = map_chr(snps$snps, ~ lm(paste0("outcome ~ ", .x, "*", exposure, "+", paste0(covariates_set1, collapse = "+")), data = figi)$coef[[2]])) %>% 
  unique(.)
lm_out
saveRDS(lm_out, file = paste0(output_dir, "lm_g_maineffects.rds"))


# helper function

# helper function
posthoc_helper <- function(x, method) {
  # pvalues dataframe
  walk(x, ~ pval_summary(figi, exposure, .x, covariates_list, method))
  # model estimates (stargazer html)
  walk(x, ~ fit_gxe_covars(figi, exposure, .x, covariates_list, method))
  # stratified odds ratio
  walk(x, ~ fit_stratified_or_q4(figi, exposure, snp = .x, hrc_version = hrc_version, covariates = covariates_set1, mod = mod, dosage = T))
  # locuszoom and functional annotation plots
  locuszoom_helper <- function(snp) {
    tmp <- gsub("chr", "", snp)
    snp_chr <- strsplit(tmp, "_")[[1]][1]
    snp_bp <- as.numeric(strsplit(tmp, '_')[[1]][2])
    paste(snp_chr, snp_bp, sep = ":")
  }
  snps_plot <- map_chr(x, locuszoom_helper)
  walk(snps_plot, ~ system(paste("bash ~/Dropbox/FIGI/FIGI_code/results/posthoc/posthoc_02_locuszoom.sh", exposure, hrc_version, .x, method)))
  walk(snps_plot, ~ system(paste("Rscript ~/Dropbox/FIGI/FIGI_code/results/posthoc/posthoc_03_functional_annotation.R", exposure, .x, method)))
}


#-----------------------------------------------------------------------------#
# suggestive GxE hits (1df)
# generate all components
#-----------------------------------------------------------------------------#
snps <- readRDS(paste0("/media/work/gwis/posthoc/gwis_sig_results_input_", exposure, ".rds")) %>% 
  filter(method == paste0("chiSqGxE_", exposure, "_clump")) %>% 
  arrange(Chromosome, Location)
snps <- paste0('chr', gsub("\\:", "\\_", snps$SNP))

posthoc_helper(snps, 'chiSqGxE')

#-----------------------------------------------------------------------------#
# findings from chiSq3df_folate_dietqc2_no_gwas_no_marginal_no_ge
#
# tophits
# - 2:135759095:C:T	
#-----------------------------------------------------------------------------#
snps <- readRDS(paste0("/media/work/gwis/posthoc/gwis_sig_results_input_", exposure, ".rds")) %>% 
  filter(method == paste0("chiSq3df_", exposure, "_no_gwas_no_marginal_no_ge"), 
         SNP %in% c("2:135759095:C:T")) %>% 
  arrange(Chromosome, Location)
snps <- paste0('chr', gsub("\\:", "\\_", snps$SNP))

posthoc_helper(snps, 'chiSq3df')


#-----------------------------------------------------------------------------#
# LD clumped results (D|G)
# 
# 8:117630683:A:C this one only for now
# 8:117596426:A:G not in input file, but can add. although not necessarily - same locus
#-----------------------------------------------------------------------------#
snps <- readRDS(paste0("/media/work/gwis/posthoc/gwis_sig_results_input_", exposure, ".rds")) %>% 
  filter(SNP %in% c("8:117630683:A:C")) %>% 
  arrange(Chromosome, Location)
snps <- unique(paste0('chr', gsub("\\:", "\\_", snps$SNP)))

posthoc_helper(snps, 'chiSqGxE')

