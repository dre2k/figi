#=============================================================================#
# FIGI GWIS Results
# Updated 04/16/2020
# Updated 08/18/2020
#
# Create plots and output dataframe of significant results
# 
# Outputs
# - QQ plots for each method 
# - Manhattan plots for each method (using LD based annotation file based on 140 GECCO GWAS loci)
# - Dataframe of significant results for each method (GxE 1df threshold is lower, 5e-6)
#
# Comments:
# - covariates - all analyses ran with a basic set of covariates
# - omit LD clumped results now. Difficult to defend without permutation testing
#=============================================================================#
library(tidyverse)
library(data.table)
library(qqman)
library(EasyStrata)
library(lmtest)
library(figifs)
library(glue)
library(ramwas)


code_dir <- "/home/rak/Dropbox/FIGI/FIGI_code/results/gwis/"
output_dir <- paste0("/media/work/gwis/posthoc/", exposure, "/")

# gwis covariates for plot title etc
if(exposure %in% c("folate_totqc2", "folate_dietqc2", "alcoholc_moderate", "alcoholc_heavy_vs_moderate", "calcium_totqc2", "calcium_dietqc2", "fiberqc2", "fruitqc2", "vegetableqc2", "procmeatqc2", "redmeatqc2")) {
  covariates <- c("age_ref_imp", "sex", "energytot_imp", "study_gxe", "pc1", "pc2", "pc3")
} else if(exposure %in% c("hrt_ref_pm2", "eo_ref_pm_gxe", "ep_ref_pm_gxe", "pure_eo_allNo", "pure_ep_allNo")) {
  covariates <- c("age_ref_imp", "study_gxe", "pc1", "pc2", "pc3")
} else if(exposure == "p_diet_std") {
  covariates <- c("age_ref_imp", "sex", "energytot_imp", "study_gxe", "pc1", "pc2", "pc3", "bmi5", "height10")
} else {
  covariates <- c("age_ref_imp", "sex", "study_gxe", "pc1", "pc2", "pc3")
}

# manhattan plot SNP annotation file (Jeroen GWAS N = 140)
annotation_file <- 'gwas_141_ld_annotation_july2020.txt'

# ----------------------------------- #
# remove GWAS hits, D|G, E|G controls for visualize top p-value SNP in each loci
# using plink LD clump output files (based on 1000 FIGI controls)
# 
# - edit - will only remove GWAS hits
# ----------------------------------- #
gxe_nogwas <- gxe %>%
  dplyr::filter(!SNP2 %in% exclude_gwas)

#-----------------------------------------------------------------------------#
# Manhattan Plots ----
#-----------------------------------------------------------------------------#
message("create Manhattan plots")
source(paste0(code_dir, "gwis_02b_manhattan_plots.R"))

#-----------------------------------------------------------------------------#
# Two-step methods
# step1:
# - using main effects statistics from exposure subset analysis, 
#   main effects models are adjusted by exposure
#-----------------------------------------------------------------------------#
message("Various two-step methods")
source(paste0(code_dir, "gwis_02c_two_step.R"))
