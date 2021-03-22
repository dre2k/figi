#=============================================================================#
# FIGI GxE results
#=============================================================================#
library(tidyverse)
library(data.table)
library(qqman)
library(EasyStrata)
library(lmtest)
library(figifs)
library(forcats)
library(car)
library(grid)
library(gridExtra)
library(stargazer)
library(nnet)
library(gtools)
library(glue)
library(ramwas)
library(flextable)
library(gtools)
library(interactionR)
library(epiR)
library(flextable)
library(jtools)
library(interactions)
library(msm)
rm(list = ls())

source("functions.R")

# input variables
exposure = 'calcium_totqc2'
hrc_version = 'v3.0'
output_dir = paste0("/media/work/gwis/posthoc/", exposure, "/")
annotation_file <- 'gwas_141_ld_annotation_july2020.txt'

covariates = sort(c('age_ref_imp', 'sex', 'energytot_imp', 'study_gxe', 'pc1', 'pc2', 'pc3'))
covariates_set1 <- covariates
covariates_list <- list(covariates_set1)
mod <- 'age_ref_imp+sex+energytot_imp+pc1+pc2+pc3+study_gxe'


#-----------------------------------------------------------------------------#
# edit code below as necessary
#-----------------------------------------------------------------------------#

# input
# have to recode calcium for posthoc analyses (protective)
figi <- posthoc_input(exposure, hrc_version, glue('gwis_sig_results_output_{exposure}.rds')) %>% 
	mutate(calcium_totqc2 = abs(calcium_totqc2 - 3))

snps <- readRDS(glue(output_dir, "gwis_sig_results_input_{exposure}.rds"))
table(snps$method)


# ------- GxE 1DF findings ------- #
snps_filter <- snps %>% 
	dplyr::filter(method == glue("manhattan_chiSqGxE_{exposure}_clump_df"),
		      Pval <= 5e-6) %>% 
dplyr::mutate(snps = paste0('chr', gsub("\\:", "\\_", SNP))) %>% 
pull(snps)

posthoc_run_models(snps_filter, 'chiSqGxE', quartile = T)
posthoc_run_reri_plot(snps_filter)
posthoc_create_plots(snps_filter, 'chiSqGxE')


# ------- Expectation based hybrid... (just do all...) -------- #
snps_filter <- readRDS(glue(output_dir, "twostep_wht_chiSqGE_{exposure}_expectation_hybrid_df.rds")) %>% 
	dplyr::mutate(snps = paste0('chr', gsub("\\:", "\\_", SNP))) %>% 
	pull(snps)

posthoc_run_models(snps_filter, 'chiSqGxE', quartile = T)
posthoc_run_reri_plot(snps_filter)
posthoc_create_plots(snps_filter, 'chiSqGxE')


# ------- 2DF/3DF clumped removing GWAS and or EG (to check 'novel' D|G results) ------- #
# as first pass, only generate locuszoom plots to gauge whether the region is already capture by GECCO meta-analysis
snps_filter <- readRDS(glue(output_dir, "manhattan_chiSq2df_{exposure}_no_gwas_clump_df.rds")) %>% 
	dplyr::mutate(snps = paste0('chr', gsub("\\:", "\\_", SNP))) %>% 
	pull(snps)
posthoc_create_plots(snps_filter, 'chiSq2df')

snps_filter <- readRDS(glue(output_dir, "manhattan_chiSq3df_{exposure}_no_gwas_no_ge_clump_df.rds")) %>% 
	dplyr::mutate(snps = paste0('chr', gsub("\\:", "\\_", SNP))) %>% 
	pull(snps)
posthoc_create_plots(snps_filter, 'chiSq3df')


