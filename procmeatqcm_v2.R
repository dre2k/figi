#=============================================================================#
# FIGI GxE redmeatqcm results
#=============================================================================#

# setup -------------------------------------------------------------------


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
library(qs)
rm(list = ls())
# source("functions.R")



# input variables
exposure = 'procmeatqcm_v2'
hrc_version = 'v3.0'
covariates <- sort(c('age_ref_imp', 'sex', 'energytot_imp', 'study_gxe', 'pc1', 'pc2', 'pc3'))
path = glue("/media/work/gwis_test/{exposure}/")

# pca
pca <- fread("/media/work/gwis_test/PCA/20210222/figi_gxe_pca_update.eigenvec") %>% 
  dplyr::rename(vcfid = IID) %>% 
  dplyr::select(-`#FID`) %>% 
  dplyr::rename_with(tolower)


# input data
esubset <- readRDS(glue("/media/work/gwis_test/{exposure}/data/FIGI_{hrc_version}_gxeset_{exposure}_basic_covars_glm.rds")) %>% 
  pull(vcfid)

input_data <- readRDS(glue("/media/work/gwis_test/data/FIGI_{hrc_version}_gxeset_analysis_data_glm.rds")) %>% 
  filter(vcfid%in% esubset) %>%
  dplyr::select(-starts_with("pc")) %>% 
  left_join(pca, 'vcfid') %>% 
  mutate(procmeatqcm_v2 = as.numeric(procmeatqcm_v2))



#-----------------------------------------------------------------------------#
# Main effects ------------------------------------------------------------
#-----------------------------------------------------------------------------#

# ------ meta-analysis ------ #
output_dir = as.character(glue("{path}/output/posthoc/"))
covariates_meta <- sort(covariates[which(!covariates %in% c(paste0(rep('pc', 20), seq(1, 20)), "study_gxe"))])

create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "all", forest_height = 15, forest_width = 9, categorical = F)
create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "proximal", forest_height = 13, forest_width = 9, categorical = F)
create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "distal", forest_height = 13, forest_width = 9, categorical = F)
create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "rectal", forest_height = 13, forest_width = 9, categorical = F)
create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "female", forest_height = 13, forest_width = 9, categorical = F)
create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "male", forest_height = 13, forest_width = 9, categorical = F)

# ------- stratified pooled analysis ------- #
pooled_analysis_glm(input_data, exposure = exposure, hrc_version = hrc_version, covariates = covariates, strata = 'sex', filename_suffix = paste0(paste0(sort(covariates), collapse = '_'), "_stratified_sex"), output_dir = glue("{path}/output/posthoc/"))

pooled_analysis_glm(input_data, exposure = exposure, hrc_version = hrc_version, covariates = covariates, strata = 'study_design', filename_suffix = paste0(paste0(sort(covariates), collapse = '_'), "_stratified_study_design"), output_dir = glue("{path}/output/posthoc/"))

pooled_analysis_multinom(input_data, exposure = exposure, hrc_version = hrc_version, covariates = covariates, filename_suffix = paste0(paste0(sort(covariates), collapse = '_'), "_stratified_cancer_site_sum2"), output_dir = glue("{path}/output/posthoc/"))



#-----------------------------------------------------------------------------#
# GxE additional analysis ---- 
#-----------------------------------------------------------------------------#


snps <- c("4:11944206:A:G", "15:92814380:C:T", "8:59807881:A:C", "1:21015685:G:A", "14:24188606:C:T", "7:100511033:A:G", "16:23076410:C:T", "9:114385930:A:G", "18:21029284:C:T", "5:71909347:C:T", "15:98258859:C:T", "2:43003484:A:C", "13:21949790:T:A")

# MAF ---------------------------------------------------------------------
walk(snps, ~ create_aaf_study_plot(data = input_data, exposure = exposure, hrc_version = hrc_version, snp = .x, path = path))


# stratified odds ratio ---------------------------------------------------
walk(snps , ~ fit_stratified_or_continuous(data_epi = input_data, exposure = exposure, snp = .x, hrc_version = hrc_version, covariates = covariates, dosage = F, path = glue("{path}/output"), flip_allele = F))

# interaction plots -------------------------------------------------------
walk(snps, ~ iplot_wrapper(data_epi = input_data, exposure = exposure, snp = .x,  covariates = covariates, hrc_version = hrc_version, path = glue("{path}/output"), flip_allele = F))


# ================================================================== #
# ======= rmarkdown reports ---- 
# ================================================================== #
main_effects_report(exposure = exposure, hrc_version = hrc_version, covariates = covariates, path = path)
gwis_report(exposure = exposure, hrc_version = hrc_version, covariates = covariates)
posthoc_report(exposure = exposure)




