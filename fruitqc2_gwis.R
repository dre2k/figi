#=============================================================================#
# FIGI GxE fruitqc2 results
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

# input variables
exposure = 'fruitqc2'
hrc_version = 'v2.3'
annotation_file <- 'gwas_200_ld_annotation_feb2021.txt'
path = glue("/media/work/gwis_test/{exposure}/")

covariates <- sort(c('age_ref_imp', 'sex', 'energytot_imp', 'study_gxe', 'pc1', 'pc2', 'pc3'))

# input data
esubset <- readRDS(glue("{path}/data/FIGI_{hrc_version}_gxeset_{exposure}_basic_covars_glm.rds")) %>% pull(vcfid)
input_data <- readRDS(glue("/media/work/gwis/data/FIGI_EpiData/FIGI_{hrc_version}_gxeset_analysis_data_glm.rds")) %>% 
  filter(vcfid%in% esubset)


#-----------------------------------------------------------------------------#
# main effects ----
#-----------------------------------------------------------------------------#

# ------ meta-analysis ------ #
output_dir = as.character(glue("{path}/output/posthoc/"))
covariates_meta <- sort(covariates[which(!covariates %in% c(paste0(rep('pc', 20), seq(1, 20)), "study_gxe"))])

create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "all", forest_height = 15, categorical = F)
create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "proximal", forest_height = 13, categorical = F)
create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "distal", forest_height = 13, categorical = F)
create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "rectal", forest_height = 13, categorical = F)
create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "female", forest_height = 13, categorical = F)
create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "male", forest_height = 13, categorical = F)


# ------- stratified pooled analysis ------- #
pooled_analysis_glm(input_data, exposure = exposure, hrc_version = hrc_version, covariates = covariates, strata = 'sex', filename_suffix = paste0(paste0(sort(covariates), collapse = '_'), "_stratified_sex"), output_dir = glue("{path}/output/posthoc/"))

pooled_analysis_glm(input_data, exposure = exposure, hrc_version = hrc_version, covariates = covariates, strata = 'study_design', filename_suffix = paste0(paste0(sort(covariates), collapse = '_'), "_stratified_study_design"), output_dir = glue("{path}/output/posthoc/"))

pooled_analysis_multinom(input_data, exposure = exposure, hrc_version = hrc_version, covariates = covariates, filename_suffix = paste0(paste0(sort(covariates), collapse = '_'), "_stratified_cancer_site_sum2"), output_dir = glue("{path}/output/posthoc/"))





#-----------------------------------------------------------------------------#
# GxE additional analysis ---- 
#-----------------------------------------------------------------------------#
source(glue("~/git/figifs/R/01_process.R"))
source(glue("~/git/figifs/R/02_plots.R"))
source(glue("~/git/figifs/R/03_posthoc.R"))
source(glue("~/git/figifs/R/03_posthoc_iplot.R"))
source(glue("~/git/figifs/R/03_posthoc_stratified_or.R"))


# output RERI plots (can't install package on CARC)
snps <- c("14:74029409:C:T", "14:74029049:G:C")
walk(snps, ~ reri_wrapper(data_epi = input_data, exposure = exposure, snp = .x, covariates = covariates, path = glue("{path}/output")))


input_data_tmp <- input_data %>% 
  mutate(fruitqc2 = abs(3 - as.numeric(fruitqc2)))
snps <- c("1:72729142:A:G")
walk(snps, ~ reri_wrapper(data_epi = input_data_tmp, exposure = exposure, snp = .x, covariates = covariates, path = glue("{path}/output")))









# output GxE models adjusted by different covariate sets
covariates_sets <- list(covariates, 
                        c(covariates, 'bmi5'), 
                        c(covariates, 'bmi5', 'smk_ever'), 
                        c(covariates, 'bmi5', 'smk_ever', 'fruitqc2', 'vegetableqc2'))
covariates_sets <- list(covariates)
walk(snps, ~ fit_gxe_covars(data_epi = input_data, exposure = exposure, snp = .x, covariates_list = covariates_sets, method = 'chiSqGxE', path = glue("{path}/output")))
suggestive_gxe <- fread(glue("{path}/data/FIGI_v2.3_gxeset_diab_chiSqGxE_ldclump.clumped"))
walk(suggestive_gxe$SNP, ~ fit_gxe_covars(data_epi = input_data, exposure = exposure, snp = .x, covariates_list = covariates_sets, method = 'chiSqGxE', path = glue("{path}/output")))



# output easystrata manhattan plot for Nikolas that excludes the right set of annotations (200 GWAS instead of 140)
gxe <- readRDS("/media/work/gwis/results/fruitqc2/processed/FIGI_v2.3_gxeset_fruitqc2_basic_covars_gxescan_results.rds")

create_manhattanplot(x = gxe, exposure = exposure, stat = 'chiSq3df', annotation_file = annotation_file, output_dir = "~/Dropbox")
gwas <- fread("/home/rak/data/Annotations/gwas_200_ld_annotation_feb2021.txt") %>% 
  mutate(SNP2 = paste0(Chr, ":", Pos))
gxe_nogwas <- dplyr::filter(gxe, !SNP2 %in% gwas$SNP2)
create_manhattanplot(x = gxe, exposure = exposure, stat = 'chiSq3df', annotation_file = annotation_file, output_dir = "~/Dropbox")
create_manhattanplot(x = gxe_nogwas, exposure = exposure, stat = 'chiSq3df', annotation_file = annotation_file, output_dir = "~/Dropbox")

# try again - remove flanking regions for better visualizaation
create_manhattanplot_ramwas(data = gxe, exposure = exposure, statistic = 'chiSqGxE', hrc_version = 'v2.3', path = "~/Dropbox/", sig_line = 5e-8)
gxe_better <- gxe_nogwas %>% 
  filter(chiSqG_p > 5e-7)

create_manhattanplot_ramwas(data = gxe, exposure = exposure, statistic = 'chiSqGxE', hrc_version = 'v2.3', path = "~/Dropbox/", sig_line = 5e-8)
create_manhattanplot_ramwas(data = gxe_better, exposure = exposure, statistic = 'chiSq3df', hrc_version = 'v2.3', path = "~/Dropbox/", sig_line = 5e-8)
create_manhattanplot_ramwas(data = gxe_better, exposure = exposure, statistic = 'chiSq2df', hrc_version = 'v2.3', path = "~/Dropbox/", sig_line = 5e-8)





# output RERI plots for suggestive findings
input_data_recode <- input_data %>% 
  mutate(fruitqc2 = abs(fruitqc2-3))
  
snps_out <- fread(glue("{path}/data/FIGI_{hrc_version}_gxeset_{exposure}_chiSqGxE_ldclump.clumped"), header = T, stringsAsFactors = F) %>%
  dplyr::filter(SNP != "2:131698423:A:G", 
                SNP != "10:38739240:C:A") %>% # for some reason these guys don't seem to generate locuszoom plots. 
  dplyr::arrange(CHR, SNP) %>%
  dplyr::pull(SNP)

walk(snps_out, ~ reri_wrapper(data_epi = input_data_recode, exposure = exposure, snp = .x, covariates = covariates, path = glue("{path}/output")))


# output AAF by study_gxe plots for suggestive findings
walk(snps_out, ~ create_aaf_study_plot(data = input_data_recode, exposure,  hrc_version, snp = .x, path = path))

# ================================================================== #
# ======= rmarkdown reports ---- 
# ================================================================== #
main_effects_report(exposure = exposure, hrc_version = hrc_version, covariates = covariates, path = path)

gwis_report(exposure = exposure, 
            hrc_version = hrc_version, 
            covariates = covariates)

posthoc_report(exposure = exposure, 
               hrc_version = hrc_version,
               covariates = covariates,
               path = path)




# ----- suggestive markers ------ # 
rmarkdown::render(glue("/home/rak/git/figi/{exposure}_posthoc_suggestive.Rmd"), 
                  params = list(exposure = exposure, hrc_version = hrc_version, 
                                covariates = covariates, path = path), 
                  output_file = glue("~/Dropbox/FIGI/Results/{exposure}_posthoc_suggestive.html"))



snps_out <- fread(glue("/media/work/gwis_test/fruitqc2/data/FIGI_v2.3_gxeset_fruitqc2_chiSqGxE_ldclump.clumped"))





