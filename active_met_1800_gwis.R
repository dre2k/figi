#=============================================================================#
# FIGI GxE smk_aveday results
#=============================================================================#

# Setup -----------------------------------------------------------------
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
exposure = 'active_met_1800'
hrc_version = 'v3.1'
annotation_file <- 'gwas_200_ld_annotation_feb2021.txt'
covariates <- sort(c('age_ref_imp', 'sex', 'study_gxe', 'energytot_imp', 'pc1', 'pc2', 'pc3'))
path = glue("/media/work/gwis_test/{exposure}/")


# input data
esubset <- readRDS(glue("/media/work/gwis_test/{exposure}/data/FIGI_{hrc_version}_gxeset_{exposure}_basic_covars_glm.rds")) %>% 
  pull(vcfid)

input_data <- readRDS(glue("/media/work/gwis_test/data/FIGI_{hrc_version}_gxeset_analysis_data_glm.rds")) %>%
  # filter(!is.na(methrswklns)) %>% 
  filter(vcfid %in% esubset) %>%
  mutate(active_met_1800 = ifelse(is.na(methrswklns),NA,
                                 ifelse(methrswklns<18, 'No','Yes')))



# input_data %>% 
#   mutate(methrswklns = as.numeric(methrswklns)) %>% 
#   filter(methrswklns != 0) %>% 
#   mutate(study = fct_reorder(study, methrswklns, .fun='median')) %>% 
#   ggplot(aes(x=study, y=methrswklns, fill = study)) + 
#   geom_boxplot(outlier.color = 'red') + 
#   theme_bw() + 
#   theme(axis.text.x = element_text(angle = 270))
# 
# input_data %>% 
#   mutate(methrswklns = as.numeric(methrswklns)) %>% 
#   filter(methrswklns != 0) %>% 
#   mutate(study = fct_reorder(study, methrswklns, .fun='median')) %>% 
#   ggplot(aes(x=study_gxe, y=methrswklns, fill = study)) + 
#   geom_boxplot(outlier.color = 'red') + 
#   theme_bw() + 
#   theme(axis.text.x = element_text(angle = 270))
# 
# 
# 
# ggsave("~/Dropbox/methrswklns_boxplot.png", width = 8, height = 6)



#-----------------------------------------------------------------------------#
# main effects ------------------------------------------------------------
#-----------------------------------------------------------------------------#


# ------ meta-analysis ------ #
output_dir = as.character(glue("{path}/output/posthoc/"))
covariates_meta <- sort(covariates[which(!covariates %in% c(paste0(rep('pc', 20), seq(1, 20)), "study_gxe"))])

create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "all", forest_height = 15, forest_width = 9, categorical = T)
create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "proximal", forest_height = 13, forest_width = 9, categorical = T)
create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "distal", forest_height = 13, forest_width = 9, categorical = T)
# create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "rectal", forest_height = 13, forest_width = 9, categorical = T)
create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "female", forest_height = 13, forest_width = 9, categorical = T)
create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "male", forest_height = 13, forest_width = 9, categorical = T)

# ------- stratified pooled analysis ------- #
pooled_analysis_glm(input_data, exposure = exposure, hrc_version = hrc_version, covariates = covariates, strata = 'sex', filename_suffix = paste0(paste0(sort(covariates), collapse = '_'), "_stratified_sex"), output_dir = glue("{path}/output/posthoc/"))

pooled_analysis_glm(input_data, exposure = exposure, hrc_version = hrc_version, covariates = covariates, strata = 'study_design', filename_suffix = paste0(paste0(sort(covariates), collapse = '_'), "_stratified_study_design"), output_dir = glue("{path}/output/posthoc/"))

pooled_analysis_multinom(input_data, exposure = exposure, hrc_version = hrc_version, covariates = covariates, filename_suffix = paste0(paste0(sort(covariates), collapse = '_'), "_stratified_cancer_site_sum2"), output_dir = glue("{path}/output/posthoc/"))

# pooled analyses stratified by BMI level 
pooled_analysis_glm(input_data, exposure = exposure, hrc_version = hrc_version, covariates = covariates, strata = 'bmic3', filename_suffix = paste0(paste0(sort(covariates), collapse = '_'), "_stratified_bmic3"), output_dir = glue("{path}/output/posthoc/"))




# # ------ meta-analysis ------ #
# # additional covariates
# # output_dir = as.character(glue("{path}/output/posthoc/"))
# covariates_meta <- sort(covariates[which(!covariates %in% c(paste0(rep('pc', 20), seq(1, 20)), "study_gxe"))])
# 
# covariates_meta <- c(covariates_meta, 'bmi')
# covariates_meta <- c(covariates_meta, 'energytot_imp' )
# covariates_meta <- c(covariates_meta, 'bmi', 'energytot_imp')
# 
# create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "all", forest_height = 15, categorical = F)
# create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "proximal", forest_height = 13, categorical = F)
# create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "distal", forest_height = 13, categorical = F)
# create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "rectal", forest_height = 13, categorical = F)
# create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "female", forest_height = 13, categorical = F)
# create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "male", forest_height = 13, categorical = F)
# 
# 
# # ------- stratified pooled analysis ------- #
# pooled_analysis_glm(input_data, exposure = exposure, hrc_version = hrc_version, covariates = covariates, strata = 'sex', filename_suffix = paste0(paste0(sort(covariates), collapse = '_'), "_stratified_sex"), output_dir = glue("{path}/output/posthoc/"))
# 
# pooled_analysis_glm(input_data, exposure = exposure, hrc_version = hrc_version, covariates = covariates, strata = 'study_design', filename_suffix = paste0(paste0(sort(covariates), collapse = '_'), "_stratified_study_design"), output_dir = glue("{path}/output/posthoc/"))
# 
# pooled_analysis_multinom(input_data, exposure = exposure, hrc_version = hrc_version, covariates = covariates, filename_suffix = paste0(paste0(sort(covariates), collapse = '_'), "_stratified_cancer_site_sum2"), output_dir = glue("{path}/output/posthoc/"))



# main effects additional ----
# 
# output_dir_dropbox = paste0("~/Dropbox/Presentations/", exposure, "/")
# 
# tmp1 <- readRDS(paste0("/media/work/gwis/results/input/FIGI_", hrc_version, "_gxeset_", exposure, "_basic_covars_glm.rds")) %>% 
#   pull(vcfid)
# xx <- filter(input_data, vcfid %in% tmp1) %>% 
#   group_by(study_gxe)
# 
# 
# results_beta <- dplyr::do(xx, broom::tidy(glm(outcome ~ folate_dietqc2 + age_ref_imp + sex + energytot_imp, data = . , family = 'binomial'))) %>% 
#   dplyr::filter(grepl("folate_dietqc2", term)) %>% 
#   dplyr::arrange(study_gxe) %>% 
#   inner_join(unique(xx[,c('study_gxe', 'study_design')]), 'study_gxe')
# 
# results_meta <- meta::metagen(estimate,
#                               std.error,
#                               data=results_beta,
#                               studlab=paste(study_gxe),
#                               comb.fixed = FALSE,
#                               comb.random = TRUE,
#                               method.tau = "SJ",
#                               hakn = TRUE,
#                               prediction=TRUE,
#                               sm="OR", 
#                               byvar=study_design)
# 
# fo <- find.outliers(results_meta)
# fo
# meta::forest(results_meta,
#              layout = "JAMA",
#              text.predict = "95% CI",
#              col.predict = "black",
#              # leftcols = c("studlab", "Control", "Case", "N", "effect", "ci", "w.random"),
#              digits.addcols=0,
#              study.results=T,
#              prediction = F,
#              col.random = 'red')
# 
# png(paste0(output_dir_dropbox, "meta_analysis_", "folate_dietqc2",  "_", "original_outliers_removed", ".png"), height = 17, width = 8.5, units = 'in', res = 150)                                                          
# forest(fo)
# dev.off()
# 
# 
# # leave on out (influence analysis)
# inf.analysis <- InfluenceAnalysis(x = results_meta,
#                                   random = TRUE)
# 
# summary(inf.analysis)
# 
# plot(inf.analysis, "influence")
# plot(inf.analysis, "baujat")
# twostep_eh_snps(gxe, 'chiSqEDGE')
# plot(inf.analysis, "es")
# plot(inf.analysis, "i2")



#-----------------------------------------------------------------------------#
# GxE additional analyses -------------------------------------------------
#-----------------------------------------------------------------------------#

# ---- MAF ---- #
snps <- c("15:33008870:T:C")
walk(snps, ~ create_aaf_study_plot(data = input_data, exposure = exposure, hrc_version = hrc_version, snp = .x, path = path))


# ---- stratified odds ratios ---- #
# stratified odds ratio table with and without smoking
input_data2 <- input_data %>% 
  mutate(active_met_1800 = as.factor(ifelse(active_met_1800 == "No", 0, 1)))

walk(snps, ~ fit_stratified_or(
  data_epi = input_data2,
  exposure = exposure,
  snp = .x,
  hrc_version = hrc_version,
  covariates = covariates,
  path = glue("{path}/output/")))

# snps <- c("6:31917540:T:C", "8:138788813:C:A")
# walk(snps, ~ fit_gxe_covars(data_epi = input_data, exposure = exposure, snp = .x, covariates_list = covariates_sets, method = 'chiSq3df', path = glue("{path}/output")))



# ---- GxE model stratified by study design, tumor, sex, bmic3 ---- #
# Need to standardize how you handle factors vs numerical in your functions

input_data <- input_data %>% 
  mutate(active_met_1800 = as.factor(ifelse(active_met_1800 == "No", 0, 1)))


# GxE findings
method = 'chiSqGxE'

walk(snps, ~ fit_gxe_stratified(data_epi = input_data, exposure = exposure, snp = .x, covariates = covariates, strata = 'sex', method = method, flip_allele = F, path = path))
walk(snps, ~ fit_gxe_stratified(data_epi = input_data, exposure = exposure, snp = .x, covariates = covariates, strata = 'study_design', method = method, flip_allele = F, path = path))
walk(snps, ~ fit_gxe_stratified(data_epi = input_data, exposure = exposure, snp = .x, covariates = covariates, strata = 'cancer_site_sum2', method = method, flip_allele = F, path = path))
walk(snps, ~ fit_gxe_stratified(data_epi = input_data, exposure = exposure, snp = .x, covariates = covariates, strata = 'bmic3', method = method, flip_allele = F, path = path))

# exclude adenomas
input_data_crc <- readRDS(glue("/media/work/gwis_test/data/FIGI_{hrc_version}_gxeset_analysis_data_glm.rds")) %>%
  filter(vcfid %in% esubset) %>%
  filter(is.na(adenoma)) %>% 
  mutate(active_met_1800 = ifelse(is.na(methrswklns),NA,
                                  ifelse(methrswklns<18, 'No','Yes')), 
         active_met_1800 = as.factor(ifelse(active_met_1800 == "No", 0, 1)))


walk(snps, ~ fit_gxe_stratified(data_epi = input_data_crc, exposure = exposure, snp = .x, covariates = covariates, strata = 'sex', method = method, flip_allele = F, path = path))
walk(snps, ~ fit_gxe_stratified(data_epi = input_data_crc, exposure = exposure, snp = .x, covariates = covariates, strata = 'study_design', method = method, flip_allele = F, path = path))
walk(snps, ~ fit_gxe_stratified(data_epi = input_data_crc, exposure = exposure, snp = .x, covariates = covariates, strata = 'cancer_site_sum2', method = method, flip_allele = F, path = path))
walk(snps, ~ fit_gxe_stratified(data_epi = input_data_crc, exposure = exposure, snp = .x, covariates = covariates, strata = 'bmic3', method = method, flip_allele = F, path = path))


# output RERI plots (can't install package on CARC yet)
walk(snps, ~ reri_wrapper(data_epi = input_data, exposure = exposure, snp = .x, covariates = covariates, path = glue("{path}/output")))





# ================================================================== #
# ======= rmarkdown reports ---- 
# ================================================================== #

main_effects_report(exposure = exposure, hrc_version = hrc_version, covariates = covariates, path = path)
gwis_report(exposure = exposure, hrc_version = hrc_version, covariates = covariates)
# posthoc_report(exposure = exposure)

rmarkdown::render(glue("/home/rak/git/figi/{exposure}_posthoc.Rmd"), 
                  params = list(exposure = exposure, 
                                covariates = covariates), 
                  output_file = glue("~/Dropbox/FIGI/Results/{exposure}_posthoc.html"))







