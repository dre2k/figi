#=============================================================================#
# FIGI GxE smk_aveday results
#=============================================================================#


# Setup -------------------------------------------------------------------
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

# input variables
exposure = 'active_met_875'
hrc_version = 'v3.1'
annotation_file <- 'gwas_200_ld_annotation_feb2021.txt'
covariates <- sort(c('age_ref_imp', 'sex', 'study_gxe', 'energytot_imp', 'pc1', 'pc2', 'pc3'))
path = glue("/media/work/gwis_test/{exposure}/")

# input data
esubset <- readRDS(glue("/media/work/gwis_test/{exposure}/data/FIGI_{hrc_version}_gxeset_{exposure}_basic_covars_glm.rds")) %>% 
  pull(vcfid)


pca <- fread("/media/work/gwis_test/PCA/20210222/figi_gxe_pca_update.eigenvec") %>% 
  rename(vcfid = IID) %>% 
  rename_with(tolower) %>% 
  dplyr::select(-`#fid`)

input_data <- readRDS(glue("/media/work/gwis_test/data/FIGI_{hrc_version}_gxeset_analysis_data_glm.rds")) %>% 
  # filter(!is.na(methrswklns)) %>% 
  filter(vcfid %in% esubset) %>% 
  dplyr::select(-starts_with('pc')) %>% 
  inner_join(pca, 'vcfid') %>% 
  mutate(active_met_875 = ifelse(is.na(methrswklns),NA,
                          ifelse(methrswklns<8.75, "No", "Yes")))




#-----------------------------------------------------------------------------#
# Main Effects ------------------------------------------------------------
# (performed separately in another script.. )
#-----------------------------------------------------------------------------#

# ------ meta-analysis ------ #
output_dir = as.character(glue("{path}/output/posthoc/"))
covariates_meta <- sort(covariates[which(!covariates %in% c(paste0(rep('pc', 20), seq(1, 20)), "study_gxe"))])

create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "all", forest_height = 15, forest_width = 9, categorical = T)
create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "proximal", forest_height = 13, forest_width = 9, categorical = T)
create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "distal", forest_height = 13, forest_width = 9, categorical = T)
create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "rectal", forest_height = 13, forest_width = 9, categorical = T)
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
# 
# 
# # main effects additional ----
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
snps <- c("15:33001734:T:A", "15:32994756:T:C", "15:32994777:C:T", "15:32995298:C:A", "15:32995764:C:G", "13:61572969:G:A")
walk(snps, ~ create_aaf_study_plot(data = input_data, exposure = exposure, hrc_version = hrc_version, snp = .x, path = path))


# ---- stratified odds ratios ---- #
# flip two of the alleles so param directions are consistent
input_data2 <- input_data %>% 
  mutate(active_met_875 = as.factor(ifelse(active_met_875 == "No", 0, 1)))

snps <- c("15:33001734:T:A", "15:32995298:C:A", "15:32995764:C:G", "13:61572969:G:A")
snps_flip <- c("15:32994756:T:C", "15:32994777:C:T")

walk(snps, ~ fit_stratified_or(
  data_epi = input_data2,
  exposure = exposure,
  snp = .x,
  hrc_version = hrc_version,
  covariates = covariates,
  path = glue("{path}/output/")))

walk(snps_flip, ~ fit_stratified_or(
  data_epi = input_data2,
  exposure = exposure,
  snp = .x,
  hrc_version = hrc_version,
  covariates = covariates,
  flip_allele = T,
  path = glue("{path}/output/")))



# ---- GxE model stratified by study design, tumor, sex, bmic3 ---- #
# Need to standardize how you handle factors vs numerical in your functions
input_data2 <- input_data %>% 
  mutate(active_met_875 = as.numeric(ifelse(active_met_875 == "No", 0, 1)))

# GxE findings
method = 'chiSqGxE'

walk(snps, ~ fit_gxe_stratified(data_epi = input_data2, exposure = exposure, snp = .x, covariates = covariates, strata = 'sex', method = method, flip_allele = F, path = path))
walk(snps, ~ fit_gxe_stratified(data_epi = input_data2, exposure = exposure, snp = .x, covariates = covariates, strata = 'study_design', method = method, flip_allele = F, path = path))
walk(snps, ~ fit_gxe_stratified(data_epi = input_data2, exposure = exposure, snp = .x, covariates = covariates, strata = 'cancer_site_sum2', method = method, flip_allele = F, path = path))
walk(snps, ~ fit_gxe_stratified(data_epi = input_data2, exposure = exposure, snp = .x, covariates = covariates, strata = 'bmic3', method = method, flip_allele = F, path = path))

walk(snps_flip, ~ fit_gxe_stratified(data_epi = input_data2, exposure = exposure, snp = .x, covariates = covariates, strata = 'sex', method = method, flip_allele = T, path = path))
walk(snps_flip, ~ fit_gxe_stratified(data_epi = input_data2, exposure = exposure, snp = .x, covariates = covariates, strata = 'study_design', method = method, flip_allele = T, path = path))
walk(snps_flip, ~ fit_gxe_stratified(data_epi = input_data2, exposure = exposure, snp = .x, covariates = covariates, strata = 'cancer_site_sum2', method = method, flip_allele = T, path = path))
walk(snps_flip, ~ fit_gxe_stratified(data_epi = input_data2, exposure = exposure, snp = .x, covariates = covariates, strata = 'bmic3', method = method, flip_allele = T, path = path))







# GE findings (different stratified model in the summary tables)
snps_chiSqGE <- c("13:61572969:G:A")
method = 'chiSqGE'

walk(snps_chiSqGE, ~ fit_gxe_stratified(data_epi = input_data2, exposure = exposure, snp = .x, covariates = covariates, strata = 'sex', method = method, flip_allele = F, path = path))
walk(snps_chiSqGE, ~ fit_gxe_stratified(data_epi = input_data2, exposure = exposure, snp = .x, covariates = covariates, strata = 'study_design', method = method, flip_allele = F, path = path))
walk(snps_chiSqGE, ~ fit_gxe_stratified(data_epi = input_data2, exposure = exposure, snp = .x, covariates = covariates, strata = 'cancer_site_sum2', method = method, flip_allele = F, path = path))
walk(snps_chiSqGE, ~ fit_gxe_stratified(data_epi = input_data2, exposure = exposure, snp = .x, covariates = covariates, strata = 'bmic3', method = method, flip_allele = F, path = path))


# all the above, but exclude adenomas (+ matched controls)
input_data_crc <- input_data %>% 
  filter(is.na(adenoma)) %>% 
  mutate(active_met_875 = as.factor(ifelse(active_met_875 == "No", 0, 1)))


# GxE findings - no adenomas
method = 'chiSqGxE'
walk(snps, ~ fit_gxe_stratified(data_epi = input_data_crc, exposure = exposure, snp = .x, covariates = covariates, strata = 'sex', method = method, flip_allele = F, path = path))
walk(snps, ~ fit_gxe_stratified(data_epi = input_data_crc, exposure = exposure, snp = .x, covariates = covariates, strata = 'study_design', method = method, flip_allele = F, path = path))
walk(snps, ~ fit_gxe_stratified(data_epi = input_data_crc, exposure = exposure, snp = .x, covariates = covariates, strata = 'cancer_site_sum2', method = method, flip_allele = F, path = path))
walk(snps, ~ fit_gxe_stratified(data_epi = input_data_crc, exposure = exposure, snp = .x, covariates = covariates, strata = 'bmic3', method = method, flip_allele = F, path = path))

walk(snps_flip, ~ fit_gxe_stratified(data_epi = input_data_crc, exposure = exposure, snp = .x, covariates = covariates, strata = 'sex', method = method, flip_allele = F, path = path))
walk(snps_flip, ~ fit_gxe_stratified(data_epi = input_data_crc, exposure = exposure, snp = .x, covariates = covariates, strata = 'study_design', method = method, flip_allele = F, path = path))
walk(snps_flip, ~ fit_gxe_stratified(data_epi = input_data_crc, exposure = exposure, snp = .x, covariates = covariates, strata = 'cancer_site_sum2', method = method, flip_allele = F, path = path))
walk(snps_flip, ~ fit_gxe_stratified(data_epi = input_data_crc, exposure = exposure, snp = .x, covariates = covariates, strata = 'bmic3', method = method, flip_allele = F, path = path))



snps_chiSqGE <- c("13:61572969:G:A")
method = 'chiSqGE'
walk(snps_chiSqGE, ~ fit_gxe_stratified(data_epi = input_data_crc, exposure = exposure, snp = .x, covariates = covariates, strata = 'sex', method = method, flip_allele = F, path = path))
walk(snps_chiSqGE, ~ fit_gxe_stratified(data_epi = input_data_crc, exposure = exposure, snp = .x, covariates = covariates, strata = 'study_design', method = method, flip_allele = F, path = path))
walk(snps_chiSqGE, ~ fit_gxe_stratified(data_epi = input_data_crc, exposure = exposure, snp = .x, covariates = covariates, strata = 'cancer_site_sum2', method = method, flip_allele = F, path = path))
walk(snps_chiSqGE, ~ fit_gxe_stratified(data_epi = input_data_crc, exposure = exposure, snp = .x, covariates = covariates, strata = 'bmic3', method = method, flip_allele = F, path = path))




# snps <- c("13:61572969:G:A")
# walk(snps, ~ fit_gxe_stratified(data_epi = input_data2, exposure = exposure, snp = .x, covariates = covariates, strata = 'sex', method = 'chiSqGE', flip_allele = F, path = glue("{path}/output")))
# walk(snps, ~ fit_gxe_stratified(data_epi = input_data2, exposure = exposure, snp = .x, covariates = covariates, strata = 'study_design', method = 'chiSqGE', flip_allele = F, path = glue("{path}/output")))
# walk(snps, ~ fit_gxe_stratified(data_epi = input_data2, exposure = exposure, snp = .x, covariates = covariates, strata = 'cancer_site_sum2', method = 'chiSqGE', flip_allele = F, path = glue("{path}/output")))


# ---- models with additional covariates ---- #
# walk(snps, ~ fit_gxe_covars(data_epi = input_data, exposure = exposure, snp = .x, covariates_list = covariates_sets, method = 'chiSq3df', path = glue("{path}/output")))

# ---- RERI table ---- #
# walk(snps, ~ reri_wrapper(data_epi = input_data, exposure = exposure, snp = .x, covariates = covariates, path = glue("{path}/output")))




# ---- conditional analyses ---- #
dosages <- lapply(list.files(path = "/media/work/gwis_test/active_met_875/output/posthoc/", pattern = "dosage_chr15", full.names = T), qread)
dosages <- reduce(.x = dosages, .f = inner_join, by = 'vcfid')


snps_flip <- c("15:32994756:T:C", "15:32994777:C:T")
epi <- input_data2 %>% 
  inner_join(dosages, 'vcfid') %>% 
  mutate(active_met_875 = ifelse(is.na(methrswklns),NA,
                                 ifelse(methrswklns<8.75, "No", "Yes"))) %>% 
  mutate(chr15_32994756_C_T_dose = abs(2-chr15_32994756_T_C_dose), 
         chr15_32994777_T_C_dose = abs(2-chr15_32994777_C_T_dose))

# top finding
model1 <- glm(outcome ~ active_met_875*chr15_32994777_T_C_dose + age_ref_imp + sex + energytot_imp + pc1 + pc2 + pc3 + study_gxe, data = epi, family ='binomial')


model1a <- glm(outcome ~ active_met_875*chr15_32994777_T_C_dose + chr15_32994756_C_T_dose + age_ref_imp + sex + energytot_imp + pc1 + pc2 + pc3 + study_gxe, data = epi, family ='binomial')
model1b <- glm(outcome ~ active_met_875*chr15_32994777_T_C_dose + active_met_875*chr15_32994756_C_T_dose + age_ref_imp + sex + energytot_imp + pc1 + pc2 + pc3 + study_gxe, data = epi, family ='binomial')
summary(model1b)


model2a <- glm(outcome ~ active_met_875*chr15_32994777_T_C_dose + chr15_32995298_C_A_dose + age_ref_imp + sex + energytot_imp + pc1 + pc2 + pc3 + study_gxe, data = epi, family ='binomial')
model2b <- glm(outcome ~ active_met_875*chr15_32994777_T_C_dose + active_met_875*chr15_32995298_C_A_dose + age_ref_imp + sex + energytot_imp + pc1 + pc2 + pc3 + study_gxe, data = epi, family ='binomial')
summary(model1b)


model3a <- glm(outcome ~ active_met_875*chr15_32994777_T_C_dose + chr15_32995764_C_G_dose + age_ref_imp + sex + energytot_imp + pc1 + pc2 + pc3 + study_gxe, data = epi, family ='binomial')
model3b <- glm(outcome ~ active_met_875*chr15_32994777_T_C_dose + active_met_875*chr15_32995764_C_G_dose + age_ref_imp + sex + energytot_imp + pc1 + pc2 + pc3 + study_gxe, data = epi, family ='binomial')
summary(model1b)


model4a <- glm(outcome ~ active_met_875*chr15_32994777_T_C_dose + chr15_33001734_T_A_dose + age_ref_imp + sex + energytot_imp + pc1 + pc2 + pc3 + study_gxe, data = epi, family ='binomial')
model4b <- glm(outcome ~ active_met_875*chr15_32994777_T_C_dose + active_met_875*chr15_33001734_T_A_dose + age_ref_imp + sex + energytot_imp + pc1 + pc2 + pc3 + study_gxe, data = epi, family ='binomial')
summary(model1b)


list_of_glms <- list(model1, model1b, model2b, model3b, model4b)
out_html <- capture.output(stargazer(list_of_glms, align = T, 
          type = 'html', 
          ci = TRUE, 
          ci.level = 0.95,
          keep.stat = "n", 
          star.cutoffs = c(0.05, 0.01, 0.001), 
          # coef = coefs, 
          omit = c("pc", "study_gxe"),
          p.auto = F, 
          single.row = T,
          column.sep.width = '20pt',
          report = ('vcsp'), 
          column.labels = c("top hit", "top hit + snp1", "top hit + snp2", "top hit + snp3", "top hit + snp4")))
filename = glue("{path}/output/posthoc/conditional_analysis_{exposure}_{hrc_version}_{deparse(substitute(data_epi))}.html")
cat(paste(out_html, collapse = "\n"), "\n", file = filename, append = F)



# ================================================================== #
# Rmarkdown reports -------------------------------------------------------
# ================================================================== #
main_effects_report(exposure = exposure, hrc_version = hrc_version, covariates = covariates, path = path)

gwis_report(exposure = exposure, hrc_version = hrc_version, covariates = covariates)

rmarkdown::render(glue("/home/rak/git/figi/{exposure}_posthoc.Rmd"), 
                  params = list(exposure = exposure), 
                  output_file = glue("~/Dropbox/FIGI/Results/{exposure}_posthoc.html"))




