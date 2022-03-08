#=============================================================================#
# FIGI GxE fruit5qcm results
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
library(qs)
rm(list = ls())

# input variables
exposure = 'fruit5qcm'
exposure_scaled = 'fruit5qcm_scaled'
hrc_version = 'v3.0'
annotation_file <- 'gwas_200_ld_annotation_feb2021.txt'
path = glue("/media/work/gwis_test/{exposure}/")

covariates <- sort(c('age_ref_imp', 'sex', 'energytot_imp', 'study_gxe', 'pc1', 'pc2', 'pc3'))

# input data
esubset <- readRDS(glue("{path}/data/FIGI_{hrc_version}_gxeset_{exposure}_basic_covars_glm.rds")) %>% pull(vcfid)
input_data <- readRDS(glue("/media/work/gwis/data/FIGI_EpiData/FIGI_{hrc_version}_gxeset_analysis_data_glm.rds")) %>% 
  filter(vcfid %in% esubset) %>% 
  mutate(fruit5qcm = as.numeric(fruit5qcm))


# ---- rescale variable by median study/sex specific IQR
# ONLY USE THIS WHEN CREATING META-ANALYSIS AND POOLED ANALYSES RESULTS

# calculate iqr
input_data <- input_data %>% 
  group_by(study_gxe, sex) %>% 
  mutate(fruit5qcm.iqr = IQR(fruit5qcm)) %>% 
  ungroup()

# test <- input_data %>% 
#   group_by(study_gxe, sex, fruit5qcm.iqr) %>% 
#   summarise(n = n())

# hist(input_data$fruit5qcm.iqr)
summary(input_data$fruit5qcm.iqr)

# scale by median sex/study specific IQR
# after scaling, center variable so that ref (0) = median value, 1 = +1 IQR
fruit5qcm.iqr.median = median(input_data$fruit5qcm.iqr)

input_data <- input_data %>% 
  mutate(fruit5qcm = fruit5qcm / fruit5qcm.iqr.median, 
         fruit5qcm_scaled = (fruit5qcm - median(fruit5qcm)) / IQR(fruit5qcm)) %>% 
  data.frame()



# test dachs_1
# dachs1 <- filter(input_data, study_gxe == "DACHS_1")
# table(dachs1$fruit5qcm_scaled)


# # let me try something
# # for meta-analyses, i'm going to rescale overall, such that ORs reflect increase in 1SD overall (based on median coded values in all studies)
# input_data <- readRDS(glue("/media/work/gwis/data/FIGI_EpiData/FIGI_{hrc_version}_gxeset_analysis_data_glm.rds")) %>% 
#   filter(vcfid%in% esubset) %>% 
#   # group_by(study_gxe) %>% 
#   mutate(fruit5qcm = as.numeric(scale(as.numeric(fruit5qcm)))) %>% data.frame()


# FYI, centering does not change estimates 
glm(outcome ~ fruit5qcm + age_ref_imp + sex, data = input_data, family = 'binomial')

glm(outcome ~ fruit5qcm_scaled + age_ref_imp + sex, data = input_data, family = 'binomial')


#-----------------------------------------------------------------------------#
# main effects ----
#-----------------------------------------------------------------------------#

# ------ meta-analysis ------ #
output_dir = as.character(glue("{path}/output/posthoc/"))
covariates_meta <- sort(covariates[which(!covariates %in% c(paste0(rep('pc', 20), seq(1, 20)), "study_gxe"))])
covariates_meta <- sort(covariates[which(!covariates %in% c("study_gxe"))])

create_forest_plot(data_epi = input_data, exposure = exposure_scaled, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "all", forest_height = 15, categorical = F)
create_forest_plot(data_epi = input_data, exposure = exposure_scaled, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "proximal", forest_height = 13, categorical = F)
create_forest_plot(data_epi = input_data, exposure = exposure_scaled, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "distal", forest_height = 13, categorical = F)
create_forest_plot(data_epi = input_data, exposure = exposure_scaled, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "rectal", forest_height = 13, categorical = F)
create_forest_plot(data_epi = input_data, exposure = exposure_scaled, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "female", forest_height = 13, categorical = F)
create_forest_plot(data_epi = input_data, exposure = exposure_scaled, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "male", forest_height = 13, categorical = F)


# ------- stratified pooled analysis ------- #
pooled_analysis_glm(input_data, exposure = exposure_scaled, hrc_version = hrc_version, covariates = covariates, strata = 'sex', filename_suffix = paste0(paste0(sort(covariates), collapse = '_'), "_stratified_sex"), output_dir = glue("{path}/output/posthoc/"))

pooled_analysis_glm(input_data, exposure = exposure_scaled, hrc_version = hrc_version, covariates = covariates, strata = 'study_design', filename_suffix = paste0(paste0(sort(covariates), collapse = '_'), "_stratified_study_design"), output_dir = glue("{path}/output/posthoc/"))

pooled_analysis_multinom(input_data, exposure = exposure_scaled, hrc_version = hrc_version, covariates = covariates, filename_suffix = paste0(paste0(sort(covariates), collapse = '_'), "_stratified_cancer_site_sum2"), output_dir = glue("{path}/output/posthoc/"))






#-----------------------------------------------------------------------------#
# main effects additional ----
#-----------------------------------------------------------------------------#

output_dir_dropbox = paste0("~/Dropbox/Presentations/", exposure, "/")

tmp1 <- readRDS(paste0("/media/work/gwis/results/input/FIGI_", hrc_version, "_gxeset_", exposure, "_basic_covars_glm.rds")) %>% 
  pull(vcfid)
xx <- filter(input_data, vcfid %in% tmp1) %>% 
  group_by(study_gxe)


results_beta <- dplyr::do(xx, broom::tidy(glm(outcome ~ folate_dietqc2 + age_ref_imp + sex + energytot_imp, data = . , family = 'binomial'))) %>% 
  dplyr::filter(grepl("folate_dietqc2", term)) %>% 
  dplyr::arrange(study_gxe) %>% 
  inner_join(unique(xx[,c('study_gxe', 'study_design')]), 'study_gxe')

results_meta <- meta::metagen(estimate,
                              std.error,
                              data=results_beta,
                              studlab=paste(study_gxe),
                              comb.fixed = FALSE,
                              comb.random = TRUE,
                              method.tau = "SJ",
                              hakn = TRUE,
                              prediction=TRUE,
                              sm="OR", 
                              byvar=study_design)

fo <- find.outliers(results_meta)
fo
meta::forest(results_meta,
             layout = "JAMA",
             text.predict = "95% CI",
             col.predict = "black",
             # leftcols = c("studlab", "Control", "Case", "N", "effect", "ci", "w.random"),
             digits.addcols=0,
             study.results=T,
             prediction = F,
             col.random = 'red')

png(paste0(output_dir_dropbox, "meta_analysis_", "folate_dietqc2",  "_", "original_outliers_removed", ".png"), height = 17, width = 8.5, units = 'in', res = 150)                                                          
forest(fo)
dev.off()


# leave on out (influence analysis)
inf.analysis <- InfluenceAnalysis(x = results_meta,
                                  random = TRUE)

summary(inf.analysis)

plot(inf.analysis, "influence")
plot(inf.analysis, "baujat")
plot(inf.analysis, "es")
plot(inf.analysis, "i2")



#-----------------------------------------------------------------------------#
# GxE additional analysis ---- 
#-----------------------------------------------------------------------------#


# ---- AAF by study_gxe ---- #
snps_out <- c("8:117689668:T:C", "8:117742957:A:C", "8:117752713:C:T", "1:72729142:A:G")
walk(snps_out, ~ create_aaf_study_plot(data = input_data, exposure,  hrc_version, snp = .x, path = path))

# ---- RERI plots (can't install package on CARC) ---- #
input_data <- input_data %>% 
  mutate(fruit5qcm = as.numeric(fruit5qcm))


snps <- c("14:74029409:C:T", "14:74029049:G:C")
snps_out <- c("8:117689668:T:C", "8:117742957:A:C", "8:117752713:C:T", "1:72729142:A:G")

walk(snps, ~ reri_wrapper(data_epi = input_data, exposure = exposure, snp = .x, covariates = covariates, path = glue("{path}/output")))


# Gxe models stratified by sex, tumor site, study_design
walk(snps_out, ~ fit_gxe_stratified(data_epi = input_data, exposure = exposure_scaled, snp = .x, covariates = covariates, strata = 'sex', method = 'chiSqGxE', path = glue("{path}/output")))
walk(snps_out, ~ fit_gxe_stratified(data_epi = input_data, exposure = exposure_scaled, snp = .x, covariates = covariates, strata = 'cancer_site_sum2', method = 'chiSqGxE', path = glue("{path}/output")))
walk(snps_out, ~ fit_gxe_stratified(data_epi = input_data, exposure = exposure_scaled, snp = .x, covariates = covariates, strata = 'study_design', method = 'chiSqGxE', path = glue("{path}/output")))


# ---- iplot ---- #
walk(snps_out, ~iplot_wrapper(data_epi = input_data, exposure = exposure_scaled, hrc_version = hrc_version, snp = .x, covariates = covariates, path = glue("{path}/output/"), flip_allele = F))


# ---- stratified odds ratios ---- #
walk(snps_out, ~ fit_stratified_or_continuous(data_epi = input_data, exposure = exposure_scaled, hrc_version = hrc_version, snp = .x, covariates = covariates, dosage = F, path = glue("{path}/output/")))

# ---- stratified odds ratios ---- #
walk(snps_out, ~ fit_stratified_or_q4(data_epi = input_data, exposure = exposure, hrc_version = hrc_version, snp = .x, covariates = covariates, dosage = F, path = glue("{path}/")))



# ---- GxE models by covariate sets ---- #
# covariates_sets <- list(covariates, 
#                         c(covariates, 'bmi5'), 
#                         c(covariates, 'bmi5', 'smk_ever'), 
#                         c(covariates, 'bmi5', 'smk_ever', 'fruitqc2', 'vegetableqc2'))
# covariates_sets <- list(covariates)
# walk(snps, ~ fit_gxe_covars(data_epi = input_data, exposure = exposure, snp = .x, covariates_list = covariates_sets, method = 'chiSqGxE', path = glue("{path}/output")))
# suggestive_gxe <- fread(glue("{path}/data/FIGI_v2.3_gxeset_diab_chiSqGxE_ldclump.clumped"))
# walk(suggestive_gxe$SNP, ~ fit_gxe_covars(data_epi = input_data, exposure = exposure, snp = .x, covariates_list = covariates_sets, method = 'chiSqGxE', path = glue("{path}/output")))


# ================================================================== #
# ======= rmarkdown reports ---- 
# ================================================================== #
main_effects_report(exposure = exposure, hrc_version = hrc_version, covariates = covariates, path = path)

rmarkdown::render("~/git/figi/main_effects/main_effects_scaled.Rmd", 
                    params = list(exposure = exposure, hrc_version = hrc_version, 
                                  covariates = covariates, path = path), output_file = glue("~/Dropbox/FIGI/Results/{exposure}_scaled_main_effects.html"))


gwis_report(exposure = exposure, 
            hrc_version = hrc_version, 
            covariates = covariates)

posthoc_report(exposure = exposure)



# ----- suggestive markers ------ # 
rmarkdown::render(glue("/home/rak/git/figi/{exposure}_posthoc_suggestive.Rmd"), 
                  params = list(exposure = exposure, hrc_version = hrc_version, 
                                covariates = covariates, path = path), 
                  output_file = glue("~/Dropbox/FIGI/Results/{exposure}_posthoc_suggestive.html"))



# ================================================================== #
# ======= explore ---- 
# ================================================================== #


input_data2 <- readRDS(glue("/media/work/gwis/data/FIGI_EpiData/FIGI_{hrc_version}_gxeset_analysis_data_glm.rds")) %>% 
  filter(vcfid%in% esubset) %>% 
  mutate(fruit5qcm = scale(as.numeric(fruit5qcm)))


# does individual study values change if i rescale the overall poppulation ... 

dachs1a <- filter(input_data, study_gxe == "DACHS_1")
dachs1b <- filter(input_data2, study_gxe == "DACHS_1") 

table(dachs1a$fruit5qcm)
table(dachs1b$fruit5qcm)

summary(input_data$fruit5qcm)

model1 <- glm(outcome ~ fruit5qcm + age_ref_imp + sex + energytot_imp + pc1 + pc2 + pc3 + study_gxe, data = input_data, family = 'binomial')
exp(coef(model1))

dachs1 <-

table(dachs1$outcome, dachs1$fruit5qcm)
table(dachs1$outcome, dachs1$test)

model1 <- glm(outcome ~ fruit5qcm + age_ref_imp + sex + energytot_imp + pc1 + pc2 + pc3, data = dachs1, family = 'binomial')
exp(coef(model1))

model1 <- glm(outcome ~ fruit5qcm + age_ref_imp + sex, data = dachs1, family = 'binomial')
tidy(model1, conf.int = T, exponentiate = T)

model1 <- glm(outcome ~ fruit5qcm + age_ref_imp + sex + energytot_imp, data = dachs1 %>% mutate(fruit5qcm = fruit5qcm * 10), family = 'binomial')
tidy(model1, conf.int = T, exponentiate = T)

model2 <- glm(outcome ~ test + age_ref_imp + sex, data = dachs1, family = 'binomial')
tidy(model2, conf.int = T, exponentiate = T)
