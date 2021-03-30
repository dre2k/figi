#=============================================================================#
# FIGI GxE vitd25 results
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
exposure = 'vitd25'
hrc_version = 'v2.3'
output_dir = paste0("/media/work/gwis/posthoc/", exposure, "/")
annotation_file <- 'gwas_141_ld_annotation_july2020.txt'

covariates <- sort(c('age_ref_imp', 'sex', 'study_gxe', 'pc1', 'pc2', 'pc3'))
covariates_set1 <- covariates
# covariates_set2 <- sort(c(covariates_set1, 'bmi', 'smk_ever', 'fruitqc2', 'vegetableqc2'))
covariates_list <- list(covariates)
mod <- 'age_ref_imp+sex+pc1+pc2+pc3+study_gxe'

# input data (need to read in another file to get the vitd25 variable)
# also keep in mind that when you eventually adjust for batch or season things will get more complicated
esubset <- readRDS(glue("/media/work/gwis/data/FIGI_EpiData/FIGI_{hrc_version}_gxeset_{exposure}_basic_covars_glm.rds")) %>% 
  pull(vcfid)
input_data <- readRDS(glue("/media/work/gwis/data/FIGI_EpiData/FIGI_{hrc_version}_gxeset_analysis_data_glm.rds")) %>% 
  filter(vcfid%in% esubset)

vitd <- readRDS("/media/work/gwis/data/FIGI_EpiData/FIGI_v2.3_gxeset_vitd25_basic_covars_glm.rds")

input_data$vitd25 <- vitd$vitd25





#-----------------------------------------------------------------------------#
# main effects ----
#-----------------------------------------------------------------------------#

# ------ meta-analysis ------ #
covariates_meta <- sort(c('age_ref_imp', 'sex'))
meta_analysis_execute(dataset = input_data, exposure = exposure, hrc_version = hrc_version, covariates = covariates_meta, output_dir = output_dir, filename_suffix2 = "")



# ------- pooled analysis ------- #
## basic covariates
covariates_pooled <- sort(covariates[which(!covariates %in% c(paste0(rep('pc', 20), seq(1, 20))))])
covariates_pooled <- covariates
pooled_analysis_glm(input_data, exposure = exposure, covariates = covariates_pooled, strata = 'sex', 
                    filename_suffix = paste0(paste0(sort(covariates), collapse = '_'), "_stratified_sex"), output_dir = output_dir)
pooled_analysis_glm(input_data, exposure = exposure, covariates = covariates_pooled, strata = 'study_design', 
                    filename_suffix = paste0(paste0(sort(covariates), collapse = '_'), "_stratified_study_design"), output_dir = output_dir)
pooled_analysis_multinom(input_data, exposure = exposure, covariates = covariates_pooled, 
                         filename_suffix = paste0(paste0(sort(covariates), collapse = '_'), "_stratified_cancer_site_sum2"), output_dir = output_dir)

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

# output GxE models adjusted by different covariate sets
covariates_sets <- list(covariates, 
                        c(covariates, 'bmi5'), 
                        c(covariates, 'bmi5', 'smk_ever'), 
                        c(covariates, 'bmi5', 'smk_ever', 'fruitqc2', 'vegetableqc2'))

fit_gxe_covars(data_epi = input_data, exposure = exposure, snp = "13:47191972:G:A", covariates_list = covariates_sets, method = 'chiSqGxE', path = glue("{path}/output"))


# additional covariates is making association more significant.. let's generate for all suggestive hits
suggestive_gxe <- fread(glue("{path}/data/FIGI_v2.3_gxeset_diab_chiSqGxE_ldclump.clumped"))
walk(suggestive_gxe$SNP, ~ fit_gxe_covars(data_epi = input_data, exposure = exposure, snp = .x, covariates_list = covariates_sets, method = 'chiSqGxE', path = glue("{path}/output")))


# output RERI plots (can't install package on CARC yet)
reri_wrapper(data_epi = input_data, exposure = exposure, snp = "13:47191972:G:A", covariates = covariates, path = glue("{path}/output"))




# ================================================================== #
# ======= rmarkdown reports ---- 
# ================================================================== #
# 
gwis_report(exposure = exposure, hrc_version = hrc_version, covariates = covariates)

posthoc_report(exposure = exposure)



