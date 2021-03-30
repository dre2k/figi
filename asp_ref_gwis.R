#=============================================================================#
# FIGI GxE bmi5 results
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
exposure = 'asp_ref'
hrc_version = 'v2.4'
annotation_file <- 'gwas_200_ld_annotation_feb2021.txt'
covariates <- sort(c('age_ref_imp', 'sex', 'study_gxe', 'pc1', 'pc2', 'pc3'))
path = glue("/media/work/gwis_test/{exposure}/")


# input data
esubset <- readRDS(glue("/media/work/gwis_test/{exposure}/data/FIGI_{hrc_version}_gxeset_{exposure}_basic_covars_glm.rds")) %>% 
  pull(vcfid)

input_data <- readRDS(glue("/media/work/gwis_test/data/FIGI_{hrc_version}_gxeset_analysis_data_glm.rds")) %>% 
  filter(vcfid%in% esubset)



#-----------------------------------------------------------------------------#
# main effects ----
#-----------------------------------------------------------------------------#

# ------ meta-analysis ------ #
output_dir = as.character(glue("{path}/output/posthoc/"))
covariates_meta <- sort(covariates[which(!covariates %in% c(paste0(rep('pc', 20), seq(1, 20)), "study_gxe"))])

create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "all", forest_height = 15)
create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "proximal", forest_height = 13)
create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "distal", forest_height = 13)
create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "rectal", forest_height = 13)
create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "female", forest_height = 13)
create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "male", forest_height = 13)


# ------- stratified pooled analysis ------- #
pooled_analysis_glm(input_data, exposure = exposure, hrc_version = hrc_version, covariates = covariates, strata = 'sex', filename_suffix = paste0(paste0(sort(covariates), collapse = '_'), "_stratified_sex"), output_dir = glue("{path}/output/posthoc/"))

pooled_analysis_glm(input_data, exposure = exposure, hrc_version = hrc_version, covariates = covariates, strata = 'study_design', filename_suffix = paste0(paste0(sort(covariates), collapse = '_'), "_stratified_study_design"), output_dir = glue("{path}/output/posthoc/"))

pooled_analysis_multinom(input_data, exposure = exposure, hrc_version = hrc_version, covariates = covariates, filename_suffix = paste0(paste0(sort(covariates), collapse = '_'), "_stratified_cancer_site_sum2"), output_dir = glue("{path}/output/posthoc/"))





# main effects additional ----
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
# functional annotation subset ---- 
# focus on pooled scores for now
#-----------------------------------------------------------------------------#

svm_pooled <- readRDS("/media/work/svm_scores/svm_pooled_filter_sd3.rds")
x1 <- gxe %>%
  filter(SNP %in% svm_pooled$SNP)


# output bin SNPs for expectation based hybrid method..  
twostep_eh_snps(x1, 'chiSqG', output_dir = glue('/media/work/gwis/twostep_expectation_hybrid/{exposure}/svm_subset/'))
twostep_eh_snps(x1, 'chiSqGE', output_dir = glue('/media/work/gwis/twostep_expectation_hybrid/{exposure}/svm_subset/'))
twostep_eh_snps(x1, 'chiSqEDGE', output_dir = glue('/media/work/gwis/twostep_expectation_hybrid/{exposure}/svm_subset/'))


# create manhattan, qq, two-step
plot_funcannot_wrap(gxe, exposure = exposure, covariates = covariates, output_dir = output_dir, filename_suffix = "_functional_subset")

# run expectation based hybrid?
simplem_wrap2(x = x1, exposure = exposure, covariates = covariates, simplem_step1_statistic = 'chiSqG', output_dir = output_dir, filename_suffix = "_functional_subset")
simplem_wrap2(x = x1, exposure = exposure, covariates = covariates, simplem_step1_statistic = 'chiSqGE', output_dir = output_dir, filename_suffix = "_functional_subset")
simplem_wrap2(x = x1, exposure = exposure, covariates = covariates, simplem_step1_statistic = 'chiSqEDGE', output_dir = output_dir, filename_suffix = "_functional_subset")









#-----------------------------------------------------------------------------#
# GxE additional analysis ---- 
#-----------------------------------------------------------------------------#

# output GxE models adjusted by different covariate sets
covariates_sets <- list(covariates)


# stratified by tumor and sex 
snps <- c("6:12577203:T:C", "6:32560631:C:T", "5:40252294:C:T")
walk(snps, ~ fit_gxe_stratified(data_epi = input_data, exposure = exposure, snp = .x, covariates = covariates, method = "chiSqGxE", strata = 'sex', path = glue("{path}/output")))
walk(snps, ~ fit_gxe_stratified(data_epi = input_data, exposure = exposure, snp = .x, covariates = covariates, method = "chiSqGxE", strata = 'cancer_site_sum2', path = glue("{path}/output")))


# output RERI plots
snps <- c("6:12577203:T:C")
walk(snps, ~ reri_wrapper(data_epi = input_data, exposure = exposure, snp = .x, covariates = covariates, path = glue("{path}/output")))

# need to flip exposure
tmp_flip <- input_data %>% 
  mutate(asp_ref = fct_relevel(asp_ref, "Yes"))
snps <- c( "6:32560631:C:T", "5:40252294:C:T")
walk(snps, ~ reri_wrapper(data_epi = tmp_flip, exposure = exposure, snp = .x, covariates = covariates, path = glue("{path}/output")))





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




# ================================================================== #
# updated PCs and how it changes results.. 
# ================================================================== #

# original GxE findings.. 
modelb <- glm(outcome ~ chr6_12577203_T_C+asp_ref + age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe, data = figi, family = 'binomial')
model1 <- glm(outcome ~ chr6_12577203_T_C*asp_ref + age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe, data = figi, family = 'binomial')
lrtest(modelb, model1)



# let's add updated PCs .. (keep in mind they're calculated using HRC v3.0)
pcs <- fread("~/figi_gxe_pca_update.eigenvec")
figi_newpc <- inner_join(figi, pcs, by = c('vcfid' = 'IID'))


modelb <- glm(outcome ~ chr6_12577203_T_C+asp_ref + age_ref_imp + sex + PC1 + PC2 + PC3 + study_gxe, data = figi_newpc, family = 'binomial')
model1 <- glm(outcome ~ chr6_12577203_T_C*asp_ref + age_ref_imp + sex + PC1 + PC2 + PC3 + study_gxe, data = figi_newpc, family = 'binomial')
summary(model1)

lrtest(modelb, model1)

modelb <- glm(outcome ~ chr6_12577203_T_C+asp_ref + age_ref_imp + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + study_gxe, data = figi_newpc, family = 'binomial')
model1 <- glm(outcome ~ chr6_12577203_T_C*asp_ref + age_ref_imp + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10+ study_gxe, data = figi_newpc, family = 'binomial')
summary(model1)

lrtest(modelb, model1)




#  how about the other hit..
# original GxE findings.. 
  modelb <- glm(outcome ~ chr6_32560631_C_T+asp_ref + age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe, data = figi, family = 'binomial')
model1 <- glm(outcome ~ chr6_32560631_C_T*asp_ref + age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe, data = figi, family = 'binomial')
lrtest(modelb, model1)

# let's add updated PCs .. (keep in mind they're calculated using HRC v3.0)
modelb <- glm(outcome ~ chr6_32560631_C_T+asp_ref + age_ref_imp + sex + PC1 + PC2 + PC3 + study_gxe, data = figi_newpc, family = 'binomial')
model1 <- glm(outcome ~ chr6_32560631_C_T*asp_ref + age_ref_imp + sex + PC1 + PC2 + PC3 + study_gxe, data = figi_newpc, family = 'binomial')
summary(model1)

lrtest(modelb, model1)




# chr5_40252294_C_T finding (two step expectation based.. not a 1:1 comparison but just to get an idea)
wtf <- filter(gxe , Location == 40252294) # EDGE chi-square = 24.32772

figi$asp_ref2 <- as.numeric(figi$asp_ref)
figi_newpc$asp_ref2 <- as.numeric(figi_newpc$asp_ref)



# marginal association
modelb <- glm(outcome ~                     asp_ref + age_ref_imp + sex + PC1 + PC2 + PC3 + study_gxe, data = figi_newpc, family = 'binomial')
model1 <- glm(outcome ~ chr5_40252294_C_T + asp_ref + age_ref_imp + sex + PC1 + PC2 + PC3 + study_gxe, data = figi_newpc, family = 'binomial')
lrtest(modelb, model1)
# 24.275 


modelb <- lm(asp_ref2 ~                     age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe, data = figi_newpc)
model1 <- lm(asp_ref2 ~ chr5_40252294_C_T + age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe, data = figi_newpc)
lrtest(modelb, model1)
# 0.5652

pchisq(24.8402, lower.tail = F, df = 2)




# GxE
modelb <- glm(outcome ~ chr5_40252294_C_T+asp_ref + age_ref_imp + sex + PC1 + PC2 + PC3 + study_gxe, data = figi_newpc, family = 'binomial')
model1 <- glm(outcome ~ chr5_40252294_C_T*asp_ref + age_ref_imp + sex + PC1 + PC2 + PC3 + study_gxe, data = figi_newpc, family = 'binomial')
lrtest(modelb, model1)
