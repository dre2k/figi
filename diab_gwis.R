#=============================================================================#
# FIGI GxE diab results
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
rm(list = ls())

# input variables
exposure = 'diab'
hrc_version = 'v2.3'
# output_dir = paste0("/media/work/gwis/posthoc/", exposure, "/")
# annotation_file <- 'gwas_141_ld_annotation_july2020.txt'
annotation_file <- 'gwas_200_ld_annotation_feb2021.txt'
path = glue("/media/work/gwis_test/{exposure}/")


covariates <- sort(c('age_ref_imp', 'sex', 'study_gxe', 'pc1', 'pc2', 'pc3'))
covariates_set1 <- covariates
# covariates_set2 <- sort(c(covariates_set1, 'bmi', 'smk_ever', 'fruitqc2', 'vegetableqc2'))
covariates_list <- list(covariates)
mod <- 'age_ref_imp+sex+pc1+pc2+pc3+study_gxe'

# input data
esubset <- readRDS(glue("{path}/data/FIGI_{hrc_version}_gxeset_{exposure}_basic_covars_glm.rds")) %>% pull(vcfid)
input_data <- readRDS(glue("/media/work/gwis/data/FIGI_EpiData/FIGI_{hrc_version}_gxeset_analysis_data_glm.rds")) %>% 
  filter(vcfid%in% esubset)




#-----------------------------------------------------------------------------#
# main effects ----
#-----------------------------------------------------------------------------#

# ------- meta-analysis (manuscript)------ # 
output_dir = as.character(glue("{path}/output/posthoc/"))
covariates_meta <- sort(covariates[which(!covariates %in% c(paste0(rep('pc', 20), seq(1, 20)), "study_gxe"))])
# covariates_meta <- sort(covariates[which(!covariates %in% c("study_gxe"))])


create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "all", forest_height = 15, forest_width = 9)
create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "proximal", forest_height = 13, forest_width = 9)
create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "distal", forest_height = 13, forest_width = 9)
create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "rectal", forest_height = 13, forest_width = 9)
create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "female", forest_height = 13, forest_width = 9)
create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "male", forest_height = 13, forest_width = 9)




# ------- meta-analysis (manuscript) + bmi------ # 
output_dir = as.character(glue("{path}/output/posthoc/"))
covariates_meta <- sort(covariates[which(!covariates %in% c(paste0(rep('pc', 20), seq(1, 20)), "study_gxe"))])
covariates_meta <- c(covariates_meta, 'bmi')

input_data2 <- input_data %>% 
  dplyr::filter(!is.na(bmi))
  
create_forest_plot(data_epi = input_data2, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "all", forest_height = 15, forest_width = 9)
create_forest_plot(data_epi = input_data2, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "proximal", forest_height = 13, forest_width = 9)
create_forest_plot(data_epi = input_data2, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "distal", forest_height = 13, forest_width = 9)
create_forest_plot(data_epi = input_data2, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "rectal", forest_height = 13, forest_width = 9)
create_forest_plot(data_epi = input_data2, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "female", forest_height = 13, forest_width = 9)
create_forest_plot(data_epi = input_data2, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "male", forest_height = 13, forest_width = 9)



# ------- meta-analysis (manuscript) + bmi + smk_ever + redmeatqc2------ # 
output_dir = as.character(glue("{path}/output/posthoc/"))
covariates_meta <- sort(covariates[which(!covariates %in% c(paste0(rep('pc', 20), seq(1, 20)), "study_gxe"))])
covariates_meta <- c(covariates_meta, 'bmi', 'smk_ever', 'redmeatqc2')

input_data2 <- input_data %>% 
  dplyr::filter(!is.na(bmi), 
                !is.na(smk_ever), 
                !is.na(redmeatqc2))

create_forest_plot(data_epi = input_data2, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "all", forest_height = 15, forest_width = 9)
create_forest_plot(data_epi = input_data2, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "proximal", forest_height = 13, forest_width = 9)
create_forest_plot(data_epi = input_data2, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "distal", forest_height = 13, forest_width = 9)
create_forest_plot(data_epi = input_data2, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "rectal", forest_height = 13, forest_width = 9)
create_forest_plot(data_epi = input_data2, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "female", forest_height = 13, forest_width = 9)
create_forest_plot(data_epi = input_data2, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "male", forest_height = 13, forest_width = 9)





# ------ meta-analysis (updated?) ------ #
output_dir = as.character(glue("{path}/posthoc/"))
covariates_meta <- sort(covariates[which(!covariates %in% c(paste0(rep('pc', 20), seq(1, 20)), "study_gxe"))])

create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("~/Dropbox/"), strata = "all", forest_height = 15, categorical = F)



# ------ meta-analysis ------ #
covariates_meta <- sort(c('age_ref_imp', 'sex'))
meta_analysis_execute(dataset = input_data, exposure = exposure, hrc_version = hrc_version, covariates = covariates_meta, output_dir = output_dir, filename_suffix2 = "")

# BMI adjustment

tmp <- input_data %>% 
  filter(!is.na(bmi))

covariates_meta <- sort(c('age_ref_imp', 'sex', 'bmi'))
meta_analysis_execute(dataset = tmp, exposure = exposure, hrc_version = hrc_version, covariates = covariates_meta, output_dir = output_dir, filename_suffix2 = "")


tmp <- input_data %>% 
  filter(!is.na(bmi), 
         !is.na(redmeatqc2), 
         !is.na(smk_ever))

covariates_meta <- sort(c('age_ref_imp', 'sex', 'bmi', 'smk_ever', 'redmeatqc2'))
meta_analysis_execute(dataset = tmp, exposure = exposure, hrc_version = hrc_version, covariates = covariates_meta, output_dir = output_dir, filename_suffix2 = "")




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
# main effects additional analysis ----
#-----------------------------------------------------------------------------#

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
# plot(inf.analysis, "es")
# plot(inf.analysis, "i2")





#-----------------------------------------------------------------------------#
# functional annotation subset ---- 
# focus on pooled scores for now
#-----------------------------------------------------------------------------#

# svm_pooled <- readRDS("/media/work/svm_scores/svm_pooled_filter_sd3.rds")
# x1 <- gxe %>%
#   filter(SNP %in% svm_pooled$SNP)
# 
# 
# # output bin SNPs for expectation based hybrid method..  
# twostep_eh_snps(x1, 'chiSqG', output_dir = glue('/media/work/gwis/twostep_expectation_hybrid/{exposure}/svm_subset/'))
# twostep_eh_snps(x1, 'chiSqGE', output_dir = glue('/media/work/gwis/twostep_expectation_hybrid/{exposure}/svm_subset/'))
# twostep_eh_snps(x1, 'chiSqEDGE', output_dir = glue('/media/work/gwis/twostep_expectation_hybrid/{exposure}/svm_subset/'))
# 
# 
# # create manhattan, qq, two-step
# plot_funcannot_wrap(gxe, exposure = exposure, covariates = covariates, output_dir = output_dir, filename_suffix = "_functional_subset")
# 
# # run expectation based hybrid?
# simplem_wrap2(x = x1, exposure = exposure, covariates = covariates, simplem_step1_statistic = 'chiSqG', output_dir = output_dir, filename_suffix = "_functional_subset")
# simplem_wrap2(x = x1, exposure = exposure, covariates = covariates, simplem_step1_statistic = 'chiSqGE', output_dir = output_dir, filename_suffix = "_functional_subset")
# simplem_wrap2(x = x1, exposure = exposure, covariates = covariates, simplem_step1_statistic = 'chiSqEDGE', output_dir = output_dir, filename_suffix = "_functional_subset")




#-----------------------------------------------------------------------------#
# GxE additional analysis ---- 
#-----------------------------------------------------------------------------#

# ---- MAF ---- #
snps <- c("13:47191972:G:A", "8:118185025:G:A")
walk(snps, ~ create_aaf_study_plot(data = input_data, exposure = exposure, hrc_version = hrc_version, snp = .x, path = path))







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
reri_wrapper(data_epi = input_data, exposure = exposure, snp = "8:118185025:G:A", covariates = covariates, path = glue("{path}/output"))


# stratified odds ratio (newer) -------------------------------------------
walk(snps, ~ fit_stratified_or(
  data_epi = input_data,
  exposure = exposure,
  snp = .x,
  hrc_version = hrc_version,
  covariates = covariates,
  dosage = T,
  path = glue("{path}/output/")))


snps <- c("8:118185025:G:A")
walk(snps, ~ fit_stratified_or(
  data_epi = input_data %>% filter(!study_gxe %in% c("PPS3", "PPS4")),
  exposure = exposure,
  snp = .x,
  hrc_version = hrc_version,
  covariates = covariates,
  dosage = T, flip_allele = T,
  path = glue("{path}/output/")))




# stratified odds ratio (older) -------------------------------------------

# stratified odds ratios
snps <- c("13:47191972:G:A")

walk(snps, ~ fit_stratified_or(
  data_epi = input_data %>% filter(!study_gxe %in% c("PPS3", "PPS4")),
  exposure = 'diab',
  snp = .x,
  hrc_version = hrc_version,
  covariates = covariates,
  path = glue("{path}/output/")))


walk(snps, ~ fit_stratified_or_old(
  data_epi = input_data,
  exposure = 'diab',
  snp = .x,
  hrc_version = hrc_version,
  covariates = covariates,
  path = glue("{path}/output/")))



# Explore removing PPS3/4 studies from manuscript -------------------------
snps <- c("13:47191972:G:A", "8:118185025:G:A")

## 13:47191972:G:A (see table 1 in Niki's manuscript)
## original p-values match, go ahead and filter PPS3/4
dos <- qread("/media/work/gwis_test/diab/output/posthoc/dosage_chr13_47191972.qs")
out <- inner_join(input_data, dos, 'vcfid') %>% 
  mutate(diab_num = as.numeric(diab)) %>% 
  dplyr::filter(!study_gxe  %in% c("PPS3", "PPS4"))

model_base <- glm(outcome ~ diab + chr13_47191972_G_A_dose + age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe, data = out, family = 'binomial')
model_gxe <- glm(outcome ~ diab * chr13_47191972_G_A_dose + age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe, data = out, family = 'binomial')
model_gwas <- glm(outcome ~ diab + age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe, data = out, family = 'binomial')
model_base_eg <- lm(diab_num ~ age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe, data = out)
model_eg <- lm(diab_num ~ chr13_47191972_G_A_dose + age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe, data = out)

lrtest(model_base, model_gwas)
lrtest(model_base_eg, model_eg)
lrtest(model_base, model_gxe)
lrtest(model_gwas, model_gxe)

model_gxe <- glm(outcome ~ diab * chr13_47191972_G_A_dose + age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe + bmi, data = out, family = 'binomial')
model_gwas <- glm(outcome ~ diab + age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe + bmi, data = out, family = 'binomial')
lrtest(model_gwas, model_gxe)







## 8:118185025:G:A (see table 1 in Niki's manuscript)
## original p-values match, go ahead and filter PPS3/4
dos <- qread("/media/work/gwis_test/diab/output/posthoc/dosage_chr8_118185025.qs")
out <- inner_join(input_data, dos, 'vcfid') %>% 
  mutate(diab_num = as.numeric(diab)) %>% 
  dplyr::filter(!study_gxe  %in% c("PPS3", "PPS4"))

model_base <- glm(outcome ~ diab + chr8_118185025_G_A_dose + age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe, data = out, family = 'binomial')
model_gxe <- glm(outcome ~ diab * chr8_118185025_G_A_dose + age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe, data = out, family = 'binomial')
model_gwas <- glm(outcome ~ diab + age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe, data = out, family = 'binomial')
model_base_eg <- lm(diab_num ~ age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe, data = out)
model_eg <- lm(diab_num ~ chr8_118185025_G_A_dose + age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe, data = out)

lrtest(model_base, model_gwas)
lrtest(model_base_eg, model_eg)
lrtest(model_base, model_gxe)

# 3df ..
a<-lrtest(model_base_eg, model_eg)
b<-lrtest(model_gwas, model_gxe)
pchisq(a$Chisq + b$Chisq, df = 3, lower.tail = F)

model_gxe <- glm(outcome ~ diab * chr8_118185025_G_A_dose + age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe + bmi, data = out, family = 'binomial')
model_gwas <- glm(outcome ~ diab + age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe + bmi, data = out, family = 'binomial')
model_base_eg <- lm(diab_num ~ age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe + bmi, data = out)
model_eg <- lm(diab_num ~ chr8_118185025_G_A_dose + age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe + bmi, data = out)
a<-lrtest(model_base_eg, model_eg)
b<-lrtest(model_gwas, model_gxe)
pchisq(a$Chisq + b$Chisq, df = 3, lower.tail = F)





# ================================================================== #
# ======= interaction stratified by diabetes status ---- 
# i forget the rationale for this analysis, maybe it's to check D|G association? 
# ================================================================== #

snps_filter <- snps %>% 
  dplyr::filter(method == glue("manhattan_chiSqCase_{exposure}_clump_df")) %>% 
  dplyr::mutate(snps = paste0('chr', gsub("\\:", "\\_", SNP))) %>% 
  dplyr::arrange(Pval) %>% 
  pull(snps)



# ================================================================== #
# ======= rmarkdown reports ---- 
# ================================================================== #

gwis_report(exposure = exposure, hrc_version = hrc_version, covariates = covariates)

posthoc_report(exposure = exposure)



main_effects_report(exposure = exposure, hrc_version = hrc_version, covariates = covariates, path = path)


gwis_report(exposure = exposure, 
            hrc_version = hrc_version, 
            covariates = covariates)

posthoc_report(exposure = exposure, 
               hrc_version = hrc_version,
               covariates = covariates,
               path = path)






# ================================================================== #
# ======= additional stuff ---- 
# ================================================================== #

# additional stuff
# 
# 2DF - chr13_47191972_G_A
table(figi$diab)

figi_nodiab <- figi %>% 
  dplyr::filter(diab == "No")

figi_diab <- figi %>% 
  dplyr::filter(diab == "Yes")

model1 <- glm(outcome ~ chr13_47191972_G_A + age_ref_imp + sex + pc1 + pc2 + pc3 , data = figi_nodiab, family = 'binomial')
model2 <- glm(outcome ~ chr13_47191972_G_A + age_ref_imp + sex + pc1 + pc2 + pc3 , data = figi_diab,   family = 'binomial')

summary(model1)
summary(model2)

mht_hits <- fread("~/mht_gxe_previously_reported.txt", stringsAsFactors = F) %>% 
  mutate(SNP2 = paste0(V1, ":", V2))

mht_hits_stats <- gxe %>% 
  filter(SNP2 %in% mht_hits$SNP2)

write.csv(mht_hits_stats, file = "~/Dropbox/mht_replication_stats.csv", quote = F, row.names = F)



# 3DF - chr8_118185025_G_A
table(figi$diab)

figi_nodiab <- figi %>% 
  dplyr::filter(diab == "No")

figi_diab <- figi %>% 
  dplyr::filter(diab == "Yes")

model1 <- glm(outcome ~ chr8_118185025_G_A + age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe, data = figi_nodiab, family = 'binomial')
model2 <- glm(outcome ~ chr8_118185025_G_A + age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe, data = figi_diab,   family = 'binomial')

summary(model1)
summary(model2)





mht_hits <- fread("~/mht_gxe_previously_reported.txt", stringsAsFactors = F) %>% 
  mutate(SNP2 = paste0(V1, ":", V2))

mht_hits_stats <- gxe %>% 
  filter(SNP2 %in% mht_hits$SNP2)

write.csv(mht_hits_stats, file = "~/Dropbox/mht_replication_stats.csv", quote = F, row.names = F)







