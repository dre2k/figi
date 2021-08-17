#=============================================================================#
# FIGI GxE smk_ever results
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
# source("functions.R")

# input variables
exposure = 'smk_ever'
hrc_version = 'v2.3'
annotation_file <- 'gwas_200_ld_annotation_feb2021.txt'
covariates <- sort(c('age_ref_imp', 'sex', 'study_gxe', 'pc1', 'pc2', 'pc3'))
path = glue("/media/work/gwis_test/{exposure}/")


# input data
esubset <- readRDS(glue("/media/work/gwis_test/{exposure}/data/FIGI_{hrc_version}_gxeset_{exposure}_basic_covars_glm.rds")) %>% 
  pull(vcfid)

input_data <- readRDS(glue("/media/work/gwis_test/data/FIGI_{hrc_version}_gxeset_analysis_data_glm.rds")) %>% 
  filter(vcfid %in% esubset)


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


# ------- main effects additional analyses ------- #
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
# output GxE models adjusted by different covariate sets (if requested)
# covariates_sets <- list(covariates)
# walk(snps, ~ fit_gxe_covars(data_epi = input_data, exposure = exposure, snp = .x, covariates_list = covariates_sets, method = 'chiSq3df', path = glue("{path}/output")))

# output RERI plots (can't install package on CARC yet)
significant_snps <- c("2:136407479:A:G", "1:151700227:G:A", "3:85461302:G:A", "4:147966066:G:A", "15:83613505:A:T") 
significant_snps <- c("3:85461302:G:A", "4:147966066:G:A") 
walk(significant_snps, ~ reri_wrapper(data_epi = input_data, exposure = exposure, snp = .x, covariates = covariates, path = glue("{path}/output")))


# get MAF plots for suggestive hits 
suggestive_hits <- fread(glue("/media/work/gwis_test/{exposure}/data/FIGI_{hrc_version}_gxeset_{exposure}_chiSqGxE_ldclump.clumped")) %>% 
  pull(SNP)

hits <- c("1:151700227:G:A", "15:83613505:A:T")
walk(hits, ~ output_aaf_plot(dat = input_data, exposure, snp = .x))









# ================================================================== # 
# investigate MECC exclusion/inclusion - some results look borderline when MECC included.. 
# answer - another case of different MAFs skewing findings
# ================================================================== #

# INCLUDE mecc
esubset <- readRDS(glue("/media/work/gwis_test/smk_ever_old/data/FIGI_{hrc_version}_gxeset_{exposure}_basic_covars_glm.rds")) %>% 
  pull(vcfid)
dat_all <- readRDS(glue("/media/work/gwis_test/data/FIGI_{hrc_version}_gxeset_analysis_data_glm.rds")) %>% 
  filter(vcfid%in% esubset)


# exclude MECC
esubset <- readRDS("/media/work/gwis_test/smk_ever/data/FIGI_v2.3_gxeset_smk_ever_basic_covars_glm.rds") %>% 
  pull(vcfid)
dat <- readRDS(glue("/media/work/gwis_test/data/FIGI_{hrc_version}_gxeset_analysis_data_glm.rds")) %>% 
  filter(vcfid%in% esubset)




# 2:136407479:A:G -- nominally high GxE p-value that disappears when excluding mecc

dose <- qread("/media/work/gwis_test/smk_ever/output/posthoc/dosage_chr2_136407479.qs")
out_all <- inner_join(dat_all, dose, 'vcfid')

model1 <- glm(outcome ~ smk_ever * chr2_136407479_A_G_dose + age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe, data = out_all, family = 'binomial')
out <- inner_join(dat, dose, 'vcfid')


model2 <- glm(outcome ~ smk_ever * chr2_136407479_A_G_dose + age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe, data = out, family = 'binomial')


summary(model1)
summary(model2) # p-value much higher (less sig)


# is it a MAF influence? 

# first, let's look at estimates for mECC only 

mecc <- filter(out_all, grepl("MECC", study_gxe))

model3 <- glm(outcome ~ smk_ever * chr2_136407479_A_G_dose + age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe, data = mecc, family = 'binomial')
summary(model3)



# compare MAF between studies, first look at all studies and then MECC vs the rest

aaf <- function(x) {
  sum(x) / nrow(x)
}

maf_compare <- out_all %>% 
  group_by(study_gxe) %>% 
  summarise(total = n(), 
            study_aaf = sum(chr2_136407479_A_G_dose) / (total*2)) %>% 
  arrange(study_aaf) %>% 
  mutate(study_gxe = fct_reorder(study_gxe, study_aaf))

ggplot(aes(x = study_gxe, y = study_aaf), data = maf_compare) + 
  geom_point() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 270)) + 
  xlab("Study") + 
  ylab("Alternate Allele Frequency")






# how about 5:134526551:C:T 
dose <- qread("/media/work/gwis_test/smk_ever/output/posthoc/dosage_chr5_134526551.qs")
out_all <- inner_join(dat_all, dose, 'vcfid')


aaf <- function(x) {
  sum(x) / nrow(x)
}

maf_compare <- out_all %>% 
  group_by(study_gxe) %>% 
  summarise(total = n(), 
            study_aaf = sum(chr5_134526551_C_T_dose) / (total*2)) %>% 
  arrange(study_aaf) %>% 
  mutate(study_gxe = fct_reorder(study_gxe, study_aaf))

ggplot(aes(x = study_gxe, y = study_aaf), data = maf_compare) + 
  geom_point() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 270)) + 
  xlab("Study") + 
  ylab("Alternate Allele Frequency")

# doesn't look amiss, how about model estimates. 



model1 <- glm(outcome ~ smk_ever * chr5_134526551_C_T_dose + age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe, data = out_all, family = 'binomial')
out <- inner_join(dat, dose, 'vcfid')
model2 <- glm(outcome ~ smk_ever * chr5_134526551_C_T_dose + age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe, data = out, family = 'binomial')
summary(model1)
summary(model2) # p-value much higher (less sig)





# how about 14:59209156:C:T

dose <- qread("/media/work/gwis_test/smk_ever/output/posthoc/dosage_chr14_59209156.qs")
out_all <- inner_join(dat_all, dose, 'vcfid')

aaf <- function(x) {
  sum(x) / nrow(x)
}

maf_compare <- out_all %>% 
  group_by(study_gxe) %>% 
  summarise(total = n(), 
            study_aaf = sum(chr14_59209156_C_T_dose) / (total*2)) %>% 
  arrange(study_aaf) %>% 
  mutate(study_gxe = fct_reorder(study_gxe, study_aaf))

ggplot(aes(x = study_gxe, y = study_aaf), data = maf_compare) + 
  geom_point() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 270)) + 
  xlab("Study") + 
  ylab("Alternate Allele Frequency")


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


posthoc_report(exposure = exposure)



