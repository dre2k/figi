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

# input variables
exposure = 'redmeatqcm_v2'
hrc_version = 'v3.0'
annotation_file <- 'gwas_200_ld_annotation_feb2021.txt'
covariates <- sort(c('age_ref_imp', 'sex', 'energytot_imp', 'study_gxe', 'pc1', 'pc2', 'pc3'))
path = glue("/media/work/gwis_test/{exposure}/")

# pca
pca <- fread("/media/work/gwis_test/PCA/20210222/figi_gxe_pca_update.eigenvec") %>% 
  dplyr::rename(vcfid = IID) %>% 
  dplyr::select(-`#FID`) %>% 
  dplyr::rename_with(tolower)

# input data
esubset <- readRDS(glue("/media/work/gwis_test/{exposure}/data/FIGI_{hrc_version}_gxeset_{exposure}_basic_covars_glm.rds")) %>% pull(vcfid)

input_data <- readRDS(glue("/media/work/gwis_test/data/FIGI_{hrc_version}_gxeset_analysis_data_glm.rds")) %>% 
  filter(vcfid%in% esubset) %>%
  dplyr::select(-starts_with("pc")) %>% 
  left_join(pca, 'vcfid') %>% 
  mutate(redmeatqcm_v2 = as.numeric(redmeatqcm))

#-----------------------------------------------------------------------------#
# Main effects ------------------------------------------------------------
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

# output GxE models adjusted by different covariate sets (if requested)
covariates_sets <- list(covariates)
snps <- c("8:122247679:C:G")
walk(snps, ~ fit_gxe_covars(data_epi = input_data, exposure = exposure, snp = .x, covariates_list = covariates_sets, method = 'chiSq2df', path = glue("{path}/output")))


# stratified odds ratios
snps <- c("8:122247679:C:G", "18:46453754:C:T")
snps <- c("18:46453754:C:T")
walk(snps , ~ fit_stratified_or_continuous(data_epi = input_data, exposure = exposure, snp = .x, hrc_version = hrc_version, covariates = covariates, dosage = F, path = glue("{path}/output"), flip_allele = T))


# interaction plots (let's rescale the cutoffs, or maybe make it binary?)
iplot_wrapper2(data_epi = input_data, exposure = exposure, snp = "8:122247679:C:G",  covariates = covariates, hrc_version = hrc_version, path = glue("{path}/output"))




# output RERI table (can't install package on CARC yet)
significant_snps <- c("8:122247679:C:G", "18:46453754:C:T")
walk(significant_snps, ~ reri_wrapper(data_epi = input_data, exposure = exposure, snp = .x, covariates = covariates, path = glue("{path}/output")))


input_data <- input_data %>% 
  mutate(energytot_imp1000 = energytot_imp / 1000)


glm(outcome ~ redmeatqcm_v2 + age_ref_imp + sex + energytot_imp + study_gxe + pc1 + pc2 + pc3, data = input_data, family = 'binomial')

glm(outcome ~ redmeatqcm_v2 + age_ref_imp + sex + energytot_imp1000 + study_gxe + pc1 + pc2 + pc3, data = input_data, family = 'binomial')








#-----------------------------------------------------------------------------#
# SNP information ---- 
#-----------------------------------------------------------------------------#


# need allele frequency by study_gxe to confirm finding (and sensitivity analysis)
tmp <- qread("/media/work/gwis_test/redmeatqcm_v2/output/posthoc/dosage_chr18_46453754.qs") %>% 
  inner_join(input_data, 'vcfid')

aaf <- function(x) {
  sum(x) / nrow(x)
}

tmp %>% 
  summarise(total = n(), 
            aaf = sum(chr18_46453754_C_T_dose / (total*2)))


out <- tmp %>% 
  group_by(study_gxe) %>% 
  summarise(total = n(), 
            study_aaf = sum(chr18_46453754_C_T_dose) / (total*2)) %>% 
  arrange(study_aaf) %>% 
  mutate(study_gxe = fct_reorder(study_gxe, study_aaf))

ggplot(aes(x = study_gxe, y = study_aaf), data = out) + 
  geom_point() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 270)) + 
  xlab("Study") + 
  ylab("Alternate Allele Frequency")

ggsave(filename = "~/Dropbox/asp_ref_chr6_32560631_C_T_AAF.png", width = 7.5, height = 4.72)


# would results change if you remove mecc (note wald statistic)
# EXCLUDE !!!!!!!!!!!!!
model1 <- glm(glue("outcome ~ {exposure}*chr6_32560631_C_T_dose + {glue_collapse(covariates, sep = '+')}"), data = chr6_325, family = 'binomial')
summary(model1)

chr6_325b <- dplyr::filter(chr6_325, !grepl("UKB", study_gxe))
model2 <- glm(glue("outcome ~ {exposure}*chr6_32560631_C_T_dose + {glue_collapse(covariates, sep = '+')}"), data = chr6_325b, family = 'binomial')
summary(model2)

















# need allele frequency by study_gxe to confirm finding (and sensitivity analysis)
tmp <- qread("/media/work/gwis_test/redmeatqcm_v2/output/posthoc/dosage_chr8_122247679.qs") %>% 
  inner_join(input_data, 'vcfid')

aaf <- function(x) {
  sum(x) / nrow(x)
}


out <- tmp %>% 
  group_by(study_gxe) %>% 
  summarise(total = n(), 
            study_aaf = sum(chr8_122247679_C_G_dose) / (total*2)) %>% 
  arrange(study_aaf) %>% 
  mutate(study_gxe = fct_reorder(study_gxe, study_aaf))

ggplot(aes(x = study_gxe, y = study_aaf), data = out) + 
  geom_point() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 270)) + 
  xlab("Study") + 
  ylab("Alternate Allele Frequency")






tmp <- qread("/media/work/gwis_test/redmeatqcm_v2/output/posthoc/dosage_chr18_46453754.qs") %>% 
  inner_join(input_data, 'vcfid')

out <- tmp %>% 
  group_by(study_gxe) %>% 
  summarise(total = n(), 
            study_aaf = sum(chr18_46453754_C_T_dose) / (total*2)) %>% 
  arrange(study_aaf) %>% 
  mutate(study_gxe = fct_reorder(study_gxe, study_aaf))

ggplot(aes(x = study_gxe, y = study_aaf), data = out) + 
  geom_point() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 270)) + 
  xlab("Study") + 
  ylab("Alternate Allele Frequency")








#-----------------------------------------------------------------------------#
# simple plots ---- 
#-----------------------------------------------------------------------------#


label <- c("TT", "CT", "CC")
mean  <- c(1.46, 1.35, 1.18)
lower <- c(1.26, 1.26, 1.11)
upper <- c(1.69, 1.44, 1.24)
df <- data.frame(label, mean, lower, upper)

ggplot(df, aes(y = mean, x = label)) +
  geom_errorbar(aes(ymax = upper, ymin = lower), size = .5, width = .1, color = "gray50") + 
  geom_point(size = 3.5, color = "orange") + 
  geom_hline(yintercept = 1.0, linetype = "dotted", size = 1) +
  scale_y_log10(breaks =  seq(0.5, 2.0, by = 0.1),  minor_breaks = NULL) +
  coord_flip(ylim = c(0.7, 2)) +
  theme_bw() + labs(y = "Odds ratio (log)", x = "") + 
  ggtitle("redmeatqcm odds ratios by chr18:46453754 genotype")


ggsave(filename = "~/Dropbox/redmeatqcm_by_chr18_46453754_C_T_genotype.png", width = 6, height = 4)




label <- c("CT vs TT, redmeatqcm = 0", 
           "CC vs TT, redmeatqcm = 0", 
           "CT vs TT, redmeatqcm = 1",
           "CC vs TT, redmeatqcm = 1")
mean  <- c(1.29, 1.54, 1.19, 1.25)
lower <- c(1.14, 1.38, 1.08, 1.13)
upper <- c(1.45, 1.73, 1.31, 1.37)
df <- data.frame(label, mean, lower, upper)


ggplot(df, aes(y = mean, x = label)) +
  geom_errorbar(aes(ymax = upper, ymin = lower), size = .5, width = .1, color = "gray50") + 
  geom_point(size = 3.5, color = "red") + 
  geom_hline(yintercept = 1.0, linetype = "dotted", size = 1) +
  # scale_y_log10(breaks = c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0),  minor_breaks = NULL) +
  scale_y_log10(breaks = seq(0.5, 2.0, by = 0.1),  minor_breaks = NULL) +
  labs(y = "Odds ratio", x = "Effect") +
  coord_flip(ylim = c(0.7, 2)) +
  theme_bw() + labs(y = "Odds ratio (log)", x = "") + 
  ggtitle("chr18:46453754 ORs by redmeatqcm intake (1 vs 0 servings/day)")

ggsave(filename = "~/Dropbox/chr18_46453754_C_T_by_redmeatqcm.png", width = 6, height = 4)




# ================================================================== #
# ======= rmarkdown reports ---- 
# ================================================================== #
main_effects_report(exposure = exposure, hrc_version = hrc_version, covariates = covariates, path = path)

gwis_report(exposure = exposure, hrc_version = hrc_version, covariates = covariates)


posthoc_report(exposure = exposure)

posthoc_report(exposure = exposure, 
               hrc_version = hrc_version,
               covariates = covariates,
               path = path)



out <- data.frame(table(input_data$study_gxe, input_data$redmeatqcm_v2, input_data$sex)) %>% 
  dplyr::filter(Freq != 0)


