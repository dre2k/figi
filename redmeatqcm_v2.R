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
# GxE additional analysis -------------------------------------------------
#-----------------------------------------------------------------------------#
GxE <- c("17:71354811:G:A", "9:120653523:T:G", "4:31160820:G:A", "9:120290042:A:G", "8:122247679:C:G", "1:21495668:G:T", "3:29246192:C:T", "1:21209766:T:A", "15:100002029:C:T", "8:112374447:C:T", "18:21029284:C:T", "20:51479759:G:T", "2:42324419:T:A", "10:29013189:C:T", "13:62734417:A:G", "1:41466550:T:A", "3:27633795:G:A", "2:73204203:G:A", "1:239473113:T:C")
df2 <- c("8:120042899:G:T", "8:122247679:C:G", "10:86469912:G:A")
df3 <- c("8:122247679:C:G")
twostep <- c("18:46453754:C:T")
eg <- c("19:49232226:G:A")


snps <- unique(c(GxE, df2, df3, twostep, eg))


# MAF ---------------------------------------------------------------------
walk(snps, ~ create_aaf_study_plot(data = input_data, exposure = exposure, hrc_version = hrc_version, snp = .x, path = path))



# stratified odds ratio ---------------------------------------------------
snps <- unique(c(df2, df3, twostep, eg))
snps <- c("8:120042899:G:T", "8:122247679:C:G", "10:86469912:G:A")
walk(snps , ~ fit_stratified_or_continuous(data_epi = input_data, exposure = exposure, snp = .x, hrc_version = hrc_version, covariates = covariates, dosage = F, path = glue("{path}/output"), flip_allele = F))

# some alt alleles are protective, consider flipping
# snps <- c('18:46453754:C:T')
# walk(snps , ~ fit_stratified_or_continuous(data_epi = input_data, exposure = exposure, snp = .x, hrc_version = hrc_version, covariates = covariates, dosage = F, path = glue("{path}/output"), flip_allele = T))



# interaction plots -------------------------------------------------------
snps <- c('18:46453754:C:T')
snps <- c("8:120042899:G:T", "8:122247679:C:G", "10:86469912:G:A")
walk(snps, ~ iplot_wrapper(data_epi = input_data, exposure = exposure, snp = .x,  covariates = covariates, hrc_version = hrc_version, path = glue("{path}/output"), flip_allele = F))



# Additional covariates ---------------------------------------------------
# output GxE models adjusted by different covariate sets (if requested)
# height, smoking, , alcohol intake, total fruit and vegetable intake, 


# physical activity - too many missing.. 

cov1 <- covariates
cov2 <- c(cov1, 'bmi5', 'height10')
cov3 <- c(cov2, 'smk_ever', 'alcoholc', 'fruitqc2', 'fiberqc2', 'vegetableqc2')


covariates_sets <- list(cov1, cov2, cov3)
# snps <- c("8:122247679:C:G")
snps <- unique(c(df2, twostep))
walk(snps, ~ fit_gxe_covars(data_epi = input_data, exposure = exposure, snp = .x, covariates_list = covariates_sets, method = 'chiSqGxE', path = glue("{path}/output")))



# Interaction plots -------------------------------------------------------
# interaction plots (let's rescale the cutoffs, or maybe make it binary?)
iplot_wrapper2(data_epi = input_data, exposure = exposure, snp = "8:122247679:C:G",  covariates = covariates, hrc_version = hrc_version, path = glue("{path}/output"))




# RERI --------------------------------------------------------------------
# output RERI table (can't install package on CARC yet)
significant_snps <- c("8:122247679:C:G", "18:46453754:C:T")
walk(significant_snps, ~ reri_wrapper(data_epi = input_data, exposure = exposure, snp = .x, covariates = covariates, path = glue("{path}/output")))




#-----------------------------------------------------------------------------#
# simple plots ---- 
#-----------------------------------------------------------------------------#


# ------ simple forest plot of the gxE interaction parameters ----- #
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

