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
rm(list = ls())
source("functions.R")

# input variables
exposure = 'asp_ref'
hrc_version = 'v2.4'
output_dir = paste0("/media/work/gwis/posthoc/", exposure, "/")
annotation_file <- 'gwas_141_ld_annotation_july2020.txt'

covariates <- sort(c('age_ref_imp', 'sex', 'study_gxe', 'pc1', 'pc2', 'pc3'))
covariates_set1 <- covariates
# covariates_set2 <- sort(c(covariates_set1, 'bmi', 'smk_ever', 'fruitqc2', 'vegetableqc2'))
covariates_list <- list(covariates)
mod <- 'age_ref_imp+sex+pc1+pc2+pc3+study_gxe'

exclude_gwas <- fread("~/data/Annotations/gwas_141_ld_annotation_july2020.txt") %>%
  mutate(SNP2 = paste0(Chr, ":", Pos)) %>%
  pull(SNP2)

# input data
esubset <- readRDS(glue("/media/work/gwis/data/FIGI_EpiData/FIGI_{hrc_version}_gxeset_{exposure}_basic_covars_glm.rds")) %>% 
  pull(vcfid)
input_data <- readRDS(glue("/media/work/gwis/data/FIGI_EpiData/FIGI_{hrc_version}_gxeset_analysis_data_glm.rds")) %>% 
  filter(vcfid%in% esubset)

# gxescan output
gxe <- readRDS(glue("/media/work/gwis/results/{exposure}/processed/FIGI_{hrc_version}_gxeset_{exposure}_basic_covars_gxescan_results.rds")) %>% 
  mutate(SNP2 = paste0(Chromosome, ":", Location))



#-----------------------------------------------------------------------------#
# main effects ----
#-----------------------------------------------------------------------------#

# ------ meta-analysis ------ #
covariates_meta <- sort(c('age_ref_imp', 'sex', 'study_gxe'))
meta_analysis_execute(dataset = input_data, exposure = exposure, hrc_version = hrc_version, covariates = covariates_meta, output_dir = output_dir, filename_suffix2 = "")

# ------- pooled analysis ------- #
## basic covariates
covariates_pooled <- sort(covariates[which(!covariates %in% c(paste0(rep('pc', 20), seq(1, 20))))])
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
# qq plots, manhattan plots, two-step plots ---- 
# (should be able to source script since variables are defined here)
#-----------------------------------------------------------------------------#
source("~/Dropbox/FIGI/FIGI_code/results/gwis/gwis_02_plots.R")


#-----------------------------------------------------------------------------#
# two-step expectation hybrid ---- 
#-----------------------------------------------------------------------------#
# ------ output list of SNPs from each bin (expectation based) to extract dosages
# exposure subset for step 1 statistics
twostep_eh_snps(gxe, 'chiSqG')
twostep_eh_snps(gxe, 'chiSqGE')
twostep_eh_snps(gxe, 'chiSqEDGE')

# gecco meta analysis for main effects step 1 statistics
gwas_results <- readRDS("/media/work/gwis/results/gecco/MarginalMeta_gecco_HRC_EUR_only.rds")
gxe_twostep_gwas_step1 <- gxe %>%
  dplyr::select(-betaG, -chiSqG, -chiSqEDGE, -chiSq3df) %>%
  inner_join(gwas_results, 'SNP') %>%
  mutate(chiSqEDGE = chiSqG + chiSqGE,
         chiSq3df = chiSqG + chiSqGxE + chiSqGE)

twostep_eh_snps(gxe_twostep_gwas_step1, 'chiSqG', step1_source = "gecco")


# ----------------------------------------------------------- #
# ------- run simpleM and generate plots + results df ------- #
simplem_wrap(x = gxe, exposure = exposure, covariates = covariates, simplem_step1_statistic = 'chiSqG', output_dir = output_dir)
simplem_wrap(x = gxe, exposure = exposure, covariates = covariates, simplem_step1_statistic = 'chiSqGE', output_dir = output_dir)
simplem_wrap(x = gxe, exposure = exposure, covariates = covariates, simplem_step1_statistic = 'chiSqEDGE', output_dir = output_dir)
simplem_wrap(x = gxe_twostep_gwas_step1, exposure = exposure, covariates = covariates, simplem_step1_statistic = 'chiSqG', output_dir = output_dir, filename_suffix = "_gwas_step1" )


simplem_wrap(x = gxe, exposure = exposure, covariates = covariates, simplem_step1_statistic = 'chiSqG', output_dir = output_dir, include_gwas = F, filename_suffix = "_no_gwas")
simplem_wrap(x = gxe, exposure = exposure, covariates = covariates, simplem_step1_statistic = 'chiSqGE', output_dir = output_dir, include_gwas = F, filename_suffix = "_no_gwas")
simplem_wrap(x = gxe, exposure = exposure, covariates = covariates, simplem_step1_statistic = 'chiSqEDGE', output_dir = output_dir, include_gwas = F, filename_suffix = "_no_gwas")
simplem_wrap(x = gxe_twostep_gwas_step1, exposure = exposure, covariates = covariates, simplem_step1_statistic = 'chiSqG', output_dir = output_dir, include_gwas = F, filename_suffix = "_gwas_step1_no_gwas" )




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
# SNP followup ---- 
#-----------------------------------------------------------------------------#
# compile significant results into data.frame
source("/home/rak/Dropbox/FIGI/FIGI_code/results/posthoc/posthoc_01_combine_results.R")

# on HPC:
# - extract dosage information on HPC
# - calculate clumped statistics - gxe, 2/3df if necessary



#-----------------------------------------------------------------------------#
# SNP followup (models and plots) ---- 
#-----------------------------------------------------------------------------#
figi <- posthoc_input(exposure, hrc_version, glue('gwis_sig_results_output_{exposure}.rds')) %>% 
  mutate(asp_ref = fct_relevel(asp_ref, "Yes", "No"))

snps <- readRDS(glue(output_dir, "gwis_sig_results_input_{exposure}.rds"))
table(snps$method)




# ------- GxE 1DF findings ------- #
snps_filter <- snps %>% 
  dplyr::filter(method == glue("manhattan_chiSqGxE_{exposure}_clump_df"),
                Pval <= 5e-8) %>% 
  dplyr::mutate(snps = paste0('chr', gsub("\\:", "\\_", SNP))) %>% 
  pull(snps)

posthoc_run_models(snps_filter, 'chiSqGxE', quartile = T)
posthoc_run_reri_plot(snps_filter)
posthoc_create_plots(snps_filter, 'chiSqGxE')


# ------- two step (E|G no gwas) ------- #
snps_filter <- snps %>% 
  dplyr::filter(method == glue("twostep_wht_chiSqGE_{exposure}_no_gwas_df")) %>% 
  dplyr::mutate(snps = paste0('chr', gsub("\\:", "\\_", SNP))) %>% 
  pull(snps)

posthoc_run_models(snps_filter, 'chiSqGxE', quartile = T)
posthoc_run_reri_plot(snps_filter)
posthoc_create_plots(snps_filter, 'chiSqGxE')





# ------- twostep eh ------- #
snps_filter <- snps %>% 
  dplyr::filter(method == glue("twostep_wht_chiSqEDGE_{exposure}_expectation_hybrid_no_gwas_df")) %>% 
  dplyr::arrange(Pval) %>% 
  dplyr::mutate(snps = paste0('chr', gsub("\\:", "\\_", SNP))) %>% 
  pull(snps)

posthoc_run_models(snps_filter, 'chiSqGxE', quartile = T)
posthoc_run_reri_plot(snps_filter)
posthoc_create_plots(snps_filter, 'chiSqGxE')




# ------- 2DF/3DF clumped removing GWAS and or EG (to check 'novel' D|G results) ------- #
# as first pass, only generate locuszoom plots to gauge whether the region is already capture by GECCO meta-analysis
snps_filter <- readRDS(glue(output_dir, "manhattan_chiSq2df_{exposure}_no_gwas_clump_df.rds")) %>% 
  dplyr::mutate(snps = paste0('chr', gsub("\\:", "\\_", SNP))) %>% 
  pull(snps)
posthoc_create_plots(snps_filter, 'chiSq2df')

snps_filter <- readRDS(glue(output_dir, "manhattan_chiSq3df_{exposure}_no_gwas_no_ge_clump_df.rds")) %>% 
  dplyr::mutate(snps = paste0('chr', gsub("\\:", "\\_", SNP))) %>% 
  pull(snps)
posthoc_create_plots(snps_filter, 'chiSq3df')










# ================================================================== #
# ======= rmarkdown reports ---- 
# ================================================================== #
gwis_report(exposure = exposure, 
            hrc_version = hrc_version, 
            covariates = covariates)

posthoc_report(exposure = exposure)





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
