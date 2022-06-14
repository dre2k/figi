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

# input variables
exposure = 'calcium_dietqc2'
hrc_version = 'v3.0'
annotation_file <- 'gwas_200_ld_annotation_feb2021.txt'
covariates <- sort(c('age_ref_imp', 'sex', 'energytot_imp', 'study_gxe', 'pc1', 'pc2', 'pc3'))
path = glue("/media/work/gwis_test/{exposure}/")

esubset <- readRDS(glue("/media/work/gwis_test/{exposure}/data/FIGI_{hrc_version}_gxeset_{exposure}_basic_covars_glm.rds")) %>% pull(vcfid)
input_data <- readRDS(glue("/media/work/gwis_test/data/FIGI_{hrc_version}_gxeset_analysis_data_glm.rds")) %>% 
  filter(vcfid%in% esubset)

#-----------------------------------------------------------------------------#
# main effects ----
#-----------------------------------------------------------------------------#


# ------ meta-analysis ------ #
output_dir = as.character(glue("{path}/output/posthoc/"))
covariates_meta <- sort(covariates[which(!covariates %in% c(paste0(rep('pc', 20), seq(1, 20)), "study_gxe"))])
# covariates_meta <- sort(covariates[which(!covariates %in% c("study_gxe"))])


create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "all", forest_height = 15, forest_width = 9)
create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "proximal", forest_height = 13, forest_width = 9)
create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "distal", forest_height = 13, forest_width = 9)
create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "rectal", forest_height = 13, forest_width = 9)
create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "female", forest_height = 13, forest_width = 9)
create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "male", forest_height = 13, forest_width = 9)




# ---- what happened below (4/20/2022) ---- #
# for meta - need to remove PLCO_4 when performing tumor stratified analysis
# or perhaps limit to sex-stratified unless there's strong evidence of site stratified result difference


# ------ meta-analysis ------ #
covariates_meta <- sort(covariates[which(!covariates %in% c(paste0(rep('pc', 20), seq(1, 20)), 'study_gxe'))])
meta_analysis_execute2(dataset = input_data, exposure = exposure, hrc_version = hrc_version, covariates = covariates_meta, output_dir = output_dir, filename_suffix2 = "")

model_form <- Reduce(paste, deparse(reformulate(c(exposure, covariates_meta), response = 'outcome')))

xx <- meta_analysis_wrapper_pwrap_tableonly(dataset = input_data, exposure = exposure, model_formula = model_form, subset_var = 'all')



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
twostep_eh_snps(gxe, 'chiSqEDGE')
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
# SNP followup (models and plots) ---- 
#-----------------------------------------------------------------------------#

# have to recode calcium for posthoc analyses (protective)
figi <- posthoc_input(exposure, hrc_version, glue('gwis_sig_results_output_{exposure}.rds')) %>% 
  mutate(calcium_dietqc2 = abs(calcium_dietqc2 - 3))

snps <- readRDS(glue(output_dir, "gwis_sig_results_input_{exposure}.rds"))
table(snps$method)


# ------- GxE 1DF suggestive findings ------- #
snps_filter <- snps %>%
  dplyr::filter(method == glue("manhattan_chiSqGxE_{exposure}_clump_df")) %>%
  dplyr::mutate(snps = paste0('chr', gsub("\\:", "\\_", SNP))) %>%
  pull(snps)

posthoc_run_models(snps_filter, 'chiSqGxE')
posthoc_run_reri_plot(snps_filter)

posthoc_create_plots(snps_filter[3], 'chiSqGxE')


# ------- 2DF/3DF clumped removing GWAS and or EG (to check 'novel' D|G results) ------- #
# as first pass, only generate locuszoom plots to gauge whether the region is already capture by GECCO meta-analysis
snps_filter <- readRDS(glue(output_dir, "manhattan_chiSq2df_{exposure}_no_gwas_clump_df.rds")) %>% 
  dplyr::mutate(snps = paste0('chr', gsub("\\:", "\\_", SNP))) %>% 
  pull(snps)
posthoc_create_plots_locuszoom(snps_filter, 'chiSq2df')

snps_filter <- readRDS(glue(output_dir, "manhattan_chiSq3df_{exposure}_no_gwas_no_ge_clump_df.rds")) %>% 
  dplyr::mutate(snps = paste0('chr', gsub("\\:", "\\_", SNP))) %>% 
  pull(snps)
posthoc_create_plots_locuszoom(snps_filter, 'chiSq3df')



# ------- 2DF/3DF clumped removing GWAS and or EG (to check 'novel' D|G results) ------- #
# - locuszoom plots
# - conditional analyses??? (more involved.. but maybe worthwhile if they request it)

snps_2df <- readRDS(glue(output_dir, "manhattan_chiSq2df_{exposure}_no_gwas_clump_df.rds")) %>% 
  dplyr::mutate(snps = paste0('chr', gsub("\\:", "\\_", SNP))) %>% 
  pull(snps)

posthoc_run_models(snps_2df, 'chiSq2df')
posthoc_run_reri_plot(snps_2df)
posthoc_create_plots(snps_2df[1], 'chiSqG')



snps_3df <- readRDS(glue(output_dir, "manhattan_chiSq3df_{exposure}_no_gwas_no_ge_clump_df.rds")) %>% 
  dplyr::mutate(snps = paste0('chr', gsub("\\:", "\\_", SNP))) %>% 
  pull(snps)

posthoc_run_models(snps_3df, 'chiSq3df')
posthoc_run_reri_plot(snps_3df)
posthoc_create_plots(snps_3df, 'chiSq3df')


# ------- Expectation based hybrid... (just do all...) -------- #
snps_filter <- snps %>% 
  dplyr::filter(method == glue("twostep_wht_chiSqG_{exposure}_expectation_hybrid_df")) %>% 
  dplyr::mutate(snps = paste0('chr', gsub("\\:", "\\_", SNP))) %>% 
  pull(snps)

posthoc_run_models(snps_filter, 'chiSqGxE')
posthoc_run_reri_plot(snps_filter)
posthoc_create_plots(snps_filter, 'chiSqGxE')




#-----------------------------------------------------------------------------#
# SNP followup (models and plots) ---- 
# updated version
#-----------------------------------------------------------------------------#

# stratified odds ratios for two-step hybrid top hit
x <- readRDS("/media/work/gwis_test/calcium_totqc2/output/posthoc/stratified_oddsratio_calcium_totqc2_v3.0_chr2_135594399_A_C_age_ref_imp_energytot_imp_pc1_pc2_pc3_sex_study_gxe.rds")

kable(x) %>% 
  kable_styling()




# get MAF plots for suggestive hits 
suggestive_hits <- fread(glue("/media/work/gwis_test/{exposure}/data/FIGI_{hrc_version}_gxeset_{exposure}_chiSqGxE_ldclump.clumped")) %>% 
  pull(SNP)

walk(suggestive_hits, ~ output_aaf_plot(dat = input_data, exposure, snp = .x))






# ================================================================== #
# ======= rmarkdown reports ---- 
# ================================================================== #
gwis_report(exposure = exposure, 
            hrc_version = hrc_version, 
            covariates = covariates)

posthoc_report(exposure = 'calcium_dietqc2')






