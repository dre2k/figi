#=============================================================================#
# FIGI GxE hrt_pm_ref2 results
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
exposure = 'hrt_ref_pm2'
hrc_version = 'v2.3'
output_dir = paste0("/media/work/gwis/posthoc/", exposure, "/")
annotation_file <- 'gwas_141_ld_annotation_july2020.txt'

covariates <- sort(c('age_ref_imp', 'study_gxe', 'pc1', 'pc2', 'pc3'))
covariates_set1 <- covariates
# covariates_set2 <- sort(c(covariates_set1, 'bmi', 'smk_ever', 'fruitqc2', 'vegetableqc2'))
covariates_list <- list(covariates)
mod <- 'age_ref_imp+pc1+pc2+pc3+study_gxe'

# input data
exposure_subset <- readRDS(glue("/media/work/gwis/data/FIGI_EpiData/FIGI_{hrc_version}_gxeset_{exposure}_basic_covars_glm.rds"))[,'vcfid']

input_data <- readRDS(glue("/media/work/gwis/data/FIGI_EpiData/FIGI_{hrc_version}_gxeset_analysis_data_glm.rds")) %>% 
  filter(vcfid %in% exposure_subset)



# gxescan output
# gxe <- readRDS(glue("/media/work/gwis/results/{exposure}/processed/FIGI_{hrc_version}_gxeset_{exposure}_basic_covars_gxescan_results.rds")) %>% 
#   mutate(SNP2 = paste0(Chromosome, ":", Location))

path = glue("/media/work/gwis_test/{exposure}/")



#-----------------------------------------------------------------------------#
# main effects ----
#-----------------------------------------------------------------------------#

# ------ meta-analysis ------ #
meta_analysis_execute(dataset = input_data, exposure = exposure, hrc_version = hrc_version, covariates = covariates, output_dir = output_dir, filename_suffix2 = "", nosex = T)


# ------- pooled analysis ------- #
## basic covariates
covariates_pooled <- sort(covariates[which(!covariates %in% c(paste0(rep('pc', 20), seq(1, 20))))])
covariates_pooled <- covariates

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

# ------- MAF Plots  ---------- #
snps_out <- c("12:13670508:G:C", "6:117823508:T:C")
walk(snps_out, ~ create_aaf_study_plot(data = input_data, exposure,  hrc_version, snp = .x, path = path))



figi <- posthoc_input(exposure, hrc_version, glue('gwis_sig_results_output_{exposure}.rds')) %>% 
  mutate(hrt_ref_pm2 = fct_relevel(hrt_ref_pm2, "Yes", "No"))

figi <- posthoc_input(exposure, hrc_version, glue('gwis_sig_results_output_{exposure}.rds')) 

snps <- readRDS(glue(output_dir, "gwis_sig_results_input_{exposure}.rds"))
table(snps$method)



# ------- GxE 1DF findings ------- #
snps_filter <- snps %>% 
  dplyr::filter(method == glue("manhattan_chiSqGxE_{exposure}_clump_df"),
                Pval <= 5e-6) %>% 
  dplyr::mutate(snps = paste0('chr', gsub("\\:", "\\_", SNP))) %>% 
  pull(snps)

posthoc_run_models(snps_filter, 'chiSqGxE')
posthoc_run_reri_plot(snps_filter)
posthoc_create_plots(snps_filter, 'chiSqGxE')


# ------- two step (E|G no gwas) ------- #
snps_filter <- snps %>% 
  dplyr::filter(method == glue("twostep_wht_chiSqGE_{exposure}_no_gwas_df")) %>% 
  dplyr::mutate(snps = paste0('chr', gsub("\\:", "\\_", SNP))) %>% 
  pull(snps)

posthoc_run_models(snps_filter, 'chiSqGxE')
posthoc_run_reri_plot(snps_filter)
posthoc_create_plots(snps_filter, 'chiSqGxE')





# ------- twostep eh ------- #
snps_filter <- snps %>% 
  dplyr::filter(method == glue("twostep_wht_chiSqEDGE_{exposure}_expectation_hybrid_no_gwas_df")) %>% 
  dplyr::arrange(Pval) %>% 
  dplyr::mutate(snps = paste0('chr', gsub("\\:", "\\_", SNP))) %>% 
  pull(snps)

snps_filter <- snps_filter[1:6]

posthoc_run_models(snps_filter, 'chiSqGxE')
posthoc_run_reri_plot(snps_filter)
posthoc_create_plots(snps_filter, 'chiSqGxE')






# ------- 2DF/3DF clumped removing GWAS and or EG (to check 'novel' D|G results) ------- #
# as first pass, only generate locuszoom plots to gauge whether the region is already capture by GECCO meta-analysis
snps_filter <- readRDS(glue(output_dir, "manhattan_chiSq2df_{exposure}_no_gwas_clump_df.rds")) %>% 
  dplyr::mutate(snps = paste0('chr', gsub("\\:", "\\_", SNP))) %>% 
  pull(snps)


posthoc_run_models(snps_filter, 'chiSq2df')
posthoc_run_reri_plot(snps_filter)
posthoc_create_plots(snps_filter, 'chiSq2df')


snps_filter <- readRDS(glue(output_dir, "twostep_wht_chiSqGE_{exposure}_df.rds")) %>% 
  dplyr::mutate(snps = paste0('chr', gsub("\\:", "\\_", SNP))) %>% 
  pull(snps)

posthoc_run_models(snps_filter, 'chiSqGxE')
posthoc_run_reri_plot(snps_filter)
posthoc_create_plots(snps_filter, 'chiSqGxE')


posthoc_create_plots(snps_filter, 'chiSq2df')

snps_filter <- readRDS(glue(output_dir, "manhattan_chiSq3df_{exposure}_no_gwas_no_ge_clump_df.rds")) %>% 
  dplyr::mutate(snps = paste0('chr', gsub("\\:", "\\_", SNP))) %>% 
  pull(snps)
posthoc_create_plots(snps_filter, 'chiSq3df')










#-----------------------------------------------------------------------------#
# Posthoc requests ---- 
#-----------------------------------------------------------------------------#

# ------- some basic checking to make sure it's right ------- #
chr6_117823508_T_C_out <- gxe %>% 
filter(Location == 117823508)
chr6_117823508_T_C_out

figi <- posthoc_input(exposure, hrc_version, glue('gwis_sig_results_output_{exposure}.rds')) %>% 
  mutate(hrt_ref_pm2 = fct_relevel(hrt_ref_pm2, "Yes", "No"))
figi <- posthoc_input(exposure, hrc_version, glue('gwis_sig_results_output_{exposure}.rds')) %>% 
  mutate(hrt_ref_pm2 = as.numeric(hrt_ref_pm2))

snps <- readRDS(glue(output_dir, "gwis_sig_results_input_{exposure}.rds"))

# remember to use likelihood ratio test
model_dg <-   glm(outcome ~ chr6_117823508_T_C_dose + hrt_ref_pm2 + age_ref_imp + pc1 + pc2 + pc3 + study_gxe, data = figi, family = 'binomial')
lrtest(model_dg, "chr6_117823508_T_C_dose")


# 2DF test.. (OK)
model_base <- glm(outcome ~                           hrt_ref_pm2 + age_ref_imp + pc1 + pc2 + pc3 + study_gxe, data = figi, family = 'binomial')
model_dg <-   glm(outcome ~ chr6_117823508_T_C_dose * hrt_ref_pm2 + age_ref_imp + pc1 + pc2 + pc3 + study_gxe, data = figi, family = 'binomial')
lrtest(model_dg, model_base)



# ------- top hits limited to COLON cancer only ------- #


# Hrt_ref_pm2 x chr6_117823508_T_C_dose -- 2DF results
figi_colon <- posthoc_input(exposure, hrc_version, glue('gwis_sig_results_output_{exposure}.rds')) %>% 
  mutate(hrt_ref_pm2 = fct_relevel(hrt_ref_pm2, "Yes", "No")) %>% 
  filter(outcome == 0 | cancer_site_sum1 == "colon")
table(figi_colon$outcome)


model_base  <- glm(outcome ~                           hrt_ref_pm2 + age_ref_imp + pc1 + pc2 + pc3 + study_gxe, data = figi_colon, family = 'binomial')
model_out   <- glm(outcome ~ chr6_117823508_T_C_dose * hrt_ref_pm2 + age_ref_imp + pc1 + pc2 + pc3 + study_gxe, data = figi_colon, family = 'binomial')
lrtest(model_out, model_base)


model_out   <- glm(outcome ~ chr6_117823508_T_C_dose + hrt_ref_pm2 + age_ref_imp + pc1 + pc2 + pc3 + study_gxe, data = figi_colon, family = 'binomial')
lrtest(model_out, model_base)


model_base  <- glm(outcome ~ chr6_117823508_T_C_dose + hrt_ref_pm2 + age_ref_imp + pc1 + pc2 + pc3 + study_gxe, data = figi_colon, family = 'binomial')
model_out   <- glm(outcome ~ chr6_117823508_T_C_dose * hrt_ref_pm2 + age_ref_imp + pc1 + pc2 + pc3 + study_gxe, data = figi_colon, family = 'binomial')
lrtest(model_out, model_base)



# Hrt_ref_pm2 x chr12_13670508_G_C_dose -- twostep results (E|G step 1 filter)
chr12_13670508_G_C_out <- gxe %>% 
  filter(Location == 13670508)
chr12_13670508_G_C_out

figi_colon <- posthoc_input(exposure, hrc_version, glue('gwis_sig_results_output_{exposure}.rds')) %>% 
  mutate(hrt_ref_pm2 = fct_relevel(hrt_ref_pm2, "Yes", "No"), 
         hrt_ref_pm2 = as.numeric(hrt_ref_pm2)) %>% 
  filter(outcome == 0 | cancer_site_sum1 == "colon")

model_base  <- lm(hrt_ref_pm2 ~                           age_ref_imp + pc1 + pc2 + pc3 + study_gxe, data = figi_colon)
model_out   <- lm(hrt_ref_pm2 ~ chr12_13670508_G_C_dose + age_ref_imp + pc1 + pc2 + pc3 + study_gxe, data = figi_colon)
lrtest(model_out, model_base)

model_base  <- glm(outcome ~ chr12_13670508_G_C_dose + hrt_ref_pm2 + age_ref_imp + pc1 + pc2 + pc3 + study_gxe, data = figi_colon, family = 'binomial')
model_out   <- glm(outcome ~ chr12_13670508_G_C_dose * hrt_ref_pm2 + age_ref_imp + pc1 + pc2 + pc3 + study_gxe, data = figi_colon, family = 'binomial')
lrtest(model_out, model_base)






# ================================================================== #
# ======= check HRT literature SNPs


# WHI - main effects

rs17724534 <- qread("/media/work/gwis_test/hrt_ref_pm2/output/posthoc/dosage_chr10_104605521.qs") %>% 
  rename(snp = chr10_104605521_C_T_dose)
rs10883782 <- qread("/media/work/gwis_test/hrt_ref_pm2/output/posthoc/dosage_chr10_104583932.qs") %>% 
  rename(snp = chr10_104583932_A_G_dose)
rs1902584 <- qread("/media/work/gwis_test/hrt_ref_pm2/output/posthoc/dosage_chr15_51611654.qs") %>% 
  rename(snp = chr15_51611654_A_T_dose)

whi <- filter(input_data, grepl("WHI", study_gxe))

# no nominal significance when using all WHI studies .. 
rs17724534_dg <- glm(outcome ~ snp + age_ref_imp + pc1 + pc2 + pc3, 
                     data = inner_join(whi, rs17724534, 'vcfid'), family = 'binomial')
summary(rs17724534_dg)

rs10883782_dg <- glm(outcome ~ snp + age_ref_imp + pc1 + pc2 + pc3, 
                     data = inner_join(whi, rs10883782, 'vcfid'), family = 'binomial')
summary(rs10883782_dg)

rs1902584_dg <- glm(outcome ~ snp + age_ref_imp + pc1 + pc2 + pc3, 
                     data = inner_join(whi, rs1902584, 'vcfid'), family = 'binomial')
summary(rs1902584_dg)







# DACHS GxE


rs1202168 <- qread("/media/work/gwis_test/hrt_ref_pm2/output/posthoc/dosage_chr6_152432902.qs") %>% 
  rename(snp = chr6_152432902_C_T_dose)

dachs <- filter(input_data, grepl("DACHS", study_gxe))


rs1202168_dg <- glm(outcome ~ snp * hrt_ref_pm2 + age_ref_imp + pc1 + pc2 + pc3, 
                     data = inner_join(dachs, rs1202168, 'vcfid'), family = 'binomial')

summary(rs1202168_dg)





rs910416 <- qread("/media/work/gwis_test/hrt_ref_pm2/output/posthoc/dosage_chr7_87195962.qs") %>% 
  rename(snp = chr7_87195962_G_A_dose)

dachs <- filter(input_data, grepl("DACHS", study_gxe))

rs910416_dg <- glm(outcome ~ snp * hrt_ref_pm2 + age_ref_imp + pc1 + pc2 + pc3, 
                    data = inner_join(dachs, rs910416, 'vcfid'), family = 'binomial')

summary(rs910416_dg)



# stratified odds ratio for garcia paper
rs964293 <- qread("/media/work/gwis_test/hrt_ref_pm2/output/posthoc/dosage_chr20_52816717.qs")


covariates <- sort(c('age_ref_imp','study_gxe', 'pc1', 'pc2', 'pc3'))
path = glue("/media/work/gwis_test/{exposure}/")
snps <- c("20:52816717:C:A")

walk(snps, ~ fit_stratified_or(
  data_epi = input_data,
  exposure = exposure,
  snp = .x,
  hrc_version = hrc_version,
  covariates = covariates,
  path = glue("{path}/output/")))


iplot_wrapper(data_epi = input_data, exposure = exposure, hrc_version = hrc_version, snp = "20:52816717:C:A", covariates = covariates, path = glue("{path}/output"), flip = F)


# ================================================================== #
# ======= rmarkdown reports ---- 
gwis_report(exposure = exposure, 
            hrc_version = hrc_version, 
            covariates = covariates)


posthoc_report(exposure = exposure)

