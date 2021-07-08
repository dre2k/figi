#=============================================================================#
# FIGI GxE alcoholc_moderate results
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
exposure = 'alcoholc_moderate'
exposure = "alcoholc_heavy_vs_moderate"
hrc_version = 'v2.3'
annotation_file <- 'gwas_200_ld_annotation_feb2021.txt'
path = glue("/media/work/gwis_test/{exposure}/")

covariates <- sort(c('age_ref_imp', 'sex', 'energytot_imp', 'study_gxe', 'pc1', 'pc2', 'pc3'))
covariates_set1 <- covariates
covariates_set2 <- c(covariates_set1, "smk_ever", "bmi", "famhx1")
covariates_list <- list(covariates_set1, covariates_set2)
mod <- 'age_ref_imp+sex+energytot_imp+pc1+pc2+pc3+study_gxe'


# input data
esubset <- readRDS(glue("{path}/data/FIGI_{hrc_version}_gxeset_{exposure}_basic_covars_glm.rds")) %>% pull(vcfid)
input_data <- readRDS(glue("/media/work/gwis/data/FIGI_EpiData/FIGI_{hrc_version}_gxeset_analysis_data_glm.rds")) %>% 
  filter(vcfid%in% esubset)



#-----------------------------------------------------------------------------#
# main effects ----
#-----------------------------------------------------------------------------#
# ------ meta-analysis ------ #
meta_analysis_execute(dataset = input_data, exposure = exposure, hrc_version = hrc_version, covariates = covariates, output_dir = output_dir, filename_suffix2 = "")
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
# 
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
# SNP followup (models and plots) ---- 
#-----------------------------------------------------------------------------#

# make sure 1-28g/d is the reference group !
figi <- posthoc_input(exposure, hrc_version, glue('gwis_sig_results_output_{exposure}.rds')) %>% 
  mutate(alcoholc_moderate = fct_relevel(alcoholc_moderate, "1-28g/d", "nondrinker"), 
         alcoholc_moderate = as.numeric(alcoholc_moderate)) # better to make it numeric.. 

figi <- posthoc_input(exposure, hrc_version, glue('gwis_sig_results_output_{exposure}.rds')) %>% 
  mutate(alcoholc_moderate = fct_relevel(alcoholc_moderate, "1-28g/d", "nondrinker"))


snps <- readRDS(glue(output_dir, "gwis_sig_results_input_{exposure}.rds"))
table(snps$method)


# ------- GxE 1DF suggestive findings ------- #
snps_filter <- snps %>% 
  dplyr::filter(method == glue("manhattan_chiSqGxE_{exposure}_clump_df")) %>% 
  dplyr::mutate(snps = paste0('chr', gsub("\\:", "\\_", SNP))) %>% 
  pull(snps)

posthoc_run_models(snps_filter, 'chiSqGxE')
posthoc_run_reri_plot(snps_filter)
posthoc_create_plots(snps_filter, 'chiSqGxE')





# ------- chr10_101476905_G_A ------- #
snps_filter <- snps %>% 
  dplyr::filter(method == glue("manhattan_chiSqGxE_{exposure}_clump_df"),
                SNP == "10:101476905:G:A") %>% 
  dplyr::mutate(snps = paste0('chr', gsub("\\:", "\\_", SNP))) %>% 
  pull(snps)

walk(snps_filter, ~ fit_stratified_or_tmp(ds = figi, exposure = exposure, snp = .x, hrc_version = hrc_version, covariates = covariates_set1, mod = mod, dosage = F, output_dir = output_dir))






# ------- 2DF/3DF clumped removing GWAS and or EG (to check 'novel' D|G results) ------- #
# - locuszoom plots
# - conditional analyses??? (more involved.. but maybe worthwhile if they request it)

snps_2df <- readRDS(glue(output_dir, "manhattan_chiSq2df_{exposure}_no_gwas_clump_df.rds")) %>% 
  dplyr::mutate(snps = paste0('chr', gsub("\\:", "\\_", SNP))) %>% 
  pull(snps)

posthoc_run_models(snps_2df, 'chiSq2df')
posthoc_run_reri_plot(snps_2df)
posthoc_create_plots(snps_2df, 'chiSqG')

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

posthoc_run_models(snps_filter[13], 'chiSqGxE')
posthoc_run_reri_plot(snps_filter[13])
posthoc_create_plots(snps_filter[13], 'chiSqGxE')


# ------- chr9_97251034_C_T ------- #
# previous gwis finding.. 

figi <- posthoc_input(exposure, hrc_version, glue('gwis_sig_results_additional_output_{exposure}.rds')) %>% 
  mutate(alcoholc_moderate = fct_relevel(alcoholc_moderate, "1-28g/d", "nondrinker"))


figi <- posthoc_input(exposure, hrc_version, glue('gwis_sig_results_additional_output_{exposure}.rds')) %>% 
  mutate(alcoholc_moderate = fct_relevel(alcoholc_moderate, "nondrinker", "1-28g/d"))


snps <- readRDS(glue(output_dir, "gwis_sig_results_additional_input_{exposure}.rds"))


snps_filter <- snps %>% 
  dplyr::mutate(snps = paste0('chr', gsub("\\:", "\\_", SNP))) %>% 
  pull(snps)


# make sure SNPs are NOT flipped
walk(snps_filter, ~ fit_gxe_covars_tmp(ds = figi, exposure, .x, covariates_list, 'chiSqGxE', output_dir = output_dir))





# repeat above but keep original alcoholc_moderate coding

figi <- posthoc_input(exposure, hrc_version, glue('gwis_sig_results_additional_output_{exposure}.rds')) %>% 
  mutate(alcoholc_moderate = fct_relevel(alcoholc_moderate, "nondrinker", "1-28g/d"))

figi <- posthoc_input(exposure, hrc_version, glue('gwis_sig_results_additional_output_{exposure}.rds')) %>% 
  mutate(alcoholc_moderate = fct_relevel(alcoholc_moderate, "1-28g/d", "nondrinker"))

snps <- readRDS(glue(output_dir, "gwis_sig_results_additional_input_{exposure}.rds"))

#chr9_97251034_C_T
snps_filter <- snps %>% 
  dplyr::mutate(snps = paste0('chr', gsub("\\:", "\\_", SNP))) %>% 
  pull(snps)

walk(snps_filter, ~ fit_gxe_covars_tmp(ds = figi, exposure, .x, covariates_list, 'chiSqGxE', output_dir = output_dir))


# 1-28g/d nondrinker 
# 35637      29300
# 






# ================================================================== #
# ======= AAF check ---- 
# ================================================================== #




# for rs2300985 chr10_101476905_G_A
tmp <- qread(glue("/media/work/gwis_test/{exposure}/output/posthoc/dosage_chr10_101476905.qs")) %>% 
  inner_join(input_data, 'vcfid')

tmp %>% 
  summarise(total = n(), 
            study_aaf = sum(chr10_101476905_G_A_dose) / (total*2))

out <- tmp %>% 
  group_by(study_gxe) %>% 
  summarise(total = n(), 
            study_aaf = sum(chr10_101476905_G_A_dose) / (total*2)) %>% 
  arrange(study_aaf) %>% 
  mutate(study_gxe = fct_reorder(study_gxe, study_aaf))

ggplot(aes(x = study_gxe, y = study_aaf), data = out) + 
  geom_point() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 270)) + 
  xlab("Study") + 
  ylab("Alternate Allele Frequency")






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









#-----------------------------------------------------------------------------#
# SNP followup (manuscript editing) ---- 
#-----------------------------------------------------------------------------#

source("~/git/figifs/R/01_process.R")
source("~/git/figifs/R/02_plots.R")
source("~/git/figifs/R/03_posthoc.R")
source("~/git/figifs/R/03_posthoc_iplot.R")
source("~/git/figifs/R/03_posthoc_stratified_or.R")
source("~/git/figifs/R/05_functional_annotation.R")


# models adjusted by BMI
covariates_sets <- list(covariates, 
                        c(covariates, 'bmi'), 
                        c(covariates, 'bmi', 'famhx1', 'diab', 'educ', 'smk_ever', 'fruitqc2', 'vegetableqc2'), 
                        c(covariates, 'bmi', 'famhx1', 'diab', 'educ', 'smk_ever', 'fruitqc2', 'vegetableqc2',  'methrswklns'))


covariates_sets <- list(covariates, 
                        c(covariates, 'bmi', 'diab', 'educ', 'smk_ever', 'redmeatqc2', 'fruitqc2', 'vegetableqc2'), 
                        c(covariates, 'bmi', 'diab', 'educ', 'smk_ever', 'redmeatqc2', 'fruitqc2', 'vegetableqc2', 'methrswklns'),
                        c(covariates, 'bmi', 'diab', 'educ', 'smk_ever', 'redmeatqc2', 'fruitqc2', 'vegetableqc2', 'hrt_ref_pm'), 
                        c(covariates, 'bmi', 'diab', 'educ', 'smk_ever', 'redmeatqc2', 'fruitqc2', 'vegetableqc2', 'hrt_ref_pm', 'methrswklns'))


tmp <- qread(glue("{path}/output/posthoc/dosage_chr10_101476905.qs")) %>% 
  inner_join(input_data, 'vcfid')

fit_gxe_covars(data_epi = input_data, exposure = exposure, snp = "10:101476905:G:A", covariates_list = covariates_sets, method = 'chiSqGxE', path = glue("{path}/output"))












# double check meta-analysis findings
tmp <- qread(glue("/media/work/gwis_test/{exposure}/output/posthoc/dosage_chr10_101476905.qs")) %>% 
  inner_join(input_data, 'vcfid')

data_epi = tmp
exposure
snp = "chr10_101476905_G_A_dose"
covariates
hrc_version
path
forest_height = 17
forest_width = 8.5
funnel_height = 8
funnel_width = 8.5
strata = 'all'
categorical = T


# flip allele coding

# tmp <- tmp %>% 
#   mutate(chr10_101476905_A_G_flip =  abs(2-chr10_101476905_G_A_dose))


create_forest_plot_gxe(data_epi = tmp, exposure = exposure, snp = 'chr10_101476905_G_A_dose', covariates = covariates, hrc_version = hrc_version, path = path)

stratified_or <- readRDS("/media/work/gwis_test/alcoholc_moderate/output/posthoc/stratified_oddsratio_alcoholc_moderate_v2.3_chr10_101476905_G_A_age_ref_imp_energytot_imp_pc1_pc2_pc3_sex_study_gxe.rds")

# for rs2300985 chr10_101476905_G_A
# make sure A is the reference allele

tmp <- qread(glue("{path}/output/posthoc/dosage_chr10_101476905.qs")) %>% 
  inner_join(input_data, 'vcfid')

model <- glm(glue("outcome ~ {exposure}*chr10_101476905_G_A_dose + {glue_collapse(covariates, sep = '+')}"), data = tmp, family = 'binomial')
summary(model)

exp(coef(model))



# fit the model as you usually do 
covariates_sets = list(covariates)
fit_gxe_covars(data_epi = input_data, exposure = exposure, snp = "10:101476905:G:A", covariates_list = covariates_sets, method = 'chiSqGxE', path = glue("{path}/output"), flip_allele = T)




fit_stratified_or(
  data_epi = input_data,
  exposure = exposure,
  snp = '10:101476905:G:A',
  hrc_version = hrc_version,
  covariates = covariates,
  path = glue("{path}/output"))


x <- readRDS("/media/work/gwis_test/alcoholc_moderate/output/posthoc/stratified_oddsratio_alcoholc_moderate_v2.3_chr10_101476905_G_A_age_ref_imp_energytot_imp_pc1_pc2_pc3_sex_study_gxe.rds")








# p-value for main effects

tmp <- qread(glue("/media/work/gwis_test/alcoholc_moderate/output/posthoc/dosage_chr10_101351704.qs")) %>% 
  inner_join(input_data, 'vcfid')

model_base <- glm(outcome ~ age_ref_imp + sex + energytot_imp + pc1 + pc2 + pc3 + study_gxe, data = tmp, family = 'binomial')
model_g <- glm(outcome ~ chr10_101351704_A_G_dose + age_ref_imp + sex + energytot_imp + pc1 + pc2 + pc3 + study_gxe, data = tmp, family = 'binomial')

lrtest(model_base, model_g)



tmp <- qread(glue("/media/work/gwis_test/alcoholc_moderate/output/posthoc/dosage_chr10_101476905.qs")) %>% 
  inner_join(input_data, 'vcfid')

model_base <- glm(outcome ~ age_ref_imp + sex + energytot_imp + pc1 + pc2 + pc3 + study_gxe, data = tmp, family = 'binomial')
model_g <- glm(outcome ~ chr10_101476905_G_A_dose + age_ref_imp + sex + energytot_imp + pc1 + pc2 + pc3 + study_gxe, data = tmp, family = 'binomial')

lrtest(model_base, model_g)





tmp <- qread(glue("/media/work/gwis_test/alcoholc_moderate/output/posthoc/dosage_chr10_101353285.qs")) %>% 
  inner_join(input_data, 'vcfid')

model_base <- glm(outcome ~ age_ref_imp + sex + energytot_imp + pc1 + pc2 + pc3 + study_gxe, data = tmp, family = 'binomial')
model_g <- glm(outcome ~ chr10_101353285_C_T_dose + age_ref_imp + sex + energytot_imp + pc1 + pc2 + pc3 + study_gxe, data = tmp, family = 'binomial')

lrtest(model_base, model_g)





#### poshtoc analysis after flipping alleles (stratified models by sex study_design tumor site, iplot)

source('/media/work/gwis_test/figifs/R/03_posthoc.R')

# for rs2300985 chr10_101476905_G_A
tmp <- qread(glue("/media/work/gwis_test/alcoholc_moderate/output/posthoc/dosage_chr10_101476905.qs")) %>% 
  inner_join(input_data, 'vcfid')


iplot_wrapper_temp(data_epi = tmp, exposure = exposure, hrc_version = hrc_version, snp = "10:101476905:G:A", covariates = covariates, path = "~/Dropbox/", flip = T)

fit_gxe_stratified_tmp(data_epi = tmp, exposure = exposure, snp = "10:101476905:G:A", covariates = covariates, strata = 'sex', method = 'chiSqGxE', path = "~/Dropbox/", flip = T)
fit_gxe_stratified_tmp(data_epi = tmp, exposure = exposure, snp = "10:101476905:G:A", covariates = covariates, strata = 'cancer_site_sum2', method = 'chiSqGxE', path = "~/Dropbox/", flip = T)
fit_gxe_stratified_tmp(data_epi = tmp, exposure = exposure, snp = "10:101476905:G:A", covariates = covariates, strata = 'study_design', method = 'chiSqGxE', path = "~/Dropbox/", flip = T)








# ------------------------------- #
# question - allelic dosages vs imputed genotypes
# they're not. i mean it must be some sort of rounding perhaps, since the odsage is liek a weigted avg of the genotype probabiltiies
# ---------------------------------# 


model1 <- glm(outcome ~ alcoholc_moderate*chr10_101476905_G_A_dose + age_ref_imp + sex + energytot_imp + pc1 + pc2 + pc3 + study_gxe, data = tmp, family = 'binomial')
summary(model1)


model2 <- glm(outcome ~ alcoholc_moderate*chr10_101476905_G_A_p1 +  alcoholc_moderate*chr10_101476905_G_A_p2 + age_ref_imp + sex + energytot_imp + pc1 + pc2 + pc3 + study_gxe, data = tmp, family = 'binomial')
summary(model2)


model3 <- glm(outcome ~ alcoholc_moderate*chr10_101476905_G_A_p1 + age_ref_imp + sex + energytot_imp + pc1 + pc2 + pc3 + study_gxe, data = tmp, family = 'binomial')
summary(model3)



x <- readRDS("/media/work/gwis_test/alcoholc_moderate/output/posthoc/stratified_oddsratio_alcoholc_moderate_v2.3_chr10_101476905_G_A_age_ref_imp_energytot_imp_pc1_pc2_pc3_sex_study_gxe.rds")


reri_wrapper(data_epi = input_data, 
            exposure = exposure, 
              snp = "10:101476905:G:A", covariates = covariates, path = glue("{path}/output/"))




# --------------------------- #
# get allele frequencies for main findings
# report for each alcoholc subset
# ---------------------------- #

# let's use moderate and heavy to obtain VCFIDs, this is because I want to ensure i have the same study inclusion/exclusion as Kristina's paper
alcoholc_moderate_vcfid <- readRDS("/media/work/gwis_test/alcoholc_moderate/data/FIGI_v2.3_gxeset_alcoholc_moderate_basic_covars_glm.rds") %>% 
  pull(vcfid)
alcoholc_heavy_vs_moderate_vcfid <- readRDS("/media/work/gwis_test/alcoholc_heavy_vs_moderate/data/FIGI_v2.3_gxeset_alcoholc_heavy_vs_moderate_basic_covars_glm.rds") %>% 
  pull(vcfid)

# N = 74099
keep <- unique(c(alcoholc_moderate_vcfid, alcoholc_heavy_vs_moderate_vcfid))

input_data <-  readRDS(glue("/media/work/gwis_test/data/FIGI_v2.3_gxeset_analysis_data_glm.rds")) %>% 
  filter(vcfid %in% keep)


# allele frequencies


tmp <- qread(glue("/media/work/gwis_test/alcoholc/output/posthoc/dosage_chr10_101353285.qs")) %>% 
  inner_join(input_data, 'vcfid')

tmp %>% 
  summarise(total = n(), 
            study_aaf = sum(chr10_101353285_C_T_dose) / (total*2))

tmp %>% 
  filter(alcoholc == "nondrinker") %>% 
  summarise(total = n(), 
            study_aaf = sum(chr10_101353285_C_T_dose) / (total*2))


tmp %>% 
  filter(alcoholc == "1-28g/d") %>% 
  summarise(total = n(), 
            study_aaf = sum(chr10_101353285_C_T_dose) / (total*2))


tmp %>% 
  filter(alcoholc == ">28g/d" ) %>% 
  summarise(total = n(), 
            study_aaf = sum(chr10_101353285_C_T_dose) / (total*2))





# --------------------------------- #
# try to estimate how many people you excluded in thsi analysis
# #---------------------------------- #


input_data <-  readRDS(glue("/media/work/gwis_test/data/FIGI_v2.3_gxeset_analysis_data_glm.rds")) %>% 
  filter(!is.na(alcoholc)) %>% 
  filter(!grepl("PPS", study_gxe))
  # filter(vcfid %in% keep)

table(input_data$study_gxe, input_data$outcome)




final <-  readRDS(glue("/media/work/gwis_test/data/FIGI_v2.3_gxeset_analysis_data_glm.rds")) %>% 
  filter(vcfid %in% keep)

table(final$study_gxe, final$outcome)



# --------------------------------- #
# how many people were excluded based on ancestry
# #---------------------------------- #

gwas <- readRDS("/media/work/gwis_test/data/FIGI_v2.3_GWAS.rds")

gxe <- readRDS("/media/work/gwis_test/data/FIGI_v2.3_gxeset_analysis_data_glm.rds")
gxe <- filter(gwas, gxe == 1)



figi_gxe2 <- figi %>% 
  dplyr::filter(gxe == 1, 
                !study_gxe %in% c("ColoCare_1", "GALEON", "MOFFITT", "NFCCR_1", "NGCCS"))

x <- figi_gxe2 %>% 
  filter(!is.na(alcoholc))

noneur <- filter(x, is.na(EUR_subset ))



# here it is - get studies that are included in the alcoholc_moderate  analysis, and then subset without  excluding EUR. might have to recreate some variables that's all

alcoholc_mod_studies <- readRDS("/media/work/gwis_test/alcoholc_moderate/data/FIGI_v2.3_gxeset_alcoholc_moderate_basic_covars_glm.rds") %>% 
  pull(study_gxe) %>% 
  unique(.) %>% 
  as.character(.)



# dataset (see untitled2 lol)


tmp <- qread(glue("/media/work/gwis_test/{exposure}/output/posthoc/dosage_chr10_101476905.qs")) %>% 
  inner_join(figi_gxe, 'vcfid') %>% 
  filter(study_gxe %in% alcoholc_mod_studies)


# main effects 
model_dg_base <- glm(outcome ~                            age_ref_imp + sex + energytot_imp + pc1 + pc2 + pc3 + study_gxe, data = tmp, family = 'binomial')
model_dg      <- glm(outcome ~ chr10_101476905_G_A_dose + age_ref_imp + sex + energytot_imp + pc1 + pc2 + pc3 + study_gxe, data = tmp, family = 'binomial')

lrtest(model_dg_base, model_dg)

# GxE
model_gxe_base <- glm(outcome ~ chr10_101476905_G_A_dose + alcoholc_moderate + age_ref_imp + sex + energytot_imp + pc1 + pc2 + pc3 + study_gxe, data = tmp, family = 'binomial')
model_gxe      <- glm(outcome ~ chr10_101476905_G_A_dose * alcoholc_moderate + age_ref_imp + sex + energytot_imp + pc1 + pc2 + pc3 + study_gxe, data = tmp, family = 'binomial')

summary(model_gxe)

lrtest(model_gxe_base, model_gxe)














# -------------------------------------------------------------------------- #
# estimate how many age and energy values were imputed ..
# -------------------------------------------------------------------------- #

# (repeated from above)
# let's use moderate and heavy to obtain VCFIDs, this is because I want to ensure i have the same study inclusion/exclusion as Kristina's paper
alcoholc_moderate_vcfid <- readRDS("/media/work/gwis_test/alcoholc_moderate/data/FIGI_v2.3_gxeset_alcoholc_moderate_basic_covars_glm.rds") %>% 
  pull(vcfid)
alcoholc_heavy_vs_moderate_vcfid <- readRDS("/media/work/gwis_test/alcoholc_heavy_vs_moderate/data/FIGI_v2.3_gxeset_alcoholc_heavy_vs_moderate_basic_covars_glm.rds") %>% 
  pull(vcfid)

# N = 74099
keep <- unique(c(alcoholc_moderate_vcfid, alcoholc_heavy_vs_moderate_vcfid))

input_data <-  readRDS(glue("/media/work/gwis_test/data/FIGI_v2.3_gxeset_analysis_data_glm.rds")) %>% 
  filter(vcfid %in% keep)



any(is.na(input_data$energytot))
summary(input_data$energytot)
summary(input_data$energytot_imp)

impute_energy <- filter(input_data, is.na(energytot))



any(is.na(input_data$age_ref))
summary(input_data$age_ref)
summary(input_data$)

impute_energy <- filter(input_data, is.na(energytot))





# -------------------------------------------------- 
# assess alcoholc ~ covariate association
# -----------------------------------------------


quick_function <- function(covar) {
  model <- glm(glue("alcoholc_moderate ~ {glue_collapse(covariates, sep = '+')} + {covar}"), data = input_data, family = 'binomial')
  print(summary(model)); print(data.frame(tidy(model, exponentiate = T)))
}

sink("~/Dropbox/Working/alcoholc_moderate_covar_association2.txt")
covar <- c('bmi', 'famhx1', 'diab', 'educ', 'smk_ever', 'redmeatqc2', 'fruitqc2', 'vegetableqc2', 'hrt_ref_pm', 'methrswklns')
walk(covar, ~ quick_function(.x))
sink()

covariates_sets <- list(covariates, 
                        c(covariates, 'bmi', 'diab', 'educ', 'smk_ever', 'redmeatqc2', 'fruitqc2', 'vegetableqc2'), 
                        c(covariates, 'bmi', 'diab', 'educ', 'smk_ever', 'redmeatqc2', 'fruitqc2', 'vegetableqc2', 'methrswklns'),
                        c(covariates, 'bmi', 'diab', 'educ', 'smk_ever', 'redmeatqc2', 'fruitqc2', 'vegetableqc2', 'hrt_ref_pm'), 
                        c(covariates, 'bmi', 'diab', 'educ', 'smk_ever', 'redmeatqc2', 'fruitqc2', 'vegetableqc2', 'hrt_ref_pm', 'methrswklns'))
                        