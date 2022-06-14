#=============================================================================#
# FIGI GxE asp_ref results
#=============================================================================#

# Setup -------------------------------------------------------------------
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
covariates <- sort(c('age_ref_imp', 'sex', 'study_gxe', 'pc1', 'pc2', 'pc3'))
path = glue("/media/work/gwis_test/{exposure}/")

# input data
esubset <- readRDS(glue("/media/work/gwis_test/{exposure}/data/FIGI_{hrc_version}_gxeset_{exposure}_basic_covars_glm.rds")) %>% 
  pull(vcfid)

input_data <- readRDS(glue("/media/work/gwis_test/data/FIGI_{hrc_version}_gxeset_analysis_data_glm.rds")) %>% 
  filter(vcfid%in% esubset) 




#-----------------------------------------------------------------------------#
# Main effects ------------------------------------------------------------
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



# having issue with study specific estimates... 
# covariates_full <- c(covariates, c('bmi', 'smk_ever', 'alcoholc', 'redmeatqc2'))
# covariates_full_meta <- sort(covariates_full[which(!covariates_full %in% c(paste0(rep('pc', 20), seq(1, 20)), "study_gxe"))])
# 
# create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_full_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "all", forest_height = 15)
# create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_full_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "proximal", forest_height = 13)
# create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_full_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "distal", forest_height = 13)
# create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_full_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "rectal", forest_height = 13)
# create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_full_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "female", forest_height = 13)
# create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_full_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "male", forest_height = 13)


# ------- stratified pooled analysis ------- #
pooled_analysis_glm(input_data, exposure = exposure, hrc_version = hrc_version, covariates = covariates, strata = 'sex', filename_suffix = paste0(paste0(sort(covariates), collapse = '_'), "_stratified_sex"), output_dir = glue("{path}/output/posthoc/"))

pooled_analysis_glm(input_data, exposure = exposure, hrc_version = hrc_version, covariates = covariates, strata = 'study_design', filename_suffix = paste0(paste0(sort(covariates), collapse = '_'), "_stratified_study_design"), output_dir = glue("{path}/output/posthoc/"))

pooled_analysis_multinom(input_data, exposure = exposure, hrc_version = hrc_version, covariates = covariates, filename_suffix = paste0(paste0(sort(covariates), collapse = '_'), "_stratified_cancer_site_sum2"), output_dir = glue("{path}/output/posthoc/"))



# with full covariates (easier to run)
covariates_full <- c(covariates, c('bmi', 'smk_ever', 'alcoholc', 'redmeatqc2'))

pooled_analysis_glm(input_data, exposure = exposure, hrc_version = hrc_version, covariates = covariates_full, strata = 'sex', filename_suffix = paste0(paste0(sort(covariates_full), collapse = '_'), "_stratified_sex"), output_dir = glue("{path}/output/posthoc/"))

pooled_analysis_glm(input_data, exposure = exposure, hrc_version = hrc_version, covariates = covariates_full, strata = 'study_design', filename_suffix = paste0(paste0(sort(covariates_full), collapse = '_'), "_stratified_study_design"), output_dir = glue("{path}/output/posthoc/"))

pooled_analysis_multinom(input_data, exposure = exposure, hrc_version = hrc_version, covariates = covariates_full, filename_suffix = paste0(paste0(sort(covariates_full), collapse = '_'), "_stratified_cancer_site_sum2"), output_dir = glue("{path}/output/posthoc/"))







# ------ let's output a cleaner meta-analysis plots of all the main effects, stratified, in a single image ----# ? 









# # main effects additional ----
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
# plot(inf.analysis, "influence")
# plot(inf.analysis, "baujat")
# plot(inf.analysis, "es")
# plot(inf.analysis, "i2")



#-----------------------------------------------------------------------------#
# main effects (SNP) ----
#-----------------------------------------------------------------------------#
tmp <- qread("/media/work/gwis_test/asp_ref/output/posthoc/dosage_chr6_12577203.qs")

tmp2 <- inner_join(input_data, tmp, 'vcfid')


model <- glm(outcome ~ asp_ref + chr6_12577203_T_C_dose + age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe, family = 'binomial', data = tmp2)

out <- tidy(model, conf.int = T, exponentiate = T)


exp(0.00285 - (1.96*0.0238))
exp(0.00285 + (1.96*0.0238))


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
# 



#-----------------------------------------------------------------------------#
# GxE additional analysis ---- 
#-----------------------------------------------------------------------------#


# Additional adjustment covariates
covariates_sets <- list(covariates,
                        c(covariates, 'bmi'), 
                        c(covariates, 'bmi', 'smk_ever'), 
                        c(covariates, 'bmi', 'smk_ever', 'alcoholc', 'redmeatqc2'))

fit_gxe_covars(data = input_data, exposure = exposure, snp = "6:12577203:T:C", covariates_list = covariates_sets, method = 'chiSqGxE', path = glue("{path}/output"))
fit_gxe_covars(data = input_data, exposure = exposure, snp = "5:40252294:C:T", covariates_list = covariates_sets, method = 'chiSqGxE', path = glue("{path}/output"))


# output RERI plots
snps <- c("6:12577203:T:C")
walk(snps, ~ reri_wrapper(data_epi = input_data, exposure = exposure, snp = .x, covariates = covariates, path = glue("{path}/output")))

# need to flip exposure
tmp_flip <- input_data %>% 
  mutate(asp_ref = fct_relevel(asp_ref, "Yes"))
snps <- c("5:40252294:C:T")
walk(snps, ~ reri_wrapper(data_epi = tmp_flip, exposure = exposure, snp = .x, covariates = covariates, path = glue("{path}/output")))


# stratified odds ratio table with additional covariates
# covariates <- c(covariates, 'bmi', 'smk_ever', 'alcoholc', 'redmeatqc2')




snps <- c( "6:12577203:T:C", "5:40252294:C:T")
walk(snps, ~ fit_stratified_or(
  data_epi = input_data,
  exposure = exposure,
  snp = .x,
  hrc_version = hrc_version,
  covariates = covariates,
  path = glue("{path}/output/")))


covariates_mult <- c(covariates, 'bmi', 'smk_ever', 'alcoholc', 'redmeatqc2')
snps <- c( "6:12577203:T:C", "5:40252294:C:T")
walk(snps, ~ fit_stratified_or(
  data_epi = input_data,
  exposure = exposure,
  snp = .x,
  hrc_version = hrc_version,
  covariates = covariates_mult,
  path = glue("{path}/output/")))


# SNP information
snp_info <- qread("/media/work/FIGI_RsqEstimate_chrALL.qs") %>% 
  filter(id %in%  c("6:12577203:T:C", "6:32560631:C:T", "5:40252294:C:T"))

ff <- do.call(c, lapply(strsplit(snps, split = ":"), function(x) paste(x[1], as.numeric(x[2])-1, x[2],  sep = ':')))

snp_mart = useMart(biomart = "ENSEMBL_MART_SNP", 
                   host    = "grch37.ensembl.org", 
                   path    = "/biomart/martservice", 
                   dataset = "hsapiens_snp")

rsid = getBM(attributes = c("refsnp_id", "allele", "chr_name", "chrom_end"),
             filters = c("chromosomal_region"),
             values = ff, mart=snp_mart)


snp_info$rsid <- rsid$refsnp_id

snp_info_out <- snp_info %>% 
  separate(SNP, into = c("chr", "bp", "REF", "ALT"), remove= F) %>%
  rename(ALT_AF = GxESet_AltAlleleFreq, 
         MAF = maf, 
         Imputation_Rsq = GxESet_Rsq) %>% 
  dplyr::select(SNP, rsid, REF, ALT, ALT_AF , MAF, Imputation_Rsq )

saveRDS(snp_info_out, glue("{path}/output/posthoc/gwis_snp_info.rds"))





# ---- MAF ---- #
snps <- c("6:12577203:T:C", "6:32560631:C:T", "9:22103341:T:G", "2:137374734:T:C", "5:40252294:C:T")
walk(snps, ~ create_aaf_study_plot(data = input_data, exposure = exposure, hrc_version = hrc_version, snp = .x, path = path))






# would results change if you remove mecc (note wald statistic)
# model1 <- glm(glue("outcome ~ {exposure}*chr6_12577203_T_C_dose + {glue_collapse(covariates, sep = '+')}"), data = chr6_125, family = 'binomial')
# summary(model1)
# 
# chr6_125b <- dplyr::filter(chr6_125, !grepl("MECC", study_gxe))
# model2 <- glm(glue("outcome ~ {exposure}*chr6_12577203_T_C_dose + {glue_collapse(covariates, sep = '+')}"), data = chr6_125b, family = 'binomial')
# summary(model2)
# 
# model3 <- glm(glue("outcome ~ chr6_12577203_T_C_dose + {glue_collapse(covariates, sep = '+')}"), data = chr6_125, family = 'binomial')
# summary(model3)



# would results change if you remove UKB (note wald statistic)
# yes, as expected
dose <- qread("/media/work/gwis_test/asp_ref/output/posthoc/dosage_chr6_32560631.qs")
out <- inner_join(dose, input_data, 'vcfid')

model1 <- glm(glue("outcome ~ {exposure}*chr6_32560631_C_T_dose + {glue_collapse(covariates, sep = '+')}"), data = out, family = 'binomial')
summary(model1)

model1 <- glm(glue("outcome ~ {exposure}*chr6_32560631_C_T_dose + {glue_collapse(covariates, sep = '+')}"), data = out %>% filter(study_gxe != "UKB_1"), family = 'binomial')
summary(model1)


model1 <- glm(glue("outcome ~ {exposure}+chr6_32560631_C_T_dose + {glue_collapse(covariates, sep = '+')}"), data = out, family = 'binomial')
summary(model1)

model1 <- glm(glue("outcome ~ {exposure}+chr6_32560631_C_T_dose + {glue_collapse(covariates, sep = '+')}"), data = out %>% filter(study_gxe != "UKB_1"), family = 'binomial')
summary(model1)

ukb <- out %>% 
  filter(study_gxe == "UKB_1") %>% 
  mutate(chr6_32560631_C_T_dose = 2-chr6_32560631_C_T_dose)
  
out2 <- filter(out, study_gxe !="UKB_1") %>% 
  bind_rows(ukb)


model1 <- glm(glue("outcome ~ {exposure}*chr6_32560631_C_T_dose + {glue_collapse(covariates, sep = '+')}"), data = out2, family = 'binomial')
summary(model1)


model1 <- glm(glue("outcome ~ {exposure}+chr6_32560631_C_T_dose + {glue_collapse(covariates, sep = '+')}"), data = out2, family = 'binomial')
summary(model1)


covariates
covariates2 <- c("age_ref_imp", "pc1", "pc2", "pc3", "sex")
model1 <- glm(glue("outcome ~ {exposure}*chr6_32560631_C_T_dose + {glue_collapse(covariates2, sep = '+')}"), data = out %>% filter(study_gxe == "UKB_1"), family = 'binomial')
summary(model1)



# is there a way UKB AF is RIGHT (is that population that different - unlikely)
# - 



# 
# chr6_325b <- dplyr::filter(chr6_325, !grepl("UKB", study_gxe))
# model2 <- glm(glue("outcome ~ {exposure}*chr6_32560631_C_T_dose + {glue_collapse(covariates, sep = '+')}"), data = chr6_325b, family = 'binomial')
# summary(model2)



# I'm not sure why I did this ... 
model1 <- glm(glue("outcome ~ {exposure}*chr5_40252294_C_T_dose + {glue_collapse(covariates, sep = '+')}"), data = chr5, family = 'binomial')
summary(model1)

chr5b <- dplyr::filter(chr5, !grepl("REACH", study_gxe))
model2 <- glm(glue("outcome ~ {exposure}*chr5_40252294_C_T_dose + {glue_collapse(covariates, sep = '+')}"), data = chr5b, family = 'binomial')
summary(model2)



# Replicate Nan et al findings --------------------------------------------

# replicate Nan et al (GxE or case-only)
snps <- c("12:17444733:A:T", "12:17488764:A:T", "15:82229999:A:C")

## MAF plots..
walk(snps, ~ create_aaf_study_plot(data = input_data, exposure = exposure, hrc_version = hrc_version, snp = .x, path = path))


## stratified odds ratio plots
walk(snps, ~ fit_stratified_or(
  data_epi = input_data,
  exposure = exposure,
  snp = .x,
  hrc_version = hrc_version,
  covariates = covariates,
  path = glue("{path}/output/")))


## forest plot of GxE term , by study 
dose <- qread("/media/work/gwis_test/asp_ref/output/posthoc/dosage_chr12_17444733.qs")
tmp <- inner_join(input_data, dose, 'vcfid')



## fit models
dose <- 





# additional covariate for the chr6 finding -------------------------------


dose <- qread("/media/work/gwis_test/asp_ref/output/posthoc/dosage_chr6_12577203.qs")
out <- inner_join(input_data, dose, 'vcfid')

mod_base <- glm(outcome ~ chr6_12577203_T_C_dose + asp_ref + age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe, data = out, family = 'binomial')
mod_alt <- glm(outcome ~ chr6_12577203_T_C_dose * asp_ref + age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe, data = out, family = 'binomial') # GxE OR = 2.776e-01

lrtest(mod_base, mod_alt)

# BMI, tobacco smoking, alcohol use, and red meat
mod_base <- glm(outcome ~ chr6_12577203_T_C_dose + asp_ref + age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe + bmi + smk_ever + alcoholc + redmeatqc2, data = out, family = 'binomial')
mod_alt <- glm(outcome ~ chr6_12577203_T_C_dose * asp_ref + age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe + bmi + smk_ever + alcoholc + redmeatqc2, data = out, family = 'binomial') # GxE OR = 2.557e-01

lrtest(mod_base, mod_alt)




2.557e-01/2.776e-01

exp(2.557e-01)/exp(2.776e-01)




# -------------------- chr5 finding GxExSex -------------------- #

dose <- qread("/media/work/gwis_test/asp_ref/output/posthoc/dosage_chr5_40252294.qs")
out <- inner_join(input_data, dose, 'vcfid')

mod_base <- glm(outcome ~ chr5_40252294_C_T_dose * asp_ref + age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe, data = out, family = 'binomial')
summary(mod_base)

mod_base <- glm(outcome ~ chr5_40252294_C_T_dose * asp_ref * sex + age_ref_imp + pc1 + pc2 + pc3 + study_gxe, data = out, family = 'binomial')
summary(mod_base)

# ================================================================== #
# ======= rmarkdown reports ---- 
# ================================================================== #
main_effects_report(exposure = exposure, hrc_version = hrc_version, covariates = covariates, path = path)

gwis_report(exposure = exposure, 
            hrc_version = hrc_version, 
            covariates = covariates)

posthoc_report(exposure = exposure)


