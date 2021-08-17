#=============================================================================#
# FIGI GxE smk_aveday results
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
exposure = 'smk_aveday'
hrc_version = 'v2.3'
annotation_file <- 'gwas_200_ld_annotation_feb2021.txt'
covariates <- sort(c('age_ref_imp', 'sex', 'study_gxe', 'pc1', 'pc2', 'pc3'))
path = glue("/media/work/gwis_test/{exposure}/")


# input data
esubset <- readRDS(glue("/media/work/gwis_test/{exposure}/data/FIGI_{hrc_version}_gxeset_{exposure}_basic_covars_glm.rds")) %>% 
  pull(vcfid)

input_data <- readRDS(glue("/media/work/gwis_test/data/FIGI_{hrc_version}_gxeset_analysis_data_glm.rds")) %>% 
  filter(vcfid %in% esubset) %>%
  mutate(smk_aveday = smk_aveday / 10)
  # mutate(smk_aveday = scale(smk_aveday))






input_data %>% 
  mutate(smk_aveday = as.numeric(smk_aveday)) %>% 
  filter(smk_aveday != 0) %>% 
  mutate(study = fct_reorder(study, smk_aveday, .fun='median')) %>% 
  ggplot(aes(x=study, y=smk_aveday, fill = study)) + 
  geom_boxplot(outlier.color = 'red') + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 270))






#-----------------------------------------------------------------------------#
# main effects ----
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
# GxE additional analysis ---- 
#-----------------------------------------------------------------------------#

# output GxE models adjusted by different covariate sets
covariates_sets <- list(covariates, 
                        c(covariates, 'bmi5'), 
                        c(covariates, 'bmi5', 'smk_ever'), 
                        c(covariates, 'bmi5', 'smk_ever', 'fruitqc2', 'vegetableqc2'))

covariates_sets <- list(covariates)



# stratified by tumor and sex 
snps <- c("6:31917540:T:C", "8:138788813:C:A")
walk(snps, ~ fit_gxe_stratified(data_epi = input_data, exposure = exposure, snp = .x, covariates = covariates, method = "chiSqGxE", strata = 'sex', path = glue("{path}/output")))
walk(snps, ~ fit_gxe_stratified(data_epi = input_data, exposure = exposure, snp = .x, covariates = covariates, method = "chiSqGxE", strata = 'cancer_site_sum2', path = glue("{path}/output")))
walk(snps, ~ fit_gxe_stratified(data_epi = input_data, exposure = exposure, snp = .x, covariates = covariates, method = "chiSqGxE", strata = 'study_design', path = glue("{path}/output")))


snps <- c("6:32006886:G:A")
walk(snps, ~ fit_gxe_stratified(data_epi = input_data, exposure = exposure, snp = .x, covariates = covariates, method = "chiSq2df", strata = 'sex', path = glue("{path}/output")))
walk(snps, ~ fit_gxe_stratified(data_epi = input_data, exposure = exposure, snp = .x, covariates = covariates, method = "chiSq2df", strata = 'cancer_site_sum2', path = glue("{path}/output")))




# additional covariates is making association more significant.. let's generate for all suggestive hits
suggestive_gxe <- fread(glue("{path}/data/FIGI_{hrc_version}_gxeset_{exposure}_chiSqGxE_ldclump.clumped"))
walk(suggestive_gxe$SNP, ~ fit_gxe_covars(data_epi = input_data, exposure = exposure, snp = .x, covariates_list = covariates_sets, method = 'chiSqGxE', path = glue("{path}/output")))


# output RERI plots (can't install package on CARC yet)
snps <- c("6:31917540:T:C", "8:138788813:C:A")
walk(snps, ~ reri_wrapper(data_epi = input_data, exposure = exposure, snp = .x, covariates = covariates, path = glue("{path}/output")))


# recreate iplot wrapper (interaction plot)
# in order to use rescaled variable
iplot_wrapper(data_epi = input_data, 
              exposure = exposure, 
              hrc_version = hrc_version ,
              covariates = covariates, 
              snp = "8:138788813:C:A", 
              path = '/media/work/gwis_test/smk_aveday/output')


# create stratified odds ratios for smk_aveday scaled to per 10 cigarettes
fit_stratified_or(data_epi = input_data, exposure = exposure, snp = '6:31917540:T:C', hrc_version = hrc_version, covariates = covariates, dosage = F, path = glue("{path}/output"))

fit_stratified_or(data_epi = input_data, exposure = exposure, snp = '8:138788813:C:A', hrc_version = hrc_version, covariates = covariates, dosage = F, flip_allele = T, path = glue("{path}/output"))





# get MAF plots for suggestive hits 
suggestive_hits <- fread(glue("/media/work/gwis_test/{exposure}/data/FIGI_{hrc_version}_gxeset_{exposure}_chiSqGxE_ldclump.clumped")) %>% 
  pull(SNP)

hits <- c("6:31917540:T:C")
hits <- c("8:138788813:C:A")

walk(hits, ~ output_aaf_plot(dat = input_data, exposure, snp = .x))







# ================================================================== #
# ======= AAF check ---- 
# ================================================================== #




# chr6_31917540_T_C_dose
tmp <- qread(glue("/media/work/gwis_test/{exposure}/output/posthoc/dosage_chr6_31917540.qs")) %>% 
  inner_join(input_data, 'vcfid')

tmp %>% 
  summarise(total = n(), 
            study_aaf = sum(chr6_31917540_T_C_dose) / (total*2))

out <- tmp %>% 
  group_by(study_gxe) %>% 
  summarise(total = n(), 
            study_aaf = sum(chr6_31917540_T_C_dose) / (total*2)) %>% 
  arrange(study_aaf) %>% 
  mutate(study_gxe = fct_reorder(study_gxe, study_aaf))

ggplot(aes(x = study_gxe, y = study_aaf), data = out) + 
  geom_point() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 270)) + 
  xlab("Study") + 
  ylab("Alternate Allele Frequency")






# chr8_138788813_C_A_dose
tmp <- qread(glue("/media/work/gwis_test/{exposure}/output/posthoc/dosage_chr8_138788813.qs")) %>% 
  inner_join(input_data, 'vcfid')

tmp %>% 
  summarise(total = n(), 
            study_aaf = sum(chr8_138788813_C_A_dose) / (total*2))

out <- tmp %>% 
  group_by(study_gxe) %>% 
  summarise(total = n(), 
            study_aaf = sum(chr8_138788813_C_A_dose) / (total*2)) %>% 
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

main_effects_report(exposure = exposure, hrc_version = hrc_version, covariates = covariates, path = path)

gwis_report(exposure = exposure, 
            hrc_version = hrc_version, 
            covariates = covariates)

posthoc_report(exposure = exposure, 
               hrc_version = hrc_version,
               covariates = covariates,
               path = path)


posthoc_report(exposure = exposure)



