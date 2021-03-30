#=============================================================================#
# FIGI GxE aspirin results
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
exposure = 'aspirin'
hrc_version = 'v2.4'
# hrc_version = 'v3.0'
output_dir = paste0("/media/work/gwis/posthoc/", exposure, "/")
# annotation_file <- 'gwas_141_ld_annotation_july2020.txt'
annotation_file <- 'gwas_200_ld_annotation_feb2021.txt'

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

gxe <- readRDS(glue("~/FIGI_{hrc_version}_gxeset_{exposure}_basic_covars_gxescan_results.rds")) %>% 
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

source(glue("/media/work/git/figifs/R/01_process.R"))
source(glue("/media/work/git/figifs/R/02_plots.R"))
source(glue("/media/work/git/figifs/R/03_posthoc.R"))
source(glue("/media/work/git/figifs/R/03_posthoc_iplot.R"))
source(glue("/media/work/git/figifs/R/03_posthoc_stratified_or.R"))

# output GxE models adjusted by different covariate sets
covariates_sets <- list(covariates, 
                        c(covariates, 'bmi5'), 
                        c(covariates, 'bmi5', 'smk_ever'), 
                        c(covariates, 'bmi5', 'smk_ever', 'fruitqc2', 'vegetableqc2'))

fit_gxe_covars(data_epi = input_data, exposure = exposure, snp = "13:47191972:G:A", covariates_list = covariates_sets, method = 'chiSqGxE', path = glue("{path}/output"))


# additional covariates is making association more significant.. let's generate for all suggestive hits
suggestive_gxe <- fread(glue("{path}/data/FIGI_{hrc_version}_gxeset_{exposure}_chiSqGxE_ldclump.clumped"))
walk(suggestive_gxe$SNP, ~ fit_gxe_covars(data_epi = input_data, exposure = exposure, snp = .x, covariates_list = covariates_sets, method = 'chiSqGxE', path = glue("{path}/output")))





# output RERI plots (can't install package on CARC yet)
significant_snps <- c("6:12577203:T:C", "5:40252294:C:T") 
significant_snps <- c("5:40252294:C:T") # not sure what's going on here.. do i need to flip the aspirin coding? -- answer is YES

walk(significant_snps, ~ reri_wrapper(data_epi = input_data, exposure = exposure, snp = .x, covariates = covariates, path = glue("{path}/output")))

input_data_flip <- input_data %>% 
  mutate(aspirin = fct_relevel(aspirin, "Yes"))

walk(significant_snps, ~ reri_wrapper(data_epi = input_data_flip, exposure = exposure, snp = .x, covariates = covariates, path = glue("{path}/output")))


# ================================================================== #
# ======= rmarkdown reports ---- 
# ================================================================== #
exposure = 'aspirin'
hrc_version = 'v3.0'
path = glue("/media/work/gwis_test/{exposure}")

gwis_report(exposure = exposure, 
            hrc_version = hrc_version, 
            covariates = covariates)

posthoc_report(exposure = exposure)














