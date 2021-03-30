#=============================================================================#
# FIGI GxE fruitqc2 results
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
exposure = 'fruitqc2'
hrc_version = 'v2.3'
annotation_file <- 'gwas_200_ld_annotation_feb2021.txt'
path = glue("/media/work/gwis_test/{exposure}/")

covariates <- sort(c('age_ref_imp', 'sex', 'energytot_imp', 'study_gxe', 'pc1', 'pc2', 'pc3'))

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





# ----- additional ---- #
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

source(glue("/media/work/gwis_test/R/03_posthoc.R"))
source(glue("/media/work/gwis_test/R/02_plots.R"))

# SNPs to generate items
snps <- c("14:74029409:C:T", "1:72729142:A:G")
snps <- c("14:74029409:C:T")
snps <- c("1:72729142:A:G")

# output GxE models adjusted by different covariate sets
covariates_sets <- list(covariates, 
                        c(covariates, 'bmi5'), 
                        c(covariates, 'bmi5', 'smk_ever'), 
                        c(covariates, 'bmi5', 'smk_ever', 'fruitqc2', 'vegetableqc2'))

covariates_sets <- list(covariates)

walk(snps, ~ fit_gxe_covars(data_epi = input_data, exposure = exposure, snp = .x, covariates_list = covariates_sets, method = 'chiSqGxE', path = glue("{path}/output")))


# additional covariates is making association more significant.. let's generate for all suggestive hits
suggestive_gxe <- fread(glue("{path}/data/FIGI_v2.3_gxeset_diab_chiSqGxE_ldclump.clumped"))
walk(suggestive_gxe$SNP, ~ fit_gxe_covars(data_epi = input_data, exposure = exposure, snp = .x, covariates_list = covariates_sets, method = 'chiSqGxE', path = glue("{path}/output")))


# output RERI plots (can't install package on CARC yet)
walk(snps, ~ reri_wrapper(data_epi = input_data, exposure = exposure, snp = .x, covariates = covariates, path = glue("{path}/output")))




# output easystrata manhattan plot for Nikolas that excludes the right set of annotations (200 GWAS instead of 140)

gxe <- readRDS("/media/work/gwis/results/fruitqc2/processed/FIGI_v2.3_gxeset_fruitqc2_basic_covars_gxescan_results.rds")

create_manhattanplot(x = gxe, exposure = exposure, stat = 'chiSq3df', annotation_file = annotation_file, output_dir = "~/Dropbox")


gwas <- fread("/home/rak/data/Annotations/gwas_200_ld_annotation_feb2021.txt") %>% 
  mutate(SNP2 = paste0(Chr, ":", Pos))
gxe_nogwas <- dplyr::filter(gxe, !SNP2 %in% gwas$SNP2)

create_manhattanplot(x = gxe, exposure = exposure, stat = 'chiSq3df', annotation_file = annotation_file, output_dir = "~/Dropbox")

create_manhattanplot(x = gxe_nogwas, exposure = exposure, stat = 'chiSq3df', annotation_file = annotation_file, output_dir = "~/Dropbox")



# try again - remove flanking regions for better visualizaation
# 
# 


create_manhattanplot_ramwas(data = gxe, exposure = exposure, statistic = 'chiSqGxE', hrc_version = 'v2.3', path = "~/Dropbox/", sig_line = 5e-8)



gxe_better <- gxe_nogwas %>% 
  filter(chiSqG_p > 5e-7)


create_manhattanplot_ramwas(data = gxe, exposure = exposure, statistic = 'chiSqGxE', hrc_version = 'v2.3', path = "~/Dropbox/", sig_line = 5e-8)
create_manhattanplot_ramwas(data = gxe_better, exposure = exposure, statistic = 'chiSq3df', hrc_version = 'v2.3', path = "~/Dropbox/", sig_line = 5e-8)
create_manhattanplot_ramwas(data = gxe_better, exposure = exposure, statistic = 'chiSq2df', hrc_version = 'v2.3', path = "~/Dropbox/", sig_line = 5e-8)






# ================================================================== #
# ======= rmarkdown reports ---- 
# ================================================================== #

gwis_report(exposure = exposure, hrc_version = hrc_version, covariates = covariates)

posthoc_report(exposure = exposure)


