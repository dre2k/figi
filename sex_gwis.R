#=============================================================================#
# FIGI GxE aspirin results
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
exposure = 'sex'
hrc_version = 'v2.3'
annotation_file <- 'gwas_200_ld_annotation_feb2021.txt'
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

create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "all", forest_height = 15, forest_width = 9, categorical = T)
create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "proximal", forest_height = 13, forest_width = 9, categorical = T)
create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "distal", forest_height = 13, forest_width = 9, categorical = T)
create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "rectal", forest_height = 13, forest_width = 9, categorical = T)

# ------- stratified pooled analysis ------- #
pooled_analysis_glm(input_data, exposure = exposure, hrc_version = hrc_version, covariates = covariates, strata = 'study_design', filename_suffix = paste0(paste0(sort(covariates), collapse = '_'), "_stratified_study_design"), output_dir = glue("{path}/output/posthoc/"))

pooled_analysis_multinom(input_data, exposure = exposure, hrc_version = hrc_version, covariates = covariates, filename_suffix = paste0(paste0(sort(covariates), collapse = '_'), "_stratified_cancer_site_sum2"), output_dir = glue("{path}/output/posthoc/"))

# pooled analyses stratified by BMI level 
pooled_analysis_glm(input_data, exposure = exposure, hrc_version = hrc_version, covariates = covariates, strata = 'bmic3', filename_suffix = paste0(paste0(sort(covariates), collapse = '_'), "_stratified_bmic3"), output_dir = glue("{path}/output/posthoc/"))



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
# 
# plot(inf.analysis, "influence")
# plot(inf.analysis, "baujat")
# plot(inf.analysis, "es")
# plot(inf.analysis, "i2")






#-----------------------------------------------------------------------------#
# GxE additional analysis ---- 
#-----------------------------------------------------------------------------#
# ---- would this be significant with this subset ---- #
dos <- qread("/media/work/gwis_test/sex/output/posthoc/dosage_chr18_46458950.qs")
out <- inner_join(input_data, dos, 'vcfid')

summary(glm(outcome ~ chr18_46458950_G_A_dose + sex + age_ref_imp + study_gxe + pc1 + pc2 + pc3, data = out, family = 'binomial'))
summary(glm(outcome ~ chr18_46458950_G_A_dose * sex + age_ref_imp + study_gxe + pc1 + pc2 + pc3, data = out, family = 'binomial'))

wtf <- readRDS("/media/work/gwis_test/sex/output/twostep_gao_0.995_sex_v2.3_chiSqG_df_eric.rds")


dos <- qread("/media/work/gwis_test/sex/output/posthoc/dosage_chr18_46458950.qs")
nat <- readRDS("/media/work/gwis_test/sex/data/cov_by_study_gxe_5pc.rds") %>% 
  pull(vcfid)
out <-  %>% 
  filter(vcfid %in% nat)


summary(glm(outcome ~ chr18_46458950_G_A_dose + sex + age_ref_imp + study_gxe + pc1 + pc2 + pc3, data = out, family = 'binomial'))
summary(glm(outcome ~ chr18_46458950_G_A_dose * sex + age_ref_imp + study_gxe + pc1 + pc2 + pc3, data = out, family = 'binomial'))




# ---- MAF ---- #
snps <- c("18:46458950:G:A", "12:116201277:T:A")
walk(snps, ~ create_aaf_study_plot(data = input_data, exposure = exposure, hrc_version = hrc_version, snp = .x, path = path))


# ---- stratified odds ratios ---- #
input_data2 <- input_data %>% 
  mutate(sex =)
walk(snps, ~ fit_stratified_or(
  data_epi = input_data2,
  exposure = exposure,
  snp = .x,
  hrc_version = hrc_version,
  covariates = covariates,
  path = glue("{path}/output/")))




# output GxE models adjusted by different covariate sets
# covariates_sets <- list(covariates, 
#                         c(covariates, 'bmi5'), 
#                         c(covariates, 'bmi5', 'smk_ever'), 
#                         c(covariates, 'bmi5', 'smk_ever', 'fruitqc2', 'vegetableqc2'))
# 
# fit_gxe_covars(data_epi = input_data, exposure = exposure, snp = "13:47191972:G:A", covariates_list = covariates_sets, method = 'chiSqGxE', path = glue("{path}/output"))
# 
# # additional covariates is making association more significant.. let's generate for all suggestive hits
# suggestive_gxe <- fread(glue("{path}/data/FIGI_{hrc_version}_gxeset_{exposure}_chiSqGxE_ldclump.clumped"))
# walk(suggestive_gxe$SNP, ~ fit_gxe_covars(data_epi = input_data, exposure = exposure, snp = .x, covariates_list = covariates_sets, method = 'chiSqGxE', path = glue("{path}/output")))




# output GxE models adjusted by different covariate sets
covariates_sets <- list(covariates,
                        c(covariates, 'bmi5'),
                        c(covariates, 'bmi5', 'smk_ever'),
                        c(covariates, 'bmi5', 'smk_ever', 'fruitqc2', 'vegetableqc2', 'fiberqc2'))
covariates_sets <- list(covariates)

fit_gxe_covars(data_epi = input_data, exposure = exposure, snp = "5:40252294:C:T", covariates_list = covariates_sets, method = 'chiSqGxE', path = glue("{path}/output"))
fit_gxe_covars(data_epi = input_data, exposure = exposure, snp = "6:12577203:T:C", covariates_list = covariates_sets, method = 'chiSqGxE', path = glue("{path}/output"))


# stratified by tumor and sex 
walk(snps, ~ fit_gxe_stratified(data_epi = input_data, exposure = exposure, snp = .x, covariates = covariates, method = "chiSqGxE", strata = 'cancer_site_sum2', path = glue("{path}")))



# output RERI plots (can't install package on CARC yet)
significant_snps <- c("6:12577203:T:C")
walk(significant_snps, ~ reri_wrapper(data_epi = input_data, exposure = exposure, snp = .x, covariates = covariates, path = glue("{path}/output")))

significant_snps <- c("5:40252294:C:T")
input_data_flip <- input_data %>% 
  mutate(aspirin = fct_relevel(aspirin, "Yes"))

walk(significant_snps, ~ reri_wrapper(data_epi = input_data_flip, exposure = exposure, snp = .x, covariates = covariates, path = glue("{path}/output")))






# SNP information
# maybe you should get info for all LD SNPs as well
chr5 <- fread(glue("/media/work/gwis_test/aspirin/output/functional_plot/functional_annotation_chr5_40252294.bed")) %>% 
  mutate(SNP = paste(gsub("chr", "", V1), V2, V3, sep = ":"))
chr6_125 <- fread(glue("/media/work/gwis_test/aspirin/output/functional_plot/functional_annotation_chr6_12577203.bed")) %>% 
  mutate(SNP = paste(gsub("chr", "", V1), V2, V3, sep = ":"))


snp_info <- qread("/media/work/FIGI_RsqEstimate_chrALL.qs") %>% 
  filter(id %in%  c("6:12577203:T:C", "5:40252294:C:T"))

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








# ================================================================== #
# Rmarkdown reports -------------------------------------------------------
# ================================================================== #

main_effects_report(exposure = exposure, hrc_version = hrc_version, covariates = covariates, path = path)

gwis_report(exposure = exposure, hrc_version = hrc_version, covariates = covariates)


rmarkdown::render(glue("/home/rak/git/figi/{exposure}_posthoc.Rmd"), 
                  params = list(exposure = exposure), 
                  output_file = glue("~/Dropbox/FIGI/Results/{exposure}_posthoc.html"))



