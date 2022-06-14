#=============================================================================#
# FIGI GxE asp_ref results
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
exposure = 'bmi5'
hrc_version = 'v2.3'
annotation_file <- 'gwas_200_ld_annotation_feb2021.txt'
covariates <- sort(c('age_ref_imp', 'sex', 'study_gxe', 'pc1', 'pc2', 'pc3'))
path = glue("/media/work/gwis_test/{exposure}/")


# input data
esubset <- readRDS(glue("/media/work/gwis_test/{exposure}/data/FIGI_{hrc_version}_gxeset_{exposure}_basic_covars_glm.rds")) %>% 
  pull(vcfid)

input_data <- readRDS(glue("/media/work/gwis_test/data/FIGI_{hrc_version}_gxeset_analysis_data_glm.rds")) %>% 
  filter(vcfid%in% esubset) %>% 
  mutate(bmic4f = case_when(bmi < 18.5 ~ "Underweight", 
                            between(bmi, 18.5, 24.99) ~ "Normal",
                            between(bmi, 25, 29.99) ~ "Overweight", 
                            bmi > 29.99 ~ "Obese"),
         bmic3f = fct_relevel(factor(bmic3), "Normal", "Overweight", "Obese"))


# mean / SD by study_gxe (bmi)
unique(input_data$study_gxe)

out <- input_data %>% 
  group_by(study_gxe, outcome) %>% 
  summarise(bmi_mean = mean(bmi),
            bmi_sd = sd(bmi))
write.csv(out, file = "~/Dropbox/Working/figi_bmi_mean_sd.csv")




#-----------------------------------------------------------------------------#
# main effects ----
#-----------------------------------------------------------------------------#

# for BMI5, 0 in model is essentially thelower end of health BMI, so it's ok to keep meta-analyses as is

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
results_gxe <- c("11:47328213:T:C", "3:100037644:A:G", "2:10961902:T:C", "11:120531547:C:G","11:46574459:T:A", "2:217410559:G:A", "6:138766021:C:T", "9:136601110:C:G", "15:33122966:C:T")
results_twostep <- c("15:33122966:C:T", "15:33120215:T:C")
results_3df <- c("1:177889480:A:G", "12:50204089:G:A", "16:28649651:C:A", "16:53800954:T:C","18:57850583:C:A", "18:57971625:A:G", "18:58083923:A:C", "19:47581242:C:T","2:41716#7:T:C", "2:558080:G:T", "2:651365:G:A", "4:45181334:A:T")
# snps <- c("20:6584196:C:G", "15:33122966:C:T", "15:32987718:A:G", "5:134486618:G:A", "10:114274269:T:C", "12:50177324:A:T")

snps <- c(results_gxe, results_twostep)
snps <- c("12:50204089:G:A")


# ---- MAF ---- #
walk(snps, ~ create_aaf_study_plot(data = input_data, exposure = exposure, hrc_version = hrc_version, snp = .x, path = path))

walk(results_3df, ~ create_aaf_study_plot(data = input_data, exposure = exposure, hrc_version = hrc_version, snp = .x, path = path))



# ---- stratified odds ratios ---- #
# stratified odds ratio table with and without smoking
walk(snps, ~ fit_stratified_or_q3(
  data_epi = input_data,
  exposure = 'bmic3f',
  snp = .x,
  hrc_version = hrc_version,
  covariates = covariates,
  path = glue("{path}/output/")))


covariates_smk <- c(covariates, 'smk_ever')
walk(snps, ~ fit_stratified_or_q3(
  data_epi = input_data,
  exposure = 'bmic3f',
  snp = .x,
  hrc_version = hrc_version,
  covariates = covariates_smk,
  path = glue("{path}/output/")))



# ---- iplot ---- #
walk(snps, ~iplot_wrapper(data_epi = input_data, exposure = 'bmic3f', hrc_version = hrc_version, snp = .x, covariates = covariates, path = glue("{path}/output/"), flip_allele = F))




# output GxE models adjusted by different covariate sets
covariates_sets <- list(covariates, 
                        c(covariates, "smk_ever"))
walk(snps, ~ fit_gxe_covars(data_epi = input_data, exposure = exposure, snp = .x, covariates_list = covariates_sets, method = 'chiSqGxE', path = glue("{path}/output")))




# stratified by tumor and sex 

walk(snps, ~ fit_gxe_stratified(data_epi = input_data, exposure = exposure, snp = .x, covariates = covariates, method = "chiSqGxE", strata = 'sex', path = glue("{path}/output")))
walk(snps, ~ fit_gxe_stratified(data_epi = input_data, exposure = exposure, snp = .x, covariates = covariates, method = "chiSqGxE", strata = 'cancer_site_sum2', path = glue("{path}/output")))
walk(snps, ~ fit_gxe_stratified(data_epi = input_data, exposure = exposure, snp = .x, covariates = covariates, method = "chiSqGxE", strata = 'study_design', path = glue("{path}/output")))




# output RERI plots
snps <- c("15:33122966:C:T", "15:33120215:T:C")
walk(snps, ~ reri_wrapper(data_epi = input_data, exposure = 'bmic3f', snp = .x, covariates = covariates, path = glue("{path}/output")))



# SNP information
snp_info <- qread("/media/work/FIGI_RsqEstimate_chrALL.qs") %>% 
  filter(id %in%  c("15:33122966:C:T"))

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









# "fine mapping" ------------------------------------------------------------


# basically - g







# --------------- #
# ----- old ----- #
# --------------- #

# stratified OR
walk(snps, ~ fit_stratified_or(
  data_epi = input_data,
  exposure = exposure,
  snp = .x,
  hrc_version = hrc_version,
  covariates = covariates,
  path = glue("{path}/output/")))

snps <- c("15:33122966:C:T", "15:33120215:T:C")
walk(snps, ~ fit_stratified_or_continuous(
  data_epi = input_data,
  exposure = exposure,
  snp = .x,
  hrc_version = hrc_version,
  covariates = covariates,
  path = glue("{path}/output/")))

# stratified odds ratio table with and without smoking
snps <- c( "15:33122966:C:T", "15:33120215:T:C")
walk(snps, ~ fit_stratified_or(
  data_epi = input_data,
  exposure = exposure,
  snp = .x,
  hrc_version = hrc_version,
  covariates = covariates,
  path = glue("{path}/output/")))


covariates_smk <- c(covariates, 'smk_ever')
walk(snps, ~ fit_stratified_or(
  data_epi = input_data,
  exposure = exposure,
  snp = .x,
  hrc_version = hrc_version,
  covariates = covariates_smk,
  path = glue("{path}/output/")))


# ================================================================== #
# ======= rmarkdown reports ---- 
# ================================================================== #
main_effects_report(exposure = exposure, hrc_version = hrc_version, covariates = covariates, path = path)


gwis_report(exposure = exposure, 
            hrc_version = hrc_version, 
            covariates = covariates)

posthoc_report(exposure = exposure)

posthoc_report(exposure = exposure, 
               hrc_version = hrc_version,
               covariates = covariates,
               path = path)

posthoc_report(exposure = exposure)




# test
rmarkdown::render(glue("/home/rak/git/figi/{exposure}_posthoc_test.Rmd"), 
                  output_file = glue("~/Dropbox/FIGI/Results/{exposure}_posthoc_test.html"))


