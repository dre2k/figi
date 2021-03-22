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
exposure = 'fruitqc2'
hrc_version = 'v2.3'
output_dir = paste0("/media/work/gwis/posthoc/", exposure, "/")
annotation_file <- 'gwas_141_ld_annotation_july2020.txt'

covariates <- sort(c('age_ref_imp', 'sex', 'energytot_imp', 'study_gxe', 'pc1', 'pc2', 'pc3'))
covariates_set1 <- covariates
# covariates_set2 <- sort(c(covariates_set1, 'bmi', 'smk_ever', 'fruitqc2', 'vegetableqc2'))
covariates_list <- list(covariates_set1)
mod <- 'age_ref_imp+sex+energytot_imp+pc1+pc2+pc3+study_gxe'

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
covariates_meta <- sort(c('age_ref_imp', 'energytot_imp','sex'))
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

# make sure 1-28g/d is the reference group !
figi <- posthoc_input(exposure, hrc_version, glue('gwis_sig_results_output_{exposure}.rds'))

snps <- readRDS(glue(output_dir, "gwis_sig_results_input_{exposure}.rds"))
table(snps$method)


# ------- GxE 1DF suggestive findings ------- #
# in this case, i'll use significant findings since we have two SNPs.. 
snps_filter <- snps %>% 
  dplyr::filter(SNP %in% c("14:74029409:C:T", "14:74029049:G:C")) %>% 
  dplyr::mutate(snps = paste0('chr', gsub("\\:", "\\_", SNP))) %>% 
  pull(snps)

posthoc_run_models(snps_filter, 'chiSqGxE')
posthoc_run_reri_plot(snps_filter)
posthoc_create_plots(snps_filter, 'chiSqGxE')



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





# ================================================================== #
# ======= rmarkdown reports ---- 
# ================================================================== #
gwis_report(exposure = exposure, 
            hrc_version = hrc_version, 
            covariates = covariates)

posthoc_report(exposure = exposure)










#-------------------------------------------- #
# output all statistics from 3df findings to answer nick email
# -------------------------------------------- #

snps_3df <- readRDS(glue(output_dir, "manhattan_chiSq3df_{exposure}_no_gwas_clump_df.rds"))

write.csv(snps_3df, file = "~/Dropbox/fruitqc2_chiSq3df_all_stats.csv", quote = F, row.names = F)



# investigate this one - 2:136608646:G:A

threedf_clump <- fread(glue("/media/work/gwis/clump_combined/FIGI_{hrc_version}_gxeset_{exposure}_chiSq3df_no_gwas_ldclump.clumped"), data.table = F)


test <- threedf_clump %>% 
  filter(grepl("136608646", SP2)) 

test <- gxe %>% 
  filter(SNP == "2:136608646:G:A")
write.csv(test, file = "~/Dropbox/fruitqc2_tmp.csv", quote = F, row.names = F)

# ================================================================== #
# ======= binarydosage output of main findings (all methods) ======= #
figi <- posthoc_input(exposure, hrc_version)
# ================================================================== #

posthoc_run_models <- function(x, method) {
  # pvalues dataframe
  # important to specify method to calculate proper p-value
  walk(x, ~ pval_summary(figi, exposure, .x, covariates_list, method, output_dir = output_dir))
  
  # model estimates (stargazer html)
  # likewise, important to specify method to calculate proper p-value
  walk(x, ~ fit_gxe_covars(figi, exposure, .x, covariates_list, method, output_dir = output_dir))
  
  # stratified odds ratio
  # walk(x, ~ fit_stratified_or_q4(figi, exposure, snp = .x, hrc_version = hrc_version, covariates = covariates_set1, mod = mod, dosage = T, output_dir = output_dir))
}


posthoc_create_plots <- function(x, statistic_to_plot) {
  # locuszoom and functional annotation plots
  # need to format SNP properly for locuszoom function call 
  locuszoom_helper <- function(snp) {
    tmp <- gsub("chr", "", snp) # have to turn snps into chr:bp format for locuszoom 
    snp_chr <- strsplit(tmp, "_")[[1]][1]
    snp_bp <- as.numeric(strsplit(tmp, '_')[[1]][2])
    paste(snp_chr, snp_bp, sep = ":")
  }
  
  snps_plot <- map_chr(x, locuszoom_helper)
  walk(snps_plot, ~ system(paste("bash ~/Dropbox/FIGI/FIGI_code/results/posthoc/posthoc_02_locuszoom.sh", exposure, hrc_version, .x, statistic_to_plot)))
  # statistic_to_plot for functional annotation is to help identify the locuszoom plot output directory (covenience)
  walk(snps_plot, ~ system(paste("Rscript ~/Dropbox/FIGI/FIGI_code/results/posthoc/posthoc_03_functional_annotation.R", exposure, .x, statistic_to_plot)))
}



# ------ Suggestive GxE hits (1df) ------ #
snps <- readRDS(paste0("/media/work/gwis/posthoc/gwis_sig_results_input_", exposure, ".rds")) %>%
  filter(method == paste0("chiSqGxE_", exposure, "_clump")) %>% 
  arrange(Chromosome, Location) %>% 
  pull(SNP)
snps <- paste0('chr', gsub("\\:", "\\_", snps))

posthoc_run_models(snps, 'chiSqGxE')
posthoc_create_plots(snps, 'chiSqGxE')


# ------ 



# ------ 3DF findings ------ #
keep = c("1:177889480:A:G", "2:558048:T:G",  "4:45181334:A:T", "10:114779013:C:T", "12:112883476:G:A", "16:28649651:C:A", 	"18:57968666:C:T")

# all findings (top hit only)
snps <- readRDS(paste0("/media/work/gwis/posthoc/gwis_sig_results_input_", exposure, ".rds")) %>%
  filter(method == "chiSq3df_bmi5_no_gwas_no_marginal_no_ge") %>% 
  arrange(Chromosome, Location) %>% 
  filter(SNP %in% keep) %>% 
  pull(SNP) %>% 
  paste0('chr', .) %>% 
  gsub("\\:", "\\_", .)

posthoc_run_models(snps, 'chiSq3df')
posthoc_create_plots(snps, 'chiSq3df')





# from aki email request (08/21/2020)
keep <- c(28561982,28565667,28566158,28649651) # chromosome 16


snps <- readRDS(paste0("/media/work/gwis/posthoc/gwis_sig_results_input_", exposure, ".rds")) %>% 
  filter(grepl("chiSq3df",method)) %>% 
  arrange(Chromosome, Location) %>% 
  filter(Chromosome == 16 & Location %in% keep)

snps <- paste0('chr', gsub("\\:", "\\_", snps$SNP))
posthoc_run_models(snps, 'chiSq3df')





