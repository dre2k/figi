#=============================================================================#
# FIGI GxE folate_totqc2 results
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
exposure = 'folate_diet400qcm'
hrc_version = 'v3.0'
annotation_file <- 'gwas_200_ld_annotation_feb2021.txt'
path = glue("/media/work/gwis_test/{exposure}/")
# path = glue("/media/work/gwis_test/{exposure}_v23/")

covariates <- sort(c('age_ref_imp', 'sex', 'energytot_imp', 'study_gxe', 'pc1', 'pc2', 'pc3'))

# input data
esubset <- readRDS(glue("{path}/data/FIGI_{hrc_version}_gxeset_{exposure}_basic_covars_glm.rds")) %>% pull(vcfid)
input_data <- readRDS(glue("/media/work/gwis/data/FIGI_EpiData/FIGI_{hrc_version}_gxeset_analysis_data_glm.rds")) %>% 
  filter(vcfid%in% esubset)


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
# SNP followup ---- 
#-----------------------------------------------------------------------------#
# compile significant results into data.frame
source("/home/rak/Dropbox/FIGI/FIGI_code/results/posthoc/posthoc_01_combine_results.R")

# on HPC:
# - extract dosage information on HPC
# - calculate clumped statistics - gxe, 2/3df if necessary








# ================================================================== #
# ======= rmarkdown reports ---- 
# ================================================================== #
main_effects_report <- function (exposure, hrc_version, covariates, path) {
  rmarkdown::render("~/git/figi/main_effects/main_effects.Rmd", 
                    params = list(exposure = exposure, hrc_version = hrc_version, 
                                  covariates = covariates, path = path), output_file = glue("~/Dropbox/FIGI/Results/{exposure}_{hrc_version}_main_effects.html"))
}

main_effects_report(exposure = exposure, hrc_version = hrc_version, covariates = covariates, path = path)

gwis_report <- function (exposure, hrc_version, covariates) {
  rmarkdown::render("~/git/figi/gwis/results.Rmd", 
                    params = list(exposure = exposure, 
                                  hrc_version = hrc_version, 
                                  covariates = covariates),
                    output_file = glue("~/Dropbox/FIGI/Results/{exposure}_{hrc_version}_gwis.html"))
}

gwis_report(exposure = exposure, 
            hrc_version = hrc_version, 
            covariates = covariates)


posthoc_report(exposure = exposure, 
               hrc_version = hrc_version,
               covariates = covariates,
               path = path)










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




# leave one out type analysis (studies, etc)

#different covariate adjustments. 




#----------------------------------------------------------------------------#
# conditional analyses ----


# conditional analysis and locuszoom plots
gwas_annotation <- fread("~/data/Annotations/gwas_141_ld_annotation_july2020.txt")

gecco_gwas <- read_tsv("/media/work/gwis/gwas_hits_140/KnownLoci/crc_gwas_indep_signals_052920.tsv") %>% 
  mutate(SNP = paste0("chr", CHROM, ":", POSITION))

# file to annotate locuszoom plots with text
# tmp <- gecco_gwas %>% 
#   mutate(snp = paste0("chr", CHROM, ":", POSITION), 
#          string = paste0(AUTHOR_FIRST_REPORTED, "_", YEAR_FIRST_REPORTED),
#          color = "black") %>% 
#   dplyr::select(snp, string, color)
# write.table(tmp, file = "/media/work/gwis/locuszoom/folate_dietqc2/folate_dietqc2_locuszoom_gwas_text_annotation.txt", sep = "\t", quote = F, row.names = F)


# first, get signifiicant hits for 2DF analysis, I want to properly create the stuff


# 2DF significant results, after removing GWAS SNPs
folate_2df_sig <- readRDS("/media/work/gwis/posthoc/folate_dietqc2/significant_results_dataframe_chiSq2df_folate_dietqc2_no_gwas.rds")


# conditional analyses? 


# look at chromosome 15 first
tmp <- filter(gecco_gwas, CHROM == 15)

gwas_140 <- readRDS("/media/work/gwis/FIGI_genotype_dosages_gwas140.rds")
folate_2df <- readRDS("/media/work/gwis/posthoc/gwis_sig_results_output_folate_dietqc2.rds")

xxx <- inner_join(gwas_140, folate_2df, 'vcfid') %>% 
  inner_join(xx, 'vcfid') %>% 
  rename(chr15_33112299_folate = X15.33112299.C.T_dose, 
         chr15_32992836_tomlinson_2011 = X15.32992836.G.A, 
         chr15_33156386_huygue_2019 = X15.33156386.G.A, 
         chr15_33010736_tomlinson_2011 = X15.33010736.G.A)

cor_matrix <- as.data.frame(cor(xxx[, c('chr15_33112299_folate', 'chr15_33010736_tomlinson_2011', 'chr15_32992836_tomlinson_2011', 'chr15_33156386_huygue_2019')]))


model_chr15_33112299 <- glm(outcome ~ folate_dietqc2 + chr15_33112299_folate + age_ref_imp + sex + energytot_imp + study_gxe + pc1 + pc2 + pc3, data = xxx, family = 'binomial')

model_chr15_32992836 <- glm(outcome ~ folate_dietqc2 + chr15_32992836_tomlinson_2011 + age_ref_imp + sex + energytot_imp + study_gxe + pc1 + pc2 + pc3, data = xxx, family = 'binomial')

model_chr15_33156386 <- glm(outcome ~ folate_dietqc2 + chr15_33156386_huygue_2019 + age_ref_imp + sex + energytot_imp + study_gxe + pc1 + pc2 + pc3, data = xxx, family = 'binomial')

model_chr15_33010736 <- glm(outcome ~ folate_dietqc2 + chr15_33010736_tomlinson_2011 + age_ref_imp + sex + energytot_imp + study_gxe + pc1 + pc2 + pc3, data = xxx, family = 'binomial')

# model individual SNPs
tab_model(model_chr15_33112299, model_chr15_32992836, model_chr15_33156386, model_chr15_33010736, 
          terms = c("(Intercept)", "folate_dietqc2", "age_ref_imp", "sex", "chr15_33112299_folate", "chr15_32992836_tomlinson_2011", "chr15_33156386_huygue_2019", "chr15_33010736_tomlinson_2011"), 
          p.style = "scientific", title = "Marginal associations for chr15 region findings")


# pairwise conditional analyses
aa <- glm(outcome ~ folate_dietqc2 + chr15_33112299_folate + chr15_32992836_tomlinson_2011 + age_ref_imp + sex + energytot_imp + study_gxe + pc1 + pc2 + pc3, data = xxx, family = 'binomial')
bb <- glm(outcome ~ folate_dietqc2 + chr15_33112299_folate + chr15_33156386_huygue_2019 + age_ref_imp + sex + energytot_imp + study_gxe + pc1 + pc2 + pc3, data = xxx, family = 'binomial')
cc <- glm(outcome ~ folate_dietqc2 + chr15_33112299_folate + chr15_33010736_tomlinson_2011 + age_ref_imp + sex + energytot_imp + study_gxe + pc1 + pc2 + pc3, data = xxx, family = 'binomial')
dd <- glm(outcome ~ folate_dietqc2 + chr15_33112299_folate + chr15_32992836_tomlinson_2011 +  chr15_33156386_huygue_2019 + chr15_33010736_tomlinson_2011 + age_ref_imp + sex + energytot_imp + study_gxe + pc1 + pc2 + pc3, data = xxx, family = 'binomial')

tab_model(aa, bb, cc, dd, 
          terms = c("(Intercept)", "folate_dietqc2", "age_ref_imp", "sex", "chr15_33112299_folate", "chr15_32992836_tomlinson_2011", "chr15_33156386_huygue_2019", "chr15_33010736_tomlinson_2011"), 
          p.style = "scientific", title = "conditional analysis chr15_33112299_folate")



# hit doesn't survive conditional analysis
model_chr15_33112299_cond_chr15_33010736 <- glm(outcome ~ folate_dietqc2 + X15.33112299.C.T_dose +X15.33010736.G.A + age_ref_imp + sex + energytot_imp + study_gxe + pc1 + pc2 + pc3, data = xx, family = 'binomial')
summary(model_chr15_33112299_cond_chr15_33010736)

cor(xx$X15.33112299.C.T_dose, xx$X15.33010736.G.A)


tab_model(model_chr15_33112299, model_chr15_33112299_cond_chr15_33010736, terms = c("(Intercept)", "folate_dietqc2", "age_ref_imp", "sex", "X15.33112299.C.T_dose", "X15.33010736.G.A"), p.style = "scientific", title = "conditional analysis chr15:33112299")

# chromosome 18

xxx <- inner_join(gwas_140, folate_2df, 'vcfid') %>% 
  inner_join(xx, 'vcfid') %>% 
  rename(chr18_46371993 = X18.46371993.C.T_dose, 
         chr18_46460385 = X18.46460385.A.G_dose, 
         chr18_46453156_broderick_2007 = X18.46453156.A.T)


tmp <- filter(gecco_gwas, CHROM == 18)

model_chr18_46371993 <- glm(outcome ~ folate_dietqc2 + chr18_46371993 + age_ref_imp + sex + energytot_imp + study_gxe + pc1 + pc2 + pc3, data = xxx, family = 'binomial')

model_chr18_46460385 <- glm(outcome ~ folate_dietqc2 + chr18_46460385 + age_ref_imp + sex + energytot_imp + study_gxe + pc1 + pc2 + pc3, data = xxx, family = 'binomial')

model_chr18_46453156_broderick_2007 <- glm(outcome ~ folate_dietqc2 + chr18_46453156_broderick_2007 + age_ref_imp + sex + energytot_imp + study_gxe + pc1 + pc2 + pc3, data = xxx, family = 'binomial')

model_chr18_all <- glm(outcome ~ folate_dietqc2 + chr18_46371993 + chr18_46460385 + chr18_46453156_broderick_2007 + age_ref_imp + sex + energytot_imp + study_gxe + pc1 + pc2 + pc3, data = xxx, family = 'binomial')

tab_model(model_chr18_46371993, model_chr18_46460385, model_chr18_46453156_broderick_2007, model_chr18_all, terms = c("(Intercept)", "folate_dietqc2", "age_ref_imp", "sex", "chr18_46371993", "chr18_46460385", "chr18_46453156_broderick_2007"), p.style = "scientific", title = "conditional analysis chr18 findings")




## first hit
model_chr18_46371993_cond_chr18_46453156 <- glm(outcome ~ folate_dietqc2 + X18.46371993.C.T_dose + X18.46453156.A.T + age_ref_imp + sex + energytot_imp + study_gxe + pc1 + pc2 + pc3, data = xx, family = 'binomial')
summary(model_chr18_46371993_cond_chr18_46453156)

tab_model(model_chr18_46371993, model_chr18_46371993_cond_chr18_46453156, terms = c("(Intercept)", "folate_dietqc2", "age_ref_imp", "sex", "X18.46371993.C.T_dose", "X18.46453156.A.T"), p.style = "scientific", title = "conditional analysis chr18:46371993")


















