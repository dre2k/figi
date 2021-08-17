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








# ----------------------------------------------------------------------------#
# simple Meff, important ---- 

# try to run simpleM on chromosome (based on functional annotation overlap....)

# ------- run simpleM ------ #
PCA_cutoff <- 0.995

# important variables, hardcoded in the simple M wrapper function 
number_of_snps <- c()
number_of_tests <- c()
simplem_step1_statistic = 'chiSqGxE'


for(chr in 1:22) {
  
  mySNP_nonmissing <- readRDS(paste0("/media/work/gwis/twostep_expectation_hybrid/", exposure, "/", exposure, "_snplist_functional_annotation_test_", simplem_step1_statistic, "_chr", chr, "_allcontrols_output.rds"))
  # mySNP_nonmissing <- readRDS("/media/work/gwis/twostep_expectation_hybrid/folate_dietqc2/folate_dietqc2_snplist_functional_annotation_test_chiSqGxE_chr22_allcontrols_output.rds")
  
  workdata <- mySNP_nonmissing %>% 
    pivot_longer(-vcfid, ) %>% 
    pivot_wider(names_from = vcfid, values_from = value) %>% 
    separate(name, into = c("chr", "bp", "ref", "alt"), remove = F) %>% 
    mutate(chr = as.numeric(gsub("X", "", chr)), 
           bp = as.numeric(bp)) %>% 
    arrange(chr, bp) %>% 
    dplyr::select(-chr, -bp, -ref, -alt, -name)
  
  numLoci <- length(pull(workdata, 1))
  
  simpleMeff <- NULL
  
  fixLength <- 150
  i <- 1
  myStart <- 1
  myStop <- 1
  iteration <- 0
  
  while(myStop < numLoci){
    myDiff <- numLoci - myStop 
    if(myDiff <= fixLength) break
    
    myStop <- myStart + i*fixLength - 1
    snpInBlk <- t(workdata[myStart:myStop, ])
    MeffBlk <- inferCutoff(snpInBlk)
    simpleMeff <- c(simpleMeff, MeffBlk)
    myStart <- myStop+1
    iteration <- iteration+1
    print(iteration)
  }
  
  snpInBlk <- t(workdata[myStart:numLoci, ])
  MeffBlk <- inferCutoff(snpInBlk)
  simpleMeff <- c(simpleMeff, MeffBlk)
  
  cat("Total number of SNPs is: ", numLoci, "\n")
  cat("Inferred Meff is: ", sum(simpleMeff), "\n")
  
  # append results to vectors
  number_of_snps <- c(number_of_snps, ncol(mySNP_nonmissing) - 1)
  number_of_tests <- c(number_of_tests, sum(simpleMeff))
}


# nbins always 18 with one million tests expectation
nbins = 18
alpha = 0.05
out <- format_data_twostep_expectation_hybrid(dat = gxescan_output, stats_step1 = simplem_step1_statistic, sizeBin0 = 5, alpha = 0.05)

# modification - bins bigger than 10 are 'fake', not plotting so not calculating (too large, not likely to yield significant results)
## if you want to plot bonferroni instead of simpleM ..
# sizeBin2 <- c(number_of_snps[1:10], c(1:8))
# alphaBin_bonferroni = alpha * 2 ^ -(1:nbins) / sizeBin2
# 
# rep_helper <- c(table(out[,grp]))
# out[ , step2p_sig_bonferroni:=rep(alphaBin_bonferroni, rep_helper)]

## simpleM
sizeBin2 <- c(number_of_tests[1:10], c(1:8))
alphaBin_simpleM = alpha * 2 ^ -(1:nbins) / sizeBin2

rep_helper <- c(table(out[,grp]))
out[ , step2p_sig_simpleM:=rep(alphaBin_simpleM, rep_helper)]

create_twostep_plot_expectation_hybrid(dat = out, exposure = exposure, covars = covariates, binsToPlot = 10, statistic = simplem_step1_statistic, output_dir = output_dir, filename_suffix = "_expectation_based_simpleM")

# output data.frame of significant SNPs
out2 <- filter(out, grp %in% c(1:10) & step2p < step2p_sig_simpleM)
saveRDS(out2, file = paste0("/media/work/gwis/posthoc/", exposure, "/significant_results_dataframe_twostep_wht_",simplem_step1_statistic, "_", exposure, "_expectation_based_simpleM.rds"), version = 2)



























######################    

dat <- gxe
stats_step1 = 'chiSqG'
sizeBin0 = 5
alpha = 0.05

format_twostep_data <- function(dat, stats_step1, sizeBin0, alpha) {
  
  # quick function to calculate p values from chisq stats depending on method
  create_pval_info <- function(dat, stats_step1, df=1) {
    tmp <- data.table(dat)
    tmp[, step1p := pchisq(tmp[, get(stats_step1)], df = df, lower.tail = F)
    ][
      , step2p := pchisq(tmp[, get('chiSqGxE')],  df = 1, lower.tail = F)
    ][
      , y := -log10(step2p)
    ][
      order(step1p)
    ][
      , MapInfo := Location
    ]
  }
  
  if(stats_step1 == 'chiSqEDGE') {
    pv <- create_pval_info(dat, stats_step1, df = 2)
  } else {
    pv <- create_pval_info(dat, stats_step1, df = 1)
  }
  
  
  # format output for plotting..
  m = nrow(pv)
  nbins = floor(log2(m/sizeBin0 + 1)) # number of bins for second-step weighted Bonferroni correction
  nbins = if (m > sizeBin0 * 2^nbins) {nbins = nbins + 1}
  sizeBin = c(sizeBin0 * 2^(0:(nbins-2)), m - sizeBin0 * (2^(nbins-1) - 1) ) # bin sizes
  endpointsBin = cumsum(sizeBin) # endpoints of the bins
  alphaBin = alpha * 2 ^ -(1:nbins) / sizeBin # alpha elevel at which a SNP landing in bin 1,2,...nbins is tested
  
  sizeBin
  
  # add grp, wt (p value threshold for each bin), rename p value to 'y'
  rk.pv <- c(1:m)
  grp <- ceiling(log(rk.pv/sizeBin0+1,base=2)) # math shortcut to put bin size to every single marker by c(1:nrow(pv))
  pv[,grp:=grp] # assigning group to the p values..
  setkey(pv,grp)
  for(i in 1:max(grp))
  {
    pv[J(i),wt:=alpha*2^(-i)/nrow(pv[J(i)])] # data.table syntax, create threshold value
  }
  
  # return the data.table
  return(pv)
}


















create_twostep_weighted_plot <- function(dat, exposure, covars, sizeBin0, alpha, binsToPlot, statistic, filename_suffix = "") {
  
  cases <- unique(data.frame(dat[, 'Cases']))
  controls <- unique(data.frame(dat[, 'Subjects'])) - unique(data.frame(dat[, 'Cases']))
  total <- cases + controls
  
  # some plot title and file names
  write_twostep_weightedHT_plot_title <- function(statistic, exposure, covars, total) {
    gxescan_tests <- c(paste0("D|G 2-step Procedure Results (N = ", total, ")\noutc ~ G+", paste0(covars, collapse = "+"),"+", exposure),
                       paste0("G|E 2-step Procedure Results (N = ", total, ")\nG ~ ", exposure, "+", paste0(covars, collapse = "+")),
                       paste0("EDGE 2-step Procedure Results (N = ", total, ")\nchiSqG + chiSqGE"))
    names(gxescan_tests) <- c("chiSqG", "chiSqGE", "chiSqEDGE")
    return(gxescan_tests[statistic])
  }
  
  # write_twostep_weightedHT_plot_title("chiSqG", exposure, covars, total)
  
  # bin information
  m = nrow(dat)
  nbins = floor(log2(m/sizeBin0 + 1)) # number of bins for second-step weighted Bonferroni correction
  nbins = if (m > sizeBin0 * 2^nbins) {nbins = nbins + 1} # add +1 bin if condition met
  sizeBin = c(sizeBin0 * 2^(0:(nbins-2)), m - sizeBin0 * (2^(nbins-1) - 1) ) # bin sizes
  endpointsBin = cumsum(sizeBin) # endpoints of the bins
  alphaBin = alpha * 2 ^ -(1:nbins) / sizeBin # alpha level for each bbin 1,2, ... N bin tested
  
  # create list and vars for plotting
  min.p = 12 # this might be plot upper limit in -log10 scale, not sure why called 'min.p'
  last.sig = alphaBin[binsToPlot]
  
  
  # create list where each component contains a bin
  # log transform bin alpha value
  # create 'x', normalized position information for each bin
  
  dat <- data.table(dat)
  
  glist<-list()
  for(i in 1:binsToPlot){
    tmp <- dat[J(i)][order(Chromosome, Location)]
    tmp[, ref := -1*log10(min(tmp[,wt]))] # -log10 of bin specific alpha
    tmp[, mapinfo_tmp := seq(Chromosome, Chromosome+1, length.out = .N), by = Chromosome]
    tmp[, x := 0.8*((tmp[,mapinfo_tmp]-min(tmp[,mapinfo_tmp])) / (max(tmp[,mapinfo_tmp])-min(tmp[,mapinfo_tmp]))) + 0.1 + i - 1]
    glist[[i]]<-tmp[order(Chromosome, Location)]
    rm(tmp)
  }
  
  significant_hits <- dat[step2p <= wt]
  
  # CREATE PLOT
  head(glist[[1]]) # for reference
  
  png(paste0("/media/work/gwis/posthoc/", exposure, "/twostep_wht_", statistic, "_", exposure, filename_suffix, ".png"), height = 720, width = 1280)
  color <- rep(c("#377EB8","#4DAF4A"),100)
  par(mar=c(6, 7, 6, 3))
  plot(glist[[1]][,x], glist[[1]][,y],
       col = ifelse(glist[[1]][,SNP] %in% significant_hits[,SNP], '#E41A1C','#377EB8'),
       pch = ifelse(glist[[1]][,SNP] %in% significant_hits[,SNP], 19, 20),
       cex = ifelse(glist[[1]][,SNP] %in% significant_hits[,SNP], 1.3, 1.7),
       xlab="Bin number for step1 p value",
       ylab="-log10(step2 chiSqGxE p value)",
       xlim=c(0,binsToPlot),
       ylim=c(0,min.p),
       axes=F,
       cex.main = 1.7,
       cex.axis = 1.7,
       cex.lab = 1.7,
       cex.sub = 1.7)
  # cex = 1.2)
  lines(glist[[1]][,x], glist[[1]][,ref], col = "black", lwd=1)
  
  # the rest of the points
  for(i in 2:binsToPlot){
    points(glist[[i]][,x], glist[[i]][,y],
           col = ifelse(glist[[i]][,SNP] %in% significant_hits$SNP, '#E41A1C', color[i]),
           pch = ifelse(glist[[i]][,SNP] %in% significant_hits$SNP, 19, 20),
           cex = ifelse(glist[[1]][,SNP] %in% significant_hits[,SNP], 1.3, 1.7),
           cex.main = 1.7,
           cex.axis = 1.7,
           cex.lab = 1.7,
           cex.sub = 1.7)
    # cex = 1.2)
    lines(glist[[i]][,x], glist[[i]][,ref], col = "black",lwd = 1)
  }
  
  ## the last bin..
  ## it's this way because the last bin has smaller number of samples compared to bins before, thus lower bar
  ## let's only plot the first 15 bins for now, so change code a bit above.
  # points(glist[[num]][,x], glist[[num]][,y],
  #        col= color[num], pch = 19, cex = 0.5)
  # lines(glist[[num]][,x], rep(last.sig, nrow(glist[[num]])) ,col="black",lwd=2) # need to fix to create last horizontal line
  
  axis(1, at = c(-1.5, seq(0.5, binsToPlot-0.5, 1)), label = c(0, seq(1, binsToPlot, 1)), cex.axis = 1.7)
  axis(2, at = c(0:floor(min.p)), label = c(0:min.p), cex.axis=1.7)
  title(main = write_twostep_weightedHT_plot_title(statistic, exposure, covars, total), sub = "iBin Size = 5, alpha = 0.05", cex.main = 2, cex.sub = 1.7)
  
  dev.off()
}






# i think you need to understand how SVM model fitting works that anna used .. 
# because it looks like jsut cuz somethign is in an open region, doesn't mean that marker itself is scored high? 


huh <- filter(anna2, grepl("38402781", SNP))

60375355

huh <- filter(anna2, grepl("60375355", SNP))



# ----------------------------------------------------------___#
# functional analysis SUBSET PLOTS ----

# let's create the same plots but for the subset of individuals?






anna <- fread("/media/work/gwis/gecco.snps.in.peaks.svm.max.abs.delta.scores.txt") %>% 
  separate(SNP, into = c("SNP", "alleles"), sep = "_") %>% 
  separate(alleles, into = c("Reference", "Alternate"), sep = "/") %>% 
  mutate(tmp1 = paste0(SNP, ":", Reference, ":", Alternate), 
         tmp2 = paste0(SNP, ":", Alternate, ":", Reference))

gxescan_output <- readRDS(paste0("/media/work/gwis/results/", exposure, "/processed/FIGI_", hrc_version, "_gxeset_", exposure, "_basic_covars_gxescan_results.rds")) %>% 
  mutate(SNP2 = paste0(Chromosome, ":", Location))



head(gxescan_output)


test1 <- inner_join(gxescan_output, anna, by = c('SNP' = 'tmp1')) # ok good 
test2 <- inner_join(gxescan_output, anna, by = c('SNP' = 'tmp2'))




# let's get overlap, and conduct testing appropriately.. 
out <- inner_join(gxescan_output, anna, by = c('SNP' = 'tmp1')) %>% 
  mutate(p = pchisq(chiSqGxE, lower.tail = F, df = 1))

hist(out$deltaScore)
head(out)


manhattan(out, chr = "Chromosome", bp = "Location", p = "p", snp = "SNP", main = "GxE Results", cex.axis = 0.8, col = c("blue4", "orange3"), annotatePval = 0.05 / nrow(out))


# one question i have is - it seems like some gxe suggestive hits are in functionally active regions, shouldn't they be in this file.. 
wtf <- filter(anna, grepl("8:", SNP)) %>% 
  separate(SNP, c("chr", "bp")) %>% 
  mutate(bp = as.numeric(bp)) %>% 
  filter(between(bp, 60375355-100000, 60375355+100000))

# i think you need to remember 1) the smaller subset not sure how it was sourced, and 2) these are predictive scores based on the SVM model, not a positional subset .. 




# first pass, we could perhaps try a 2sd window for influential points? e.g. greater or less than 2 filter
anna2 <- fread("/media/work/gwis/gecco.snps.flora.svm.delta.scores.txt.gz") %>% 
  filter(abs(deltaScore) > 2)

anna3 <- filter(anna2, abs(deltaScore) > 3)

gxescan_output <- readRDS(paste0("/media/work/gwis/results/", exposure, "/processed/FIGI_", hrc_version, "_gxeset_", exposure, "_basic_covars_gxescan_results.rds"))


out <- inner_join(gxescan_output, anna2, by = c('SNP')) %>% 
  mutate(p = pchisq(chiSqGxE, lower.tail = F, df = 1),
         SNP2 = paste0(Chromosome, ":", Location))

annotatePval = -log10(0.05 / nrow(out))

manhattan(out, chr = "Chromosome", bp = "Location", p = "p", snp = "SNP", main = "GxE Results", cex.axis = 0.8, col = c("blue4", "orange3"), suggestiveline = F, genomewideline = annotatePval)

# multiple testing still stringent, we don't get significant results. need to use number of tests or something, but would be annoying to perform with 800k SNPs. 
table(out$Chromosome)

for(x in 1:22) {
  tmp <- filter(out , Chromosome == x) %>% pull(SNP)
  saveRDS(tmp, file = paste0("/media/work/gwis/twostep_expectation_hybrid/", exposure, "/", exposure, "_snplist_functional_annotation_test_", 'chiSqGxE', "_chr", x, ".rds"))
}


output_helper('chiSqGxE')




### i ran meff but the number of effective tests is still pretty high ...
# > number_of_tests
# [1] 30916 33838 27218 25460 25102 24615 22950 21862 18287 20489 19956 18949 14057 13064 12022 14489 12009 11589  9988 10281  5553  6259

meff = 398953


annotatePval = -log10(0.05 / meff)

manhattan(out, chr = "Chromosome", bp = "Location", p = "p", snp = "SNP", main = "GxE Results", cex.axis = 0.8, col = c("blue4", "orange3"), suggestiveline = -log10(5e-6), genomewideline = annotatePval)



# two step plot

tmp <- format_twostep_data(out, stats_step1 = 'chiSqG', sizeBin0 = 5, alpha = 0.05)
create_twostep_weighted_plot(tmp, exposure = 'folate_dietqc2', covars = covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = "chiSqG", filename_suffix = "_TESTING")

tmp <- format_twostep_data(out, stats_step1 = 'chiSqGE', sizeBin0 = 5, alpha = 0.05)
create_twostep_weighted_plot(tmp, exposure = 'folate_dietqc2', covars = covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = "chiSqGE", filename_suffix = "_TESTING")

tmp <- format_twostep_data(out, stats_step1 = 'chiSqEDGE', sizeBin0 = 5, alpha = 0.05)
create_twostep_weighted_plot(tmp, exposure = 'folate_dietqc2', covars = covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = "chiSqEDGE", filename_suffix = "_TESTING")





















create_qqplot <- function(dat, exposure, covars, stat, df, filename_suffix = "") {
  
  # calculate p value + lambda
  pvals <- calculate_pval(dat[, stat], df)
  lambda <- round( (median(qchisq(1-pvals, df)) / qchisq(0.5, df)), 4)
  
  # calculate lambda1000
  cases <- unique(dat[, 'Cases'])
  controls <- unique(dat[, 'Subjects']) - unique(dat[, 'Cases'])
  total <- cases + controls
  cases1000 <- (cases/total) * 1000
  controls1000 <- (controls/total) * 1000
  lambda1000 <- 1 + (lambda - 1) * ( (1/cases + 1/controls) / (1/(2*cases1000) + 1/(2*controls1000)))
  
  # plotting function
  png(paste0("/media/work/gwis/posthoc/", exposure, "/qq_", stat, "_", exposure, filename_suffix, ".png"), height = 720, width = 1280)
  qqman::qq(pvals,
            xlab = "Expected -log10(p)",
            ylab = "Observed -log10(p)",
            main = write_plot_title(stat, exposure, covars, total),
            cex.main = 1.6,
            cex.axis = 1.3,
            cex.lab = 1.3,
            cex.sub = 1.3,
            cex = 1.4,
            pch = 1,
            col = 'blue4')
  par(adj = 1)
  title(sub = bquote(lambda ~ '=' ~ .(signif(lambda, 4)) ~~ lambda[1000] ~ '=' ~.(signif(lambda1000, 4))), cex.sub = 1.3) # FYI ~~ adds spaces when using signif
  dev.off()
}







create_manhattanplot <- function(dat, exposure, stat, df, annotation_file, output_dir, filename_suffix = "") {
  
  dat <- dat %>%
    mutate(P = calculate_pval(.data[[stat]], df)) %>%
    filter(!(P > 0.01)) %>%
    dplyr::rename(CHR = Chromosome,
                  BP = Location) %>%
    dplyr::select(SNP, CHR, BP, P)
  
  write.table(dat, file = paste0("/media/work/tmp/manhattan_", stat, "_", exposure, filename_suffix), quote = F, row.names = F, sep = '\t')
  
  # create ecf file
  ecf1 <- paste0(output_dir)
  ecf2 <- "SNP;CHR;BP;P"
  ecf3 <- "character;numeric;numeric;numeric"
  ecf4 <- paste0("/media/work/tmp/manhattan_", stat, "_", exposure, filename_suffix)
  ecf_file_name <- paste0(output_dir,"EasyStrata_", stat, "_", exposure, filename_suffix, ".ecf")
  
  cat(paste0("DEFINE	--pathOut ", ecf1, "
      --acolIn ", ecf2, "
      --acolInClasses ", ecf3, "
      --strMissing NA
      --strSeparator TAB

      EASYIN	--fileIn ", ecf4, "

      START EASYX

      ################
      ## MHplot
      ################

      MHPLOT
      --colMHPlot P
      --colInChr CHR
      --colInPos BP
      --numWidth 1280
      --numHeight 720
      #--astrDefaultColourChr gray51;gray66
      --astrDefaultColourChr gray70;gray80
      --blnYAxisBreak 1
      --numYAxisBreak 22
      #--numPvalOffset 0.01
      # Annotation
      --fileAnnot /home/rak/data/Annotations/", annotation_file, "
      --numAnnotPosLim 1
      # Horizontal lines
      --anumAddPvalLine 1.25328e-07;5.698441e-08
      --anumAddPvalLineLty 6;6
      --astrAddPvalLineCol blue;red
      # Other Graphical Params
      --anumParMar 7;5;7;3
      --numDefaultSymbol 1
      --numDefaultCex 0.6
      --numCexAxis 1.7
      --numCexLab 1.7
      --arcdSymbolCrit  (P>5.698441e-08 & P<1.25328e-07);P<5.698441e-08
      --anumSymbol 20;19
      --arcdColourCrit (P>5.698441e-08 & P<1.25328e-07);P<5.698441e-08
      --astrColour gray55;black
      --arcdCexCrit (P>55.698441e-08 & P<1.25328e-07);P<5.698441e-08
      --anumCex 0.8;0.6

      STOP EASYX"), file = ecf_file_name, append = F)
  
  # run EasyStrata
  EasyStrata(ecf_file_name)
}



# SNP subset while also using hybrid expectation based.

# does 3df finding overlap with the functinal snps 
# using plink LD clump output, create vector of SNPs to remove (GWAS, D|G, or E|G + SNPs in LD)
format_clump <- function(x) {
  # x should be the plink output *.clumped read as dataframe
  tmp1 <- as.list(x$SP2)
  tmp2 <- sapply(tmp1, function(y) strsplit(y, split = ','))
  tmp3 <- do.call(c, tmp2) %>% 
    gsub("\\(1\\)", "", .) %>% 
    c(x$SNP) # add tag SNPs back into vector
}


check <- filter(out, grepl("77205974", Location)) # NOPE

exclude_gwas <- fread("~/data/Annotations/gwas_141_ld_annotation_july2020.txt") %>% 
  mutate(SNP2 = paste0(Chr, ":", Pos)) %>% 
  pull(SNP2)

clump_g_file <- paste0("/media/work/gwis/clump_combined/FIGI_", hrc_version, "_gxeset_", exposure, "_chiSqG_ldclump.clumped")
if(file.exists(clump_g_file)) {
  clump_g <- fread(clump_g_file)
} else {
  clump_g <- data.frame()
}

exclude_g    <- format_clump(clump_g)


out_nogwas_nomain <- out %>%
  dplyr::filter(!SNP2 %in% exclude_gwas,
                !SNP %in% exclude_g) %>% 
  mutate(p = pchisq(chiSq2df, lower.tail = F, df = 2)) %>% 
  filter(p < 0.05)
png(paste0("~/Dropbox/Presentations/folate_presentation/snp_subset_2df.png"), height = 600, width = 800)
manhattan(out_nogwas_nomain, chr = "Chromosome", bp = "Location", p = "p", snp = "SNP", main = "2DF Results - exclude GWAS and D|G", cex.axis = 0.8, col = c("blue4", "orange3"), suggestiveline = -log10(5e-6), genomewideline = -log10(5e-8))
dev.off()


out_nogwas_nomain <- out %>%
  dplyr::filter(!SNP2 %in% exclude_gwas,
                !SNP %in% exclude_g) %>% 
  mutate(p = pchisq(chiSq3df, lower.tail = F, df = 3)) %>% 
  filter(p < 0.05)

png(paste0("~/Dropbox/Presentations/folate_presentation/snp_subset_3df.png"), height = 600, width = 800)
manhattan(out_nogwas_nomain, chr = "Chromosome", bp = "Location", p = "p", snp = "SNP", main = "3DF Results - exclude GWAS and D|G", cex.axis = 0.8, col = c("blue4", "orange3"), suggestiveline = -log10(5e-6), genomewideline = -log10(5e-8))
dev.off()


out_nogwas_nomain <- out %>%
  # dplyr::filter(!SNP2 %in% exclude_gwas,
  #               !SNP %in% exclude_g) %>% 
  mutate(p = pchisq(chiSqGxE, lower.tail = F, df = 1)) %>% 
  filter(p < 0.05)
png(paste0("~/Dropbox/Presentations/folate_presentation/snp_subset_gxe.png"), height = 600, width = 800)
manhattan(out_nogwas_nomain, chr = "Chromosome", bp = "Location", p = "p", snp = "SNP", main = "GxE Results", cex.axis = 0.8, col = c("blue4", "orange3"), suggestiveline = -log10(5e-6), genomewideline = -log10(5e-8))
dev.off()





# let me make sure the 3df locus survives or something
out_nogwas_nomain <- out %>%
  dplyr::filter(!SNP2 %in% exclude_gwas,
                !SNP %in% exclude_g) %>% 
  mutate(p = pchisq(chiSq3df, lower.tail = F, df = 3)) %>% 
  filter(p < 5e-6)



# maybe some of these 3df are significant using meff cutoff
meff = 398953
annotatePval = -log10(0.05 / meff)



out_nogwas_nomain <- out %>%
  dplyr::filter(!SNP2 %in% exclude_gwas,
                !SNP %in% exclude_g) %>% 
  mutate(p = pchisq(chiSq2df, lower.tail = F, df = 2)) %>% 
  filter(p < 0.05)
png(paste0("~/Dropbox/Presentations/folate_presentation/snp_subset_2df_meff.png"), height = 600, width = 800)
manhattan(out_nogwas_nomain, chr = "Chromosome", bp = "Location", p = "p", snp = "SNP", main = "2DF Results - exclude GWAS and D|G", cex.axis = 0.8, col = c("blue4", "orange3"), suggestiveline = -log10(5e-6), genomewideline = annotatePval)
dev.off()


out_nogwas_nomain <- out %>%
  dplyr::filter(!SNP2 %in% exclude_gwas,
                !SNP %in% exclude_g) %>% 
  mutate(p = pchisq(chiSq3df, lower.tail = F, df = 3)) %>% 
  filter(p < 0.05)

png(paste0("~/Dropbox/Presentations/folate_presentation/snp_subset_3df_meff.png"), height = 600, width = 800)
manhattan(out_nogwas_nomain, chr = "Chromosome", bp = "Location", p = "p", snp = "SNP", main = "3DF Results - exclude GWAS and D|G", cex.axis = 0.8, col = c("blue4", "orange3"), suggestiveline = -log10(5e-6), genomewideline = annotatePval)
dev.off()


# one SNP
siggg <- filter(out_nogwas_nomain, p < 0.05/ meff)

# to properly plot, maybe you should only plot this subset
locuszoom_wrapper <- function(statistic, df) {
  # gxe defined above
  locuszoom <- gxe %>%
    dplyr::mutate(`P-value` = pchisq(.data[[statistic]], df = df, lower.tail = F),
                  MarkerName = paste0("chr", Chromosome, ":", Location)) %>%
    dplyr::select(MarkerName, `P-value`)
  
  write.table(locuszoom, file = paste0("/media/work/gwis/locuszoom/FIGI_", 
                                       hrc_version, "_gxeset_", exposure, 
                                       paste0("_basic_covars_gxescan_", statistic, "_locuszoom.txt")), 
              quote = F, row.names = F, sep = "\t")
}
locuszoom_wrapper('chiSq3df', df = 3)


annotatePval

#### previous folate x CRC findings ----
filter(gxescan_output, grepl("111171709", SNP))

previous_findings <- filter(gxescan_output, SNP %in% c("12:50880216:T:C", "11:111171709:C:A"))
saveRDS(previous_findings, file = paste0("/media/work/gwis/posthoc/gwis_sig_results_input_", exposure, "_additional_snps.rds"), version = 2)


# analysis
gwas_140 <- readRDS("/media/work/gwis/FIGI_genotype_dosages_gwas140.rds")
folate_2df <- readRDS("/media/work/gwis/posthoc/gwis_sig_results_output_folate_dietqc2_additional_snps.rds")

xxx <- inner_join(gwas_140, folate_2df, 'vcfid') %>% 
  inner_join(xx, 'vcfid') %>% 
  rename(rs7136702 = X12.50880216.T.C_dose, 
         rs3802842 = X11.111171709.C.A_dose) %>% 
  mutate(folate_diet400qcm = as.numeric(folate_diet400qcm))

# [1] "folate_diet"          "folate_sup"           "folate_diet400qcm"    "folate_dietqc2"       "folate_sup_yn"       
# [6] "folate_tot"           "folate_totqc2"        "folate_diet400qcm.ca" "folate_dietqc2.ca"    "folate_totqc2.ca" 

model_rs7136702 <- glm(outcome ~ folate_dietqc2 * rs7136702 + age_ref_imp + sex + energytot_imp + study_gxe + pc1 + pc2 + pc3, data = xxx, family = 'binomial')

tab_model(model_rs7136702, terms = c("(Intercept)", "folate_dietqc2", "age_ref_imp", "sex", "rs7136702", "folate_dietqc2:rs7136702"),
          p.style = "scientific", title = "GxE rs7136702")



model_rs3802842 <- glm(outcome ~ folate_dietqc2 * rs3802842 + age_ref_imp + sex + energytot_imp + study_gxe + pc1 + pc2 + pc3, data = xxx, family = 'binomial')
tab_model(model_rs3802842, terms = c("(Intercept)", "folate_dietqc2", "age_ref_imp", "sex", "rs3802842", "folate_dietqc2:rs3802842"),
          p.style = "scientific", title = "GxE rs3802842")


model_rs3802842 <- glm(outcome ~ folate_diet400qcm * rs3802842 + age_ref_imp + sex + energytot_imp + study_gxe + pc1 + pc2 + pc3, data = xxx, family = 'binomial')
tab_model(model_rs3802842, terms = c("(Intercept)", "folate_diet400qcm", "age_ref_imp", "sex", "rs3802842", "folate_diet400qcm:rs3802842"),
          p.style = "scientific", title = "GxE rs3802842")


head(anna2)
