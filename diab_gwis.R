#=============================================================================#
# FIGI GxE diab results
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
exposure = 'diab'
hrc_version = 'v2.3'
# output_dir = paste0("/media/work/gwis/posthoc/", exposure, "/")
# annotation_file <- 'gwas_141_ld_annotation_july2020.txt'
annotation_file <- 'gwas_200_ld_annotation_feb2021.txt'
path = glue("/media/work/gwis_test/{exposure}/")


covariates <- sort(c('age_ref_imp', 'sex', 'study_gxe', 'pc1', 'pc2', 'pc3'))
covariates_set1 <- covariates
# covariates_set2 <- sort(c(covariates_set1, 'bmi', 'smk_ever', 'fruitqc2', 'vegetableqc2'))
covariates_list <- list(covariates)
mod <- 'age_ref_imp+sex+pc1+pc2+pc3+study_gxe'

# input data
esubset <- readRDS(glue("{path}/data/FIGI_{hrc_version}_gxeset_{exposure}_basic_covars_glm.rds")) %>% pull(vcfid)
input_data <- readRDS(glue("/media/work/gwis/data/FIGI_EpiData/FIGI_{hrc_version}_gxeset_analysis_data_glm.rds")) %>% 
  filter(vcfid%in% esubset)



input_data_bmi <- filter(input_data, 
                         !is.na(bmi))
table(input_data$outcome)
table(input_data_bmi$outcome)

#-----------------------------------------------------------------------------#
# main effects ----
#-----------------------------------------------------------------------------#

# ------ meta-analysis ------ #
covariates_meta <- sort(c('age_ref_imp', 'sex'))
meta_analysis_execute(dataset = input_data, exposure = exposure, hrc_version = hrc_version, covariates = covariates_meta, output_dir = output_dir, filename_suffix2 = "")

# BMI adjustment

tmp <- input_data %>% 
  filter(!is.na(bmi))

covariates_meta <- sort(c('age_ref_imp', 'sex', 'bmi'))
meta_analysis_execute(dataset = tmp, exposure = exposure, hrc_version = hrc_version, covariates = covariates_meta, output_dir = output_dir, filename_suffix2 = "")


tmp <- input_data %>% 
  filter(!is.na(bmi), 
         !is.na(redmeatqc2), 
         !is.na(smk_ever))

covariates_meta <- sort(c('age_ref_imp', 'sex', 'bmi', 'smk_ever', 'redmeatqc2'))
meta_analysis_execute(dataset = tmp, exposure = exposure, hrc_version = hrc_version, covariates = covariates_meta, output_dir = output_dir, filename_suffix2 = "")




# ------- pooled analysis ------- #
## basic covariates
covariates_pooled <- sort(covariates[which(!covariates %in% c(paste0(rep('pc', 20), seq(1, 20))))])
covariates_pooled <- covariates
pooled_analysis_glm(input_data, exposure = exposure, covariates = covariates_pooled, strata = 'sex', 
                    filename_suffix = paste0(paste0(sort(covariates), collapse = '_'), "_stratified_sex"), output_dir = output_dir)
pooled_analysis_glm(input_data, exposure = exposure, covariates = covariates_pooled, strata = 'study_design', 
                    filename_suffix = paste0(paste0(sort(covariates), collapse = '_'), "_stratified_study_design"), output_dir = output_dir)
pooled_analysis_multinom(input_data, exposure = exposure, covariates = covariates_pooled, 
                         filename_suffix = paste0(paste0(sort(covariates), collapse = '_'), "_stratified_cancer_site_sum2"), output_dir = output_dir)



#-----------------------------------------------------------------------------#
# main effects additional analysis ----
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




#-----------------------------------------------------------------------------#
# GxE additional analysis ---- 
#-----------------------------------------------------------------------------#

# output GxE models adjusted by different covariate sets
covariates_sets <- list(covariates, 
                        c(covariates, 'bmi5'), 
                        c(covariates, 'bmi5', 'smk_ever'), 
                        c(covariates, 'bmi5', 'smk_ever', 'fruitqc2', 'vegetableqc2'))

fit_gxe_covars(data_epi = input_data, exposure = exposure, snp = "13:47191972:G:A", covariates_list = covariates_sets, method = 'chiSqGxE', path = glue("{path}/output"))


# additional covariates is making association more significant.. let's generate for all suggestive hits
suggestive_gxe <- fread(glue("{path}/data/FIGI_v2.3_gxeset_diab_chiSqGxE_ldclump.clumped"))
walk(suggestive_gxe$SNP, ~ fit_gxe_covars(data_epi = input_data, exposure = exposure, snp = .x, covariates_list = covariates_sets, method = 'chiSqGxE', path = glue("{path}/output")))


# output RERI plots (can't install package on CARC yet)
reri_wrapper(data_epi = input_data, exposure = exposure, snp = "13:47191972:G:A", covariates = covariates, path = glue("{path}/output"))
reri_wrapper(data_epi = input_data, exposure = exposure, snp = "8:118185025:G:A", covariates = covariates, path = glue("{path}/output"))



# ---- calculate 2DF LR TEST ---- #
# ---- chr13_47191972_G_A_dose ---- #

# need allele frequency by study_gxe to confirm finding (and sensitivity analysis)
tmp <- qread("/media/work/gwis_test/diab/output/posthoc/dosage_chr13_47191972.qs") %>% 
  inner_join(input_data, 'vcfid')


# calculate 2DF original to confirm
covariates_original = covariates
model_base <- glm(glue("outcome ~ diab + {glue_collapse(covariates_original, sep = '+')}"), data = tmp, family = 'binomial')
model_gxe <- glm(glue("outcome ~ diab * chr13_47191972_G_A_dose + {glue_collapse(covariates_original, sep = '+')}"), data = tmp, family = 'binomial')
lrtest(model_base, model_gxe)

covariates_bmi = c(covariates_original, 'bmi')
model_base <- glm(glue("outcome ~ diab + {glue_collapse(covariates_bmi, sep = '+')}"), data = tmp, family = 'binomial')
model_gxe <- glm(glue("outcome ~ diab * chr13_47191972_G_A_dose + {glue_collapse(covariates_bmi, sep = '+')}"), data = tmp, family = 'binomial')
lrtest(model_base, model_gxe)

tmp %>% 
  summarise(total = n(), 
            study_aaf = sum(chr13_47191972_G_A_dose) / (total*2))

aaf <- function(x) {
  sum(x) / nrow(x)
}

out <- tmp %>% 
  group_by(study_gxe) %>% 
  summarise(total = n(), 
            study_aaf = sum(chr13_47191972_G_A_dose) / (total*2)) %>% 
  arrange(study_aaf) %>% 
  mutate(study_gxe = fct_reorder(study_gxe, study_aaf))

ggplot(aes(x = study_gxe, y = study_aaf), data = out) + 
  geom_point() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 270)) + 
  xlab("Study") + 
  ylab("Alternate Allele Frequency")










# ---- calculate 3DF LR TEST ---- #
# ---- chr13_47191972_G_A_dose ---- #

# need allele frequency by study_gxe to confirm finding (and sensitivity analysis)
tmp <- qread("/media/work/gwis_test/diab/output/posthoc/dosage_chr8_118185025.qs") %>% 
  inner_join(input_data, 'vcfid')

tmp %>% 
  summarise(total = n(), 
            study_aaf = sum(chr8_118185025_G_A_dose) / (total*2))

out <- tmp %>% 
  group_by(study_gxe) %>% 
  summarise(total = n(), 
            study_aaf = sum(chr8_118185025_G_A_dose) / (total*2)) %>% 
  arrange(study_aaf) %>% 
  mutate(study_gxe = fct_reorder(study_gxe, study_aaf))

ggplot(aes(x = study_gxe, y = study_aaf), data = out) + 
  geom_point() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 270)) + 
  xlab("Study") + 
  ylab("Alternate Allele Frequency")


# calculate 3DF original to confirm

covariates_original = covariates

model_base <- glm(glue("outcome ~ diab + {glue_collapse(covariates_original, sep = '+')}"), data = tmp, family = 'binomial')
model_g <- glm(glue("outcome ~ diab + chr8_118185025_G_A_dose + {glue_collapse(covariates_original, sep = '+')}"), data = tmp, family = 'binomial')
lrtest(model_base, model_g)

model_base_gxe <- glm(glue("outcome ~ diab + chr8_118185025_G_A_dose + {glue_collapse(covariates_original, sep = '+')}"), data = tmp, family = 'binomial')
model_gxe <- glm(glue("outcome ~ diab * chr8_118185025_G_A_dose + {glue_collapse(covariates_original, sep = '+')}"), data = tmp, family = 'binomial')
lrtest(model_base_gxe, model_gxe)

model_base_eg <- lm(glue("chr8_118185025_G_A_dose ~ {glue_collapse(covariates_original, sep = '+')}"), data = tmp)
model_eg <- lm(glue("chr8_118185025_G_A_dose ~ diab + {glue_collapse(covariates_original, sep = '+')}"), data = tmp)
lrtest(model_base_eg, model_eg)


covariates_bmi = c(covariates_original, 'bmi')
model_base <- glm(glue("outcome ~ diab + {glue_collapse(covariates_bmi, sep = '+')}"), data = tmp, family = 'binomial')
model_g <- glm(glue("outcome ~ diab + chr8_118185025_G_A_dose + {glue_collapse(covariates_bmi, sep = '+')}"), data = tmp, family = 'binomial')
lrtest(model_base, model_g)

model_base_gxe <- glm(glue("outcome ~ diab + chr8_118185025_G_A_dose + {glue_collapse(covariates_bmi, sep = '+')}"), data = tmp, family = 'binomial')
model_gxe <- glm(glue("outcome ~ diab * chr8_118185025_G_A_dose + {glue_collapse(covariates_bmi, sep = '+')}"), data = tmp, family = 'binomial')
lrtest(model_base_gxe, model_gxe)

model_base_eg <- lm(glue("chr8_118185025_G_A_dose ~ {glue_collapse(covariates_bmi, sep = '+')}"), data = tmp)
model_eg <- lm(glue("chr8_118185025_G_A_dose ~ diab + {glue_collapse(covariates_bmi, sep = '+')}"), data = tmp)
lrtest(model_base_eg, model_eg)

pchisq(0.4015+9.696+39.23, df = 3, lower.tail = F)



# AAF for suggestive hits
aaf_wrapper <- function(chr, snp) {
  tmp <- qread(glue("/media/work/gwis_test/diab/output/posthoc/dosage_chr{chr}_{snp}.qs")) %>% 
    inner_join(input_data, 'vcfid')
  
  tmp %>% 
    summarise(total = n(), 
              study_aaf = sum(.[, 1]) / (total*2))
}

dosage_chr13_47191972.qs
chr8_118185025_G_A_dose
aaf_wrapper(13, 47191972)
aaf_wrapper(8, 118185025)

aaf_wrapper(2, 3731978)
aaf_wrapper(3, 193408847)
aaf_wrapper(4, 31909491)
aaf_wrapper(5, 60920832)
aaf_wrapper(13, 47170118)
# ================================================================== #
# ======= interaction stratified by diabetes status ---- 
# i forget the rationale for this analysis, maybe it's to check D|G association? 
# ================================================================== #

snps_filter <- snps %>% 
  dplyr::filter(method == glue("manhattan_chiSqCase_{exposure}_clump_df")) %>% 
  dplyr::mutate(snps = paste0('chr', gsub("\\:", "\\_", SNP))) %>% 
  dplyr::arrange(Pval) %>% 
  pull(snps)



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






# ================================================================== #
# ======= additional stuff ---- 
# ================================================================== #

# additional stuff
# 
# 2DF - chr13_47191972_G_A
table(figi$diab)

figi_nodiab <- figi %>% 
  dplyr::filter(diab == "No")

figi_diab <- figi %>% 
  dplyr::filter(diab == "Yes")

model1 <- glm(outcome ~ chr13_47191972_G_A + age_ref_imp + sex + pc1 + pc2 + pc3 , data = figi_nodiab, family = 'binomial')
model2 <- glm(outcome ~ chr13_47191972_G_A + age_ref_imp + sex + pc1 + pc2 + pc3 , data = figi_diab,   family = 'binomial')

summary(model1)
summary(model2)

mht_hits <- fread("~/mht_gxe_previously_reported.txt", stringsAsFactors = F) %>% 
  mutate(SNP2 = paste0(V1, ":", V2))

mht_hits_stats <- gxe %>% 
  filter(SNP2 %in% mht_hits$SNP2)

write.csv(mht_hits_stats, file = "~/Dropbox/mht_replication_stats.csv", quote = F, row.names = F)



# 3DF - chr8_118185025_G_A
table(figi$diab)

figi_nodiab <- figi %>% 
  dplyr::filter(diab == "No")

figi_diab <- figi %>% 
  dplyr::filter(diab == "Yes")

model1 <- glm(outcome ~ chr8_118185025_G_A + age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe, data = figi_nodiab, family = 'binomial')
model2 <- glm(outcome ~ chr8_118185025_G_A + age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe, data = figi_diab,   family = 'binomial')

summary(model1)
summary(model2)





mht_hits <- fread("~/mht_gxe_previously_reported.txt", stringsAsFactors = F) %>% 
  mutate(SNP2 = paste0(V1, ":", V2))

mht_hits_stats <- gxe %>% 
  filter(SNP2 %in% mht_hits$SNP2)

write.csv(mht_hits_stats, file = "~/Dropbox/mht_replication_stats.csv", quote = F, row.names = F)







