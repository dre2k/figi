#=============================================================================#
# FIGI GxE asp_ref results
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
exposure = 'comp2'
hrc_version = 'v2.3'
annotation_file <- 'gwas_200_ld_annotation_feb2021.txt'
covariates <- sort(c('age_ref_imp', 'sex', 'study_gxe', 'pc1', 'pc2', 'pc3'))
path = glue("/media/work/gwis_test/{exposure}/")


# input data
esubset <- readRDS(glue("/media/work/gwis_test/{exposure}/data/FIGI_{hrc_version}_gxeset_{exposure}_basic_covars_gxescan.rds")) %>% 
  pull(vcfid)

comp2var <- readRDS(glue("/media/work/gwis_test/{exposure}/data/FIGI_{hrc_version}_gxeset_{exposure}_basic_covars_gxescan.rds")) %>% 
  dplyr::select(vcfid, Comp.2)

input_data <- readRDS(glue("/media/work/gwis_test/data/FIGI_{hrc_version}_gxeset_analysis_data_glm.rds")) %>% 
  inner_join(comp2var, 'vcfid') %>% 
  mutate(comp2 = Comp.2)

# saveRDS(input_data, file = glue("/media/work/gwis_test/{exposure}/data/FIGI_{hrc_version}_gxeset_{exposure}_basic_covars_glm.rds"))




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
# GxE additional analysis ---- 
#-----------------------------------------------------------------------------#


# need allele frequency by study_gxe to confirm finding (and sensitivity analysis)
tmp <- qread("/media/work/gwis_test/comp2/output/posthoc/dosage_chr7_13816598.qs") %>% 
  inner_join(input_data, 'vcfid')

aaf <- function(x) {
  sum(x) / nrow(x)
}

out <- tmp %>% 
  group_by(study_gxe) %>% 
  summarise(total = n(), 
            study_aaf = sum(chr7_13816598_A_C_dose) / (total*2)) %>% 
  arrange(study_aaf) %>% 
  mutate(study_gxe = fct_reorder(study_gxe, study_aaf))

ggplot(aes(x = study_gxe, y = study_aaf), data = out) + 
  geom_point() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 270)) + 
  xlab("Study") + 
  ylab("Alternate Allele Frequency")


# would results change if you remove mecc (note wald statistic)
model1 <- glm(glue("outcome ~ {exposure}*chr7_13816598_A_C_dose + {glue_collapse(covariates, sep = '+')}"), data = tmp, family = 'binomial')
summary(model1)

tmp2 <- dplyr::filter(tmp, !grepl("MECC_2", study_gxe))
model2 <- glm(glue("outcome ~ {exposure}*chr7_13816598_A_C_dose + {glue_collapse(covariates, sep = '+')}"), data = tmp2, family = 'binomial')
summary(model2)

tmp3 <- dplyr::filter(tmp, grepl("UKB_1", study_gxe))
covariates_nostudy <- covariates[which(covariates != "study_gxe")]
model3 <- glm(glue("outcome ~ {exposure}*chr14_74029409_C_T_dose + {glue_collapse(covariates_nostudy, sep = '+')}"), data = tmp3, family = 'binomial')
summary(model3)








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





# ================================================================== #
# ======= interaction test with gwas hits ---- 
# ================================================================== #

x <- qread("/media/work/gwis_test/gwas/output/dosage_chr1_110365045.qs")

out <- inner_join(x, input_data, 'vcfid')

model_formula <- glue("outcome ~ {names(out)[1]}*comp6 + age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe")

model <- glm(model_formula, data = out, family = 'binomial')

tidy(model) %>% 
  dplyr::filter(grepl(":", term))


# wraper

snp_list <- list.files("/media/work/gwis_test/gwas/output/", full.names = T)


wrap <- function(x) {
  dose <- qread(x)
  out <- inner_join(dose, input_data, 'vcfid')
  model_formula <- glue("outcome ~ {names(out)[1]}*comp2 + age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe")
  
  tidy(glm(model_formula, data = out, family = 'binomial')) %>% 
    dplyr::filter(grepl(":", term)) 
}


output <- map(snp_list, ~ wrap(.x))
outputf <- do.call(bind_rows, output)
write.csv(outputf, file = "~/Dropbox/FIGI_comp2_gxe_knownhits.csv", quote= F, row.names = F)

