rm(list = ls())
# source("functions.R")

# input variables
exposure = 'smk_ever'
hrc_version = 'v2.3'
annotation_file <- 'gwas_200_ld_annotation_feb2021.txt'
covariates <- sort(c('age_ref_imp', 'sex', 'study_gxe', 'pc1', 'pc2', 'pc3'))
path = glue("/media/work/gwis_test/{exposure}/")


# input data
esubset <- readRDS(glue("/media/work/gwis_test/{exposure}/data/FIGI_{hrc_version}_gxeset_{exposure}_basic_covars_glm.rds")) %>% 
  pull(vcfid)

input_data <- readRDS(glue("/media/work/gwis_test/data/FIGI_{hrc_version}_gxeset_analysis_data_glm.rds")) %>% 
  filter(vcfid %in% esubset)



dose <- qread("/media/work/gwis_test/smk_ever/output/posthoc/dosage_chr3_85461302.qs")
x <- inner_join(input_data, dose, 'vcfid')
model <- glm(outcome ~ smk_ever*chr3_85461302_G_A_dose + age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe, data = x, family = 'binomial')
summary(model)



dose <- qread("/media/work/gwis_test/smk_ever/output/posthoc/dosage_chr15_83613505.qs")
x <- inner_join(input_data, dose, 'vcfid')
model <- glm(outcome ~ smk_ever*chr15_83613505_A_T_dose + age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe, data = x, family = 'binomial')
summary(model)






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





dose <- qread("/media/work/gwis_test/smk_aveday/output/posthoc/dosage_chr8_138788813.qs")
x <- inner_join(input_data, dose, 'vcfid')
model <- glm(outcome ~ smk_aveday*chr8_138788813_C_A_dose + age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe, data = x, family = 'binomial')
summary(model)

