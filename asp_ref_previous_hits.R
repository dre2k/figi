# ================================================================== #
# ======= check previously published hits ---- 
# ======= asp_ref 
# ================================================================== #

# input variables
exposure = 'asp_ref'
hrc_version = 'v2.4'
covariates <- sort(c('age_ref_imp', 'sex', 'study_gxe', 'pc1', 'pc2', 'pc3'))
path = glue("/media/work/gwis_test/{exposure}/")
esubset <- readRDS(glue("/media/work/gwis_test/{exposure}/data/FIGI_{hrc_version}_gxeset_{exposure}_basic_covars_glm.rds")) %>% 
  pull(vcfid)
input_data <- readRDS(glue("/media/work/gwis_test/data/FIGI_{hrc_version}_gxeset_analysis_data_glm.rds")) %>% 
  filter(vcfid%in% esubset)

sink('~/Dropbox/asp_ref_previous_gxe.txt')

tmp <- qread("/media/work/gwis_test/asp_ref/output/posthoc/dosage_chr12_17444733.qs")
out <- inner_join(input_data, tmp, 'vcfid')
model <- glm(outcome ~ asp_ref*chr12_17444733_A_T_dose + sex + age_ref_imp + pc1 + pc2 + pc3 + study_gxe, data = out, family = 'binomial')
summary(model)

tmp <- qread("/media/work/gwis_test/asp_ref/output/posthoc/dosage_chr12_17488764.qs")
out <- inner_join(input_data, tmp, 'vcfid')
model <- glm(outcome ~ asp_ref*chr12_17488764_A_T_dose + sex + age_ref_imp + pc1 + pc2 + pc3 + study_gxe, data = out, family = 'binomial')
summary(model)

tmp <- qread("/media/work/gwis_test/asp_ref/output/posthoc/dosage_chr15_82229999.qs")
out <- inner_join(input_data, tmp, 'vcfid')
model <- glm(outcome ~ asp_ref*chr15_82229999_A_C_dose + sex + age_ref_imp + pc1 + pc2 + pc3 + study_gxe, data = out, family = 'binomial')
summary(model)

sink()

# ================================================================== #
# ======= check previously published hits ---- 
# ======= aspirin 
# ================================================================== #


# input variables
exposure = 'aspirin'
hrc_version = 'v2.4'
covariates <- sort(c('age_ref_imp', 'sex', 'study_gxe', 'pc1', 'pc2', 'pc3'))
path = glue("/media/work/gwis_test/{exposure}/")
esubset <- readRDS(glue("/media/work/gwis_test/{exposure}/data/FIGI_{hrc_version}_gxeset_{exposure}_basic_covars_glm.rds")) %>% 
  pull(vcfid)
input_data <- readRDS(glue("/media/work/gwis_test/data/FIGI_{hrc_version}_gxeset_analysis_data_glm.rds")) %>% 
  filter(vcfid%in% esubset)


sink('~/Dropbox/aspirin_previous_gxe.txt')

tmp <- qread("/media/work/gwis_test/asp_ref/output/posthoc/dosage_chr12_17444733.qs")
out <- inner_join(input_data, tmp, 'vcfid')
model <- glm(outcome ~ aspirin*chr12_17444733_A_T_dose + sex + age_ref_imp + pc1 + pc2 + pc3 + study_gxe, data = out, family = 'binomial')
summary(model)

tmp <- qread("/media/work/gwis_test/asp_ref/output/posthoc/dosage_chr12_17488764.qs")
out <- inner_join(input_data, tmp, 'vcfid')
model <- glm(outcome ~ aspirin*chr12_17488764_A_T_dose + sex + age_ref_imp + pc1 + pc2 + pc3 + study_gxe, data = out, family = 'binomial')
summary(model)

tmp <- qread("/media/work/gwis_test/asp_ref/output/posthoc/dosage_chr15_82229999.qs")
out <- inner_join(input_data, tmp, 'vcfid')
model <- glm(outcome ~ aspirin*chr15_82229999_A_C_dose + sex + age_ref_imp + pc1 + pc2 + pc3 + study_gxe, data = out, family = 'binomial')
summary(model)

sink()





