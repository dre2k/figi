library(tidyverse)
library(glue)
library(qs)
library(lmtest)

# the results you sent jim include dachs!

# fruit5qcm 3df finding information

dose <- qread("/media/work/gwis_test/fruit5qcm/output/posthoc/dosage_chr1_72729142.qs")


# input variables
exposure = 'fruit5qcm'
hrc_version = 'v3.0'
covariates <- sort(c('age_ref_imp', 'sex', 'energytot_imp', 'study_gxe', 'pc1', 'pc2', 'pc3'))
path = glue("/media/work/gwis_test/{exposure}/")


# input data
esubset <- readRDS(glue("{path}/data/FIGI_{hrc_version}_gxeset_{exposure}_basic_covars_glm.rds")) %>% pull(vcfid)
input_data <- readRDS(glue("/media/work/gwis/data/FIGI_EpiData/FIGI_{hrc_version}_gxeset_analysis_data_glm.rds")) %>% 
  # filter(vcfid %in% esubset) %>% 
  mutate(fruit5qcm = as.numeric(fruit5qcm))


input_data <- readRDS("/media/work/gwis_test/fruit5qcm/data/FIGI_v3.0_gxeset_fruit5qcm_basic_covars_glm.rds")
input_data <- readRDS("/media/work/gwis_test/fruit5qcm_includeDACHS/data/FIGI_v3.0_gxeset_fruit5qcm_basic_covars_glm.rds")



out <- inner_join(input_data, dose, 'vcfid')



head(out)


# gxe model (1df)

model1 <- glm(glue("outcome ~ {exposure}*chr1_72729142_A_G_dose + {glue_collapse(covariates, sep = '+')}"), data = out, family = 'binomial')
summary(model1)

model_base <- glm(glue("outcome ~ {exposure}+chr1_72729142_A_G_dose + {glue_collapse(covariates, sep = '+')}"), data = out, family = 'binomial')
pval <- lrtest(model1, model_base)
pval$`Pr(>Chisq)`


# GxE model (case only)
model1 <- lm(glue("{exposure} ~ chr1_72729142_A_G_dose + {glue_collapse(covariates, sep = '+')}"), data = out %>% filter(outcome == 1))
summary(model1)

model_base <- lm(glue("{exposure} ~ {glue_collapse(covariates, sep = '+')}"), data = out %>% filter(outcome == 1))
pval <- lrtest(model1, model_base)
pval$`Pr(>Chisq)`


# GxE model (control only)
model1 <- lm(glue("{exposure} ~ chr1_72729142_A_G_dose + {glue_collapse(covariates, sep = '+')}"), data = out %>% filter(outcome == 0))
summary(model1)

model_base <- lm(glue("{exposure} ~ {glue_collapse(covariates, sep = '+')}"), data = out %>% filter(outcome == 0))
pval <- lrtest(model1, model_base)
pval$`Pr(>Chisq)`






# additional request
# Basically same structure as your case only and control only models (linear regression with E as the outcome)…but
# apply to all subjects and include a disease indicator (1=case, 0=control) and an interaction between G and D. 
# 
# A slightly different way to look at GxE interaction that would be roughly equivalent to standard 1 df GxE if G and E were binary (and we used logistic reg for the case-only/control only modeling)…but may yield something different when G and E are continuous.

# case only
model1 <- lm(glue("{exposure} ~ outcome*chr1_72729142_A_G_dose + {glue_collapse(covariates, sep = '+')}"), data = out)
summary(model1)


