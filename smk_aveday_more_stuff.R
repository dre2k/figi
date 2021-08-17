#=============================================================================#
# stratified odds ratio but more generalized
# (allow single unit increment for continuous variable)
#=============================================================================#

# goddamn i did this at one point but can't find it. get your shit together

#=============================================================================#
# FIGI GxE alcoholc_moderate results
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
exposure = 'smk_aveday'
hrc_version = 'v2.3'
annotation_file <- 'gwas_200_ld_annotation_feb2021.txt'
path = glue("/media/work/gwis_test/{exposure}/")

covariates <- sort(c('age_ref_imp', 'sex', 'study_gxe', 'pc1', 'pc2', 'pc3'))
covariates_set1 <- covariates
covariates_set2 <- c(covariates_set1, "smk_ever", "bmi", "famhx1")
covariates_list <- list(covariates_set1, covariates_set2)
mod <- 'age_ref_imp+sex+pc1+pc2+pc3+study_gxe'


# input data
esubset <- readRDS(glue("{path}/data/FIGI_{hrc_version}_gxeset_{exposure}_basic_covars_glm.rds")) %>% pull(vcfid)
input_data <- readRDS(glue("/media/work/gwis/data/FIGI_EpiData/FIGI_{hrc_version}_gxeset_analysis_data_glm.rds")) %>% 
  filter(vcfid%in% esubset)


# for rs2300985 chr10_101476905_G_A
tmp2 <- qread(glue("/media/work/gwis_test/{exposure}/output/posthoc/dosage_chr6_31917540.qs")) %>% 
  inner_join(input_data, 'vcfid') %>% 
  mutate(p0 = round(chr6_31917540_T_C_p0,0),
         p1 = round(chr6_31917540_T_C_p1,0),
         p2 = round(chr6_31917540_T_C_p2,0))


mod <- 'age_ref_imp+sex+pc1+pc2+pc3+study_gxe'
mModel <- paste0('outcome~',exposure,'*chr6_31917540_T_C_p1+',exposure,'*chr6_31917540_T_C_p2+', mod)
model <- glm(mModel, data = tmp2, family = 'binomial')


tmp <- summary(model)
COV <- tmp$cov.unscaled


# stratified by GE
# e1 <- 'alcoholc1-28g/d'
e1 <- 'smk_aveday'
p1 <- 'chr6_31917540_T_C_p1'
p2 <- 'chr6_31917540_T_C_p2'


e1g0 <- tmp$coef[e1,c(1,2)] # est + SE

beta.e1g1 <- tmp$coef[e1,1]+
  tmp$coef[p1,1]+
  tmp$coef[glue(e1, ':', p1),1]

beta.e1g2 <- tmp$coef[e1,1]+
  tmp$coef[p2,1]+
  tmp$coef[glue(e1, ':', p2),1]

e0g <- tmp$coef[c(p1, p2), c(1,2)] # getting both e0g1 and e0g2 estimates (matrix)

I3 =  rep(1,3)
covs = COV[c(e1, p1 , glue(e1, ':', p1)), c(e1, p1, glue(e1, ':', p1))]
se.e1g1 = sqrt(t(I3)%*%covs%*%I3)

covs = COV[c(e1, p2 , glue(e1, ':', p2)), c(e1, p2, glue(e1, ':', p2))]
se.e1g2 = sqrt(t(I3)%*%covs%*%I3)

# calculate ODDS ratios

rnd2 <- function(x) {round(x, 2)}

OR_e1g0 = glue("{rnd2(exp(e1g0[1]))} ({rnd2(exp(e1g0[1] - qnorm(0.975)*e1g0[2]))} - {rnd2(exp(e1g0[1] + qnorm(0.975)*e1g0[2]))})")
OR_e0g1 = glue("{rnd2(exp(e0g[1,1]))} ({rnd2(exp(e0g[1,1] - qnorm(0.975)*e0g[1,2]))} - {rnd2(exp(e0g[1,1] + qnorm(0.975)*e0g[1,2]))})")
OR_e0g2 = glue("{rnd2(exp(e0g[2,1]))} ({rnd2(exp(e0g[2,1] - qnorm(0.975)*e0g[2,2]))} - {rnd2(exp(e0g[2,1] + qnorm(0.975)*e0g[2,2]))})")

OR_e1g1 = glue("{rnd2(exp(beta.e1g1))} ({rnd2(exp(beta.e1g1 - qnorm(0.975)*se.e1g1))} - {rnd2(exp(beta.e1g1 + qnorm(0.975)*se.e1g1))})")
OR_e1g2 = glue("{rnd2(exp(beta.e1g2))} ({rnd2(exp(beta.e1g2 - qnorm(0.975)*se.e1g2))} - {rnd2(exp(beta.e1g2 + qnorm(0.975)*se.e1g2))})")



# E stratified by G
beta.e1g1 <- tmp$coef[e1,1]+
  tmp$coef[glue(e1, ':', p1),1]

beta.e1g2 <- tmp$coef[e1,1]+
  tmp$coef[glue(e1, ':', p2),1]


I3 = rep(1,2)

covs = COV[c(e1, glue(e1, ':', p1)), c(e1, glue(e1, ':', p1))]
se.e1g1 = sqrt(t(I3)%*%covs%*%I3)

covs = COV[c(e1, glue(e1, ':', p2)), c(e1, glue(e1, ':', p2))]
se.e1g2 = sqrt(t(I3)%*%covs%*%I3)



OR_e1g1 = glue("{rnd2(exp(beta.e1g1))} ({rnd2(exp(beta.e1g1 - qnorm(0.975)*se.e1g1))} - {rnd2(exp(beta.e1g1 + qnorm(0.975)*se.e1g1))})")
OR_e1g2 = glue("{rnd2(exp(beta.e1g2))} ({rnd2(exp(beta.e1g2 - qnorm(0.975)*se.e1g2))} - {rnd2(exp(beta.e1g2 + qnorm(0.975)*se.e1g2))})")


OR_e1g1
OR_e1g2


# G stratified by E

beta.e1g1 <- tmp$coef[p1,1]+
  tmp$coef[glue(e1, ':', p1),1]

beta.e1g2 <- tmp$coef[p2,1]+
  tmp$coef[glue(e1, ':', p2),1]


I3 = rep(1,2)

covs = COV[c(p1, glue(e1, ':', p1)), c(p1, glue(e1, ':', p1))]
se.e1g1 = sqrt(t(I3)%*%covs%*%I3)

covs = COV[c(p2, glue(e1, ':', p2)), c(p2, glue(e1, ':', p2))]
se.e1g2 = sqrt(t(I3)%*%covs%*%I3)


OR_e1g1 = glue("{rnd2(exp(beta.e1g1))} ({rnd2(exp(beta.e1g1 - qnorm(0.975)*se.e1g1))} - {rnd2(exp(beta.e1g1 + qnorm(0.975)*se.e1g1))})")
OR_e1g2 = glue("{rnd2(exp(beta.e1g2))} ({rnd2(exp(beta.e1g2 - qnorm(0.975)*se.e1g2))} - {rnd2(exp(beta.e1g2 + qnorm(0.975)*se.e1g2))})")

OR_e1g1
OR_e1g2



