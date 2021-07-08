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
exposure = 'alcoholc_moderate'
hrc_version = 'v2.3'
annotation_file <- 'gwas_200_ld_annotation_feb2021.txt'
path = glue("/media/work/gwis_test/{exposure}/")

covariates <- sort(c('age_ref_imp', 'sex', 'energytot_imp', 'study_gxe', 'pc1', 'pc2', 'pc3'))
covariates_set1 <- covariates
covariates_set2 <- c(covariates_set1, "smk_ever", "bmi", "famhx1")
covariates_list <- list(covariates_set1, covariates_set2)
mod <- 'age_ref_imp+sex+energytot_imp+pc1+pc2+pc3+study_gxe'


# input data
esubset <- readRDS(glue("{path}/data/FIGI_{hrc_version}_gxeset_{exposure}_basic_covars_glm.rds")) %>% pull(vcfid)
input_data <- readRDS(glue("/media/work/gwis/data/FIGI_EpiData/FIGI_{hrc_version}_gxeset_analysis_data_glm.rds")) %>% 
  filter(vcfid%in% esubset)











#-----------------------------------------------------------------------------#
# SNP followup (models and plots) ---- 
#-----------------------------------------------------------------------------#

# make sure 1-28g/d is the reference group !
figi <- posthoc_input(exposure, hrc_version, glue('gwis_sig_results_output_{exposure}.rds')) %>% 
  mutate(alcoholc_moderate = fct_relevel(alcoholc_moderate, "1-28g/d", "nondrinker"), 
         alcoholc_moderate = as.numeric(alcoholc_moderate)) # better to make it numeric.. 

figi <- posthoc_input(exposure, hrc_version, glue('gwis_sig_results_output_{exposure}.rds')) %>% 
  mutate(alcoholc_moderate = fct_relevel(alcoholc_moderate, "1-28g/d", "nondrinker"))


snps <- readRDS(glue(output_dir, "gwis_sig_results_input_{exposure}.rds"))
table(snps$method)





# ================================================================== #
# ======= AAF check ---- 
# ================================================================== #



# for rs2300985 chr10_101476905_G_A
tmp <- qread(glue("/media/work/gwis_test/{exposure}/output/posthoc/dosage_chr10_101476905.qs")) %>% 
  inner_join(input_data, 'vcfid')

tmp %>% 
  summarise(total = n(), 
            study_aaf = sum(chr10_101476905_G_A_dose) / (total*2))

out <- tmp %>% 
  group_by(study_gxe) %>% 
  summarise(total = n(), 
            study_aaf = sum(chr10_101476905_G_A_dose) / (total*2)) %>% 
  arrange(study_aaf) %>% 
  mutate(study_gxe = fct_reorder(study_gxe, study_aaf))

ggplot(aes(x = study_gxe, y = study_aaf), data = out) + 
  geom_point() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 270)) + 
  xlab("Study") + 
  ylab("Alternate Allele Frequency")




# ------------------ #
# stuff for duncan
# -------------------#




# ----- alcoholc_moderate ------ #

# input variables
exposure = 'alcoholc_moderate'
hrc_version = 'v2.3'
annotation_file <- 'gwas_200_ld_annotation_feb2021.txt'
path = glue("/media/work/gwis_test/{exposure}/")
covariates <- sort(c('age_ref_imp', 'sex', 'energytot_imp', 'study_gxe', 'pc1', 'pc2', 'pc3'))

# input data
esubset <- readRDS(glue("{path}/data/FIGI_{hrc_version}_gxeset_{exposure}_basic_covars_glm.rds")) %>% pull(vcfid)
input_data <- readRDS(glue("/media/work/gwis/data/FIGI_EpiData/FIGI_{hrc_version}_gxeset_analysis_data_glm.rds")) %>% 
  filter(vcfid%in% esubset)

# stratified odds ratios table (redo, make sure you're replicating Kristina's results)
fit_stratified_or(data_epi = input_data, exposure = exposure, snp = "10:101476905:G:A", hrc_version = 'v2.3', covariates = covariates, dosage = F)

input_data <- input_data %>% 
  mutate(alcoholc_moderate_flip = fct_relevel(alcoholc_moderate, "nondrinker"))
fit_stratified_or(data_epi = input_data, exposure = "alcoholc_moderate_flip", snp = "10:101476905:G:A", hrc_version = 'v2.3', covariates = covariates, dosage = F)


test <- qread("/media/work/tmp_duncan/alcoholc_moderate/dosage_chr10_101353285.qs")




# ----- alcoholc_heavy_vs_moderate ------ #

# input variables
exposure = 'alcoholc_heavy_vs_moderate'
hrc_version = 'v2.3'
annotation_file <- 'gwas_200_ld_annotation_feb2021.txt'
path = glue("/media/work/gwis_test/{exposure}/")
covariates <- sort(c('age_ref_imp', 'sex', 'energytot_imp', 'study_gxe', 'pc1', 'pc2', 'pc3'))

# input data
esubset <- readRDS(glue("{path}/data/FIGI_{hrc_version}_gxeset_{exposure}_basic_covars_glm.rds")) %>% pull(vcfid)
input_data <- readRDS(glue("/media/work/gwis/data/FIGI_EpiData/FIGI_{hrc_version}_gxeset_analysis_data_glm.rds")) %>% 
  filter(vcfid%in% esubset)

# stratified odds ratios table (redo, make sure you're replicating Kristina's results)
fit_stratified_or(data_epi = input_data, exposure = exposure, snp = "10:101476905:G:A", hrc_version = 'v2.3', covariates = covariates, dosage = F)

test <- qread("/media/work/tmp_duncan/alcoholc_heavy_vs_moderate/dosage_chr10_101353285.qs")








# ----- alcoholc heavy vs nondrinkers ----- #

# input variables
exposure = 'alcoholc_heavy'
hrc_version = 'v2.3'
annotation_file <- 'gwas_200_ld_annotation_feb2021.txt'
path = glue("/media/work/gwis_test/{exposure}/")
covariates <- sort(c('age_ref_imp', 'sex', 'energytot_imp', 'study_gxe', 'pc1', 'pc2', 'pc3'))


# let's use moderate and heavy to obtain VCFIDs, this is because I want to ensure i have the same study inclusion/exclusion as Kristina's paper
alcoholc_moderate_vcfid <- readRDS("/media/work/gwis_test/alcoholc_moderate/data/FIGI_v2.3_gxeset_alcoholc_moderate_basic_covars_glm.rds") %>% 
  pull(vcfid)
alcoholc_heavy_vs_moderate_vcfid <- readRDS("/media/work/gwis_test/alcoholc_heavy_vs_moderate/data/FIGI_v2.3_gxeset_alcoholc_heavy_vs_moderate_basic_covars_glm.rds") %>% 
  pull(vcfid)

# N = 74099
keep <- unique(c(alcoholc_moderate_vcfid, alcoholc_heavy_vs_moderate_vcfid))

input_data <-  readRDS(glue("/media/work/gwis_test/data/FIGI_v2.3_gxeset_analysis_data_glm.rds")) %>% 
  filter(vcfid %in% keep) %>% 
  filter(!is.na(alcoholc_heavy)) %>% 
  dplyr::select(vcfid, alcoholc_heavy)
# output, send to END so you can dosages for this subset
saveRDS(input_data, file = glue("/media/work/gwis_test/alcoholc_heavy/data/FIGI_{hrc_version}_gxeset_{exposure}_basic_covars_glm.rds"))


# input data
esubset <- readRDS(glue("{path}/data/FIGI_{hrc_version}_gxeset_{exposure}_basic_covars_glm.rds")) %>% pull(vcfid)
input_data <- readRDS(glue("/media/work/gwis/data/FIGI_EpiData/FIGI_{hrc_version}_gxeset_analysis_data_glm.rds")) %>% 
  filter(vcfid%in% esubset)


fit_stratified_or(data_epi = input_data, exposure = exposure, snp = "10:101476905:G:A", hrc_version = 'v2.3', covariates = covariates, dosage = F)




test <- qread("/media/work/tmp_duncan/alcoholc_heavy/dosage_chr10_101353285.qs")









# ============================================ #
# test something - fitting both alcohol genotype in same model.. 
# do this for duncan and kristina as welll 
# ============================================= #


# input variables
exposure = 'alcoholc'
hrc_version = 'v2.3'
annotation_file <- 'gwas_200_ld_annotation_feb2021.txt'
path = glue("/media/work/gwis_test/{exposure}/")
covariates <- sort(c('age_ref_imp', 'sex', 'energytot_imp', 'study_gxe', 'pc1', 'pc2', 'pc3'))

# let's use moderate and heavy to obtain VCFIDs, this is because I want to ensure i have the same study inclusion/exclusion as Kristina's paper
alcoholc_moderate_vcfid <- readRDS("/media/work/gwis_test/alcoholc_moderate/data/FIGI_v2.3_gxeset_alcoholc_moderate_basic_covars_glm.rds") %>% 
  pull(vcfid)
alcoholc_heavy_vs_moderate_vcfid <- readRDS("/media/work/gwis_test/alcoholc_heavy_vs_moderate/data/FIGI_v2.3_gxeset_alcoholc_heavy_vs_moderate_basic_covars_glm.rds") %>% 
  pull(vcfid)

# N = 74099
keep <- unique(c(alcoholc_moderate_vcfid, alcoholc_heavy_vs_moderate_vcfid))

input_data <-  readRDS(glue("/media/work/gwis_test/data/FIGI_v2.3_gxeset_analysis_data_glm.rds")) %>% 
  filter(vcfid %in% keep)




# for rs2300985 chr10_101476905_G_A
tmp2 <- qread(glue("/media/work/gwis_test/{exposure}/output/posthoc/dosage_chr10_101476905.qs")) %>% 
  inner_join(input_data, 'vcfid') %>% 
  mutate(p0 = round(chr10_101476905_G_A_p0,0),
         p1 = round(chr10_101476905_G_A_p1,0),
         p2 = round(chr10_101476905_G_A_p2,0))


mod <- 'age_ref_imp+sex+energytot_imp+pc1+pc2+pc3+study_gxe'
mModel <- paste0('outcome~',exposure,'*chr10_101476905_G_A_p1+',exposure,'*chr10_101476905_G_A_p2+', mod)
model <- glm(mModel, data = tmp2, family = 'binomial')
       


# model with all levels of alcoholc, nondrinkers as reference 
# tmp2 <- qread(glue("/media/work/gwis_test/{exposure}/output/posthoc/dosage_chr10_101476905.qs")) %>% 
#   inner_join(input_data, 'vcfid') %>% 
#   mutate(alcoholc = fct_relevel(alcoholc, "nondrinker"))
# 
# mModel <- paste0('outcome~alcoholc*chr10_101476905_G_A_p1 + alcoholc*chr10_101476905_G_A_p2+', mod)
# model <- glm(mModel, data = tmp2, family = 'binomial')


# now just need to obtain different estimates for duncan. keep the layout of the tables the same. 

tmp <- summary(model)
COV <- tmp$cov.unscaled


# stratified by GE
# e1 <- 'alcoholc1-28g/d'
e1 <- 'alcoholcnondrinker'
e2 <- 'alcoholc>28g/d'
p1 <- 'chr10_101476905_G_A_p1'
p2 <- 'chr10_101476905_G_A_p2'


e1g0 <- tmp$coef[e1,c(1,2)] # est + SE
e2g0 <- tmp$coef[e2,c(1,2)]

beta.e1g1 <- tmp$coef[e1,1]+
  tmp$coef[p1,1]+
  tmp$coef[glue(e1, ':', p1),1]
beta.e2g1 <- tmp$coef[e2,1]+
  tmp$coef[p1,1]+
  tmp$coef[glue(e2, ':', p1),1]

beta.e1g2 <- tmp$coef[e1,1]+
  tmp$coef[p2,1]+
  tmp$coef[glue(e1, ':', p2),1]
beta.e2g2 <- tmp$coef[e2,1]+
  tmp$coef[p2,1]+
  tmp$coef[glue(e2, ':', p2),1]

e0g <- tmp$coef[c(p1, p2), c(1,2)] # getting both e0g1 and e0g2 estimates (matrix)

# calculate standard errors
I3 =  rep(1,3)
# se.eg1 = se.eg2 = NULL
# for(i in 1:(nlevels(data[,env])-1))
# {
#   covs = COV[c(env.idx[i],'p1',paste0(env.idx[i],':p1')),c(env.idx[i],'p1',paste0(env.idx[i],':p1'))]		
#   se.eg1 = c(se.eg1,sqrt(t(I3)%*%covs%*%I3))
#   covs = COV[c(env.idx[i],'p2',paste0(env.idx[i],':p2')),c(env.idx[i],'p2',paste0(env.idx[i],':p2'))]		
#   se.eg2 = c(se.eg2,sqrt(t(I3)%*%covs%*%I3))
# }
# GE<- c(eg0,e0g,beta.eg1,se.eg1,beta.eg2,se.eg2)


covs = COV[c(e1, p1 , glue(e1, ':', p1)), c(e1, p1, glue(e1, ':', p1))]
se.e1g1 = sqrt(t(I3)%*%covs%*%I3)

covs = COV[c(e1, p2 , glue(e1, ':', p2)), c(e1, p2, glue(e1, ':', p2))]
se.e1g2 = sqrt(t(I3)%*%covs%*%I3)

covs = COV[c(e2, p1 , glue(e2, ':', p1)), c(e2, p1, glue(e2, ':', p1))]
se.e2g1 = sqrt(t(I3)%*%covs%*%I3)

covs = COV[c(e2, p2 , glue(e2, ':', p2)), c(e2, p2, glue(e2, ':', p2))]
se.e2g2 = sqrt(t(I3)%*%covs%*%I3)



# calculate ODDS ratios
OR_e1g0 = glue("{rnd2(exp(e1g0[1]))} ({rnd2(exp(e1g0[1] - qnorm(0.975)*e1g0[2]))} - {rnd2(exp(e1g0[1] + qnorm(0.975)*e1g0[2]))})")
OR_e2g0 = glue("{rnd2(exp(e2g0[1]))} ({rnd2(exp(e2g0[1] - qnorm(0.975)*e2g0[2]))} - {rnd2(exp(e2g0[1] + qnorm(0.975)*e2g0[2]))})")
OR_e0g1 = glue("{rnd2(exp(e0g[1,1]))} ({rnd2(exp(e0g[1,1] - qnorm(0.975)*e0g[1,2]))} - {rnd2(exp(e0g[1,1] + qnorm(0.975)*e0g[1,2]))})")
OR_e0g2 = glue("{rnd2(exp(e0g[2,1]))} ({rnd2(exp(e0g[2,1] - qnorm(0.975)*e0g[2,2]))} - {rnd2(exp(e0g[2,1] + qnorm(0.975)*e0g[2,2]))})")

OR_e1g1 = glue("{rnd2(exp(beta.e1g1))} ({rnd2(exp(beta.e1g1 - qnorm(0.975)*se.e1g1))} - {rnd2(exp(beta.e1g1 + qnorm(0.975)*se.e1g1))})")
OR_e2g1 = glue("{rnd2(exp(beta.e2g1))} ({rnd2(exp(beta.e2g1 - qnorm(0.975)*se.e2g1))} - {rnd2(exp(beta.e2g1 + qnorm(0.975)*se.e2g1))})")
OR_e1g2 = glue("{rnd2(exp(beta.e1g2))} ({rnd2(exp(beta.e1g2 - qnorm(0.975)*se.e1g2))} - {rnd2(exp(beta.e1g2 + qnorm(0.975)*se.e1g2))})")
OR_e2g2 = glue("{rnd2(exp(beta.e2g2))} ({rnd2(exp(beta.e2g2 - qnorm(0.975)*se.e2g2))} - {rnd2(exp(beta.e2g2 + qnorm(0.975)*se.e2g2))})")





OR_e0g1
OR_e0g2

OR_e1g0
OR_e1g1
OR_e1g2

OR_e2g0
OR_e2g1
OR_e2g2





# E stratified by G
beta.e1g1 <- tmp$coef[e1,1]+
  tmp$coef[glue(e1, ':', p1),1]

beta.e2g1 <- tmp$coef[e2,1]+
  tmp$coef[glue(e2, ':', p1),1]

beta.e1g2 <- tmp$coef[e1,1]+
  tmp$coef[glue(e1, ':', p2),1]

beta.e2g2 <- tmp$coef[e2,1]+
  tmp$coef[glue(e2, ':', p2),1]



I3 = rep(1,2)

covs = COV[c(e1, glue(e1, ':', p1)), c(e1, glue(e1, ':', p1))]
se.e1g1 = sqrt(t(I3)%*%covs%*%I3)

covs = COV[c(e2, glue(e2, ':', p1)), c(e2, glue(e2, ':', p1))]
se.e2g1 = sqrt(t(I3)%*%covs%*%I3)

covs = COV[c(e1, glue(e1, ':', p2)), c(e1, glue(e1, ':', p2))]
se.e1g2 = sqrt(t(I3)%*%covs%*%I3)

covs = COV[c(e2, glue(e2, ':', p2)), c(e2, glue(e2, ':', p2))]
se.e2g2 = sqrt(t(I3)%*%covs%*%I3)
  
  

OR_e1g1 = glue("{rnd2(exp(beta.e1g1))} ({rnd2(exp(beta.e1g1 - qnorm(0.975)*se.e1g1))} - {rnd2(exp(beta.e1g1 + qnorm(0.975)*se.e1g1))})")
OR_e2g1 = glue("{rnd2(exp(beta.e2g1))} ({rnd2(exp(beta.e2g1 - qnorm(0.975)*se.e2g1))} - {rnd2(exp(beta.e2g1 + qnorm(0.975)*se.e2g1))})")
OR_e1g2 = glue("{rnd2(exp(beta.e1g2))} ({rnd2(exp(beta.e1g2 - qnorm(0.975)*se.e1g2))} - {rnd2(exp(beta.e1g2 + qnorm(0.975)*se.e1g2))})")
OR_e2g2 = glue("{rnd2(exp(beta.e2g2))} ({rnd2(exp(beta.e2g2 - qnorm(0.975)*se.e2g2))} - {rnd2(exp(beta.e2g2 + qnorm(0.975)*se.e2g2))})")

  
OR_e1g1
OR_e2g1
OR_e1g2
OR_e2g2


# G stratified by E

beta.e1g1 <- tmp$coef[p1,1]+
  tmp$coef[glue(e1, ':', p1),1]

beta.e1g2 <- tmp$coef[p2,1]+
  tmp$coef[glue(e1, ':', p2),1]

beta.e2g1 <- tmp$coef[p1,1]+
  tmp$coef[glue(e2, ':', p1),1]

beta.e2g2 <- tmp$coef[p2,1]+
  tmp$coef[glue(e2, ':', p2),1]


I3 = rep(1,2)

covs = COV[c(p1, glue(e1, ':', p1)), c(p1, glue(e1, ':', p1))]
se.e1g1 = sqrt(t(I3)%*%covs%*%I3)

covs = COV[c(p2, glue(e1, ':', p2)), c(p2, glue(e1, ':', p2))]
se.e1g2 = sqrt(t(I3)%*%covs%*%I3)

covs = COV[c(p1, glue(e2, ':', p1)), c(p1, glue(e2, ':', p1))]
se.e2g1 = sqrt(t(I3)%*%covs%*%I3)

covs = COV[c(p2, glue(e2, ':', p2)), c(p2, glue(e2, ':', p2))]
se.e2g2 = sqrt(t(I3)%*%covs%*%I3)
  
  
OR_e1g1 = glue("{rnd2(exp(beta.e1g1))} ({rnd2(exp(beta.e1g1 - qnorm(0.975)*se.e1g1))} - {rnd2(exp(beta.e1g1 + qnorm(0.975)*se.e1g1))})")
OR_e1g2 = glue("{rnd2(exp(beta.e1g2))} ({rnd2(exp(beta.e1g2 - qnorm(0.975)*se.e1g2))} - {rnd2(exp(beta.e1g2 + qnorm(0.975)*se.e1g2))})")
OR_e2g1 = glue("{rnd2(exp(beta.e2g1))} ({rnd2(exp(beta.e2g1 - qnorm(0.975)*se.e2g1))} - {rnd2(exp(beta.e2g1 + qnorm(0.975)*se.e2g1))})")
OR_e2g2 = glue("{rnd2(exp(beta.e2g2))} ({rnd2(exp(beta.e2g2 - qnorm(0.975)*se.e2g2))} - {rnd2(exp(beta.e2g2 + qnorm(0.975)*se.e2g2))})")

OR_e1g1
OR_e1g2
OR_e2g1
OR_e2g2


# --------------------------------------- #
# - same as above, but flip allele coding... 
# --------------------------------------- #
# 
# 
# model with all levels of alcoholc, nondrinkers as reference 
tmp2 <- qread(glue("/media/work/gwis_test/{exposure}/output/posthoc/dosage_chr10_101476905.qs")) %>% 
  inner_join(input_data, 'vcfid') %>% 
  mutate(alcoholc = fct_relevel(alcoholc, "nondrinker"), 
         chr10_101476905_A_G_p0 = chr10_101476905_G_A_p2, 
         chr10_101476905_A_G_p1 = chr10_101476905_G_A_p1, 
         chr10_101476905_A_G_p2 = chr10_101476905_G_A_p0)

mModel <- paste0('outcome~alcoholc*chr10_101476905_A_G_p1 + alcoholc*chr10_101476905_A_G_p2+', mod)
model <- glm(mModel, data = tmp2, family = 'binomial')


# now just need to obtain different estimates for duncan. keep the layout of the tables the same. 

tmp <- summary(model)
COV <- tmp$cov.unscaled
env.idx = env
if(is.factor(data[,env])) env.idx = paste0(env,levels(data[,env])[-1])


# stratified by GE

e1 <- 'alcoholc1-28g/d'
e2 <- 'alcoholc>28g/d'
p1 <- 'chr10_101476905_A_G_p1'
p2 <- 'chr10_101476905_A_G_p2'



e1g0 <- tmp$coef[e1,c(1,2)] # est + SE
e2g0 <- tmp$coef[e2,c(1,2)]

beta.e1g1 <- tmp$coef[e1,1]+
  tmp$coef[p1,1]+
  tmp$coef[glue(e1, ':', p1),1]
beta.e2g1 <- tmp$coef[e2,1]+
  tmp$coef[p1,1]+
  tmp$coef[glue(e2, ':', p1),1]

beta.e1g2 <- tmp$coef[e1,1]+
  tmp$coef[p2,1]+
  tmp$coef[glue(e1, ':', p2),1]
beta.e2g2 <- tmp$coef[e2,1]+
  tmp$coef[p2,1]+
  tmp$coef[glue(e2, ':', p2),1]

e0g <- tmp$coef[c(p1, p2), c(1,2)] # getting both e0g1 and e0g2 estimates (matrix)

# calculate standard errors
I3 =  rep(1,3)
se.eg1 = se.eg2 = NULL
for(i in 1:(nlevels(data[,env])-1))
{
  covs = COV[c(env.idx[i],'p1',paste0(env.idx[i],':p1')),c(env.idx[i],'p1',paste0(env.idx[i],':p1'))]		
  se.eg1 = c(se.eg1,sqrt(t(I3)%*%covs%*%I3))
  covs = COV[c(env.idx[i],'p2',paste0(env.idx[i],':p2')),c(env.idx[i],'p2',paste0(env.idx[i],':p2'))]		
  se.eg2 = c(se.eg2,sqrt(t(I3)%*%covs%*%I3))
}
GE<- c(eg0,e0g,beta.eg1,se.eg1,beta.eg2,se.eg2)



covs = COV[c(e1, p1 , glue(e1, ':', p1)), c(e1, p1, glue(e1, ':', p1))]
se.e1g1 = sqrt(t(I3)%*%covs%*%I3)

covs = COV[c(e1, p2 , glue(e1, ':', p2)), c(e1, p2, glue(e1, ':', p2))]
se.e1g2 = sqrt(t(I3)%*%covs%*%I3)

covs = COV[c(e2, p1 , glue(e2, ':', p1)), c(e2, p1, glue(e2, ':', p1))]
se.e2g1 = sqrt(t(I3)%*%covs%*%I3)

covs = COV[c(e2, p2 , glue(e2, ':', p2)), c(e2, p2, glue(e2, ':', p2))]
se.e2g2 = sqrt(t(I3)%*%covs%*%I3)



# calculate ODDS ratios
OR_e1g0 = glue("{rnd2(exp(e1g0[1]))} ({rnd2(exp(e1g0[1] - qnorm(0.975)*e1g0[2]))} - {rnd2(exp(e1g0[1] + qnorm(0.975)*e1g0[2]))})")
OR_e2g0 = glue("{rnd2(exp(e2g0[1]))} ({rnd2(exp(e2g0[1] - qnorm(0.975)*e2g0[2]))} - {rnd2(exp(e2g0[1] + qnorm(0.975)*e2g0[2]))})")
OR_e0g1 = glue("{rnd2(exp(e0g[1,1]))} ({rnd2(exp(e0g[1,1] - qnorm(0.975)*e0g[1,2]))} - {rnd2(exp(e0g[1,1] + qnorm(0.975)*e0g[1,2]))})")
OR_e0g2 = glue("{rnd2(exp(e0g[2,1]))} ({rnd2(exp(e0g[2,1] - qnorm(0.975)*e0g[2,2]))} - {rnd2(exp(e0g[2,1] + qnorm(0.975)*e0g[2,2]))})")

OR_e1g1 = glue("{rnd2(exp(beta.e1g1))} ({rnd2(exp(beta.e1g1 - qnorm(0.975)*se.e1g1))} - {rnd2(exp(beta.e1g1 + qnorm(0.975)*se.e1g1))})")
OR_e2g1 = glue("{rnd2(exp(beta.e2g1))} ({rnd2(exp(beta.e2g1 - qnorm(0.975)*se.e2g1))} - {rnd2(exp(beta.e2g1 + qnorm(0.975)*se.e2g1))})")
OR_e1g2 = glue("{rnd2(exp(beta.e1g2))} ({rnd2(exp(beta.e1g2 - qnorm(0.975)*se.e1g2))} - {rnd2(exp(beta.e1g2 + qnorm(0.975)*se.e1g2))})")
OR_e2g2 = glue("{rnd2(exp(beta.e2g2))} ({rnd2(exp(beta.e2g2 - qnorm(0.975)*se.e2g2))} - {rnd2(exp(beta.e2g2 + qnorm(0.975)*se.e2g2))})")

OR_e1g0
OR_e2g0

OR_e0g1
OR_e0g2
OR_e1g1
OR_e2g1
OR_e1g2
OR_e2g2


OR_e1g1
OR_e1g2

OR_e0g1
OR_e0g2

OR_e2g1
OR_e2g2







# 
# 
# 
# 
# # stratified by G 
# g0 <-tmp$coef[env.idx,1:2]
# if(length(env.idx)==1){
#   idx.g1 = c(env.idx,paste0(env.idx,':p1'))
#   idx.g2 = c(env.idx,paste0(env.idx,':p2'))
#   g1 = sum(tmp$coef[idx.g1,1],na.rm=T)
#   g2 = sum(tmp$coef[idx.g2,1],na.rm=T)
#   Is1 =  rep(1,length(idx.g1))
#   Is2 =  rep(1,length(idx.g2))
# }else{
#   est.g1 <- data.frame(tmp$coef[env.idx,1],tmp$coef[paste0(env.idx,':p1'),1])
#   est.g2 <- data.frame(tmp$coef[env.idx,1],tmp$coef[paste0(env.idx,':p2'),1])
#   g1 = rowSums(est.g1,na.rm=T)
#   g2 = rowSums(est.g2,na.rm=T)
#   Is1 =  rep(1,ncol(est.g1))
#   Is2 =  rep(1,ncol(est.g2))
# }
# 
# se.g1 = se.g2 = NULL
# for(i in 1:(nlevels(data[,env])-1))
# {
#   covs1 = COV[c(env.idx[i],paste0(env.idx[i],':p1')),c(env.idx[i],paste0(env.idx[i],':p1'))]  
#   covs2 = COV[c(env.idx[i],paste0(env.idx[i],':p2')),c(env.idx[i],paste0(env.idx[i],':p2'))]  
#   se.g1 <- c(se.g1,sqrt(t(Is1)%*%covs1%*%Is1))
#   se.g2 <- c(se.g2,sqrt(t(Is2)%*%covs2%*%Is2))
# }
# 
# G<- c(g0,g1,se.g1,g2,se.g2)
# 
# # stratified by E
# 
# e0 <-tmp$coef[c('p1','p2'),1:2]
# if(length(env.idx)==1){
#   idx.e11 = c('p1',paste0(env.idx,':p1'))
#   idx.e12 = c('p2',paste0(env.idx,':p2'))
#   e11 = sum(tmp$coef[idx.e11,1],na.rm=T)
#   e12 = sum(tmp$coef[idx.e12,1],na.rm=T)
#   Is1 =  rep(1,length(idx.e11))
#   Is2 =  rep(1,length(idx.e12))
# }else{
#   est.e11 <- data.frame(tmp$coef['p1',1],tmp$coef[paste0(env.idx,':p1'),1])
#   est.e12 <- data.frame(tmp$coef['p2',1],tmp$coef[paste0(env.idx,':p2'),1])
#   e11 = rowSums(est.e11,na.rm=T)
#   e12 = rowSums(est.e12,na.rm=T)
#   Is1 =  rep(1,ncol(est.e11))
#   Is2 =  rep(1,ncol(est.e12))
# }
# 
# se.e11 = se.e12 = NULL
# for(i in 1:(nlevels(data[,env])-1))
# {
#   covs1 = COV[c('p1',paste0(env.idx[i],':p1')),c('p1',paste0(env.idx[i],':p1'))]  
#   covs2 = COV[c('p2',paste0(env.idx[i],':p2')),c('p2',paste0(env.idx[i],':p2'))]  
#   se.e11 <- c(se.e11,sqrt(t(Is1)%*%covs1%*%Is1))
#   se.e12 <- c(se.e12,sqrt(t(Is2)%*%covs2%*%Is2))
# }
# E <- c(e0[1,1],e11,e0[1,2],se.e11,e0[2,1],e12,e0[2,2],se.e12)
# }
# res<-list(GE=GE,G=G,E=E)
# }
# 




# =========================================================================== #
# more analysis for kristina
# =========================================================================== #


# pooled analysis with all alcoholc categories



model <- glm(glue("outcome ~ alcoholc:chr10_101476905_G_A_dose + {glue_collapse(covariates, sep = '+')}"), data = tmp2, family = 'binomial')
summary(model)




