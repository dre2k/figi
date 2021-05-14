#=============================================================================#
# FIGI GxE
# PCA analysis of multi E (jim)
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
library(purrr)
rm(list = ls())


hrc_version = 'v2.3'
annotation_file <- 'gwas_200_ld_annotation_feb2021.txt'
covariates <- sort(c('age_ref_imp', 'sex', 'study_gxe', 'pc1', 'pc2', 'pc3'))
input_data <- readRDS(glue("/media/work/gwis_test/data/FIGI_{hrc_version}_gxeset_analysis_data_glm.rds")) 

# filtered gxe == 1, EUR_subset == 1
input_data <- readRDS(glue("/media/work/gwis_test/data/FIGI_v2.3_GWAS.rds")) %>% 
  filter(gxe == 1) %>% 
  filter(EUR_subset == 1) %>% 
  mutate(outcome = ifelse(outcome == "Case", 1, ifelse(outcome == "Control", 0, NA)))
names(input_data) <- tolower(names(input_data))



# set alcoholc NA to 0
gxe <- input_data %>% 
  dplyr::filter(! (is.na(alcoholc_heavy) & is.na(alcoholc_moderate))) %>% 
  # alcoholc_heavy levels = nondrinker, >28g/d
  # alcoholc_moderate levels = nondrink, 1-28g/d
  mutate(alcoholc_heavy = ifelse(is.na(alcoholc_heavy), 1, alcoholc_heavy), 
         alcoholc_moderate = ifelse(is.na(alcoholc_moderate), 2, alcoholc_moderate),
         asp_ref = as.numeric(asp_ref), 
         smk_ever = as.numeric(smk_ever), 
         calcium_totqc2 = as.numeric(calcium_totqc2), 
         fruitqc2 = as.numeric(fruitqc2),  
         vegetableqc2 = as.numeric(vegetableqc2),
         heightcm = as.numeric(heightcm), 
         bmi10 = as.numeric(bmi10),
         redmeatqc2 = as.numeric(redmeatqc2), 
         procmeatqc2 = as.numeric(procmeatqc2)) %>% 
  dplyr::select(vcfid, outcome, age_ref_imp, sex, study_gxe, asp_ref, smk_ever, alcoholc_heavy, alcoholc_moderate, calcium_totqc2, redmeatqc2, procmeatqc2, fruitqc2, vegetableqc2, heightcm, bmi10) %>% 
  filter(complete.cases(.))

# tabulations
table(gxe$alcoholc_heavy)
table(gxe$alcoholc_moderate)
table(gxe$alcoholc_heavy, gxe$alcoholc_moderate, useNA = 'ifany')
table(gxe$asp_ref)
table(gxe$calcium_totqc2)


# verify directions of association (not changing directions for PCA)
glm(outcome ~ asp_ref + sex + age_ref_imp + study_gxe, data = gxe, family = 'binomial')
glm(outcome ~ smk_ever + sex + age_ref_imp + study_gxe, data = gxe, family = 'binomial')
glm(outcome ~ alcoholc_heavy + sex + age_ref_imp + study_gxe, data = gxe, family = 'binomial')
glm(outcome ~ alcoholc_moderate + sex + age_ref_imp + study_gxe, data = gxe, family = 'binomial')
# glm(outcome ~ alcoholc_heavy_vs_moderate + sex + age_ref_imp + study_gxe, data = gxe, family = 'binomial')
glm(outcome ~ calcium_totqc2 + sex + age_ref_imp + study_gxe, data = gxe, family = 'binomial')
glm(outcome ~ redmeatqc2 + sex + age_ref_imp + study_gxe, data = gxe, family = 'binomial')
glm(outcome ~ procmeatqc2 + sex + age_ref_imp + study_gxe, data = gxe, family = 'binomial')
glm(outcome ~ fruitqc2 + sex + age_ref_imp + study_gxe, data = gxe, family = 'binomial')
glm(outcome ~ vegetableqc2 + sex + age_ref_imp + study_gxe, data = gxe, family = 'binomial')
glm(outcome ~ heightcm + sex + age_ref_imp + study_gxe, data = gxe, family = 'binomial')
glm(outcome ~ bmi10 + sex + age_ref_imp + study_gxe, data = gxe, family = 'binomial')


pca_in <- gxe %>%
  dplyr::select(asp_ref, smk_ever, alcoholc_heavy, alcoholc_moderate, calcium_totqc2, redmeatqc2, procmeatqc2, fruitqc2, vegetableqc2, heightcm, bmi10)

pca_incor <- cor(pca_in)



# PCA
eigen(pca_incor)$values
eigen(pca_incor)$vectors
prcomp(pca_in, scale = T)$sd^2
prcomp(pca_in, scale = T)
out <- prcomp(pca_in, scale = T)
summary(out)

princomp(pca_in, cor = T)
princomp(pca_in, cor = T)$sd^2

out2 <- princomp(pca_in, cor = T)
summary(out2)
out2$sdev
out2$loadings
out2$center
out2$scale
out2$scores

pcE <- data.frame(out2$scores)

# regression analysis

tmp1 <- gxe %>%
  dplyr::select(vcfid, asp_ref, smk_ever, alcoholc_heavy, alcoholc_moderate, calcium_totqc2, redmeatqc2, procmeatqc2, fruitqc2, vegetableqc2, heightcm, bmi10)

gxe_out <- inner_join(input_data[, c("vcfid", "outcome", "sex", "age_ref_imp", "study_gxe", "pc1", "pc2", "pc3")], tmp1, by = 'vcfid') %>% 
  bind_cols(pcE)

gxe_out <- inner_join(input_data[, c("vcfid", "outcome", "sex", "age_ref_imp", "studyname", "pc1", "pc2", "pc3")], tmp1, by = 'vcfid') %>% 
  bind_cols(pcE)


varlist <- c("Comp.1", "Comp.2", "Comp.3", "Comp.4", 
             "Comp.5", "Comp.6", "Comp.7", "Comp.8", "Comp.9", "Comp.10", 
             "Comp.11")

run_glm <- function(x) {
  tidy(glm(glue("outcome ~ {x} + age_ref_imp + sex + pc1 + pc2 + pc3 + studyname"), data = gxe_out, family = 'binomial'))
}

results <- map_dfr(varlist, ~ run_glm(.x))
results_out <- filter(results, grepl("Comp", term))

write.table(results_out , file = "~/Dropbox/polyE_PCA_maineffects.txt", quote = F, row.names = F)




# is there significant interaction between PC6 and rs13086367
tmp <- qread("/media/work/gwis_test/gwas/output/dosage_chr3_112903888.qs")
gxe_out_geno <- inner_join(gxe_out, tmp, 'vcfid')

glm_out <- tidy(glm(outcome ~ chr3_112903888_A_G_dose*Comp.6 + age_ref_imp + sex + pc1 + pc2 + pc3 + studyname, data = gxe_out_geno, family='binomial'))
