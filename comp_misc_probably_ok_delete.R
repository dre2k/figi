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
hrc_version = 'v3.0'





input_data <- readRDS(glue("/media/work/gwis_test/data/FIGI_{hrc_version}_GWAS.rds")) %>% 
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

pca_in <- gxe %>%
  dplyr::select(asp_ref, smk_ever, alcoholc_heavy, alcoholc_moderate, calcium_totqc2, redmeatqc2, procmeatqc2, fruitqc2, vegetableqc2, heightcm, bmi10)
pca_incor <- cor(pca_in)

out2 <- princomp(pca_in, cor = T)
pcE <- data.frame(out2$scores)

tmp1 <- gxe %>%
  dplyr::select(vcfid, asp_ref, smk_ever, alcoholc_heavy, alcoholc_moderate, calcium_totqc2, redmeatqc2, procmeatqc2, fruitqc2, vegetableqc2, heightcm, bmi10)

gxe_out <- inner_join(input_data[, c("vcfid", "outcome", "sex", "age_ref_imp", "study_gxe", "pc1", "pc2", "pc3")], tmp1, by = 'vcfid') %>% 
  bind_cols(pcE)


# is there significant interaction between PC6 and rs13086367
tmp <- qread("/media/work/gwis_test/gwas/output/dosage_chr3_112903888.qs")
gxe_out_geno <- inner_join(gxe_out, tmp, 'vcfid')
glm_out <- tidy(glm(outcome ~ chr3_112903888_A_G_dose*Comp.6 + age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe, data = gxe_out_geno, family='binomial'))



gxescan <- gxe_out_geno %>% 
  dplyr::select()

ref_study <- as.character(unique(gxe_out_geno$study_gxe)[1])
for(t in unique(gxe_out_geno$study_gxe)) {
  gxe_out_geno[paste0(t)] <- ifelse(gxe_out_geno$study_gxe == t, 1, 0)
}

gxe_out_geno <- dplyr::select(gxe_out_geno, -ref_study, -study_gxe, -Comp.6, Comp.6)

format_data_gxescan <- function(d, exposure) {
  
  tmp <- d
  ref_study <- as.character(unique(d[, 'study_gxe'])[1])
  
  for(t in unique(tmp$study_gxe)) {
    tmp[paste0(t)] <- ifelse(tmp$study_gxe==t,1,0)
  }
  
  # # tmp <- dplyr::select(tmp, -ref_study, -study_gxe, -exposure, exposure)
  tmp <- dplyr::select(tmp, -ref_study, -study_gxe, -exposure, exposure)
  
}





