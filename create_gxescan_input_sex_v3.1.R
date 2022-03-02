#=============================================================================#
# FIGI Analysis 12/12/2019
# FIGI Analysis 10/14/2020
#
# Create covariate tables for GxEScanR
# 
# ------ Notes ------ #
# use gxe set as determined in `gxe` variable (in addition to drop == 0)
# 	- gxe == 1 removes nonEUR studies
#		- c("HispanicCCS", "Taiwan", "SMHS", "SWHS", "HCES-CRC", "PURIFICAR")
# 	- gxe == 1 removes outc == "other"!
# study_gxe = adj. for study + platform
# remove case-only or control-only studies (causes gxescanr to crash)
#=============================================================================#
library(tidyverse)
library(data.table)
library(figifs)
library(glue)
rm(list = ls())
output_dir <- "/media/work/gwis_test/data/"

figi_gwas <- readRDS(glue("{output_dir}/FIGI_v3.1_GWAS.rds")) 
# FYI - Yi's summaries she shared with Mireia et al did not filter EUR == 1
figi_gxe <- dplyr::filter(figi_gwas, gxe == 1, EUR_subset == 1)
# figi_gxe <- dplyr::filter(figi_gwas, gxe == 1)


# other than EUR they all match ok 
# test <- filter(gxe, EUR_subset == 1)
# test2 <- inner_join(figi_gxe, test, 'vcfid')

#-----------------------------------------------------------------------------#
# create exposure phenotype files for gxescanr
#
# ------ Notes ------ #
# I check counts against old dataset. Gwas and Gxe match in numbers
# (minus the single exclusion from CCFR_4)
#
# remember that you do NOT need to have a reference study indicator variable
# when running models. But I am using a reference study just in case
# GxEScanR complains. The reference study will always be the first
# study in a e main effect dataset (alphabetical, typically ATBC)
#
# Updates:
#   - HawaiiCCS_AD N = 1264
#   - LCCS N = 2085
#   - nsaids coding changed a little bit
#   - NOTE! --- still waiting for one more study.. not sure what to do?
#
# IMPORTANT
# - this step requires some thought, make sure it conforms to analysis plan 
#   decisions
# - be mindful of covariates
#-----------------------------------------------------------------------------#
wrap <- function(d, exposure, is_e_categorical = T, min_cell_size = 0, vars_to_exclude = c("energytot_imp"), vars_to_include = c(), studies_to_exclude) {
  cov <- format_data_glm(d, exposure, is_e_categorical, min_cell_size, vars_to_exclude, vars_to_include, eur_only = T) %>% 
    filter(!study_gxe %in% studies_to_exclude) %>% 
    mutate(study_gxe = fct_drop(study_gxe))
  # saveRDS(cov, file = glue("/media/work/gwis/results/input/FIGI_v3.0_gxeset_{exposure}_basic_covars_glm.rds"), version = 2)
  saveRDS(cov, file = glue(output_dir, "/FIGI_v3.0_gxeset_{exposure}_basic_covars_glm.rds"), version = 2)
  
  cov_gxescan <- format_data_gxescan(cov, exposure)
  # saveRDS(cov_gxescan, file = glue("/media/work/gwis/results/input/FIGI_v3.0_gxeset_{exposure}_basic_covars_gxescan.rds"), version = 2)
  saveRDS(cov_gxescan, file = glue(output_dir, "/FIGI_v3.0_gxeset_{exposure}_basic_covars_gxescan.rds"), version = 2)
  
}


# ---- update from 7/27/2021 ---- #
# this is the format_data_glm / gxescan functions

# ---- functions
format_data_glm <- function(d, 
                            exposure, 
                            is_e_categorical, 
                            min_cell_size = 0, 
                            vars_to_exclude = c(), 
                            vars_to_include = c(), 
                            eur_only=T) {
  
  vars_to_keep <- c("vcfid", "outcome", exposure, "age_ref_imp", "sex", "energytot_imp", "study_gxe", "PC1", "PC2", "PC3", vars_to_include)
  vars_to_keep <- vars_to_keep[!vars_to_keep %in% vars_to_exclude]
  
  # note that in gxe set, outcome+age_ref_imp+sex+study_gxe+energytot_imp do NOT have missing values
  # OK to subset simply by using is.na(exposure)
  if(eur_only == T) {
    tmp <- d %>%
      dplyr::filter(gxe == 1,
                    EUR_subset == 1,
                    !is.na(get(exposure))) # %>%
    # dplyr::mutate(outcome = ifelse(outcome == "Control", 0, 1),
    #               sex = ifelse(sex == "Female", 0, 1))
  } else {
    tmp <- d %>%
      dplyr::filter(gxe == 1,
                    !is.na(get(exposure))) # %>%
    # dplyr::mutate(outcome = ifelse(outcome == "Control", 0, 1),
    #               sex = ifelse(sex == "Female", 0, 1))
  }
  
  
  # drop zero cells, keep vars_to_keep
  if (is_e_categorical == T) {
    tmp <- mutate(tmp, {{exposure}} := as.numeric(get(exposure)) - 1)
    
    drops <- data.frame(table(tmp$outcome, tmp[, exposure], tmp$study_gxe)) %>%
      filter(Freq <= min_cell_size)
    
    tmp <- filter(tmp, !study_gxe %in% unique(drops$Var3)) %>%
      dplyr::mutate(study_gxe = fct_drop(study_gxe)) %>%
      dplyr::select(vars_to_keep) %>%
      filter(complete.cases(.))}
  else {
    drops <- data.frame(table(tmp$outcome, tmp$study_gxe)) %>%
      filter(Freq <= min_cell_size)
    tmp <- filter(tmp, !study_gxe %in% unique(drops$Var2)) %>%
      dplyr::mutate(study_gxe = fct_drop(study_gxe)) %>%
      dplyr::select(vars_to_keep) %>%
      filter(complete.cases(.))}
  
  return(tmp)
  
}



format_data_gxescan <- function(d, exposure) {
  
  tmp <- d
  ref_study <- as.character(unique(d[, 'study_gxe'])[1])
  
  for(t in unique(tmp$study_gxe)) {
    tmp[paste0(t)] <- ifelse(tmp$study_gxe==t,1,0)
  }
  
  # # tmp <- dplyr::select(tmp, -ref_study, -study_gxe, -exposure, exposure)
  tmp <- dplyr::select(tmp, -ref_study, -study_gxe, -exposure, exposure)
  
}



# ------ methrswklns N = 42602 ------ #
# input file has methrswklns scaled by 10

pca <- fread("/media/work/gwis_test/PCA/20210222/figi_gxe_pca_update.eigenvec") %>% 
  rename(vcfid = `#FID`) %>% 
  dplyr::select(vcfid, PC1, PC2, PC3)




exclude <- c("CRCGEN", "CLUEII", "NCCCSII")
exclude <- c("")


vars_to_keep <- c("vcfid", "outcome", "age_ref_imp", "sex", "study_gxe", "PC1", "PC2", "PC3")
tmp <- figi_gxe %>%
  dplyr::filter(gxe == 1,
                EUR_subset == 1) %>% 
  mutate(outcome = fct_drop(outcome))


drops <- data.frame(table(tmp$outcome, tmp$study_gxe)) %>%
  filter(Freq <= 0)

drops_e <- data.frame(table(tmp$sex, tmp$study_gxe)) %>% 
  filter(Freq <= 0)

tmp <- tmp %>% 
  filter(!study_gxe %in% unique(drops$Var2), 
         !study_gxe %in% unique(drops_e$Var2)) %>%
  dplyr::select(all_of(vars_to_keep)) %>%
  filter(complete.cases(.)) %>% 
  dplyr::filter(!study_gxe %in% exclude) %>% 
  dplyr::mutate(study_gxe = fct_drop(study_gxe), 
                outcome = ifelse(outcome == "Control", 0, 1),
                sex_v2 = ifelse(sex == "Female", 0, 1)) %>% 
  dplyr::select(-sex)



# should i incorporate new PCs ... 
tmp <- tmp %>% 
  dplyr::select(-starts_with("PC")) %>% 
  inner_join(pca, 'vcfid')
table(tmp$study_gxe, tmp$outcome)


saveRDS(tmp, file = glue("/media/work/gwis_test/sex_v2/data/FIGI_v3.1_gxeset_sex_v2_basic_covars_glm.rds"))

cov_gxescan <- format_data_gxescan(tmp, 'sex_v2')
saveRDS(cov_gxescan, file = glue("/media/work/gwis_test/sex_v2/data/FIGI_v3.1_gxeset_sex_v2_basic_covars_gxescan.rds"), version = 2)

x <- readRDS(glue("/media/work/gwis_test/sex_v2/data/FIGI_v3.1_gxeset_sex_v2_basic_covars_glm.rds"))
table(x$study_gxe, x$outcome)
y <- readRDS(glue("/media/work/gwis_test/sex_v2/data/FIGI_v3.1_gxeset_sex_v2_basic_covars_gxescan.rds"))
head(y)







# ------ active_met_875 N = 42602 ------ #

# active_met_875
load("/media/work/gwis_test/data/physact_additional_vars/PA_Dataset.Rdata")

tmp <- figi_gxe %>%
  dplyr::select(-methrswklns) %>% 
  inner_join(gxe[,c('vcfid', 'methrswklns')], 'vcfid') %>%
  dplyr::filter(gxe == 1,
                EUR_subset == 1,
                !is.na(methrswklns)) %>% 
  mutate(outcome = fct_drop(outcome), 
         active_met_875 = ifelse(is.na(methrswklns),NA,
                                 ifelse(methrswklns>=8.75,'Yes','No')))


exclude <- c("CRCGEN", "CLUEII", "NCCCSII")
pca <- fread("/media/work/gwis_test/PCA/20210222/figi_gxe_pca_update.eigenvec") %>% 
  rename(vcfid = `#FID`) %>% 
  dplyr::select(vcfid, PC1, PC2, PC3)

vars_to_keep <- c("vcfid", "outcome", "active_met_875", "age_ref_imp", "sex", "study_gxe", "energytot_imp", "PC1", "PC2", "PC3")

drops <- data.frame(table(tmp$outcome, tmp$study_gxe)) %>%
  filter(Freq <= 0)

tmp <- filter(tmp, !study_gxe %in% unique(drops$Var2)) %>%
  dplyr::select(all_of(vars_to_keep)) %>%
  filter(complete.cases(.)) %>% 
  dplyr::filter(!study_gxe %in% exclude) %>% 
  dplyr::mutate(study_gxe = fct_drop(study_gxe), 
                active_met_875 = as.numeric(ifelse(active_met_875 == "No", 0, 1)), 
                outcome = ifelse(outcome == "Control", 0, 1),
                sex = ifelse(sex == "Female", 0, 1))


table(tmp$active_met_875, tmp$study_gxe)
drops_e <- data.frame(table(tmp$active_met_875, tmp$study_gxe)) %>% 
  filter(Freq <= 0)


tmp <- tmp %>% 
  filter(!study_gxe %in% drops_e$Var2)



# should i incorporate new PCs ... 
# (already did for other variables)
tmp <- tmp %>% 
  dplyr::select(-starts_with("PC")) %>% 
  inner_join(pca, 'vcfid')
table(tmp$study_gxe, tmp$outcome)

saveRDS(tmp, file = glue("/media/work/gwis_test/active_met_875/data/FIGI_v3.1_gxeset_active_met_875_basic_covars_glm.rds"))

cov_gxescan <- format_data_gxescan(tmp, 'active_met_875')
saveRDS(cov_gxescan, file = glue("/media/work/gwis_test/active_met_875/data/FIGI_v3.1_gxeset_active_met_875_basic_covars_gxescan.rds"), version = 2)

x <- readRDS(glue("/media/work/gwis_test/active_met_875/data/FIGI_v3.1_gxeset_active_met_875_basic_covars_glm.rds"))
table(x$study_gxe, x$outcome)
table(x$outcome, x$active_met_875)
y <- readRDS(glue("/media/work/gwis_test/active_met_875/data/FIGI_v3.1_gxeset_active_met_875_basic_covars_gxescan.rds"))
head(y)















