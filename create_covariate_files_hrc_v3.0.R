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
output_dir <- "/media/work/gwis/data/"
setwd(glue("{output_dir}/FIGI_samplefile_epi-201014/"))


# figi_gwas <- readRDS("~/data/FIGI_EpiData_rdata/FIGI_v2.4_GWAS.rds")
# figi_gwas <- readRDS("~/data/FIGI_EpiData_rdata/FIGI_v2.4_GWAS.rds")
figi_gwas <- readRDS(glue("{output_dir}/FIGI_EpiData/FIGI_v3.0_GWAS.rds")) 
figi_gxe <- dplyr::filter(figi_gwas, gxe == 1, EUR_subset == 1)


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
  saveRDS(cov, file = glue(output_dir, "FIGI_EpiData/FIGI_v3.0_gxeset_{exposure}_basic_covars_glm.rds"), version = 2)
  
  cov_gxescan <- format_data_gxescan(cov, exposure)
  # saveRDS(cov_gxescan, file = glue("/media/work/gwis/results/input/FIGI_v3.0_gxeset_{exposure}_basic_covars_gxescan.rds"), version = 2)
  saveRDS(cov_gxescan, file = glue(output_dir, "FIGI_EpiData/FIGI_v3.0_gxeset_{exposure}_basic_covars_gxescan.rds"), version = 2)
  
}





# ------- aspirin ------- # 
# i want to check the numbers HRC v3 vs v2.4 to make sure the variables aren't different
exclude <- c()
wrap(figi_gxe, 'aspirin', studies_to_exclude = exclude)
x <- readRDS("~/data/gwis_v2.4/input/FIGI_v2.4_gxeset_aspirin_basic_covars_glm.rds")
y <- readRDS("/media/work/gwis/results/input/FIGI_v3.0_gxeset_aspirin_basic_covars_glm.rds")
table(x$study_gxe, x$outcome)
table(y$study_gxe, y$outcome)

# look OK, few added cases where the intersection showed. 




# ------ Calcium ------
# 
# calcium_totqc2
# recreate this variable now that it has been reimputed. 
# imputation was simply carried out, no checks performed yet (maybe it'll be an issue for reviewers for not sure yet)

exclude = c("ATBC", "PPS3", "PPS4")
wrap(figi_gwas, 'calcium_totqc2', vars_to_exclude = c(""), studies_to_exclude = exclude)
x <- readRDS(glue(output_dir, "FIGI_EpiData/FIGI_v3.0_gxeset_calcium_totqc2_basic_covars_glm.rds"))
table(x$study_gxe, x$outcome)

ccfr <- filter(figi_gwas, 
               !(is.na(calcium_totqc2)), 
               EUR_subset == 1,
               grepl("CCFR", study_gxe))
table(ccfr$outcome, ccfr$calcium_totqc2)

reach <- filter(figi_gwas, 
                !(is.na(calcium_totqc2)), 
                EUR_subset == 1,
                grepl("REACH", study_gxe))
table(reach$outcome, reach$calcium_totqc2)

sms_ad <- filter(figi_gwas, 
                 !(is.na(calcium_totqc2)), 
                 EUR_subset == 1,
                 grepl("SMS_AD", study_gxe))
table(sms_ad$outcome, sms_ad$calcium_totqc2)

ukb <- filter(figi_gwas, 
              !(is.na(calcium_totqc2)), 
              EUR_subset == 1,
              grepl("UKB", study_gxe))
table(ukb$outcome, ukb$calcium_totqc2)


# to add the above studies, let me just fix it here directly with rbinds
tmp <- figi_gxe %>% 
  dplyr::filter(!is.na(calcium_totqc2), 
                EUR_subset == 1,
                study_gxe %in% c("CCFR_1", "CCFR_4", "CCFR_3", "REACH_AD", "UKB_1")) %>% 
  dplyr::mutate(outcome = ifelse(outcome == "Control", 0, 1),
                sex = ifelse(sex == "Female", 0, 1),
                calcium_totqc2 = as.numeric(calcium_totqc2) - 1 ) %>% 
  dplyr::select(vcfid, outcome, calcium_totqc2, age_ref_imp, sex, energytot_imp, study_gxe, PC1, PC2, PC3)

out <- rbind(x, tmp)
saveRDS(out, glue(output_dir, "FIGI_EpiData/FIGI_v3.0_gxeset_calcium_totqc2_basic_covars_glm.rds"))

cov_gxescan <- format_data_gxescan(out, 'calcium_totqc2')
saveRDS(cov_gxescan, file = paste0("/media/work/gwis/results/input/FIGI_v3.0_gxeset_", "calcium_totqc2", "_basic_covars_gxescan.rds"), version = 2)

y <- readRDS(glue(output_dir, "FIGI_EpiData/FIGI_v3.0_gxeset_calcium_totqc2_basic_covars_glm.rds"))
table(y$study_gxe, y$outcome)









# ------ Calcium ------
# calcium_dietqc2 (we will be updating results to version HRC v3.0 to match total calcium variable)

exclude = c("ATBC", "PPS3", "PPS4")
wrap(figi_gwas, 'calcium_dietqc2', vars_to_exclude = c(""), studies_to_exclude = exclude)
x <- readRDS(glue(output_dir, "FIGI_EpiData/FIGI_v3.0_gxeset_calcium_dietqc2_basic_covars_glm.rds"))
table(x$study_gxe, x$outcome)
table(x$study_gxe, x$calcium_dietqc2)


cov_gxescan <- format_data_gxescan(x, 'calcium_dietqc2')
saveRDS(cov_gxescan, file = paste0("/media/work/gwis/results/input/FIGI_v3.0_gxeset_", "calcium_dietqc2", "_basic_covars_gxescan.rds"), version = 2)
saveRDS(cov_gxescan, file = paste0(output_dir, "FIGI_EpiData/FIGI_v3.0_gxeset_", "calcium_dietqc2", "_basic_covars_gxescan.rds"), version = 2)


y <- readRDS(glue(output_dir, "FIGI_EpiData/FIGI_v3.0_gxeset_calcium_dietqc2_basic_covars_glm.rds"))
table(y$study_gxe, y$outcome)





# ------ Folate ------

# folate_totqc2
exclude = c("UKB_1", "SMS_AD")
wrap(figi_gwas, 'folate_totqc2', vars_to_exclude = c(""), vars_to_include = c(""), studies_to_exclude = exclude)
x <- readRDS(glue(output_dir, "FIGI_EpiData/FIGI_v3.0_gxeset_folate_totqc2_basic_covars_glm.rds"))
table(x$study_gxe, x$outcome)

wtf <- filter(figi_gwas, study_gxe == "ASTERISK")
wtf <- filter(figi_gwas, grepl("DACH", study_gxe))
wtf <- filter(figi_gwas, grepl("PHS", study_gxe)) ## **
wtf <- filter(figi_gwas, grepl("REACH_AD", study_gxe)) ## **
wtf <- filter(figi_gwas, grepl("SMS", study_gxe)) ## **
table(wtf$folate_totqc2)

tmp <- figi_gxe %>% 
  dplyr::filter(!is.na(folate_totqc2), 
                EUR_subset == 1,
                study_gxe %in% c("PHS", "REACH_AD")) %>% 
  dplyr::mutate(outcome = ifelse(outcome == "Control", 0, 1),
                sex = ifelse(sex == "Female", 0, 1),
                folate_totqc2 = as.numeric(folate_totqc2) - 1 ) %>% 
  dplyr::select(vcfid, outcome, folate_totqc2, age_ref_imp, sex, energytot_imp, study_gxe, PC1, PC2, PC3)

out <- rbind(x, tmp)
saveRDS(out, glue(output_dir, "FIGI_EpiData/FIGI_v3.0_gxeset_folate_totqc2_basic_covars_glm.rds"))

cov_gxescan <- format_data_gxescan(out, 'folate_totqc2')
saveRDS(cov_gxescan, file = paste0("/media/work/gwis/results/input/FIGI_v3.0_gxeset_", "folate_totqc2", "_basic_covars_gxescan.rds"), version = 2)

y <- readRDS(glue(output_dir, "FIGI_EpiData/FIGI_v3.0_gxeset_folate_totqc2_basic_covars_glm.rds"))
table(y$study_gxe, y$outcome)
table(y$study_gxe, y$folate_totqc2)









# folate_dietqc2
exclude = c("UKB_1", "SMS_AD")
wrap(figi_gwas, 'folate_dietqc2', vars_to_exclude = c(""), vars_to_include = c(""), studies_to_exclude = exclude)
x <- readRDS(glue(output_dir, "FIGI_EpiData/FIGI_v3.0_gxeset_folate_dietqc2_basic_covars_glm.rds"))
# x <- readRDS(glue(output_dir, "FIGI_EpiData/FIGI_v3.0_gxeset_folate_dietqc2_basic_covars_gxescan.rds"))

table(x$study_gxe, x$outcome)
table(x$study_gxe, x$folate_dietqc2)





saveRDS(x, glue(output_dir, "FIGI_EpiData/FIGI_v3.0_gxeset_folate_dietqc2_basic_covars_glm.rds"))

cov_gxescan <- format_data_gxescan(x, 'folate_dietqc2')
saveRDS(cov_gxescan, file = paste0("/media/work/gwis/results/input/FIGI_v3.0_gxeset_", "folate_dietqc2", "_basic_covars_gxescan.rds"), version = 2)

y <- readRDS(glue(output_dir, "FIGI_EpiData/FIGI_v3.0_gxeset_folate_dietqc2_basic_covars_glm.rds"))
table(y$study_gxe, y$outcome)
table(y$study_gxe, y$folate_totqc2)














