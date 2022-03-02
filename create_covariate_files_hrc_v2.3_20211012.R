#=============================================================================#
# FIGI Analysis 12/12/2019
# updated 2021/10/13 to create bmi5 var (just making sure things look right..)
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
rm(list = ls())

figi_gwas <- readRDS("~/data/FIGI_EpiData_rdata/FIGI_v2.3_GWAS.rds")
figi_gxe <- dplyr::filter(figi_gwas, gxe == 1, EUR_subset == 1) %>% 
  mutate(outcome = fct_drop(outcome))


#-----------------------------------------------------------------------------#
# create exposure phenotype files for gxescan
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
# looks like additional epi data from lccs and hawaiiccs, confirm with Yi
# (this is specific for NSAIDs since I know from previous runs)
#   - HawaiiCCS_AD N = 1264
#   - LCCS N = 2085
#
# IMPORANT
# - this step requires some thought, make sure it conforms to analysis plan 
#   decisions
# - be mindful of covariates
#-----------------------------------------------------------------------------#


# ---- functions
format_data_glm <- function(d, 
                            exposure, 
                            is_e_categorical, 
                            min_cell_size = 0, 
                            vars_to_exclude = c(), 
                            vars_to_include = c(), 
                            eur_only=T) {
  
  vars_to_keep <- c("vcfid", "outcome", exposure, "age_ref_imp", "sex", "energytot_imp", "study_gxe", "PC1", "PC2", "PC3")
  vars_to_keep <- vars_to_keep[!vars_to_keep %in% vars_to_exclude]
  
  # note that in gxe set, outcome+age_ref_imp+sex+study_gxe+energytot_imp do NOT have missing values
  # OK to subset simply by using is.na(exposure)
  if(eur_only == T) {
    tmp <- d %>%
      dplyr::filter(gxe == 1,
                    EUR_subset == 1,
                    !is.na(get(exposure)))  %>%
    dplyr::mutate(outcome = ifelse(outcome == "Control", 0, 1),
                  sex = ifelse(sex == "Female", 0, 1))
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
      filter(complete.cases(.))
  }
  
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





wrap <- function(d, exposure, is_e_categorical = T, min_cell_size = 0, vars_to_exclude = c("energytot_imp"), vars_to_include = c(), studies_to_exclude) {
  cov <- format_data_glm(d, exposure, is_e_categorical, min_cell_size, vars_to_exclude, vars_to_include, eur_only = T) %>% 
    filter(!study_gxe %in% studies_to_exclude) %>% 
    mutate(study_gxe = fct_drop(study_gxe))
  saveRDS(cov, file = paste0("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_", exposure, "_basic_covars_glm.rds"), version = 2)
  
  cov_gxescan <- format_data_gxescan(cov, exposure)
  saveRDS(cov_gxescan, file = paste0("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_", exposure, "_basic_covars_gxescan.rds"), version = 2)
}






# ----- BMI5 ------ 

# bmi5
exclude = c("")
wrap(figi_gxe, exposure = 'bmi5', is_e_categorical = F, vars_to_exclude = c("energytot_imp"), studies_to_exclude = exclude)
x <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_bmi5_basic_covars_glm.rds")
table(x$study_gxe, x$outcome)
table(x$study_gxe, x$bmi5)


# ----- BMI5 female -----
x <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_bmi5_basic_covars_glm.rds")
y <- readRDS("/media/work/gwis_test/bmi5/data/FIGI_v2.3_gxeset_bmi5_basic_covars_glm.rds")
all(x==y)
wtf <- setdiff(x,y)
identical(x,y)


cov <- filter(x, sex == 0) %>% 
  mutate(bmi5_female = bmi5) %>%
  dplyr::select(-sex, -bmi5)
cov_gxescan <- format_data_gxescan(cov, 'bmi5_female')
saveRDS(cov, file = paste0("/media/work/gwis_test/bmi5_female/data/FIGI_v2.3_gxeset_bmi5_female_basic_covars_glm.rds"), version = 2)
saveRDS(cov_gxescan, file = paste0("/media/work/gwis_test/bmi5_female/data/FIGI_v2.3_gxeset_bmi5_female_basic_covars_gxescan.rds"), version = 2)


# ----- BMI5 male -----
cov <- filter(x, sex == 1) %>% 
  mutate(bmi5_male = bmi5) %>%
  dplyr::select(-sex, -bmi5)
cov_gxescan <- format_data_gxescan(cov, 'bmi5_male')
saveRDS(cov, file = paste0("/media/work/gwis_test/bmi5_male/data/FIGI_v2.3_gxeset_bmi5_male_basic_covars_glm.rds"), version = 2)
saveRDS(cov_gxescan, file = paste0("/media/work/gwis_test/bmi5_male/data/FIGI_v2.3_gxeset_bmi5_male_basic_covars_gxescan.rds"), version = 2)



