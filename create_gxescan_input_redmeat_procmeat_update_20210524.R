#=============================================================================#
# Create redmeat/procmeat gxescan inputs
#=============================================================================#

# ---- explore inclusion/exclusion ---- 
library(tidyverse)
library(data.table)

hrc_version = 'v3.0'
# hrc_version = 'v2.3'
input_data <- readRDS(glue("/media/work/gwis/data/FIGI_EpiData/FIGI_{hrc_version}_gxeset_analysis_data_glm.rds"))

meatonly <- filter(input_data, !is.na(redmeatqcm ))

meatonly <- filter(input_data, !is.na(procmeatqc2 ))


table(meatonly$study)
table(meatonly$adenoma, useNA = 'ifany')

table(meatonly$study, meatonly$adenoma, useNA = 'ifany')
table(meatonly$study, meatonly$crc, useNA = 'ifany')


adenomas <- filter(meatonly, !is.na(adenoma))



# SMS
# - ALL ADENOMAS
sms <- filter(meatonly, grepl("SMS", study))
table(sms$adenoma)
table(sms$crc)
table(sms$outc, sms$adenoma, useNA = 'ifany')
table(sms$outc, sms$crc, sms$adenoma, useNA = 'ifany')


# HPFS
# - remove adenomas, final N ~ 1165
hpfs <- filter(meatonly, grepl("HPFS", study))
table(hpfs$outc, hpfs$crc, useNA = 'ifany')

table(hpfs$adenoma, useNA = 'ifany')
table(hpfs$crc, useNA = 'ifany')
table(hpfs$crc, hpfs$adenoma, useNA = 'ifany')
table(hpfs$outc, hpfs$crc, hpfs$adenoma, useNA = 'ifany')


hpfs <- filter(meatonly, grepl("HPFS", study), 
               !is.na(crc))
table(hpfs$crc)
table(hpfs$adenoma)
table(hpfs$crc, hpfs$adenoma, useNA = 'ifany')


# NHS
# - remove adenomas, final N ~ 2084
nhs <- filter(meatonly, grepl("NHS", study))
table(nhs$outc, nhs$crc, useNA = 'ifany')

table(nhs$adenoma, useNA = 'ifany')
table(nhs$crc, useNA = 'ifany')
table(nhs$crc, nhs$adenoma, useNA = 'ifany')
table(nhs$outc, nhs$crc, nhs$adenoma, useNA = 'ifany')



# UKB
# - remove adenomas, final N ~ 13644
ukb <- filter(meatonly, grepl("UKB", study))
table(ukb$outc, ukb$crc, useNA = 'ifany')

table(ukb$adenoma, useNA = 'ifany')
table(ukb$crc, useNA = 'ifany')
table(ukb$crc, ukb$adenoma, useNA = 'ifany')
table(ukb$outc, ukb$crc, ukb$adenoma, useNA = 'ifany')



# PLCO
# - remove adenomas, final N ~ 4757
plco <- filter(meatonly, grepl("PLCO", study))
table(plco$outc, plco$crc, useNA = 'ifany')

table(plco$adenoma, useNA = 'ifany')
table(plco$crc, useNA = 'ifany')
table(plco$crc, plco$adenoma, useNA = 'ifany')
table(plco$outc, plco$crc, plco$adenoma, useNA = 'ifany')


# NHS2
# - remove adenomas, final N ~ 4757
nhs2 <- filter(meatonly, grepl("NHSII", study))
table(nhs2$outc, nhs2$crc, useNA = 'ifany')

table(nhs2$adenoma, useNA = 'ifany')
table(nhs2$crc, useNA = 'ifany')
table(nhs2$crc, nhs2$adenoma, useNA = 'ifany')
table(nhs2$outc, nhs2$crc, nhs2$adenoma, useNA = 'ifany')



# hawaiiccs
# - exclude
hawaii <- filter(meatonly, grepl("Hawai", study))
table(hawaii$outc, hawaii$crc, useNA = 'ifany')

table(hawaii$adenoma, useNA = 'ifany')
table(hawaii$crc, useNA = 'ifany')
table(hawaii$crc, hawaii$adenoma, useNA = 'ifany')
table(hawaii$outc, hawaii$crc, hawaii$adenoma, useNA = 'ifany')

# reach
# - exclude
reach <- filter(meatonly, grepl("REACH", study))
table(reach$outc, reach$crc, useNA = 'ifany')

table(reach$adenoma, useNA = 'ifany')
table(reach$crc, useNA = 'ifany')
table(reach$crc, reach$adenoma, useNA = 'ifany')
table(reach$outc, reach$crc, reach$adenoma, useNA = 'ifany')






#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#

# ---- create gxescan input file ---- 
# comments:
# - to keep analysis consistent with previous runs, should use v2.3 + old PCs
# - vs .. updating analysis to latest dataset + updated PCs
# - DECISION - 


# ---- functions
format_data_glm <- function(d, 
                            exposure, 
                            is_e_categorical, 
                            min_cell_size = 0, 
                            vars_to_exclude = c(), 
                            vars_to_include = c(), 
                            eur_only=T) {
  
  vars_to_keep <- c("vcfid", "outcome", exposure, "age_ref_imp", "sex", "energytot_imp", "study_gxe", "pc1", "pc2", "pc3", vars_to_include)
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





wrap <- function(d, 
                 exposure, 
                 hrc_version,
                 is_e_categorical = T, 
                 min_cell_size = 0, 
                 vars_to_exclude = c(), 
                 vars_to_include = c(), 
                 studies_to_exclude, 
                 path) {
  
  cov <- format_data_glm(d, exposure, is_e_categorical, min_cell_size, vars_to_exclude, vars_to_include, eur_only = T) %>% 
    dplyr::filter(!study_gxe %in% studies_to_exclude) %>% 
    mutate(study_gxe = fct_drop(study_gxe))
  saveRDS(cov, file = glue(path, "/FIGI_{hrc_version}_gxeset_{exposure}_basic_covars_glm.rds"), version = 2)
  
  cov_gxescan <- format_data_gxescan(cov, exposure)
  saveRDS(cov_gxescan, file = glue(path, "FIGI_{hrc_version}_gxeset_{exposure}_basic_covars_gxescan.rds"), version = 2)
  
}






#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
# ---- redmeatqcm ----
hrc_version = 'v3.0'

# remove old PC values
# remove adenoma cases (exclude list + subsets in NHS, HPFS, PLCO, UKB)
# (best way might be to remove *_AD studies in the first place. 
# there are some left over cases that are 'not' AD, but don't want to collapse them to other studies anyway.. )
# keep if !is.na(redmeatqcm)
exclude = c("HawaiiCCS_AD", "PPS3", "PPS4", "REACH_AD", "SMS_AD", "NFCCR_2", "HPFS_3_AD", "NHS_3_AD", "NHS_5_AD", "PLCO_4_AD", "HPFS_5_AD")
input_data <- readRDS(glue("/media/work/gwis/data/FIGI_EpiData/FIGI_{hrc_version}_gxeset_analysis_data_glm.rds")) %>% 
  dplyr::select(-starts_with("pc")) %>% 
  dplyr::filter(!is.na(redmeatqcm)) %>% 
  dplyr::filter(! study_gxe %in% exclude, 
                is.na(adenoma))

table(input_data$study_gxe, input_data$redmeatqc2)
table(input_data$study_gxe, input_data$adenoma, useNA = 'ifany')
table(input_data$study_gxe, input_data$crc, useNA = 'ifany')
table(input_data$study_gxe, input_data$outc, useNA = 'ifany')

pca <- fread("/media/work/gwis_test/PCA/20210222/figi_gxe_pca_update.eigenvec") %>% 
  dplyr::rename(vcfid = IID) %>% 
  dplyr::select(-`#FID`) %>% 
  dplyr::rename_with(tolower)

# creating redmeatqcm_v2 for convenience for workflow creation
gxe <- inner_join(input_data, pca, 'vcfid') %>% 
  mutate(redmeatqcm_v2 = as.numeric(redmeatqcm))



wrap(d = gxe, exposure = 'redmeatqcm_v2', hrc_version = 'v3.0', is_e_categorical = F, path = '/media/work/gwis_test/redmeatqcm_v2/data/', studies_to_exclude = "")
check <- readRDS("/media/work/gwis_test/redmeatqcm_v2/data/FIGI_v3.0_gxeset_redmeatqcm_v2_basic_covars_glm.rds")
table(check$study_gxe, check$outcome, useNA = 'ifany')
model <- glm(outcome ~ redmeatqcm_v2 + age_ref_imp + sex + study_gxe, data = check, family = 'binomial')
summary(model)
exp(coef(model))

check <- readRDS("/media/work/gwis_test/redmeatqcm_v2/data/FIGI_v3.0_gxeset_redmeatqcm_v2_basic_covars_gxescan.rds")






#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
# ---- procmeatqcm ----

exclude = c("HawaiiCCS_AD", "PPS3", "PPS4", "REACH_AD", "SMS_AD", "NFCCR_2", "HPFS_3_AD", "NHS_3_AD", "NHS_5_AD", "PLCO_4_AD", "HPFS_5_AD")
input_data <- readRDS(glue("/media/work/gwis/data/FIGI_EpiData/FIGI_{hrc_version}_gxeset_analysis_data_glm.rds")) %>% 
  dplyr::select(-starts_with("pc")) %>% 
  dplyr::filter(!is.na(procmeatqcm)) %>% 
  dplyr::filter(! study_gxe %in% exclude, 
                is.na(adenoma))

table(input_data$study_gxe, input_data$procmeatqc2)
table(input_data$study_gxe, input_data$adenoma, useNA = 'ifany')
table(input_data$study_gxe, input_data$crc, useNA = 'ifany')
table(input_data$study_gxe, input_data$outc, useNA = 'ifany')

# mecc <- filter(input_data, grepl("MECC", study_gxe))
# table(mecc$procmeatqc2, mecc$procmeatqcm)

pca <- fread("/media/work/gwis_test/PCA/20210222/figi_gxe_pca_update.eigenvec") %>% 
  dplyr::rename(vcfid = IID) %>% 
  dplyr::select(-`#FID`) %>% 
  dplyr::rename_with(tolower)

# creating redmeatqcm_v2 for convenience for workflow creation
gxe <- inner_join(input_data, pca, 'vcfid') %>% 
  mutate(procmeatqcm_v2 = as.numeric(procmeatqcm))


wrap(d = gxe, exposure = 'procmeatqcm_v2', hrc_version = 'v3.0', is_e_categorical = F, path = '/media/work/gwis_test/procmeatqcm_v2/data/', studies_to_exclude = "")
check <- readRDS("/media/work/gwis_test/procmeatqcm_v2/data/FIGI_v3.0_gxeset_procmeatqcm_v2_basic_covars_glm.rds")
table(check$study_gxe, check$outcome, useNA = 'ifany')
model <- glm(outcome ~ procmeatqcm_v2 + age_ref_imp + sex + study_gxe, data = check, family = 'binomial')
summary(model)
exp(coef(model))

check <- readRDS("/media/work/gwis_test/procmeatqcm_v2/data/FIGI_v3.0_gxeset_procmeatqcm_v2_basic_covars_gxescan.rds")






