#=============================================================================#
# Create redmeat/procmeat gxescan inputs
#=============================================================================#
# ---- explore inclusion/exclusion ---- 
library(tidyverse)
library(data.table)
rm(list = ls())


hrc_version = 'v3.0'
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


# ---- folate_diet400qcm ----
hrc_version = 'v3.0'
input_data <- readRDS(glue("/media/work/gwis/data/FIGI_EpiData/FIGI_{hrc_version}_gxeset_analysis_data_glm.rds"))

# make sure sample sizes are the same?
sum(!is.na(input_data$folate_dietqc2))
sum(!is.na(input_data$folate_diet400qcm))


# remove old PC values (... should I?)
exclude = c("UKB_1", "SMS_AD")

input_data <- readRDS(glue("/media/work/gwis/data/FIGI_EpiData/FIGI_{hrc_version}_gxeset_analysis_data_glm.rds")) %>% 
  dplyr::select(-starts_with("pc")) %>% 
  dplyr::filter(!is.na(folate_diet400qcm)) %>% 
  dplyr::filter(!study_gxe %in% exclude)

x <- data.frame(table(input_data$study_gxe, input_data$folate_diet400qcm, input_data$sex)) %>% 
  dplyr::filter(Freq != 0) %>% 
  arrange(Var1, Var3)

pca <- fread("/media/work/gwis_test/PCA/20210222/figi_gxe_pca_update.eigenvec") %>% 
  dplyr::rename(vcfid = IID) %>% 
  dplyr::select(-`#FID`) %>% 
  dplyr::rename_with(tolower)

# creating redmeatqcm_v2 for convenience for workflow creation
gxe <- inner_join(input_data, pca, 'vcfid') %>% 
  mutate(folate_diet400qcm = as.numeric(folate_diet400qcm))


wrap(d = gxe, exposure = 'folate_diet400qcm', hrc_version = 'v3.0', is_e_categorical = F, path = '/media/work/gwis_test/folate_diet400qcm/data/', studies_to_exclude = "")
check <- readRDS("/media/work/gwis_test/folate_diet400qcm/data/FIGI_v3.0_gxeset_folate_diet400qcm_basic_covars_glm.rds")
table(check$study_gxe, check$outcome, useNA = 'ifany')
model <- glm(outcome ~ folate_diet400qcm + age_ref_imp + energytot_imp + sex + study_gxe, data = check, family = 'binomial')
summary(model)
exp(coef(model))

check <- readRDS("/media/work/gwis_test/folate_diet400qcm/data/FIGI_v3.0_gxeset_folate_diet400qcm_basic_covars_gxescan.rds")









#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#


# ---- folate_sup ----
# for supplemental folate, i will include UKB in analysis
hrc_version = 'v3.0'
input_data <- readRDS(glue("/media/work/gwis/data/FIGI_EpiData/FIGI_{hrc_version}_gxeset_analysis_data_glm.rds"))

# make sure sample sizes are the same?
sum(!is.na(input_data$folate_sup))


# remove old PC values (... should I?)
exclude = c("SMS_AD")
input_data <- readRDS(glue("/media/work/gwis/data/FIGI_EpiData/FIGI_{hrc_version}_gxeset_analysis_data_glm.rds")) %>% 
  dplyr::select(-starts_with("pc")) %>% 
  dplyr::filter(!is.na(folate_sup)) %>% 
  dplyr::filter(!study_gxe %in% exclude)

# x <- data.frame(table(input_data$study_gxe, input_data$folate_sup, input_data$sex)) %>% 
#   dplyr::filter(Freq != 0) %>% 
#   arrange(Var1, Var3)

pca <- fread("/media/work/gwis_test/PCA/20210222/figi_gxe_pca_update.eigenvec") %>% 
  dplyr::rename(vcfid = IID) %>% 
  dplyr::select(-`#FID`) %>% 
  dplyr::rename_with(tolower)

# creating redmeatqcm_v2 for convenience for workflow creation
gxe <- inner_join(input_data, pca, 'vcfid') %>% 
  mutate(folate_sup = as.numeric(folate_sup) / 400)


wrap(d = gxe, exposure = 'folate_sup', hrc_version = 'v3.0', is_e_categorical = F, path = '/media/work/gwis_test/folate_sup/data/', studies_to_exclude = "")
check <- readRDS("/media/work/gwis_test/folate_sup/data/FIGI_v3.0_gxeset_folate_sup_basic_covars_glm.rds")
table(check$study_gxe, check$outcome, useNA = 'ifany')
model <- glm(outcome ~ folate_sup + age_ref_imp + energytot_imp + sex + study_gxe, data = check, family = 'binomial')
summary(model)
exp(coef(model))

check <- readRDS("/media/work/gwis_test/folate_sup/data/FIGI_v3.0_gxeset_folate_sup_basic_covars_gxescan.rds")




#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#

# ---- folate_sup_yn ----
# for supplemental folate, i will include UKB in analysis
hrc_version = 'v3.0'
input_data <- readRDS(glue("/media/work/gwis/data/FIGI_EpiData/FIGI_{hrc_version}_gxeset_analysis_data_glm.rds"))

# make sure sample sizes are the same?
sum(!is.na(input_data$folate_sup_yn))


# remove old PC values (... should I?)
exclude = c("SMS_AD")
input_data <- readRDS(glue("/media/work/gwis/data/FIGI_EpiData/FIGI_{hrc_version}_gxeset_analysis_data_glm.rds")) %>% 
  dplyr::select(-starts_with("pc")) %>% 
  dplyr::filter(!is.na(folate_sup_yn)) %>% 
  dplyr::filter(!study_gxe %in% exclude)

# x <- data.frame(table(input_data$study_gxe, input_data$folate_sup, input_data$sex)) %>% 
#   dplyr::filter(Freq != 0) %>% 
#   arrange(Var1, Var3)

pca <- fread("/media/work/gwis_test/PCA/20210222/figi_gxe_pca_update.eigenvec") %>% 
  dplyr::rename(vcfid = IID) %>% 
  dplyr::select(-`#FID`) %>% 
  dplyr::rename_with(tolower)

# creating redmeatqcm_v2 for convenience for workflow creation
gxe <- inner_join(input_data, pca, 'vcfid') %>% 
  mutate(folate_sup_yn = ifelse(folate_sup_yn == "No", 0, 1))


wrap(d = gxe, exposure = 'folate_sup_yn', hrc_version = 'v3.0', is_e_categorical = F, path = '/media/work/gwis_test/folate_sup_yn/data/', studies_to_exclude = "")
check <- readRDS("/media/work/gwis_test/folate_sup_yn/data/FIGI_v3.0_gxeset_folate_sup_yn_basic_covars_glm.rds")
table(check$study_gxe, check$outcome, useNA = 'ifany')
model <- glm(outcome ~ folate_sup_yn + age_ref_imp + energytot_imp + sex + study_gxe, data = check, family = 'binomial')
summary(model)
exp(coef(model))

check <- readRDS("/media/work/gwis_test/folate_sup_yn/data/FIGI_v3.0_gxeset_folate_sup_yn_basic_covars_gxescan.rds")




#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#

# ---- folate_tot ----
# for supplemental folate, i will include UKB in analysis
hrc_version = 'v3.0'
input_data <- readRDS(glue("/media/work/gwis/data/FIGI_EpiData/FIGI_{hrc_version}_gxeset_analysis_data_glm.rds"))

# make sure sample sizes are the same?
sum(!is.na(input_data$folate_tot))


# remove old PC values (... should I?)
exclude = c("UKB_1", "SMS_AD")
input_data <- readRDS(glue("/media/work/gwis/data/FIGI_EpiData/FIGI_{hrc_version}_gxeset_analysis_data_glm.rds")) %>% 
  dplyr::select(-starts_with("pc")) %>% 
  dplyr::filter(!is.na(folate_sup_yn)) %>% 
  dplyr::filter(!study_gxe %in% exclude)

# x <- data.frame(table(input_data$study_gxe, input_data$folate_sup, input_data$sex)) %>% 
#   dplyr::filter(Freq != 0) %>% 
#   arrange(Var1, Var3)

pca <- fread("/media/work/gwis_test/PCA/20210222/figi_gxe_pca_update.eigenvec") %>% 
  dplyr::rename(vcfid = IID) %>% 
  dplyr::select(-`#FID`) %>% 
  dplyr::rename_with(tolower)

gxe <- inner_join(input_data, pca, 'vcfid') %>% 
  mutate(folate_tot = as.numeric(folate_tot) / 400)

wrap(d = gxe, exposure = 'folate_tot', hrc_version = 'v3.0', is_e_categorical = F, path = '/media/work/gwis_test/folate_tot/data/', studies_to_exclude = "")
check <- readRDS("/media/work/gwis_test/folate_tot/data/FIGI_v3.0_gxeset_folate_tot_basic_covars_glm.rds")
# table(check$study_gxe, check$outcome, useNA = 'ifany')
model <- glm(outcome ~ folate_tot + age_ref_imp + energytot_imp + sex + study_gxe, data = check, family = 'binomial')
summary(model)
exp(coef(model))

check <- readRDS("/media/work/gwis_test/folate_tot/data/FIGI_v3.0_gxeset_folate_tot_basic_covars_gxescan.rds")
