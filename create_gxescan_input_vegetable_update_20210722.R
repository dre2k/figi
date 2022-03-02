#=============================================================================#
# Create fruit/veg/fiber gxescan inputs
# remove adenomas (and matched controls)
#
# use updated PCs
#=============================================================================#
library(tidyverse)
library(data.table)
library(glue)
library(qs)
rm(list = ls())

# if you're going to exclude adenomas, you should use v3.0? seems like some cases were defined in future releases
hrc_version = 'v3.0'

# this is the HRC v3.0  GxE set (before subsetting for presence of exposure info. N = 97391)
pca <- fread("/media/work/gwis_test/PCA/20210222/figi_gxe_pca_update.eigenvec") %>% 
  rename(vcfid = `#FID`) %>% 
  dplyr::select(-IID)


# vcfid_fruit <- readRDS("/media/work/gwis_test/fruitqc2/data/FIGI_v2.3_gxeset_fruitqc2_basic_covars_glm.rds") %>% 
#   pull(vcfid) # original N = 76014
# input_data <- readRDS(glue("/media/work/gwis/data/FIGI_EpiData/FIGI_v2.3_gxeset_analysis_data_glm.rds")) %>% 
#   filter(vcfid %in% vcfid_fruit) %>% 
#   rename_with(tolower, .cols = starts_with("PC")) 


# since samples change between v3.0 and 2.3, i should just subset studies and go from there (to make sure i'm not losing individuals...)
vcfid_veg <- readRDS("/media/work/gwis_test/vegetableqc2/data/FIGI_v2.3_gxeset_vegetableqc2_basic_covars_glm.rds") %>% 
  pull(vcfid) # original N = 76157
study_veg <- readRDS(glue("/media/work/gwis/data/FIGI_EpiData/FIGI_v2.3_gxeset_analysis_data_glm.rds")) %>% 
  pull(study) %>% unique(.)

input_data <- readRDS(glue("/media/work/gwis/data/FIGI_EpiData/FIGI_{hrc_version}_gxeset_analysis_data_glm.rds")) %>% 
  filter(study %in% study_veg, 
         !is.na(vegetable5qcm)) %>% 
  dplyr::select(!starts_with("PC")) %>% 
  inner_join(pca, 'vcfid') %>% 
  rename_with(tolower, .cols = starts_with("PC"))

table(input_data$adenoma, useNA = 'ifany')
table(input_data$crc, input_data$adenoma, useNA = 'ifany')
table(input_data$outcome, input_data$adenoma, useNA = 'ifany')
table(input_data$study, input_data$adenoma, useNA = 'ifany')

# check on NA in both CRC and adenoma
# These are controls - these are blank when using 'crc' variable because if you use crc variable for analysis, these studies would be limited to only controls
# for meat analysis - i removed these studies entirely since the cases were comprised only of adenomas
# clean analysis would exclude. but I'll include because 'outc' is what we've using throughout 
xx <- filter(input_data, is.na(crc) & is.na(adenoma))
table(xx$study_gxe, xx$outcome)

# reach, sms, and hawaii will have to be removed in analysis (controls only after removing adenomas)


#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#

# this should remove adenomas + matched controls (assuming they're matched)
# exclude SMS, REACH, and Hawaii because once AD are removed, only controls left when analyzing using 'outc' variable
# also need to exclude HPFS_3_AD, NHS_3_AD, NHS_5_AD, PLCO_4_AD (sparse)
# note that some additional participants added when using v3.0 compared to 2.3 (~ 600 individuals)
out <- dplyr::filter(input_data, is.na(adenoma))
table(out$study, out$adenoma, useNA = 'ifany')
table(out$study, out$outc, useNA = 'ifany')
table(out$study_gxe, out$outc, useNA = 'ifany')


out <- dplyr::filter(input_data, is.na(adenoma),
                     !study_gxe %in% c("SMS_AD", "REACH_AD", "HawaiiCCS_AD", "HPFS_3_AD", "NHS_3_AD", "NHS_5_AD", "PLCO_4_AD"))
table(out$study_gxe, out$outc, useNA = 'ifany')

out <- out %>% 
  filter(!grepl("DACHS", study_gxe))
table(out$study, out$outc, useNA = 'ifany')




# which exposure to use? 
# either fruitqc2 (original quartile) or fruit5qcm (median 5servings/day, sex study specific)
# let's try fruit5qcm for this first pass 

exposure = 'vegetable5qcm'

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
      dplyr::select(all_of(vars_to_keep)) %>%
      filter(complete.cases(.))}
  else {
    drops <- data.frame(table(tmp$outcome, tmp$study_gxe)) %>%
      filter(Freq <= min_cell_size)
    tmp <- filter(tmp, !study_gxe %in% unique(drops$Var2)) %>%
      dplyr::mutate(study_gxe = fct_drop(study_gxe),
                    !! exposure := as.numeric(.data[[exposure]])) %>%
      dplyr::select(all_of(vars_to_keep)) %>%
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
  tmp <- dplyr::select(tmp, -all_of(ref_study), -study_gxe, -exposure, exposure)
  
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
  saveRDS(cov, file = glue(path, "FIGI_{hrc_version}_gxeset_{exposure}_basic_covars_glm.rds"), version = 2)
  
  cov_gxescan <- format_data_gxescan(cov, exposure)
  saveRDS(cov_gxescan, file = glue(path, "FIGI_{hrc_version}_gxeset_{exposure}_basic_covars_gxescan.rds"), version = 2)
  
}



wrap(d = out, exposure = exposure, hrc_version = hrc_version, is_e_categorical = F, path = glue('/media/work/gwis_test/{exposure}/data/'), studies_to_exclude = "")

# 69599
check <- readRDS(glue("/media/work/gwis_test/{exposure}/data/FIGI_{hrc_version}_gxeset_{exposure}_basic_covars_glm.rds"))
check_gxe <- readRDS(glue("/media/work/gwis_test/{exposure}/data/FIGI_{hrc_version}_gxeset_{exposure}_basic_covars_gxescan.rds"))

table(check$study_gxe, check$outcome, useNA = 'ifany')
model <- glm(outcome ~ vegetable5qcm + age_ref_imp + energytot_imp + sex + study_gxe, data = check, family = 'binomial')
summary(model)
exp(coef(model))






