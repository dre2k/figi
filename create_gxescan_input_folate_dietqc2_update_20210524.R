#=============================================================================#
# Create folate_dietqc2 gxescan inputs
# (just to recreate results markdown - can't find the glm input file)
#
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
vcfid_folate_diet <- readRDS("/media/work/gwis_test/folate_dietqc2/data/FIGI_v3.0_gxeset_folate_dietqc2_basic_covars_gxescan.rds") %>% 
  pull(vcfid) # original N = 51649

study_folate_diet <- readRDS(glue("/media/work/gwis/data/FIGI_EpiData/FIGI_v3.0_gxeset_analysis_data_glm.rds")) %>% 
  pull(study) %>% unique(.)

input_data <- readRDS(glue("/media/work/gwis/data/FIGI_EpiData/FIGI_{hrc_version}_gxeset_analysis_data_glm.rds")) %>% 
  filter(study %in% study_folate_diet, 
         !is.na(folate_dietqc2))

exposure = 'folate_dietqc2'

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



wrap(d = input_data, exposure = exposure, hrc_version = hrc_version, is_e_categorical = F, path = glue('/media/work/gwis_test/{exposure}/data/'), studies_to_exclude = "")

# 51649
check <- readRDS(glue("/media/work/gwis_test/{exposure}/data/FIGI_{hrc_version}_gxeset_{exposure}_basic_covars_glm.rds"))
check_gxe <- readRDS(glue("/media/work/gwis_test/{exposure}/data/FIGI_{hrc_version}_gxeset_{exposure}_basic_covars_gxescan.rds"))





