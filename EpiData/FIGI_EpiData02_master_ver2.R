#=============================================================================#
# FIGI Analysis 04/24/2020
# (redo more thoroughly)
# (edited to incorporate polyE variables)
#
# create master analytical dataset for GWAS and GWIS
#
# Use principal components calculated 190729
# - excludes 557 TCGA samples (genotype data N/A)
#=============================================================================#
library(tidyverse)
library(data.table)
library(figifs)
library(forcats)
library(table1)
rm(list = ls())
output_dir <- "/media/work/gwis/data/"
# setwd(glue("{output_dir}/FIGI_samplefile_epi-201014/"))

# hrc_version = 'v2.3'
# hrc_version = 'v2.4'
hrc_version = 'v3.0'

pca <- "/home/rak/data/PCA/190729/FIGI_GwasSet_190729.eigenvec"
pc <- fread(pca, skip = 1, col.names = c("FID", "IID", paste0(rep("PC", 20), seq(1,20))))

if(hrc_version == "v2.3") {
  load("/media/work/gwis/data/FIGI_EpiData/FIGI_Genotype_Epi_191205.RData")
} else if(hrc_version == "v2.4") {
  load("/media/work/gwis/data/FIGI_EpiData/FIGI_Genotype_Epi_200309.RData")
} else if(hrc_version == "v3.0") { 
  load("/media/work/gwis/data/FIGI_EpiData/FIGI_Genotype_Epi_201014.RData")
}


#-----------------------------------------------------------------------------#
# GWAS SET
#
# Questions for Yi (important for manuscript)
# - remind me again criteria for determinig gxe == 1. some studies are
#   gxe == 0 but still have some exposure information like aspirin
# - calcium supplemental - which variable to use if they want to analyze this var
#
# Comments
# - Q4 vars like redmeatqc2 are categorical (Q4) but modeled continuously. 
#   Code as factors for plots, numerical for meta-analysis/pooled analysis
#-----------------------------------------------------------------------------#
figi_gwas <- inner_join(figi, pc, by = c("vcfid" = "IID")) %>% 
  dplyr::filter(drop == 0) %>% 
  dplyr::mutate(outcome = factor(outc),
                age_ref_imp = as.numeric(age_ref_imp),
                sex = factor(sex, labels = c("Female", "Male")), 
                famhx1 = factor(famhx1), 
                energytot = as.numeric(energytot), 
                energytot_imp = as.numeric(energytot_imp), 
                educ = factor(educ, labels = c("Less than High School", "High School/GED", "Some College", "College/Graduate School")), 
                studyname = ifelse(studyname == "Colo2&3", "Colo23", studyname), 
                study_gxe = ifelse(study_gxe == "Colo2&3", "Colo23", study_gxe), # removing special character
                # create factors for study name variables
                studyname_fct = factor(studyname), # use in gwas analysis
                study_gxe_fct = factor(study_gxe), # use in gxe  analysis
                # studyname_fct = fct_recode(studyname_fct, Colo23 = "Colo2&3"), 
                # study_gxe_fct = fct_recode(study_gxe_fct, Colo23 = "Colo2&3"), 
                # nsaids
                asp_ref = factor(asp_ref),
                nsaids = factor(nsaids),
                aspirin = factor(aspirin),
                # bmi
                bmi = as.numeric(bmi), 
                bmi5 = as.numeric(bmi5), 
                # height
                heightcm = as.numeric(heightcm),
                height10 = as.numeric(height10),
                # smoking
                smoke = factor(smoke),
                smoke = fct_relevel(smoke, "Never smoker", "Former smoker", "Smoker"), 
                smk_ever = factor(smk_ever), 
                smk_pkyrqc2 = factor(smk_pkyrqc2), 
                smk_aveday = as.numeric(smk_aveday),
                smk_pkyr = as.numeric(smk_pkyr),
                # alcohol
                # modified so moderate drinking is ref in alcoholc_moderate
                alcoholc = factor(alcoholc), # keep original labels
                alcoholc = fct_relevel(alcoholc, "nondrinker", "1-28g/d", ">28g/d"),
                alcoholc_heavy = factor(ifelse(alcoholc == "1-28g/d", NA, alcoholc), labels = c("nondrinker", ">28g/d")),
                alcoholc_heavy = fct_drop(alcoholc_heavy),
                alcoholc_moderate = factor(ifelse(alcoholc == ">28g/d", NA, alcoholc), labels = c("nondrinker", "1-28g/d")),
                alcoholc_moderate = fct_drop(alcoholc_moderate), 
                alcoholc_moderate = fct_relevel(alcoholc_moderate, "1-28g/d", "nondrinker"),
                alcoholc_heavy_vs_moderate = ifelse(alcoholc == "1-28g/d", 0,
                                                    ifelse(alcoholc == ">28g/d", 1, NA)),
                alcoholc_heavy_vs_moderate = factor(alcoholc_heavy_vs_moderate, labels = c("1-28g/d", ">28g/d")),
                # HRT
                hrt_ref_pm = factor(hrt_ref_pm), 
                hrt_ref_pm2 = factor(hrt_ref_pm2), 
                eo_ref_pm = factor(eo_ref_pm),
                ep_ref_pm = factor(ep_ref_pm),
                eo_ref_pm_gxe = ifelse(eo_ref_pm == "Yes" & hrt_ref_pm2 == "Yes", "Yes", 
                                       ifelse(hrt_ref_pm2 == "No", "No", NA)), 
                ep_ref_pm_gxe = ifelse(ep_ref_pm == "Yes" & hrt_ref_pm2 == "Yes", "Yes", 
                                       ifelse(hrt_ref_pm2 == "No", "No", NA)),
                eo_ref_pm_gxe = factor(eo_ref_pm_gxe),
                ep_ref_pm_gxe = factor(ep_ref_pm_gxe),
                # pure_eo_allNo = ifelse(eo_ref_pm == "Yes" & ep_ref_pm == "No", "Yes",
                #                        ifelse(eo_ref_pm == "No" & hrt_ref_pm2 == "No" & ep_ref_pm == "No", "No", NA)), 
                # pure_eo_allNo = factor(pure_eo_allNo), 
                # pure_ep_allNo = ifelse(ep_ref_pm == "Yes" & eo_ref_pm == "No", "Yes",
                #                        ifelse(ep_ref_pm == "No" & hrt_ref_pm2 == "No" & eo_ref_pm == "No", "No", NA)),
                # pure_ep_allNo = factor(pure_ep_allNo), 
                # diabetes
                diab = factor(diab), 
                # calcium
                calcium_totqc2 = factor(calcium_totqc2),
                calcium_dietqc2 = factor(calcium_dietqc2),
                calcium_supp = factor(calcium_supp),
                # folate
                folate_totqc2 = factor(folate_totqc2),
                folate_dietqc2 = factor(folate_dietqc2),
                folate_sup_yn = factor(folate_sup_yn), 
                # fiber
                fiberqc2 = factor(fiberqc2),
                # vegetables
                vegetableqc2 = factor(vegetableqc2),
                # fruits
                fruitqc2 = factor(fruitqc2), 
                # meat
                redmeatqc2 = factor(redmeatqc2),
                procmeatqc2 = factor(procmeatqc2), 
                # helper variables to generate table1 with t-test/chisq p values
                table1_outcome = factor(outcome, exclude = NULL, levels = c("Case", "Control", "P-value")),
                table1_outcome = fct_relevel(table1_outcome, "Control", "Case"), # ok to ignore "other" level, gets excluded with gxe == 1
                table1_asp_ref = factor(asp_ref, exclude = NULL, levels = c("No", "Yes", "P-value")))


# saveRDS(figi_gwas, "~/data/FIGI_EpiData_rdata/FIGI_v2.3_GWAS.rds", version = 2)
saveRDS(figi_gwas, glue("{output_dir}/FIGI_EpiData/FIGI_{hrc_version}_GWAS.rds"), version = 2)


#-----------------------------------------------------------------------------#
# GxE SET for analysis (MODEL FITTING)
#
# Notes:
# - gxe == 1 has binary outc "Case/Control"
# - drop these studies (case-only)
#   c("COLOCARE_1", "GALEON", "MOFFITT", "NFCCR_1", "NGCCS")
# - fyi - glm treats character as factor (ordered alphabetically)
#
#
#
# - need to edit study_info for HRC v3.0 (Yi split some of the studies)
#-----------------------------------------------------------------------------#
# study_info <- readRDS("~/data/Annotations/FIGI_studygxe_info.rds")

if(hrc_version == "v2.3") {
  study_info <- readRDS("~/data/Annotations/FIGI_studygxe_info.rds") %>% 
    mutate(study_gxe = as.character(study_gxe))
} else if(hrc_version == "v2.4") {
  study_info <- readRDS("~/data/Annotations/FIGI_studygxe_info.rds") %>% 
    mutate(study_gxe = as.character(study_gxe), 
           study_gxe=replace(study_gxe, study_gxe=="WHI_1", "WHI_1_Rematch"))
} else if(hrc_version == "v3.0") {
  study_info <- readRDS("~/data/Annotations/FIGI_studygxe_info_v3.0.rds") %>% 
    mutate(study_gxe = as.character(study_gxe))
}












#-----------------------------------------------------------------------------#
# create the polyE variable
# (a single score summarizing exposure)
# components need to be scaled, and recoded such that risk directions are consistent
polyE <- inner_join(figi, pc, by = c("vcfid" = "IID")) %>% 
  filter(gxe == 1, EUR_subset == 1) %>% 
  mutate(outcome = factor(outc),
         age_ref_imp = as.numeric(age_ref_imp),
         sex = factor(sex, labels = c("Female", "Male")), 
         energytot_imp = as.numeric(energytot_imp), 
         bmi5 = as.numeric(bmi5), 
         height10 = as.numeric(height10),
         studyname = ifelse(studyname == "Colo2&3", "Colo23", studyname), 
         study_gxe = ifelse(study_gxe == "Colo2&3", "Colo23", study_gxe), # removing special character
         # create factors for study name variables
         studyname_fct = factor(studyname), # use in gwas analysis
         study_gxe_fct = factor(study_gxe), # use in gxe  analysis
         p_redmeatqc2 = (as.numeric(redmeatqc2) - 1) / 3,
         p_procmeatqc2 = (as.numeric(procmeatqc2) - 1) / 3,
         p_calcium_totqc2 = 1 - (as.numeric(calcium_totqc2) - 1) / 3,
         p_fruitqc2 = 1 - (as.numeric(fruitqc2) - 1) / 3,
         p_vegetableqc2 = 1 - (as.numeric(vegetableqc2) - 1) / 3,
         alcoholc_heavy = factor(ifelse(alcoholc == "1-28g/d", NA, alcoholc), labels = c("nondrinker", ">28g/d")),
         alcoholc_heavy = fct_drop(alcoholc_heavy),
         alcoholc_moderate = factor(ifelse(alcoholc == ">28g/d", NA, alcoholc), labels = c("nondrinker", "1-28g/d")),
         alcoholc_moderate = fct_drop(alcoholc_moderate),
         p_alcoholc_moderate = 1 - (as.numeric(alcoholc_moderate) - 1),
         p_alcoholc_heavy = (as.numeric(alcoholc_heavy) - 1), 
         p_nsaids = 1 - (as.numeric(factor(asp_ref)) - 1),
         p_smk_ever = as.numeric(factor(smk_ever)) - 1) %>% 
  mutate(p_diet = p_redmeatqc2 + p_procmeatqc2 + p_calcium_totqc2 + p_fruitqc2 + p_vegetableqc2,
         p_nas = p_nsaids + p_smk_ever + p_alcoholc_moderate + p_alcoholc_heavy,
         p_all = p_diet + p_nas) %>% 
  mutate(p_diet_std = p_diet / max(p_diet,na.rm = T),
         p_nas_std = p_nas / max(p_nas, na.rm = T),
         p_all_std = p_all / max(p_all, na.rm = T)) %>% 
  dplyr::select(vcfid, p_diet_std, p_nas_std, p_all_std)


figi_gxe <- figi %>% 
  dplyr::filter(gxe == 1, EUR_subset == 1,
                !study_gxe %in% c("ColoCare_1", "GALEON", "MOFFITT", "NFCCR_1", "NGCCS")) %>% 
  dplyr::mutate(studyname = ifelse(studyname == "Colo2&3", "Colo23", studyname), 
                study_gxe = ifelse(study_gxe == "Colo2&3", "Colo23", study_gxe)) %>% 
  inner_join(pc, by = c("vcfid" = "IID")) %>% 
  inner_join(study_info, 'study_gxe') %>% 
  inner_join(polyE, 'vcfid') %>% 
  # A E S T H E T I C S
  dplyr::rename_all(tolower) %>% 
  rename(EUR_subset = eur_subset) %>% 
  # variables
  dplyr::mutate(study_design = fct_drop(study_design), 
                study_design = fct_relevel(study_design, "Cohort"),
                outcome = ifelse(outc == "Control", 0, 1),
                age_ref_imp = as.numeric(age_ref_imp),
                sex = ifelse(sex == "Female", 0, 1),
                famhx1 = factor(famhx1), 
                energytot = as.numeric(energytot), 
                energytot_imp = as.numeric(energytot_imp), 
                educ = factor(educ, labels = c("Less than High School", "High School/GED", "Some College", "College/Graduate School")), 
                # nsaids
                asp_ref = factor(asp_ref),
                nsaids = factor(nsaids),
                aspirin = factor(aspirin),
                # bmi
                bmi = as.numeric(bmi), 
                bmi5 = as.numeric(bmi5),
                # height
                heightcm = as.numeric(heightcm),
                height10 = as.numeric(height10),
                # smoking
                smoke = factor(smoke),
                smoke = fct_relevel(smoke, "Never smoker", "Former smoker", "Smoker"), 
                smk_ever = factor(smk_ever), 
                smk_pkyrqc2 = factor(smk_pkyrqc2), 
                smk_aveday = as.numeric(smk_aveday),
                smk_pkyr = as.numeric(smk_pkyr),
                # alcohol (moderate drinking is ref group)
                alcoholc = factor(alcoholc),
                alcoholc = fct_relevel(alcoholc, "1-28g/d", "nondrinker", ">28g/d"),
                alcoholc_heavy = fct_drop(na_if(alcoholc, "1-28g/d")),
                alcoholc_moderate = fct_drop(na_if(alcoholc, ">28g/d")),
                alcoholc_heavy_vs_moderate = fct_drop(na_if(alcoholc, "nondrinker")),
                # HRT
                hrt_ref_pm = factor(hrt_ref_pm), 
                hrt_ref_pm2 = factor(hrt_ref_pm2), 
                eo_ref_pm = factor(eo_ref_pm),
                ep_ref_pm = factor(ep_ref_pm),
                eo_ref_pm_gxe = factor(ifelse(eo_ref_pm == "Yes" & hrt_ref_pm2 == "Yes", "Yes", 
                                       ifelse(hrt_ref_pm2 == "No", "No", NA))), 
                ep_ref_pm_gxe = factor(ifelse(ep_ref_pm == "Yes" & hrt_ref_pm2 == "Yes", "Yes", 
                                       ifelse(hrt_ref_pm2 == "No", "No", NA))),
                eo_ref_pm_gxe = factor(eo_ref_pm_gxe),
                ep_ref_pm_gxe = factor(ep_ref_pm_gxe),
                # pure_eo_allNo = ifelse(eo_ref_pm == "Yes" & ep_ref_pm == "No", "Yes",
                #                        ifelse(eo_ref_pm == "No" & hrt_ref_pm2 == "No" & ep_ref_pm == "No", "No", NA)), 
                # pure_eo_allNo = factor(pure_eo_allNo), 
                # pure_ep_allNo = ifelse(ep_ref_pm == "Yes" & eo_ref_pm == "No", "Yes",
                #                        ifelse(ep_ref_pm == "No" & hrt_ref_pm2 == "No" & eo_ref_pm == "No", "No", NA)),
                # pure_ep_allNo = factor(pure_ep_allNo), 
                # diabetes
                diab = factor(diab), 
                # calcium
                calcium_totqc2 = as.numeric(calcium_totqc2)-1,
                calcium_dietqc2 = as.numeric(calcium_dietqc2)-1,
                calcium_supp = factor(calcium_supp, labels = c("No","Yes")),
                # folate
                folate_totqc2 = as.numeric(folate_totqc2)-1,
                folate_dietqc2 = as.numeric(folate_dietqc2)-1,
                folate_sup_yn = factor(folate_sup_yn, labels = c("No","Yes")),
                # fiber
                fiberqc2 = as.numeric(fiberqc2)-1,
                # vegetables
                vegetableqc2 = as.numeric(vegetableqc2)-1,
                # fruits
                fruitqc2 = as.numeric(fruitqc2)-1, 
                # meat
                redmeatqc2 = as.numeric(redmeatqc2)-1,
                procmeatqc2 = as.numeric(procmeatqc2)-1)
                
for(x in c("asp_ref", "nsaids", "aspirin", "smk_ever", "alcoholc", "alcoholc_heavy", "alcoholc_moderate", "alcoholc_heavy_vs_moderate", "hrt_ref_pm2", "eo_ref_pm_gxe", "ep_ref_pm_gxe", "diab", "calcium_totqc2", "calcium_dietqc2", "calcium_supp", "folate_totqc2", "folate_dietqc2", "folate_sup_yn", "fiberqc2", "vegetableqc2", "fruitqc2", "redmeatqc2", "procmeatqc2")) {
  print(x)
  print(class(figi_gxe[, x]))
  print(table(figi_gxe[, x]))
}

# saveRDS(figi_gxe, file = paste0("/media/work/gwis/results/input/FIGI_", hrc_version, "_gxeset_analysis_data_glm.rds"), version = 2)
saveRDS(figi_gxe, file = glue("{output_dir}/FIGI_EpiData/FIGI_", hrc_version, "_gxeset_analysis_data_glm.rds"), version = 2)




#-----------------------------------------------------------------------------#
# GxE SET for creating descriptive table (table1 package)
#
# Notes:
# - binary/quartile variables should be factors
# - binary/quartile create helper variable to tabulate vars by exposures
#-----------------------------------------------------------------------------#

binary_to_factor_table1 <- function(x) {
  factor(x, exclude = NA, levels = c("No", "Yes", "P-value"))
}
quartile_to_factor_table1 <- function(x) {
  factor(x, exclude = NA, levels = c("0", "1", "2", "3", "P-value"), labels = c("Q1", "Q2", "Q3", "Q4", "P-value"))
}

binary_vars <- c("asp_ref", "nsaids", "aspirin", "smk_ever", "hrt_ref_pm2", "eo_ref_pm_gxe", "ep_ref_pm_gxe", "diab", "calcium_supp", "folate_sup_yn")
quartile_vars <- c("calcium_totqc2", "calcium_dietqc2", "folate_totqc2", "folate_dietqc2", "fiberqc2", "vegetableqc2", "fruitqc2", "redmeatqc2", "procmeatqc2")

figi_gxe_table1 <- figi_gxe %>% 
  dplyr::mutate(outcome_table1 = factor(outcome, exclude = NA, levels = c("0", "1", "P-value"), labels = c("Controls", "Cases", "P-value")), 
                alcoholc_table1 = factor(alcoholc, exclude = NA, levels = c("1-28g/d", "nondrinker", ">28g/d" , "P-value")),
                alcoholc_moderate_table1 = factor(alcoholc_moderate, exclude = NA, levels = c("1-28g/d", "nondrinker", "P-value")),
                alcoholc_heavy_vs_moderate_table1 = factor(alcoholc_heavy_vs_moderate, exclude = NA, levels = c("1-28g/d", ">28g/d", "P-value"))) %>% 
  mutate_at(binary_vars, .funs = list(table1 = ~binary_to_factor_table1(.))) %>% 
  mutate_at(quartile_vars, .funs = list(table1 = ~quartile_to_factor_table1(.)))


label(figi_gxe_table1$age_ref_imp) <- "Age (mean imputed)"
label(figi_gxe_table1$sex) <- "Sex"
label(figi_gxe_table1$outcome) <- "Outcome"
label(figi_gxe_table1$educ) <- "Education (highest completed)"
label(figi_gxe_table1$famhx1) <- "Family History of CRC"

label(figi_gxe_table1$alcoholc) <- "Alcohol intake - 3 categories"
label(figi_gxe_table1$alcoholc_moderate) <- "Moderate alcohol intake"
label(figi_gxe_table1$alcoholc_heavy) <- "Heavy alcohol intake"
label(figi_gxe_table1$alcoholc_heavy_vs_moderate) <- "Heavy vs. moderate alcohol intake"

label(figi_gxe_table1$bmi) <- "BMI"
label(figi_gxe_table1$bmi5) <- "BMI/5, set to NA if BMI<18.5"

label(figi_gxe_table1$calcium_totqc2) <- "Total calcium intake (diet + supp)"
label(figi_gxe_table1$calcium_dietqc2) <- "Dietary calcium intake"

label(figi_gxe_table1$fiberqc2) <- "Total fiber intake"
label(figi_gxe_table1$fruitqc2) <- "Total fruit intake"
label(figi_gxe_table1$vegetableqc2) <- "Total vegetable intake"

label(figi_gxe_table1$folate_totqc2) <- "Total folate intake"
label(figi_gxe_table1$folate_dietqc2) <- "Dietary folate intake"

label(figi_gxe_table1$heightcm) <- "Height"
units(figi_gxe_table1$heightcm) <- "cm"

label(figi_gxe_table1$hrt_ref_pm2) <- "Any Post-Menopausal HRT Use"
label(figi_gxe_table1$eo_ref_pm_gxe) <- "Post-menopausal HRT (estrogen only) use"
label(figi_gxe_table1$ep_ref_pm_gxe) <- "Post-menopausal HRT (estrogen+progesterone) use"
# label(figi_gxe_table1$pure_eo_allNo) <- "Post-menopausal HRT (estrogen only) use - pure def"
# label(figi_gxe_table1$pure_ep_allNo) <- "Post-menopausal HRT (estrogen+progesterone) use - pure def"


label(figi_gxe_table1$asp_ref) <- "Reg aspirin or NSAID use"
label(figi_gxe_table1$aspirin) <- "Reg aspirin use"
label(figi_gxe_table1$nsaids) <- "Ref NSAID use"

label(figi_gxe_table1$procmeatqc2) <- "Total processed meat intake"
label(figi_gxe_table1$redmeatqc2) <- "Total red meat intake"

label(figi_gxe_table1$smk_ever) <- "Smoking, never/ever"
label(figi_gxe_table1$smk_aveday) <- "Smoking, avg cig/day"
label(figi_gxe_table1$smk_pkyr) <- "Smoking, avg pks/yr"

label(figi_gxe_table1$diab) <- "T2D (ever diagnosed)"

# saveRDS(figi_gxe_table1, file = paste0("/media/work/gwis/results/input/FIGI_", hrc_version, "_gxeset_analysis_data_table1.rds"), version = 2)
saveRDS(figi_gxe_table1, file = glue("{output_dir}/FIGI_EpiData/FIGI_", hrc_version, "_gxeset_analysis_data_table1.rds"), version = 2)


for(x in paste0(binary_vars, "_table1")) {
  print(x)
  print(class(figi_gxe_table1[, x]))
  print(table(figi_gxe_table1[, x]))
}



#-----------------------------------------------------------------------------#
# GxE SET for creating ggplot (summary stats of exposures) 
#
# Notes:
# - setting NA as first factor level (for binary/quartile variables)
#-----------------------------------------------------------------------------#

binary_factor_with_na <- function(x) {
  # binary variables are already factors, good to go
  fct_relevel(fct_explicit_na(x, na_level = "NA"), "NA", "Yes", "No")
}

binary_factor_with_na_counts_table <- function(x) {
  fct_relevel(fct_explicit_na(x, na_level = "NA"), "No", "Yes", "NA")
}

quartile_factor_with_na <- function(x) {
  # first, set as factor with labels Q1-Q4
  tmp <- factor(x, levels = c("0", "1", "2", "3"), labels = c("Q1", "Q2", "Q3", "Q4"))
  fct_relevel(fct_explicit_na(tmp, na_level = "NA"), "NA", "Q1", "Q2", "Q3", "Q4")
}

quartile_factor_with_na_counts_table <- function(x) {
  tmp <- factor(x, levels = c("0", "1", "2", "3"), labels = c("Q1", "Q2", "Q3", "Q4"))
  fct_relevel(fct_explicit_na(tmp, na_level = "NA"), "Q1", "Q2", "Q3", "Q4", "NA")
}

binary_vars <- c("asp_ref", "nsaids", "aspirin", "smk_ever", "hrt_ref_pm2", "eo_ref_pm_gxe", "ep_ref_pm_gxe", "diab", "calcium_supp", "folate_sup_yn")
quartile_vars <- c("calcium_totqc2", "calcium_dietqc2", "folate_totqc2", "folate_dietqc2", "fiberqc2", "vegetableqc2", "fruitqc2", "redmeatqc2", "procmeatqc2")

# data for ggplot (bar and box plots). Setting NA as first factor level
# not excluding subsets of individuals since I want to show missing values by studies
# data also includes variables for counts table
figi_gxe_ggplot <- figi_gxe %>% 
  mutate(alcoholc_ggplot = fct_relevel(fct_explicit_na(alcoholc, na_level = "NA"), "NA", ">28g/d", "nondrinker", "1-28g/d"),
         alcoholc_moderate_ggplot = fct_relevel(fct_explicit_na(alcoholc_moderate, na_level = "NA"), "NA", "nondrinker", "1-28g/d"),
         alcoholc_heavy_vs_moderate_ggplot = fct_relevel(fct_explicit_na(alcoholc_heavy_vs_moderate, na_level = "NA"), "NA", ">28g/d", "1-28g/d"),
         alcoholc_count = fct_relevel(fct_explicit_na(alcoholc, na_level = "NA"), "1-28g/d", "nondrinker", ">28g/d", "NA"),
         alcoholc_moderate_count = fct_relevel(fct_explicit_na(alcoholc_moderate, na_level = "NA"), "1-28g/d", "nondrinker", "NA"),
         alcoholc_heavy_vs_moderate_count = fct_relevel(fct_explicit_na(alcoholc_heavy_vs_moderate, na_level = "NA"), "1-28g/d", ">28g/d", "NA")) %>% 
  mutate_at(binary_vars, .funs = list(ggplot = ~binary_factor_with_na(.))) %>% 
  mutate_at(binary_vars, .funs = list(count = ~binary_factor_with_na_counts_table(.))) %>% 
  mutate_at(quartile_vars, .funs = list(ggplot = ~quartile_factor_with_na(.))) %>% 
  mutate_at(quartile_vars, .funs = list(count = ~quartile_factor_with_na_counts_table(.)))

# for(x in paste0(binary_vars, "_ggplot")) {
#   print(x)
#   print(class(figi_gxe_ggplot[, x]))
#   print(table(figi_gxe_ggplot[, x]))
# }
# 
# for(x in paste0(quartile_vars, "_ggplot")) {
#   print(x)
#   print(class(figi_gxe_ggplot[, x]))
#   print(table(figi_gxe_ggplot[, x]))
# }

# saveRDS(figi_gxe_ggplot, file = paste0("/media/work/gwis/results/input/FIGI_", hrc_version, "_gxeset_analysis_data_ggplot.rds"), version = 2)
saveRDS(figi_gxe_ggplot, file = glue("{output_dir}/FIGI_EpiData/FIGI_", hrc_version, "_gxeset_analysis_data_ggplot.rds"), version = 2)




#-----------------------------------------------------------------------------#
# csv file for jim
#-----------------------------------------------------------------------------#
# rm(list = ls())
# figi_gwas <- readRDS("~/data/FIGI_EpiData_rdata/FIGI_HRC_v2.3_GWAS.rds")
# 
# write.csv(figi_gwas, file = "/media/work/tmp/FIGI_v2.3_gwas_set_N138014_20191211.csv", quote = T, row.names = F)



