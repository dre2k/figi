#=============================================================================#
# FIGI Poly E Analysis
# using hrc v2.4
#
# 1) submit G x PolyE scan
# 2) clustering?
#=============================================================================#
library(tidyverse)
library(data.table)
library(figifs)
library(forcats)
rm(list = ls())
pca <- "/home/rak/data/PCA/190729/FIGI_GwasSet_190729.eigenvec"
pc <- fread(pca, skip = 1, col.names = c("FID", "IID", paste0(rep("PC", 20), seq(1,20))))
load("~/data/FIGI_EpiData_rdata/FIGI_Genotype_Epi_200309.RData")

# we should limit to GxE set
# create 3 sets of poly E variable
# 1) smoke/drink
# 2) diet
# 3) all the above

gxe <- inner_join(figi, pc, by = c("vcfid" = "IID")) %>% 
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
         p_all_std = p_all / max(p_all, na.rm = T))

table(gxe$alcoholc, gxe$p_alcoholc)
table(gxe$asp_ref, gxe$p_nsaids)
table(gxe$alcoholc_moderate, gxe$p_alcoholc_moderate)
table(gxe$alcoholc_heavy, gxe$p_alcoholc_heavy)
table(gxe$smk_ever, gxe$p_smk_ever)

hist(gxe$p_diet)
hist(gxe$p_diet_std)
hist(gxe$p_nas_std)
hist(gxe$p_all_std)


# format data for gxescan
gxe_glm <- format_data_glm(gxe, exposure = 'p_diet_std', is_e_categorical = F, vars_to_exclude = '', vars_to_include = c('bmi5', 'height10'), eur_only = T)
saveRDS(gxe_glm , file = "~/data/results/input/FIGI_v2.4_gxeset_p_diet_std_basic_covars_bmi5_height10_glm.rds")
gxe_gxescan <- format_data_gxescan(gxe_glm, exposure = 'p_diet_std')
saveRDS(gxe_gxescan , file = "~/data/results/input/FIGI_v2.4_gxeset_p_diet_std_basic_covars_bmi5_height10_gxescan.rds")





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

# # test glm
# test_char <- glm(outcome ~ asp_ref + age_ref_imp, data = figi, family = 'binomial')
# summary(test)
# 
# test_factor <- glm(outcome ~ asp_ref + age_ref_imp, data = figi_gwas, family = 'binomial')
# summary(test)
# 
# # test table1
# table1(~ asp_ref + age_ref_imp | outcome, data=figi_gwas)

saveRDS(figi_gwas, "~/data/FIGI_EpiData_rdata/FIGI_HRC_v2.4_GWAS.rds", version = 2)



#-----------------------------------------------------------------------------#
# csv file for jim
#-----------------------------------------------------------------------------#
# rm(list = ls())
# figi_gwas <- readRDS("~/data/FIGI_EpiData_rdata/FIGI_HRC_v2.4_GWAS.rds")
# 
# write.csv(figi_gwas, file = "/media/work/tmp/FIGI_v2.4_gwas_set_N138014_20191211.csv", quote = T, row.names = F)



#-----------------------------------------------------------------------------#
# data set for E main effects analysis (Rmarkdown results report)
#
#
# EDIT - change the names to distinguish from v2.3
# (you'll need to run this again when you generate rmarkdown results report for NSAIDs)
#-----------------------------------------------------------------------------#
figi_gwas <- readRDS("~/data/FIGI_EpiData_rdata/FIGI_HRC_v2.4_GWAS.rds")

# exposure_subset <- readRDS(paste0("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_", params$exposure, "_basic_covars_glm.rds"))$vcfid

# data for meta-analysis and pooled analysis
# categorical variables (e.g. Q4) are modeled CONTINUOUSLY
factor_to_numeric <- function(x) {
  as.numeric(x) - 1
}

gxe <- figi_gwas %>%
  # dplyr::filter(vcfid %in% exposure_subset) %>% # subset this specifically for each E. Use vcfids included in gxe scans
  dplyr::mutate(outcome = ifelse(outcome == "Control", 0, 1),
                sex = ifelse(sex == "Female", 0, 1)) %>% 
  dplyr::mutate_at(c("calcium_totqc2", "calcium_dietqc2","folate_totqc2", "folate_dietqc2", "fiberqc2", "vegetableqc2", "fruitqc2", "redmeatqc2", "procmeatqc2"), factor_to_numeric)

saveRDS(gxe, "~/data/FIGI_EpiData_rdata/FIGI_results_exposure_main_effects.rds")


# data for descriptive statistics (tables and plots)
# if something looks weird, likely there were errors in this step
# (data used with table1 package)
gxe_table1 <- figi_gwas %>%
  # dplyr::filter(vcfid %in% exposure_subset) %>%
  dplyr::filter(outcome != "Other") %>% 
  dplyr::mutate(outcome = droplevels.factor(outcome), 
                table1_asp_ref = factor(asp_ref, exclude = NA, levels = c("No", "Yes", "P-value")),
                asp_ref = fct_relevel(asp_ref, "No", "Yes"),
                table1_aspirin = factor(aspirin, exclude = NA, levels = c("No", "Yes", "P-value")),
                aspirin = fct_relevel(aspirin, "No", "Yes"),
                table1_nsaids = factor(nsaids, exclude = NA, levels = c("No", "Yes", "P-value")),
                nsaids = fct_relevel(nsaids, "No", "Yes"),
                table1_folate_totqc2 = factor(folate_totqc2, exclude = NA, levels = c("1", "2", "3", "4", "P-value")),
                folate_totqc2 = fct_relevel(folate_totqc2, "1", "2", "3", "4"),
                table1_folate_dietqc2 = factor(folate_dietqc2, exclude = NA, levels = c("1", "2", "3", "4", "P-value")),
                folate_dietqc2 = fct_relevel(folate_dietqc2, "1", "2", "3", "4"),
                table1_calcium_totqc2 = factor(calcium_totqc2, exclude = NA, levels = c("1", "2", "3", "4", "P-value")),
                calcium_totqc2 = fct_relevel(calcium_totqc2, "1", "2", "3", "4"),
                table1_calcium_dietqc2 = factor(calcium_dietqc2, exclude = NA, levels = c("1", "2", "3", "4", "P-value")),
                calcium_dietqc2 = fct_relevel(calcium_dietqc2, "1", "2", "3", "4"),
                table1_alcoholc_heavy = factor(alcoholc_heavy, exclude = NA, levels = c("nondrinker", ">28g/d", "P-value")),
                alcoholc_heavy = fct_relevel(alcoholc_heavy, "nondrinker", ">28g/d"),
                table1_alcoholc_moderate = factor(alcoholc_moderate, exclude = NA, levels = c("nondrinker", "1-28g/d", "P-value")),
                alcoholc_moderate = fct_relevel(alcoholc_moderate, "nondrinker", "1-28g/d"),
                table1_alcoholc_heavy_vs_moderate = factor(alcoholc_heavy_vs_moderate, exclude = NA, levels = c("1-28g/d", ">28g/d", "P-value")),
                alcoholc_heavy_vs_moderate = fct_relevel(alcoholc_heavy_vs_moderate, "1-28g/d", ">28g/d"),
                table1_hrt_ref_pm = factor(hrt_ref_pm, exclude = NA, levels = c("No", "Yes", "P-value")),
                hrt_ref_pm = fct_relevel(hrt_ref_pm, "No", "Yes"),
                table1_hrt_ref_pm2 = factor(hrt_ref_pm2, exclude = NA, levels = c("No", "Yes", "P-value")),
                hrt_ref_pm2 = fct_relevel(hrt_ref_pm2, "No", "Yes"),
                table1_eo_ref_pm_gxe = factor(eo_ref_pm_gxe, exclude = NA, levels = c("No", "Yes", "P-value")),
                eo_ref_pm_gxe = fct_relevel(eo_ref_pm_gxe, "No", "Yes"),
                table1_ep_ref_pm_gxe = factor(ep_ref_pm_gxe, exclude = NA, levels = c("No", "Yes", "P-value")),
                ep_ref_pm_gxe = fct_relevel(ep_ref_pm_gxe, "No", "Yes"),
                table1_diab = factor(diab, exclude = NA, levels = c("No", "Yes", "P-value")),
                diab = fct_relevel(diab, "No", "Yes"),
                table1_smk_ever = factor(smk_ever, exclude = NA, levels = c("No", "Yes", "P-value")),
                smk_ever = fct_relevel(smk_ever, "No", "Yes"),
                table1_fiberqc2 = factor(fiberqc2, exclude = NA, levels = c("1", "2", "3", "4", "P-value")),
                fiberqc2 = fct_relevel(fiberqc2, "1", "2", "3", "4"),
                table1_fruitqc2 = factor(fruitqc2, exclude = NA, levels = c("1", "2", "3", "4", "P-value")),
                fruitqc2 = fct_relevel(fruitqc2, "1", "2", "3", "4"), 
                table1_vegetableqc2 = factor(vegetableqc2, exclude = NA, levels = c("1", "2", "3", "4", "P-value")),
                vegetableqc2 = fct_relevel(vegetableqc2, "1", "2", "3", "4"),
                table1_redmeatqc2 = factor(redmeatqc2, exclude = NA, levels = c("1", "2", "3", "4", "P-value")),
                redmeatqc2 = fct_relevel(redmeatqc2, "1", "2", "3", "4"),
                table1_procmeatqc2 = factor(procmeatqc2, exclude = NA, levels = c("1", "2", "3", "4", "P-value")),
                procmeatqc2 = fct_relevel(procmeatqc2, "1", "2", "3", "4"))

label(gxe_table1$age_ref_imp) <- "Age (mean imputed)"
label(gxe_table1$sex) <- "Sex"
label(gxe_table1$outcome) <- "Outcome"
label(gxe_table1$educ) <- "Education (highest completed)"
label(gxe_table1$famhx1) <- "Family History of CRC"

label(gxe_table1$alcoholc) <- "Alcohol intake - 3 categories"
label(gxe_table1$alcoholc_moderate) <- "Noderate alcohol intake"
label(gxe_table1$alcoholc_heavy) <- "Heavy alcohol intake"

label(gxe_table1$bmi) <- "BMI"
label(gxe_table1$bmi5) <- "BMI/5, set to NA if BMI<18.5"

label(gxe_table1$calcium_totqc2) <- "Total calcium intake (diet + supp)"
label(gxe_table1$calcium_dietqc2) <- "Dietary calcium intake"

label(gxe_table1$fiberqc2) <- "Total fiber intake"
label(gxe_table1$fruitqc2) <- "Total fruit intake"
label(gxe_table1$vegetableqc2) <- "Total vegetable intake"

label(gxe_table1$folate_totqc2) <- "Total folate intake"
label(gxe_table1$folate_dietqc2) <- "Dietary folate intake"

label(gxe_table1$heightcm) <- "Height"
units(gxe_table1$heightcm) <- "cm"

label(gxe_table1$hrt_ref_pm2) <- "Any Post-Menopausal HRT Use"
label(gxe_table1$eo_ref_pm_gxe) <- "Post-menopausal HRT (estrogen only) use"
label(gxe_table1$ep_ref_pm_gxe) <- "Post-menopausal HRT (estrogen+progesterone) use"

label(gxe_table1$asp_ref) <- "Reg aspirin or NSAID use"
label(gxe_table1$aspirin) <- "Reg aspirin use"
label(gxe_table1$nsaids) <- "Ref NSAID use"

label(gxe_table1$procmeatqc2) <- "Total processed meat intake"
label(gxe_table1$redmeatqc2) <- "Total red meat intake"

label(gxe_table1$smk_ever) <- "Smoking, never/ever"
label(gxe_table1$smk_aveday) <- "Smoking, avg cig/day"
label(gxe_table1$smk_pkyr) <- "Smoking, avg pks/yr"

label(gxe_table1$diab) <- "T2D (ever diagnosed)"

saveRDS(gxe_table1, "~/data/FIGI_EpiData_rdata/FIGI_results_exposure_table1.rds")


# data for ggplot (bar and box plots). Setting NA as first factor level
# not excluding subsets of individuals since I want to show missing values by studies
# data also includes variables for counts table
gxe_ggplot <- figi_gwas %>% 
  filter(gxe == 1, EUR_subset == 1) %>%  # important - note it's EUR and gxe set
  mutate(asp_ref = fct_explicit_na(asp_ref, na_level = "NA"), 
         asp_ref = fct_relevel(asp_ref, "NA", "Yes", "No"),
         asp_ref_count = fct_relevel(asp_ref, "No", "Yes", "NA"),
         aspirin = fct_explicit_na(aspirin, na_level = "NA"), 
         aspirin = fct_relevel(aspirin, "NA", "Yes", "No"),
         aspirin_count = fct_relevel(aspirin, "No", "Yes", "NA"),
         nsaids = fct_explicit_na(nsaids, na_level = "NA"), 
         nsaids = fct_relevel(nsaids, "NA", "Yes", "No"),
         nsaids_count = fct_relevel(nsaids, "No", "Yes", "NA"),
         folate_totqc2 = fct_explicit_na(folate_totqc2, na_level = "NA"), 
         folate_totqc2 = fct_relevel(folate_totqc2, "NA", "1", "2", "3", "4"),
         folate_totqc2_count = fct_relevel(folate_totqc2, "1", "2", "3", "4", "NA"),
         folate_dietqc2 = fct_explicit_na(folate_dietqc2, na_level = "NA"), 
         folate_dietqc2 = fct_relevel(folate_dietqc2, "NA", "1", "2", "3", "4"),
         folate_dietqc2_count = fct_relevel(folate_dietqc2, "1", "2", "3", "4", "NA"),
         calcium_totqc2 = fct_explicit_na(calcium_totqc2, na_level = "NA"), 
         calcium_totqc2 = fct_relevel(calcium_totqc2, "NA", "1", "2", "3", "4"),
         calcium_totqc2_count = fct_relevel(calcium_totqc2, "1", "2", "3", "4", "NA"),
         calcium_dietqc2 = fct_explicit_na(calcium_dietqc2, na_level = "NA"), 
         calcium_dietqc2 = fct_relevel(calcium_dietqc2, "NA", "1", "2", "3", "4"),
         calcium_dietqc2_count = fct_relevel(calcium_dietqc2, "1", "2", "3", "4", "NA"),
         alcoholc_heavy = fct_explicit_na(alcoholc_heavy, na_level = "NA"), 
         alcoholc_heavy = fct_relevel(alcoholc_heavy, "NA", ">28g/d", "nondrinker"),
         alcoholc_heavy_count = fct_relevel(alcoholc_heavy, "nondrinker", ">28g/d", "NA"),
         alcoholc_moderate = fct_explicit_na(alcoholc_moderate, na_level = "NA"), 
         alcoholc_moderate = fct_relevel(alcoholc_moderate, "NA", "nondrinker", "1-28g/d"),
         alcoholc_moderate_count = fct_relevel(alcoholc_moderate, "1-28g/d", "nondrinker", "NA"),
         alcoholc_heavy_vs_moderate = fct_explicit_na(alcoholc_heavy_vs_moderate, na_level = "NA"), 
         alcoholc_heavy_vs_moderate = fct_relevel(alcoholc_heavy_vs_moderate, "NA", ">28g/d", "1-28g/d"),
         alcoholc_heavy_vs_moderate_count = fct_relevel(alcoholc_heavy_vs_moderate, "1-28g/d", ">28g/d", "NA"),
         hrt_ref_pm = fct_explicit_na(hrt_ref_pm, na_level = "NA"), 
         hrt_ref_pm = fct_relevel(hrt_ref_pm, "NA", "Yes", "No"),
         hrt_ref_pm_count = fct_relevel(hrt_ref_pm, "No", "Yes", "NA"),
         hrt_ref_pm2 = fct_explicit_na(hrt_ref_pm2, na_level = "NA"), 
         hrt_ref_pm2 = fct_relevel(hrt_ref_pm2, "NA", "Yes", "No"),
         hrt_ref_pm2_count = fct_relevel(hrt_ref_pm2, "No", "Yes", "NA"),
         eo_ref_pm_gxe = fct_explicit_na(eo_ref_pm_gxe, na_level = "NA"), 
         eo_ref_pm_gxe = fct_relevel(eo_ref_pm_gxe, "NA", "Yes", "No"),
         eo_ref_pm_gxe_count = fct_relevel(eo_ref_pm_gxe, "No", "Yes", "NA"),
         ep_ref_pm_gxe = fct_explicit_na(ep_ref_pm_gxe, na_level = "NA"), 
         ep_ref_pm_gxe = fct_relevel(ep_ref_pm_gxe, "NA", "Yes", "No"),
         ep_ref_pm_gxe_count = fct_relevel(ep_ref_pm_gxe, "No", "Yes", "NA"),
         diab = fct_explicit_na(diab, na_level = "NA"), 
         diab = fct_relevel(diab, "NA", "Yes", "No"),
         diab_count = fct_relevel(diab, "No", "Yes", "NA"),
         smk_ever = fct_explicit_na(smk_ever, na_level = "NA"), 
         smk_ever = fct_relevel(smk_ever, "NA", "Yes", "No"),
         smk_ever_count = fct_relevel(smk_ever, "No", "Yes", "NA"),
         fiberqc2 = fct_explicit_na(fiberqc2, na_level = "NA"), 
         fiberqc2 = fct_relevel(fiberqc2, "NA", "1", "2", "3", "4"),
         fiberqc2_count = fct_relevel(fiberqc2, "1", "2", "3", "4", "NA"),
         fruitqc2 = fct_explicit_na(fruitqc2, na_level = "NA"), 
         fruitqc2 = fct_relevel(fruitqc2, "NA", "1", "2", "3", "4"),
         fruitqc2_count = fct_relevel(fruitqc2, "1", "2", "3", "4", "NA"),
         vegetableqc2 = fct_explicit_na(vegetableqc2, na_level = "NA"), 
         vegetableqc2 = fct_relevel(vegetableqc2, "NA", "1", "2", "3", "4"),
         vegetableqc2_count = fct_relevel(vegetableqc2, "1", "2", "3", "4", "NA"),
         redmeatqc2 = fct_explicit_na(redmeatqc2, na_level = "NA"), 
         redmeatqc2 = fct_relevel(redmeatqc2, "NA", "1", "2", "3", "4"),
         redmeatqc2_count = fct_relevel(redmeatqc2, "1", "2", "3", "4", "NA"),
         procmeatqc2 = fct_explicit_na(procmeatqc2, na_level = "NA"), 
         procmeatqc2 = fct_relevel(procmeatqc2, "NA", "1", "2", "3", "4"),
         procmeatqc2_count = fct_relevel(procmeatqc2, "1", "2", "3", "4", "NA"))

saveRDS(gxe_ggplot, "~/data/FIGI_EpiData_rdata/FIGI_results_exposure_ggplot.rds")


# complete case for second plot
# gxe_ggplot_subset <- gxe_ggplot %>%
#   filter(vcfid %in% exposure_subset)
