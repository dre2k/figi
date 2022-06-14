#=============================================================================#
# FIGI GxE physact
# 04/01/2022
#
# Sensitivity analyses for main findings 
# spans multiple exposures so better to have a separate script
#=============================================================================#


# Setup -------------------------------------------------------------------
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
library(broom)
rm(list = ls())

# ---------------------------------------------- #
# variable          finding           method
# ---------------   ---------------   ---------
# active_met_875    15:32994756:T:C   gao 99.5
# active_met_1800   15:33008870:T:C   gao 97.5
# methrswklnsqc2    20:49693755:T:C   GxE
# methrswklns       21:44227260:T:C   case-only
# ---------------------------------------------- #


# input files
hrc_version <- 'v3.1'
covariates <- sort(c('age_ref_imp', 'sex', 'study_gxe', 'energytot_imp', 'pc1', 'pc2', 'pc3'))
pca <- fread("/media/work/gwis_test/PCA/20210222/figi_gxe_pca_update.eigenvec") %>% 
  rename(vcfid = IID) %>% 
  rename_with(tolower) %>% 
  dplyr::select(-`#fid`)

active_met_875_ids <- readRDS(glue("/media/work/gwis_test/active_met_875/data/FIGI_{hrc_version}_gxeset_active_met_875_basic_covars_glm.rds")) %>%
  pull(vcfid)
active_met_1800_ids <- readRDS(glue("/media/work/gwis_test/active_met_1800/data/FIGI_{hrc_version}_gxeset_active_met_1800_basic_covars_glm.rds")) %>%
  pull(vcfid)
methrswklnsqc2_ids <- readRDS(glue("/media/work/gwis_test/methrswklnsqc2/data/FIGI_{hrc_version}_gxeset_methrswklnsqc2_basic_covars_glm.rds")) %>%
  pull(vcfid)
methrswklns_ids <- readRDS(glue("/media/work/gwis_test/methrswklns/data/FIGI_{hrc_version}_gxeset_methrswklns_basic_covars_glm.rds")) %>%
  pull(vcfid)


d1 <- qread("/media/work/gwis_test/active_met_875/output/posthoc/dosage_chr15_32994756.qs")
d2 <- qread("/media/work/gwis_test/active_met_1800/output/posthoc/dosage_chr15_33008870.qs")
d3 <- qread("/media/work/gwis_test/methrswklnsqc2/output/posthoc/dosage_chr20_49693755.qs")
d4 <- qread("/media/work/gwis_test/methrswklns/output/posthoc/dosage_chr21_44227260.qs")

doses <- Reduce(function(x,y) {merge(x,y)}, list(d1,d2,d3,d4))

gxe <- readRDS(glue("/media/work/gwis_test/data/FIGI_{hrc_version}_gxeset_analysis_data_glm.rds")) %>%
  dplyr::select(-starts_with('pc')) %>%
  inner_join(pca, 'vcfid') %>% 
  inner_join(doses, 'vcfid') %>% 
  # change coding for some covariates
  mutate(across(.cols = c(redmeatqc2, procmeatqc2, fiberqc2, fruitqc2, vegetableqc2, calcium_totqc2, folate_totqc2), ~ as.numeric(.x)))





#   mutate(active_met_875 = ifelse(is.na(methrswklns), NA,
#                                  ifelse(methrswklns < 8.75, "No", "Yes")))
# 
# filter(vcfid %in% esubset) %>%



# generate datasets why not -----------------------------------------------

active_met_875 <- gxe %>% 
  filter(vcfid %in% active_met_875_ids) %>% 
  mutate(active_met_875 = ifelse(is.na(methrswklns), NA, 
                                 ifelse(methrswklns < 8.75, "No", "Yes")))

active_met_1800 <- gxe %>% 
  filter(vcfid %in% active_met_1800_ids) %>% 
  mutate(active_met_1800 = ifelse(is.na(methrswklns), NA, 
                                 ifelse(methrswklns < 18.00, "No", "Yes")))

methrswklns <- gxe %>% 
  filter(vcfid %in% methrswklns_ids) 


methrswklnsqc2 <- gxe %>% 
  filter(vcfid %in% methrswklnsqc2_ids) 




# additional adjustment covariates -------------------------------------------------------
# output results as tibbles (whole and filtered)

# quick function
mcovariates = covariates

runmod <- function(exposure, snp, mcovariates, model_label = 'age_sex_energy_study_pc') {
  # out <- tidy(glm(glue("outcome ~ {exposure}*{snp} + {glue_collapse(mcovariates, sep = '+')}"),
  #                 data = get(exposure), family = 'binomial')) %>%
  #   mutate(model = model_label) %>% 
  #   filter(grepl(glue("{exposure}|{snp}"), term))
  
  model <- glm(glue("outcome ~ {exposure}*{snp} + {glue_collapse(mcovariates, sep = '+')}"),
                  data = get(exposure), family = 'binomial')
  nobs <- glance(model) %>% pull(nobs)
  out <- tidy(model) %>% 
    mutate(model = model_label, 
           nobs = nobs) %>% 
    filter(grepl(glue("{exposure}|{snp}"), term))
}





# baseline models 
out1 <- do.call(bind_rows, (list(runmod('active_met_875', 'chr15_32994756_T_C_dose', mcovariates), 
             runmod('active_met_1800', 'chr15_33008870_T_C_dose', mcovariates), 
             runmod('methrswklnsqc2', 'chr20_49693755_T_C_dose', mcovariates))))
# write.table(out1, file = "~/Dropbox/physact_tmp/physact_gxe_results_baseline.txt", quote = F, row.names = F, sep = '\t')

# bmi + smk_ever
covar2 = c(mcovariates, 'smk_ever', 'bmi5')
out2 <- do.call(bind_rows, (list(runmod('active_met_875', 'chr15_32994756_T_C_dose', covar2, 'age_sex_energy_study_pc_bmi_smk'), 
                                 runmod('active_met_1800', 'chr15_33008870_T_C_dose', covar2, 'age_sex_energy_study_pc_bmi_smk'), 
                                 runmod('methrswklnsqc2', 'chr20_49693755_T_C_dose', covar2, 'age_sex_energy_study_pc_bmi_smk'))))
# write.table(out2, file = "~/Dropbox/physact_tmp/physact_gxe_results_baseline_smk_bmi.txt", quote = F, row.names = F, sep = '\t')

# famhx
covar3 = c(mcovariates, 'famhx1')
out3 <- do.call(bind_rows, (list(runmod('active_met_875', 'chr15_32994756_T_C_dose', covar3, 'age_sex_energy_study_pc_famhx1'), 
                                 runmod('active_met_1800', 'chr15_33008870_T_C_dose', covar3, 'age_sex_energy_study_pc_famhx1'), 
                                 runmod('methrswklnsqc2', 'chr20_49693755_T_C_dose', covar3, 'age_sex_energy_study_pc_famhx1'))))
# write.table(out3, file = "~/Dropbox/physact_tmp/physact_gxe_results_baseline_famhx1.txt", quote = F, row.names = F, sep = '\t')


# bmi + smk_ever + famhx
covar4 = c(mcovariates, 'smk_ever', 'bmi5', 'famhx1')
out4 <- do.call(bind_rows, (list(runmod('active_met_875', 'chr15_32994756_T_C_dose', covar4, 'age_sex_energy_study_pc_bmi_smk_famhx1'), 
                                 runmod('active_met_1800', 'chr15_33008870_T_C_dose', covar4, 'age_sex_energy_study_pc_bmi_smk_famhx1'), 
                                 runmod('methrswklnsqc2', 'chr20_49693755_T_C_dose', covar4, 'age_sex_energy_study_pc_bmi_smk_famhx1'))))
# write.table(out4, file = "~/Dropbox/physact_tmp/physact_gxe_results_baseline_smk_bmi_famhx1.txt", quote = F, row.names = F, sep = '\t')


# 'full' covar (let's include all the covariates that jim explored)
covar5 <- c(mcovariates, 'smk_ever', 'bmi5', 'famhx1', 'redmeatqc2', 'procmeatqc2', 'folate_totqc2', 'fiberqc2', 'fruitqc2', 'vegetableqc2', 'calcium_totqc2', 'asp_ref', 'alcoholc')
out5 <- do.call(bind_rows, (list(runmod('active_met_875', 'chr15_32994756_T_C_dose', covar5, 'full'), 
                                 runmod('active_met_1800', 'chr15_33008870_T_C_dose', covar5, 'full'), 
                                 runmod('methrswklnsqc2', 'chr20_49693755_T_C_dose', covar5, 'full'))))
write.table(out5, file = "~/Dropbox/physact_tmp/physact_gxe_results_baseline_smk_bmi_famhx1.txt", quote = F, row.names = F, sep = '\t')

out_final <- bind_rows(out1, out2, out3, out4, out5)
write.csv(out_final,  file = "~/Dropbox/physact_tmp/physact_gxe_results_sensitivity.csv", quote = F, row.names = F)



# exclude clinical trials (WHI_2, WHI_3) ----------------------------------
# just do all the above, excluding the appropriate studies

# quick function
mcovariates = covariates

# exclude studies.. easier to do thatin the function than individually 
runmod2 <- function(exposure, snp, mcovariates, model_label = 'age_sex_energy_study_pc') {
  model <- glm(glue("outcome ~ {exposure}*{snp} + {glue_collapse(mcovariates, sep = '+')}"),
               data = get(exposure) %>% 
                 # important part here
                 filter(!study_gxe %in% c("WHI_2", "WHI_3")), family = 'binomial')
  nobs <- glance(model) %>% pull(nobs)
  out <- tidy(model) %>% 
    mutate(model = model_label, 
           nobs = nobs) %>% 
    filter(grepl(glue("{exposure}|{snp}"), term))
}



# baseline models 
out1 <- do.call(bind_rows, (list(runmod2('active_met_875', 'chr15_32994756_T_C_dose', mcovariates), 
                                 runmod2('active_met_1800', 'chr15_33008870_T_C_dose', mcovariates), 
                                 runmod2('methrswklnsqc2', 'chr20_49693755_T_C_dose', mcovariates))))

# bmi + smk_ever
covar2 = c(mcovariates, 'smk_ever', 'bmi5')
out2 <- do.call(bind_rows, (list(runmod2('active_met_875', 'chr15_32994756_T_C_dose', covar2, 'age_sex_energy_study_pc_bmi_smk'), 
                                 runmod2('active_met_1800', 'chr15_33008870_T_C_dose', covar2, 'age_sex_energy_study_pc_bmi_smk'), 
                                 runmod2('methrswklnsqc2', 'chr20_49693755_T_C_dose', covar2, 'age_sex_energy_study_pc_bmi_smk'))))

# famhx
covar3 = c(mcovariates, 'famhx1')
out3 <- do.call(bind_rows, (list(runmod2('active_met_875', 'chr15_32994756_T_C_dose', covar3, 'age_sex_energy_study_pc_famhx1'), 
                                 runmod2('active_met_1800', 'chr15_33008870_T_C_dose', covar3, 'age_sex_energy_study_pc_famhx1'), 
                                 runmod2('methrswklnsqc2', 'chr20_49693755_T_C_dose', covar3, 'age_sex_energy_study_pc_famhx1'))))


# bmi + smk_ever + famhx
covar4 = c(mcovariates, 'smk_ever', 'bmi5', 'famhx1')
out4 <- do.call(bind_rows, (list(runmod2('active_met_875', 'chr15_32994756_T_C_dose', covar4, 'age_sex_energy_study_pc_bmi_smk_famhx1'), 
                                 runmod2('active_met_1800', 'chr15_33008870_T_C_dose', covar4, 'age_sex_energy_study_pc_bmi_smk_famhx1'), 
                                 runmod2('methrswklnsqc2', 'chr20_49693755_T_C_dose', covar4, 'age_sex_energy_study_pc_bmi_smk_famhx1'))))


# 'full' covar (let's include all the covariates that jim explored)

covar5 <- c(mcovariates, 'smk_ever', 'bmi5', 'famhx1', 'redmeatqc2', 'procmeatqc2', 'folate_totqc2', 'fiberqc2', 'fruitqc2', 'vegetableqc2', 'calcium_totqc2', 'asp_ref', 'alcoholc')
out5 <- do.call(bind_rows, (list(runmod2('active_met_875', 'chr15_32994756_T_C_dose', covar5, 'full'), 
                                 runmod2('active_met_1800', 'chr15_33008870_T_C_dose', covar5, 'full'), 
                                 runmod2('methrswklnsqc2', 'chr20_49693755_T_C_dose', covar5, 'full'))))


out_final <- bind_rows(out1, out2, out3, out4, out5)
write.csv(out_final,  file = "~/Dropbox/physact_tmp/physact_gxe_results_sensitivity_exclude_whi2_whi3.csv", quote = F, row.names = F)




# repeat above, but exclude studies primarily composed of adenomas --------
# remove HPFS_3_AD, HPFS_5_AD, NHS_3_AD, NHS_5_AD
# other studies included, meat intake exclude all adenomas individually (not by study)
# i'll remove IN ADDITION to WHI here

# quick function
mcovariates = covariates

# exclude studies.. easier to do thatin the function than individually 
runmod2 <- function(exposure, snp, mcovariates, model_label = 'age_sex_energy_study_pc') {
  model <- glm(glue("outcome ~ {exposure}*{snp} + {glue_collapse(mcovariates, sep = '+')}"),
               data = get(exposure) %>% 
                 # important part here
                 filter(!study_gxe %in% c("WHI_2", "WHI_3", "HPFS_3_AD", "HPFS_5_AD", "NHS_3_AD", "NHS_5_AD")), family = 'binomial')
  nobs <- glance(model) %>% pull(nobs)
  out <- tidy(model) %>% 
    mutate(model = model_label, 
           nobs = nobs) %>% 
    filter(grepl(glue("{exposure}|{snp}"), term))
}

# baseline models 
out1 <- do.call(bind_rows, (list(runmod2('active_met_875', 'chr15_32994756_T_C_dose', mcovariates), 
                                 runmod2('active_met_1800', 'chr15_33008870_T_C_dose', mcovariates), 
                                 runmod2('methrswklnsqc2', 'chr20_49693755_T_C_dose', mcovariates))))

# bmi + smk_ever
covar2 = c(mcovariates, 'smk_ever', 'bmi5')
out2 <- do.call(bind_rows, (list(runmod2('active_met_875', 'chr15_32994756_T_C_dose', covar2, 'age_sex_energy_study_pc_bmi_smk'), 
                                 runmod2('active_met_1800', 'chr15_33008870_T_C_dose', covar2, 'age_sex_energy_study_pc_bmi_smk'), 
                                 runmod2('methrswklnsqc2', 'chr20_49693755_T_C_dose', covar2, 'age_sex_energy_study_pc_bmi_smk'))))

# famhx
covar3 = c(mcovariates, 'famhx1')
out3 <- do.call(bind_rows, (list(runmod2('active_met_875', 'chr15_32994756_T_C_dose', covar3, 'age_sex_energy_study_pc_famhx1'), 
                                 runmod2('active_met_1800', 'chr15_33008870_T_C_dose', covar3, 'age_sex_energy_study_pc_famhx1'), 
                                 runmod2('methrswklnsqc2', 'chr20_49693755_T_C_dose', covar3, 'age_sex_energy_study_pc_famhx1'))))


# bmi + smk_ever + famhx
covar4 = c(mcovariates, 'smk_ever', 'bmi5', 'famhx1')
out4 <- do.call(bind_rows, (list(runmod2('active_met_875', 'chr15_32994756_T_C_dose', covar4, 'age_sex_energy_study_pc_bmi_smk_famhx1'), 
                                 runmod2('active_met_1800', 'chr15_33008870_T_C_dose', covar4, 'age_sex_energy_study_pc_bmi_smk_famhx1'), 
                                 runmod2('methrswklnsqc2', 'chr20_49693755_T_C_dose', covar4, 'age_sex_energy_study_pc_bmi_smk_famhx1'))))


# 'full' covar (let's include all the covariates that jim explored)

covar5 <- c(mcovariates, 'smk_ever', 'bmi5', 'famhx1', 'redmeatqc2', 'procmeatqc2', 'folate_totqc2', 'fiberqc2', 'fruitqc2', 'vegetableqc2', 'calcium_totqc2', 'asp_ref', 'alcoholc')
out5 <- do.call(bind_rows, (list(runmod2('active_met_875', 'chr15_32994756_T_C_dose', covar5, 'full'), 
                                 runmod2('active_met_1800', 'chr15_33008870_T_C_dose', covar5, 'full'), 
                                 runmod2('methrswklnsqc2', 'chr20_49693755_T_C_dose', covar5, 'full'))))


out_final <- bind_rows(out1, out2, out3, out4, out5)
write.csv(out_final,  file = "~/Dropbox/physact_tmp/physact_gxe_results_sensitivity_exclude_whi2_whi3_hpfs_nhs.csv", quote = F, row.names = F)








# repeat above, but exclude studies primarily composed of adenomas --------
# remove HPFS_3_AD, HPFS_5_AD, NHS_3_AD, NHS_5_AD
# other studies included, meat intake exclude all adenomas individually (not by study)
# i'll remove IN ADDITION to WHI here

# quick function
mcovariates = covariates

# exclude individuals
runmod2 <- function(exposure, snp, mcovariates, model_label = 'age_sex_energy_study_pc') {
  model <- glm(glue("outcome ~ {exposure}*{snp} + {glue_collapse(mcovariates, sep = '+')}"),
               data = get(exposure) %>% 
                 # important part here
                 filter(is.na(adenoma_adv)), family = 'binomial')
  nobs <- glance(model) %>% pull(nobs)
  out <- tidy(model) %>% 
    mutate(model = model_label, 
           nobs = nobs) %>% 
  filter(grepl(glue("{exposure}|{snp}"), term))
}



# baseline models 
out1 <- do.call(bind_rows, (list(runmod2('active_met_875', 'chr15_32994756_T_C_dose', mcovariates), 
                                 runmod2('active_met_1800', 'chr15_33008870_T_C_dose', mcovariates), 
                                 runmod2('methrswklnsqc2', 'chr20_49693755_T_C_dose', mcovariates))))

# bmi + smk_ever
covar2 = c(mcovariates, 'smk_ever', 'bmi5')
out2 <- do.call(bind_rows, (list(runmod2('active_met_875', 'chr15_32994756_T_C_dose', covar2, 'age_sex_energy_study_pc_bmi_smk'), 
                                 runmod2('active_met_1800', 'chr15_33008870_T_C_dose', covar2, 'age_sex_energy_study_pc_bmi_smk'), 
                                 runmod2('methrswklnsqc2', 'chr20_49693755_T_C_dose', covar2, 'age_sex_energy_study_pc_bmi_smk'))))

# famhx
covar3 = c(mcovariates, 'famhx1')
out3 <- do.call(bind_rows, (list(runmod2('active_met_875', 'chr15_32994756_T_C_dose', covar3, 'age_sex_energy_study_pc_famhx1'), 
                                 runmod2('active_met_1800', 'chr15_33008870_T_C_dose', covar3, 'age_sex_energy_study_pc_famhx1'), 
                                 runmod2('methrswklnsqc2', 'chr20_49693755_T_C_dose', covar3, 'age_sex_energy_study_pc_famhx1'))))


# bmi + smk_ever + famhx
covar4 = c(mcovariates, 'smk_ever', 'bmi5', 'famhx1')
out4 <- do.call(bind_rows, (list(runmod2('active_met_875', 'chr15_32994756_T_C_dose', covar4, 'age_sex_energy_study_pc_bmi_smk_famhx1'), 
                                 runmod2('active_met_1800', 'chr15_33008870_T_C_dose', covar4, 'age_sex_energy_study_pc_bmi_smk_famhx1'), 
                                 runmod2('methrswklnsqc2', 'chr20_49693755_T_C_dose', covar4, 'age_sex_energy_study_pc_bmi_smk_famhx1'))))


# 'full' covar (let's include all the covariates that jim explored)
covar5 <- c(mcovariates, 'smk_ever', 'bmi5', 'famhx1', 'redmeatqc2', 'procmeatqc2', 'folate_totqc2', 'fiberqc2', 'fruitqc2', 'vegetableqc2', 'calcium_totqc2', 'asp_ref', 'alcoholc')
out5 <- do.call(bind_rows, (list(runmod2('active_met_875', 'chr15_32994756_T_C_dose', covar5, 'full'), 
                                 runmod2('active_met_1800', 'chr15_33008870_T_C_dose', covar5, 'full'), 
                                 runmod2('methrswklnsqc2', 'chr20_49693755_T_C_dose', covar5, 'full'))))


out_final <- bind_rows(out1, out2, out3, out4, out5)
write.csv(out_final,  file = "~/Dropbox/physact_tmp/physact_gxe_results_sensitivity_exclude_adenoma_samples.csv", quote = F, row.names = F)





# output active_met_875 dataset for convenience ---------------------------


saveRDS(active_met_875, "/media/work/gwis_test/active_met_875/data/tmp.rds")









# =========================================================================== #
# =========================================================================== #
# =========================================================================== #
# =========================================================================== #
# =========================================================================== #
# =========================================================================== #


# additional requests 4/13/2022 ------------------------------------------
# (see notion page)


# ---- exclude height from covar5, rerun
# ---- final model should be age, sex, energy, pc, study, bmi, smk, famhx, 
#      red, proc, folate, fiber, fruit, veg, calcium, NSAID, alcoholc
mcovariates
covar6 <- c(mcovariates, 'smk_ever', 'bmi5', 'famhx1', 'redmeatqc2', 'procmeatqc2', 'folate_totqc2', 'fiberqc2', 'fruitqc2', 'vegetableqc2', 'calcium_totqc2', 'asp_ref', 'alcoholc')
out6 <- do.call(bind_rows, (list(runmod('active_met_875', 'chr15_32994756_T_C_dose', covar6, 'covar6'), 
                                 runmod('active_met_1800', 'chr15_33008870_T_C_dose', covar6, 'covar6'), 
                                 runmod('methrswklnsqc2', 'chr20_49693755_T_C_dose', covar6, 'covar6'))))

out_final <- bind_rows(out1, out2, out3, out4, out5)
write.csv(out6,  file = "~/Dropbox/physact_tmp/physact_gxe_results_sensitivity_20220413_set6.csv", quote = F, row.names = F)




