#=============================================================================#
# create gwis results rmarkdown summaries
#=============================================================================#

# args <- commandArgs(trailingOnly=T)
# exposure <- args[1] # asp_ref
# hrc_version = args[2]
# covariates <- sort(c(args[3:length(args)]))

# rmarkdown::render("~/Dropbox/FIGI/FIGI_code/results/results_report/results.Rmd",
#                   params = list(exposure = exposure,
#                                 hrc_version = hrc_version,
#                                 covariates = covariates),
#                   output_file = paste0("~/Dropbox/FIGI/Results/rmarkdown_reports/", exposure, "_results.html"))


gwis_helper <- function(exposure, hrc_version, covariates) {
  rmarkdown::render('~/git/figi/results.Rmd',
                    params = list(exposure = exposure, 
                                  hrc_version = hrc_version,
                                  covariates = covariates), 
                    output_file = glue('~/Dropbox/FIGI/Results/gwis_{exposure}.html'))
}

posthoc_helper <- function(exposure) {
  rmarkdown::render(glue('/home/rak/git/figi/posthoc_{exposure}.Rmd'),
                    params = list(exposure = exposure), 
                    output_file = glue('~/Dropbox/FIGI/Results/posthoc_{exposure}.html'))
}



# ---- asp_ref ----
gwis_helper(exposure = 'asp_ref', 
            hrc_version = 'v2.4', 
            covariates = c('age_ref_imp', 'sex' , 'pc1', 'pc2', 'pc3', 'study_gxe'))

# ---- aspirin ----
gwis_helper(exposure = 'aspirin', 
            hrc_version = 'v2.4', 
            covariates = c('age_ref_imp', 'sex' , 'pc1', 'pc2', 'pc3', 'study_gxe'))



# ---- alcoholc_moderate ----
gwis_helper(exposure = 'alcoholc_moderate', 
            hrc_version = 'v2.3', 
            covariates = c('age_ref_imp', 'energytot_imp', 'sex' , 'pc1', 'pc2', 'pc3', 'study_gxe'))

# ---- alcoholc_heavy_vs_moderate ----
gwis_helper(exposure = 'alcoholc_heavy_vs_moderate', 
            hrc_version = 'v2.3', 
            covariates = c('age_ref_imp', 'energytot_imp', 'sex' , 'pc1', 'pc2', 'pc3', 'study_gxe'))

# ---- redmeatqc2 ----
gwis_helper(exposure = 'redmeatqc2', 
            hrc_version = 'v2.3', 
            covariates = c('age_ref_imp', 'energytot_imp', 'sex' , 'pc1', 'pc2', 'pc3', 'study_gxe'))

# ---- procmeatqc2 ----
gwis_helper(exposure = 'procmeatqc2', 
            hrc_version = 'v2.3', 
            covariates = c('age_ref_imp', 'energytot_imp', 'sex' , 'pc1', 'pc2', 'pc3', 'study_gxe'))

# ---- calcium_totqc2 ----
gwis_helper(exposure = 'calcium_totqc2', 
            hrc_version = 'v3.0', 
            covariates = c('age_ref_imp', 'energytot_imp', 'sex' , 'pc1', 'pc2', 'pc3', 'study_gxe'))

# ---- calcium_dietqc2 ----
gwis_helper(exposure = 'calcium_dietqc2', 
            hrc_version = 'v3.0', 
            covariates = c('age_ref_imp', 'energytot_imp', 'sex' , 'pc1', 'pc2', 'pc3', 'study_gxe'))
posthoc_report(exposure = 'calcium_dietqc2')



# ---- folate_totqc2 ----
gwis_helper(exposure = 'folate_totqc2', 
            hrc_version = 'v3.0', 
            covariates = c('age_ref_imp', 'energytot_imp', 'sex' , 'pc1', 'pc2', 'pc3', 'study_gxe'))

# ---- folate_dietqc2 ----
gwis_helper(exposure = 'folate_dietqc2', 
            hrc_version = 'v3.0', 
            covariates = c('age_ref_imp', 'energytot_imp', 'sex' , 'pc1', 'pc2', 'pc3', 'study_gxe'))

# ---- fruitqc2 ----
gwis_helper(exposure = 'fruitqc2', 
            hrc_version = 'v2.3', 
            covariates = c('age_ref_imp', 'energytot_imp', 'sex' , 'pc1', 'pc2', 'pc3', 'study_gxe'))

# ---- fiberqc2 ----
gwis_helper(exposure = 'fiberqc2', 
            hrc_version = 'v2.3', 
            covariates = c('age_ref_imp', 'energytot_imp', 'sex' , 'pc1', 'pc2', 'pc3', 'study_gxe'))

# ---- vegetableqc2 ----
gwis_helper(exposure = 'vegetableqc2', 
            hrc_version = 'v2.3', 
            covariates = c('age_ref_imp', 'energytot_imp', 'sex' , 'pc1', 'pc2', 'pc3', 'study_gxe'))

# ---- diab ----
gwis_helper(exposure = 'diab', 
            hrc_version = 'v2.3', 
            covariates = c('age_ref_imp', 'sex' , 'pc1', 'pc2', 'pc3', 'study_gxe'))

# ---- hrt_ref_pm2 ----
gwis_helper(exposure = 'hrt_ref_pm2', 
            hrc_version = 'v2.3', 
            covariates = c('age_ref_imp', 'pc1', 'pc2', 'pc3', 'study_gxe'))


# ---- pure_eo_allNo ----
gwis_helper(exposure = 'pure_eo_allNo', 
            hrc_version = 'v2.3', 
            covariates = c('age_ref_imp', 'pc1', 'pc2', 'pc3', 'study_gxe'))

# ---- pure_eo_all ----
gwis_helper(exposure = 'pure_eo_all', 
            hrc_version = 'v2.3', 
            covariates = c('age_ref_imp', 'pc1', 'pc2', 'pc3', 'study_gxe'))


# ---- pure_ep_allNo ----
gwis_helper(exposure = 'pure_ep_allNo', 
            hrc_version = 'v2.3', 
            covariates = c('age_ref_imp', 'pc1', 'pc2', 'pc3', 'study_gxe'))


# ---- pure_ep_all ----
gwis_helper(exposure = 'pure_ep_all', 
            hrc_version = 'v2.3', 
            covariates = c('age_ref_imp', 'pc1', 'pc2', 'pc3', 'study_gxe'))





# ---- bmi5 ----
gwis_helper(exposure = 'bmi5', 
            hrc_version = 'v2.3', 
            covariates = c('age_ref_imp', 'sex' , 'pc1', 'pc2', 'pc3', 'study_gxe'))

# ---- asp_ref ----
gwis_helper(exposure = 'asp_ref', 
            hrc_version = 'v2.4', 
            covariates = c('age_ref_imp', 'sex' , 'pc1', 'pc2', 'pc3', 'study_gxe'))

# ---- aspirin ---- 
gwis_helper(exposure = 'aspirin', 
            hrc_version = 'v2.4', 
            covariates = c('age_ref_imp', 'sex' , 'pc1', 'pc2', 'pc3', 'study_gxe'))

# ---- hrt_ref_pm2 ---- #
gwis_helper(exposure = 'hrt_ref_pm2', 
            hrc_version = 'v2.3', 
            covariates = c('age_ref_imp', 'pc1', 'pc2', 'pc3', 'study_gxe'))

# ---- ep_ref_pm_gxe ---- #
gwis_helper(exposure = 'ep_ref_pm_gxe', 
            hrc_version = 'v2.3', 
            covariates = c('age_ref_imp', 'pc1', 'pc2', 'pc3', 'study_gxe'))

# ---- eo_ref_pm_gxe ---- #
gwis_helper(exposure = 'eo_ref_pm_gxe', 
            hrc_version = 'v2.3', 
            covariates = c('age_ref_imp', 'pc1', 'pc2', 'pc3', 'study_gxe'))


# ---- folate_totqc2 ---- #
gwis_helper(exposure = 'folate_totqc2', 
            hrc_version = 'v2.3', 
            covariates = c('age_ref_imp', 'sex', 'energytot_imp', 'pc1', 'pc2', 'pc3', 'study_gxe'))

# ---- folate_dietqc2 ---- #
gwis_helper(exposure = 'folate_dietqc2', 
            hrc_version = 'v2.3', 
            covariates = c('age_ref_imp', 'sex', 'energytot_imp', 'pc1', 'pc2', 'pc3', 'study_gxe'))


# ---- calcium_totqc2 ---- #
gwis_helper(exposure = 'calcium_totqc2', 
            hrc_version = 'v2.3', 
            covariates = c('age_ref_imp', 'sex', 'energytot_imp', 'pc1', 'pc2', 'pc3', 'study_gxe'))

# ---- calcium_dietqc2 ---- #
gwis_helper(exposure = 'calcium_dietqc2', 
            hrc_version = 'v2.3', 
            covariates = c('age_ref_imp', 'sex', 'energytot_imp', 'pc1', 'pc2', 'pc3', 'study_gxe'))


# ---- height10 ---- #
gwis_helper(exposure = 'height10', 
            hrc_version = 'v2.3', 
            covariates = c('age_ref_imp', 'sex', 'pc1', 'pc2', 'pc3', 'study_gxe'))


# ---- height10_cohort ---- #
gwis_helper(exposure = 'height10_cohort', 
            hrc_version = 'v2.3', 
            covariates = c('age_ref_imp', 'sex', 'pc1', 'pc2', 'pc3', 'study_gxe'))


# ---- height10_casecontrol ---- #
gwis_helper(exposure = 'height10_casecontrol', 
            hrc_version = 'v2.3', 
            covariates = c('age_ref_imp', 'sex', 'pc1', 'pc2', 'pc3', 'study_gxe'))





# =========================================================================== #
# main effects only report ----
# =========================================================================== #

gwis_helper <- function(exposure, hrc_version, covariates) {
  rmarkdown::render('~/Dropbox/FIGI/FIGI_code/results/gwis/results_main_effects.Rmd',
                    params = list(exposure = exposure, 
                                  hrc_version = hrc_version,
                                  covariates = covariates), 
                    output_file = paste0('~/Dropbox/FIGI/FIGI_code/results/', exposure, '/', 'main_effects_', exposure, '.html'))
}

# ---- pure_eo_allNo ----
gwis_helper(exposure = 'pure_eo_allNo', 
            hrc_version = 'v2.3', 
            covariates = c('age_ref_imp', 'pc1', 'pc2', 'pc3', 'study_gxe'))

# ---- pure_eo_allNo ----
gwis_helper(exposure = 'pure_ep_allNo', 
            hrc_version = 'v2.3', 
            covariates = c('age_ref_imp', 'pc1', 'pc2', 'pc3', 'study_gxe'))

# ---- redmeatqc2 ----
gwis_helper(exposure = 'redmeatqc2', 
            hrc_version = 'v2.3', 
            covariates = c('age_ref_imp', 'pc1', 'pc2', 'pc3', 'sex', 'study_gxe', "energytot_imp"))

# ---- procmeatqc2 ----
gwis_helper(exposure = 'procmeatqc2', 
            hrc_version = 'v2.3', 
            covariates = c('age_ref_imp', 'pc1', 'pc2', 'pc3', 'sex', 'study_gxe', "energytot_imp"))

























# rmarkdown_wrapper <- function(exposure, covariates, hrc_version = "v2.3") {
#   rmarkdown::render("~/git/FIGI_code/results/results_report/results_report_parent.Rmd", 
#                     params = list(exposure = exposure,
#                                   covariates = covariates), 
#                     output_file = paste0("~/Dropbox/FIGI/Results/rmarkdown_reports/", exposure, "_results.html"))
# }
# covariates <- c('age_ref_imp', 'sex', 'study_gxe', 'pc1', 'pc2', 'pc3')
# 
# 
# rmarkdown_wrapper(exposure = 'diab',
#                   covariates = covariates)
# 
# 
# rmarkdown_wrapper(exposure = 'alcoholc_heavy_vs_moderate',
#                   covariates = covariates)
# 
# 
# 
# 
# rmarkdown_wrapper(exposure_var = 'alcoholc_moderate')
# rmarkdown_wrapper(exposure_var = 'asp_ref')
# rmarkdown_wrapper(exposure_var = 'aspirin')
# rmarkdown_wrapper(exposure_var = 'bmi5')
# rmarkdown_wrapper(exposure_var = 'diab')
# rmarkdown_wrapper(exposure_var = 'folate_totqc2')
# rmarkdown_wrapper(exposure_var = 'folate_dietqc2')
# 
# rmarkdown_wrapper(exposure_var = 'calcium_totqc2')
# rmarkdown_wrapper(exposure_var = 'calcium_dietqc2')
# 
# rmarkdown_wrapper(exposure_var = 'fruitqc2')
# rmarkdown_wrapper(exposure_var = 'fiberqc2')
# rmarkdown_wrapper(exposure_var = 'vegetableqc2')
# 
# rmarkdown_wrapper(exposure_var = 'redmeatqc2')
# rmarkdown_wrapper(exposure_var = 'procmeatqc2')
# 
# # rmarkdown_wrapper(exposure_var = 'hrt_ref_pm2')
# 
# rmarkdown_wrapper(exposure_var = 'smk_ever')
# rmarkdown_wrapper(exposure_var = 'smk_pkyr')
# rmarkdown_wrapper(exposure_var = 'smk_aveday')
# 
# rmarkdown_wrapper(exposure_var = 'height10')
# rmarkdown_wrapper(exposure_var = 'eo_ref_pm_gxe')
# rmarkdown_wrapper(exposure_var = 'ep_ref_pm_gxe')
# 
# 
# # not yet (generates error/warning not sure what's going on )
# rmarkdown_wrapper(exposure_var = 'height10')
# 
# 
# 
# 
# # ---- NSAIDS ---- 
# # asp_ref
# rm(list = ls())
# rmarkdown::render("~/git/FIGI_code/results/results_report_parent.Rmd", params = list(
#   exposure = 'asp_ref',
#   additional_analyses = F,
#   covariates_suffix = 'age_ref_imp_sex_study_gxe_pc1_pc2_pc3'
# ), output_file = "~/Dropbox/FIGI/Results/asp_ref_results.html")
# 
# 
# # aspirin
# rm(list = ls())
# rmarkdown::render("~/git/FIGI_code/results/results_report_parent.Rmd", params = list(
#   exposure = 'aspirin',
#   additional_analyses = F,
#   covariates_suffix = 'age_ref_imp_sex_study_gxe_pc1_pc2_pc3'
# ), output_file = "~/Dropbox/FIGI/Results/aspirin_results.html")
# 
# 
# # nsaids
# rm(list = ls())
# rmarkdown::render("~/git/FIGI_code/results/results_report_parent.Rmd", params = list(
#   exposure = 'nsaids',
#   is_exposure_categorical = T,
#   energy_adj = F,
#   table1_by_e = T, 
#   additional_analyses = F,
#   covariates_suffix = 'age_ref_imp_sex_study_gxe_pc1_pc2_pc3'
# ), output_file = "~/Dropbox/FIGI/Results/nsaids_results.html")
# 
# 
# 
# 
# 
# # ---- Smoking ---- 
# # smk_ever
# rm(list = ls())
# rmarkdown::render("~/git/FIGI_code/results/results_report_parent.Rmd", params = list(
#   exposure = 'smk_ever',
#   is_exposure_categorical = T,
#   energy_adj = F,
#   table1_by_e = T, 
#   additional_analyses = F,
#   covariates_suffix = 'age_ref_imp_sex_study_gxe_pc1_pc2_pc3'
# ), output_file = "~/Dropbox/FIGI/Results/smk_ever_results.html")
# 
# 
# # smk_aveday
# rm(list = ls())
# rmarkdown::render("~/git/FIGI_code/results/results_report_parent.Rmd", params = list(
#   exposure = 'smk_aveday',
#   is_exposure_categorical = F,
#   energy_adj = F,
#   table1_by_e = F, 
#   additional_analyses = F,
#   covariates_suffix = 'age_ref_imp_sex_study_gxe_pc1_pc2_pc3'
# ), output_file = "~/Dropbox/FIGI/Results/smk_aveday_results.html")
# 
# 
# # smk_pkyr
# rm(list = ls())
# rmarkdown::render("~/git/FIGI_code/results/results_report_parent.Rmd", params = list(
#   exposure = 'smk_pkyr',
#   is_exposure_categorical = F,
#   energy_adj = F,
#   table1_by_e = F, 
#   additional_analyses = F,
#   covariates_suffix = 'age_ref_imp_sex_study_gxe_pc1_pc2_pc3'
# ), output_file = "~/Dropbox/FIGI/Results/smk_pkyr_results.html")
# 
# 
# 
# # ---- Red and Processed meat ----
# # redmeatqc2
# rm(list = ls())
# rmarkdown::render("~/git/FIGI_code/results/results_report_parent.Rmd", params = list(
#   exposure = 'redmeatqc2',
#   is_exposure_categorical = T,
#   energy_adj = T,
#   table1_by_e = F, 
#   additional_analyses = F,
#   covariates_suffix = 'age_ref_imp_sex_energytot_imp_study_gxe_pc1_pc2_pc3'
# ), output_file = "~/Dropbox/FIGI/Results/redmeatqc2_results.html")
# 
# 
# # procmeatqc2
# rm(list = ls())
# rmarkdown::render("~/git/FIGI_code/results/results_report_parent.Rmd", params = list(
#   exposure = 'procmeatqc2',
#   is_exposure_categorical = T,
#   energy_adj = T,
#   table1_by_e = F, 
#   additional_analyses = F,
#   covariates_suffix = 'age_ref_imp_sex_energytot_imp_study_gxe_pc1_pc2_pc3'
# ), output_file = "~/Dropbox/FIGI/Results/procmeatqc2_results.html")
# 
# 
# 
# 
# # ---- Fruits, vegetables, fiber ----
# 
# # fruitqc2
# rm(list = ls())
# rmarkdown::render("~/git/FIGI_code/results/results_report_parent.Rmd", params = list(
#   exposure = 'fruitqc2',
#   is_exposure_categorical = T,
#   energy_adj = T,
#   table1_by_e = F, 
#   additional_analyses = F,
#   covariates_suffix = 'age_ref_imp_sex_energytot_imp_study_gxe_pc1_pc2_pc3'
# ), output_file = "~/Dropbox/FIGI/Results/fruitqc2_results.html")
# 
# 
# # fiberqc2
# rm(list = ls())
# rmarkdown::render("~/git/FIGI_code/results/results_report_parent.Rmd", params = list(
#   exposure = 'fiberqc2',
#   is_exposure_categorical = T,
#   energy_adj = T,
#   table1_by_e = F, 
#   additional_analyses = F,
#   covariates_suffix = 'age_ref_imp_sex_energytot_imp_study_gxe_pc1_pc2_pc3'
# ), output_file = "~/Dropbox/FIGI/Results/fiberqc2_results.html")
# 
# 
# # vegetableqc2
# rm(list = ls())
# rmarkdown::render("~/git/FIGI_code/results/results_report_parent.Rmd", params = list(
#   exposure = 'vegetableqc2',
#   is_exposure_categorical = T,
#   energy_adj = T,
#   table1_by_e = F, 
#   additional_analyses = F,
#   covariates_suffix = 'age_ref_imp_sex_energytot_imp_study_gxe_pc1_pc2_pc3'
# ), output_file = "~/Dropbox/FIGI/Results/vegetableqc2_results.html")
# 
# 
# 
# # ---- HRT ---- 
# # NO SEX STRATIFIED
# # hrt_ref_pm
# rmarkdown::render("~/git/FIGI_code/results/results_report_parent.Rmd", params = list(
#   exposure = 'hrt_ref_pm',
#   is_exposure_categorical = T,
#   energy_adj = F,
#   covariates_suffix = 'age_ref_imp_sex_study_gxe_pc1_pc2_pc3'
# ), output_file = "~/Dropbox/FIGI/Results/hrt_ref_pm_results.html")
# 
# 
# # hrt_ref_pm2
# rm(list = ls())
# rmarkdown::render("~/git/FIGI_code/results/results_report_parent.Rmd", params = list(
#   exposure = 'hrt_ref_pm2',
#   is_exposure_categorical = T,
#   energy_adj = F,
#   table1_by_e = T, 
#   additional_analyses = F,
#   covariates_suffix = 'age_ref_imp_study_gxe_pc1_pc2_pc3'
# ), output_file = "~/Dropbox/FIGI/Results/hrt_ref_pm2_results.html")
# 
# 
# # eo_ref_pm_gxe
# rmarkdown::render("~/git/FIGI_code/results/results_report_parent.Rmd", params = list(
#   exposure = 'eo_ref_pm_gxe',
#   is_exposure_categorical = T,
#   energy_adj = F,
#   covariates_suffix = 'age_ref_imp_sex_study_gxe_pc1_pc2_pc3'
# ), output_file = "~/Dropbox/FIGI/Results/eo_ref_pm_gxe_results.html")
# 
# # ep_ref_pm_gxe
# rmarkdown::render("~/git/FIGI_code/results/results_report_parent.Rmd", params = list(
#   exposure = 'ep_ref_pm_gxe',
#   is_exposure_categorical = T,
#   energy_adj = F,
#   covariates_suffix = 'age_ref_imp_sex_study_gxe_pc1_pc2_pc3'
# ), output_file = "~/Dropbox/FIGI/Results/ep_ref_pm_gxe_results.html")
# 
# 
# # ---- T2D ---- #
# # diab
# rmarkdown::render("~/git/FIGI_code/results/results_report_parent.Rmd", params = list(
#   exposure = 'diab',
#   is_exposure_categorical = T,
#   energy_adj = F,
#   covariates_suffix = 'age_ref_imp_sex_study_gxe_pc1_pc2_pc3'
# ), output_file = "~/Dropbox/FIGI/Results/diab_results.html")
# 
# 
# # ---- Folate ----
# # folate_totqc2
# rm(list = ls())
# rmarkdown::render("~/git/FIGI_code/results/results_report_parent.Rmd", 
#            params = list(exposure = 'folate_totqc2'), 
#            output_file = "~/Dropbox/FIGI/Results/folate_totqc2_results.html")
# 
# 
# # folate_dietqc2
# rm(list = ls())
# rmarkdown::render("~/git/FIGI_code/results/results_report_parent.Rmd", params = list(
#   exposure = 'folate_dietqc2',
#   is_exposure_categorical = T,
#   energy_adj = T,
#   table1_by_e = F, 
#   covariates_suffix = 'age_ref_imp_sex_energytot_imp_study_gxe_pc1_pc2_pc3'
# ), output_file = "~/Dropbox/FIGI/Results/folate_dietqc2_results.html")
