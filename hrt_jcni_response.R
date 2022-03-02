#=============================================================================#
# FIGI HRT
# Reviewer comment responses
#
# 10/28/2021
#=============================================================================#



#-----------------------------------------------------------------------------#
# ------ MHT ------ 
# Create forest plots of GxE term for main findings
#
# MHT x SNPs
# 
# rs117868593 - chr12_13670508_G_C
# rs10782186 -  chr6_117823508_T_C
#-----------------------------------------------------------------------------#

rm(list = ls())
exposure = 'hrt_ref_pm2'
hrc_version = 'v2.3'
covariates <- sort(c('age_ref_imp', 'study_gxe', 'pc1', 'pc2', 'pc3'))

# exposure subset
# re-do continent so there's only North America, Australia, Europe
# try keeping CCFR_1 with Australia for now - might more fine tuned separation? 
# combined Israel with Europe .. 
exposure_subset <- readRDS(glue("/media/work/gwis/data/FIGI_EpiData/FIGI_{hrc_version}_gxeset_{exposure}_basic_covars_glm.rds"))[,'vcfid']
input_data <- readRDS(glue("/media/work/gwis/data/FIGI_EpiData/FIGI_{hrc_version}_gxeset_analysis_data_glm.rds")) %>% 
  filter(vcfid %in% exposure_subset) %>% 
  mutate(
    study_continent = fct_drop(fct_recode(study_continent, Europe = "Israel", Australia = "Australia/America", "North America" = "America")), 
    study_continent = replace(study_continent, study_gxe == "CCFR_1", "North America")) # sometimes base is best

study_design <- input_data %>% 
  dplyr::select(study_gxe, study_design) %>% 
  filter(!duplicated(.)) %>% 
  arrange(study_gxe)

study_continent <- input_data %>% 
  dplyr::select(study_gxe, study_continent) %>%
  filter(!duplicated(.)) %>% 
  arrange(study_gxe)







# ----- MHT (main effects by region) -----

# forest plot
covariates <- sort(c('age_ref_imp', 'pc1', 'pc2', 'pc3'))
model_formula <- glue("outcome ~ hrt_ref_pm2 + {glue_collapse(covariates, sep = '+')}")

glm_out <- input_data %>% 
  tidyr::nest(data = -study_gxe) %>% 
  dplyr::mutate(fit = purrr::map(data, ~glm(model_formula, data = .x, family = "binomial")), 
                tidied = purrr::map(fit, ~tidy(.x)), quality = purrr::map(fit, ~glance(.x))) %>% 
  dplyr::select(-data, -fit) %>% 
  tidyr::unnest(tidied) %>% 
  tidyr::unnest(quality) %>% 
  # dplyr::filter(grepl(":", term),
  #               null.deviance > quantile(null.deviance, 0.05)) %>%
  dplyr::filter(grepl(exposure, term)) %>%
  dplyr::arrange(study_gxe) %>% 
  left_join(study_continent, 'study_gxe')


# --- sanity check - make sure that the study specific values are correct --- #
# dachs1 <- filter(out, study_gxe == "DACHS_1")
# sanity <- glm(model_formula, data = dachs1, family = 'binomial')
# summary(sanity)

results_meta <- meta::metagen(estimate, std.error, data = glm_out, 
                              studlab = paste(study_gxe), 
                              method.tau = 'SJ',
                              hakn = T,
                              sm = "OR",
                              fixed = TRUE, random = FALSE,
                              subgroup = study_continent)

png(glue("~/Dropbox/Working/forest_plot_{exposure}_{hrc_version}_basic_covar_region.png"), 
    height = 13, 
    width = 10, 
    units = "in", res = 300)
meta::forest(results_meta, 
             leftcols = c('studlab', 'TE', 'seTE', 'effect', 'pval', 'w.random'),
             xlim = c(0.25, 4),
             layout = "JAMA",
             test.overall.fixed = TRUE,
             test.overall.random = FALSE)
grid.text("Any MHT by study region", 0.5, 0.98, gp = gpar(cex = 2))
dev.off()



# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

# ----- MHT x rs117868593 (original) -----
# (do you have to recode anything? SNP, or exposure?)
# (according to study tables, not for this SNP)

# Replicate table statistics for these two SNPs (MHT)
dose <- qread("/media/work/gwis_test/hrt_ref_pm2/output/posthoc/dosage_chr12_13670508.qs")
out <- inner_join(input_data, dose, 'vcfid')

# GxE
model_base <- glm(glue("outcome ~ hrt_ref_pm2 + chr12_13670508_G_C_dose + {glue_collapse(covariates, sep = '+')}"), data = out, family = 'binomial')
model_int <-  glm(glue("outcome ~ hrt_ref_pm2 * chr12_13670508_G_C_dose + {glue_collapse(covariates, sep = '+')}"), data = out, family = 'binomial')
lrtest(model_base, model_int) # OK


# forest plot
covariates_forest <- sort(c('age_ref_imp', 'pc1', 'pc2', 'pc3'))
model_formula <- glue("outcome ~ hrt_ref_pm2 * chr12_13670508_G_C_dose + {glue_collapse(covariates_forest, sep = '+')}")

glm_out <- out %>% 
  tidyr::nest(data = -study_gxe) %>% 
  dplyr::mutate(fit = purrr::map(data, ~glm(model_formula, data = .x, family = "binomial")), 
                tidied = purrr::map(fit, ~tidy(.x)), quality = purrr::map(fit, ~glance(.x))) %>% 
  dplyr::select(-data, -fit) %>% 
  tidyr::unnest(tidied) %>% 
  tidyr::unnest(quality) %>% 
  # dplyr::filter(grepl(":", term),
  #               null.deviance > quantile(null.deviance, 0.05)) %>%
  dplyr::filter(grepl(":", term)) %>%
  dplyr::arrange(study_gxe) %>% 
  left_join(study_design, 'study_gxe')

results_meta <- meta::metagen(estimate, std.error, data = glm_out, 
                              studlab = paste(study_gxe), 
                              method.tau = 'SJ',
                              hakn = T,
                              sm = "OR",
                              fixed = TRUE, random = FALSE,
                              subgroup = study_design)

png(glue("~/Dropbox/Working/forest_plot_{exposure}_x_rs117868593_{hrc_version}_basic_covar.png"), 
    height = 13, 
    width = 10, 
    units = "in", res = 300)
meta::forest(results_meta, 
             leftcols = c('studlab', 'TE', 'seTE', 'effect', 'pval', 'w.random'),
             xlim = c(0.1, 10),
             layout = "JAMA",
             test.overall.fixed = TRUE,
             test.overall.random = FALSE)
grid.text("Any MHT x rs117868593", 0.5, 0.98, gp = gpar(cex = 2))
dev.off()




# ----- MHT x rs117868593 (by study COUNTRY) -----


# we know there's no heterogeneity by study design, so we can investigate country now
glm_out <- out %>% 
  tidyr::nest(data = -study_gxe) %>% 
  dplyr::mutate(fit = purrr::map(data, ~glm(model_formula, data = .x, family = "binomial")), 
                tidied = purrr::map(fit, ~tidy(.x)), quality = purrr::map(fit, ~glance(.x))) %>% 
  dplyr::select(-data, -fit) %>% 
  tidyr::unnest(tidied) %>% 
  tidyr::unnest(quality) %>% 
  dplyr::filter(grepl(":", term)) %>%
  dplyr::arrange(study_gxe) %>% 
  left_join(study_continent, 'study_gxe')

results_meta <- meta::metagen(estimate, std.error, data = glm_out, 
                              studlab = paste(study_gxe), 
                              method.tau = 'SJ',
                              hakn = T,
                              sm = "OR",
                              fixed = TRUE, random = FALSE,
                              subgroup = study_continent)

png(glue("~/Dropbox/Working/forest_plot_{exposure}_x_rs117868593_{hrc_version}_basic_covar_region.png"), 
    height = 13, 
    width = 10, 
    units = "in", res = 300)
meta::forest(results_meta, 
             leftcols = c('studlab', 'TE', 'seTE', 'effect', 'pval', 'w.random'),
             xlim = c(0.1, 10),
             layout = "JAMA",
             test.overall.fixed = TRUE,
             test.overall.random = FALSE)
grid.text("Any MHT x rs117868593 by study continent", 0.5, 0.98, gp = gpar(cex = 2))
dev.off()



# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #



# ----- MHT x rs10782186 (pending) -----
# need to ask collaborators (this Locus is newly identified as GWAS in the latest AsianUK Gwas ... )
# (Proceed per Riki email)

# Replicate table statistics for these two SNPs (MHT)
dose <- qread("/media/work/gwis_test/hrt_ref_pm2/output/posthoc/dosage_chr6_117823508.qs")
out <- inner_join(input_data, dose, 'vcfid')

# 2DF identified (replicate statistic here for check)
model_base <- glm(glue("outcome ~ hrt_ref_pm2 +                           {glue_collapse(covariates, sep = '+')}"), data = out, family = 'binomial')
model_int <-  glm(glue("outcome ~ hrt_ref_pm2 * chr6_117823508_T_C_dose + {glue_collapse(covariates, sep = '+')}"), data = out, family = 'binomial')
lrtest(model_base, model_int) # OK

# GxE component .. identified (replicate statistic here for check)
model_base <- glm(glue("outcome ~ hrt_ref_pm2 + chr6_117823508_T_C_dose + {glue_collapse(covariates, sep = '+')}"), data = out, family = 'binomial')
model_int <-  glm(glue("outcome ~ hrt_ref_pm2 * chr6_117823508_T_C_dose + {glue_collapse(covariates, sep = '+')}"), data = out, family = 'binomial')
lrtest(model_base, model_int) # OK



# ----- MHT x rs10782186 (original) -----
# (do you have to recode anything? SNP, or exposure?)
# (according to study tables, not for this SNP)

# forest plot
covariates_forest <- sort(c('age_ref_imp', 'pc1', 'pc2', 'pc3'))
model_formula <- glue("outcome ~ hrt_ref_pm2 * chr6_117823508_T_C_dose + {glue_collapse(covariates_forest, sep = '+')}")

glm_out <- out %>% 
  tidyr::nest(data = -study_gxe) %>% 
  dplyr::mutate(fit = purrr::map(data, ~glm(model_formula, data = .x, family = "binomial")), 
                tidied = purrr::map(fit, ~tidy(.x)), quality = purrr::map(fit, ~glance(.x))) %>% 
  dplyr::select(-data, -fit) %>% 
  tidyr::unnest(tidied) %>% 
  tidyr::unnest(quality) %>% 
  # dplyr::filter(grepl(":", term),
  #               null.deviance > quantile(null.deviance, 0.05)) %>%
  dplyr::filter(grepl(":", term)) %>%
  dplyr::arrange(study_gxe) %>% 
  left_join(study_design, 'study_gxe')

results_meta <- meta::metagen(estimate, std.error, data = glm_out, 
                              studlab = paste(study_gxe), 
                              method.tau = 'SJ',
                              hakn = T,
                              sm = "OR",
                              fixed = TRUE, random = FALSE,
                              subgroup = study_design)

png(glue("~/Dropbox/Working/forest_plot_{exposure}_x_rs10782186_{hrc_version}_basic_covar.png"), 
    height = 13, 
    width = 10, 
    units = "in", res = 300)
meta::forest(results_meta, 
             leftcols = c('studlab', 'TE', 'seTE', 'effect', 'pval', 'w.random'),
             xlim = c(0.1, 10),
             layout = "JAMA",
             test.overall.fixed = TRUE,
             test.overall.random = FALSE)
grid.text("Any MHT x rs10782186", 0.5, 0.98, gp = gpar(cex = 2))
dev.off()




# ----- MHT x rs10782186 (by study COUNTRY) -----


# we know there's no heterogeneity by study design, so we can investigate country now
glm_out <- out %>% 
  tidyr::nest(data = -study_gxe) %>% 
  dplyr::mutate(fit = purrr::map(data, ~glm(model_formula, data = .x, family = "binomial")), 
                tidied = purrr::map(fit, ~tidy(.x)), quality = purrr::map(fit, ~glance(.x))) %>% 
  dplyr::select(-data, -fit) %>% 
  tidyr::unnest(tidied) %>% 
  tidyr::unnest(quality) %>% 
  dplyr::filter(grepl(":", term)) %>%
  dplyr::arrange(study_gxe) %>% 
  left_join(study_continent, 'study_gxe')

results_meta <- meta::metagen(estimate, std.error, data = glm_out, 
                              studlab = paste(study_gxe), 
                              method.tau = 'SJ',
                              hakn = T,
                              sm = "OR",
                              fixed = TRUE, random = FALSE,
                              subgroup = study_continent)

png(glue("~/Dropbox/Working/forest_plot_{exposure}_x_rs10782186_{hrc_version}_basic_covar_region.png"), 
    height = 13, 
    width = 10, 
    units = "in", res = 300)
meta::forest(results_meta, 
             leftcols = c('studlab', 'TE', 'seTE', 'effect', 'pval', 'w.random'),
             xlim = c(0.1, 10),
             layout = "JAMA",
             test.overall.fixed = TRUE,
             test.overall.random = FALSE)
grid.text("Any MHT x rs10782186 by study continent", 0.5, 0.98, gp = gpar(cex = 2))
dev.off()









#-----------------------------------------------------------------------------#
# ------ ESTROGEN + PROGESTERONE (NEW DEF) ------ 
# Create forest plots of GxE term for main findings
#
# pure_ep_all x SNPs
# 
# rs117868593 - chr12_13670508_G_C
# rs10782186 -  chr6_117823508_T_C **
#-----------------------------------------------------------------------------#

exposure = 'pure_ep_all'
# exposure = 'pure_eo_allNo'

hrc_version = 'v2.3'
output_dir = paste0("/media/work/gwis/posthoc/", exposure, "/")
covariates <- sort(c('age_ref_imp', 'study_gxe', 'pc1', 'pc2', 'pc3'))

# input data
# use Yi's data subsets for pure_eo/ep variables
figi_gxe <- readRDS(glue("/media/work/gwis/data/FIGI_EpiData/FIGI_{hrc_version}_gxeset_analysis_data_glm.rds"))
load("/media/work/gwis/data/FIGI_EpiData/Data_HRT_USC_epi_v2.3-20201204.Rdata")

table(data.usc$pure_ep_all)
table(data.usc$pure_ep_allNo)
table(data.usc$pure_ep_all, data.usc$pure_ep_allNo) # flipped coding

pure_eo_all_yi <- data.usc %>% 
  filter(!is.na(pure_ep_all))

# discrepancy is due to case-only studies being removed from GxE subset
input_data <- inner_join(figi_gxe, pure_eo_all_yi[, c('vcfid', 'pure_ep_all', 'pure_ep_allNo')], 'vcfid') %>% 
  mutate(pure_ep_allNo = as.numeric(pure_ep_allNo), 
         pure_ep_all = as.numeric(pure_ep_all))

# drop studies that have zero cell counts when tabulating outcome/pure_eo_allNo
drops <- data.frame(table(input_data$outcome, input_data[, exposure], 
                          input_data$study_gxe)) %>% 
  filter(Freq <= 1)

# sample sizes OK (N = 6887)
input_data <- filter(input_data, !study_gxe %in% unique(drops$Var3)) %>% 
  mutate(study_gxe = fct_drop(study_gxe))







# forest plot
covariates_forest <- sort(c('age_ref_imp', 'pc1', 'pc2', 'pc3'))
model_formula <- glue("outcome ~ pure_ep_all + {glue_collapse(covariates_forest, sep = '+')}")

glm_out <- input_data %>% 
  tidyr::nest(data = -study_gxe) %>% 
  dplyr::mutate(fit = purrr::map(data, ~glm(model_formula, data = .x, family = "binomial")), 
                tidied = purrr::map(fit, ~tidy(.x)), quality = purrr::map(fit, ~glance(.x))) %>% 
  dplyr::select(-data, -fit) %>% 
  tidyr::unnest(tidied) %>% 
  tidyr::unnest(quality) %>% 
  dplyr::filter(grepl(exposure, term)) %>%
  dplyr::arrange(study_gxe) %>% 
  left_join(study_continent, 'study_gxe')

results_meta <- meta::metagen(estimate, std.error, data = glm_out, 
                              studlab = paste(study_gxe), 
                              method.tau = 'SJ',
                              hakn = T,
                              sm = "OR",
                              fixed = TRUE, random = FALSE,
                              subgroup = study_continent)

png(glue("~/Dropbox/Working/forest_plot_{exposure}_{hrc_version}_basic_covar_region.png"), 
    height = 13, 
    width = 10, 
    units = "in", res = 300)
meta::forest(results_meta, 
             leftcols = c('studlab', 'TE', 'seTE', 'effect', 'pval', 'w.random'),
             xlim = c(0.25, 4),
             layout = "JAMA",
             test.overall.fixed = TRUE,
             test.overall.random = FALSE)
grid.text("pure_ep_all by study region", 0.5, 0.98, gp = gpar(cex = 2))
dev.off()





# ----- pure_ep_all x rs117868593 (original) -----

# super sanity check - recreate forest plot from publication
# it's right, use pure_ep_all
# move onto real results now

# Replicate table statistics for these two SNPs (MHT)
dose <- qread("/media/work/gwis_test/hrt_ref_pm2/output/posthoc/dosage_chr12_13670508.qs")
out <- inner_join(input_data, dose, 'vcfid')

# GxE
# covariates <- sort(c('age_ref_imp', 'study_gxe', 'pc1', 'pc2', 'pc3'))
# model_base <- glm(glue("outcome ~ pure_ep_all + chr12_13670508_G_C_dose + {glue_collapse(covariates, sep = '+')}"), data = out, family = 'binomial')
# model_int <-  glm(glue("outcome ~ pure_ep_all * chr12_13670508_G_C_dose + {glue_collapse(covariates, sep = '+')}"), data = out, family = 'binomial')
# lrtest(model_base, model_int) # OK



# forest plot
covariates_forest <- sort(c('age_ref_imp', 'pc1', 'pc2', 'pc3'))
model_formula <- glue("outcome ~ pure_ep_all * chr12_13670508_G_C_dose + {glue_collapse(covariates_forest, sep = '+')}")

glm_out <- out %>% 
  tidyr::nest(data = -study_gxe) %>% 
  dplyr::mutate(fit = purrr::map(data, ~glm(model_formula, data = .x, family = "binomial")), 
                tidied = purrr::map(fit, ~tidy(.x)), quality = purrr::map(fit, ~glance(.x))) %>% 
  dplyr::select(-data, -fit) %>% 
  tidyr::unnest(tidied) %>% 
  tidyr::unnest(quality) %>% 
  dplyr::filter(grepl(':', term)) %>%
  dplyr::arrange(study_gxe) %>% 
  left_join(study_design, 'study_gxe')

results_meta <- meta::metagen(estimate, std.error, data = glm_out, 
                              studlab = paste(study_gxe), 
                              method.tau = 'SJ',
                              hakn = T,
                              sm = "OR",
                              fixed = TRUE, random = FALSE,
                              subgroup = study_design)

png(glue("~/Dropbox/Working/forest_plot_pure_ep_all_x_rs117868593_{hrc_version}_basic_covar.png"), 
    height = 13, 
    width = 10, 
    units = "in", res = 300)
meta::forest(results_meta, 
             leftcols = c('studlab', 'TE', 'seTE', 'effect', 'pval', 'w.random'),
             xlim = c(0.1, 10),
             layout = "JAMA",
             test.overall.fixed = TRUE,
             test.overall.random = FALSE)
grid.text("pure_ep_all x rs117868593", 0.5, 0.98, gp = gpar(cex = 2))
dev.off()








# ----- pure_ep_all x rs117868593 (region) -----

# forest plot
glm_out <- out %>% 
  tidyr::nest(data = -study_gxe) %>% 
  dplyr::mutate(fit = purrr::map(data, ~glm(model_formula, data = .x, family = "binomial")), 
                tidied = purrr::map(fit, ~tidy(.x)), quality = purrr::map(fit, ~glance(.x))) %>% 
  dplyr::select(-data, -fit) %>% 
  tidyr::unnest(tidied) %>% 
  tidyr::unnest(quality) %>% 
  dplyr::filter(grepl(':', term)) %>%
  dplyr::arrange(study_gxe) %>% 
  left_join(study_continent, 'study_gxe')

results_meta <- meta::metagen(estimate, std.error, data = glm_out, 
                              studlab = paste(study_gxe), 
                              method.tau = 'SJ',
                              hakn = T,
                              sm = "OR",
                              fixed = TRUE, random = FALSE,
                              subgroup = study_continent)

png(glue("~/Dropbox/Working/forest_plot_pure_ep_all_x_rs117868593_{hrc_version}_basic_covar_region.png"), 
    height = 7, 
    width = 10, 
    units = "in", res = 300)
meta::forest(results_meta, 
             leftcols = c('studlab', 'TE', 'seTE', 'effect', 'pval', 'w.random'),
             xlim = c(0.1, 10),
             layout = "JAMA",
             test.overall.fixed = TRUE,
             test.overall.random = FALSE)
grid.text("pure_ep_all x rs117868593 by region", 0.5, 0.98, gp = gpar(cex = 2))
dev.off()






# ----- pure_ep_all x rs10782186 (original) -----
# (do you have to recode anything? SNP, or exposure?)
# (according to study tables, not for this SNP)
dose <- qread("/media/work/gwis_test/hrt_ref_pm2/output/posthoc/dosage_chr6_117823508.qs")
out <- inner_join(input_data, dose, 'vcfid')

model_formula <- glue("outcome ~ pure_ep_all * chr6_117823508_T_C_dose + {glue_collapse(covariates_forest, sep = '+')}")

glm_out <- out %>% 
  tidyr::nest(data = -study_gxe) %>% 
  dplyr::mutate(fit = purrr::map(data, ~glm(model_formula, data = .x, family = "binomial")), 
                tidied = purrr::map(fit, ~tidy(.x)), quality = purrr::map(fit, ~glance(.x))) %>% 
  dplyr::select(-data, -fit) %>% 
  tidyr::unnest(tidied) %>% 
  tidyr::unnest(quality) %>% 
  dplyr::filter(grepl(":", term)) %>%
  dplyr::arrange(study_gxe) %>% 
  left_join(study_design, 'study_gxe')

results_meta <- meta::metagen(estimate, std.error, data = glm_out, 
                              studlab = paste(study_gxe), 
                              method.tau = 'SJ',
                              hakn = T,
                              sm = "OR",
                              fixed = TRUE, random = FALSE,
                              subgroup = study_design)

png(glue("~/Dropbox/Working/forest_plot_{exposure}_x_rs10782186_{hrc_version}_basic_covar.png"), 
    height = 8, 
    width = 10, 
    units = "in", res = 300)
meta::forest(results_meta, 
             leftcols = c('studlab', 'TE', 'seTE', 'effect', 'pval', 'w.random'),
             xlim = c(0.1, 10),
             layout = "JAMA",
             test.overall.fixed = TRUE,
             test.overall.random = FALSE)
grid.text(glue("{exposure} x rs10782186"), 0.5, 0.98, gp = gpar(cex = 2))
dev.off()




# ----- pure_ep_all x rs10782186 (by study COUNTRY) -----


# we know there's no heterogeneity by study design, so we can investigate country now
glm_out <- out %>% 
  tidyr::nest(data = -study_gxe) %>% 
  dplyr::mutate(fit = purrr::map(data, ~glm(model_formula, data = .x, family = "binomial")), 
                tidied = purrr::map(fit, ~tidy(.x)), quality = purrr::map(fit, ~glance(.x))) %>% 
  dplyr::select(-data, -fit) %>% 
  tidyr::unnest(tidied) %>% 
  tidyr::unnest(quality) %>% 
  dplyr::filter(grepl(":", term)) %>%
  dplyr::arrange(study_gxe) %>% 
  left_join(study_continent, 'study_gxe')

results_meta <- meta::metagen(estimate, std.error, data = glm_out, 
                              studlab = paste(study_gxe), 
                              method.tau = 'SJ',
                              hakn = T,
                              sm = "OR",
                              fixed = TRUE, random = FALSE,
                              subgroup = study_continent)

png(glue("~/Dropbox/Working/forest_plot_{exposure}_x_rs10782186_{hrc_version}_basic_covar_region.png"), 
    height = 8, 
    width = 10, 
    units = "in", res = 300)
meta::forest(results_meta, 
             leftcols = c('studlab', 'TE', 'seTE', 'effect', 'pval', 'w.random'),
             xlim = c(0.1, 10),
             layout = "JAMA",
             test.overall.fixed = TRUE,
             test.overall.random = FALSE)
grid.text(glue("{exposure} x rs10782186 by study continent"), 0.5, 0.98, gp = gpar(cex = 2))
dev.off()






















#-----------------------------------------------------------------------------#
# ------ ESTROGEN + PROGESTERONE (NEW DEF) ------ 
# Create forest plots of GxE term for main findings
#
# pure_ep_all x SNPs
# 
# rs117868593 - chr12_13670508_G_C
# rs10782186 -  chr6_117823508_T_C **
#-----------------------------------------------------------------------------#

exposure = 'pure_eo_all'

hrc_version = 'v2.3'
output_dir = paste0("/media/work/gwis/posthoc/", exposure, "/")
covariates <- sort(c('age_ref_imp', 'study_gxe', 'pc1', 'pc2', 'pc3'))

# input data
# use Yi's data subsets for pure_eo/ep variables
figi_gxe <- readRDS(glue("/media/work/gwis/data/FIGI_EpiData/FIGI_{hrc_version}_gxeset_analysis_data_glm.rds"))
load("/media/work/gwis/data/FIGI_EpiData/Data_HRT_USC_epi_v2.3-20201204.Rdata")

pure_eo_all_yi <- data.usc %>% 
  filter(!is.na(pure_eo_all))

# discrepancy is due to case-only studies being removed from GxE subset
input_data <- inner_join(figi_gxe, pure_eo_all_yi[, c('vcfid', 'pure_eo_all', 'pure_eo_allNo')], 'vcfid') %>% 
  mutate(pure_eo_allNo = as.numeric(pure_eo_allNo), 
         pure_eo_all = as.numeric(pure_eo_all))

# drop studies that have zero cell counts when tabulating outcome/pure_eo_allNo
drops <- data.frame(table(input_data$outcome, input_data[, exposure], 
                          input_data$study_gxe)) %>% 
  filter(Freq <= 1)

# sample sizes OK (N = 7637)
input_data <- filter(input_data, !study_gxe %in% unique(drops$Var3)) %>% 
  mutate(study_gxe = fct_drop(study_gxe))










# forest plot
covariates_forest <- sort(c('age_ref_imp', 'pc1', 'pc2', 'pc3'))
model_formula <- glue("outcome ~ pure_eo_all + {glue_collapse(covariates_forest, sep = '+')}")

glm_out <- input_data %>% 
  tidyr::nest(data = -study_gxe) %>% 
  dplyr::mutate(fit = purrr::map(data, ~glm(model_formula, data = .x, family = "binomial")), 
                tidied = purrr::map(fit, ~tidy(.x)), quality = purrr::map(fit, ~glance(.x))) %>% 
  dplyr::select(-data, -fit) %>% 
  tidyr::unnest(tidied) %>% 
  tidyr::unnest(quality) %>% 
  dplyr::filter(grepl(exposure, term)) %>%
  dplyr::arrange(study_gxe) %>% 
  left_join(study_continent, 'study_gxe') %>% 
  filter(study_gxe != "MEC_2")

results_meta <- meta::metagen(estimate, std.error, data = glm_out, 
                              studlab = paste(study_gxe), 
                              method.tau = 'SJ',
                              hakn = T,
                              sm = "OR",
                              fixed = TRUE, random = FALSE,
                              subgroup = study_continent)

png(glue("~/Dropbox/Working/forest_plot_{exposure}_{hrc_version}_basic_covar_region.png"), 
    height = 8, 
    width = 10, 
    units = "in", res = 300)
meta::forest(results_meta, 
             leftcols = c('studlab', 'TE', 'seTE', 'effect', 'pval', 'w.random'),
             xlim = c(0.1, 10),
             layout = "JAMA",
             test.overall.fixed = TRUE,
             test.overall.random = FALSE)
grid.text(glue("{exposure} by study region"), 0.5, 0.98, gp = gpar(cex = 2))
dev.off()





# ----- pure_eo_all x rs117868593 (original) -----

# super sanity check - recreate forest plot from publication
# it's right, use pure_ep_all
# move onto real results now

# Replicate table statistics for these two SNPs (MHT)
dose <- qread("/media/work/gwis_test/hrt_ref_pm2/output/posthoc/dosage_chr12_13670508.qs")
out <- inner_join(input_data, dose, 'vcfid')

# GxE
# covariates <- sort(c('age_ref_imp', 'study_gxe', 'pc1', 'pc2', 'pc3'))
# model_base <- glm(glue("outcome ~ pure_ep_all + chr12_13670508_G_C_dose + {glue_collapse(covariates, sep = '+')}"), data = out, family = 'binomial')
# model_int <-  glm(glue("outcome ~ pure_ep_all * chr12_13670508_G_C_dose + {glue_collapse(covariates, sep = '+')}"), data = out, family = 'binomial')
# lrtest(model_base, model_int) # OK



# forest plot
covariates_forest <- sort(c('age_ref_imp', 'pc1', 'pc2', 'pc3'))
model_formula <- glue("outcome ~ pure_eo_all * chr12_13670508_G_C_dose + {glue_collapse(covariates_forest, sep = '+')}")

glm_out <- out %>% 
  tidyr::nest(data = -study_gxe) %>% 
  dplyr::mutate(fit = purrr::map(data, ~glm(model_formula, data = .x, family = "binomial")), 
                tidied = purrr::map(fit, ~tidy(.x)), quality = purrr::map(fit, ~glance(.x))) %>% 
  dplyr::select(-data, -fit) %>% 
  tidyr::unnest(tidied) %>% 
  tidyr::unnest(quality) %>% 
  dplyr::filter(grepl(':', term)) %>%
  dplyr::arrange(study_gxe) %>% 
  left_join(study_design, 'study_gxe')

results_meta <- meta::metagen(estimate, std.error, data = glm_out, 
                              studlab = paste(study_gxe), 
                              method.tau = 'SJ',
                              hakn = T,
                              sm = "OR",
                              fixed = TRUE, random = FALSE,
                              subgroup = study_design)

png(glue("~/Dropbox/Working/forest_plot_pure_eo_all_x_rs117868593_{hrc_version}_basic_covar.png"), 
    height = 7, 
    width = 10, 
    units = "in", res = 300)
meta::forest(results_meta, 
             leftcols = c('studlab', 'TE', 'seTE', 'effect', 'pval', 'w.random'),
             xlim = c(0.1, 10),
             layout = "JAMA",
             test.overall.fixed = TRUE,
             test.overall.random = FALSE)
grid.text("pure_eo_all x rs117868593", 0.5, 0.98, gp = gpar(cex = 2))
dev.off()








# ----- pure_eo_all x rs117868593 (region) -----

# forest plot
glm_out <- out %>% 
  tidyr::nest(data = -study_gxe) %>% 
  dplyr::mutate(fit = purrr::map(data, ~glm(model_formula, data = .x, family = "binomial")), 
                tidied = purrr::map(fit, ~tidy(.x)), quality = purrr::map(fit, ~glance(.x))) %>% 
  dplyr::select(-data, -fit) %>% 
  tidyr::unnest(tidied) %>% 
  tidyr::unnest(quality) %>% 
  dplyr::filter(grepl(':', term)) %>%
  dplyr::arrange(study_gxe) %>% 
  left_join(study_continent, 'study_gxe')

results_meta <- meta::metagen(estimate, std.error, data = glm_out, 
                              studlab = paste(study_gxe), 
                              method.tau = 'SJ',
                              hakn = T,
                              sm = "OR",
                              fixed = TRUE, random = FALSE,
                              subgroup = study_continent)

png(glue("~/Dropbox/Working/forest_plot_pure_eo_all_x_rs117868593_{hrc_version}_basic_covar_region.png"), 
    height = 7, 
    width = 10, 
    units = "in", res = 300)
meta::forest(results_meta, 
             leftcols = c('studlab', 'TE', 'seTE', 'effect', 'pval', 'w.random'),
             xlim = c(0.1, 10),
             layout = "JAMA",
             test.overall.fixed = TRUE,
             test.overall.random = FALSE)
grid.text("pure_eo_all x rs117868593 by region", 0.5, 0.98, gp = gpar(cex = 2))
dev.off()






# ----- pure_eo_all x rs10782186 (original) -----
# (do you have to recode anything? SNP, or exposure?)
# (according to study tables, not for this SNP)
dose <- qread("/media/work/gwis_test/hrt_ref_pm2/output/posthoc/dosage_chr6_117823508.qs")
out <- inner_join(input_data, dose, 'vcfid')

model_formula <- glue("outcome ~ pure_eo_all * chr6_117823508_T_C_dose + {glue_collapse(covariates_forest, sep = '+')}")

glm_out <- out %>% 
  tidyr::nest(data = -study_gxe) %>% 
  dplyr::mutate(fit = purrr::map(data, ~glm(model_formula, data = .x, family = "binomial")), 
                tidied = purrr::map(fit, ~tidy(.x)), quality = purrr::map(fit, ~glance(.x))) %>% 
  dplyr::select(-data, -fit) %>% 
  tidyr::unnest(tidied) %>% 
  tidyr::unnest(quality) %>% 
  dplyr::filter(grepl(":", term)) %>%
  dplyr::arrange(study_gxe) %>% 
  left_join(study_design, 'study_gxe')

results_meta <- meta::metagen(estimate, std.error, data = glm_out, 
                              studlab = paste(study_gxe), 
                              method.tau = 'SJ',
                              hakn = T,
                              sm = "OR",
                              fixed = TRUE, random = FALSE,
                              subgroup = study_design)

png(glue("~/Dropbox/Working/forest_plot_{exposure}_x_rs10782186_{hrc_version}_basic_covar.png"), 
    height = 8, 
    width = 10, 
    units = "in", res = 300)
meta::forest(results_meta, 
             leftcols = c('studlab', 'TE', 'seTE', 'effect', 'pval', 'w.random'),
             xlim = c(0.1, 10),
             layout = "JAMA",
             test.overall.fixed = TRUE,
             test.overall.random = FALSE)
grid.text(glue("{exposure} x rs10782186"), 0.5, 0.98, gp = gpar(cex = 2))
dev.off()




# ----- pure_eo_all x rs10782186 (by study COUNTRY) -----


# we know there's no heterogeneity by study design, so we can investigate country now
glm_out <- out %>% 
  tidyr::nest(data = -study_gxe) %>% 
  dplyr::mutate(fit = purrr::map(data, ~glm(model_formula, data = .x, family = "binomial")), 
                tidied = purrr::map(fit, ~tidy(.x)), quality = purrr::map(fit, ~glance(.x))) %>% 
  dplyr::select(-data, -fit) %>% 
  tidyr::unnest(tidied) %>% 
  tidyr::unnest(quality) %>% 
  dplyr::filter(grepl(":", term)) %>%
  dplyr::arrange(study_gxe) %>% 
  left_join(study_continent, 'study_gxe')

results_meta <- meta::metagen(estimate, std.error, data = glm_out, 
                              studlab = paste(study_gxe), 
                              method.tau = 'SJ',
                              hakn = T,
                              sm = "OR",
                              fixed = TRUE, random = FALSE,
                              subgroup = study_continent)

png(glue("~/Dropbox/Working/forest_plot_{exposure}_x_rs10782186_{hrc_version}_basic_covar_region.png"), 
    height = 8, 
    width = 10, 
    units = "in", res = 300)
meta::forest(results_meta, 
             leftcols = c('studlab', 'TE', 'seTE', 'effect', 'pval', 'w.random'),
             xlim = c(0.1, 10),
             layout = "JAMA",
             test.overall.fixed = TRUE,
             test.overall.random = FALSE)
grid.text(glue("{exposure} x rs10782186 by study continent"), 0.5, 0.98, gp = gpar(cex = 2))
dev.off()




















