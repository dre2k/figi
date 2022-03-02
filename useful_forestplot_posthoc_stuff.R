tmp <- qread(glue("/media/work/gwis_test/{exposure}/output/posthoc/dosage_chr10_101476905.qs")) %>% 
  inner_join(input_data, 'vcfid')

data_epi = input_data
exposure
# snp = "chr10_101476905_G_A_dose"
covariates <- sort(c('age_ref_imp', 'sex', 'pc1', 'pc2', 'pc3'))
covariates <- sort(c('age_ref_imp', 'sex'))
# covariates = c(covariates, 'bmi', 'diab')
# covariates = c(covariates, 'bmi', 'redmeatqc2', 'fruitqc2', 'vegetableqc2')
# covariates = c(covariates, 'bmi', 'diab', 'educ', 'smk_ever', 'redmeatqc2', 'fruitqc2', 'vegetableqc2')
hrc_version
path
forest_height = 15
forest_width = 8.5
funnel_height = 8
funnel_width = 8.5
strata = 'all'
categorical = F
group_var = 'study_gxe'

strata = "all"



# just make sure study isn't included in the covariate list please
# group_var should be either study or study_gxe ... 
create_forest_plot_ver2 <- function(data_epi, exposure, covariates, hrc_version, path, forest_height = 17, forest_width = 8.5, funnel_height = 8, 
                                    funnel_width = 8.5, strata = 'all', categorical = T, group_var) {
  

  if (strata == 'female') {
    data_epi <- dplyr::filter(data_epi, sex == 0)
    covariates <- covariates[which(!covariates %in% c('sex'))]
  } else if (strata == 'male') {
    data_epi <- dplyr::filter(data_epi, sex == 1)
    covariates <- covariates[which(!covariates %in% c('sex'))]
  } else if (strata == 'proximal') {
    data_epi <- dplyr::filter(data_epi, cancer_site_sum2 == 'proximal' | data_epi$outcome == 0)
  } else if (strata == 'distal') {
    data_epi <- dplyr::filter(data_epi, cancer_site_sum2 == 'distal' | data_epi$outcome == 0)
  } else if (strata == 'rectal') {
    data_epi <- dplyr::filter(data_epi, cancer_site_sum2 == 'rectal' | data_epi$outcome == 0)
  }
  
  # some subsets generate empty cells, need to remove them
  if(categorical == T) {
    drops <- data.frame(table(data_epi[,group_var], data_epi$outcome, data_epi[, exposure])) %>%
      dplyr::filter(Freq == 0) %>%
      dplyr::pull(Var1) %>%
      unique(.)
    
    data_epi <- data_epi %>%
      dplyr::filter(!study_gxe %in% drops)
  } else if(categorical == F) {
    drops <- data.frame(table(data_epi[, group_var], data_epi$outcome)) %>%
      dplyr::filter(Freq < 5) %>%
      dplyr::pull(Var1) %>%
      unique(.)
    data_epi <- data_epi %>%
      dplyr::filter(! (!!sym(group_var) %in% drops))

  }
  
  
  
  # create model term
  # model_formula <- Reduce(paste, deparse(reformulate(c(exposure, sort(covariates)), response = 'outcome')))
  # model_formula <- deparse(reformulate(c(exposure, sort(covariates)), response = 'outcome'))
  model_formula <- glue("outcome ~ {exposure} + {glue_collapse(sort(covariates), sep = '+')}")
  
  
  
  
  
  # create study_design data.frame
  study_design <- data_epi %>%
    dplyr::select(!!sym(group_var), study_design) %>%
    filter(!duplicated(.)) %>%
    arrange(!!sym(group_var))
  
  
  
  
  
  
  glm_out <- data_epi %>%
    tidyr::nest(data = -(!!sym(group_var))) %>%
    dplyr::mutate(
      fit = purrr::map(data, ~ glm(model_formula, data = .x, family = 'binomial')),
      tidied = purrr::map(fit, ~ tidy(.x)),
      quality = purrr::map(fit, ~ glance(.x))
    ) %>%
    dplyr::select(-data, -fit) %>%
    tidyr::unnest(tidied) %>% tidyr::unnest(quality) %>%
    # dplyr::filter(grepl(":", term),
    #               null.deviance > quantile(null.deviance, 0.01)) %>%
    dplyr::filter(grepl(exposure, term)) %>%
    dplyr::arrange(!!sym(group_var))
  
  
  
  
  
  
  
  
  # include study sample sizes and study design information
  meta_input <- get_counts_outcome_by_group(data_epi, 'outcome', group_var) %>%
    mutate(study = as.character(!!sym(group_var))) %>%
    inner_join(study_design, group_var) %>%
    full_join(glm_out, group_var)
  
  results_meta <- meta::metagen(estimate,
                                std.error,
                                data=meta_input,
                                studlab=paste(study),
                                comb.fixed = FALSE,
                                comb.random = TRUE,
                                method.tau = "SJ",
                                hakn = TRUE,
                                prediction=TRUE,
                                sm="OR",
                                byvar = study_design)
  
  # results_meta
  # nrow(data_epi)
  # table(data_epi$outcome)
  

  # plot_title = glue(strata, "\n{model_formula}", .na = "All")
  plot_title = glue(strata, " (N=", sum(glm_out$nobs), ")", "\n{model_formula}")

  # png(glue("~/Dropbox/Working/forest_plot_{exposure}_{hrc_version}_", glue_collapse(covariates, "_"), "_{strata}.png"), height = forest_height, width = forest_width, units = 'in', res = 150)
  png(glue("{path}/forest_plot_{exposure}_{hrc_version}_", glue_collapse(covariates, "_"), "_{strata}.png"), height = forest_height, width = forest_width, units = 'in', res = 150)
  
  meta::forest(results_meta,
               layout = "JAMA",
               # text.predict = "95% CI",
               # col.predict = "black",
               leftcols = c("studlab", "Control", "Case", "N", "effect", "ci", "w.random"),
               digits.addcols=0,
               study.results=T,
               prediction = F,
               col.random = 'red',
               xlim = c(0.25, 4))
  grid.text(plot_title, 0.5, .98, gp=gpar(cex=1))
  dev.off()

  # png(glue("~/Dropbox/Working/funnel_plot_{exposure}_{hrc_version}_", glue_collapse(covariates, "_"), "_{strata}.png"), height = funnel_height, width = funnel_width, units = 'in', res = 150)
  png(glue("{path}/funnel_plot_{exposure}_{hrc_version}_", glue_collapse(covariates, "_"), "_{strata}.png"), height = funnel_height, width = funnel_width, units = 'in', res = 150)
  meta::funnel(results_meta, sm="OR", studlab = T, pos = 4, col.random = 'red')
  dev.off()

  return(glue("{path}/forest_plot_{exposure}_{hrc_version}_", glue_collapse(covariates, "_"), "_{strata}.png"))
}

covariates <- sort(c('age_ref_imp', 'sex', 'energytot_imp'))
create_forest_plot_ver2(data_epi = input_data, exposure = exposure, covariates = covariates, hrc_version = hrc_version, strata = 'all', path = path, categorical = F, group_var = 'study')
create_forest_plot_ver2(data_epi = input_data, exposure = exposure, covariates = covariates, hrc_version = hrc_version, strata = 'male', path = path, categorical = F, group_var = 'study_gxe')
create_forest_plot_ver2(data_epi = input_data, exposure = exposure, covariates = covariates, hrc_version = hrc_version, strata = 'female', path = path, categorical = F, group_var = 'study_gxe')
create_forest_plot_ver2(data_epi = input_data, exposure = exposure, covariates = covariates, hrc_version = hrc_version, strata = 'proximal', path = path, categorical = F, group_var = 'study_gxe')
create_forest_plot_ver2(data_epi = input_data, exposure = exposure, covariates = covariates, hrc_version = hrc_version, strata = 'distal', path = path, categorical = F, group_var = 'study_gxe')
create_forest_plot_ver2(data_epi = input_data, exposure = exposure, covariates = covariates, hrc_version = hrc_version, strata = 'rectal', path = path, categorical = F, group_var = 'study_gxe')













# just make sure study isn't included in the covariate list please
# group_var should be either study or study_gxe ... 
create_forest_plot_ver3 <- function(data_epi, exposure, covariates, hrc_version, path, strata = 'all', categorical = T, group_var) {
  
  # some subsets generate empty cells, need to remove them
  if(categorical == T) {
    drops <- data.frame(table(data_epi[,group_var], data_epi$outcome, data_epi[, exposure])) %>%
      dplyr::filter(Freq == 0) %>%
      dplyr::pull(Var1) %>%
      unique(.)
    
    data_epi <- data_epi %>%
      dplyr::filter(!study_gxe %in% drops)
  } else if(categorical == F) {
    drops <- data.frame(table(data_epi[, group_var], data_epi$outcome)) %>%
      dplyr::filter(Freq < 5) %>%
      dplyr::pull(Var1) %>%
      unique(.)
    data_epi <- data_epi %>%
      dplyr::filter(! (!!sym(group_var) %in% drops))
    
  }
  
  # create model term
  # model_formula <- Reduce(paste, deparse(reformulate(c(exposure, sort(covariates)), response = 'outcome')))
  # model_formula <- deparse(reformulate(c(exposure, sort(covariates)), response = 'outcome'))
  model_formula <- glue("outcome ~ {exposure} + {glue_collapse(sort(covariates), sep = '+')}")
  
  # create study_design data.frame
  
  glm_out <- data_epi %>%
    tidyr::nest(data = -(!!sym(group_var))) %>%
    dplyr::mutate(
      fit = purrr::map(data, ~ glm(model_formula, data = .x, family = 'binomial')),
      tidied = purrr::map(fit, ~ tidy(.x)),
      quality = purrr::map(fit, ~ glance(.x))
    ) %>%
    dplyr::select(-data, -fit) %>%
    tidyr::unnest(tidied) %>% tidyr::unnest(quality) %>%
    # dplyr::filter(grepl(":", term),
    #               null.deviance > quantile(null.deviance, 0.01)) %>%
    dplyr::filter(grepl(exposure, term)) %>%
    dplyr::arrange(!!sym(group_var))
  
  # include study sample sizes and study design information
  meta_input <- get_counts_outcome_by_group(data_epi, 'outcome', group_var) %>%
    mutate(study = as.character(!!sym(group_var))) %>%
    full_join(glm_out, group_var)
  
  results_meta <- meta::metagen(estimate,
                                std.error,
                                data=meta_input,
                                studlab=paste(study),
                                comb.fixed = FALSE,
                                comb.random = TRUE,
                                method.tau = "SJ",
                                hakn = TRUE,
                                prediction=TRUE,
                                sm="OR")
  
  results_meta
  
}


cc <- filter(study_design, study_design == "Case-Control") %>% pull(study_gxe)
create_forest_plot_ver3(data_epi = input_data[which(input_data$study_gxe %in% cc),], exposure = exposure, covariates = covariates, hrc_version = hrc_version, strata = 'all', path = path, categorical = F, group_var = 'study_gxe')

cc <- filter(study_design, study_design == "Cohort") %>% pull(study_gxe)
create_forest_plot_ver3(data_epi = input_data[which(input_data$study_gxe %in% cc),], exposure = exposure, covariates = covariates, hrc_version = hrc_version, strata = 'all', path = path, categorical = F, group_var = 'study_gxe')
