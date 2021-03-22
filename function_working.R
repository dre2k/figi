# =========================================================================== #
# additional one-time modificatio to functions, as needed 
# =========================================================================== #



# --------------------------------------------------------------------------- #
# meta analysis ---- 
# --------------------------------------------------------------------------- #



dataset <- input_data
exposure
model_formula = model_form

meta_analysis_wrapper_pwrap_tableonly <- function(dataset, exposure, model_formula, subset_var) {
  
  # analysis subset (matches GxEScanR input)
  exposure_subset <- readRDS(paste0("/media/work/gwis/results/input/FIGI_", hrc_version, "_gxeset_", exposure, "_basic_covars_glm.rds"))[, 'vcfid']
  
  # create study_design data.frame
  study_design <- dataset %>%
    dplyr::select(study_gxe, study_design) %>%
    filter(!duplicated(.)) %>%
    arrange(study_gxe)
  
  # create temporary analysis dataset subset
  tmp <- dataset %>% 
    filter(vcfid %in% exposure_subset) %>% 
    filter(case_when(subset_var == "female" ~ sex == 0, 
                     subset_var == "male" ~ sex == 1,
                     subset_var == "proximal" ~ cancer_site_sum2 == "proximal" | outcome == 0,
                     subset_var == "distal" ~ cancer_site_sum2 == "distal" | outcome == 0,
                     subset_var == "rectal" ~ cancer_site_sum2 == "rectal" | outcome == 0,
                     TRUE ~ TRUE))
  dataset_tmp <- remove_low_cell_count(tmp, exposure)
  
  # get study sample sizes
  count_gxe_subset <- get_counts_outcome_by_group(dataset_tmp, 'outcome', 'study_gxe') %>%
    mutate(study_gxe = as.character(study_gxe))
  count_gxe_all <- full_join(study_design, count_gxe_subset, 'study_gxe')
  
  # run GLM by study_gxe, merge with counts table for meta-analysis
  gxe_glm <- get_estimates_e_by_group_pwalk(dataset_tmp, exposure, model_formula) %>%
    mutate(study_gxe = as.character(study_gxe))
  gxe_meta <- dplyr::full_join(count_gxe_all, gxe_glm, 'study_gxe')
  
  # don't display warnings
  oldw <- getOption("warn")
  options(warn = -1)
  
  results_meta <- meta::metagen(estimate,
                                std.error,
                                data=gxe_meta,
                                studlab=paste(study_gxe),
                                comb.fixed = FALSE,
                                comb.random = TRUE,
                                method.tau = "SJ",
                                hakn = TRUE,
                                prediction=TRUE,
                                sm="OR",
                                byvar = study_design)
  return(results_meta)
}














meta_analysis_execute2 <- function(dataset, exposure, hrc_version, covariates, output_dir, filename_suffix2 = "", nosex=F) {
  
  # assemble data.frame of arguments for meta_analysis_wrapper function
  
  ## stratifying variables
  stratifying_vars_all <- c("all")
  stratifying_vars_sex <- c("female", "male")
  
  if(nosex) {
    stratifying_vars <- c(stratifying_vars_all)
  } else {
    stratifying_vars <- c(stratifying_vars_all, stratifying_vars_sex)
  }
  
  ## model formulas - check for and remove study_gxe and PCs variables
  covariates_sorted <- sort(covariates[which(!covariates %in% c("study_gxe", paste0(rep('pc', 20), seq(1, 20))))])
  covariates_sorted <- sort(covariates[which(!covariates %in% c("study_gxe"))])
  
  model_formula_all <- Reduce(paste, deparse(reformulate(c(exposure, covariates_sorted), response = 'outcome')))
  model_formula_sex <- Reduce(paste, deparse(reformulate(c(exposure, covariates_sorted[!covariates_sorted %in% c("sex")]), response = 'outcome')))
  
  if(nosex) {
    model_formula <- c(rep(model_formula_all, length(stratifying_vars_all)))
  } else {
    model_formula <- c(rep(model_formula_all, length(stratifying_vars_all)), rep(model_formula_sex, length(stratifying_vars_sex)))
  }
  
  ## figure titles
  forest_plot_title <- purrr::map2_chr(.x = stratifying_vars, .y = model_formula, ~ paste0(.x, "\n", .y))
  
  ## filename
  filename_all <- paste0(rep(paste0(covariates_sorted, collapse = "_"), length(stratifying_vars_all)), "_", stratifying_vars_all, filename_suffix2)
  filename_sex <- paste0(rep(paste0(covariates_sorted[!covariates_sorted %in% c("sex")], collapse = "_"), length(stratifying_vars_sex)), "_", stratifying_vars_sex, filename_suffix2)
  
  if(nosex) {
    filename <- c(filename_all)
  } else {
    filename <- c(filename_all, filename_sex)
  }
  
  ## argument data.frame
  meta_analysis_args <- data.frame(
    exposure = exposure,
    model_formula = model_formula,
    subset_var = stratifying_vars,
    output_dir = output_dir,
    forest_plot_title = forest_plot_title,
    filename_suffix = filename,
    forest_height = 17, 
    forest_width = 8.5, 
    funnel_height = 8, 
    funnel_width = 8.5, 
    stringsAsFactors = F)
  
  # execute
  purrr::pwalk(meta_analysis_args, meta_analysis_wrapper_pwrap, dataset = dataset)
}







# --------------------------------------------------------------------------- #
# meta analysis, study ---- 
# --------------------------------------------------------------------------- #







meta_analysis_execute2 <- function(dataset, exposure, hrc_version, covariates, output_dir, filename_suffix2 = "", nosex=F) {
  
  # assemble data.frame of arguments for meta_analysis_wrapper function
  
  ## stratifying variables
  stratifying_vars_all <- c("all")
  stratifying_vars_sex <- c("female", "male")
  
  if(nosex) {
    stratifying_vars <- c(stratifying_vars_all)
  } else {
    stratifying_vars <- c(stratifying_vars_all, stratifying_vars_sex)
  }
  
  ## model formulas - check for and remove study_gxe and PCs variables
  covariates_sorted <- sort(covariates[which(!covariates %in% c("study_gxe", paste0(rep('pc', 20), seq(1, 20))))])
  covariates_sorted <- sort(covariates[which(!covariates %in% c("study_gxe"))])
  
  model_formula_all <- Reduce(paste, deparse(reformulate(c(exposure, covariates_sorted), response = 'outcome')))
  model_formula_sex <- Reduce(paste, deparse(reformulate(c(exposure, covariates_sorted[!covariates_sorted %in% c("sex")]), response = 'outcome')))
  
  if(nosex) {
    model_formula <- c(rep(model_formula_all, length(stratifying_vars_all)))
  } else {
    model_formula <- c(rep(model_formula_all, length(stratifying_vars_all)), rep(model_formula_sex, length(stratifying_vars_sex)))
  }
  
  ## figure titles
  forest_plot_title <- purrr::map2_chr(.x = stratifying_vars, .y = model_formula, ~ paste0(.x, "\n", .y))
  
  ## filename
  filename_all <- paste0(rep(paste0(covariates_sorted, collapse = "_"), length(stratifying_vars_all)), "_", stratifying_vars_all, filename_suffix2)
  filename_sex <- paste0(rep(paste0(covariates_sorted[!covariates_sorted %in% c("sex")], collapse = "_"), length(stratifying_vars_sex)), "_", stratifying_vars_sex, filename_suffix2)
  
  if(nosex) {
    filename <- c(filename_all)
  } else {
    filename <- c(filename_all, filename_sex)
  }
  
  ## argument data.frame
  meta_analysis_args <- data.frame(
    exposure = exposure,
    model_formula = model_formula,
    subset_var = stratifying_vars,
    output_dir = output_dir,
    forest_plot_title = forest_plot_title,
    filename_suffix = filename,
    forest_height = 10, 
    forest_width = 8.5, 
    funnel_height = 8, 
    funnel_width = 8.5, 
    stringsAsFactors = F)
  
  # execute
  purrr::pwalk(meta_analysis_args, meta_analysis_wrapper_pwrap2, dataset = dataset)
}



meta_analysis_wrapper_pwrap2 <- function(dataset, exposure, model_formula, subset_var, output_dir, forest_plot_title, filename_suffix, forest_height, forest_width, funnel_height, funnel_width) {
  
  # analysis subset (matches GxEScanR input)
  exposure_subset <- readRDS(paste0("/media/work/gwis/results/input/FIGI_", hrc_version, "_gxeset_", exposure, "_basic_covars_glm.rds"))[, 'vcfid']
  
  # create study_design data.frame
  study_design <- dataset %>%
    dplyr::select(study, study_design) %>%
    filter(!duplicated(.)) %>%
    arrange(study)
  
  # create temporary analysis dataset subset
  tmp <- dataset %>% 
    filter(vcfid %in% exposure_subset) %>% 
    filter(case_when(subset_var == "female" ~ sex == 0, 
                     subset_var == "male" ~ sex == 1,
                     subset_var == "proximal" ~ cancer_site_sum2 == "proximal" | outcome == 0,
                     subset_var == "distal" ~ cancer_site_sum2 == "distal" | outcome == 0,
                     subset_var == "rectal" ~ cancer_site_sum2 == "rectal" | outcome == 0,
                     TRUE ~ TRUE))
  dataset_tmp <- remove_low_cell_count(tmp, exposure)
  
  # get study sample sizes
  count_gxe_subset <- get_counts_outcome_by_group(dataset_tmp, 'outcome', 'study') %>%
    mutate(study = as.character(study))
  count_gxe_all <- full_join(study_design, count_gxe_subset, 'study')
  
  # run GLM by study, merge with counts table for meta-analysis
  gxe_glm <- get_estimates_e_by_group_pwalk2(dataset_tmp, exposure, model_formula) %>%
    mutate(study = as.character(study))
  gxe_meta <- dplyr::full_join(count_gxe_all, gxe_glm, 'study')
  
  # don't display warnings
  oldw <- getOption("warn")
  options(warn = -1)
  
  results_meta <- meta::metagen(estimate,
                                std.error,
                                data=gxe_meta,
                                studlab=paste(study),
                                comb.fixed = FALSE,
                                comb.random = TRUE,
                                method.tau = "SJ",
                                hakn = TRUE,
                                prediction=TRUE,
                                sm="OR",
                                byvar = study_design)
  options(warn = oldw)
  
  png(paste0(output_dir, "meta_analysis_", exposure,  "_", filename_suffix, ".png"), height = forest_height, width = forest_width, units = 'in', res = 150)
  meta::forest(results_meta,
               layout = "JAMA",
               # text.predict = "95% CI",
               # col.predict = "black",
               leftcols = c("studlab", "Control", "Case", "N", "effect", "ci", "w.random"),
               digits.addcols=0,
               study.results=T,
               prediction = F,
               col.random = 'red')
  grid.text(forest_plot_title, 0.5, .98, gp=gpar(cex=1))
  dev.off()
  
  png(paste0(output_dir, "funnel_plot_", exposure,  "_", filename_suffix, ".png"), height = funnel_height, width = funnel_width, units = 'in', res = 150)
  oldw <- getOption("warn")
  options(warn = -1)
  meta::funnel(results_meta, sm="OR", studlab = T, pos = 4, col.random = 'red')
  options(warn = oldw)
  dev.off()
}





get_estimates_e_by_group_pwalk2 <- function(dataset, exposure, model_formula) {
  
  # group by study_gxe
  dataset <- dataset %>%
    group_by(study)
  
  # run glm with model_formula
  results_beta <- dplyr::do(dataset, broom::tidy(glm(model_formula, data = . , family = 'binomial'), conf.int = T))
  
  # clean up output, return exposure main effect
  results <- results_beta %>%
    dplyr::ungroup() %>%
    dplyr::filter(grepl(exposure, term)) %>% 
    dplyr::arrange(study)
  return(results)
}





## new script to get LD snps list with appropriate statistics for leads to followup on ##


# for example, take one of the hits and follow-up here, generalize
# use the locuszoom outputs for convenience















