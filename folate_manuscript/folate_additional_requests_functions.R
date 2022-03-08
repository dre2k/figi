
# data_epi <- input_data
# exposure = 'folate_diet400qcm'
# covariates
# hrc_version = 'v3.0'
# strata = 'all'
# categorical = F



#' Title
#'
#' @param data_epi 
#' @param exposure 
#' @param covariates 
#' @param hrc_version 
#' @param strata 
#' @param categorical 
#'
#' @return
#' @export
#'
#' @examples
create_forest_plot_rmarkdown <- function(data_epi,
                                         exposure,
                                         covariates,
                                         hrc_version,
                                         strata = 'all',
                                         categorical = T) {
  
  covariates <- covariates[which(covariates != 'study_gxe')]
  
  # specify stratified analyses
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
  # do we? don't they just output warnings now ..
  if(categorical == T) {
    drops <- data.frame(table(data_epi$study_gxe, data_epi$outcome, data_epi[, exposure])) %>% 
      dplyr::filter(Freq == 0) %>% 
      dplyr::pull(Var1) %>% 
      unique(.)
    
    data_epi <- data_epi %>% 
      dplyr::filter(!study_gxe %in% drops)
  } else if(categorical == F) {
    drops <- data.frame(table(data_epi$study_gxe, data_epi$outcome)) %>% 
      dplyr::filter(Freq == 0) %>% 
      dplyr::pull(Var1) %>% 
      unique(.)
    
    data_epi <- data_epi %>% 
      dplyr::filter(!study_gxe %in% drops)
  }
  
  
  
  # create model term
  # model_formula <- deparse(reformulate(c(exposure, sort(covariates)), response = 'outcome'))
  model_formula <- glue("outcome ~ {exposure} + {glue_collapse(covariates, sep = '+')}")
  
  # create study_design data.frame
  study_design <- data_epi %>%
    dplyr::select(study_gxe, study_design) %>%
    filter(!duplicated(.)) %>%
    arrange(study_gxe)
  
  glm_out <- data_epi %>%
    tidyr::nest(data = -study_gxe) %>% 
    dplyr::mutate(
      fit = purrr::map(data, ~ glm(model_formula, data = .x, family = 'binomial')),
      tidied = purrr::map(fit, ~ tidy(.x)), 
      quality = purrr::map(fit, ~ glance(.x))) %>% 
    dplyr::select(-data, -fit) %>% 
    tidyr::unnest(tidied) %>% 
    tidyr::unnest(quality) %>% 
    dplyr::filter(grepl(exposure, term)) %>% 
    # null.deviance > quantile(null.deviance, 0.01)) %>% 
    dplyr::arrange(study_gxe)
  
  # include study sample sizes and study design information
  meta_input <- get_counts_outcome_by_group(data_epi, 'outcome', 'study_gxe') %>%
    mutate(study_gxe = as.character(study_gxe)) %>% 
    inner_join(study_design, 'study_gxe') %>% 
    full_join(glm_out, 'study_gxe')
  
  results_meta <- meta::metagen(estimate,
                                std.error,
                                data=meta_input,
                                studlab=paste(study_gxe),
                                fixed = TRUE,
                                # random = TRUE,
                                method.tau = "SJ",
                                hakn = TRUE,
                                prediction=TRUE,
                                sm="OR",
                                subgroup = study_design)
  
  
  # plot_title = glue(strata, "\n{model_formula}", .na = "All")
  plot_title = glue(strata, " (N=", sum(glm_out$nobs), ")", "\n{model_formula}")
  
  meta::forest(results_meta,
               layout = "JAMA",
               # text.predict = "95% CI",
               # col.predict = "black",
               leftcols = c("studlab", "Control", "Case", "N", "effect", "ci", "p.value"),
               # digits.addcols=0,
               study.results=T,
               prediction = F,
               col.random = 'red',
               xlim = c(0.25, 4))
  grid.text(plot_title, 0.5, .98, gp=gpar(cex=1))
  
}











# 
# 
# 
# data_epi <- input_data
# exposure = 'folate_sup_yn'
# snp = '3:12041456:T:A'
# snp = '6:23445253:T:A'
# covariates
# strata = 'alcohol_ref'
# method = 'chiSqGxE'
# flip_allele = F
# path = '/media/work/gwis_test/folate_sup_yn/output/'
# 
# wtf <- create_stratified_gxe_rmarkdown(data_epi = input_data, 
#                                 exposure = exposure, 
#                                 snp = snp, 
#                                 covariates = covariates, 
#                                 strata = "alcohol_ref", 
#                                 method = "chiSqGxE", 
#                                 flip_allele = F)





#' Title
#'
#' @param data_epi 
#' @param exposure 
#' @param snp 
#' @param covariates 
#' @param strata 
#' @param method 
#' @param flip_allele 
#'
#' @return
#' @export
#'
#' @examples
create_stratified_gxe_rmarkdown <- function(data_epi, 
                                            exposure, 
                                            snp, 
                                            covariates, 
                                            strata = c("sex", "study_design", "cancer_site_sum2", "alcohol_ref"), 
                                            method = c("chiSqGxE", "two-step", "chiSqCase", "chiSq2df", "chiSq3df"), 
                                            flip_allele = F) {
  
  snp_info <- unlist(strsplit(snp, split = ":"))
  snpname_clean <- function(x) {
    tmp <- gsub("\\:", "\\_", x)
    tmp <- glue("chr{tmp}_dose")
    return(tmp)
  }
  
  snpfix <- snpname_clean(snp)
  snpfix_short <- paste0("chr", gsub("\\:", "\\_", snp))
  data_dose <- qread(glue("/media/work/gwis_test/{exposure}/output/posthoc/dosage_chr{snp_info[1]}_{snp_info[2]}.qs"))
  data <- inner_join(data_epi, data_dose, "vcfid")
  
  if (flip_allele == T) {
    snp_old <- snpfix
    snp_tmp <- unlist(strsplit(snpfix, split = "_"))
    chr <- snp_tmp[1]
    bp <- snp_tmp[2]
    a1 <- snp_tmp[3]
    a2 <- snp_tmp[4]
    snp_new <- glue("{chr}_{bp}_{a2}_{a1}_dose_flipped")
    data[[snp_new]] <- abs(2 - data[, snp_old])
  }
  else {
    snp_new <- snpfix
  }
  
  strata <- match.arg(strata)
  method <- match.arg(method)
  
  covariates_nostrata <- paste0(covariates[!covariates %in% strata], collapse = " + ")
  
  data[, "strata_num"] <- as.numeric(factor(data[, strata])) - 1
  data[, "exposure_num"] = as.numeric(data[, exposure])
  
  out <- list()
  out_all <- fit_gxe(data, exposure, snp_new, covariates)
  out[["all"]] <- out_all
  number_of_levels <- nlevels(factor(data[, strata]))
  
  
  for (level in seq(number_of_levels) - 1) {
    if (strata == "cancer_site_sum2") {
      index_vector <- which(data[, "strata_num"] == level | 
                              data[, "outcome"] == 0)
    }
    else {
      index_vector <- which(data[, "strata_num"] == level)
    }
    out_level <- fit_gxe(data[index_vector, ], exposure, 
                         snp_new, covariates_nostrata)
    out[[paste0(strata, "_", as.character(level))]] <- out_level
  }
  
  
  
  list_of_glms <- lapply(out, function(x) x[[1]])
  list_of_samplesizes <- lapply(list_of_glms, function(x) paste0(c("Ca=", "Co="), rev(as.character(table(x$model$outcome))), collapse = ","))
  coefs <- lapply(list_of_glms, function(x) (exp(coef(x))))
  
  if (strata == "sex") {
    col_label = paste0(c("All", "Female", "Male"), " (", 
                       list_of_samplesizes, ")")
  }
  else if (strata == "study_design") {
    col_label = paste0(c("All", "Cohort", "Case-Control"), 
                       " (", list_of_samplesizes, ")")
  }
  else if (strata == "cancer_site_sum2") {
    col_label = paste0(c("All", "Proximal", "Distal", "Rectal"), 
                       " (", list_of_samplesizes, ")")
  } 
  else if (strata == "alcohol_ref") {
    col_label = paste0(c("All", "Yes", "No"), 
                       " (", list_of_samplesizes, ")")
  }
  
  
  
  list_of_chisq <- lapply(out, function(x) x[[2]])
  if (method %in% c("chiSqGxE", "two-step", "chiSqCase")) {
    gxe_pvalues <- do.call(c, lapply(list_of_chisq, function(x) formatC(pchisq(x[[1]], 
                                                                               df = 1, lower.tail = F), format = "e", digits = 4)))
    notes <- c("(PC and Study estimates omitted from table)", 
               paste0("GxE term LRtest p = ", gxe_pvalues))
  }
  else if (method == "chiSq2df") {
    gxe_pvalues <- do.call(c, lapply(list_of_chisq, function(x) formatC(pchisq(x[[2]], 
                                                                               df = 2, lower.tail = F), format = "e", digits = 4)))
    notes <- c("(PC and Study estimates omitted from table)", 
               paste0("2DF LRtest p = ", gxe_pvalues))
  }
  else if (method == "chiSq3df") {
    gxe_pvalues <- do.call(c, lapply(list_of_chisq, function(x) formatC(pchisq((x[[2]] + 
                                                                                  x[[3]]), df = 3, lower.tail = F), format = "e", digits = 4)))
    notes <- c("(PC and Study estimates omitted from table)", 
               paste0("3DF LRtest p = ", gxe_pvalues))
  }
  
  # output
  stargazer_helper(list_of_glms, title = paste0(gsub("\\_", "\\\\_", strata), " stratified ", gsub("\\_", "\\\\_", snp_new), " x ", gsub("\\_", "\\\\_", exposure)), column.labels = col_label, coef = coefs, notes = notes, single.row = T)

  # cat(paste(out_html, collapse = "\n"), "\n", file = glue("~/Dropbox/gxe_models_{exposure}_{hrc_version}_{snpfix}_{glue_collapse(sort(covariates), sep = '_')}_stratified_by_{strata}.html"), append = F)
  # return(glue("~/Dropbox/gxe_models_{exposure}_{hrc_version}_{snpfix}_{glue_collapse(sort(covariates), sep = '_')}_stratified_by_{strata}.html"))
}


