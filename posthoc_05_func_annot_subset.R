#=============================================================================#
# FUNCTIONAL ANNOTATION SUBSET ANALYSIS
# Scores are SVM predictions based on Scacheri data
# (Anna Scherbina)
#
# create plots, including expectation based hybrid
#
# make sure to source this script in the results code 
# (to access global environments...)
#=============================================================================#
library(tidyverse)
library(data.table)
library(qqman)
library(EasyStrata)
library(lmtest)
library(figifs)

#------------------------------------------------------------------------------#
# ------- process SVM scores (filter and save RDS, only need run once) ------- #
# svm_pooled <- fread("/media/work/svm_scores/gecco.snps.flora.svm.delta.scores.txt.gz") %>%
#   filter(abs(deltaScore) > 3)
# saveRDS(svm_pooled, file = glue("/media/work/svm_scores/svm_pooled_filter_sd3.rds"))
# 
# svm_tumor <- fread("/media/work/svm_scores/gecco.snps.flora.svm.scores.scaled.txt.gz") %>%
#   dplyr::rename(SNP = V1) %>%
#   dplyr::select(SNP, dnase_v.ALT.REF) %>%
#   filter(abs(dnase_v.ALT.REF) > 3)
# saveRDS(svm_tumor, file = glue("/media/work/svm_scores/svm_tumor_filter_sd3.rds"))
# 
# svm_control <- fread("/media/work/svm_scores/gecco.snps.flora.svm.scores.scaled.txt.gz") %>%
#   dplyr::rename(SNP = V1) %>%
#   dplyr::select(SNP, dnase_c.ALT.REF) %>%
#   filter(abs(dnase_c.ALT.REF) > 3)
# saveRDS(svm_control, file = glue("/media/work/svm_scores/svm_control_filter_sd3.rds"))



#-----------------------------------------------------------------------------------#
# ------- output SNPs for binarydosage (for expectation based hybrid plots) ------- #
# (not doing this anymore, no expectation based hybrid for functional subset)

# input - gxescan output and step1 statistic
# output - writes RDS files to output_dir for binarydosage
# output_bin_snps <- function(x, stats_step1, output_suffix) {
#   out <- data_twostep_eh(x, stats_step1)
# 
#   # only plot bins 1-10 (other bins too large)
#   for(y in 1:10) {
#     tmp <- filter(out, bin_number == y) %>% pull(SNP)
#     saveRDS(tmp, file = glue("{output_dir}expectation_hybrid/{exposure}_snplist_twostep_func_annot_subset_{output_suffix}_{stats_step1}_bin{y}.rds"))
#   }
# }
# 
# 
# output_bin_snps_wrap <- function(x) {
# 
#   # svm_pooled
#   svm_pooled <- readRDS("/media/work/svm_scores/svm_pooled_filter_sd3.rds")
#   xx <- x %>%
#     filter(SNP %in% svm_pooled$SNP)
# 
#   output_bin_snps(x = xx, 'chiSqG', "pooled")
#   output_bin_snps(x = xx, 'chiSqGE', "pooled")
#   output_bin_snps(x = xx, 'chiSqEDGE', "pooled")
# 
#   # svm_tumor
#   svm_tumor <- readRDS("/media/work/svm_scores/svm_tumor_filter_sd3.rds")
#   xx <- x %>%
#     filter(SNP %in% svm_tumor$SNP)
# 
#   output_bin_snps(x = xx, 'chiSqG', "tumor")
#   output_bin_snps(x = xx, 'chiSqGE', "tumor")
#   output_bin_snps(x = xx, 'chiSqEDGE', "tumor")
# 
#   # svm_control
#   svm_control <- readRDS("/media/work/svm_scores/svm_control_filter_sd3.rds")
#   xx <- x %>%
#     filter(SNP %in% svm_control$SNP)
# 
#   output_bin_snps(x = xx, 'chiSqG', "control")
#   output_bin_snps(x = xx, 'chiSqGE', "control")
#   output_bin_snps(x = xx, 'chiSqEDGE', "control")
# }
# 
# # this is fine, move on
# output_bin_snps_wrap(gxe)





#--------------------------------------------------------------------#
# ------- generate qq, manhattan, and regular two-step plots ------- #
# this has to be done 3x (pooled, tumor, control...)


# D|G p < 5e-8
clump_g_file <- paste0("/media/work/gwis/clump_combined/FIGI_", hrc_version, "_gxeset_", exposure, "_chiSqG_ldclump.clumped")
if(file.exists(clump_g_file)) {
  clump_g <- fread(clump_g_file)
} else {
  clump_g <- data.frame()
}

# E|G p < 5e-8 (CONTROLS ONLY)
# (note most exposures don't have genomewide sig association with a SNP)
clump_control_file <- paste0("/media/work/gwis/clump_combined/FIGI_", hrc_version, "_gxeset_", exposure, "_chiSqControl_ldclump.clumped")
if(file.exists(clump_control_file)) {
  clump_control <- fread(clump_control_file)
} else {
  clump_control <- data.frame()
}

# using plink LD clump output, create vector of SNPs to remove (GWAS, D|G, or E|G + SNPs in LD)
format_clump <- function(x) {
  # x should be the plink output *.clumped read as dataframe
  tmp1 <- as.list(x$SP2)
  tmp2 <- sapply(tmp1, function(y) strsplit(y, split = ','))
  tmp3 <- do.call(c, tmp2) %>% 
    gsub("\\(1\\)", "", .) %>% 
    c(x$SNP) # add tag SNPs back into vector
}

# these are vectors of SNPs to exclude (chr:bp:ref:alt)
exclude_gwas <- fread("~/data/Annotations/gwas_141_ld_annotation_july2020.txt") %>% 
  mutate(SNP2 = paste0(Chr, ":", Pos)) %>% 
  pull(SNP2)
exclude_g    <- format_clump(clump_g)
exclude_eg_control <- format_clump(clump_control)












plot_funcannot <- function(x, exposure, covariates, output_dir, filename_suffix) {
  
  # # qqplots
  create_qqplot_ramwas(x = x, stat = 'chiSqG',    df = 1, filename_suffix = filename_suffix, output_dir = output_dir)
  create_qqplot_ramwas(x = x, stat = 'chiSqGxE',  df = 1, filename_suffix = filename_suffix, output_dir = output_dir)
  create_qqplot_ramwas(x = x, stat = 'chiSq2df',  df = 2, filename_suffix = filename_suffix, output_dir = output_dir)
  create_qqplot_ramwas(x = x, stat = 'chiSq3df',  df = 3, filename_suffix = filename_suffix, output_dir = output_dir)
  create_qqplot_ramwas(x = x, stat = 'chiSqCase', df = 1, filename_suffix = filename_suffix, output_dir = output_dir)

  
  # # manhattan plots
  # sig_line_value = formatC( (0.05 / nrow(x)), format = "e", digits = 2)
  sig_line_value = 0.05 / nrow(x)
  create_manhattanplot(x = x, exposure, stat = 'chiSqG',    output_dir = output_dir, annotation_file = annotation_file, filename_suffix = filename_suffix, sig_line = sig_line_value)
  create_manhattanplot(x = x, exposure, stat = 'chiSqGxE',  output_dir = output_dir, annotation_file = annotation_file, filename_suffix = filename_suffix, sig_line = sig_line_value)
  create_manhattanplot(x = x, exposure, stat = 'chiSq2df',  output_dir = output_dir, annotation_file = annotation_file, filename_suffix = filename_suffix, sig_line = sig_line_value)
  create_manhattanplot(x = x, exposure, stat = 'chiSq3df',  output_dir = output_dir, annotation_file = annotation_file, filename_suffix = filename_suffix, sig_line = sig_line_value)
  create_manhattanplot(x = x, exposure, stat = 'chiSqCase', output_dir = output_dir, annotation_file = annotation_file, filename_suffix = filename_suffix, sig_line = sig_line_value)
  
  # additional manhattan plots
  # for 2df/3df in particular, only generate the following:
  # - 2df - exclude gwas
  # - 3df - exclude gwas and E|G 
  # - everything else - exclude gwas only 
  
  # manhattan plots without GWAS loci
  x_nogwas <- x %>%
    filter(!SNP2 %in% exclude_gwas)

  # sig_line_value = formatC( (0.05 / nrow(x_nogwas)), format = "e", digits = 2)
  sig_line_value = 0.05 / nrow(x_nogwas)
  create_manhattanplot(x = x_nogwas, exposure, stat = 'chiSqG',   output_dir = output_dir, annotation_file = annotation_file,
                       filename_suffix = glue("{filename_suffix}_no_gwas"), sig_line = sig_line_value)
  create_manhattanplot(x = x_nogwas, exposure, stat = 'chiSqGxE', output_dir = output_dir, annotation_file = annotation_file,
                       filename_suffix = glue("{filename_suffix}_no_gwas"), sig_line = sig_line_value)
  create_manhattanplot(x = x_nogwas, exposure, stat = 'chiSqCase', output_dir = output_dir, annotation_file = annotation_file,
                       filename_suffix = glue("{filename_suffix}_no_gwas"), sig_line = sig_line_value)
  create_manhattanplot(x = x_nogwas, exposure, stat = 'chiSq2df',  output_dir = output_dir, annotation_file = annotation_file,
                       filename_suffix = glue("{filename_suffix}_no_gwas"), sig_line = sig_line_value)
  
  # manhattan plots for 3df (exclude gwas and ge)
  xx_nogwas_noge <- x %>%
    dplyr::filter(!SNP2 %in% exclude_gwas,
                  !SNP %in% exclude_eg_control)

  # sig_line_value = formatC( (0.05 / nrow(xx_nogwas_noge)), format = "e", digits = 2)
  sig_line_value = 0.05 / nrow(xx_nogwas_noge)
  create_manhattanplot(x = xx_nogwas_noge, exposure, stat = 'chiSq3df',  output_dir = output_dir, annotation_file = annotation_file,
                       filename_suffix = glue("{filename_suffix}_no_gwas_no_ge"), sig_line = sig_line_value)
  
  # # original two-step plots
  plot_twostep(x = x, exposure = exposure, covars = covariates, binsToPlot = 8, stats_step1 = 'chiSqG',    sizeBin0 = 5, alpha = 0.05, output_dir = output_dir, filename_suffix = filename_suffix)
  plot_twostep(x = x, exposure = exposure, covars = covariates, binsToPlot = 8, stats_step1 = 'chiSqGE',   sizeBin0 = 5, alpha = 0.05, output_dir = output_dir, filename_suffix = filename_suffix)
  plot_twostep(x = x, exposure = exposure, covars = covariates, binsToPlot = 8, stats_step1 = 'chiSqEDGE', sizeBin0 = 5, alpha = 0.05, output_dir = output_dir, filename_suffix = filename_suffix)

  plot_twostep(x = x_nogwas, exposure = exposure, covars = covariates, binsToPlot = 8, stats_step1 = 'chiSqG',    sizeBin0 = 5, alpha = 0.05, output_dir = output_dir,
               filename_suffix = glue("{filename_suffix}_no_gwas"))
  plot_twostep(x = x_nogwas, exposure = exposure, covars = covariates, binsToPlot = 8, stats_step1 = 'chiSqGE',   sizeBin0 = 5, alpha = 0.05, output_dir = output_dir,
               filename_suffix = glue("{filename_suffix}_no_gwas"))
  plot_twostep(x = x_nogwas, exposure = exposure, covars = covariates, binsToPlot = 8, stats_step1 = 'chiSqEDGE', sizeBin0 = 5, alpha = 0.05, output_dir = output_dir,
               filename_suffix = glue("{filename_suffix}_no_gwas"))
  
  
}

# call the function above for each score subset
plot_funcannot_wrap <- function(x, exposure, covariates, output_dir, filename_suffix) {
  
  # svm_pooled
  svm_pooled <- readRDS("/media/work/svm_scores/svm_pooled_filter_sd3.rds")
  x1 <- x %>%
    filter(SNP %in% svm_pooled$SNP)
  plot_funcannot(x = x1, exposure = exposure, covariates = covariates, output_dir = output_dir, filename_suffix = glue("{filename_suffix}_pooled"))

  # svm_tumor
  svm_tumor <- readRDS("/media/work/svm_scores/svm_tumor_filter_sd3.rds")
  x2 <- x %>%
    filter(SNP %in% svm_tumor$SNP)
  plot_funcannot(x = x2, exposure = exposure, covariates = covariates, output_dir = output_dir, filename_suffix = glue("{filename_suffix}_tumor"))

  # svm_control
  svm_control <- readRDS("/media/work/svm_scores/svm_control_filter_sd3.rds")
  x3 <- x %>%
    filter(SNP %in% svm_control$SNP)
  plot_funcannot(x = x3, exposure = exposure, covariates = covariates, output_dir = output_dir, filename_suffix = glue("{filename_suffix}_control"))
}



# some lazy stuff with filename suffix but it's ok
plot_funcannot_wrap(gxe, exposure = exposure, covariates = covariates, output_dir = output_dir, filename_suffix = "_functional_subset")





















# 
# 
# #------------------------------------------#
# # ------- expectation based hybrid ------- #
# # again.. do this for each of the SVM scores
# 
# plot_funcannot_eh <- function(x, exposure, covariates, output_dir, filename_suffix = "", input_helper) {
#   
#   # for(simplem_step1_statistic in c("chiSqG")) {
#   for(simplem_step1_statistic in c("chiSqG", "chiSqGE", "chiSqEDGE")) {
# 
#     files_input <- mixedsort(list.files(glue("/media/work/gwis/posthoc/{exposure}/expectation_hybrid/"), 
#                                         pattern = glue("annot_subset_{input_helper}_{simplem_step1_statistic}(?:.+)output.rds"), full.names = T))[1:8]
#     files_list <- map(files_input, ~ readRDS(.x))
#     files_list <- map(files_list, function(x) if (is.null(x)) data.frame(x) else x) 
#     number_of_snps <- map_int(files_list, ~ length(names(.x)[!names(.x) %in% "vcfid"])) # -1 to remove vcfid column from count
#     names(number_of_snps) <- seq(1,8)
#     number_of_tests <- map_int(files_list, function(x) if(ncol(x) < 2) ncol(x) else meff_r(dat = x, PCA_cutoff = 0.995, fixLength = 150) ) 
#     names(number_of_tests) <- seq(1,8)
#     
#     output <- plot_twostep_eh(x,
#                          exposure = exposure,
#                          covars = covariates, 
#                          binsToPlot = 8, 
#                          stats_step1 = simplem_step1_statistic, 
#                          sizeBin0 = 5, 
#                          alpha = 0.05, 
#                          output_dir = output_dir, 
#                          filename_suffix = filename_suffix, 
#                          number_of_snps = number_of_snps, 
#                          number_of_tests = number_of_tests)
#     
#     saveRDS(output, file = glue("{output_dir}twostep_wht_{simplem_step1_statistic}_{exposure}_{filename_suffix}_df.rds"))
#     }
# }
# 
# 
# plot_funcannot_eh_wrap <- function(x, exposure, covariates, output_dir, filename_suffix = "") {
#   
#   # svm_pooled
#   svm_pooled <- readRDS("/media/work/svm_scores/svm_pooled_filter_sd3.rds")
#   xx <- x %>% 
#     filter(SNP %in% svm_pooled$SNP)
#   plot_funcannot_eh(xx, exposure = exposure, covariates = covariates, output_dir = output_dir, filename_suffix = glue("{filename_suffix}_pooled"), input_helper = "pooled")
# 
#   # svm_tumor
#   svm_tumor <- readRDS("/media/work/svm_scores/svm_tumor_filter_sd3.rds")
#   xx <- x %>%
#     filter(SNP %in% svm_tumor$SNP)
#   plot_funcannot_eh(xx, exposure = exposure, covariates = covariates, output_dir = output_dir, filename_suffix = glue("{filename_suffix}_tumor"), input_helper = "tumor")
# 
#   # svm_control
#   svm_control <- readRDS("/media/work/svm_scores/svm_control_filter_sd3.rds")
#   xx <- x %>%
#     filter(SNP %in% svm_control$SNP)
#   plot_funcannot_eh(xx, exposure = exposure, covariates = covariates, output_dir = output_dir, filename_suffix = glue("{filename_suffix}_control"), input_helper = "control")
# }
# 
# plot_funcannot_eh_wrap(gxescan_output, exposure, covariates, output_dir, filename_suffix = "expectation_hybrid_func_annot_subset")
# 
#   
