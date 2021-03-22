# =========================================================================== #
# functions not worth putting in figifs package 
# because I change these around so much with each exposure
# =========================================================================== #



# =========================================================================== #
# posthoc function wrappers
# =========================================================================== #

# pvalue table, gxe models
posthoc_run_models <- function(x, method, quartile = F) {
  walk(x, ~ pval_summary(ds = figi, exposure, .x, covariates_list, method, output_dir = output_dir))
  walk(x, ~ fit_gxe_covars(ds = figi, exposure, .x, covariates_list, method, output_dir = output_dir))
  
  if(quartile == F) {
    walk(x, ~ fit_stratified_or(ds = figi, exposure = exposure, snp = .x, hrc_version = hrc_version, covariates = covariates_set1, mod = mod, dosage = F, output_dir = output_dir))
  }
}


# reri and interaction plot
posthoc_run_reri_plot <- function(y) {
  # RERI
  walk(y, ~ reri_wrapper(x = figi, exposure = exposure, snp = .x, covariates = covariates, output_dir = output_dir))
  # interaction plot
  walk(y, ~ iplot_wrapper(x = figi, exposure = exposure, snp = .x, covariates = covariates, output_dir = output_dir)) 
}


# locuszoom and functional annotation plots
posthoc_create_plots <- function(x, statistic_to_plot) {
  
  locuszoom_helper <- function(snp) {
    tmp <- gsub("chr", "", snp) # have to turn snps into chr:bp format for locuszoom 
    snp_chr <- strsplit(tmp, "_")[[1]][1]
    snp_bp <- as.numeric(strsplit(tmp, '_')[[1]][2])
    paste(snp_chr, snp_bp, sep = ":")
  }
  
  # --- going to cat empty annotation to the denote-markers-file to avoid errors --- #
  snps_plot <- map_chr(x, locuszoom_helper)
  
  tmp <- read.table("/media/work/gwis/locuszoom/locuszoom_gecco_gwas_annotation.txt", header = T)
  
  locuszoom_marker_helper <- function(x) {
    tmp_xx <- c(paste0("chr", x), "", "white")
    tmp_out <- rbind(tmp, tmp_xx)
    write.table(tmp_out, file = paste0("/media/work/gwis/locuszoom/locuszoom_annotation_chr", x, "_", statistic_to_plot, ".txt"), quote = F, row.names = F, sep = '\t')
  }
  
  walk(snps_plot, locuszoom_marker_helper)
  #walk(snps_plot, ~ system(paste("bash ~/Dropbox/FIGI/FIGI_code/results/posthoc/posthoc_02_locuszoom.sh", exposure, hrc_version, .x, statistic_to_plot)))
  walk(snps_plot, ~ system(paste("bash ~/git/figi/posthoc_02_locuszoom.sh", exposure, hrc_version, .x, statistic_to_plot)))
  walk(snps_plot, ~ system(paste("Rscript ~/git/figi/posthoc_03_functional_annotation.R", exposure, .x, statistic_to_plot)))
}

# locuszoom and functional annotation plots
posthoc_create_plots_locuszoom <- function(x, statistic_to_plot) {
  
  locuszoom_helper <- function(snp) {
    tmp <- gsub("chr", "", snp) # have to turn snps into chr:bp format for locuszoom 
    snp_chr <- strsplit(tmp, "_")[[1]][1]
    snp_bp <- as.numeric(strsplit(tmp, '_')[[1]][2])
    paste(snp_chr, snp_bp, sep = ":")
  }
  
  # --- going to cat empty annotation to the denote-markers-file to avoid errors --- #
  snps_plot <- map_chr(x, locuszoom_helper)
  
  tmp <- read.table("/media/work/gwis/locuszoom/locuszoom_gecco_gwas_annotation.txt", header = T)
  
  locuszoom_marker_helper <- function(x) {
    tmp_xx <- c(paste0("chr", x), "", "white")
    tmp_out <- rbind(tmp, tmp_xx)
    write.table(tmp_out, file = paste0("/media/work/gwis/locuszoom/locuszoom_annotation_chr", x, "_", statistic_to_plot, ".txt"), quote = F, row.names = F, sep = '\t')
  }
  
  walk(snps_plot, locuszoom_marker_helper)
  walk(snps_plot, ~ system(paste("bash ~/git/figi/posthoc_02_locuszoom.sh", exposure, hrc_version, .x, statistic_to_plot)))
  walk(snps_plot, ~ system(paste("Rscript ~/git/figi/posthoc_03_functional_annotation_convertpng.R", exposure, .x, statistic_to_plot)))
}

# # pvalue table, gxe models
# posthoc_run_models <- function(x, method) {
#   walk(x, ~ pval_summary(ds = figi, exposure, .x, covariates_list, method, output_dir = output_dir))
#   walk(x, ~ fit_gxe_covars(ds = figi, exposure, .x, covariates_list, method, output_dir = output_dir))
#   walk(x, ~ fit_stratified_or_q4(ds = figi, exposure = exposure, snp = .x, hrc_version = hrc_version, covariates = covariates_set1, mod = mod, dosage = F, output_dir = output_dir))
# }










# =========================================================================== #
# functional annotation subset analysis 
# 
# simplify results - generate only:
# - qq/manhattan plots
# - two-step edge (normal)
# - two-step expectation based hybrid 
# =========================================================================== #



plot_funcannot <- function(x, exposure, covariates, output_dir, filename_suffix) {
  
  # # qqplots
  create_qqplot_ramwas(x = x, exposure = exposure, stat = 'chiSqG',    df = 1, filename_suffix = filename_suffix, output_dir = output_dir)
  create_qqplot_ramwas(x = x, exposure = exposure, stat = 'chiSqGxE',  df = 1, filename_suffix = filename_suffix, output_dir = output_dir)
  create_qqplot_ramwas(x = x, exposure = exposure, stat = 'chiSq2df',  df = 2, filename_suffix = filename_suffix, output_dir = output_dir)
  create_qqplot_ramwas(x = x, exposure = exposure, stat = 'chiSq3df',  df = 3, filename_suffix = filename_suffix, output_dir = output_dir)
  create_qqplot_ramwas(x = x, exposure = exposure, stat = 'chiSqControl',  df = 1, filename_suffix = filename_suffix, output_dir = output_dir)
  create_qqplot_ramwas(x = x, exposure = exposure, stat = 'chiSqCase', df = 1, filename_suffix = filename_suffix, output_dir = output_dir)
  
  # # manhattan plots
  # sig_line_value = formatC( (0.05 / nrow(x)), format = "e", digits = 2)
  sig_line_value = 0.05 / nrow(x)
  create_manhattanplot_ramwas(x = x, exposure, stat = 'chiSqG',    output_dir = output_dir, annotation_file = annotation_file, filename_suffix = filename_suffix, sig_line = sig_line_value)
  create_manhattanplot_ramwas(x = x, exposure, stat = 'chiSqGxE',  output_dir = output_dir, annotation_file = annotation_file, filename_suffix = filename_suffix, sig_line = sig_line_value)
  create_manhattanplot_ramwas(x = x, exposure, stat = 'chiSq2df',  output_dir = output_dir, annotation_file = annotation_file, filename_suffix = filename_suffix, sig_line = sig_line_value)
  create_manhattanplot_ramwas(x = x, exposure, stat = 'chiSq3df',  output_dir = output_dir, annotation_file = annotation_file, filename_suffix = filename_suffix, sig_line = sig_line_value)
  create_manhattanplot_ramwas(x = x, exposure, stat = 'chiSqControl', output_dir = output_dir, annotation_file = annotation_file, filename_suffix = filename_suffix, sig_line = sig_line_value)
  create_manhattanplot_ramwas(x = x, exposure, stat = 'chiSqCase', output_dir = output_dir, annotation_file = annotation_file, filename_suffix = filename_suffix, sig_line = sig_line_value)
  
  # manhattan plots without GWAS loci
  x_nogwas <- x %>%
    filter(!SNP2 %in% exclude_gwas)
  
  # sig_line_value = formatC( (0.05 / nrow(x_nogwas)), format = "e", digits = 2)
  sig_line_value = 0.05 / nrow(x_nogwas)
  create_manhattanplot_ramwas(x = x_nogwas, exposure, stat = 'chiSqG',   output_dir = output_dir, annotation_file = annotation_file,
                       filename_suffix = glue("{filename_suffix}_no_gwas"), sig_line = sig_line_value)
  create_manhattanplot_ramwas(x = x_nogwas, exposure, stat = 'chiSqGxE', output_dir = output_dir, annotation_file = annotation_file,
                       filename_suffix = glue("{filename_suffix}_no_gwas"), sig_line = sig_line_value)
  create_manhattanplot_ramwas(x = x_nogwas, exposure, stat = 'chiSqCase', output_dir = output_dir, annotation_file = annotation_file,
                       filename_suffix = glue("{filename_suffix}_no_gwas"), sig_line = sig_line_value)
  create_manhattanplot_ramwas(x = x_nogwas, exposure, stat = 'chiSq2df',  output_dir = output_dir, annotation_file = annotation_file,
                       filename_suffix = glue("{filename_suffix}_no_gwas"), sig_line = sig_line_value)
  create_manhattanplot_ramwas(x = x_nogwas, exposure, stat = 'chiSq3df',  output_dir = output_dir, annotation_file = annotation_file,
                       filename_suffix = glue("{filename_suffix}_no_gwas"), sig_line = sig_line_value)
  
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
  
  # # svm_tumor
  # svm_tumor <- readRDS("/media/work/svm_scores/svm_tumor_filter_sd3.rds")
  # x2 <- x %>%
  #   filter(SNP %in% svm_tumor$SNP)
  # plot_funcannot(x = x2, exposure = exposure, covariates = covariates, output_dir = output_dir, filename_suffix = glue("{filename_suffix}_tumor"))
  # 
  # # svm_control
  # svm_control <- readRDS("/media/work/svm_scores/svm_control_filter_sd3.rds")
  # x3 <- x %>%
  #   filter(SNP %in% svm_control$SNP)
  # plot_funcannot(x = x3, exposure = exposure, covariates = covariates, output_dir = output_dir, filename_suffix = glue("{filename_suffix}_control"))
}



# some lazy stuff with filename suffix but it's ok
# plot_funcannot_wrap(gxe, exposure = exposure, covariates = covariates, output_dir = output_dir, filename_suffix = "_functional_subset")



#-----------------
##------------------------
###-----------------------------

# svm_pooled <- readRDS("/media/work/svm_scores/svm_pooled_filter_sd3.rds")
# x1 <- gxe %>%
#   filter(SNP %in% svm_pooled$SNP)
# 
# twostep_eh_snps(x1, 'chiSqG', output_dir = glue('/media/work/gwis/twostep_expectation_hybrid/{exposure}/svm_subset/'))
# twostep_eh_snps(x1, 'chiSqGE', output_dir = glue('/media/work/gwis/twostep_expectation_hybrid/{exposure}/svm_subset/'))
# twostep_eh_snps(x1, 'chiSqEDGE', output_dir = glue('/media/work/gwis/twostep_expectation_hybrid/{exposure}/svm_subset/'))
# 
# twostep_eh_snps(gxe_twostep_gwas_step1, 'chiSqG', step1_source = "gecco")
# 










simplem_wrap2 <- function (x, exposure, covariates, simplem_step1_statistic, output_dir, 
          filename_suffix = "", include_gwas = T) 
{
  files_input <- mixedsort(list.files(glue("/media/work/gwis/twostep_expectation_hybrid/{exposure}/svm_subset"), 
                                      pattern = paste0(paste0("twostep_", simplem_step1_statistic, 
                                                              "_bin"), "(?:.+)", "output.rds"), full.names = T))
  files_list <- map(files_input, ~readRDS(.x))
  exclude_gwas_snps <- fread("~/data/Annotations/gwas_141_ld_annotation_july2020.txt") %>% 
    mutate(snps = paste0("X", Chr, ".", Pos)) %>% pull(snps)
  tmp_function <- function(zz) {
    zznames <- substr(names(zz), 1, nchar(names(zz)) - 4)
    zz_index <- !zznames %in% exclude_gwas_snps
    zz_out <- zz[, zz_index]
    return(zz_out)
  }
  if (include_gwas == F) {
    files_list <- map(files_list, ~tmp_function(.x))
  }
  number_of_snps <- map_int(files_list, ~ncol(.x)) - 1
  number_of_tests <- map_int(files_list, ~meff_r(dat = .x, 
                                                 PCA_cutoff = 0.995, fixLength = 150))
  plot_twostep_eh(x, exposure = exposure, covars = covariates, 
                  binsToPlot = 8, stats_step1 = simplem_step1_statistic, 
                  sizeBin0 = 5, alpha = 0.05, output_dir = output_dir, 
                  filename_suffix = glue("_expectation_hybrid{filename_suffix}"), 
                  number_of_snps = number_of_snps, number_of_tests = number_of_tests)
}


