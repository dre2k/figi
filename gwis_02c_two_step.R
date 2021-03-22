# --------------------------------------------------------------------------- #
# generate two-step method plots
# multiple iterations
# - original
# - gwis statistic (from fred hutch)
#
# wrapper function also outputs data.frame of significant/suggestive results
# it's better to do it in this step because the function to process data also generates bins/pvalues
#
# expectation based results will be in separate script (making modifications to it)
# --------------------------------------------------------------------------- #


#=============================================================================#
#=============================================================================#
# USING EXPOSURE SUBSET D|G AS STEP 1 STATISTIC
#=============================================================================#
#=============================================================================#


#-----------------------------------------------------------------------------#
# Two-step methods
# step1:
# - using main effects statistics from exposure subset
# - main effects models are adjusted by exposure
message("two-step methods")
#-----------------------------------------------------------------------------#
plot_twostep(x = gxe, exposure = exposure, covars = covariates, binsToPlot = 10, stats_step1 = 'chiSqG', sizeBin0 = 5, alpha = 0.05, output_dir = output_dir, filename_suffix = "")
plot_twostep(x = gxe, exposure = exposure, covars = covariates, binsToPlot = 10, stats_step1 = 'chiSqGE', sizeBin0 = 5, alpha = 0.05, output_dir = output_dir, filename_suffix = "")
plot_twostep(x = gxe, exposure = exposure, covars = covariates, binsToPlot = 10, stats_step1 = 'chiSqEDGE', sizeBin0 = 5, alpha = 0.05, output_dir = output_dir, filename_suffix = "")

#-----------------------------------------------------------------------------#
# Two-step methods, no gwas
message("two-step methods, no gwas")
#-----------------------------------------------------------------------------#
gxe_nogwas <- gxe %>%
  filter(!SNP2 %in% exclude_gwas)

plot_twostep(x = gxe_nogwas, exposure = exposure, covars = covariates, binsToPlot = 10, stats_step1 = 'chiSqG', sizeBin0 = 5, alpha = 0.05, output_dir = output_dir, filename_suffix = "_no_gwas")
plot_twostep(x = gxe_nogwas, exposure = exposure, covars = covariates, binsToPlot = 10, stats_step1 = 'chiSqGE', sizeBin0 = 5, alpha = 0.05, output_dir = output_dir, filename_suffix = "_no_gwas")
plot_twostep(x = gxe_nogwas, exposure = exposure, covars = covariates, binsToPlot = 10, stats_step1 = 'chiSqEDGE', sizeBin0 = 5, alpha = 0.05, output_dir = output_dir, filename_suffix = "_no_gwas")

rm(gxe_nogwas)
gc()








#=============================================================================#
#=============================================================================#
# USING FRED HUTCH GWAS AS STEP 1 STATISTIC
#=============================================================================#
#=============================================================================#

#-----------------------------------------------------------------------------#
# Two-step methods, gwas step 1
# step1:
# - using FIGI GWAS set marginal statistics (N ~ 120K) - I ran this one
# - OR.... using meta-analysis p value (new update 04/17/2020)
message("two step methods, gwas step 1")
#-----------------------------------------------------------------------------#
# incorporate GWAS set main effects statistics to the two-step methods
# gwas_results <- readRDS("~/data/results/gwas/processed/FIGI_v2.3_gwasset_basic_covars_gxescan_results_short.rds")
gwas_results <- readRDS("/media/work/gwis/results/gecco/MarginalMeta_gecco_HRC_EUR_only.rds")

gxe_twostep_gwas_step1 <- gxe %>%
  dplyr::select(-betaG, -chiSqG, -chiSqEDGE, -chiSq3df) %>%
  inner_join(gwas_results, 'SNP') %>%
  mutate(chiSqEDGE = chiSqG + chiSqGE,
         chiSq3df = chiSqG + chiSqGxE + chiSqGE)

plot_twostep(x = gxe_twostep_gwas_step1, exposure = exposure, covars = covariates, binsToPlot = 10, stats_step1 = 'chiSqG', sizeBin0 = 5, alpha = 0.05, output_dir = output_dir, filename_suffix = "_gwas_step1")
plot_twostep(x = gxe_twostep_gwas_step1, exposure = exposure, covars = covariates, binsToPlot = 10, stats_step1 = 'chiSqGE', sizeBin0 = 5, alpha = 0.05, output_dir = output_dir, filename_suffix = "_gwas_step1")
plot_twostep(x = gxe_twostep_gwas_step1, exposure = exposure, covars = covariates, binsToPlot = 10, stats_step1 = 'chiSqEDGE', sizeBin0 = 5, alpha = 0.05, output_dir = output_dir, filename_suffix = "_gwas_step1")

#-----------------------------------------------------------------------------#
# Two-step methods, gwas step 1, no gwas loci
message("two step methods, gwas step 1, no gwas loci")
#-----------------------------------------------------------------------------#
gxe_twostep_gwas_step1_nogwas <- gxe_twostep_gwas_step1 %>%
  filter(!SNP2 %in% exclude_gwas)

plot_twostep(x = gxe_twostep_gwas_step1_nogwas, exposure = exposure, covars = covariates, binsToPlot = 10, stats_step1 = 'chiSqG', sizeBin0 = 5, alpha = 0.05, output_dir = output_dir, filename_suffix = "_gwas_step1_no_gwas")
plot_twostep(x = gxe_twostep_gwas_step1_nogwas, exposure = exposure, covars = covariates, binsToPlot = 10, stats_step1 = 'chiSqGE', sizeBin0 = 5, alpha = 0.05, output_dir = output_dir, filename_suffix = "_gwas_step1_no_gwas")
plot_twostep(x = gxe_twostep_gwas_step1_nogwas, exposure = exposure, covars = covariates, binsToPlot = 10, stats_step1 = 'chiSqEDGE', sizeBin0 = 5, alpha = 0.05, output_dir = output_dir, filename_suffix = "_gwas_step1_no_gwas")



