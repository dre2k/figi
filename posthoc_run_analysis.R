#=============================================================================#
# FIGI GWIS Post hoc analyses
# 03/16/2020
#
#=============================================================================#
library(tidyverse)
library(data.table)
library(ggplot2)
library(qqman)
library(table1)
library(meta)
library(rlang)
library(broom)
library(effects)
library(figifs)
library(DT)
library(grid)
library(sjPlot)
library(stargazer)
library(forcats)
library(kableExtra)
rm(list = ls())

# exposure = 'hrt_ref_pm2'
# covariates = c('age_ref_imp', 'sex', 'study_gxe', 'pc1', 'pc2', 'pc3')



args <- commandArgs(trailingOnly=T)
exposure <- args[1] # ex: asp_ref
covariates <- c(args[2:length(args)]) # space delim list of covars, turns into vector. Useful for plot titles

dir.create(paste0("/media/work/tmp/posthoc/", exposure), recursive = T, showWarnings = F)

# dat <- hrt
# exposure <- 'hrt_ref_pm2'
# snp <- 'X18.10396692.A.G'
# covariates = c("age_ref_imp", "study_gxe", "pc1", "pc2", "pc3")
# covariates = c("age_ref_imp")


# GLM wrap function. Outcome is ALWAYS ALWAYS ALWAYS 'outcome'
glm_wrapper <- function(dat, exposure, snp, covariates = c()) {
  
  # check if exposure variable is binary, fit models accordingly (categorical vs numeric)
  factor_or_numeric <- function(exposure) {
    if(length(table(exposure)) <= 4) {
      return(as.factor)
    } else {
      return(as.numeric)
    }}
  
  snp_rename <- (gsub('\\.', '_', snp) %>%  gsub('X', 'chr', .))

  tmp <- dat %>% 
    mutate(!!snp_rename := as.numeric(!!sym(snp))) %>% 
    mutate_at(vars(contains(exposure)), factor_or_numeric(exposure))

  # fit model, output plots and odds ratios data.frame
  mod_formula <- as.formula(paste0('outcome ~ ', exposure , "*", snp_rename, "+", paste0(covariates, collapse = "+")))
  mod <- glm(mod_formula, data = tmp, family = 'binomial')

  mod_list <- list(x=c(0,1,2))
  names(mod_list) <- snp_rename
  model_eff <- effect(paste0(exposure,":", snp_rename),
                      mod,
                      se = TRUE,
                      confidence.level = 0.95,
                      typical = mean,
                      family = 'binomial',
                      xlevels = mod_list)
  
  # output plots
  oldw <- getOption("warn")
  options(warn = -1)
  png(paste0("/media/work/tmp/posthoc/", exposure, "/effects_plot_", exposure, "_", snp_rename, "_", paste0(sort(covariates), collapse = "_"), "_e.png"), height = 4, width = 6, units = 'in', res = 150)
  print({
  plot(model_eff, snp_rename, multiline = T, confint = list(style = 'bands'), type = 'link', ylab = 'predicted outcome log-odds')
  })
  dev.off()

  png(paste0("/media/work/tmp/posthoc/", exposure, "/effects_plot_", exposure, "_", snp_rename, "_", paste0(sort(covariates), collapse = "_"), "_g.png"), height = 4, width = 6, units = 'in', res = 150)
  print({
  plot(model_eff, exposure, multiline = T, confint = list(style = 'auto'), type = 'link', ylab = 'predicted outcome log-odds')
  })
  dev.off()

  options(warn = oldw)
  
  # output data.frame
  # out <- as.data.frame(model_eff, type = 'link') %>%
  #   dplyr::select(-se) %>%
  #   dplyr::mutate(fit = fit / (1 - fit),
  #                 lower = lower / (1 - lower),
  #                 upper = upper / (1 - upper))
  
  ## you shouldn't need to calculate odds, since using 'link' gives you logodds 
  ## output odds ratio table .. 
  out <- as.data.frame(model_eff, type = 'link') %>%
    dplyr::select(-se) 
  
  e1ve0_g0 <- c("E=1 vs Ref", "G=0", (exp(out[2,3:5]) / (exp(out[1, 3:5]))))
  e1ve0_g1 <- c("E=1 vs Ref", "G=1", (exp(out[4,3:5]) / (exp(out[3, 3:5]))))
  e1ve0_g2 <- c("E=1 vs Ref", "G=2", (exp(out[6,3:5]) / (exp(out[5, 3:5]))))
  e0_g1vg0 <- c("E=0", "G=1 vs Ref", (exp(out[3,3:5]) / (exp(out[1, 3:5]))))
  e0_g2vg0 <- c("E=0", "G=2 vs Ref", (exp(out[5,3:5]) / (exp(out[1, 3:5]))))
  e1_g1vg0 <- c("E=1", "G=1 vs Ref", (exp(out[4,3:5]) / (exp(out[2, 3:5]))))
  e1_g2vg0 <- c("E=1", "G=2 vs Ref", (exp(out[6,3:5]) / (exp(out[2, 3:5]))))
  
  out2 <- data.frame(rbind(e1ve0_g0,
                           e1ve0_g1,
                           e1ve0_g2,
                           e0_g1vg0,
                           e0_g2vg0,
                           e1_g1vg0,
                           e1_g2vg0))
  names(out2) <- c("E", "G", "OR", "Lower", "Upper")
  
  saveRDS(out2, file = paste0("/media/work/tmp/posthoc/", exposure, "/effects_df_", exposure, "_", snp_rename, "_", paste0(sort(covariates), collapse = "_"), ".rds"), version = 2)
  }


# # maybe you should calculate ORs for them.. 
# # this will differ whether the variable is binary, Q4, or continuous. Get the results for binary asap, adapt for other exposures. 
# e1ve0_g0 <- c("E=1 vs Ref", "G=0", (exp(out[2,3:5]) / (exp(out[1, 3:5]))))
# e1ve0_g1 <- c("E=1 vs Ref", "G=1", (exp(out[4,3:5]) / (exp(out[3, 3:5]))))
# e1ve0_g2 <- c("E=1 vs Ref", "G=2", (exp(out[6,3:5]) / (exp(out[5, 3:5]))))
# e0_g1vg0 <- c("E=0", "G=1 vs Ref", (exp(out[3,3:5]) / (exp(out[1, 3:5]))))
# e0_g2vg0 <- c("E=0", "G=2 vs Ref", (exp(out[5,3:5]) / (exp(out[1, 3:5]))))
# e1_g1vg0 <- c("E=1", "G=1 vs Ref", (exp(out[4,3:5]) / (exp(out[2, 3:5]))))
# e1_g2vg0 <- c("E=1", "G=2 vs Ref", (exp(out[6,3:5]) / (exp(out[2, 3:5]))))
# 
# out2 <- data.frame(rbind(e1ve0_g0,
#       e1ve0_g1,
#       e1ve0_g2,
#       e0_g1vg0,
#       e0_g2vg0,
#       e1_g1vg0,
#       e1_g2vg0))
# 
# names(out2) <- c("E", "G", "OR", "Lower", "Upper")
# out2
# 
# 
# 
# # one more test - get odds ratios of predicted values 
# newdata = data.frame(age_ref_imp = mean(dat$age_ref_imp), hrt_ref_pm2 = "0", chr18_10396692_A_G = 1)
# predict(mod, newdata, type = 'link')
# 
# 
# # I just want to make sure that I can get odds ratio from logits.. 
# # yes I can
# tmp <- dat %>% 
#   mutate(!!snp_rename := as.numeric(!!sym(snp))) %>% 
#   mutate_at(vars(contains(exposure)), factor_or_numeric(exposure))
# 
# mod_formula <- as.formula(paste0('outcome ~ ', exposure , "+", snp_rename))
# mod <- glm(mod_formula, data = tmp, family = 'binomial')
# mod
# exp(coef(mod))
# 
# #0.9463954  for genotype
# 
# newdata = data.frame(hrt_ref_pm2 = "0", chr18_10396692_A_G = 0)
# predict(mod, newdata, type = 'link')
# 
# -0.338497
# 
# newdata = data.frame(hrt_ref_pm2 = "0", chr18_10396692_A_G = 1)
# predict(mod, newdata, type = 'link')
# 
# -0.3935919
# 
# 
# exp(-0.3935919) / exp(-0.338497)
#-----------------------------------------------------------------------------#
# apply glm_wrapper over SNPs
#
# be aware that for now the wrapper outputs logits
# it'll up to the leads to calculate odds ratios and 95% CI if they want
#-----------------------------------------------------------------------------#

# input data
geno <- readRDS(paste0("/media/work/tmp/posthoc/gwis_sig_results_output_", exposure, ".rds"))
dat <- readRDS(paste0("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_", exposure, "_basic_covars_glm.rds")) %>% 
  inner_join(geno, 'vcfid') %>% 
  rename_at(vars(c("PC1", "PC2", "PC3")), tolower) # for neatness

# create vector of snp variables to apply wrapper over
snp_list <- names(dat)[grepl("X\\d{1,2}", names(dat))]
snp_list

sapply(snp_list, function(x) glm_wrapper(dat, exposure = exposure, snp = x, covariates = covariates))


# --------- one off tests ---------- #
# glm_wrapper(dat, exposure = exposure, snp = 'X2.44319250.A.T', covariates = c("age_ref_imp", "study_gxe", "pc1", "pc2", "pc3"))

# make sure it's odds ratios.. 
# outout <- glm(outcome ~ hrt_ref_pm2 + age_ref_imp, data = dat, family = 'binomial')

# newdata = data.frame(age_ref_imp = mean(dat$age_ref_imp), hrt_ref_pm2 = 0, X18.10396692.A.G = 0)
# exp(predict(outout, newdata, type = 'link'))





